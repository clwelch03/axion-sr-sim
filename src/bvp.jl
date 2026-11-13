"""
Discretisation of the coupled radial boundary-value problem and its Newton
solver.

The unknowns are the node values of every channel's profile `R` and tortoise
slope `Q`, followed by a parameter block holding the complex mode frequencies,
the leading horizon coefficients of every channel, and the leading outer
coefficients of the channels not fixed by the prescribed amplitudes. MIRK4
relations on each mesh interval and series-matching conditions at both
endpoints form a square nonlinear system solved by Newton's method with a
structured sparse Jacobian and a backtracking line search.

All profile-like unknowns (`R`, `Q`, leading coefficients) are stored in units
of a reference `scale` (normally the largest amplitude) so that the Newton
tolerances are relative to the cloud size.
"""

# --- layout -------------------------------------------------------------------

"""
    UnknownLayout(problem)

Index bookkeeping for the Newton unknown vector: `n_nodes` node blocks of
`4 n_channels` reals (`Re R, Im R, Re Q, Im Q` per channel) followed by the
parameter block (mode frequencies, horizon leading coefficients, free outer
leading coefficients; each complex number as two reals).
"""
struct UnknownLayout
    n_nodes::Int
    n_channels::Int
    n_modes::Int
    fundamental_channels::Vector{Int}
    free_outer_channels::Vector{Int}
end

function UnknownLayout(problem::RadialProblem)
    fundamentals = fundamental_channels(problem.mode_set)
    free_outer = setdiff(1:n_channels(problem), fundamentals)
    UnknownLayout(n_nodes(problem.mesh), n_channels(problem), n_modes(problem), fundamentals, free_outer)
end

node_block_size(layout::UnknownLayout) = 4layout.n_channels
n_parameters(layout::UnknownLayout) = 2(layout.n_modes + layout.n_channels + length(layout.free_outer_channels))
n_unknowns(layout::UnknownLayout) = layout.n_nodes * node_block_size(layout) + n_parameters(layout)
node_range(layout::UnknownLayout, node::Int) = (node - 1) * node_block_size(layout) .+ (1:node_block_size(layout))
parameter_range(layout::UnknownLayout) = layout.n_nodes * node_block_size(layout) .+ (1:n_parameters(layout))
frequency_range(layout::UnknownLayout) = 1:2layout.n_modes   # within the parameter block

# Residual rows: one block per interval, then the horizon block, then the outer block.
interval_rows(layout::UnknownLayout, interval::Int) = node_range(layout, interval)
horizon_rows(layout::UnknownLayout) = node_range(layout, layout.n_nodes)
outer_rows(layout::UnknownLayout) = node_range(layout, layout.n_nodes + 1)

# --- packing between real vectors and complex quantities ----------------------

"""Interleaved real pairs → complex vector."""
_complexify(reals::AbstractVector) = [complex(reals[2i - 1], reals[2i]) for i in 1:length(reals) ÷ 2]

"""Complex vector → interleaved real pairs."""
_realify(complexes::AbstractVector) = collect(Iterators.flatten((real(z), imag(z)) for z in complexes))

"""
    NodeState(state_block)

Profiles `R` and slopes `Q` of all channels at one node, decoded from the
`4 n_channels` reals of a node block.
"""
struct NodeState{T<:Number}
    profiles::Vector{T}
    slopes::Vector{T}
end

"""Node state from profiles and slopes of possibly different numeric types (promoted)."""
function NodeState(profiles::AbstractVector{<:Number}, slopes::AbstractVector{<:Number})
    T = promote_type(eltype(profiles), eltype(slopes))
    NodeState{T}(convert(Vector{T}, profiles), convert(Vector{T}, slopes))
end

function NodeState(block::AbstractVector)
    n_chan = length(block) ÷ 4
    NodeState([complex(block[4c - 3], block[4c - 2]) for c in 1:n_chan],
              [complex(block[4c - 1], block[4c]) for c in 1:n_chan])
end

node_block(state::NodeState) =
    collect(Iterators.flatten((real(R), imag(R), real(Q), imag(Q)) for (R, Q) in zip(state.profiles, state.slopes)))

"""
    RadialParameters(layout, parameter_block, amplitudes, scale)

Decoded parameter block: complex mode frequencies, leading horizon coefficient
of every channel, and leading outer coefficient of every channel with the
prescribed amplitudes (in units of `scale`) filled into the fundamental slots.
"""
struct RadialParameters{T<:Number}
    mode_frequencies::Vector{T}
    horizon_leading::Vector{T}
    outer_leading::Vector{T}
end

function RadialParameters(layout::UnknownLayout, block::AbstractVector, amplitudes::AbstractVector, scale)
    n_modes, n_chan = layout.n_modes, layout.n_channels
    mode_frequencies = _complexify(block[1:2n_modes])
    horizon_leading = _complexify(block[2n_modes .+ (1:2n_chan)])
    free_outer = _complexify(block[2n_modes + 2n_chan + 1:end])
    outer_leading = zeros(Complex{promote_type(real(eltype(mode_frequencies)), typeof(scale), eltype(amplitudes))}, n_chan)
    outer_leading[layout.fundamental_channels] = amplitudes ./ scale
    outer_leading[layout.free_outer_channels] = free_outer
    T = promote_type(eltype(mode_frequencies), eltype(outer_leading))
    RadialParameters(Vector{T}(mode_frequencies), Vector{T}(horizon_leading), Vector{T}(outer_leading))
end

# --- workspaces and threading -------------------------------------------------

"""
    NodeWorkspace{T}(problem)

Scratch arrays for evaluating the right-hand side and its derivatives at one
radial point and one interval's MIRK4 relation: an `RhsWorkspace`, the decoded
complex profiles and slopes and their derivatives, the Hermite midpoint state,
its right-hand side and derivatives, and the chain-rule factors of the
midpoint with respect to the left node, the right node and the frequencies.
The element type `T` may be a dual number (residual only). One workspace
serves one task at a time.
"""
struct NodeWorkspace{T<:Number, R<:RhsWorkspace{T}}
    rhs::R
    profiles::Vector{Complex{T}}
    slopes::Vector{Complex{T}}
    profile_rates::Vector{Complex{T}}
    slope_rates::Vector{Complex{T}}
    midpoint::Vector{T}
    midpoint_rhs::Vector{T}
    midpoint_wrt_state::Matrix{T}
    midpoint_wrt_frequencies::Matrix{T}
    chain_state::Matrix{T}
    chain_frequencies::Matrix{T}
end

function NodeWorkspace{T}(problem::RadialProblem) where {T}
    n_chan, n_frequency_reals = n_channels(problem), 2n_modes(problem)
    n_state = 4n_chan
    complex_buffer() = zeros(Complex{T}, n_chan)
    rhs = RhsWorkspace{T}(problem)
    NodeWorkspace{T, typeof(rhs)}(rhs, complex_buffer(), complex_buffer(), complex_buffer(), complex_buffer(),
                                  zeros(T, n_state), zeros(T, n_state), zeros(T, n_state, n_state),
                                  zeros(T, n_state, n_frequency_reals), zeros(T, n_state, n_state),
                                  zeros(T, n_state, n_frequency_reals))
end

"""
    _foreach_in_chunks(body, indices, make_workspace; chunks_per_thread=4)

Call `body(workspace, index)` for every index. The indices are split into
contiguous chunks (`chunks_per_thread` per Julia thread, so that threads busy
with other tasks do not hold up the rest), run as parallel tasks; each task
gets its own workspace from `make_workspace()`. `body` must write only to
storage owned by its index, so the result does not depend on the number of
threads.
"""
function _foreach_in_chunks(body, indices, make_workspace; chunks_per_thread=4)
    n_tasks = Threads.nthreads() == 1 ? 1 : min(chunks_per_thread * Threads.nthreads(), length(indices))
    if n_tasks ≤ 1
        _run_chunk(body, indices, make_workspace)
    else
        chunks = Iterators.partition(indices, cld(length(indices), n_tasks))
        foreach(wait, [Threads.@spawn(_run_chunk(body, chunk, make_workspace)) for chunk in chunks])
    end
    return nothing
end

"""Run `body(workspace, index)` over one chunk of indices with a fresh workspace (a function of its own, so that tasks share no local variables)."""
function _run_chunk(body, chunk, make_workspace)
    workspace = make_workspace()
    for index in chunk
        body(workspace, index)
    end
    return nothing
end

"""Call `body(index)` for every index, one task each when several threads are available (for few, uneven jobs)."""
function _foreach_task(body, indices)
    if Threads.nthreads() == 1
        foreach(body, indices)
    else
        foreach(wait, [Threads.@spawn(body(index)) for index in indices])
    end
    return nothing
end

# --- residual pieces ----------------------------------------------------------

"""Number type of the assembly arithmetic for these inputs (dual numbers propagate)."""
_assembly_eltype(problem::RadialProblem, arguments...) =
    promote_type(map(argument -> real(eltype(argument)), arguments)..., typeof(problem.mass_coupling))

"""Decode a node block into the complex profiles and slopes of `workspace`."""
function _decode_node!(workspace::NodeWorkspace, block::AbstractVector)
    for channel in eachindex(workspace.profiles)
        workspace.profiles[channel] = complex(block[4channel - 3], block[4channel - 2])
        workspace.slopes[channel] = complex(block[4channel - 1], block[4channel])
    end
    return workspace
end

"""
    rhs_block!(result, workspace, problem, s, block, component_frequencies, scale)

Right-hand side at offset `s` as a real node block, written into `result`;
`component_frequencies` are the complex `ω̃_q` of every component.
"""
function rhs_block!(result::AbstractVector, workspace::NodeWorkspace, problem::RadialProblem, s, block::AbstractVector,
                    component_frequencies::AbstractVector, scale)
    _decode_node!(workspace, block)
    radial_rhs!(workspace.profile_rates, workspace.slope_rates, workspace.rhs, problem, s, workspace.profiles,
                workspace.slopes, component_frequencies, scale)
    for channel in eachindex(workspace.profiles)
        profile_rate, slope_rate = workspace.profile_rates[channel], workspace.slope_rates[channel]
        result[4channel - 3], result[4channel - 2] = real(profile_rate), imag(profile_rate)
        result[4channel - 1], result[4channel] = real(slope_rate), imag(slope_rate)
    end
    return result
end

"""Right-hand side as a real node block, with the mode frequencies given as interleaved reals."""
function rhs_block(problem::RadialProblem, s, block::AbstractVector, frequency_reals::AbstractVector, scale)
    T = _assembly_eltype(problem, block, frequency_reals, scale, s)
    frequencies = component_frequencies(problem.mode_set, _complexify(frequency_reals))
    return rhs_block!(zeros(T, length(block)), NodeWorkspace{T}(problem), problem, s, block, frequencies, scale)
end

"""
    interval_residual!(result, workspace, problem, interval, left, right, f_left, f_right, component_frequencies, scale)

MIRK4 (Simpson with Hermite midpoint) relation on one mesh interval of width `h`,
given the node states and the right-hand sides already evaluated there,
written into `result`:

    y_R − y_L − h/6 (f_L + 4 f_mid + f_R),   y_mid = (y_L + y_R)/2 − h/8 (f_R − f_L).
"""
function interval_residual!(result::AbstractVector, workspace::NodeWorkspace, problem::RadialProblem, interval::Int,
                            left::AbstractVector, right::AbstractVector, f_left::AbstractVector, f_right::AbstractVector,
                            component_frequencies::AbstractVector, scale)
    mesh = problem.mesh
    width = mesh.rstar[interval + 1] - mesh.rstar[interval]
    workspace.midpoint .= (left .+ right) ./ 2 .- (width / 8) .* (f_right .- f_left)
    f_mid = rhs_block!(workspace.midpoint_rhs, workspace, problem, mesh.midpoint_offsets[interval], workspace.midpoint,
                       component_frequencies, scale)
    result .= right .- left .- (width / 6) .* (f_left .+ 4 .* f_mid .+ f_right)
    return result
end

"""
    expected_boundary_values(problem, side, component_index, frequency, leading, offset)

Asymptotic profiles `R` and tortoise slopes `Q` of one component's channels at
the endpoint offset, from its complex frequency and leading coefficients: the
Frobenius series at the horizon, and at the outer end the truncated series or,
for propagating components with `outer_boundary.data == :integrated`, the far
series integrated inward (`integrated_outer_values`). `side` is `:horizon` or
`:outer`.
"""
function expected_boundary_values(problem::RadialProblem, side::Symbol, component_index::Int, frequency,
                                  leading::AbstractVector, offset)
    background, mass_coupling = problem.background, problem.mass_coupling
    component, mixing = problem.mode_set.components[component_index], problem.mixing[component_index]
    if side == :horizon
        series = HorizonSeries(background, mass_coupling, frequency, component.m, component.ls, mixing, leading,
                               problem.horizon_order)
        return endpoint_values(series, background, offset)
    elseif uses_integrated_data(problem.outer_boundary, frequency, mass_coupling)
        integrated = integrated_outer_values(background, mass_coupling, frequency, component.m, component.ls, mixing,
                                             leading, problem.outer_boundary, offset)
        return integrated.profiles, integrated.slopes
    else
        series = OuterSeries(background, mass_coupling, frequency, component.m, component.ls, mixing, leading,
                             problem.outer_order)
        return endpoint_values(series, background, offset)
    end
end

"""Mesh offset of an endpoint (`:horizon` or `:outer`)."""
endpoint_offset(problem::RadialProblem, side::Symbol) = side == :horizon ? problem.mesh.offsets[1] : problem.mesh.offsets[end]

"""
    boundary_residual(problem, layout, side, block, parameter_block, amplitudes, scale)

Mismatch between the node values of `R`, `Q` at an endpoint and the asymptotic
data evaluated there with the current frequencies and leading coefficients
(`expected_boundary_values`, one task per component). `side` is `:horizon` or
`:outer`.
"""
function boundary_residual(problem::RadialProblem, layout::UnknownLayout, side::Symbol, block::AbstractVector,
                           parameter_block::AbstractVector, amplitudes::AbstractVector, scale)
    mode_set = problem.mode_set
    parameters = RadialParameters(layout, parameter_block, amplitudes, scale)
    frequencies = component_frequencies(mode_set, parameters.mode_frequencies)
    leading = side == :horizon ? parameters.horizon_leading : parameters.outer_leading
    state = NodeState(block)
    offset = endpoint_offset(problem, side)

    value_type = promote_type(eltype(state.profiles), eltype(frequencies), Complex{typeof(problem.mass_coupling)})
    expected = NodeState(similar(state.profiles, value_type), similar(state.slopes, value_type))
    _foreach_task(eachindex(mode_set.components)) do component_index
        channels = mode_set.component_channels[component_index]
        expected.profiles[channels], expected.slopes[channels] =
            expected_boundary_values(problem, side, component_index, frequencies[component_index], leading[channels], offset)
    end
    return node_block(NodeState(state.profiles .- expected.profiles, state.slopes .- expected.slopes))
end

"""
    residual!(result, problem, layout, unknowns, amplitudes, scale)

Full residual vector of the discretised boundary-value problem, written into
`result`. Node right-hand sides and interval relations are evaluated in
parallel over the Julia threads (`_foreach_in_chunks`), concurrently with the
boundary rows; every entry is computed by the same arithmetic whatever the
thread count. Works for dual-number
unknowns, amplitudes or problems (background, mass coupling, mesh offsets).
"""
function residual!(result::AbstractVector, problem::RadialProblem, layout::UnknownLayout, unknowns::AbstractVector,
                   amplitudes::AbstractVector, scale)
    T = eltype(result)
    parameter_block = unknowns[parameter_range(layout)]
    frequencies = component_frequencies(problem.mode_set, _complexify(parameter_block[frequency_range(layout)]))
    node_rhs = zeros(T, node_block_size(layout), layout.n_nodes)
    # The boundary rows run as a separate task, concurrently with the node and interval passes.
    boundary_task = Threads.@spawn for (side, node, rows) in ((:horizon, 1, horizon_rows(layout)),
                                                              (:outer, layout.n_nodes, outer_rows(layout)))
        result[rows] = boundary_residual(problem, layout, side, unknowns[node_range(layout, node)], parameter_block,
                                         amplitudes, scale)
    end
    make_workspace() = NodeWorkspace{T}(problem)
    _foreach_in_chunks(1:layout.n_nodes, make_workspace) do workspace, node
        rhs_block!(view(node_rhs, :, node), workspace, problem, problem.mesh.offsets[node],
                   view(unknowns, node_range(layout, node)), frequencies, scale)
    end
    _foreach_in_chunks(1:layout.n_nodes-1, make_workspace) do workspace, interval
        interval_residual!(view(result, interval_rows(layout, interval)), workspace, problem, interval,
                           view(unknowns, node_range(layout, interval)), view(unknowns, node_range(layout, interval + 1)),
                           view(node_rhs, :, interval), view(node_rhs, :, interval + 1), frequencies, scale)
    end
    wait(boundary_task)
    return result
end

"""
    residual(problem, layout, unknowns, amplitudes, scale)

Full residual vector of the discretised boundary-value problem (see `residual!`).
"""
function residual(problem::RadialProblem, layout::UnknownLayout, unknowns::AbstractVector,
                  amplitudes::AbstractVector, scale)
    result = zeros(promote_type(eltype(unknowns), eltype(amplitudes), typeof(problem.mass_coupling)), n_unknowns(layout))
    return residual!(result, problem, layout, unknowns, amplitudes, scale)
end

# --- Jacobian -----------------------------------------------------------------

"""
    RhsDerivatives

Jacobian of the right-hand side at one radial point with respect to the node
block and to the (interleaved real) mode frequencies.
"""
struct RhsDerivatives{T<:Real}
    wrt_state::Matrix{T}
    wrt_frequencies::Matrix{T}
end

"""Reference implementation of `rhs_derivatives` by forward-mode automatic differentiation (used in tests)."""
function rhs_derivatives_autodiff(problem::RadialProblem, s, block::AbstractVector, frequency_reals::AbstractVector, scale)
    n_state = length(block)
    combined = ForwardDiff.jacobian(vcat(block, frequency_reals)) do inputs
        rhs_block(problem, s, inputs[1:n_state], inputs[n_state+1:end], scale)
    end
    RhsDerivatives(combined[:, 1:n_state], combined[:, n_state+1:end])
end

"""Write the real `2×2` block of the complex-linear map `z ↦ factor·z` at `(rows, columns)` of `matrix`, accumulating."""
function _add_complex_block!(matrix::AbstractMatrix, row::Int, column::Int, factor::Number)
    matrix[row, column] += real(factor)
    matrix[row, column + 1] -= imag(factor)
    matrix[row + 1, column] += imag(factor)
    matrix[row + 1, column + 1] += real(factor)
    return matrix
end

"""
    rhs_derivatives!(wrt_state, wrt_frequencies, workspace, problem, s, block, component_frequencies, scale)

Analytic Jacobian of `rhs_block` at offset `s`, written into `wrt_state`
(node block) and `wrt_frequencies` (interleaved real mode frequencies). The
right-hand side is complex-linear in the node state except for the nonlinear
source, whose real Jacobian comes from `nonlinear_source_jacobian!`; the
frequency dependence enters only through the linear coefficients `ΔB_{l,q}`
and `a²(ω̃_q² − α²)`.
"""
function rhs_derivatives!(wrt_state::AbstractMatrix, wrt_frequencies::AbstractMatrix, workspace::NodeWorkspace,
                          problem::RadialProblem, s, block::AbstractVector, component_frequencies::AbstractVector, scale)
    background, mode_set, mass_coupling = problem.background, problem.mode_set, problem.mass_coupling
    _decode_node!(workspace, block)
    profiles = workspace.profiles
    a = background.a
    r = radius(background, s)
    horizon_function, sum_r2a2 = delta(background, s), r2a2(background, s)
    n_chan = n_channels(mode_set)
    fill!(wrt_state, 0)
    fill!(wrt_frequencies, 0)

    slope_coefficient = -2r * horizon_function / sum_r2a2^2
    for channel in 1:n_chan
        profile_row, slope_row = 4channel - 3, 4channel - 1
        wrt_state[profile_row, slope_row] = 1
        wrt_state[profile_row + 1, slope_row + 1] = 1
        wrt_state[slope_row, slope_row] = slope_coefficient
        wrt_state[slope_row + 1, slope_row + 1] = slope_coefficient
    end

    for (component_index, component) in enumerate(mode_set.components)
        channels = mode_set.component_channels[component_index]
        frequency = component_frequencies[component_index]
        m = component.m
        mixing = problem.mixing[component_index]
        coupling_strength = a^2 * (frequency^2 - mass_coupling^2)
        for (local_index, channel) in enumerate(channels)
            l = component.ls[local_index]
            coupled = sum(mixing[local_index, other_local] * profiles[other_channel]
                          for (other_local, other_channel) in enumerate(channels))
            delta_B = (sum_r2a2 * frequency - a * m)^2 +
                      horizon_function * (2a * m * frequency - a^2 * frequency^2 - mass_coupling^2 * r^2 - l * (l + 1))
            for (other_local, other_channel) in enumerate(channels)
                factor = horizon_function * coupling_strength * mixing[local_index, other_local]
                other_channel == channel && (factor += delta_B)
                _add_complex_block!(wrt_state, 4channel - 1, 4other_channel - 3, -factor / sum_r2a2^2)
            end
            # ∂(forcing)/∂ω̃_q, then chain through ω̃_q = Σ_i q_i ω̃_i.
            forcing_wrt_frequency = (2sum_r2a2 * (sum_r2a2 * frequency - a * m) + horizon_function * (2a * m - 2a^2 * frequency)) *
                                    profiles[channel] + horizon_function * 2a^2 * frequency * coupled
            for mode_index in eachindex(component.q)
                wrt_real_part = -forcing_wrt_frequency * component.q[mode_index] / sum_r2a2^2
                wrt_imag_part = im * wrt_real_part
                wrt_frequencies[4channel - 1, 2mode_index - 1] = real(wrt_real_part)
                wrt_frequencies[4channel, 2mode_index - 1] = imag(wrt_real_part)
                wrt_frequencies[4channel - 1, 2mode_index] = real(wrt_imag_part)
                wrt_frequencies[4channel, 2mode_index] = imag(wrt_imag_part)
            end
        end
    end

    rhs = workspace.rhs
    rhs.physical_profiles .= scale .* profiles
    source_jacobian = nonlinear_source_jacobian!(rhs.source_jacobian, rhs.field, problem.quadrature, problem.potential,
                                                 background, mass_coupling, s, rhs.physical_profiles)
    source_factor = -horizon_function / sum_r2a2^2
    for channel in 1:n_chan, other_channel in 1:n_chan
        wrt_real_part = source_factor * source_jacobian.wrt_real[channel, other_channel]
        wrt_imag_part = source_factor * source_jacobian.wrt_imag[channel, other_channel]
        wrt_state[4channel - 1, 4other_channel - 3] += real(wrt_real_part)
        wrt_state[4channel, 4other_channel - 3] += imag(wrt_real_part)
        wrt_state[4channel - 1, 4other_channel - 2] += real(wrt_imag_part)
        wrt_state[4channel, 4other_channel - 2] += imag(wrt_imag_part)
    end
    return wrt_state, wrt_frequencies
end

"""
    rhs_derivatives(problem, s, block, frequency_reals, scale)

`RhsDerivatives` of `rhs_block` (see `rhs_derivatives!`), with the mode
frequencies given as interleaved reals.
"""
function rhs_derivatives(problem::RadialProblem, s, block::AbstractVector, frequency_reals::AbstractVector, scale)
    T = _assembly_eltype(problem, block, frequency_reals, scale, s)
    workspace = NodeWorkspace{T}(problem)
    frequencies = component_frequencies(problem.mode_set, _complexify(frequency_reals))
    wrt_state, wrt_frequencies = rhs_derivatives!(zeros(T, length(block), length(block)),
                                                  zeros(T, length(block), length(frequency_reals)),
                                                  workspace, problem, s, block, frequencies, scale)
    return RhsDerivatives(wrt_state, wrt_frequencies)
end

"""Add `value` to the diagonal of a square matrix."""
function _add_to_diagonal!(matrix::AbstractMatrix, value)
    for index in axes(matrix, 1)
        matrix[index, index] += value
    end
    return matrix
end

"""
    interval_jacobian!(wrt_left, wrt_right, wrt_frequencies, workspace, problem, interval,
                       left_derivatives, right_derivatives, left, right, f_left, f_right, component_frequencies, scale)

Blocks of one interval's MIRK4 residual with respect to the left node, the
right node and the mode frequencies, written into the three matrices, obtained
by chaining the right-hand-side derivatives through the Hermite midpoint.
`left_derivatives` and `right_derivatives` are `(wrt_state, wrt_frequencies)`
pairs of the right-hand side at the two nodes (`rhs_derivatives!`).
"""
function interval_jacobian!(wrt_left::AbstractMatrix, wrt_right::AbstractMatrix, wrt_frequencies::AbstractMatrix,
                            workspace::NodeWorkspace, problem::RadialProblem, interval::Int,
                            left_derivatives::Tuple, right_derivatives::Tuple, left::AbstractVector, right::AbstractVector,
                            f_left::AbstractVector, f_right::AbstractVector, component_frequencies::AbstractVector, scale)
    mesh = problem.mesh
    width = mesh.rstar[interval + 1] - mesh.rstar[interval]
    left_wrt_state, left_wrt_frequencies = left_derivatives
    right_wrt_state, right_wrt_frequencies = right_derivatives
    workspace.midpoint .= (left .+ right) ./ 2 .- (width / 8) .* (f_right .- f_left)
    mid_wrt_state, mid_wrt_frequencies = rhs_derivatives!(workspace.midpoint_wrt_state, workspace.midpoint_wrt_frequencies,
                                                          workspace, problem, mesh.midpoint_offsets[interval],
                                                          workspace.midpoint, component_frequencies, scale)
    chain = workspace.chain_state

    # ∂y_mid/∂y_L = I/2 + (h/8) ∂f_L/∂y_L
    chain .= (width / 8) .* left_wrt_state
    _add_to_diagonal!(chain, 1 / 2)
    mul!(wrt_left, mid_wrt_state, chain)
    wrt_left .= .-(width / 6) .* (left_wrt_state .+ 4 .* wrt_left)
    _add_to_diagonal!(wrt_left, -1)

    # ∂y_mid/∂y_R = I/2 − (h/8) ∂f_R/∂y_R
    chain .= .-(width / 8) .* right_wrt_state
    _add_to_diagonal!(chain, 1 / 2)
    mul!(wrt_right, mid_wrt_state, chain)
    wrt_right .= .-(width / 6) .* (right_wrt_state .+ 4 .* wrt_right)
    _add_to_diagonal!(wrt_right, 1)

    # ∂y_mid/∂ω = −(h/8) (∂f_R/∂ω − ∂f_L/∂ω)
    workspace.chain_frequencies .= .-(width / 8) .* (right_wrt_frequencies .- left_wrt_frequencies)
    mul!(wrt_frequencies, mid_wrt_state, workspace.chain_frequencies)
    wrt_frequencies .= .-(width / 6) .* (left_wrt_frequencies .+ right_wrt_frequencies .+ 4 .* (mid_wrt_frequencies .+ wrt_frequencies))
    return wrt_left, wrt_right, wrt_frequencies
end

"""
    boundary_jacobian!(wrt_parameters, problem, layout, side, parameter_block, amplitudes, scale)

Derivative of `boundary_residual` with respect to the parameter block, written
into `wrt_parameters` (`4 n_channels × n_parameters`); the derivative with
respect to the endpoint node block is the identity. Each component's expected
values depend only on its frequency `ω̃_q = q·ω̃` and its own leading
coefficients, so they are differentiated one component at a time (in parallel
tasks) with ForwardDiff over `(Re ω̃_q, Im ω̃_q)` and the real and imaginary
parts of those coefficients, and chained to the parameter columns. The
leading outer coefficients of fundamental channels are fixed by the
amplitudes and have no column.
"""
function boundary_jacobian!(wrt_parameters::AbstractMatrix, problem::RadialProblem, layout::UnknownLayout, side::Symbol,
                            parameter_block::AbstractVector, amplitudes::AbstractVector, scale)
    mode_set = problem.mode_set
    parameters = RadialParameters(layout, parameter_block, amplitudes, scale)
    frequencies = component_frequencies(mode_set, parameters.mode_frequencies)
    leading = side == :horizon ? parameters.horizon_leading : parameters.outer_leading
    offset = endpoint_offset(problem, side)
    leading_offset = side == :horizon ? 2layout.n_modes : 2layout.n_modes + 2layout.n_channels
    leading_index = side == :horizon ? collect(1:layout.n_channels) :
                    [something(findfirst(==(channel), layout.free_outer_channels), 0) for channel in 1:layout.n_channels]
    fill!(wrt_parameters, 0)
    _foreach_task(eachindex(mode_set.components)) do component_index
        channels = mode_set.component_channels[component_index]
        frequency = frequencies[component_index]
        inputs = vcat(real(frequency), imag(frequency), _realify(leading[channels]))
        derivatives = ForwardDiff.jacobian(inputs) do component_inputs
            profiles, slopes = expected_boundary_values(problem, side, component_index,
                                                        complex(component_inputs[1], component_inputs[2]),
                                                        _complexify(component_inputs[3:end]), offset)
            node_block(NodeState(profiles, slopes))
        end
        rows = 4first(channels) - 3:4last(channels)
        for (mode_index, weight) in enumerate(mode_set.components[component_index].q)
            wrt_parameters[rows, 2mode_index - 1] .-= weight .* view(derivatives, :, 1)
            wrt_parameters[rows, 2mode_index] .-= weight .* view(derivatives, :, 2)
        end
        for (local_index, channel) in enumerate(channels)
            leading_index[channel] == 0 && continue
            column = leading_offset + 2leading_index[channel] - 1
            wrt_parameters[rows, column] .-= view(derivatives, :, 2local_index + 1)
            wrt_parameters[rows, column + 1] .-= view(derivatives, :, 2local_index + 2)
        end
    end
    return wrt_parameters
end

"""
    JacobianBlocks

Dense blocks of the sparse Jacobian; every other entry is structurally zero.
The rows of interval `i` (`interval_rows(layout, i)`) have three blocks:
`wrt_left[:, :, i]` (columns of node `i`), `wrt_right[:, :, i]` (node `i + 1`)
and `wrt_frequencies[:, :, i]` (the `2 n_modes` mode-frequency columns at the
start of the parameter block). The horizon rows have `horizon_wrt_node`
(node 1, the identity) and `horizon_wrt_parameters` (all parameter columns);
the outer rows `outer_wrt_node` (last node, the identity) and
`outer_wrt_parameters`. With `T = Int` the same struct holds the position of
every entry in the `nzval` of the sparse matrix (`JacobianPattern`).
"""
struct JacobianBlocks{T<:Real}
    wrt_left::Array{T,3}
    wrt_right::Array{T,3}
    wrt_frequencies::Array{T,3}
    horizon_wrt_node::Matrix{T}
    horizon_wrt_parameters::Matrix{T}
    outer_wrt_node::Matrix{T}
    outer_wrt_parameters::Matrix{T}
end

_block_arrays(blocks::JacobianBlocks) = (blocks.wrt_left, blocks.wrt_right, blocks.wrt_frequencies, blocks.horizon_wrt_node,
                                         blocks.horizon_wrt_parameters, blocks.outer_wrt_node, blocks.outer_wrt_parameters)

"""Array shapes of the seven `JacobianBlocks` fields, in field order."""
function _block_shapes(layout::UnknownLayout)
    n_state, n_frequency_reals, n_param = node_block_size(layout), 2layout.n_modes, n_parameters(layout)
    n_intervals = layout.n_nodes - 1
    return ((n_state, n_state, n_intervals), (n_state, n_state, n_intervals), (n_state, n_frequency_reals, n_intervals),
            (n_state, n_state), (n_state, n_param), (n_state, n_state), (n_state, n_param))
end

"""
    JacobianWorkspace(problem, layout=UnknownLayout(problem))

Preallocated storage for Jacobian assembly: the `JacobianBlocks` and the
right-hand sides and their derivatives at every node. Used by
`jacobian_blocks!`; not safe to share between concurrent calls.
"""
struct JacobianWorkspace
    layout::UnknownLayout
    blocks::JacobianBlocks{Float64}
    node_rhs::Matrix{Float64}
    node_wrt_state::Array{Float64,3}
    node_wrt_frequencies::Array{Float64,3}
end

function JacobianWorkspace(problem::RadialProblem, layout::UnknownLayout=UnknownLayout(problem))
    n_state, n_frequency_reals, n_nodes = node_block_size(layout), 2layout.n_modes, layout.n_nodes
    blocks = JacobianBlocks{Float64}((zeros(shape) for shape in _block_shapes(layout))...)
    JacobianWorkspace(layout, blocks, zeros(n_state, n_nodes), zeros(n_state, n_state, n_nodes),
                      zeros(n_state, n_frequency_reals, n_nodes))
end

"""
    JacobianPattern(layout)

Sparsity pattern of the Jacobian, fixed by the layout: the column pointers and
row indices of a `SparseMatrixCSC` holding exactly the `JacobianBlocks`
entries, and the `nzval` position of every block entry (`JacobianBlocks{Int}`).
Used by `sparse_jacobian`; read-only after construction.
"""
struct JacobianPattern
    size::Int
    colptr::Vector{Int}
    rowval::Vector{Int}
    positions::JacobianBlocks{Int}
end

function JacobianPattern(layout::UnknownLayout)
    n_intervals, n_nodes = layout.n_nodes - 1, layout.n_nodes
    shapes = _block_shapes(layout)
    parameter_columns = parameter_range(layout)
    frequency_columns = parameter_columns[frequency_range(layout)]
    # Row and column of every block entry, in the memory order of the block arrays (column-major, interval last).
    rows, columns = Int[], Int[]
    function add_block!(row_indices, column_indices)
        for column in column_indices, row in row_indices
            push!(rows, row)
            push!(columns, column)
        end
    end
    foreach(interval -> add_block!(interval_rows(layout, interval), node_range(layout, interval)), 1:n_intervals)
    foreach(interval -> add_block!(interval_rows(layout, interval), node_range(layout, interval + 1)), 1:n_intervals)
    foreach(interval -> add_block!(interval_rows(layout, interval), frequency_columns), 1:n_intervals)
    add_block!(horizon_rows(layout), node_range(layout, 1))
    add_block!(horizon_rows(layout), parameter_columns)
    add_block!(outer_rows(layout), node_range(layout, n_nodes))
    add_block!(outer_rows(layout), parameter_columns)

    total, n_entries = n_unknowns(layout), length(rows)
    numbered = sparse(rows, columns, 1:n_entries, total, total)   # nzval[slot] = entry number
    nnz(numbered) == n_entries || error("JacobianPattern: overlapping Jacobian blocks")
    position = zeros(Int, n_entries)
    position[numbered.nzval] = 1:n_entries
    ends = cumsum(collect(prod.(shapes)))
    starts = ends .- prod.(shapes) .+ 1
    positions = JacobianBlocks{Int}((reshape(position[starts[index]:ends[index]], shapes[index]) for index in eachindex(shapes))...)
    JacobianPattern(total, numbered.colptr, numbered.rowval, positions)
end

"""
    jacobian_blocks!(workspace, problem, layout, unknowns, amplitudes, scale)

Fill `workspace.blocks` (`JacobianBlocks`) with the Jacobian of `residual` at
`unknowns` and return them. Node derivatives and interval blocks are computed
in parallel over the Julia threads, concurrently with the boundary blocks (one
task per component); each block entry is computed by the same arithmetic
whatever the thread count.
"""
function jacobian_blocks!(workspace::JacobianWorkspace, problem::RadialProblem, layout::UnknownLayout,
                          unknowns::AbstractVector, amplitudes::AbstractVector, scale)
    mesh, blocks = problem.mesh, workspace.blocks
    parameter_block = unknowns[parameter_range(layout)]
    frequencies = component_frequencies(problem.mode_set, _complexify(parameter_block[frequency_range(layout)]))
    node_rhs, node_wrt_state, node_wrt_frequencies = workspace.node_rhs, workspace.node_wrt_state, workspace.node_wrt_frequencies
    # The boundary blocks run as a separate task, concurrently with the node and interval passes.
    boundary_task = Threads.@spawn for (side, wrt_node, wrt_parameters) in
                                       ((:horizon, blocks.horizon_wrt_node, blocks.horizon_wrt_parameters),
                                        (:outer, blocks.outer_wrt_node, blocks.outer_wrt_parameters))
        fill!(wrt_node, 0)
        _add_to_diagonal!(wrt_node, 1)
        boundary_jacobian!(wrt_parameters, problem, layout, side, parameter_block, amplitudes, scale)
    end
    make_workspace() = NodeWorkspace{Float64}(problem)
    _foreach_in_chunks(1:layout.n_nodes, make_workspace) do node_workspace, node
        block = view(unknowns, node_range(layout, node))
        rhs_block!(view(node_rhs, :, node), node_workspace, problem, mesh.offsets[node], block, frequencies, scale)
        rhs_derivatives!(view(node_wrt_state, :, :, node), view(node_wrt_frequencies, :, :, node), node_workspace, problem,
                         mesh.offsets[node], block, frequencies, scale)
    end
    _foreach_in_chunks(1:layout.n_nodes-1, make_workspace) do node_workspace, interval
        left, right = interval, interval + 1
        interval_jacobian!(view(blocks.wrt_left, :, :, interval), view(blocks.wrt_right, :, :, interval),
                           view(blocks.wrt_frequencies, :, :, interval), node_workspace, problem, interval,
                           (view(node_wrt_state, :, :, left), view(node_wrt_frequencies, :, :, left)),
                           (view(node_wrt_state, :, :, right), view(node_wrt_frequencies, :, :, right)),
                           view(unknowns, node_range(layout, left)), view(unknowns, node_range(layout, right)),
                           view(node_rhs, :, left), view(node_rhs, :, right), frequencies, scale)
    end
    wait(boundary_task)
    return blocks
end

"""
    sparse_jacobian(blocks, pattern)

New `SparseMatrixCSC` with the pattern `pattern` (`JacobianPattern`) holding the
entries of `blocks` (`JacobianBlocks`).
"""
function sparse_jacobian(blocks::JacobianBlocks, pattern::JacobianPattern)
    values = zeros(length(pattern.rowval))
    for (block_values, block_positions) in zip(_block_arrays(blocks), _block_arrays(pattern.positions))
        for index in eachindex(block_values, block_positions)
            values[block_positions[index]] = block_values[index]
        end
    end
    return SparseMatrixCSC(pattern.size, pattern.size, copy(pattern.colptr), copy(pattern.rowval), values)
end

"""
    jacobian(problem, layout, unknowns, amplitudes, scale; workspace=JacobianWorkspace(problem, layout),
             pattern=JacobianPattern(layout))

Sparse Jacobian of `residual` with respect to `unknowns`: block-bidiagonal in
the node blocks, bordered by the dense parameter columns (see `JacobianBlocks`).
Returns a new matrix; pass a `workspace` and a `pattern` to reuse them between
calls.
"""
jacobian(problem::RadialProblem, layout::UnknownLayout, unknowns::AbstractVector, amplitudes::AbstractVector, scale;
         workspace::JacobianWorkspace=JacobianWorkspace(problem, layout), pattern::JacobianPattern=JacobianPattern(layout)) =
    sparse_jacobian(jacobian_blocks!(workspace, problem, layout, unknowns, amplitudes, scale), pattern)

# --- solution and Newton iteration -------------------------------------------

"""
    RadialSolution

Converged (or final) state of the radial solve: the problem, the prescribed
amplitudes and reference `scale`, the unknown vector, its decoded parameters,
the factorised Jacobian at the solution (reused for sensitivities), and
convergence diagnostics.
"""
struct RadialSolution{T<:Real, P<:RadialProblem, F}
    problem::P
    layout::UnknownLayout
    amplitudes::Vector{T}
    scale::T
    unknowns::Vector{T}
    parameters::RadialParameters{Complex{T}}
    jacobian_factorization::F
    iterations::Int
    residual_norm::T
    converged::Bool
end

"""Physical profiles `R_{l,q}` at every node, as a matrix `[node, channel]`."""
profiles(solution::RadialSolution) = _node_matrix(solution, :profiles)

"""Physical tortoise slopes `Q_{l,q} = dR/dr*` at every node, as a matrix `[node, channel]`."""
slopes(solution::RadialSolution) = _node_matrix(solution, :slopes)

function _node_matrix(solution::RadialSolution, field::Symbol)
    layout = solution.layout
    stack = [getfield(NodeState(solution.unknowns[node_range(layout, node)]), field) for node in 1:layout.n_nodes]
    return solution.scale .* permutedims(reduce(hcat, stack))
end

"""Complex mode frequencies `ω̃_i = ω_i + iγ_i`."""
mode_frequencies(solution::RadialSolution) = solution.parameters.mode_frequencies

"""Growth rates `γ_i = Im ω̃_i` of the bound modes."""
growth_rates(solution::RadialSolution) = imag.(mode_frequencies(solution))

# --- envelope scaling of the Newton system ----------------------------------------

"""
    channel_envelopes(problem, layout, mode_frequencies)

Matrix `[node, channel]` of the expected size of each channel along the mesh,
used to scale the Newton system. Bound components (`Re ω_q < α`) decay like
`r^ν e^{−κr}` with `κ = √(α² − ω_q²)` and `ν = −1 + (2ω_q² − α²)/κ`; a driven
component whose nonlinear source decays more slowly (rate `Σ|q_i|κ_i` over the
fundamentals) follows the source instead. The envelope is one up to the
radius where `r^ν e^{−κr}` peaks and follows it beyond; propagating components
get one everywhere. The frequencies are those at the start of a solve and are
not differentiated.

The amplitude of each bound mode is imposed at the outer endpoint, where its
profile can be many orders of magnitude below its peak. Without scaling, that
row, the only one that sets the size of the solution, changes the residual by
less than the Newton tolerance on domains that contain the cloud, and the
solve drifts towards `R ≈ 0` (`validation/logs/domain_outer_debug_1000.log`).
"""
function channel_envelopes(problem::RadialProblem, layout::UnknownLayout, mode_frequencies::AbstractVector)
    mode_set, background, mass_coupling = problem.mode_set, problem.background, problem.mass_coupling
    real_frequencies = real.(component_frequencies(mode_set, mode_frequencies))
    decay_rate(frequency) = frequency < mass_coupling ? sqrt(mass_coupling^2 - frequency^2) : 0.0
    long_range_power(frequency) = -1 + (2frequency^2 - mass_coupling^2) / decay_rate(frequency)
    fundamental_decay = [decay_rate(real(frequency)) for frequency in mode_frequencies]
    fundamental_power = [long_range_power(real(frequency)) for frequency in mode_frequencies]
    radii = [radius(background, offset) for offset in problem.mesh.offsets]
    envelopes = ones(layout.n_nodes, layout.n_channels)
    for (component_index, component) in enumerate(mode_set.components)
        frequency = real_frequencies[component_index]
        frequency < mass_coupling || continue
        decay, power = decay_rate(frequency), long_range_power(frequency)
        source_decay = sum(abs.(component.q) .* fundamental_decay)
        if source_decay < decay
            decay, power = source_decay, sum(abs.(component.q) .* fundamental_power)
        end
        envelope_peak = max(power / decay, radii[1])
        log_envelope = [r > envelope_peak ? power * log(r / envelope_peak) - decay * (r - envelope_peak) : 0.0
                        for r in radii]
        # keep the weights representable; such tails are below underflow anyway
        envelopes[:, mode_set.component_channels[component_index]] .= exp.(max.(log_envelope, -600.0))
    end
    return envelopes
end

"""
    NewtonScaling(problem, layout, mode_frequencies)

Row and column weights of the scaled Newton system `(D_r J D_c) y = −D_r F`,
step `D_c y`. Each node block of a channel is measured in units of its envelope
(`channel_envelopes`); the residual rows of an interval and of the outer
boundary are divided by the envelope at the interval's inner node and at the
endpoint. Parameters and horizon rows are unscaled. The discrete solution is
unchanged; only the conditioning and the convergence test change.
"""
struct NewtonScaling
    row::Vector{Float64}
    column::Vector{Float64}
    layout::UnknownLayout
end

function NewtonScaling(problem::RadialProblem, layout::UnknownLayout, mode_frequencies::AbstractVector)
    envelopes = channel_envelopes(problem, layout, mode_frequencies)
    row, column = ones(n_unknowns(layout)), ones(n_unknowns(layout))
    for channel in 1:layout.n_channels
        reals = 4channel-3:4channel   # Re R, Im R, Re Q, Im Q of this channel within a node block
        for node in 1:layout.n_nodes
            column[node_range(layout, node)[reals]] .= envelopes[node, channel]
        end
        for interval in 1:layout.n_nodes-1
            row[interval_rows(layout, interval)[reals]] .= 1 / envelopes[interval, channel]
        end
        row[outer_rows(layout)[reals]] .= 1 / envelopes[end, channel]
    end
    return NewtonScaling(row, column, layout)
end

"""Factorization method used by `ScaledFactorization` when none is given (`:block` or `:sparse_lu`)."""
const DEFAULT_FACTORIZATION_METHOD = Ref(:block)

"""
    ScaledFactorization

LU factorization of the scaled Jacobian `D_r J D_c` with its weights;
`factorization \\ b` returns `J⁻¹ b` in unscaled variables, so the implicit-
function solves (`amplitude_response`, `background_charge_derivatives`) use the
well-conditioned system without knowing about the scaling.

`ScaledFactorization(jacobian, scaling; method)` factorises with the block
elimination of `BlockFactorization` (`method=:block`, which needs the Jacobian structure of
`scaling.layout`) or with the general sparse LU (`method=:sparse_lu`, UMFPACK, kept as the reference).
The default `method` is `DEFAULT_FACTORIZATION_METHOD[]`, which validation scripts set to compare the two.
"""
struct ScaledFactorization{F}
    factorization::F
    scaling::NewtonScaling
end

function ScaledFactorization(jacobian::AbstractMatrix, scaling::NewtonScaling;
                            method::Symbol=DEFAULT_FACTORIZATION_METHOD[])
    if method == :block
        sparse_matrix = jacobian isa SparseMatrixCSC ? jacobian : sparse(jacobian)   # `sparse` of a sparse matrix copies it
        factors = BlockFactorization(sparse_matrix, scaling.layout, scaling.row, scaling.column)
    elseif method == :sparse_lu
        factors = lu(Diagonal(scaling.row) * jacobian * Diagonal(scaling.column))
    else
        throw(ArgumentError("unknown factorization method $method; use :block or :sparse_lu"))
    end
    return ScaledFactorization(factors, scaling)
end

Base.:\(scaled::ScaledFactorization, right_hand_side::AbstractVecOrMat) =
    scaled.scaling.column .* (scaled.factorization \ (scaled.scaling.row .* right_hand_side))

"""
    solve_radial(problem, amplitudes; guess=nothing, guess_change=nothing, tol=1e-10, floor_tol=tol, maxiter=30,
                 verbose=false, reuse_factorization=true, chord_threshold=1e-3)

Solve the coupled radial boundary-value problem at the prescribed amplitudes.
`guess` may be a previous `RadialSolution` (warm start, rescaled to the new
amplitudes) or `nothing` for the hydrogenic initial guess. `guess_change`, in
units of `guess.scale`, is added to the guess's unknowns before rescaling, for
example a first-order prediction (`predicted_change`); if the predicted start
has a residual above `chord_threshold` and the plain warm start a smaller one,
the plain warm start is used instead. The Newton system
is solved with envelope scaling (`NewtonScaling`, built once per solve from the
starting frequencies), and iterations stop when both the scaled residual and
the scaled step fall below `tol` in the max norm, or when the residual is below
`tol` and a step below `√tol` no longer reduces it (its roundoff floor). On
domains that do not reach far into the tails the scaling is close to the identity.
A floor above `tol` stops the solve unconverged as soon as a Newton step below
`tol` fails to reduce the residual (`_floor_action`); with `floor_tol > tol`
(opt-in) a Newton step at a floor below `floor_tol` counts as converged. The floor is set by the outer-boundary rows, whose roundoff is
amplified by the envelope row weights: about 1e-14 to 1e-10 on the 690-node
α = 0.3 meshes, 5e-10 on the 1355-node α = 0.187 mesh.

With `reuse_factorization` and a warm start, the iteration begins as a chord
method: steps use the guess's factorised Jacobian (and later the most recent
one) as long as each step reduces the residual by at least a factor of two,
and the Jacobian is only re-assembled and refactorised when that fails; chord
steps are attempted only while the residual is below `chord_threshold`. Along
a time evolution, where successive amplitudes are close, this replaces most
of the factorisations by residual evaluations. The returned solution always
carries a factorisation taken at (or within `tol/10` of) the converged point,
as the amplitude sensitivities require; a solve at the guess's own amplitudes
and problem keeps the guess's factorisation.
"""
function solve_radial(problem::RadialProblem, amplitudes::AbstractVector; guess=nothing, guess_change=nothing,
                      tol=1e-10, floor_tol=tol, maxiter=30, verbose=false, reuse_factorization=true, chord_threshold=1e-3)
    layout = UnknownLayout(problem)
    scale = maximum(amplitudes)
    warm = guess isa RadialSolution && length(guess.unknowns) == n_unknowns(layout)
    predicted = warm && guess_change !== nothing && any(!iszero, guess_change)
    unknowns = warm ? warm_start_unknowns(guess, scale; change=guess_change) : initial_unknowns(problem, layout, amplitudes, scale)
    scaling = NewtonScaling(problem, layout, _complexify(unknowns[parameter_range(layout)][frequency_range(layout)]))
    assembly = nothing   # Jacobian workspace and pattern, built when the first Jacobian is needed and reused
    scaled_norm(residual_vector) = norm(scaling.row .* residual_vector, Inf)
    scaled_step_norm(step_vector) = norm(step_vector ./ scaling.column, Inf)
    current = residual(problem, layout, unknowns, amplitudes, scale)
    current_norm = scaled_norm(current)
    if predicted && current_norm ≥ chord_threshold
        # The prediction is not clearly close: start from the plain warm start if that is closer.
        plain = warm_start_unknowns(guess, scale)
        plain_residual = residual(problem, layout, plain, amplitudes, scale)
        if scaled_norm(plain_residual) < current_norm
            verbose && println("  prediction |F| = $current_norm worse than the plain warm start; using the latter")
            unknowns, current, current_norm = plain, plain_residual, scaled_norm(plain_residual)
        end
    end

    # Chord steps are only worthwhile near a solution; far from it every rejected step costs a residual.
    factorization = warm && reuse_factorization && current_norm < chord_threshold ? guess.jacobian_factorization : nothing
    # Max-norm distance from the current point to where the factorisation was taken. The guess's factorisation
    # belongs to the current system only if problem and amplitudes are the guess's own.
    same_system = warm && !predicted && amplitudes == guess.amplitudes && _same_discretisation(problem, guess.problem)
    factorization_distance = same_system ? 0.0 : Inf
    # Reference scale of the unknowns when the factorisation was taken (see `_rescale_frequency_step!`).
    factorization_scale = warm ? float(guess.scale) : float(scale)
    converged = false
    iterations = 0
    stalled = 0
    previous_step_norm = Inf
    for iteration in 1:maxiter
        iterations = iteration
        chord = factorization !== nothing
        if !chord
            assembly === nothing && (assembly = (JacobianWorkspace(problem, layout), JacobianPattern(layout)))
            factorization = _factorised_jacobian(assembly, problem, layout, unknowns, amplitudes, scale, scaling)
            factorization_distance = 0.0
            factorization_scale = float(scale)
        end
        step = _rescale_frequency_step!(-(factorization \ current), layout, scale / factorization_scale)
        # Backtracking line search on the residual norm (Armijo condition); a
        # chord step must at least halve the residual or it is abandoned.
        # A small full step that does not reduce the residual may mean the residual is at its roundoff floor
        # (`_floor_action`); backtracking would then only shrink the step (up to 14 residual evaluations).
        damping = 1.0
        trial, trial_residual, trial_norm = unknowns, current, current_norm
        floor_action = nothing
        while damping ≥ 1e-4
            trial = unknowns .+ damping .* step
            trial_residual = residual(problem, layout, trial, amplitudes, scale)
            trial_norm = scaled_norm(trial_residual)
            trial_norm ≤ (1 - 1e-4damping) * current_norm && break
            damping == 1 && (floor_action = _floor_action(chord, current_norm, scaled_step_norm(step), tol, floor_tol))
            floor_action === nothing || break
            damping /= 2
        end
        if floor_action == :refactorise
            verbose && println("  chord step $iteration does not reduce |F| = $current_norm; refactorising")
            factorization = nothing
            continue
        elseif floor_action == :converged
            verbose && println("  $(chord ? "chord" : "newton") $iteration: |F| = $current_norm is at its roundoff floor; converged")
            converged = true
            break
        elseif floor_action == :stuck
            verbose && println("  newton $iteration: |F| = $current_norm is at its roundoff floor, above tol; stopping")
            break
        end
        if chord && trial_norm > max(current_norm / 2, tol)
            verbose && println("  chord step $iteration rejected (|F| $current_norm → $trial_norm); refactorising")
            factorization = nothing
            continue
        end
        unknowns, current, current_norm = trial, trial_residual, trial_norm
        reuse_factorization && current_norm < chord_threshold || (factorization = nothing)
        factorization_distance += damping * norm(step, Inf)
        step_norm = damping * scaled_step_norm(step)
        verbose && println("  $(chord ? "chord" : "newton") $iteration: |step| = $step_norm  |F| = $current_norm  damping = $damping")
        if current_norm < tol && step_norm < tol
            converged = true
            break
        end
        # Residual below tol while the step stops shrinking. With envelope scaling the scaled step can sit at a
        # roundoff floor of ~1e-9 (tail unknowns measured in units of their envelopes): accept that as converged.
        # A step that stalls above √tol still means a (nearly) singular system. (Chord steps near convergence
        # shrink geometrically, which is fine.)
        if current_norm < tol && step_norm > previous_step_norm / 2
            stalled += 1
            if stalled ≥ 3
                if step_norm < sqrt(tol)
                    converged = true
                else
                    @warn "solve_radial: |F| = $current_norm is below tol but the Newton step is $step_norm; " *
                          "the system is ill-conditioned"
                end
                break
            end
        else
            stalled = 0
        end
        previous_step_norm = step_norm
    end
    converged || @warn "solve_radial did not converge: |F| = $current_norm after $iterations iterations"
    # Sensitivities need the Jacobian at the converged point: the relative rate errors are 0.1 to 6 × the distance
    # (validation/newton_jacobian_accuracy.jl), so refactorise unless the factorisation was taken within tol/10,
    # where they are at the noise floor of the converged solution itself.
    if factorization === nothing || factorization_distance > tol / 10
        assembly === nothing && (assembly = (JacobianWorkspace(problem, layout), JacobianPattern(layout)))
        factorization = _factorised_jacobian(assembly, problem, layout, unknowns, amplitudes, scale, scaling)
    end
    parameters = RadialParameters(layout, unknowns[parameter_range(layout)], amplitudes, scale)
    RadialSolution(problem, layout, collect(float.(amplitudes)), float(scale), unknowns, parameters, factorization,
                   iterations, current_norm, converged)
end

"""Scaled factorisation of the Jacobian at `unknowns`, assembled with `assembly = (workspace, pattern)`."""
function _factorised_jacobian(assembly::Tuple{JacobianWorkspace, JacobianPattern}, problem::RadialProblem, layout::UnknownLayout,
                              unknowns::AbstractVector, amplitudes::AbstractVector, scale, scaling::NewtonScaling)
    workspace, pattern = assembly
    blocks = jacobian_blocks!(workspace, problem, layout, unknowns, amplitudes, scale)
    return ScaledFactorization(sparse_jacobian(blocks, pattern), scaling)
end

"""
    _rescale_frequency_step!(step, layout, scale_ratio)

Correct a step computed with a Jacobian factorised when the unknowns had another reference scale. Profile-like
unknowns and residuals both scale like `1/scale` while the frequencies do not, so the Jacobian in units of
`scale` equals the one in units of `scale₀` with its frequency columns multiplied by `scale₀/scale`; the
frequency entries of the step are therefore multiplied by `scale_ratio = scale/scale₀`. Uncorrected, the
frequency columns are off by `|scale/scale₀ − 1|`, 20–30% for a nonlinear-regime spindown step, where the chord
iteration then reduced the residual by 0.38 per step instead of 0.11–0.18.
"""
function _rescale_frequency_step!(step::AbstractVector, layout::UnknownLayout, scale_ratio)
    scale_ratio == 1 || (step[parameter_range(layout)[frequency_range(layout)]] .*= scale_ratio)
    return step
end

"""
    _floor_action(chord, residual_norm, step_norm, tol, floor_tol)

What to do when a full step of scaled size `step_norm` does not reduce the residual, or `nothing` to backtrack
as usual. With the residual below `tol` and the step below `√tol` the residual is at its roundoff floor:
`:converged` (for a chord step only if the step is also below `tol`; a larger one may come from an out-of-date
Jacobian, so `:refactorise` and let a Newton step decide). A Newton step at a residual below the opt-in
`floor_tol` likewise converges. A Newton step below `tol` that does not reduce a residual above `tol` means a
floor above `tol`: `:stuck`, since repeating the step from the same point gives the same result.
"""
function _floor_action(chord::Bool, residual_norm, step_norm, tol, floor_tol)
    step_norm < sqrt(tol) || return nothing
    if chord
        residual_norm < tol && step_norm < tol && return :converged
        return :refactorise
    end
    (residual_norm < tol || residual_norm < floor_tol) && return :converged
    step_norm < tol && return :stuck
    return nothing
end

"""
Whether two problems give the same discrete system (same background, mass coupling, potential, modes, mesh,
quadrature and boundary treatment). Shared components are compared by identity, so problems built separately
from equal inputs may compare unequal; `solve_radial` then only refactorises once more than necessary.
"""
_same_discretisation(problem::RadialProblem, other::RadialProblem) =
    problem === other ||
    (problem.background == other.background && problem.mass_coupling == other.mass_coupling &&
     problem.potential == other.potential && problem.mode_set === other.mode_set &&
     problem.quadrature === other.quadrature && problem.mixing === other.mixing && problem.mesh.rstar == other.mesh.rstar &&
     problem.horizon_order == other.horizon_order && problem.outer_order == other.outer_order &&
     problem.outer_boundary == other.outer_boundary)

"""
    warm_start_unknowns(previous, scale; change=nothing)

A previous solution's unknowns, plus `change` if given (in units of `previous.scale`), rescaled to a new
reference `scale`; frequencies are not rescaled.
"""
function warm_start_unknowns(previous::RadialSolution, scale; change=nothing)
    moved = change === nothing ? previous.unknowns : previous.unknowns .+ change
    unknowns = moved .* (previous.scale / scale)
    layout = previous.layout
    frequency_indices = parameter_range(layout)[frequency_range(layout)]
    unknowns[frequency_indices] = moved[frequency_indices]
    return unknowns
end
