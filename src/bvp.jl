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

# --- residual pieces ----------------------------------------------------------

"""Right-hand side as a real node block, with the mode frequencies given as interleaved reals."""
function rhs_block(problem::RadialProblem, s, block::AbstractVector, frequency_reals::AbstractVector, scale)
    state = NodeState(block)
    frequencies = component_frequencies(problem.mode_set, _complexify(frequency_reals))
    dR, dQ = radial_rhs(problem, s, state.profiles, state.slopes, frequencies, scale)
    return node_block(NodeState(dR, dQ))
end

"""Right-hand side at every mesh node, as a vector of node blocks."""
node_rhs_blocks(problem::RadialProblem, layout::UnknownLayout, unknowns::AbstractVector, frequency_reals::AbstractVector, scale) =
    [rhs_block(problem, problem.mesh.offsets[node], unknowns[node_range(layout, node)], frequency_reals, scale)
     for node in 1:layout.n_nodes]

"""
    interval_residual(problem, interval, left, right, f_left, f_right, frequency_reals, scale)

MIRK4 (Simpson with Hermite midpoint) relation on one mesh interval of width `h`,
given the node states and the right-hand sides already evaluated there:

    y_R − y_L − h/6 (f_L + 4 f_mid + f_R),   y_mid = (y_L + y_R)/2 − h/8 (f_R − f_L).
"""
function interval_residual(problem::RadialProblem, interval::Int, left::AbstractVector, right::AbstractVector,
                           f_left::AbstractVector, f_right::AbstractVector, frequency_reals::AbstractVector, scale)
    mesh = problem.mesh
    width = mesh.rstar[interval + 1] - mesh.rstar[interval]
    midpoint = (left .+ right) ./ 2 .- (width / 8) .* (f_right .- f_left)
    f_mid = rhs_block(problem, mesh.midpoint_offsets[interval], midpoint, frequency_reals, scale)
    return right .- left .- (width / 6) .* (f_left .+ 4 .* f_mid .+ f_right)
end

"""
    boundary_residual(problem, layout, side, block, parameter_block, amplitudes, scale)

Mismatch between the node values of `R`, `Q` at an endpoint and the asymptotic
data evaluated there with the current frequencies and leading coefficients:
the Frobenius series at the horizon, and at the outer end the truncated series
or, for propagating components with `outer_boundary.data == :integrated`, the
far series integrated inward (`integrated_outer_values`). `side` is `:horizon`
or `:outer`.
"""
function boundary_residual(problem::RadialProblem, layout::UnknownLayout, side::Symbol, block::AbstractVector,
                           parameter_block::AbstractVector, amplitudes::AbstractVector, scale)
    mode_set, background, mass_coupling = problem.mode_set, problem.background, problem.mass_coupling
    parameters = RadialParameters(layout, parameter_block, amplitudes, scale)
    frequencies = component_frequencies(mode_set, parameters.mode_frequencies)
    state = NodeState(block)
    offset = side == :horizon ? problem.mesh.offsets[1] : problem.mesh.offsets[end]

    value_type = promote_type(eltype(state.profiles), eltype(frequencies), Complex{typeof(mass_coupling)})
    expected = NodeState(similar(state.profiles, value_type), similar(state.slopes, value_type))
    for (component_index, component) in enumerate(mode_set.components)
        channels = mode_set.component_channels[component_index]
        mixing = problem.mixing[component_index]
        frequency = frequencies[component_index]
        if side == :horizon
            series = HorizonSeries(background, mass_coupling, frequency, component.m, component.ls,
                                   mixing, parameters.horizon_leading[channels], problem.horizon_order)
            expected.profiles[channels], expected.slopes[channels] = endpoint_values(series, background, offset)
        elseif uses_integrated_data(problem.outer_boundary, frequency, mass_coupling)
            integrated = integrated_outer_values(background, mass_coupling, frequency, component.m, component.ls,
                                                 mixing, parameters.outer_leading[channels], problem.outer_boundary, offset)
            expected.profiles[channels], expected.slopes[channels] = integrated.profiles, integrated.slopes
        else
            series = OuterSeries(background, mass_coupling, frequency, component.m, component.ls,
                                 mixing, parameters.outer_leading[channels], problem.outer_order)
            expected.profiles[channels], expected.slopes[channels] = endpoint_values(series, background, offset)
        end
    end
    return node_block(NodeState(state.profiles .- expected.profiles, state.slopes .- expected.slopes))
end

"""
    residual(problem, layout, unknowns, amplitudes, scale)

Full residual vector of the discretised boundary-value problem.
"""
function residual(problem::RadialProblem, layout::UnknownLayout, unknowns::AbstractVector,
                  amplitudes::AbstractVector, scale)
    parameter_block = unknowns[parameter_range(layout)]
    frequency_reals = parameter_block[frequency_range(layout)]
    result = zeros(promote_type(eltype(unknowns), eltype(amplitudes), typeof(problem.mass_coupling)), n_unknowns(layout))
    node_rhs = node_rhs_blocks(problem, layout, unknowns, frequency_reals, scale)
    for interval in 1:layout.n_nodes-1
        result[interval_rows(layout, interval)] = interval_residual(
            problem, interval, unknowns[node_range(layout, interval)], unknowns[node_range(layout, interval + 1)],
            node_rhs[interval], node_rhs[interval + 1], frequency_reals, scale)
    end
    result[horizon_rows(layout)] = boundary_residual(problem, layout, :horizon, unknowns[node_range(layout, 1)],
                                                     parameter_block, amplitudes, scale)
    result[outer_rows(layout)] = boundary_residual(problem, layout, :outer, unknowns[node_range(layout, layout.n_nodes)],
                                                   parameter_block, amplitudes, scale)
    return result
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
    rhs_derivatives(problem, s, block, frequency_reals, scale)

Analytic Jacobian of `rhs_block`. The right-hand side is complex-linear in the
node state except for the nonlinear source, whose real Jacobian comes from
`nonlinear_source_jacobian`; the frequency dependence enters only through the
linear coefficients `ΔB_{l,q}` and `a²(ω̃_q² − α²)`.
"""
function rhs_derivatives(problem::RadialProblem, s, block::AbstractVector, frequency_reals::AbstractVector, scale)
    background, mode_set, mass_coupling = problem.background, problem.mode_set, problem.mass_coupling
    state = NodeState(block)
    mode_frequencies = _complexify(frequency_reals)
    frequencies = component_frequencies(mode_set, mode_frequencies)
    a = background.a
    r = radius(background, s)
    horizon_function, sum_r2a2 = delta(background, s), r2a2(background, s)
    n_chan = n_channels(mode_set)
    n_state = 4n_chan
    wrt_state = zeros(n_state, n_state)
    wrt_frequencies = zeros(n_state, length(frequency_reals))

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
        frequency = frequencies[component_index]
        m = component.m
        mixing = problem.mixing[component_index]
        coupled = mixing * state.profiles[channels]
        coupling_strength = a^2 * (frequency^2 - mass_coupling^2)
        for (local_index, channel) in enumerate(channels)
            l = component.ls[local_index]
            delta_B = (sum_r2a2 * frequency - a * m)^2 +
                      horizon_function * (2a * m * frequency - a^2 * frequency^2 - mass_coupling^2 * r^2 - l * (l + 1))
            for (other_local, other_channel) in enumerate(channels)
                factor = horizon_function * coupling_strength * mixing[local_index, other_local]
                other_channel == channel && (factor += delta_B)
                _add_complex_block!(wrt_state, 4channel - 1, 4other_channel - 3, -factor / sum_r2a2^2)
            end
            # ∂(forcing)/∂ω̃_q, then chain through ω̃_q = Σ_i q_i ω̃_i.
            forcing_wrt_frequency = (2sum_r2a2 * (sum_r2a2 * frequency - a * m) + horizon_function * (2a * m - 2a^2 * frequency)) *
                                    state.profiles[channel] + horizon_function * 2a^2 * frequency * coupled[local_index]
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

    source_jacobian = nonlinear_source_jacobian(problem.quadrature, problem.potential, background, mass_coupling, s,
                                                scale .* state.profiles)
    source_factor = -horizon_function / sum_r2a2^2
    for channel in 1:n_chan, other_channel in 1:n_chan
        wrt_real_part = source_factor * source_jacobian.wrt_real[channel, other_channel]
        wrt_imag_part = source_factor * source_jacobian.wrt_imag[channel, other_channel]
        wrt_state[4channel - 1, 4other_channel - 3] += real(wrt_real_part)
        wrt_state[4channel, 4other_channel - 3] += imag(wrt_real_part)
        wrt_state[4channel - 1, 4other_channel - 2] += real(wrt_imag_part)
        wrt_state[4channel, 4other_channel - 2] += imag(wrt_imag_part)
    end
    RhsDerivatives(wrt_state, wrt_frequencies)
end

"""
    IntervalJacobian

Blocks of one interval's MIRK4 residual with respect to the left node, the
right node, and the mode frequencies, obtained by chaining the right-hand-side
derivatives through the Hermite midpoint.
"""
struct IntervalJacobian{T<:Real}
    wrt_left::Matrix{T}
    wrt_right::Matrix{T}
    wrt_frequencies::Matrix{T}
end

function interval_jacobian(problem::RadialProblem, interval::Int, left_derivs::RhsDerivatives, right_derivs::RhsDerivatives,
                           left::AbstractVector, right::AbstractVector, f_left::AbstractVector, f_right::AbstractVector,
                           frequency_reals::AbstractVector, scale)
    mesh = problem.mesh
    width = mesh.rstar[interval + 1] - mesh.rstar[interval]
    midpoint = (left .+ right) ./ 2 .- (width / 8) .* (f_right .- f_left)
    mid_derivs = rhs_derivatives(problem, mesh.midpoint_offsets[interval], midpoint, frequency_reals, scale)

    identity = I(length(left))
    midpoint_wrt_left = identity / 2 .+ (width / 8) .* left_derivs.wrt_state
    midpoint_wrt_right = identity / 2 .- (width / 8) .* right_derivs.wrt_state
    midpoint_wrt_frequencies = -(width / 8) .* (right_derivs.wrt_frequencies .- left_derivs.wrt_frequencies)

    wrt_left = -identity .- (width / 6) .* (left_derivs.wrt_state .+ 4 .* mid_derivs.wrt_state * midpoint_wrt_left)
    wrt_right = identity .- (width / 6) .* (right_derivs.wrt_state .+ 4 .* mid_derivs.wrt_state * midpoint_wrt_right)
    wrt_frequencies = -(width / 6) .* (left_derivs.wrt_frequencies .+ right_derivs.wrt_frequencies .+
                                       4 .* (mid_derivs.wrt_frequencies .+ mid_derivs.wrt_state * midpoint_wrt_frequencies))
    IntervalJacobian(wrt_left, wrt_right, wrt_frequencies)
end

"""
    jacobian(problem, layout, unknowns, amplitudes, scale)

Sparse Jacobian of `residual` with respect to `unknowns`: block-bidiagonal in
the node blocks, bordered by the dense parameter columns.
"""
function jacobian(problem::RadialProblem, layout::UnknownLayout, unknowns::AbstractVector,
                  amplitudes::AbstractVector, scale)
    mesh = problem.mesh
    parameter_block = unknowns[parameter_range(layout)]
    frequency_reals = parameter_block[frequency_range(layout)]
    parameter_columns = parameter_range(layout)
    frequency_columns = parameter_columns[frequency_range(layout)]

    rows, cols, vals = Int[], Int[], Float64[]
    function stamp!(row_indices, col_indices, block)
        for (j, col) in enumerate(col_indices), (i, row) in enumerate(row_indices)
            push!(rows, row); push!(cols, col); push!(vals, block[i, j])
        end
    end

    node_derivs = [rhs_derivatives(problem, mesh.offsets[node], unknowns[node_range(layout, node)], frequency_reals, scale)
                   for node in 1:layout.n_nodes]
    node_rhs = node_rhs_blocks(problem, layout, unknowns, frequency_reals, scale)
    for interval in 1:layout.n_nodes-1
        left, right = node_range(layout, interval), node_range(layout, interval + 1)
        blocks = interval_jacobian(problem, interval, node_derivs[interval], node_derivs[interval + 1],
                                   unknowns[left], unknowns[right], node_rhs[interval], node_rhs[interval + 1],
                                   frequency_reals, scale)
        stamp!(interval_rows(layout, interval), left, blocks.wrt_left)
        stamp!(interval_rows(layout, interval), right, blocks.wrt_right)
        stamp!(interval_rows(layout, interval), frequency_columns, blocks.wrt_frequencies)
    end

    n_state = node_block_size(layout)
    for (side, node, row_indices) in ((:horizon, 1, horizon_rows(layout)), (:outer, layout.n_nodes, outer_rows(layout)))
        block = unknowns[node_range(layout, node)]
        combined = ForwardDiff.jacobian(vcat(block, parameter_block)) do inputs
            boundary_residual(problem, layout, side, inputs[1:n_state], inputs[n_state+1:end], amplitudes, scale)
        end
        stamp!(row_indices, node_range(layout, node), combined[:, 1:n_state])
        stamp!(row_indices, parameter_columns, combined[:, n_state+1:end])
    end
    total = n_unknowns(layout)
    return sparse(rows, cols, vals, total, total)
end

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

"""
    solve_radial(problem, amplitudes; guess=nothing, tol=1e-10, maxiter=30, verbose=false,
                 reuse_factorization=true, chord_threshold=1e-3)

Solve the coupled radial boundary-value problem at the prescribed amplitudes.
`guess` may be a previous `RadialSolution` (warm start, rescaled to the new
amplitudes) or `nothing` for the hydrogenic initial guess. Newton iterations
stop when both the residual and the step fall below `tol` in the max norm.

With `reuse_factorization` and a warm start, the iteration begins as a chord
method: steps use the guess's factorised Jacobian (and later the most recent
one) as long as each step reduces the residual by at least a factor of two,
and the Jacobian is only re-assembled and refactorised when that fails; chord
steps are attempted only while the residual is below `chord_threshold`. Along
a time evolution, where successive amplitudes are close, this replaces most
of the factorisations by residual evaluations. The returned solution always
carries a factorisation taken at (or within `√tol` of) the converged point,
as the amplitude sensitivities require.
"""
function solve_radial(problem::RadialProblem, amplitudes::AbstractVector; guess=nothing,
                      tol=1e-10, maxiter=30, verbose=false, reuse_factorization=true, chord_threshold=1e-3)
    layout = UnknownLayout(problem)
    scale = maximum(amplitudes)
    warm = guess isa RadialSolution && length(guess.unknowns) == n_unknowns(layout)
    unknowns = warm ? warm_start_unknowns(guess, scale) : initial_unknowns(problem, layout, amplitudes, scale)
    current = residual(problem, layout, unknowns, amplitudes, scale)
    current_norm = norm(current, Inf)

    # Chord steps are only worthwhile near a solution; far from it every rejected step costs a residual.
    factorization = warm && reuse_factorization && current_norm < chord_threshold ? guess.jacobian_factorization : nothing
    factorization_distance = Inf   # max-norm distance from the current point to where the factorisation was taken
    converged = false
    iterations = 0
    stalled = 0
    previous_step_norm = Inf
    for iteration in 1:maxiter
        iterations = iteration
        chord = factorization !== nothing
        if !chord
            factorization = lu(jacobian(problem, layout, unknowns, amplitudes, scale))
            factorization_distance = 0.0
        end
        step = -(factorization \ current)
        # Backtracking line search on the residual norm (Armijo condition); a
        # chord step must at least halve the residual or it is abandoned.
        damping = 1.0
        trial, trial_residual, trial_norm = unknowns, current, current_norm
        while damping ≥ 1e-4
            trial = unknowns .+ damping .* step
            trial_residual = residual(problem, layout, trial, amplitudes, scale)
            trial_norm = norm(trial_residual, Inf)
            trial_norm ≤ (1 - 1e-4damping) * current_norm && break
            damping /= 2
        end
        if chord && trial_norm > max(current_norm / 2, tol)
            verbose && println("  chord step $iteration rejected (|F| $current_norm → $trial_norm); refactorising")
            factorization = nothing
            continue
        end
        unknowns, current, current_norm = trial, trial_residual, trial_norm
        reuse_factorization && current_norm < chord_threshold || (factorization = nothing)
        factorization_distance += damping * norm(step, Inf)
        verbose && println("  $(chord ? "chord" : "newton") $iteration: |step| = $(norm(step, Inf))  |F| = $current_norm  damping = $damping")
        if current_norm < tol && damping * norm(step, Inf) < tol
            converged = true
            break
        end
        # Residual at its floor while the step neither vanishes nor shrinks: the system is (nearly)
        # singular, typically because an amplitude is imposed where its profile is negligible. Stop
        # rather than spin. (Chord steps near convergence shrink geometrically, which is fine.)
        step_norm = damping * norm(step, Inf)
        if current_norm < tol && step_norm > previous_step_norm / 2
            stalled += 1
            if stalled ≥ 3
                @warn "solve_radial: |F| = $current_norm is below tol but the Newton step is $(norm(step, Inf)); " *
                      "the system is ill-conditioned (is the outer boundary too far out for the amplitude normalization?)"
                break
            end
        else
            stalled = 0
        end
        previous_step_norm = step_norm
    end
    converged || @warn "solve_radial did not converge: |F| = $current_norm after $iterations iterations"
    # Sensitivities need the Jacobian at the converged point: refactorise unless the
    # current factorisation was taken within √tol of it (Jacobian error O(distance)).
    if factorization === nothing || factorization_distance > sqrt(tol)
        factorization = lu(jacobian(problem, layout, unknowns, amplitudes, scale))
    end
    parameters = RadialParameters(layout, unknowns[parameter_range(layout)], amplitudes, scale)
    RadialSolution(problem, layout, collect(float.(amplitudes)), float(scale), unknowns, parameters, factorization,
                   iterations, current_norm, converged)
end

"""Rescale a previous solution's unknowns to a new reference `scale`; frequencies are untouched."""
function warm_start_unknowns(previous::RadialSolution, scale)
    unknowns = previous.unknowns .* (previous.scale / scale)
    layout = previous.layout
    frequency_indices = parameter_range(layout)[frequency_range(layout)]
    unknowns[frequency_indices] = previous.unknowns[frequency_indices]
    return unknowns
end
