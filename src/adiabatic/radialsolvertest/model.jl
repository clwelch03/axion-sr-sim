# Experimental, type-generic physics adapters; production files are unchanged.
using LinearAlgebra, SparseArrays, Serialization, ForwardDiff
for filename in ("angular", "nonlinear", "tortoise", "radial", "boundary_series",
                 "boundaries", "integration", "matching", "continuation_predictor", "energy")
    include(joinpath(@__DIR__, "..", filename * ".jl"))
end

"""
    benchmark_context(saved)

Build reusable angular tables for the frozen six-mode cosine problem.

Preserve the production angular quadrature and boundary conventions. Precompute
harmonics for a type-generic RHS that supports ForwardDiff real dual numbers.
"""
function benchmark_context(saved)
    grid = angular_grid(saved.n_polar, saved.n_azimuth)
    polar = build_polar_values(5, grid.cosine_theta_nodes)
    azimuthal = build_azimuthal_values(5, grid.azimuth_nodes)
    couplings = build_cosine_squared_couplings(5, polar, grid)
    harmonics = hcat([vec(polar[mode] * transpose(azimuthal[mode[2]])) for mode in saved.modes]...)
    weights = repeat(grid.polar_weights, length(grid.azimuth_nodes)) .* grid.azimuth_weight
    cosine_squared = repeat(grid.cosine_theta_nodes.^2, length(grid.azimuth_nodes))
    coupling_matrix = [m == other_m ? couplings[(l, other_l, m)] : 0.0
        for (l,m) in saved.modes, (other_l,other_m) in saved.modes]
    return (; spin=saved.spin, alpha=saved.alpha, modes=saved.modes, options=saved.options,
        grid, polar, azimuthal, couplings, harmonics, weights, cosine_squared, coupling_matrix)
end

"""
    generic_boundary(side, frequency, amplitudes, context)

Assemble boundary data without discarding automatic derivatives.

Reuse the production asymptotic series; promote the output element type from
frequency and amplitudes instead of allocating a ComplexF64 array.
"""
function generic_boundary(side, frequency, amplitudes, context)
    count = length(context.modes)
    state = zeros(promote_type(typeof(frequency), eltype(amplitudes)), 2count)
    position = side === :horizon ? context.options.inner : context.options.outer
    for harmonic in unique(last.(context.modes))
        indices = findall(mode -> mode[2] == harmonic, context.modes)
        degrees = first.(context.modes[indices])
        mode_frequency = harmonic * frequency
        data = if side === :horizon
            horizon_boundary_data(position, context.spin, context.alpha, mode_frequency,
                degrees, harmonic, context.couplings, amplitudes[indices];
                order=get(context.options, :horizon_order, 3))
        else
            kind = real(mode_frequency) < context.alpha ? :bound : :outgoing
            outer_boundary_data(position, context.spin, context.alpha, mode_frequency,
                degrees, harmonic, context.couplings, amplitudes[indices];
                kind, order=get(context.options, :outer_order, 7))
        end
        state[indices] = data.values
        state[count .+ indices] = data.tortoise_slopes
    end
    return state
end

"""
    radial_rhs(state, frequency, geometry, context)

Evaluate the real-packed nonlinear radial equations with cached angular tables.

Use the same sine remainder and quadrature as the production cosine RHS.
`state` packs real then imaginary parts of physical `[R; Q]`.
"""
function radial_rhs(state, frequency, geometry, context)
    values = unpack_complex(state)
    count = length(context.modes)
    derivative = similar(values)
    source = zeros(eltype(values), count)
    for point in eachindex(context.weights)
        field = zero(eltype(state))
        for mode in 1:count
            field += 2real(values[mode] * context.harmonics[point, mode])
        end
        weighted = context.alpha^2 * (geometry.radius^2 +
            context.spin^2 * context.cosine_squared[point]) *
            sine_remainder(field) * context.weights[point]
        for mode in 1:count
            source[mode] += weighted * conj(context.harmonics[point, mode])
        end
    end
    for (index, (degree, harmonic)) in enumerate(context.modes)
        coefficients = tortoise_operator_coefficients(geometry, context.spin,
            context.alpha, harmonic * frequency, degree, harmonic)
        coupled = sum(context.coupling_matrix[index, other] * values[other] for other in 1:count)
        derivative[index] = values[count+index]
        derivative[count+index] = -coefficients.slope_coefficient * values[count+index] -
            coefficients.local_coefficient * values[index] - coefficients.coupling_scale * coupled -
            coefficients.source_scale * source[index]
    end
    return pack_complex(derivative)
end

"""
    integration_mesh(start, stop, step, context; count_override=nothing)

Cache geometry at the endpoints and midpoint of every fixed RK4 step.

Use `count_override` to preserve the exact step count of an aligned subinterval.
"""
function integration_mesh(start, stop, step, context; count_override=nothing)
    count = isnothing(count_override) ? max(1, ceil(Int, abs(stop-start)/step)) : count_override
    spacing = (stop-start)/count
    geometry = [tortoise_geometry(start + index*spacing/2, context.spin) for index in 0:2count]
    return (; spacing, geometry, count)
end

"""
    integrate_cached(initial, frequency, mesh, context; save=false)

Advance the real-packed state with RK4 and cached geometry.
"""
function integrate_cached(initial, frequency, mesh, context; save=false)
    state = copy(initial)
    profiles = save ? [copy(state)] : nothing
    for index in 1:mesh.count
        spacing = mesh.spacing
        left, middle, right = mesh.geometry[2index-1:2index+1]
        first = radial_rhs(state, frequency, left, context)
        second = radial_rhs(state + spacing/2*first, frequency, middle, context)
        third = radial_rhs(state + spacing/2*second, frequency, middle, context)
        fourth = radial_rhs(state + spacing*third, frequency, right, context)
        state += spacing/6 * (first + 2second + 2third + fourth)
        save && push!(profiles, copy(state))
    end
    return save ? hcat(profiles...) : state
end

"""
    shooting_problem(context, amplitude, prediction)

Construct the same normalized two-sided matching system as production.
"""
function shooting_problem(context, amplitude, prediction)
    count = length(context.modes)
    initial = pack_complex(vcat(prediction.frequency, prediction.horizon_amplitudes/amplitude,
        prediction.outer_amplitudes[2:end]/amplitude))
    inner_mesh = integration_mesh(context.options.inner, context.options.match, context.options.step, context)
    outer_mesh = integration_mesh(context.options.outer, context.options.match, context.options.step, context)
    boundaries = function(parameters)
        unknowns = unpack_complex(parameters)
        frequency = unknowns[1]
        horizon = generic_boundary(:horizon, frequency, amplitude*unknowns[2:count+1], context)
        outer = generic_boundary(:outer, frequency,
            amplitude*vcat(one(frequency), unknowns[count+2:end]), context)
        return frequency, pack_complex(horizon), pack_complex(outer)
    end
    residual = function(parameters)
        frequency, horizon, outer = boundaries(parameters)
        return (integrate_cached(horizon, frequency, inner_mesh, context) -
            integrate_cached(outer, frequency, outer_mesh, context))/amplitude
    end
    return (; initial, residual, boundaries, inner_mesh, outer_mesh, amplitude, context)
end

"""
    physical_solution(parameters, amplitude, iterations, residual_norm)

Recover the physical frequency and amplitudes from the real Newton vector.
"""
function physical_solution(parameters, amplitude, iterations, residual_norm)
    unknowns = unpack_complex(parameters)
    count = length(unknowns) ÷ 2
    return (; frequency=unknowns[1], horizon_amplitudes=amplitude*unknowns[2:count+1],
        outer_amplitudes=amplitude*vcat(one(unknowns[1]), unknowns[count+2:end]),
        iterations, scaled_residual_norm=residual_norm)
end
