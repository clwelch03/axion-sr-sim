# julia --project=. src/adiabatic/drivers/branch_energy_derivatives.jl \
#     checkpoint.jls output_prefix [shard_index=0] [shard_count=1]
# Recompute energy, signed fluxes, and their amplitude derivatives at saved points.
using Serialization, LinearAlgebra, SHA
include(joinpath(@__DIR__, "../load.jl"))

function branch_energy_derivatives(
    checkpoint_path, output_prefix; shard_index=0, shard_count=1,
)
    0 <= shard_index < shard_count || throw(ArgumentError("invalid shard selection"))

    # Reuse the seed's physics and spatial discretization.
    saved = deserialize(checkpoint_path)
    input_sha256 = bytes2hex(sha256(read(checkpoint_path)))
    spin, alpha, modes = saved.spin, saved.alpha, saved.modes
    potential = Symbol(saved.potential)
    fundamental_m = modes[1][2]
    options = merge(saved.options, (; potential, fundamental_m))

    maximum_l = maximum(first, modes)
    grid = angular_grid(saved.n_polar, saved.n_azimuth)
    polar = build_polar_values(maximum_l, grid.cosine_theta_nodes)
    azimuthal = build_azimuthal_values(maximum_l, grid.azimuth_nodes)
    couplings = build_cosine_squared_couplings(maximum_l, polar, grid)
    tolerance = 1e-10

    metadata = (;
        input_sha256, checkpoint_path=abspath(checkpoint_path),
        spin, alpha, modes, options,
        n_polar=saved.n_polar, n_azimuth=saved.n_azimuth,
        tolerance, shard_index, shard_count, flux_derivatives=true,
    )

    # Both callbacks use the same normalized matching state and settings.
    residual = (packed, amplitude) -> pack_complex(
        scaled_nonlinear_matching_residual(
            unpack_complex(packed), amplitude, spin, alpha, modes, couplings,
            polar, azimuthal, grid; options...,
        ),
    )
    observables = (packed, amplitude) -> energy_and_fluxes(
        packed, amplitude, spin, alpha, modes, couplings,
        polar, azimuthal, grid; options...,
    )

    # Resume only compatible outputs; older tables lack the flux slopes.
    mkpath(dirname(output_prefix))
    rows, solutions = NamedTuple[], NamedTuple[]
    if isfile(output_prefix * ".jls")
        previous = deserialize(output_prefix * ".jls")
        previous.metadata == metadata ||
            error("checkpoint metadata differs; use a new output prefix")
        append!(rows, previous.rows)
        append!(solutions, previous.solutions)
    end

    completed = Set(row.branch_index for row in rows)
    indices = collect(shard_index+1:shard_count:length(saved.continuation))
    println("SETTINGS ", metadata, " point_count=", length(indices))
    flush(stdout)

    # Save after every accepted point, replacing each file via a temporary file.
    function checkpoint()
        serialize(output_prefix * ".jls.tmp", (; metadata, rows, solutions))
        mv(output_prefix * ".jls.tmp", output_prefix * ".jls"; force=true)

        open(output_prefix * ".csv.tmp", "w") do io
            println(io, join(string.(keys(first(rows))), ','))
            for row in rows
                println(io, join(values(row), ','))
            end
        end
        mv(output_prefix * ".csv.tmp", output_prefix * ".csv"; force=true)
    end

    for branch_index in indices
        branch_index in completed && continue

        entry = saved.continuation[branch_index]
        amplitude, guess = entry.amplitude, entry.solution
        started = time()
        println("START index=", branch_index, " amplitude=", amplitude)
        flush(stdout)

        # Correct the seed before differentiating the stationary branch.
        solution = solve_nonlinear_mode(
            guess.frequency, amplitude, spin, alpha, modes, couplings,
            polar, azimuthal, grid;
            horizon_guess=guess.horizon_amplitudes,
            outer_tail_guess=guess.outer_amplitudes[2:end],
            tolerance, max_iterations=12, jacobian_backend=:forwarddiff,
            options...,
        )
        solution.scaled_residual_norm <= tolerance ||
            error("unconverged stationary solution")

        packed = pack_complex(vcat(
            solution.frequency,
            solution.horizon_amplitudes ./ amplitude,
            solution.outer_amplitudes[2:end] ./ amplitude,
        ))

        # One state tangent supplies all three observable slopes.
        derivative_started = time()
        derivative = implicit_observable_derivatives(residual, observables, packed, amplitude)
        derivative_seconds = time() - derivative_started

        all(isfinite, derivative.state_derivative) &&
            all(isfinite, derivative.observable_derivatives) ||
            error("nonfinite implicit derivative")

        cloud_energy, horizon_flux, outer_flux = observables(packed, amplitude)
        energy_derivative, horizon_flux_derivative, outer_flux_derivative =
            derivative.observable_derivatives
        total_flux = horizon_flux + outer_flux
        total_flux_derivative = horizon_flux_derivative + outer_flux_derivative

        amplitude_rate = iszero(energy_derivative) ? NaN : -total_flux / energy_derivative
        frequency_derivative = first(unpack_complex(derivative.state_derivative))

        row = (;
            branch_index, amplitude, energy=cloud_energy, energy_derivative,
            horizon_flux, outer_flux, total_flux, amplitude_rate,
            horizon_flux_derivative, outer_flux_derivative, total_flux_derivative,
            fractional_growth=amplitude_rate/amplitude,

            omega_real=real(solution.frequency), omega_imag=imag(solution.frequency),
            omega_real_derivative=real(frequency_derivative),
            omega_imag_derivative=imag(frequency_derivative),

            residual=solution.scaled_residual_norm, iterations=solution.iterations,
            derivative_seconds, elapsed_seconds=time()-started,
        )
        push!(rows, row)
        push!(solutions, (;
            branch_index, amplitude, solution, state_derivative=derivative.state_derivative,
        ))
        checkpoint()

        println("RESULT ", row)
        flush(stdout)
    end
    println("COMPLETED ", length(rows), " / ", length(indices), " points")
end


if abspath(PROGRAM_FILE) == @__FILE__
    2 <= length(ARGS) <= 4 ||
        error("require checkpoint.jls output_prefix [shard_index] [shard_count]")

    branch_energy_derivatives(
        ARGS[1], ARGS[2];
        shard_index=length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 0,
        shard_count=length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 1,
    )
end
