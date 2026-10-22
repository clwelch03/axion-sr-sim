using Serialization, Printf, LinearAlgebra, Test
include(joinpath(@__DIR__, "../../../../src/adiabatic/multiphase/load.jl"))
include(joinpath(@__DIR__, "../../../../src/adiabatic/single_cloud/load.jl"))

# Run from a frozen checkout; the output prefix may point outside that checkout.
# Every accepted endpoint is checkpointed so a stopped job can resume it.
"""
    save_checkpoint(path, value)

TODO: describe atomic replacement of a campaign checkpoint.
"""
function save_checkpoint(path, value)
    open(path * ".tmp", "w") do io
        serialize(io, value)
    end
    mv(path * ".tmp", path; force=true)
end

"""
    print_sample(label, sample)

TODO: describe the progress diagnostics and their units.
"""
function print_sample(label, sample)
    @printf("%s: A=(%.8e, %.8e), E=%.8e, J=%.8e, residual=%.3e, cond=%.3e\n",
        label, sample.solution.amplitudes..., sample.charges...,
        sample.solution.residual_norm, sample.balance.scaled_condition)
    println("  logarithmic amplitude rates: ", sample.balance.log_rates)
    println("  total losses: ", sample.losses)
    flush(stdout)
end

"""
    validate_neighbors(problem, initial, prefix)

TODO: describe independently corrected neighbors and derivative comparisons.
"""
function validate_neighbors(problem::TwoPhaseProblem, initial, prefix)
    path = prefix * "-neighbors.csv"
    marker = prefix * "-neighbors.jls"
    isfile(marker) && return deserialize(marker)
    rows = NamedTuple[]
    reference = initial.solution
    derivative = initial.response.observable_derivatives
    for fraction in (0.01, 0.005), column in 1:2
        displacement = fraction * reference.amplitudes[column]
        offset = zeros(2)
        offset[column] = displacement
        neighbors = map((-1, 1)) do sign
            amplitudes = reference.amplitudes + sign .* offset
            guess = reference.state + initial.response.state_derivative * (sign .* offset)
            solution = correct_two_phase_profile(problem, amplitudes, guess; tolerance=1e-9)
            observables = two_phase_trial_observables(solution.state, amplitudes, problem)
            (; solution, charges=[observables.energy, observables.angular_momentum])
        end
        finite_charges = (neighbors[2].charges-neighbors[1].charges) / (2displacement)
        finite_state = (neighbors[2].solution.state-neighbors[1].solution.state) / (2displacement)
        charge_error = norm(finite_charges-derivative[:,column]) / norm(derivative[:,column])
        state_error = norm(finite_state-initial.response.state_derivative[:,column]) /
            norm(initial.response.state_derivative[:,column])
        push!(rows, (; fraction, column, charge_error, state_error,
            energy_derivative=derivative[1,column], energy_difference=finite_charges[1],
            angular_derivative=derivative[2,column], angular_difference=finite_charges[2],
            minus_residual=neighbors[1].solution.residual_norm,
            plus_residual=neighbors[2].solution.residual_norm))
        println("Neighbor check: ", rows[end])
        flush(stdout)
        @test charge_error < 3e-4
        @test state_error < 3e-4
    end
    open(path, "w") do io
        println(io, join(string.(keys(first(rows))), ','))
        for row in rows
            println(io, join(values(row), ','))
        end
    end
    save_checkpoint(marker, rows)
    return rows
end

"""
    write_trajectory(path, history)

TODO: document trajectory columns and atomic table replacement.
"""
function write_trajectory(path, history)
    open(path * ".tmp", "w") do io
        println(io, "time,A1,A2,energy,angular_momentum,omega1_real,omega2_real," *
            "omega1_imag,omega2_imag,horizon_energy_loss,outer_energy_loss," *
            "horizon_angular_loss,outer_angular_loss,log_rate1,log_rate2," *
            "scaled_condition,matching_condition,matching_residual,response_residual," *
            "integrated_energy_loss,integrated_angular_loss,energy_balance_error,angular_balance_error")
        for index in eachindex(history.times)
            sample = history.samples[index]
            frequencies = unpack_complex(sample.solution.state)[1:2]
            observables = sample.response.observables
            row = [history.times[index]; sample.solution.amplitudes; sample.charges;
                real.(frequencies); imag.(frequencies);
                observables.horizon_energy_flux; observables.outer_energy_flux;
                observables.horizon_angular_momentum_flux; observables.outer_angular_momentum_flux;
                sample.balance.log_rates; sample.balance.scaled_condition;
                sample.response.state_jacobian_condition; sample.solution.residual_norm;
                sample.response.response_residual; history.integrated_losses[index];
                history.balance_errors[index]]
            println(io, join(row, ','))
        end
    end
    mv(path * ".tmp", path; force=true)
end

"""
    run_trajectory(problem, initial, prefix, final_time, step_count)

TODO: describe accepted-endpoint checkpoints and continuation after interruption.
"""
function run_trajectory(problem::TwoPhaseProblem, initial, prefix, final_time, step_count)
    path = prefix * "-steps$(step_count)"
    if isfile(path * ".jls")
        history = deserialize(path * ".jls")
        history.times[end] == final_time && return history
    else
        history = (; times=Float64[], samples=typeof(initial)[],
            integrated_losses=Vector{Float64}[], balance_errors=Vector{Float64}[])
    end
    offset_time = isempty(history.times) ? 0.0 : history.times[end]
    offset_losses = isempty(history.times) ? zeros(2) : history.integrated_losses[end]
    starting_sample = isempty(history.times) ? initial : history.samples[end]
    callback = (time, sample, losses, errors) -> begin
        absolute_time = offset_time + time
        if !isempty(history.times) && absolute_time == history.times[end]
            return
        end
        cumulative = offset_losses + losses
        push!(history.times, absolute_time)
        push!(history.samples, sample)
        push!(history.integrated_losses, cumulative)
        push!(history.balance_errors, sample.charges - initial.charges + cumulative)
        save_checkpoint(path * ".jls", history)
        write_trajectory(path * ".csv", history)
        print_sample("steps=$step_count, t=$absolute_time", sample)
    end
    evaluate = (amplitudes, guess) -> two_phase_evolution_sample(problem,
        amplitudes, guess; tolerance=1e-9)
    evolve_two_phase(evaluate, starting_sample; final_time=final_time-offset_time,
        time_step=final_time/step_count, callback)
    return history
end

"""
    main(prefix)

TODO: describe the fixed-background cosine pilot and its validation sequence.
"""
function main(prefix)
    mkpath(dirname(prefix))
    # A directory lock prevents accidental simultaneous writers to one prefix.
    mkdir(prefix * ".lock")
    try
        # A prefix belongs to one source snapshot and one spatial configuration.
        # Reject a changed resume instead of silently reusing incompatible results.
        manifest_path = joinpath(@__DIR__, "../../../../provenance.json")
        isfile(manifest_path) || error("run from a source snapshot made by freeze_source.py")
        configuration = (; source_manifest=read(manifest_path, String),
            radial_step=parse(Float64, get(ENV, "AXION_RADIAL_STEP", "0.4")))
        configuration_path = prefix * "-configuration.jls"
        if isfile(configuration_path)
            deserialize(configuration_path) == configuration ||
                error("resume configuration changed; use a distinct output prefix")
        else
            isfile(prefix * "-initial.jls") && error("existing prefix has no configuration record")
            save_checkpoint(configuration_path, configuration)
        end
        initial_path = prefix * "-initial.jls"
        if isfile(initial_path)
            problem, initial = deserialize(initial_path)
        else
            grid = angular_grid(8, 9)
            polar = build_polar_values(3, grid.cosine_theta_nodes)
            couplings = build_cosine_squared_couplings(3, polar, grid)
            cloud_m, spin, alpha = [1, 2], 0.99, 0.42
            components = [fourier_component(q, [l], cloud_m) for (l,q) in
                [(1,(1,0)), (2,(0,1)), (0,(2,-1)), (3,(-1,2))]]
            options = (; inner=-30.0, outer=80.0, match=5.0,
                step=parse(Float64, get(ENV, "AXION_RADIAL_STEP", "0.4")))
            amplitudes = [1e-4, 1e-5]
            first_mode = solve_linear_mode(0.4088+1e-7im, spin, alpha, [1], 1, couplings;
                options..., tolerance=1e-9)
            second_mode = solve_linear_mode(0.4159+1e-9im, spin, alpha, [2], 2, couplings;
                options..., tolerance=1e-9)
            frequencies = [first_mode.frequency, second_mode.frequency]
            horizon = ComplexF64[amplitudes[1]*first_mode.horizon_amplitudes[1],
                amplitudes[2]*second_mode.horizon_amplitudes[1], 0, 0]
            outer = ComplexF64[amplitudes; 0; 0]
            # Set up the fixed response settings alongside the reference profile.
            length(frequencies) == length(amplitudes) == length(cloud_m) == 2 ||
                throw(DimensionMismatch("require two independent clouds"))
            integration = (; options..., horizon_order=3, outer_order=7)
            all(isfinite, (integration.inner, integration.match, integration.outer, integration.step)) &&
                integration.inner < integration.match < integration.outer && integration.step > 0 ||
                throw(ArgumentError("require ordered endpoints and a positive finite radial step"))

            # Fix the first outer coefficient of each fundamental component.
            # All others are free, including additional fundamental angular harmonics.
            (; free_outer_coefficient_indices) = two_phase_coefficient_indices(components)

            coefficient_scale = minimum(amplitudes)
            check_cloud_amplitudes_and_coefficient_scale_validity(amplitudes, coefficient_scale)

            # Record the decaying/outgoing branch at the outer boundary for each component.
            ω1, ω2 = real.(frequencies)
            component_frequencies = [component.q[1] * ω1 + component.q[2] * ω2
                for component in components]
            initial_bound = component_frequencies .< alpha
            problem = TwoPhaseProblem(;
                # Background and retained modes
                spin,
                alpha,
                components,
                cloud_m,
                # Angular projection and phase sampling
                couplings,
                polar_values = polar,
                grid,
                phase1_nodes = grid.azimuth_nodes,
                phase2_nodes = grid.azimuth_nodes,
                # Potential and fixed coefficient/residual scale
                potential = :cosine,
                coefficient_scale,
                free_outer_coefficient_indices,
                # Boundary branches and radial integration settings
                initial_bound,
                integration,
            )
            profile = TwoPhaseProfile(frequencies, horizon, outer)
            guess = pack_two_phase_profile(profile, problem)
            two_phase_profile_admissible(guess, problem) ||
                throw(ArgumentError("reference frequencies are outside the supported boundary branches"))
            initial = two_phase_evolution_sample(problem, amplitudes, guess; tolerance=1e-9)
            save_checkpoint(initial_path, (problem, initial))
        end
        print_sample("Initial matched cosine profile", initial)
        @test initial.response.response_residual < 1e-12
        validate_neighbors(problem, initial, prefix)

        # Changing the numerical coefficient scale must leave physical charges
        # and their derivatives unchanged. It may change the raw Newton condition.
        scale_path = prefix * "-scale-check.jls"
        if !isfile(scale_path)
            scaled_problem = TwoPhaseProblem(;
                spin = problem.spin,
                alpha = problem.alpha,
                components = problem.components,
                cloud_m = problem.cloud_m,
                couplings = problem.couplings,
                polar_values = problem.polar_values,
                grid = problem.grid,
                phase1_nodes = problem.phase1_nodes,
                phase2_nodes = problem.phase2_nodes,
                potential = problem.potential,
                coefficient_scale = 7problem.coefficient_scale,
                free_outer_coefficient_indices = problem.free_outer_coefficient_indices,
                initial_bound = problem.initial_bound,
                integration = problem.integration,
            )
            trial = unpack_complex(initial.solution.state)
            trial[3:end] ./= 7
            scaled_solution = correct_two_phase_profile(scaled_problem,
                initial.solution.amplitudes, pack_complex(trial); tolerance=1e-9)
            scaled_response = two_phase_profile_response(scaled_problem, scaled_solution; tolerance=1e-9)
            scale_error = norm(scaled_response.observable_derivatives-initial.response.observable_derivatives) /
                norm(initial.response.observable_derivatives)
            @test scale_error < 1e-7
            println("Fixed-scale invariance relative error: ", scale_error)
            save_checkpoint(scale_path, (; scale_error))
        end

        # Aim for a 3% log-amplitude change in the faster cloud. This is a
        # numerical demonstration on this finite domain, not a long physical history.
        final_time = 0.03 / maximum(abs.(initial.balance.log_rates))
        coarse = run_trajectory(problem, initial, prefix, final_time, 4)
        fine = run_trajectory(problem, initial, prefix, final_time, 8)
        endpoint_difference = (fine.samples[end].solution.amplitudes -
            coarse.samples[end].solution.amplitudes) ./ fine.samples[end].solution.amplitudes
        coarse_budget = coarse.balance_errors[end] ./ initial.charges
        fine_budget = fine.balance_errors[end] ./ initial.charges
        println("Endpoint relative amplitude difference, 4 vs 8 steps: ", endpoint_difference)
        println("Relative E/J balance errors, 4 steps: ", coarse_budget)
        println("Relative E/J balance errors, 8 steps: ", fine_budget)
        @test maximum(abs.(endpoint_difference)) < 1e-5
        @test maximum(abs.(fine_budget)) < 1e-5
        @test norm(fine_budget) < norm(coarse_budget)
        save_checkpoint(prefix * "-summary.jls", (; final_time, endpoint_difference,
            coarse_budget, fine_budget, radial_step=problem.integration.step))
        println("COMPLETED: matched cosine profiles, two-amplitude response, neighbor checks, and evolution.")
    finally
        rm(prefix * ".lock")
    end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    length(ARGS) == 1 || error("usage: response_and_evolution.jl OUTPUT_PREFIX")
    main(abspath(only(ARGS)))
end
