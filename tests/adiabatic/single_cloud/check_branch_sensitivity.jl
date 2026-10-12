# Physical branch check; run on a compute node from the repository root:
# julia --project=. tests/adiabatic/single_cloud/check_branch_sensitivity.jl [checkpoint.jls]
# Reads a saved branch without modifying it. Reports results to stdout.
using Test, Serialization, LinearAlgebra
include(joinpath(@__DIR__, "../../../src/adiabatic/load.jl"))

function check_branch_sensitivity(checkpoint_path)
    saved = deserialize(checkpoint_path)
    target = 2.564697905478641
    index = argmin(abs.([entry.amplitude for entry in saved.continuation] .- target))
    entry = saved.continuation[index]
    amplitude = entry.amplitude
    spin, alpha, modes = saved.spin, saved.alpha, saved.modes
    potential = Symbol(saved.potential)
    options = merge(saved.options, (; potential))
    mode_count = length(modes)
    maximum_l = maximum(first, modes)
    grid = angular_grid(saved.n_polar, saved.n_azimuth)
    polar = build_polar_values(maximum_l, grid.cosine_theta_nodes)
    azimuthal = build_azimuthal_values(maximum_l, grid.azimuth_nodes)
    couplings = build_cosine_squared_couplings(maximum_l, polar, grid)
    tolerance = 1e-10

    function packed_solution(solution, trial_amplitude)
        return pack_complex(vcat(
            solution.frequency,
            solution.horizon_amplitudes ./ trial_amplitude,
            solution.outer_amplitudes[2:end] ./ trial_amplitude,
        ))
    end

    function correct(trial_amplitude, guess)
        correction_state = unpack_complex(guess)
        solution = solve_nonlinear_mode(
            correction_state[1], trial_amplitude, spin, alpha, modes, couplings,
            polar, azimuthal, grid;
            horizon_guess=trial_amplitude .* correction_state[2:mode_count+1],
            outer_tail_guess=trial_amplitude .* correction_state[mode_count+2:end],
            tolerance, max_iterations=12, jacobian_backend=:forwarddiff,
            progress_callback=(iteration, residual_norm) -> begin
                println("CORRECT amplitude=", trial_amplitude,
                        " iteration=", iteration, " residual=", residual_norm)
                flush(stdout)
            end,
            options...,
        )
        @test solution.scaled_residual_norm <= tolerance
        return packed_solution(solution, trial_amplitude)
    end

    # Differentiate real and imaginary parts independently: the residual
    # contains conjugation, so complex-analytic differentiation is invalid.
    residual = (packed, trial_amplitude) -> pack_complex(
        scaled_nonlinear_matching_residual(
            unpack_complex(packed), trial_amplitude, spin, alpha, modes,
            couplings, polar, azimuthal, grid; options...,
        ),
    )

    function energy(packed, trial_amplitude)
        energy_state = unpack_complex(packed)
        profiles = nonlinear_matching_residual(
            energy_state[1], trial_amplitude .* energy_state[2:mode_count+1],
            trial_amplitude .* vcat(one(energy_state[1]), energy_state[mode_count+2:end]),
            spin, alpha, modes, couplings, polar, azimuthal, grid;
            save_profile=true, options...,
        )
        return sum(
            integrate_energy_branch(
                branch, energy_state[1], spin, alpha, modes, polar, azimuthal, grid;
                potential,
            ) for branch in (profiles.from_horizon, profiles.from_outer)
        )
    end

    println("SETTINGS ", (; checkpoint_path, amplitude, spin, alpha, modes,
        options, n_polar=saved.n_polar, n_azimuth=saved.n_azimuth, tolerance))
    flush(stdout)
    state = correct(amplitude, packed_solution(entry.solution, amplitude))
    base_energy = energy(state, amplitude)
    println("BASE ", (; energy=base_energy, residual=norm(residual(state, amplitude))))
    flush(stdout)

    derivative = implicit_energy_derivative(residual, energy, state, amplitude)
    tangent = derivative.state_derivative
    energy_derivative = derivative.energy_derivative
    @test all(isfinite, tangent)
    @test isfinite(energy_derivative)
    frequency_derivative = first(unpack_complex(tangent))
    println("IMPLICIT ", (; energy_derivative, frequency_derivative,
        state_derivative_norm=norm(tangent)))
    flush(stdout)

    # Nearby solutions are corrected from the same baseline state, without
    # using the implicit tangent as a predictor.
    energy_errors = Float64[]
    state_errors = Float64[]
    for fraction in (0.01, 0.005, 0.0025)
        spacing = fraction * amplitude
        above = correct(amplitude + spacing, state)
        below = correct(amplitude - spacing, state)
        difference_tangent = (above .- below) ./ (2spacing)
        difference_energy = (energy(above, amplitude + spacing) -
                             energy(below, amplitude - spacing)) / (2spacing)
        state_error = norm(difference_tangent - tangent) / norm(tangent)
        energy_error = abs(difference_energy - energy_derivative) / abs(energy_derivative)
        push!(state_errors, state_error)
        push!(energy_errors, energy_error)
        println("COMPARISON ", (; spacing, difference_energy,
            energy_relative_error=energy_error, state_relative_error=state_error,
            frequency_derivative=first(unpack_complex(difference_tangent))))
        flush(stdout)
    end
    # 0.01% agreement is a check of differentiation on this fixed grid,
    # not an error bound on the continuum energy derivative.
    @test last(energy_errors) < 1e-4
    @test last(state_errors) < 1e-4
    @test last(energy_errors) < first(energy_errors)
    @test last(state_errors) < first(state_errors)
end

if abspath(PROGRAM_FILE) == @__FILE__
    default_checkpoint = normpath(joinpath(@__DIR__,
        "../../../runs/adiabatic/fig2_20260922/results/cosine_branch.jls"))
    @testset "Physical branch: implicit versus corrected finite differences" begin
        check_branch_sensitivity(isempty(ARGS) ? default_checkpoint : only(ARGS))
    end
end
