# julia --project=. tests/adiabatic/check_secondary_step_convergence.jl
# Radial-step convergence of the averaged secondary at selected primary points.
# Re-solve the primary at steps 0.1, 0.05, 0.025 (Newton from the saved solution),
# then solve the secondary on each at twice that step. If RK4 error dominates,
# successive differences shrink by about 2^4 = 16.
using Serialization, Printf, LinearAlgebra
include(joinpath(@__DIR__, "../../src/adiabatic/drivers/continue_secondary.jl"))

primary_path = joinpath(@__DIR__, "../../runs/adiabatic/fig2_20260922/results/cosine_branch.jls")
secondary_path = joinpath(@__DIR__, "../../runs/adiabatic/multimode_secondary/results/cosine_secondary.jls")

setup = load_secondary_setup(primary_path)
(; primary, secondary_modes, grid, polar, azimuthal, couplings) = setup
secondary_branch = deserialize(secondary_path)
count = length(secondary_modes)
(; inner, outer, match) = primary.options
potential = Symbol(primary.potential)

for index in (50, 232)
    entry = primary.continuation[index]
    amplitude = entry.amplitude
    saved_primary = entry.solution
    saved_secondary = secondary_branch.rows[index].solution
    secondary_guess = vcat(saved_secondary.frequency, saved_secondary.horizon_amplitudes,
                           saved_secondary.outer_amplitudes[2:end])
    @printf("\nA1 = %.6g (primary index %d)\n", amplitude, index)
    @printf("%-8s %-36s %-36s %s\n", "step", "omega1", "omega2", "iterations (1, 2)")

    previous = nothing
    differences = ComplexF64[]
    for primary_step in (0.1, 0.05, 0.025)
        primary_solution = solve_nonlinear_mode(
            saved_primary.frequency, amplitude, primary.spin, primary.alpha,
            primary.modes, couplings, polar, azimuthal, grid;
            horizon_guess=saved_primary.horizon_amplitudes,
            outer_tail_guess=saved_primary.outer_amplitudes[2:end],
            jacobian_backend=:forwarddiff, potential,
            step=primary_step, inner, outer, match,
        )
        profiles = nonlinear_matching_residual(
            primary_solution.frequency, primary_solution.horizon_amplitudes,
            primary_solution.outer_amplitudes, primary.spin, primary.alpha,
            primary.modes, couplings, polar, azimuthal, grid;
            potential, save_profile=true, step=primary_step, inner, outer, match,
        )
        secondary = solve_secondary_point(
            setup, profiles, secondary_guess;
            tolerance=1e-8, jacobian_backend=:forwarddiff, reading=:averaged,
        )
        omega1, omega2 = primary_solution.frequency, secondary.frequency
        @printf("%-8.3g %.12f %+.6ei   %.12f %+.6ei   (%d, %d)\n",
                primary_step, real(omega1), imag(omega1), real(omega2), imag(omega2),
                primary_solution.iterations, secondary.iterations)
        flush(stdout)
        previous === nothing || push!(differences, omega2 - previous)
        previous = omega2
    end

    first_change, second_change = differences
    @printf("omega2 change 0.1->0.05: %+.3e %+.3ei   0.05->0.025: %+.3e %+.3ei\n",
            real(first_change), imag(first_change), real(second_change), imag(second_change))
    @printf("ratio of changes (Re, Im): %.2f, %.2f   (16 if RK4 error dominates)\n",
            real(first_change) / real(second_change), imag(first_change) / imag(second_change))
end
