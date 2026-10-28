# julia --project=. tests/adiabatic/primary_secondary/check_radiative_convergence.jl
# Convergence of the m = 3 flux to infinity per unit secondary energy squared,
# F_inf / Ê2² (piece 5e), on the Mμ = 0.42 cosine branch.
#   1. Radial step at A1 = 2.05 (primary index 50): re-solve the primary at
#      steps 0.1, 0.05, 0.025 (Newton from the saved solution), the secondary
#      and the m = 3 field at twice that. If RK4 error dominates, successive
#      changes shrink by about 16. Same procedure as
#      check_dissipative_convergence.jl.
#   2. l-truncation: add (9, 3) to the m = 3 field at the saved solutions.
#   3. Far outer leg: far start 2500 and 10000 instead of 5000, and far step
#      0.125 instead of 0.25.
using Printf, Serialization, LinearAlgebra
include(joinpath(@__DIR__, "../../../src/adiabatic/primary_secondary/load.jl"))

results = joinpath(@__DIR__, "../../../runs/adiabatic/multimode_secondary/results")
primary_path = joinpath(@__DIR__, "../../../runs/adiabatic/fig2_20260922/results/cosine_branch.jls")

setup = load_primary_secondary_setup(primary_path)
(; primary, secondary_modes, grid) = setup
# Tables up to l = 9 for the truncation test.
polar = build_polar_values(9, grid.cosine_theta_nodes)
azimuthal = build_azimuthal_values(9, grid.azimuth_nodes)
couplings = build_cosine_squared_couplings(9, polar, grid)
setup = merge(setup, (; polar, azimuthal, couplings))
secondary_branch = deserialize(joinpath(results, "cosine_secondary.jls"))
(; inner, outer, match) = primary.options
potential = Symbol(primary.potential)
radiative_modes = [(3, 3), (5, 3), (7, 3)]

# F_inf / Ê2² for a primary solution's profiles and a secondary solved on them.
function flux_per_energy_squared(profiles, secondary, omega1;
                                 modes=radiative_modes, series_options...)
    frequency = 2real(secondary.frequency) - real(omega1)
    solution = solve_radiative_mode(
        frequency, secondary, profiles, primary.modes, primary.spin, primary.alpha,
        secondary_modes, modes, couplings, polar, azimuthal, grid; series_options...,
    )
    (; energy) = secondary_energy_and_fluxes(
        secondary, profiles, primary.modes, primary.spin, primary.alpha,
        secondary_modes, couplings, polar, azimuthal, grid,
    )
    return solution.infinity_flux / energy^2
end

println("1. Radial step at A1 = 2.05")
let index = 50
    entry = primary.continuation[index]
    saved_primary = entry.solution
    saved_secondary = secondary_branch.rows[index].solution
    secondary_guess = vcat(saved_secondary.frequency, saved_secondary.horizon_amplitudes,
                           saved_secondary.outer_amplitudes[2:end])
    @printf("%-14s %-22s %s\n", "primary step", "F_inf / E2hat^2", "change")
    values = Float64[]
    for primary_step in (0.1, 0.05, 0.025)
        primary_solution = solve_nonlinear_mode(
            saved_primary.frequency, entry.amplitude, primary.spin, primary.alpha,
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
        push!(values, flux_per_energy_squared(profiles, secondary, primary_solution.frequency))
        change = length(values) > 1 ? @sprintf("%+.3e", values[end] - values[end-1]) : ""
        @printf("%-14.3g %-22.15e %s\n", primary_step, values[end], change)
        flush(stdout)
    end
    @printf("ratio of changes: %.2f (16 if RK4 error dominates); relative change at step 0.1: %.2e\n",
            (values[2] - values[1]) / (values[3] - values[2]),
            abs(values[2] - values[1]) / abs(values[2]))
end

println("\n2. l-truncation and 3. far outer leg, relative change of F_inf / E2hat^2")
for index in (1, 50, 232)
    (; profiles) = primary_profiles_at(setup, index)
    secondary = secondary_branch.rows[index].solution
    omega1 = primary.continuation[index].solution.frequency
    base = flux_per_energy_squared(profiles, secondary, omega1)
    @printf("A1=%-8.4g F_inf/E2hat^2 = %.12e\n", primary.continuation[index].amplitude, base)
    variants = (("with l = 9", (; modes=[radiative_modes; (9, 3)])),
                ("far start 2500", (; far_position=2500.0)),
                ("far start 10000", (; far_position=10000.0)),
                ("far step 0.125", (; far_step=0.125)))
    for (label, options) in variants
        value = flux_per_energy_squared(profiles, secondary, omega1; options...)
        @printf("    %-16s %.12e  relative change %.2e\n", label, value, abs(value - base) / base)
        flush(stdout)
    end
end
