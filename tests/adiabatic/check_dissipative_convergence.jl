# julia --project=. tests/adiabatic/check_dissipative_convergence.jl
# Convergence of the m = 0 horizon flux per unit secondary energy, F_H / Ê2
# (piece 4e, docs/multimode-piece4-guide.md), on the Mμ = 0.42 cosine branch.
#   1. Radial step at A1 = 2.05 (primary index 50): re-solve the primary at
#      steps 0.1, 0.05, 0.025 (Newton from the saved solution), the secondary
#      and the m = 0 field at twice that. If RK4 error dominates, successive
#      changes shrink by about 16. Same procedure as
#      check_secondary_step_convergence.jl.
#   2. l-truncation: add (6, 0) to the m = 0 field at the saved solutions.
using Printf, Serialization, LinearAlgebra
include(joinpath(@__DIR__, "../../src/adiabatic/drivers/dissipative_fluxes.jl"))

results = joinpath(@__DIR__, "../../runs/adiabatic/multimode_secondary/results")
primary_path = joinpath(@__DIR__, "../../runs/adiabatic/fig2_20260922/results/cosine_branch.jls")

setup = load_secondary_setup(primary_path)
(; primary, secondary_modes, grid, polar, azimuthal, couplings) = setup
secondary_branch = deserialize(joinpath(results, "cosine_secondary.jls"))
(; inner, outer, match) = primary.options
potential = Symbol(primary.potential)

# F_H / Ê2 for a primary solution's profiles and a secondary solved on them.
function flux_per_energy(profiles, secondary, omega1; dissipative_modes=DISSIPATIVE_MODES)
    frequency = 2real(omega1) - real(secondary.frequency)
    solution = solve_dissipative_mode(
        frequency, secondary, profiles, primary.modes, primary.spin, primary.alpha,
        secondary_modes, dissipative_modes, couplings, polar, azimuthal, grid,
    )
    (; energy) = secondary_energy_and_fluxes(
        secondary, profiles, primary.modes, primary.spin, primary.alpha,
        secondary_modes, couplings, polar, azimuthal, grid,
    )
    return solution.horizon_flux / energy
end

println("1. Radial step at A1 = 2.05")
let index = 50
    entry = primary.continuation[index]
    saved_primary = entry.solution
    saved_secondary = secondary_branch.rows[index].solution
    secondary_guess = vcat(saved_secondary.frequency, saved_secondary.horizon_amplitudes,
                           saved_secondary.outer_amplitudes[2:end])
    @printf("%-14s %-22s %s\n", "primary step", "F_H / E2hat", "change")
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
        push!(values, flux_per_energy(profiles, secondary, primary_solution.frequency))
        change = length(values) > 1 ? @sprintf("%+.3e", values[end] - values[end-1]) : ""
        @printf("%-14.3g %-22.15e %s\n", primary_step, values[end], change)
        flush(stdout)
    end
    @printf("ratio of changes: %.2f (16 if RK4 error dominates); relative change at step 0.1: %.2e\n",
            (values[2] - values[1]) / (values[3] - values[2]),
            abs(values[2] - values[1]) / abs(values[2]))
end

println("\n2. l-truncation: (0,0), (2,0), (4,0) against adding (6,0)")
for index in (1, 50, 232)
    (; profiles) = primary_profiles_at(setup, index)
    secondary = secondary_branch.rows[index].solution
    omega1 = primary.continuation[index].solution.frequency
    base = flux_per_energy(profiles, secondary, omega1)
    extended = flux_per_energy(profiles, secondary, omega1;
                               dissipative_modes=[(0, 0), (2, 0), (4, 0), (6, 0)])
    @printf("A1=%-8.4g F_H/E2hat = %.12e  with l=6: %.12e  relative change %.2e\n",
            primary.continuation[index].amplitude, base, extended, abs(extended - base) / base)
    flush(stdout)
end
