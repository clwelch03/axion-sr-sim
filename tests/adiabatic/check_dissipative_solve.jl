# julia --project=. tests/adiabatic/check_dissipative_solve.jl
# Checks of solve_dissipative_mode (piece 4d, docs/multimode-piece4-guide.md) at
# saved points of the Mμ = 0.42 cosine branch and its averaged secondary:
#   - the matching residual after the solve is at round-off (relative to the
#     source-only residual r0);
#   - power identity: F_H + F_out = -4π ω0 ∫ dr Im Σ_l conj(R_l) S_l, the work
#     done by the source (ratio 1 up to RK4 and quadrature error);
#   - the column-scaled condition number of M (large would mean a resonance);
#   - the share of the horizon flux carried by l = 0, 2, 4.
# F_H / Ê2 uses the secondary energy of piece 3 (cosine_secondary_energy.csv).
using Printf, Serialization, LinearAlgebra
include(joinpath(@__DIR__, "../../src/adiabatic/drivers/continue_secondary.jl"))

results = joinpath(@__DIR__, "../../runs/adiabatic/multimode_secondary/results")
primary_path = joinpath(@__DIR__, "../../runs/adiabatic/fig2_20260922/results/cosine_branch.jls")

setup = load_secondary_setup(primary_path)
(; primary, secondary_modes, grid, polar, azimuthal, couplings) = setup
secondary_branch = deserialize(joinpath(results, "cosine_secondary.jls"))
dissipative_modes = [(0, 0), (2, 0), (4, 0)]
spin, alpha = primary.spin, primary.alpha
count = length(dissipative_modes)
dissipative_part = 2length(secondary_modes)+1:2length(secondary_modes)+2count

secondary_energy = Dict{Int,Float64}()
for line in Iterators.drop(eachline(joinpath(results, "cosine_secondary_energy.csv")), 1)
    fields = split(line, ",")
    secondary_energy[parse(Int, fields[1])] = parse(Float64, fields[5])
end

# Work done by the source per unit time: -4π ω0 ∫ dr Im Σ_l conj(R_l) S_l.
function source_work(solution, primary_profiles, frequency)
    total = 0.0
    branch_pairs = ((solution.profiles.from_horizon, primary_profiles.from_horizon),
                    (solution.profiles.from_outer, primary_profiles.from_outer))
    for (branch, primary_branch) in branch_pairs
        positions = branch.positions
        # Integrate toward increasing tortoise radius.
        order = issorted(positions) ? eachindex(positions) : reverse(eachindex(positions))
        integrand = map(order) do sample
            state = branch.states[:, sample]
            geometry = tortoise_geometry(positions[sample], spin)
            secondary_values = Dict(
                mode => state[index] for (index, mode) in enumerate(secondary_modes))
            source = dissipative_source(
                primary_values_at(primary_branch, positions[sample], primary.modes),
                secondary_values, geometry.radius, spin, alpha, dissipative_modes,
                polar, azimuthal, grid,
            )
            values = state[dissipative_part][1:count]
            # dr = (dr/dr_*) dr_*.
            -4π * frequency * geometry.dr_drstar * imag(sum(
                conj(values[index]) * source[mode]
                for (index, mode) in enumerate(dissipative_modes)))
        end
        total += integrate_uniform(integrand, abs(positions[2] - positions[1]))
    end
    return total
end

@printf("%4s %-8s %-9s %-7s %-9s %-11s %-10s %-12s %s\n", "idx", "A1", "omega0", "cond",
        "res/r0", "F_H/E2hat", "F_out", "power ratio", "share of F_H in l = 0, 2, 4")
for index in (1, 30, 60, 100, 150, 232)
    (; profiles) = primary_profiles_at(setup, index)
    secondary_solution = secondary_branch.rows[index].solution
    frequency = 2real(primary.continuation[index].solution.frequency) -
                real(secondary_solution.frequency)
    solution = solve_dissipative_mode(
        frequency, secondary_solution, profiles, primary.modes, spin, alpha,
        secondary_modes, dissipative_modes, couplings, polar, azimuthal, grid,
    )
    power_ratio = (solution.horizon_flux + solution.outer_flux) /
                  source_work(solution, profiles, frequency)

    # Each l's part of the horizon flux, from the inner endpoint.
    inner_state = solution.profiles.from_horizon.states[dissipative_part, 1]
    radial_factor = tortoise_geometry(
        first(solution.profiles.from_horizon.positions), spin).radial_factor
    shares = [-4π * frequency * radial_factor *
              imag(conj(inner_state[i]) * inner_state[count+i]) / solution.horizon_flux
              for i in 1:count]

    @printf("%4d %-8.4g %-9.6f %-7.1f %-9.1e %-11.4e %-10.1e %-12.9f %s\n",
            index, primary.continuation[index].amplitude, frequency,
            solution.condition_number, solution.residual_norm / solution.source_residual_norm,
            solution.horizon_flux / secondary_energy[index], solution.outer_flux, power_ratio,
            join([@sprintf("%.2e", share) for share in shares], " "))
    flush(stdout)
end
