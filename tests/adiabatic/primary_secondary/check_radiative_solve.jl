# julia --project=. tests/adiabatic/primary_secondary/check_radiative_solve.jl
# Checks of solve_radiative_mode (piece 5d) at saved points of the Mμ = 0.42
# cosine branch and its averaged secondary:
#   - the matching residual after the solve is at round-off (relative to r0);
#   - power identity: F_H + F_inf = -4π ω3 ∫ dr Im Σ_l conj(R_l) S_l, the work
#     done by the source (ratio 1 up to RK4 and quadrature error);
#   - the column-scaled condition number of M;
#   - the share of F_inf carried by l = 3, 5, 7;
#   - the small-A1 law F_inf / Ê2² ∝ A1² (E1 E2² of Eq. (39)): log-slope
#     between successive points tends to 2.
# Fluxes are divided by Ê2² (piece 3 energy per unit outer amplitude squared),
# giving flux per unit secondary energy squared, the E2² of Eq. (39).
using Printf, Serialization, LinearAlgebra
include(joinpath(@__DIR__, "../../../src/adiabatic/primary_secondary/load.jl"))

results = joinpath(@__DIR__, "../../../runs/adiabatic/multimode_secondary/results")
primary_path = joinpath(@__DIR__, "../../../runs/adiabatic/fig2_20260922/results/cosine_branch.jls")

setup = load_secondary_setup(primary_path)
(; primary, secondary_modes, grid) = setup
polar = build_polar_values(7, grid.cosine_theta_nodes)
azimuthal = build_azimuthal_values(7, grid.azimuth_nodes)
couplings = build_cosine_squared_couplings(7, polar, grid)
secondary_branch = deserialize(joinpath(results, "cosine_secondary.jls"))
radiative_modes = [(3, 3), (5, 3), (7, 3)]
spin, alpha = primary.spin, primary.alpha
count = length(radiative_modes)
radiative_part = 2length(secondary_modes)+1:2length(secondary_modes)+2count

secondary_energy = Dict{Int,Float64}()
for line in Iterators.drop(eachline(joinpath(results, "cosine_secondary_energy.csv")), 1)
    fields = split(line, ",")
    secondary_energy[parse(Int, fields[1])] = parse(Float64, fields[5])
end

# Work done by the source per unit time: -4π ω3 ∫ dr Im Σ_l conj(R_l) S_l.
function source_work(solution, primary_profiles, frequency)
    total = 0.0
    branch_pairs = ((solution.profiles.from_horizon, primary_profiles.from_horizon),
                    (solution.profiles.from_outer, primary_profiles.from_outer))
    for (branch, primary_branch) in branch_pairs
        positions = branch.positions
        order = issorted(positions) ? eachindex(positions) : reverse(eachindex(positions))
        integrand = map(order) do sample
            state = branch.states[:, sample]
            geometry = tortoise_geometry(positions[sample], spin)
            secondary_values = Dict(
                mode => state[index] for (index, mode) in enumerate(secondary_modes))
            source = radiative_source(
                primary_values_at(primary_branch, positions[sample], primary.modes),
                secondary_values, geometry.radius, spin, alpha, radiative_modes,
                polar, azimuthal, grid,
            )
            values = state[radiative_part][1:count]
            -4π * frequency * geometry.dr_drstar * imag(sum(
                conj(values[index]) * source[mode]
                for (index, mode) in enumerate(radiative_modes)))
        end
        total += integrate_uniform(integrand, abs(positions[2] - positions[1]))
    end
    return total
end

@printf("%4s %-8s %-9s %-7s %-9s %-11s %-11s %-12s %s\n", "idx", "A1", "omega3", "cond",
        "res/r0", "F_inf/E2^2", "F_H/E2^2", "power ratio", "share of F_inf in l = 3, 5, 7")
fluxes = Dict{Int,Float64}()
for index in (1:6..., 8, 10, 30, 60, 100, 150, 232)
    time = @elapsed begin
        (; profiles) = primary_profiles_at(setup, index)
        secondary_solution = secondary_branch.rows[index].solution
        frequency = 2real(secondary_solution.frequency) -
                    real(primary.continuation[index].solution.frequency)
        solution = solve_radiative_mode(
            frequency, secondary_solution, profiles, primary.modes, spin, alpha,
            secondary_modes, radiative_modes, couplings, polar, azimuthal, grid,
        )
        power_ratio = (solution.horizon_flux + solution.infinity_flux) /
                      source_work(solution, profiles, frequency)
        outer_state = solution.profiles.from_outer.states[radiative_part, 1]
        radial_factor = tortoise_geometry(
            first(solution.profiles.from_outer.positions), spin).radial_factor
        shares = [4π * frequency * radial_factor *
                  imag(conj(outer_state[i]) * outer_state[count+i]) / solution.infinity_flux
                  for i in 1:count]
    end
    energy_squared = secondary_energy[index]^2
    fluxes[index] = solution.infinity_flux / energy_squared
    @printf("%4d %-8.4g %-9.6f %-7.1f %-9.1e %-11.4e %-11.4e %-12.9f %s  (%.0f s)\n",
            index, primary.continuation[index].amplitude, frequency,
            solution.condition_number, solution.residual_norm / solution.source_residual_norm,
            solution.infinity_flux / energy_squared, solution.horizon_flux / energy_squared,
            power_ratio, join([@sprintf("%.2e", share) for share in shares], " "), time)
    flush(stdout)
end

println("
log-slope of F_inf / Ê2² against A1 between successive small-A1 points")
small = (1:6..., 8, 10)
for (previous, index) in zip(small[1:end-1], small[2:end])
    amplitude(i) = primary.continuation[i].amplitude
    @printf("A1 %.4g -> %.4g: %.4f
", amplitude(previous), amplitude(index),
            log(fluxes[index] / fluxes[previous]) / log(amplitude(index) / amplitude(previous)))
end
