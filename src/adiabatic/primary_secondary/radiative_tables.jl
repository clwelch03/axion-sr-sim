using Serialization, Printf, LinearAlgebra

"""
    radiative_source_work(solution, primary_profiles, primary_modes, frequency,
                          spin, alpha, secondary_modes, radiative_modes,
                          polar_values, azimuthal_values, grid)

Return the work done per unit time by the source on the m = 3 field,
`-4π ω3 ∫ dr Im Σ_l conj(R_l) S_l`, over the saved profiles of a
[`solve_radiative_mode`](@ref) solution. The power identity says this equals
`horizon_flux + infinity_flux`.
"""
function radiative_source_work(solution, primary_profiles, primary_modes, frequency,
                               spin, alpha, secondary_modes, radiative_modes,
                               polar_values, azimuthal_values, grid)
    offset = 2length(secondary_modes)
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
            source = radiative_source(
                primary_values_at(primary_branch, positions[sample], primary_modes),
                secondary_values, geometry.radius, spin, alpha, radiative_modes,
                polar_values, azimuthal_values, grid,
            )
            # dr = (dr/dr_*) dr_*.
            -4π * frequency * geometry.dr_drstar * imag(sum(
                conj(state[offset+index]) * source[mode]
                for (index, mode) in enumerate(radiative_modes)))
        end
        total += integrate_uniform(integrand, abs(positions[2] - positions[1]))
    end
    return total
end

"""
    infinity_flux_shares(solution, frequency, spin, secondary_modes, radiative_modes)

Return each m = 3 mode's fraction of the flux to infinity, from the current at
the primary's outer boundary.
"""
function infinity_flux_shares(solution, frequency, spin, secondary_modes, radiative_modes)
    count = length(radiative_modes)
    offset = 2length(secondary_modes)
    outer_state = solution.profiles.from_outer.states[:, 1]
    radial_factor = tortoise_geometry(
        first(solution.profiles.from_outer.positions), spin).radial_factor
    return [4π * frequency * radial_factor *
            imag(conj(outer_state[offset+i]) * outer_state[offset+count+i]) /
            solution.infinity_flux
            for i in 1:count]
end

function radiative_fluxes(primary_path, secondary_path, energy_path, output_path;
                          maximum_l=7, maximum_points=nothing)
    isodd(maximum_l) && maximum_l >= 3 || error("maximum_l must be odd and at least 3")
    radiative_modes = [(l, 3) for l in 3:2:maximum_l]
    setup = load_secondary_setup(primary_path)
    (; primary, secondary_modes, grid) = setup
    # Angular tables covering the m = 3 modes as well as the primary and secondary.
    table_l = max(maximum_l, maximum(first, primary.modes), maximum(first, secondary_modes))
    polar = build_polar_values(table_l, grid.cosine_theta_nodes)
    azimuthal = build_azimuthal_values(table_l, grid.azimuth_nodes)
    couplings = build_cosine_squared_couplings(table_l, polar, grid)

    secondary_branch = deserialize(secondary_path)
    get(secondary_branch, :reading, "fixed_time") == "averaged" ||
        error("the m = 3 field is defined for the averaged reading only")
    secondary_energies = read_secondary_energies(energy_path)

    rows = secondary_branch.rows
    last_index = maximum_points === nothing ? length(rows) : min(length(rows), maximum_points)
    share_columns = join(["share_l$(l)" for (l, _) in radiative_modes], ",")

    temporary = output_path * ".tmp"
    open(temporary, "w") do io
        println(io, "primary_index,amplitude,omega1_real,omega2_real,omega3,",
                    "infinity_flux,horizon_flux,secondary_energy,",
                    "infinity_flux_per_energy_squared,horizon_flux_per_energy_squared,",
                    "power_ratio,condition_number,relative_residual,", share_columns)
        for row in rows[1:last_index]
            (; profiles) = primary_profiles_at(setup, row.primary_index)
            omega1 = real(primary.continuation[row.primary_index].solution.frequency)
            omega2 = real(row.solution.frequency)
            # Footnote 2: real parts only.
            frequency = 2omega2 - omega1

            elapsed = @elapsed solution = solve_radiative_mode(
                frequency, row.solution, profiles, primary.modes, primary.spin,
                primary.alpha, secondary_modes, radiative_modes, couplings,
                polar, azimuthal, grid,
            )
            work = radiative_source_work(
                solution, profiles, primary.modes, frequency, primary.spin, primary.alpha,
                secondary_modes, radiative_modes, polar, azimuthal, grid,
            )
            power_ratio = (solution.horizon_flux + solution.infinity_flux) / work
            shares = infinity_flux_shares(
                solution, frequency, primary.spin, secondary_modes, radiative_modes)
            secondary_energy = secondary_energies[row.primary_index]

            @printf(io, "%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.6g,%.6g,%s\n",
                    row.primary_index, row.amplitude, omega1, omega2, frequency,
                    solution.infinity_flux, solution.horizon_flux, secondary_energy,
                    solution.infinity_flux / secondary_energy^2,
                    solution.horizon_flux / secondary_energy^2, power_ratio,
                    solution.condition_number,
                    solution.residual_norm / solution.source_residual_norm,
                    join([@sprintf("%.6g", share) for share in shares], ","))
            flush(io)
            @printf("%3d A1=%-9.5g omega3=%.6f F_inf/E2^2=%.6e F_H/F_inf=%.1e power=%.9f cond=%.1f t=%.1fs\n",
                    row.primary_index, row.amplitude, frequency,
                    solution.infinity_flux / secondary_energy^2,
                    solution.horizon_flux / solution.infinity_flux, power_ratio,
                    solution.condition_number, elapsed)
            solution.condition_number > 1e8 &&
                @warn "matching matrix nearly singular" row.primary_index solution.condition_number
            flush(stdout)
        end
    end
    mv(temporary, output_path; force=true)
end
