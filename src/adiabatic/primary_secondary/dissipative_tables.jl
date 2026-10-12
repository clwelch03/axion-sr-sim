using Serialization, Printf, LinearAlgebra

const DISSIPATIVE_MODES = [(0, 0), (2, 0), (4, 0)]

"""
    dissipative_source_work(solution, primary_profiles, primary_modes, frequency,
                            spin, alpha, secondary_modes, dissipative_modes,
                            polar_values, azimuthal_values, grid)

Return the work done per unit time by the source on the m = 0 field,
`-4π ω0 ∫ dr Im Σ_l conj(R_l) S_l`, over the saved profiles of a
[`solve_dissipative_mode`](@ref) solution. The power identity says this equals
`horizon_flux + outer_flux`.
"""
function dissipative_source_work(solution, primary_profiles, primary_modes, frequency,
                                 spin, alpha, secondary_modes, dissipative_modes,
                                 polar_values, azimuthal_values, grid)
    count = length(dissipative_modes)
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
            source = dissipative_source(
                primary_values_at(primary_branch, positions[sample], primary_modes),
                secondary_values, geometry.radius, spin, alpha, dissipative_modes,
                polar_values, azimuthal_values, grid,
            )
            # dr = (dr/dr_*) dr_*.
            -4π * frequency * geometry.dr_drstar * imag(sum(
                conj(state[offset+index]) * source[mode]
                for (index, mode) in enumerate(dissipative_modes)))
        end
        total += integrate_uniform(integrand, abs(positions[2] - positions[1]))
    end
    return total
end

"""
    horizon_flux_shares(solution, frequency, spin, secondary_modes, dissipative_modes)

Return each m = 0 mode's fraction of the horizon flux, from the inner endpoint.
"""
function horizon_flux_shares(solution, frequency, spin, secondary_modes, dissipative_modes)
    count = length(dissipative_modes)
    offset = 2length(secondary_modes)
    inner_state = solution.profiles.from_horizon.states[:, 1]
    radial_factor = tortoise_geometry(
        first(solution.profiles.from_horizon.positions), spin).radial_factor
    return [-4π * frequency * radial_factor *
            imag(conj(inner_state[offset+i]) * inner_state[offset+count+i]) /
            solution.horizon_flux
            for i in 1:count]
end

"""
    read_secondary_energies(path)

Read `Ê2` per primary index from the CSV written by `secondary_energies.jl`.
"""
function read_secondary_energies(path)
    energies = Dict{Int,Float64}()
    for line in Iterators.drop(eachline(path), 1)
        fields = split(line, ",")
        energies[parse(Int, fields[1])] = parse(Float64, fields[5])
    end
    return energies
end

function dissipative_fluxes(primary_path, secondary_path, energy_path, output_path;
                            dissipative_modes=DISSIPATIVE_MODES, maximum_points=nothing)
    setup = load_secondary_setup(primary_path)
    (; primary, secondary_modes, grid, polar, azimuthal, couplings) = setup
    secondary_branch = deserialize(secondary_path)
    get(secondary_branch, :reading, "fixed_time") == "averaged" ||
        error("the m = 0 field is defined for the averaged reading only")
    secondary_energies = read_secondary_energies(energy_path)

    rows = secondary_branch.rows
    last_index = maximum_points === nothing ? length(rows) : min(length(rows), maximum_points)
    share_columns = join(["share_l$(l)" for (l, _) in dissipative_modes], ",")

    temporary = output_path * ".tmp"
    open(temporary, "w") do io
        println(io, "primary_index,amplitude,omega1_real,omega2_real,omega0,",
                    "horizon_flux,outer_flux,secondary_energy,horizon_flux_per_energy,",
                    "power_ratio,condition_number,relative_residual,", share_columns)
        for row in rows[1:last_index]
            (; profiles) = primary_profiles_at(setup, row.primary_index)
            omega1 = real(primary.continuation[row.primary_index].solution.frequency)
            omega2 = real(row.solution.frequency)
            # Footnote 2: real parts only.
            frequency = 2omega1 - omega2

            elapsed = @elapsed solution = solve_dissipative_mode(
                frequency, row.solution, profiles, primary.modes, primary.spin,
                primary.alpha, secondary_modes, dissipative_modes, couplings,
                polar, azimuthal, grid,
            )
            work = dissipative_source_work(
                solution, profiles, primary.modes, frequency, primary.spin, primary.alpha,
                secondary_modes, dissipative_modes, polar, azimuthal, grid,
            )
            power_ratio = (solution.horizon_flux + solution.outer_flux) / work
            shares = horizon_flux_shares(
                solution, frequency, primary.spin, secondary_modes, dissipative_modes)
            secondary_energy = secondary_energies[row.primary_index]

            @printf(io, "%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.6g,%.6g,%s\n",
                    row.primary_index, row.amplitude, omega1, omega2, frequency,
                    solution.horizon_flux, solution.outer_flux, secondary_energy,
                    solution.horizon_flux / secondary_energy, power_ratio,
                    solution.condition_number,
                    solution.residual_norm / solution.source_residual_norm,
                    join([@sprintf("%.6g", share) for share in shares], ","))
            flush(io)
            @printf("%3d A1=%-9.5g omega0=%.6f F_H/E2=%.6e power=%.8f cond=%.1f t=%.1fs\n",
                    row.primary_index, row.amplitude, frequency,
                    solution.horizon_flux / secondary_energy, power_ratio,
                    solution.condition_number, elapsed)
            # Shooting lines up the outer l = 2, 4 solutions as they cross the
            # cloud (docs/multimode-piece4-guide.md, Section 8); 4e3 at the end
            # of the Mμ = 0.3 branch. Far above that, re-orthogonalize them.
            solution.condition_number > 1e8 &&
                @warn "matching matrix nearly singular" row.primary_index solution.condition_number
            flush(stdout)
        end
    end
    mv(temporary, output_path; force=true)
end
