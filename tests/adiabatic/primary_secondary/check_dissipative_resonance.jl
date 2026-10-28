# julia --project=. tests/adiabatic/primary_secondary/check_dissipative_resonance.jl
# Is the growing condition number of the m = 0 matching matrix at Mμ = 0.3 a
# resonance at ω0? (piece 4e, docs/multimode-piece4-guide.md, Section 8.)
#   1. Scan the drive frequency at fixed primary, secondary, and source. A
#      resonance shows as a peak in F_H and a dip in σ_min/σ_max of the
#      column-scaled matrix at the same frequency.
#   2. Print the right singular vector of σ_min: the combination of the six
#      column-scaled amplitudes (horizon l = 0, 2, 4 | outer l = 0, 2, 4) that
#      the matrix nearly maps to zero. A true near-resonance needs horizon and
#      outer parts together.
using Printf, Serialization, LinearAlgebra
include(joinpath(@__DIR__, "../../../src/adiabatic/primary_secondary/load.jl"))

results = joinpath(@__DIR__, "../../../runs/adiabatic/multimode_secondary/results")
setup = load_primary_secondary_setup(joinpath(results, "alpha030_cosine_branch.jls"))
(; primary, secondary_modes, grid, polar, azimuthal, couplings) = setup
secondary_branch = deserialize(joinpath(results, "alpha030_cosine_secondary.jls"))
energies = read_secondary_energies(joinpath(results, "alpha030_cosine_secondary_energy.csv"))
count = length(DISSIPATIVE_MODES)

# SVD of the matching matrix with columns scaled to unit norm.
function scaled_matching_svd(frequency, secondary, profiles)
    columns = [
        dissipative_matching_residual(
            unit[1:count], unit[count+1:end], frequency, secondary, profiles,
            primary.modes, primary.spin, primary.alpha, secondary_modes,
            DISSIPATIVE_MODES, couplings, polar, azimuthal, grid; include_source=false,
        )
        for unit in ([i == j ? 1.0 + 0im : 0im for j in 1:2count] for i in 1:2count)
    ]
    matrix = hcat(columns...)
    return svd(matrix ./ transpose([norm(column) for column in eachcol(matrix)]))
end

for index in (1, 83)
    (; profiles) = primary_profiles_at(setup, index)
    secondary = secondary_branch.rows[index].solution
    actual = 2real(primary.continuation[index].solution.frequency) - real(secondary.frequency)
    @printf("\nA1 = %.4g, actual omega0 = %.6f\n", primary.continuation[index].amplitude, actual)
    @printf("%-10s %-14s %s\n", "omega", "F_H/E2hat", "sigma_min/sigma_max")
    for frequency in sort(vcat(collect(0.280:0.001:0.2995), [actual]))
        solution = solve_dissipative_mode(
            frequency, secondary, profiles, primary.modes, primary.spin, primary.alpha,
            secondary_modes, DISSIPATIVE_MODES, couplings, polar, azimuthal, grid,
        )
        singular_values = scaled_matching_svd(frequency, secondary, profiles).S
        @printf("%-10.6f %-14.6e %.3e%s\n", frequency,
                solution.horizon_flux / energies[index],
                minimum(singular_values) / maximum(singular_values),
                frequency == actual ? "   <- actual omega0" : "")
        flush(stdout)
    end
    decomposition = scaled_matching_svd(actual, secondary, profiles)
    @printf("near-null combination at omega0, |entries| (H l=0,2,4 | O l=0,2,4): %s\n",
            join([@sprintf("%.3f", abs(entry)) for entry in decomposition.V[:, end]], " "))
end
