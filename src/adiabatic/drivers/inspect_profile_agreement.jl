using Serialization
using LinearAlgebra: norm

include(joinpath(@__DIR__, "../load.jl"))

saved = deserialize(
    joinpath(@__DIR__, "../../../runs/adiabatic/small_amplitude/results/continuation.jls")
)
entry = last(saved.continuation)
solution = entry.solution

# Reconstruct the angular tables used for this saved continuation.
grid = angular_grid(16, 31)
polar_values = build_polar_values(5, grid.cosine_theta_nodes)
azimuthal_values = build_azimuthal_values(5, grid.azimuth_nodes)
couplings = build_cosine_squared_couplings(5, polar_values, grid)

profiles = nonlinear_matching_residual(
    solution.frequency,
    solution.horizon_amplitudes,
    solution.outer_amplitudes,
    saved.spin, saved.alpha, saved.modes, couplings,
    polar_values, azimuthal_values, grid;
    inner=-100.0, outer=100.0, match=5.0, step=0.05,
    save_profile=true,
)

# Both integrations END at the matching point.
horizon_match = profiles.from_horizon.states[:, end]
outer_match = profiles.from_outer.states[:, end]

mode_count = length(saved.modes)
field_jump = horizon_match[1:mode_count] - outer_match[1:mode_count]
slope_jump = horizon_match[mode_count+1:end] - outer_match[mode_count+1:end]

println("A0 = ", entry.amplitude)
println("norm(ΔR)/A0 = ", norm(field_jump) / entry.amplitude)
println("norm(ΔQ)/A0 = ", norm(slope_jump) / entry.amplitude)
println(
    "Reproduces saved residual: ",
    profiles.residual == solution.residual,
)

# Investigate energies
inner_energy = integrate_energy_branch(
    profiles.from_horizon,
    solution.frequency, saved.spin, saved.alpha, saved.modes,
    polar_values, azimuthal_values, grid,
)

outer_energy = integrate_energy_branch(
    profiles.from_outer,
    solution.frequency, saved.spin, saved.alpha, saved.modes,
    polar_values, azimuthal_values, grid,
)

total_energy = inner_energy + outer_energy

println("Inner-branch E/Fa² = ", inner_energy)
println("Outer-branch E/Fa² = ", outer_energy)
println("Total E/Fa² = ", total_energy)
println("E/(Fa² A0²) = ", total_energy / entry.amplitude^2)

# Export branches for plotting
output_path = joinpath(@__DIR__, "../../runs/adiabatic/small_amplitude/results/profiles.csv")

open(output_path, "w") do output
    println(output, "branch,rstar,l,m,R_real,R_imag,Q_real,Q_imag")

    for (label, branch) in (
        ("horizon", profiles.from_horizon),
        ("outer", profiles.from_outer),
    )
        for column in eachindex(branch.positions)
            for (mode_index, (l, m)) in enumerate(saved.modes)
                radial_value = branch.states[mode_index, column]
                radial_slope = branch.states[mode_count + mode_index, column]

                println(output, join((
                    label, branch.positions[column], l, m,
                    real(radial_value), imag(radial_value),
                    real(radial_slope), imag(radial_slope),
                ), ','))
            end
        end
    end
end

println("Saved profiles to ", output_path)