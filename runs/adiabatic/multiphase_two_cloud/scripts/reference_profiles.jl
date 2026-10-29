# Usage: julia --project=. reference_profiles.jl [CAMPAIGN_DIRECTORY [OUTPUT_PREFIX]]
# Run on a compute node. Match each linear cloud's energy to the saved Omiya
# reference, as in Fig. 7. Existing results are preserved.
using Serialization, LinearAlgebra, Printf, SHA, TOML, Dates
include(joinpath(@__DIR__, "../../../../src/adiabatic/primary_secondary/load.jl"))
include(joinpath(@__DIR__, "../../../../src/adiabatic/multiphase/load.jl"))

length(ARGS) <= 2 || error("usage: reference_profiles.jl [CAMPAIGN_DIRECTORY [OUTPUT_PREFIX]]")
campaign_directory = isempty(ARGS) ? normpath(joinpath(@__DIR__, "..")) : abspath(ARGS[1])
output_prefix = length(ARGS) == 2 ? abspath(ARGS[2]) :
    joinpath(campaign_directory, "results/reference_curves_20261028_v2")
for suffix in (".jls", ".csv", ".toml", ".jls.tmp", ".csv.tmp", ".toml.tmp")
    ispath(output_prefix * suffix) && error("Preserve existing output; choose a new prefix: $output_prefix")
end
mkpath(dirname(output_prefix))

# Both historical branches now expose the same points, profiles, and diagnostics.
reference_directory = joinpath(campaign_directory, "work/direct_fig7_20261026_v5/inputs")
primary_path = joinpath(reference_directory, "alpha030_cosine_branch.jls")
secondary_path = joinpath(reference_directory, "alpha030_cosine_secondary.jls")
joint_input_path = joinpath(campaign_directory, "results/direct_fig7_20261026_v5-input.jls")
joint_run_input = deserialize(joint_input_path)
primary_branch = read_branch(primary_path;
    checkpoint_adapter=single_cloud_continuation_checkpoint_adapter)
secondary_branch = read_branch(secondary_path;
    checkpoint_adapter=secondary_continuation_checkpoint_adapter)

primary_point = branch_point(primary_branch, joint_run_input.reference.primary_index)
secondary_point = branch_point(secondary_branch, joint_run_input.reference.secondary_index)
primary_profile, secondary_profile = primary_point.profile, secondary_point.profile
primary_outer_amplitude, secondary_outer_amplitude = joint_run_input.cloud_amplitudes
primary_settings, secondary_settings = primary_branch.settings, secondary_branch.settings
joint_problem = joint_run_input.problem

# A secondary solution belongs to a particular primary background. Its historical
# path can differ from this frozen copy, so check the physical reference as well.
@assert secondary_point.backgrounds.primary.point_id == primary_point.point_id
@assert secondary_point.parameter_values.primary_amplitude == primary_point.parameter_values.amplitude
@assert secondary_point.diagnostics.primary_frequency == only(primary_profile.frequencies)
@assert primary_point.parameter_values.amplitude == primary_outer_amplitude
@assert secondary_point.parameter_values.amplitude == 1
@assert secondary_settings.reading == "averaged"
@assert primary_settings.spin == secondary_settings.spin == joint_problem.spin
@assert primary_settings.alpha == secondary_settings.alpha == joint_problem.alpha
@assert primary_settings.options == secondary_settings.primary_options
@assert Symbol(primary_settings.potential) == joint_problem.potential
@assert secondary_point.diagnostics.residual_norm <= secondary_point.diagnostics.tolerance

# Build angular tables once, covering both clouds.
primary_modes, secondary_modes = primary_settings.modes, secondary_settings.secondary_modes
maximum_l = max(maximum(first, primary_modes), maximum(first, secondary_modes))
grid = angular_grid(primary_settings.n_polar, primary_settings.n_azimuth)
polar_values = build_polar_values(maximum_l, grid.cosine_theta_nodes)
azimuthal_values = build_azimuthal_values(maximum_l, grid.azimuth_nodes)
couplings = build_cosine_squared_couplings(maximum_l, polar_values, grid)
spin, alpha = primary_settings.spin, primary_settings.alpha
(; step, inner, outer, match) = primary_settings.options
@assert (; step, inner, outer, match) ==
    (; joint_problem.integration.step, joint_problem.integration.inner,
       joint_problem.integration.outer, joint_problem.integration.match)
@assert joint_problem.integration.horizon_order == 3 && joint_problem.integration.outer_order == 7

# Omiya primary: reconstruct the saved nonlinear profile at its own amplitude.
omiya_primary_profiles = nonlinear_matching_residual(
    only(primary_profile.frequencies), primary_profile.horizon_coefficients,
    primary_profile.outer_coefficients, spin, alpha, primary_modes,
    couplings, polar_values, azimuthal_values, grid;
    potential=Symbol(primary_settings.potential), save_profile=true,
    primary_settings.options...,
)
@assert maximum(abs, omiya_primary_profiles.residual - primary_point.diagnostics.residual) /
    primary_outer_amplitude <= 1e-10

# Omiya secondary: reconstruct at unit amplitude, then scale the samples once.
# Its equation is linear on the fixed primary; the nonlinear primary is unchanged.
unit_secondary_profiles = secondary_matching_residual(
    only(secondary_profile.frequencies), secondary_profile.horizon_coefficients,
    secondary_profile.outer_coefficients, omiya_primary_profiles, primary_modes,
    spin, alpha, secondary_modes, couplings, polar_values, azimuthal_values, grid;
    save_profile=true, reading=:averaged,
)
@assert unit_secondary_profiles.residual == secondary_point.diagnostics.residual
omiya_secondary_profiles = (
    from_horizon=(positions=unit_secondary_profiles.from_horizon.positions,
        states=secondary_outer_amplitude .* unit_secondary_profiles.from_horizon.states),
    from_outer=(positions=unit_secondary_profiles.from_outer.positions,
        states=secondary_outer_amplitude .* unit_secondary_profiles.from_outer.states),
    residual=secondary_outer_amplitude .* unit_secondary_profiles.residual,
)

# Linear references: start with the previously tested A1/A2 reconstruction.
# After checking the archives below, normalize the linears by ENERGY, not A1/A2.
linear_primary_modes = filter(mode -> last(mode) == 1, primary_modes)
linear_secondary_modes = secondary_modes
linear_primary_solution = solve_linear_mode(
    complex(alpha * (1 - alpha^2 / 8)), spin, alpha,
    first.(linear_primary_modes), 1, couplings;
    step, inner, outer, match, tolerance=1e-8,
)
linear_secondary_solution = solve_linear_mode(
    complex(alpha * (1 - alpha^2 / 18)), spin, alpha,
    first.(linear_secondary_modes), 2, couplings;
    step=2step, inner, outer, match, tolerance=1e-8,
)
@assert linear_primary_solution.residual_norm <= 1e-8
@assert linear_secondary_solution.residual_norm <= 1e-8

# m is the fundamental azimuthal number: use omega2, not 2omega2, for the secondary.
# Retained samples now have the same format as both Omiya reconstruction results.
linear_primary_profiles = matching_residual(
    linear_primary_solution.frequency,
    primary_outer_amplitude .* linear_primary_solution.horizon_amplitudes,
    primary_outer_amplitude .* linear_primary_solution.outer_amplitudes,
    spin, alpha, first.(linear_primary_modes), 1, couplings;
    step, inner, outer, match, save_profile=true,
)
linear_secondary_profiles = matching_residual(
    linear_secondary_solution.frequency,
    secondary_outer_amplitude .* linear_secondary_solution.horizon_amplitudes,
    secondary_outer_amplitude .* linear_secondary_solution.outer_amplitudes,
    spin, alpha, first.(linear_secondary_modes), 2, couplings;
    step=2step, inner, outer, match, save_profile=true,
)

# Compare all radial states against the independently tested reference archives.
# This validates the interface migration, not spatial convergence.
omiya_archive_path = joinpath(campaign_directory, "results/reference_secondary_20261027_v2.jls")
linear_archive_path = joinpath(campaign_directory, "results/reference_linears_20261027_v1.jls")
omiya_archive, linear_archive = deserialize(omiya_archive_path), deserialize(linear_archive_path)
input_sha256 = bytes2hex(sha256(read(joint_input_path)))
@assert omiya_archive.input_sha256 == linear_archive.input_sha256 == input_sha256
@assert omiya_archive.cloud_amplitudes == linear_archive.cloud_amplitudes == joint_run_input.cloud_amplitudes
references = (
    (; name="omiya_primary", modes=primary_modes, archive_amplitude=primary_outer_amplitude,
       radial_result=omiya_primary_profiles, saved_radial_result=omiya_archive.omiya_primary_profiles),
    (; name="omiya_secondary", modes=secondary_modes, archive_amplitude=secondary_outer_amplitude,
       radial_result=omiya_secondary_profiles, saved_radial_result=omiya_archive.omiya_secondary_profiles),
    (; name="linear_primary", modes=linear_primary_modes, archive_amplitude=primary_outer_amplitude,
       radial_result=linear_primary_profiles, saved_radial_result=linear_archive.linear_primary_profiles),
    (; name="linear_secondary", modes=linear_secondary_modes, archive_amplitude=secondary_outer_amplitude,
       radial_result=linear_secondary_profiles, saved_radial_result=linear_archive.linear_secondary_profiles),
)
matching_norms, reference_differences = Dict{String,Float64}(), Dict{String,Float64}()
for reference in references
    radial_result = reference.radial_result
    matching_norm = norm(radial_result.residual) / reference.archive_amplitude
    @assert matching_norm <= joint_run_input.matching_tolerance
    endpoint_mismatch = radial_result.from_horizon.states[:, end] - radial_result.from_outer.states[:, end]
    @assert norm(endpoint_mismatch) / reference.archive_amplitude <= joint_run_input.matching_tolerance
    maximum_difference = 0.0
    for boundary in (:from_horizon, :from_outer)
        radial_branch = getproperty(radial_result, boundary)
        saved_radial_branch = getproperty(reference.saved_radial_result, boundary)
        @assert all(isfinite, radial_branch.states)
        @assert radial_branch.positions == saved_radial_branch.positions
        @assert size(radial_branch.states) == size(saved_radial_branch.states)
        relative_difference = maximum(abs, radial_branch.states - saved_radial_branch.states) /
            maximum(abs, saved_radial_branch.states)
        @assert relative_difference <= 1e-10
        maximum_difference = max(maximum_difference, relative_difference)
    end
    matching_norms[reference.name] = matching_norm
    reference_differences[reference.name] = maximum_difference
    println(reference.name, ": matching=", matching_norm, ", archive difference=", maximum_difference)
end

# Fig. 7 compares equal cloud energies. The saved primary is the available
# nonlinear solution nearest 660; do not rescale a nonlinear solution by hand.
reference_energy = primary_point.diagnostics.observables.energy
@assert reference_energy == joint_run_input.reference.reference_energy
omiya_primary_energy = sum(
    integrate_energy_branch(radial_branch, only(primary_profile.frequencies),
        spin, alpha, primary_modes, polar_values, azimuthal_values, grid;
        potential=Symbol(primary_settings.potential))
    for radial_branch in (omiya_primary_profiles.from_horizon, omiya_primary_profiles.from_outer)
)

# The Omiya secondary uses its quadratic energy on the fixed primary background,
# including the averaged effective mass. Measure the already-scaled samples.
secondary_branch_energies = Float64[]
for boundary in (:from_horizon, :from_outer)
    secondary_radial_branch = getproperty(omiya_secondary_profiles, boundary)
    primary_radial_branch = getproperty(omiya_primary_profiles, boundary)
    positions = secondary_radial_branch.positions
    sample_order = issorted(positions) ? eachindex(positions) : reverse(eachindex(positions))
    densities = [secondary_energy_per_tortoise_radius(
        view(secondary_radial_branch.states, :, sample_index), positions[sample_index],
        primary_values_at(primary_radial_branch, positions[sample_index], primary_modes),
        only(secondary_profile.frequencies), spin, alpha, secondary_modes,
        polar_values, azimuthal_values, grid,
    ) for sample_index in sample_order]
    push!(secondary_branch_energies,
        integrate_uniform(densities, abs(positions[2] - positions[1])))
end
omiya_secondary_energy = sum(secondary_branch_energies)
@assert isapprox(omiya_primary_energy, reference_energy; rtol=1e-9)
@assert isapprox(omiya_secondary_energy, reference_energy; rtol=1e-9)

# With self-interaction absent, U(phi)=phi²/2 and E(c*R)=c²*E(R).
# Scale R, its radial derivative Q, and the matching residual together.
energies_by_curve = Dict("omiya_primary" => omiya_primary_energy,
    "omiya_secondary" => omiya_secondary_energy)
outer_amplitudes_by_curve = Dict("omiya_primary" => primary_outer_amplitude,
    "omiya_secondary" => secondary_outer_amplitude)
linear_energies_before_scaling = Dict{String,Float64}()
linear_scale_factors = Dict{String,Float64}()
for (name, radial_result, solution, modes, original_amplitude) in (
    ("linear_primary", linear_primary_profiles, linear_primary_solution,
        linear_primary_modes, primary_outer_amplitude),
    ("linear_secondary", linear_secondary_profiles, linear_secondary_solution,
        linear_secondary_modes, secondary_outer_amplitude),
)
    fundamental_m = last(first(modes))
    energy_before_scaling = sum(integrate_energy_branch(
        radial_branch, solution.frequency, spin, alpha, modes,
        polar_values, azimuthal_values, grid; fundamental_m, potential=:quadratic,
    ) for radial_branch in (radial_result.from_horizon, radial_result.from_outer))
    @assert isfinite(energy_before_scaling) && energy_before_scaling > 0
    amplitude_factor = sqrt(reference_energy / energy_before_scaling)
    radial_result.from_horizon.states .*= amplitude_factor
    radial_result.from_outer.states .*= amplitude_factor
    radial_result.residual .*= amplitude_factor

    # Reintegrate the energy density after scaling to check the exported data.
    measured_energy = sum(integrate_energy_branch(
        radial_branch, solution.frequency, spin, alpha, modes,
        polar_values, azimuthal_values, grid; fundamental_m, potential=:quadratic,
    ) for radial_branch in (radial_result.from_horizon, radial_result.from_outer))
    @assert isapprox(measured_energy, reference_energy; rtol=1e-12)
    normalized_mismatch = norm(radial_result.residual) / (original_amplitude * amplitude_factor)
    @assert normalized_mismatch <= joint_run_input.matching_tolerance
    energies_by_curve[name] = measured_energy
    outer_amplitudes_by_curve[name] = original_amplitude * amplitude_factor
    linear_energies_before_scaling[name] = energy_before_scaling
    linear_scale_factors[name] = amplitude_factor
    matching_norms[name] = normalized_mismatch
    println(name, ": energy before=", energy_before_scaling, ", amplitude factor=",
        amplitude_factor, ", energy after=", measured_energy)
end

# Export the fundamental modes. Keep BOTH matching-point samples so the small
# mismatch remains visible. The JLS below retains every mode and radial slope.
open(output_prefix * ".csv.tmp", "w") do stream
    println(stream, "curve,ell,m,r_star_over_M,R_real,R_imag,R_abs,radial_branch")
    for reference in references
        ell, m = first(reference.modes)
        for boundary in (:from_horizon, :from_outer)
            radial_branch = getproperty(reference.radial_result, boundary)
            sample_indices = boundary == :from_horizon ? eachindex(radial_branch.positions) :
                reverse(eachindex(radial_branch.positions))
            for sample_index in sample_indices
                radial_value = radial_branch.states[1, sample_index]
                @printf(stream, "%s,%d,%d,%.17g,%.17g,%.17g,%.17g,%s\n",
                    reference.name, ell, m, radial_branch.positions[sample_index],
                    real(radial_value), imag(radial_value), abs(radial_value), string(boundary))
            end
        end
    end
end
mv(output_prefix * ".csv.tmp", output_prefix * ".csv")

# Record input and source hashes. The batch job uses a frozen source copy.
repository_directory = normpath(joinpath(@__DIR__, "../../../.."))
source_hashes = Dict{String,String}()
for (directory, _, filenames) in walkdir(joinpath(repository_directory, "src/adiabatic"))
    for filename in filenames
        endswith(filename, ".jl") || continue
        source_path = joinpath(directory, filename)
        source_hashes[relpath(source_path, repository_directory)] = bytes2hex(sha256(read(source_path)))
    end
end
for source_path in (joinpath(repository_directory, "Project.toml"),
        joinpath(repository_directory, "Manifest.toml"), @__FILE__)
    source_hashes[relpath(source_path, repository_directory)] = bytes2hex(sha256(read(source_path)))
end
metadata = Dict(
    "normalization" => "equal_cloud_energies",
    "reference_energy" => reference_energy,
    "energies_by_curve" => energies_by_curve,
    "outer_amplitudes_by_curve" => outer_amplitudes_by_curve,
    "linear_energies_before_scaling" => linear_energies_before_scaling,
    "linear_scale_factors" => linear_scale_factors,
    "omiya_outer_amplitudes" => joint_run_input.cloud_amplitudes,
    "spin" => spin, "alpha" => alpha,
    "primary_point_id" => primary_point.point_id, "secondary_point_id" => secondary_point.point_id,
    "inner" => inner, "outer" => outer, "match" => match,
    "primary_step" => step, "secondary_step" => 2step,
    "matching_tolerance" => joint_run_input.matching_tolerance,
    "matching_norms" => matching_norms,
    "pre_normalization_archive_relative_differences" => reference_differences,
    "created_utc" => string(now(UTC)), "slurm_job_id" => get(ENV, "SLURM_JOB_ID", "local"),
    "input_sha256" => Dict(basename(path) => bytes2hex(sha256(read(path)))
        for path in (joint_input_path, primary_path, secondary_path, omiya_archive_path, linear_archive_path)),
    "source_sha256" => source_hashes,
)
open(output_prefix * ".toml.tmp", "w") do stream
    TOML.print(stream, metadata; sorted=true)
end
mv(output_prefix * ".toml.tmp", output_prefix * ".toml")
serialize(output_prefix * ".jls.tmp", (;
    omiya_primary_profiles, omiya_secondary_profiles, linear_primary_profiles, linear_secondary_profiles,
    primary_modes, secondary_modes, linear_primary_modes, linear_secondary_modes,
    linear_primary_solution, linear_secondary_solution, metadata,
))
mv(output_prefix * ".jls.tmp", output_prefix * ".jls")
println("Saved four reference curves: ", output_prefix, ".{csv,jls,toml}")
