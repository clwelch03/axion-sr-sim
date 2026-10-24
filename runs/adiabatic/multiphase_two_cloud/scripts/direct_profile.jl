using Serialization, LinearAlgebra
include(joinpath(@__DIR__, "../../../../src/adiabatic/multiphase/load.jl"))

# Time each expensive residual call, including the passes used for derivatives.
# ForwardDiff sends Dual numbers through the same residual to build the Jacobian.
# Only the printed norm discards derivatives; Newton receives the original result.
function logged_profile_newton(residual, initial; tolerance, max_iterations,
        admissible, log_io=stdout)
    solve_start = time_ns()
    current_iteration = -1
    jacobian_evaluations = 0
    trial_evaluations = 0
    elapsed_seconds(start) = round((time_ns() - start) / 1e9; digits=2)
    function report(message)
        println(log_io, "[", elapsed_seconds(solve_start), " s] ", message)
        flush(log_io)
    end

    function logged_residual(state)
        if eltype(state) <: ForwardDiff.Dual
            jacobian_evaluations += 1
            stage = "Correction $(current_iteration + 1), Jacobian evaluation $jacobian_evaluations"
        elseif current_iteration < 0
            stage = "Initial residual"
        else
            trial_evaluations += 1
            stage = "Correction $(current_iteration + 1), trial residual $trial_evaluations"
        end
        report("$stage: starting")
        evaluation_start = time_ns()
        try
            mismatch = residual(state)
            printed_norm = norm(ForwardDiff.value.(mismatch))
            report("$stage: finished in $(elapsed_seconds(evaluation_start)) s; scaled residual = $printed_norm")
            return mismatch
        catch
            report("$stage: failed after $(elapsed_seconds(evaluation_start)) s")
            rethrow()
        end
    end

    # This callback distinguishes accepted iterates from the trial evaluations.
    function progress_callback(iteration, residual_norm)
        current_iteration = iteration
        jacobian_evaluations = 0
        trial_evaluations = 0
        report("Newton iteration $iteration: scaled residual = $residual_norm")
    end

    solution = newton_solve(logged_residual, initial; tolerance, max_iterations,
        admissible, progress_callback, jacobian_backend=:forwarddiff)
    report("Newton converged after $(solution.iterations) corrections")
    return solution
end

# Choose the reference directory and optional solve output
# No arguments: check the setup. Two paths: solve using those inputs and output prefix.
# Heavy solves use the existing Slurm script and a fresh frozen source snapshot.
if isempty(ARGS)
    reference_directory = normpath(joinpath(@__DIR__, "../../multimode_secondary/results"))
    output_prefix = nothing
elseif length(ARGS) == 2 && ARGS[1] == "--check"
    reference_directory = abspath(ARGS[2])
    output_prefix = nothing
elseif length(ARGS) == 2 && !startswith(ARGS[1], "--")
    reference_directory = abspath(ARGS[1])
    output_prefix = abspath(ARGS[2])
else
    error("usage: direct_profile.jl [--check INPUT_DIRECTORY | INPUT_DIRECTORY OUTPUT_PREFIX]")
end

# Load the saved primary and secondary branches
primary_branch = deserialize(joinpath(reference_directory, "alpha030_cosine_branch.jls"))
secondary_branch = deserialize(joinpath(reference_directory, "alpha030_cosine_secondary.jls"))

@show length(primary_branch.continuation) length(secondary_branch.rows)

# Select the saved primary solution closest to 660*F_a^2*M
target_energy = 660.0
primary_energies = [row.energy for row in primary_branch.rows]
primary_index = argmin(abs.(primary_energies .- target_energy))

primary_entry = primary_branch.continuation[primary_index]
reference_energy = primary_energies[primary_index]
primary_amplitude = primary_entry.amplitude

@show primary_index reference_energy primary_amplitude

# Select the secondary solution computed on that primary
secondary_index = findfirst(
    row -> row.primary_index == primary_index,
    secondary_branch.rows,
)
secondary_index === nothing && error("No secondary solution for the selected primary")

secondary_entry = secondary_branch.rows[secondary_index]
@assert secondary_entry.amplitude == primary_amplitude

@show secondary_index secondary_entry.amplitude

# Scale the secondary's energy
secondary_energy_lines = readlines(joinpath(
    reference_directory, "alpha030_cosine_secondary_energy.csv",
))
secondary_energy_header = split(first(secondary_energy_lines), ',')
secondary_energy_rows = [split(line, ',') for line in secondary_energy_lines[2:end]]

primary_index_column = findfirst(==("primary_index"), secondary_energy_header)
energy_column = findfirst(==("energy"), secondary_energy_header)

secondary_energy_row = only(filter(
    row -> parse(Int, row[primary_index_column]) == primary_index,
    secondary_energy_rows,
))
secondary_energy_per_amplitude_squared = parse(Float64, secondary_energy_row[energy_column])
secondary_amplitude = sqrt(reference_energy / secondary_energy_per_amplitude_squared)

@show secondary_energy_per_amplitude_squared secondary_amplitude

# Specify modes
cloud_m = [1, 2]

components = [
    fourier_component((1, 0), [1, 3, 5], cloud_m),
    fourier_component((3, 0), [3, 5], cloud_m),
    fourier_component((5, 0), [5], cloud_m),
    fourier_component((0, 1), [2, 4, 6], cloud_m),
    fourier_component((2, -1), [0, 2, 4], cloud_m),
    fourier_component((-1, 2), [3, 5], cloud_m),
]

@show sum(length(component.l_values) for component in components)

# Construct initial guess for the joint profile
primary_solution = primary_entry.solution
secondary_solution = secondary_entry.solution
frequency_guess = [primary_solution.frequency, secondary_solution.frequency]

horizon_guess = vcat(
    primary_solution.horizon_amplitudes,
    secondary_amplitude .* secondary_solution.horizon_amplitudes,
    zeros(ComplexF64, 5),
)
outer_guess = vcat(
    primary_solution.outer_amplitudes,
    secondary_amplitude .* secondary_solution.outer_amplitudes,
    zeros(ComplexF64, 5),
)

profile_guess = TwoPhaseProfile(frequency_guess, horizon_guess, outer_guess)

@show length(horizon_guess) length(outer_guess)
@show outer_guess[[1, 7]]

# Build grids and angular tables
polar_node_count = primary_branch.n_polar
phase_node_count = primary_branch.n_azimuth
angular_quadrature = angular_grid(polar_node_count, phase_node_count)

phase1_nodes = copy(angular_quadrature.azimuth_nodes)
phase2_nodes = copy(angular_quadrature.azimuth_nodes)

maximum_l = maximum(maximum(component.l_values) for component in components)
polar_values = build_polar_values(maximum_l, angular_quadrature.cosine_theta_nodes)
couplings = build_cosine_squared_couplings(maximum_l, polar_values, angular_quadrature)

@show polar_node_count length(phase1_nodes) length(phase2_nodes) maximum_l

# Create the TwoPhaseProblem instance from the above fixed settings
cloud_amplitudes = [primary_amplitude, secondary_amplitude]
coefficient_indices = two_phase_coefficient_indices(components)

initial_bound = BitVector([
    real(sum(component.q .* frequency_guess)) < primary_branch.alpha
    for component in components
])

integration_settings = (; primary_branch.options..., horizon_order=3, outer_order=7)

problem = TwoPhaseProblem(;
    spin=primary_branch.spin, alpha=primary_branch.alpha,
    components, cloud_m, couplings, polar_values,
    grid=angular_quadrature, phase1_nodes, phase2_nodes,
    potential=:cosine,
    coefficient_scale=minimum(cloud_amplitudes),
    free_outer_coefficient_indices=coefficient_indices.free_outer_coefficient_indices,
    initial_bound,
    integration=integration_settings,
)

@show coefficient_indices.fixed_outer_coefficient_indices initial_bound problem.integration

# Check that the concatenated boundary guesses follow the chosen mode ordering
radial_modes = [(l, component.q) for component in components for l in component.l_values]
@assert radial_modes[1:6] == [(l, (m, 0)) for (l, m) in primary_branch.modes]
@assert radial_modes[7:9] == [(l, (0, 1)) for (l, m) in secondary_branch.secondary_modes]
@assert length(horizon_guess) == length(outer_guess) == length(radial_modes)
@assert outer_guess[coefficient_indices.fixed_outer_coefficient_indices] == cloud_amplitudes

# Pack the physical guess into Newton's real coordinates and check the conversion
packed_guess = pack_two_phase_profile(profile_guess, problem)
recovered_profile = unpack_two_phase_profile(packed_guess, cloud_amplitudes, problem)

@assert recovered_profile.frequencies ≈ profile_guess.frequencies
@assert recovered_profile.horizon_coefficients ≈ profile_guess.horizon_coefficients
@assert recovered_profile.outer_coefficients ≈ profile_guess.outer_coefficients
@assert two_phase_profile_admissible(packed_guess, problem)

@show length(packed_guess)

# Set the Newton stopping criteria and describe the normalization being compared
# The tolerance applies to the matching residual divided by coefficient_scale.
# The saved reference secondary has a residual of about 3.2e-6 per unit outer
# coefficient; 1e-5 is the requested tolerance for this preliminary joint attempt.
# This is a matching tolerance, not an estimate of spatial or mode-truncation error.
matching_tolerance = 1e-5
maximum_newton_iterations = 12
normalization_description = "Fixed outer coefficients from the equal-energy Omiya reference; joint cloud energies are not constrained."

@show matching_tolerance maximum_newton_iterations
println(normalization_description)
flush(stdout)

if output_prefix === nothing
    println("Setup checks passed. Supply INPUT_DIRECTORY OUTPUT_PREFIX to run the solve from a frozen snapshot.")
else
    # Read the frozen-source record and reserve a fresh output prefix
    source_root = normpath(joinpath(@__DIR__, "../../../.."))
    provenance = read(joinpath(source_root, "provenance.json"), String)
    output_suffixes = ("-input.jls", "-input.jls.tmp", "-solution.jls", "-solution.jls.tmp",
        "-profiles.csv", "-profiles.csv.tmp", "-failure.txt", ".lock")
    any(ispath(output_prefix * suffix) for suffix in output_suffixes) &&
        error("Output prefix already exists; choose a new run name")
    mkpath(dirname(output_prefix))
    mkdir(output_prefix * ".lock")

    try
        # Save the starting problem, guess, reference normalization, and source record
        println("Saving the initial matching problem...")
        flush(stdout)
        reference = (; primary_index, secondary_index, reference_energy,
            secondary_energy_per_amplitude_squared, normalization_description)
        serialize(output_prefix * "-input.jls.tmp", (; problem, profile_guess,
            cloud_amplitudes, packed_guess, reference, matching_tolerance,
            maximum_newton_iterations, provenance))
        mv(output_prefix * "-input.jls.tmp", output_prefix * "-input.jls")

        # Use the same residual, admissibility check, and Newton settings as
        # correct_two_phase_profile, with timing around each residual evaluation.
        matching_residual = state -> two_phase_real_residual(state, cloud_amplitudes, problem)
        admissible = state -> two_phase_profile_admissible(state, problem)
        newton_solution = logged_profile_newton(matching_residual, packed_guess;
            tolerance=matching_tolerance, max_iterations=maximum_newton_iterations,
            admissible)
        solution = (; state=newton_solution.state, amplitudes=copy(cloud_amplitudes),
            residual_norm=norm(newton_solution.residual), iterations=newton_solution.iterations)

        # Unpack the converged coefficients and integrate again to retain radial samples
        println("Integrating the matched profiles again to save radial samples...")
        flush(stdout)
        matched_profile = unpack_two_phase_profile(solution.state, cloud_amplitudes, problem)
        radial_result = two_phase_matching_residual(matched_profile, problem;
            save_radial_result=true)
        saved_residual_norm = norm(radial_result.residual) / problem.coefficient_scale
        @assert isfinite(saved_residual_norm) && saved_residual_norm <= matching_tolerance

        @show matched_profile.frequencies solution.iterations saved_residual_norm
        flush(stdout)

        # Save the complete joint solution and both radial integration histories
        serialize(output_prefix * "-solution.jls.tmp", (; problem, reference, solution,
            profile=matched_profile, radial_result, cloud_amplitudes,
            matching_tolerance, saved_residual_norm, provenance))
        mv(output_prefix * "-solution.jls.tmp", output_prefix * "-solution.jls")

        # Export the primary and secondary fundamental profiles for the overlay
        primary_radial_index = findfirst(==((1, (1, 0))), radial_modes)
        secondary_radial_index = findfirst(==((2, (0, 1))), radial_modes)
        open(output_prefix * "-profiles.csv.tmp", "w") do io
            println(io, "branch,tortoise_radius,primary_real_R11,primary_imag_R11,primary_abs_R11,secondary_real_R22,secondary_imag_R22,secondary_abs_R22")
            for (branch_name, branch) in
                    (("horizon", radial_result.from_horizon), ("outer", radial_result.from_outer))
                for sample_index in eachindex(branch.positions)
                    primary_value = branch.states[primary_radial_index, sample_index]
                    secondary_value = branch.states[secondary_radial_index, sample_index]
                    println(io, join((branch_name, branch.positions[sample_index],
                        real(primary_value), imag(primary_value), abs(primary_value),
                        real(secondary_value), imag(secondary_value), abs(secondary_value)), ','))
                end
            end
        end
        mv(output_prefix * "-profiles.csv.tmp", output_prefix * "-profiles.csv")
        println("Saved ", output_prefix, "-solution.jls and -profiles.csv")
    catch exception
        # Preserve the error alongside the saved starting state if this attempt fails
        open(output_prefix * "-failure.txt", "w") do io
            showerror(io, exception, catch_backtrace())
        end
        rethrow()
    finally
        # Release the one-writer lock; all starting data and results remain on disk
        rm(output_prefix * ".lock")
    end
end
