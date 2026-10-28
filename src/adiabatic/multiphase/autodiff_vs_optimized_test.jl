# Standalone verification experiment; not included by the library loader.
# Run on Slurm: initial|seed|matched INPUT_JLS OUTPUT_PREFIX [SEED_JLS]
using Serialization, LinearAlgebra, SHA
include(joinpath(@__DIR__, "load.jl"))

struct OrdinarySourceDirection end
function nonlinear_two_phase_source(
        parameters::AbstractVector{ForwardDiff.Dual{OrdinarySourceDirection,Float64,1}},
        source_basis::TwoPhaseSourceBasis; potential=:cosine)
    return invoke(nonlinear_two_phase_source, Tuple{Any,TwoPhaseSourceBasis},
        parameters, source_basis; potential)
end

function verification_setup(direct_input; secondary_fraction=1e-3)
    problem = direct_input.problem
    cloud_amplitudes = direct_input.cloud_amplitudes .* [1.0, secondary_fraction]
    mode_labels = problem.source_basis.mode_labels
    coefficient_factors = [secondary_fraction^abs(q[2]) for (l, q) in mode_labels]
    weak_profile = TwoPhaseProfile(copy(direct_input.profile_guess.frequencies),
        coefficient_factors .* direct_input.profile_guess.horizon_coefficients,
        coefficient_factors .* direct_input.profile_guess.outer_coefficients)
    initial_state = pack_two_phase_profile(weak_profile, problem)
    residual_weights = repeat([q[2] == 0 ? 1.0 : secondary_fraction
        for (l, q) in mode_labels], 4)
    matching_residual = state -> two_phase_real_matching_residual(
        state, cloud_amplitudes, problem) ./ residual_weights
    return (; problem, cloud_amplitudes, secondary_fraction, initial_state,
        residual_weights, matching_residual)
end

function ordinary_ad_column(residual, state, column_index)
    @assert 1 <= column_index <= length(state)
    tagged_state = [ForwardDiff.Dual{OrdinarySourceDirection}(
        state[index], (index == column_index ? 1.0 : 0.0,)) for index in eachindex(state)]
    tagged_residual = residual(tagged_state)
    return [ForwardDiff.partials(value)[1] for value in tagged_residual]
end

# Each file is replaced only by the single job holding this prefix's lock.
function save_verification(path, result)
    serialize(path * ".tmp", result)
    mv(path * ".tmp", path; force=true)
end

function compare_ad_jacobians(residual, state; relative_tolerance=1e-8, output_path=nothing)
    @assert isfinite(relative_tolerance) && relative_tolerance > 0
    println("Building the optimized Jacobian..."); flush(stdout)
    optimized_seconds = @elapsed optimized_jacobian = residual_jacobian(residual, state;
        backend=:forwarddiff, chunk_size=length(state))
    println("Optimized Jacobian seconds: ", optimized_seconds); flush(stdout)
    ordinary_jacobian = fill(NaN, size(optimized_jacobian))
    absolute_errors, relative_errors = fill(NaN, length(state)), fill(NaN, length(state))
    column_seconds = zeros(length(state))
    result = nothing
    for column_index in eachindex(state)
        println("Ordinary AD column ", column_index, "/", length(state)); flush(stdout)
        column_seconds[column_index] = @elapsed ordinary_column =
            ordinary_ad_column(residual, state, column_index)
        ordinary_jacobian[:, column_index] = ordinary_column
        optimized_column = optimized_jacobian[:, column_index]
        absolute_errors[column_index] = norm(ordinary_column - optimized_column)
        column_scale = max(norm(ordinary_column), norm(optimized_column))
        relative_errors[column_index] = column_scale == 0 ? 0.0 :
            absolute_errors[column_index] / column_scale
        passed = all(relative_errors .<= relative_tolerance)
        result = (; state=copy(state), optimized_jacobian, ordinary_jacobian,
            absolute_errors, relative_errors, relative_tolerance, passed,
            completed_columns=column_index, optimized_seconds, column_seconds)
        output_path === nothing || save_verification(output_path, result)
        println("  absolute=", absolute_errors[column_index], " relative=",
            relative_errors[column_index], " seconds=", column_seconds[column_index]); flush(stdout)
    end
    return result
end

function check_finite_differences(residual, state, jacobian, admissible;
        steps=10.0 .^ (-9:-1:-15), relative_tolerance=1e-5, output_path=nothing)
    checks = []
    for pattern in (cos.(1:length(state)), sin.(1:length(state)))
        direction = max.(1.0, abs.(state)) .* (pattern ./ norm(pattern))
        predicted = jacobian * direction
        differences = []
        for step in steps
            plus, minus = state + step * direction, state - step * direction
            if !(admissible(plus) && admissible(minus))
                push!(differences, (; step, admissible=false, relative_error=NaN,
                    realized_relative_error=NaN, direction_rounding=NaN))
                continue
            end
            difference = (residual(plus) - residual(minus)) / (2step)
            # Also compare against the displacement actually representable in Float64.
            realized_direction = (plus - minus) / (2step)
            realized_prediction = jacobian * realized_direction
            relative_error = norm(difference - predicted) / norm(predicted)
            realized_relative_error = norm(difference - realized_prediction) / norm(realized_prediction)
            direction_rounding = norm(realized_prediction - predicted) / norm(predicted)
            push!(differences, (; step, admissible=true, relative_error,
                realized_relative_error, direction_rounding))
            println("FD direction ", length(checks)+1, " step=", step,
                " relative=", relative_error, " realized=", realized_relative_error,
                " rounding=", direction_rounding); flush(stdout)
        end
        # Nominal centered differences retain the original acceptance criterion.
        passed = any(row.admissible && row.relative_error <= relative_tolerance for row in differences)
        push!(checks, (; direction, differences, passed))
        output_path === nothing || save_verification(output_path,
            (; state, checks, relative_tolerance, completed_directions=length(checks),
                passed=length(checks)==2 && all(check.passed for check in checks)))
    end
    return (; state, checks, relative_tolerance, completed_directions=length(checks),
        passed=all(check.passed for check in checks))
end

function run_weak_seed(setup, output_prefix)
    (; problem, cloud_amplitudes, initial_state, matching_residual, residual_weights) = setup
    admissible = state -> two_phase_profile_admissible(state, problem)
    accepted_state_callback = (iteration, state, mismatch) -> begin
        save_verification(output_prefix * "-accepted.jls", (; iteration, state, mismatch))
        println("Newton ", iteration, " weighted residual=", norm(mismatch)); flush(stdout)
    end
    solution = newton_solve(matching_residual, initial_state; tolerance=1e-5,
        max_iterations=20, admissible, accepted_state_callback,
        jacobian_backend=:forwarddiff, jacobian_chunk_size=length(initial_state))
    profile = unpack_two_phase_profile(solution.state, cloud_amplitudes, problem)
    radial_result = two_phase_matching_residual(profile, problem; save_radial_result=true)
    mismatch = pack_complex(radial_result.residual) ./ problem.coefficient_scale ./ residual_weights
    @assert norm(mismatch) <= 1e-5
    save_verification(output_prefix * "-seed.jls", (; problem, cloud_amplitudes,
        residual_weights, solution, profile, radial_result, weighted_residual_norm=norm(mismatch)))
    println("Saved matched seed, reintegrated weighted residual=", norm(mismatch)); flush(stdout)
end

function run_verification(arguments)
    length(arguments) in (3, 4) || error("usage: initial|seed|matched INPUT_JLS OUTPUT_PREFIX [SEED_JLS]")
    stage, input_path, output_prefix = arguments[1:3]
    stage in ("initial", "seed", "matched") || error("unknown verification stage")
    (stage == "matched") == (length(arguments) == 4) || error("matched stage requires SEED_JLS")
    any(ispath(output_prefix * suffix) for suffix in
        ("-setup.jls", "-ad.jls", "-fd.jls", "-accepted.jls", "-seed.jls", ".lock")) &&
        error("Output prefix already exists; preserve it and choose a fresh prefix")
    setup = verification_setup(deserialize(input_path))
    (; problem, cloud_amplitudes, secondary_fraction, residual_weights, matching_residual) = setup
    state = setup.initial_state
    seed_path = stage == "matched" ? abspath(arguments[4]) : nothing
    if seed_path !== nothing
        seed = deserialize(seed_path)
        @assert seed.cloud_amplitudes == cloud_amplitudes && seed.residual_weights == residual_weights
        @assert seed.problem.integration == problem.integration
        @assert seed.problem.source_basis.mode_labels == problem.source_basis.mode_labels
        state = seed.solution.state
    end
    @assert two_phase_profile_admissible(state, problem)
    mkpath(dirname(output_prefix)); mkdir(output_prefix * ".lock")
    try
        source_root = normpath(joinpath(@__DIR__, "../../.."))
        provenance_path = joinpath(source_root, "provenance.json")
        provenance = isfile(provenance_path) ? read(provenance_path, String) : nothing
        save_verification(output_prefix * "-setup.jls", (; stage, problem, cloud_amplitudes,
            secondary_fraction, residual_weights, state, input_path=abspath(input_path),
            input_sha256=bytes2hex(sha256(read(input_path))), seed_path,
            seed_sha256=seed_path === nothing ? nothing : bytes2hex(sha256(read(seed_path))), provenance))
        @show stage cloud_amplitudes problem.integration length(state)
        if stage == "seed"
            run_weak_seed(setup, output_prefix)
            return
        end
        mismatch_norm = norm(matching_residual(state))
        @show mismatch_norm
        stage == "matched" && @assert mismatch_norm <= 1e-5
        comparison = compare_ad_jacobians(matching_residual, state;
            output_path=output_prefix * "-ad.jls")
        println("AD all columns passed=", comparison.passed,
            " maximum relative error=", maximum(comparison.relative_errors)); flush(stdout)
        finite_differences = check_finite_differences(matching_residual, state,
            comparison.optimized_jacobian, state -> two_phase_profile_admissible(state, problem);
            output_path=output_prefix * "-fd.jls")
        println("Finite differences passed=", finite_differences.passed); flush(stdout)
        comparison.passed && finite_differences.passed || error("Verification incomplete; inspect saved comparisons")
    finally
        rm(output_prefix * ".lock")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_verification(ARGS)
end
