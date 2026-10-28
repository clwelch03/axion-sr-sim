using Serialization, LinearAlgebra
include(joinpath(@__DIR__, "../../../../src/adiabatic/multiphase/load.jl"))

# Only this diagnostic tag bypasses the optimized source derivative. It lets us
# differentiate the ordinary matrix arithmetic independently, one direction at
# a time, without changing the library or Newton's dispatch.
struct OrdinarySourceDirection end
function nonlinear_two_phase_source(
        parameters::AbstractVector{ForwardDiff.Dual{OrdinarySourceDirection,Float64,1}},
        basis::TwoPhaseSourceBasis; potential=:cosine)
    return invoke(nonlinear_two_phase_source, Tuple{Any,TwoPhaseSourceBasis},
        parameters, basis; potential)
end

function weak_secondary_setup(input; fraction=1e-3)
    0 < fraction <= 1 || error("secondary fraction must lie in (0, 1]")
    problem = input.problem
    amplitudes = input.cloud_amplitudes .* [1.0, fraction]
    labels = problem.source_basis.mode_labels
    # A component with second-phase index q2 begins at order A2^abs(q2).
    factors = [fraction^abs(q[2]) for (l, q) in labels]
    profile = TwoPhaseProfile(input.profile_guess.frequencies,
        factors .* input.profile_guess.horizon_coefficients,
        factors .* input.profile_guess.outer_coefficients)
    state = pack_two_phase_profile(profile, problem)
    # Keep the original coordinate scale. Tighten every secondary-dependent
    # residual by 1/fraction, so a weak secondary cannot be ignored by Newton.
    weights = repeat([q[2] == 0 ? 1.0 : fraction for (l, q) in labels], 4)
    residual = state -> two_phase_real_matching_residual(state, amplitudes, problem) ./ weights
    return (; problem, amplitudes, fraction, state, weights, residual)
end

function check_matching_jacobian(residual, state, admissible, output_path)
    println("Building the optimized matching Jacobian at the supplied state...")
    flush(stdout)
    jacobian = residual_jacobian(residual, state;
        backend=:forwarddiff, chunk_size=length(state))
    checks = []
    for pattern in (cos.(1:length(state)), sin.(1:length(state)))
        direction = max.(1.0, abs.(state)) .* (pattern ./ norm(pattern))
        predicted = jacobian * direction
        ordinary_state = [ForwardDiff.Dual{OrdinarySourceDirection}(
            state[i], (direction[i],)) for i in eachindex(state)]
        ordinary_residual = residual(ordinary_state)
        ordinary_derivative = [ForwardDiff.partials(value)[1] for value in ordinary_residual]
        automatic_error = norm(ordinary_derivative - predicted) / norm(predicted)
        println("Direction ", length(checks)+1, ": ordinary/optimized AD error = ", automatic_error)
        flush(stdout)
        differences = []
        for step in (1e-6, 1e-7, 1e-8)
            plus, minus = state + step * direction, state - step * direction
            admissible(plus) && admissible(minus) || error("inadmissible derivative probe")
            difference = (residual(plus) - residual(minus)) / (2step)
            relative_error = norm(difference - predicted) / norm(predicted)
            push!(differences, (; step, relative_error))
            println("  finite difference step=", step, ": relative error=", relative_error)
            flush(stdout)
        end
        push!(checks, (; direction, automatic_error, differences))
        # Preserve failed checks as well as successful ones.
        serialize(output_path, (; state, jacobian, checks))
        automatic_error <= 1e-8 || error("ordinary and optimized AD disagree")
        minimum(row.relative_error for row in differences) <= 1e-5 ||
            error("directional finite differences did not validate the Jacobian")
    end
    println("Matching Jacobian checks passed.")
    flush(stdout)
    return checks
end

function main(arguments)
    check_only = length(arguments) == 2 && arguments[1] == "--check"
    check_only || length(arguments) == 2 ||
        error("usage: weak_secondary_profile.jl [--check] INPUT_JLS [OUTPUT_PREFIX]")
    input_path = abspath(arguments[check_only ? 2 : 1])
    input = deserialize(input_path)
    setup = weak_secondary_setup(input)
    (; problem, amplitudes, fraction, state, residual, weights) = setup
    admissible = state -> two_phase_profile_admissible(state, problem)
    @assert admissible(state)
    unpacked = unpack_two_phase_profile(state, amplitudes, problem)
    fixed = two_phase_coefficient_indices(problem.components).fixed_outer_coefficient_indices
    @assert unpacked.outer_coefficients[fixed] == amplitudes
    @show fraction amplitudes problem.integration length(state)
    println("Residual tolerance 1e-5; secondary-dependent rows use the smaller amplitude.")
    flush(stdout)
    check_only && return

    output_prefix = abspath(arguments[2])
    any(ispath(output_prefix * suffix) for suffix in (".jls", "-jacobian.jls", ".lock")) &&
        error("Output prefix already exists")
    mkpath(dirname(output_prefix))
    mkdir(output_prefix * ".lock")
    try
        started = time()
        progress = (iteration, mismatch) -> begin
            println("[", round(time()-started; digits=1), " s] Newton ", iteration,
                ": weighted residual=", mismatch)
            flush(stdout)
        end
        solution = newton_solve(residual, state; tolerance=1e-5, max_iterations=12,
            admissible, progress_callback=progress, jacobian_backend=:forwarddiff,
            jacobian_chunk_size=length(state))
        profile = unpack_two_phase_profile(solution.state, amplitudes, problem)
        radial_result = two_phase_matching_residual(profile, problem; save_radial_result=true)
        recomputed = pack_complex(radial_result.residual) ./ problem.coefficient_scale ./ weights
        @assert norm(recomputed) <= 1e-5
        provenance = read(joinpath(@__DIR__, "../../../../provenance.json"), String)
        serialize(output_prefix * ".jls", (; problem, amplitudes, fraction, profile,
            solution, radial_result, weighted_residual_norm=norm(recomputed),
            target_amplitudes=input.cloud_amplitudes, provenance))
        println("Saved converged weak-secondary seed: ", output_prefix, ".jls")
        flush(stdout)
        check_matching_jacobian(residual, solution.state, admissible,
            output_prefix * "-jacobian.jls")
    finally
        rm(output_prefix * ".lock")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
