# Continue the cosine branch from a saved checkpoint.
using Serialization, LinearAlgebra

# Load the solver.
root = normpath(joinpath(@__DIR__, "../../../.."))
for filename in (
    "angular", "nonlinear", "tortoise", "radial", "boundary_series",
    "boundaries", "integration", "matching", "continuation_predictor", "energy",
)
    include(joinpath(root, "src/adiabatic", filename * ".jl"))
end

# Resume this branch unless an explicit seed is supplied.
const output = get(ENV, "AXION_OUTPUT_PREFIX", joinpath(@__DIR__, "../results/cosine_branch"))
const default_seed = isfile(output * ".jls") ? output * ".jls" :
    joinpath(root, "runs/adiabatic/fig2_20260922/seeds/preliminary_cosine_20260922.jls")
const seed_path = get(ENV, "AXION_SEED_PATH", default_seed)
const saved = deserialize(seed_path)

# Preserve the checkpoint's physics and integration settings.
const spin, alpha = saved.spin, saved.alpha
const modes = saved.modes
const options = saved.options

const grid = angular_grid(saved.n_polar, saved.n_azimuth)
const polar = build_polar_values(5, grid.cosine_theta_nodes)
const azimuthal = build_azimuthal_values(5, grid.azimuth_nodes)
const couplings = build_cosine_squared_couplings(5, polar, grid)

# Predictor and correction settings. difference_step applies only to finite differences.
const predictor_kind = Symbol(get(ENV, "AXION_PREDICTOR", "secant"))
const history_spacing = parse(Float64, get(ENV, "AXION_HISTORY_SPACING", "0.025"))
const jacobian_backend = Symbol(get(ENV, "AXION_JACOBIAN", "forwarddiff"))
const difference_step = parse(Float64, get(ENV, "AXION_DIFFERENCE_STEP", "1e-8"))
const newton_max_iterations = parse(Int, get(ENV, "AXION_NEWTON_MAX_ITERATIONS", "12"))

predictor_kind in (:secant, :cubic) ||
    throw(ArgumentError("predictor must be secant or cubic"))
isfinite(history_spacing) && history_spacing > 0 ||
    throw(ArgumentError("require positive finite history spacing"))
isfinite(difference_step) && difference_step > 0 ||
    throw(ArgumentError("difference_step must be positive and finite"))
newton_max_iterations > 0 ||
    throw(ArgumentError("require a positive correction limit"))


# Finite-radius outward current, including the conjugate-mode normalization.
function outward_flux(state, position, frequency)
    count = length(modes)
    geometry = tortoise_geometry(position, spin)

    return 4π * geometry.radial_factor * sum(
        real(im * m * frequency * state[index] * conj(state[count + index]))
        for (index, (_, m)) in enumerate(modes)
    )
end


function diagnostics(amplitude, solution)
    # Reintegrate the accepted solution, retaining profiles for the energy integral.
    profiles = nonlinear_matching_residual(
        solution.frequency, solution.horizon_amplitudes, solution.outer_amplitudes,
        spin, alpha, modes, couplings, polar, azimuthal, grid;
        save_profile=true, options...,
    )
    branches = (profiles.from_horizon, profiles.from_outer)

    energy = sum(
        integrate_energy_branch(
            branch, solution.frequency, spin, alpha, modes, polar, azimuthal, grid,
        )
        for branch in branches
    )

    # Loss into the horizon has the opposite sign to the outward current there.
    horizon_flux = -outward_flux(
        profiles.from_horizon.states[:, 1], options.inner, solution.frequency,
    )
    outer_flux = outward_flux(
        profiles.from_outer.states[:, 1], options.outer, solution.frequency,
    )
    total_flux = horizon_flux + outer_flux

    # Estimate the peak field using every tenth radial sample.
    peak_field = 0.0
    for branch in branches, index in 1:10:length(branch.positions)
        radial_values = Dict(
            mode => branch.states[k, index] for (k, mode) in enumerate(modes)
        )
        field = reconstruct_field(radial_values, polar, azimuthal, grid)
        peak_field = max(peak_field, maximum(abs, field))
    end

    return (;
        amplitude, energy, normalized_energy=energy / amplitude^2,
        omega_real=real(solution.frequency), omega_imag=imag(solution.frequency),
        horizon_flux, outer_flux, total_flux,
        balance_ratio=-total_flux / (2 * imag(solution.frequency) * energy),
        peak_field, residual=solution.scaled_residual_norm,
        iterations=solution.iterations,
    )
end


function checkpoint(continuation, rows)
    # Save full solutions for restarting, plus the settings used to obtain them.
    data = (;
        spin, alpha, modes, options,
        n_polar=saved.n_polar, n_azimuth=saved.n_azimuth,
        potential="cosine", tolerance=1e-9,
        predictor=String(predictor_kind), history_spacing,
        jacobian_backend=String(jacobian_backend),
        newton_max_iterations, difference_step, continuation, rows,
    )

    # Replace each previous output only after its temporary file is complete.
    serialize(output * ".jls.tmp", data)
    mv(output * ".jls.tmp", output * ".jls"; force=true)

    open(output * ".csv.tmp", "w") do stream
        println(stream, join(string.(keys(first(rows))), ','))
        for row in rows
            println(stream, join(values(row), ','))
        end
    end
    mv(output * ".csv.tmp", output * ".csv"; force=true)
end


function main()
    started = time()
    runtime_seconds = parse(Float64, get(ENV, "AXION_RUN_SECONDS", "1000"))
    target_amplitude = parse(Float64, get(ENV, "AXION_TARGET_AMPLITUDE", "13"))

    isfinite(runtime_seconds) && runtime_seconds > 0 ||
        throw(ArgumentError("runtime must be positive and finite"))
    isfinite(target_amplitude) && target_amplitude > 0 ||
        throw(ArgumentError("target amplitude must be positive and finite"))

    continuation = copy(saved.continuation)
    rows = copy(saved.rows)
    amplitude = last(continuation).amplitude

    # Respect small accepted steps when restarting a difficult branch.
    amplitude_step = 0.1 * amplitude
    if length(continuation) >= 2
        previous_step = amplitude - continuation[end - 1].amplitude
        amplitude_step = min(1.5 * previous_step, amplitude_step)
    end

    if haskey(ENV, "AXION_INITIAL_STEP")
        amplitude_step = parse(Float64, ENV["AXION_INITIAL_STEP"])
        isfinite(amplitude_step) && amplitude_step > 0 ||
            throw(ArgumentError("require positive finite initial step"))
    end

    checkpoint(continuation, rows)

    while time() - started < runtime_seconds && amplitude < target_amplitude
        # Predict the next point without advancing beyond the requested target.
        target = min(target_amplitude, amplitude + amplitude_step)
        prediction = if predictor_kind === :cubic
            predict_nonlinear_mode_cubic(continuation, target; history_spacing)
        else
            predict_nonlinear_mode(continuation, target)
        end

        valid_frequency = isfinite(prediction.frequency) &&
            0 < real(prediction.frequency) < alpha
        valid_amplitudes = all(isfinite, prediction.horizon_amplitudes) &&
            all(isfinite, prediction.outer_amplitudes)

        if !valid_frequency || !valid_amplitudes
            amplitude_step /= 2
            amplitude_step < 1e-5 * amplitude &&
                error("No admissible continuation prediction")
            continue
        end

        println("Starting A0=$target (elapsed $(round(time() - started)) s)")
        flush(stdout)

        progress = function(iteration, residual_norm)
            println("  Newton iteration=$iteration scaled_residual=$residual_norm")
            flush(stdout)
        end

        # Correct the guess. On failure, retry from the last accepted point.
        solution = try
            solve_nonlinear_mode(
                prediction.frequency, target,
                spin, alpha, modes, couplings, polar, azimuthal, grid;
                horizon_guess=prediction.horizon_amplitudes,
                outer_tail_guess=prediction.outer_amplitudes[2:end],
                tolerance=1e-9, max_iterations=newton_max_iterations,
                difference_step, jacobian_backend,
                progress_callback=progress, options...,
            )
        catch exception
            exception isa ErrorException || rethrow()
            println("FAILED A0=$target: ", sprint(showerror, exception))
            flush(stdout)

            amplitude_step /= 2
            amplitude_step < 1e-5 * amplitude &&
                error("Continuation step exhausted; last accepted A0=$amplitude")
            continue
        end

        # Measure and save the accepted point before attempting another.
        row = diagnostics(target, solution)
        push!(continuation, (; amplitude=target, solution))
        push!(rows, row)
        checkpoint(continuation, rows)

        println("SAVED ", row, " elapsed=", round(time() - started))
        flush(stdout)

        # Cheap corrections allow a larger next step; expensive ones reduce it.
        accepted_step = target - amplitude
        amplitude = target

        growth = if solution.iterations <= 4
            1.6
        elseif solution.iterations <= 8
            1.3
        elseif solution.iterations <= 16
            1.0
        else
            0.75
        end

        amplitude_step = min(growth * accepted_step, 0.25 * amplitude)
    end

    println("Finished at A0=", amplitude, " with ", length(rows), " saved points.")
end


# Allow diagnostic scripts to include this file without starting continuation.
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end