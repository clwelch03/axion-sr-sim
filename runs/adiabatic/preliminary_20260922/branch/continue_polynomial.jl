# Preliminary run recipe. Run from repository root with julia --project=.
using Serialization, LinearAlgebra
root = normpath(joinpath(@__DIR__, "../../../.."))
for filename in ("angular", "nonlinear", "tortoise", "radial", "boundary_series", "boundaries", "integration", "matching", "continuation_predictor", "energy")
    include(joinpath(root, "src/adiabatic", filename * ".jl"))
end
const potential_kind = Symbol(get(ENV, "AXION_POTENTIAL", "quartic"))
potential_kind in (:quartic, :sextic) ||
    throw(ArgumentError("AXION_POTENTIAL must be quartic or sextic"))
const output = joinpath(@__DIR__, "$(potential_kind)_branch")
const restarting = isfile(output * ".jls")
const seed_path = get(ENV, "AXION_SEED_PATH", restarting ? output * ".jls" :
    joinpath(root, "runs/adiabatic/preliminary_cosine_20260922.jls"))
const saved = deserialize(seed_path)
if restarting
    Symbol(saved.potential) == potential_kind ||
        error("Checkpoint potential does not match the requested potential")
end
const grid = angular_grid(saved.n_polar, saved.n_azimuth)
const polar = build_polar_values(5, grid.cosine_theta_nodes)
const azimuthal = build_azimuthal_values(5, grid.azimuth_nodes)
const couplings = build_cosine_squared_couplings(5, polar, grid)
const options = saved.options
const modes = saved.modes
const spin, alpha = saved.spin, saved.alpha
# Select the derivative method; difference_step applies only to finite differences.
const jacobian_backend = Symbol(get(ENV, "AXION_JACOBIAN", "finite"))
const difference_step = parse(Float64, get(ENV, "AXION_DIFFERENCE_STEP", "1e-8"))
isfinite(difference_step) && difference_step > 0 ||
    throw(ArgumentError("AXION_DIFFERENCE_STEP must be positive and finite"))
const newton_max_iterations = parse(Int, get(ENV, "AXION_NEWTON_MAX_ITERATIONS", "12"))
newton_max_iterations > 0 || throw(ArgumentError("require a positive correction limit"))
const predictor_kind = Symbol(get(ENV, "AXION_PREDICTOR", "secant"))
predictor_kind in (:secant, :cubic) || throw(ArgumentError("predictor must be secant or cubic"))
const history_spacing = parse(Float64, get(ENV, "AXION_HISTORY_SPACING", "0.05"))
isfinite(history_spacing) && history_spacing > 0 ||
    throw(ArgumentError("require positive finite history spacing"))

# Finite-radius outward current, with positive-m and conjugate normalization.
function outward_flux(state, position, frequency)
    count = length(modes)
    geometry = tortoise_geometry(position, spin)
    return 4π * geometry.radial_factor * sum(
        real(im * m * frequency * state[index] * conj(state[count+index]))
        for (index, (l, m)) in enumerate(modes))
end

function diagnostics(amplitude, solution)
    profiles = nonlinear_matching_residual(solution.frequency,
        solution.horizon_amplitudes, solution.outer_amplitudes,
        spin, alpha, modes, couplings, polar, azimuthal, grid;
        potential=potential_kind, save_profile=true, options...)
    energy = sum(integrate_energy_branch(branch, solution.frequency,
        spin, alpha, modes, polar, azimuthal, grid; potential=potential_kind)
        for branch in (profiles.from_horizon, profiles.from_outer))
    horizon_flux = -outward_flux(profiles.from_horizon.states[:,1], options.inner, solution.frequency)
    outer_flux = outward_flux(profiles.from_outer.states[:,1], options.outer, solution.frequency)
    peak_field = 0.0
    for branch in (profiles.from_horizon, profiles.from_outer), index in 1:10:length(branch.positions)
        values = Dict(mode=>branch.states[k,index] for (k,mode) in enumerate(modes))
        peak_field = max(peak_field, maximum(abs, reconstruct_field(values, polar, azimuthal, grid)))
    end
    return (; amplitude, energy, normalized_energy=energy/amplitude^2,
        omega_real=real(solution.frequency), omega_imag=imag(solution.frequency),
        horizon_flux, outer_flux, total_flux=horizon_flux+outer_flux,
        balance_ratio=-(horizon_flux+outer_flux)/(2imag(solution.frequency)*energy),
        peak_field, residual=solution.scaled_residual_norm, iterations=solution.iterations)
end

function checkpoint(continuation, rows)
    serialize(output * ".jls.tmp", (; spin, alpha, modes, options, n_polar=saved.n_polar,
        n_azimuth=saved.n_azimuth, potential=String(potential_kind), tolerance=1e-9,
        predictor=String(predictor_kind), history_spacing,
        jacobian_backend=String(jacobian_backend),
        newton_max_iterations, difference_step, continuation, rows))
    mv(output * ".jls.tmp", output * ".jls"; force=true)
    open(output * ".csv.tmp", "w") do io
        println(io, join(string.(keys(first(rows))), ','))
        for row in rows
            println(io, join(values(row), ','))
        end
    end
    mv(output * ".csv.tmp", output * ".csv"; force=true)
end

function main()
    started = time()
    runtime_seconds = parse(Float64, get(ENV, "AXION_RUN_SECONDS", "1000"))
    isfinite(runtime_seconds) && runtime_seconds > 0 ||
        throw(ArgumentError("AXION_RUN_SECONDS must be positive and finite"))
    target_amplitude = parse(Float64, get(ENV, "AXION_TARGET_AMPLITUDE", "0.1"))
    isfinite(target_amplitude) && target_amplitude > 0 ||
        throw(ArgumentError("AXION_TARGET_AMPLITUDE must be positive and finite"))
    println("Potential=", potential_kind, " restarting=", restarting, " output=", output)
    println("Jacobian=", jacobian_backend, " predictor=", predictor_kind)
    flush(stdout)
    if restarting
        continuation = copy(saved.continuation)
        rows = copy(saved.rows)
    else
        # Use the cosine seed only as a guess; solve and diagnose this potential.
        seed_entry = first(saved.continuation)
        seed_solution = seed_entry.solution
        solution = solve_nonlinear_mode(
            seed_solution.frequency, seed_entry.amplitude,
            spin, alpha, modes, couplings, polar, azimuthal, grid;
            horizon_guess=seed_solution.horizon_amplitudes,
            outer_tail_guess=seed_solution.outer_amplitudes[2:end],
            potential=potential_kind,
            tolerance=1e-9, max_iterations=newton_max_iterations,
            difference_step, jacobian_backend, options...,
        )
        continuation = [(amplitude=seed_entry.amplitude, solution)]
        rows = [diagnostics(seed_entry.amplitude, solution)]
        println("SEED ", first(rows))
        flush(stdout)
    end
    amplitude = last(continuation).amplitude
    # Respect the accepted spacing; do not force a large restart jump.
    amplitude_step = length(continuation) >= 2 ?
        min(1.5 * (amplitude - continuation[end-1].amplitude), 0.1amplitude) :
        0.1amplitude
    if haskey(ENV, "AXION_INITIAL_STEP")
        amplitude_step = parse(Float64, ENV["AXION_INITIAL_STEP"])
        isfinite(amplitude_step) && amplitude_step > 0 ||
            throw(ArgumentError("AXION_INITIAL_STEP must be positive and finite"))
    end
    checkpoint(continuation, rows)
    while time()-started < runtime_seconds && amplitude < target_amplitude
        target = min(target_amplitude, amplitude + amplitude_step)
        prediction = predictor_kind === :cubic ?
            predict_nonlinear_mode_cubic(continuation, target; history_spacing) :
            predict_nonlinear_mode(continuation, target)
        if !(isfinite(prediction.frequency) && 0 < real(prediction.frequency) < alpha &&
             all(isfinite, prediction.horizon_amplitudes) && all(isfinite, prediction.outer_amplitudes))
            amplitude_step /= 2
            amplitude_step < 1e-5amplitude && error("No admissible continuation prediction")
            continue
        end
        println("Starting A0=$target (elapsed $(round(time()-started)) s)")
        flush(stdout)
        progress = (iteration, residual_norm) -> begin
            println("  Newton iteration=$iteration scaled_residual=$residual_norm")
            flush(stdout)
        end
        solution = try
            solve_nonlinear_mode(prediction.frequency, target,
                spin, alpha, modes, couplings, polar, azimuthal, grid;
                horizon_guess=prediction.horizon_amplitudes,
                outer_tail_guess=prediction.outer_amplitudes[2:end],
                potential=potential_kind,
                tolerance=1e-9, max_iterations=newton_max_iterations,
                difference_step, jacobian_backend, progress_callback=progress, options...)
        catch exception
            exception isa ErrorException || rethrow()
            println("FAILED A0=$target: ", sprint(showerror, exception))
            flush(stdout)
            amplitude_step /= 2
            amplitude_step < 1e-5amplitude &&
                error("Continuation step exhausted; last accepted A0=$amplitude")
            continue
        end
        row = diagnostics(target, solution)
        push!(continuation, (; amplitude=target, solution))
        push!(rows, row)
        checkpoint(continuation, rows)
        println("SAVED ", row, " elapsed=", round(time()-started))
        flush(stdout)
        accepted_step = target - amplitude
        amplitude = target
        growth = solution.iterations <= 4 ? 1.6 : solution.iterations <= 8 ? 1.3 :
            solution.iterations <= 16 ? 1.0 : 0.75
        amplitude_step = min(growth * accepted_step, 0.25amplitude)
    end
    println("Finished at A0=", amplitude, " with ", length(rows), " saved points.")
end
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
