include("model.jl")
include("correction.jl")
include("spatial.jl")

const cases = [
    (label="shooting_fd_1e8", method=:shooting, derivative=:finite, increment=1e-8, strategy=:residual, resolution=0.0),
    (label="shooting_fd_1e10", method=:shooting, derivative=:finite, increment=1e-10, strategy=:residual, resolution=0.0),
    (label="shooting_ad", method=:shooting, derivative=:automatic, increment=0.0, strategy=:residual, resolution=0.0),
    (label="shooting_ad_deuflhard", method=:shooting, derivative=:automatic, increment=0.0, strategy=:deuflhard, resolution=0.0),
    (label="multiple_4", method=:multiple, derivative=:automatic, increment=0.0, strategy=:residual, resolution=2.0),
    (label="multiple_8", method=:multiple, derivative=:automatic, increment=0.0, strategy=:residual, resolution=4.0),
    (label="collocation_050", method=:collocation, derivative=:automatic, increment=0.0, strategy=:residual, resolution=0.5),
    (label="collocation_025", method=:collocation, derivative=:automatic, increment=0.0, strategy=:residual, resolution=0.25),
    (label="collocation_0125", method=:collocation, derivative=:automatic, increment=0.0, strategy=:residual, resolution=0.125),
]

"""
    evaluate_solution(solution, amplitude, context)

Independently replay the original production integrator and compute observables.

Report both matching residual and profile diagnostics. A converged collocation
defect is not the same quantity as a converged shooting mismatch.
"""
function evaluate_solution(solution, amplitude, context)
    profiles = nonlinear_matching_residual(solution.frequency, solution.horizon_amplitudes,
        solution.outer_amplitudes, context.spin, context.alpha, context.modes,
        context.couplings, context.polar, context.azimuthal, context.grid;
        save_profile=true, context.options...)
    energy = sum(integrate_energy_branch(branch, solution.frequency, context.spin,
        context.alpha, context.modes, context.polar, context.azimuthal, context.grid)
        for branch in (profiles.from_horizon, profiles.from_outer))
    flux = function(branch)
        state = branch.states[:,1]
        geometry = tortoise_geometry(first(branch.positions), context.spin)
        return 4π*geometry.radial_factor*sum(real(im*harmonic*solution.frequency*state[index]*
            conj(state[6+index])) for (index,(_,harmonic)) in enumerate(context.modes))
    end
    return (; replay_residual=norm(profiles.residual)/amplitude, energy,
        total_flux=flux(profiles.from_outer)-flux(profiles.from_horizon),
        omega_real=real(solution.frequency), omega_imag=imag(solution.frequency))
end

function run_case(case_index)
    configuration = cases[case_index]
    configuration = merge(configuration, (; label=configuration.label * get(ENV, "AXION_LABEL_SUFFIX", "")))
    tolerance = parse(Float64, get(ENV, "AXION_NATIVE_TOLERANCE", "1e-9"))
    isfinite(tolerance) && tolerance > 0 || error("Require a positive native tolerance")
    if haskey(ENV, "AXION_RADIAL_RESOLUTION")
        configuration = merge(configuration, (; resolution=parse(Float64, ENV["AXION_RADIAL_RESOLUTION"])))
    end
    root = normpath(joinpath(@__DIR__, "../../.."))
    directory = joinpath(root, "runs/adiabatic/radialsolvertest")
    saved = deserialize(joinpath(directory,"seed_cosine.jls"))
    context = benchmark_context(saved)
    if haskey(ENV, "AXION_RK_STEP")
        context = merge(context, (; options=merge(context.options, (; step=parse(Float64, ENV["AXION_RK_STEP"])))))
    end
    amplitude = last(saved.continuation).amplitude + 0.1
    prediction = predict_nonlinear_mode_cubic(saved.continuation, amplitude; history_spacing=0.025)
    started = time()
    println("CASE ", configuration, " A0=", amplitude)
    flush(stdout)
    shooting = shooting_problem(context, amplitude, prediction)
    parameter_count = length(shooting.initial)
    spatial = configuration.method === :shooting ? nothing :
        spatial_problem(shooting; method=configuration.method, resolution=configuration.resolution,
            align_mesh=get(ENV, "AXION_ALIGN_MESH", "false") == "true")
    initial = isnothing(spatial) ? shooting.initial : spatial.initial
    residual = isnothing(spatial) ? shooting.residual : spatial.residual
    jacobian = if !isnothing(spatial)
        spatial.jacobian
    elseif configuration.derivative === :finite
        state -> difference_jacobian(residual, state, configuration.increment)
    else
        state -> ForwardDiff.jacobian(residual, state)
    end
    # Warm residual and derivative specialization before measuring corrections.
    initial_residual = norm(residual(initial))
    first_jacobian = jacobian(initial)
    setup_seconds = time()-started
    println("SETUP seconds=", setup_seconds, " unknowns=", length(initial),
        " residual=", initial_residual, " jacobian_nnz=", count(!iszero, first_jacobian))
    flush(stdout)
    result_path = joinpath(directory, configuration.label * ".jls")
    callback = iterate -> serialize(result_path * ".iterate", iterate)
    admissible = state -> all(isfinite,state) && 0 < state[1] < context.alpha
    started = time()
    result = correct_system(residual, jacobian, initial; strategy=configuration.strategy,
        admissible, callback, tolerance, max_iterations=20)
    solve_seconds = time()-started
    solution = physical_solution(result.state[1:parameter_count], amplitude,
        length(result.history)-1, result.residual)
    diagnostics = evaluate_solution(solution, amplitude, context)
    native_energy = diagnostics.energy
    if configuration.method === :collocation
        states = reshape(result.state[parameter_count+1:end], parameter_count, :)
        complex_states = hcat([amplitude*unpack_complex(states[:,index]) for index in axes(states,2)]...)
        match_index = findfirst(==(context.options.match), spatial.positions)
        branches = ((positions=spatial.positions[1:match_index], states=complex_states[:,1:match_index]),
            (positions=spatial.positions[end:-1:match_index], states=complex_states[:,end:-1:match_index]))
        native_energy = sum(integrate_energy_branch(branch, solution.frequency, context.spin,
            context.alpha, context.modes, context.polar, context.azimuthal, context.grid) for branch in branches)
    end
    serialize(result_path, (; configuration, amplitude, tolerance, options=context.options, setup_seconds, solve_seconds,
        initial_residual, result, solution, diagnostics, native_energy,
        positions=isnothing(spatial) ? nothing : spatial.positions))
    summary = (; method=configuration.label, amplitude, tolerance, rk_step=context.options.step,
        resolution=configuration.resolution, status=result.status, setup_seconds,
        solve_seconds, iterations=length(result.history)-1, residual=result.residual,
        diagnostics..., native_energy)
    open(joinpath(directory, configuration.label * ".csv"), "w") do stream
        println(stream, join(string.(keys(summary)), ','))
        println(stream, join(values(summary), ','))
    end
    println("SUMMARY ", summary)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_case(parse(Int, get(ENV,"SLURM_ARRAY_TASK_ID","0"))+1)
end
