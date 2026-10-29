# Solve the existing multimode primary/secondary formalism at the sweep energies.
# Usage: julia --project=. omiya_energy_references.jl INPUT_DIRECTORY OUTPUT_DIRECTORY
using Serialization, LinearAlgebra, Printf, SHA, TOML, Dates
include(joinpath(@__DIR__, "../../../../src/adiabatic/primary_secondary/load.jl"))

function reconstruct_primary(setup, solution)
    p = setup.primary
    nonlinear_matching_residual(solution.frequency, solution.horizon_amplitudes,
        solution.outer_amplitudes, p.spin, p.alpha, p.modes, setup.couplings,
        setup.polar, setup.azimuthal, setup.grid; potential=Symbol(p.potential),
        p.options..., save_profile=true)
end

function primary_energy(setup, solution, profiles)
    p = setup.primary
    sum(integrate_energy_branch(branch, solution.frequency, p.spin, p.alpha,
        p.modes, setup.polar, setup.azimuthal, setup.grid; potential=Symbol(p.potential))
        for branch in (profiles.from_horizon, profiles.from_outer))
end

function primary_target(setup, target, prefix)
    p = setup.primary
    energies = [row.energy for row in p.rows]
    all(diff(energies) .> 0) || error("Require a monotone seed branch")
    first(energies) <= target <= last(energies) || error("Target outside seed branch")
    nearest = argmin(abs.(energies .- target))
    if abs(energies[nearest] / target - 1) < 1e-8
        entry = p.continuation[nearest]
        profiles = reconstruct_primary(setup, entry.solution)
        energy = primary_energy(setup, entry.solution, profiles)
        @assert norm(profiles.residual) / entry.amplitude < 1e-8
        @assert abs(energy / target - 1) < 1e-6
        result = (; entry.amplitude, entry.solution, profiles, energy, seed_indices=[nearest])
        serialize(prefix * "-primary.jls", result)
        return result
    end
    hi = searchsortedfirst(energies, target)
    lo = hi - 1
    left = merge(p.continuation[lo], (; energy=energies[lo]))
    right = merge(p.continuation[hi], (; energy=energies[hi]))
    for iteration in 1:8
        fraction = log(target / left.energy) / log(right.energy / left.energy)
        amplitude = exp((1-fraction)*log(left.amplitude) + fraction*log(right.amplitude))
        # Interpolation supplies a Newton guess; every new primary is solved.
        mix(key) = (1-fraction)*getproperty(left.solution, key) + fraction*getproperty(right.solution, key)
        guess_amplitude = (1-fraction)*left.amplitude + fraction*right.amplitude
        scale = amplitude / guess_amplitude
        progress = (i, residual) -> begin
            println("PRIMARY target=", target, " step=", iteration, " Newton=", i,
                    " residual=", residual); flush(stdout)
        end
        solution = solve_nonlinear_mode(mix(:frequency), amplitude, p.spin, p.alpha,
            p.modes, setup.couplings, setup.polar, setup.azimuthal, setup.grid;
            horizon_guess=scale .* mix(:horizon_amplitudes),
            outer_tail_guess=scale .* mix(:outer_amplitudes)[2:end],
            potential=Symbol(p.potential), p.options..., tolerance=1e-9,
            jacobian_backend=:forwarddiff, max_iterations=15, progress_callback=progress)
        profiles = reconstruct_primary(setup, solution)
        energy = primary_energy(setup, solution, profiles)
        @assert all(isfinite, profiles.from_horizon.states) && all(isfinite, profiles.from_outer.states)
        @assert norm(profiles.residual) / amplitude < 1e-8
        current = (; amplitude, solution, profiles, energy, seed_indices=[lo, hi])
        serialize(prefix * "-primary-step$iteration.jls", current)
        println("PRIMARY ENERGY target=", target, " measured=", energy,
                " relative_error=", energy/target-1); flush(stdout)
        if abs(energy / target - 1) < 1e-6
            serialize(prefix * "-primary.jls", current)
            return current
        end
        energy < target ? (left=current) : (right=current)
    end
    error("Primary energy targeting exhausted eight solves")
end

function secondary_target(setup, saved_secondary, primary, prefix)
    p = setup.primary
    rows = saved_secondary.rows
    amplitudes = [row.amplitude for row in rows]
    hi = clamp(searchsortedfirst(amplitudes, primary.amplitude), 2, length(rows))
    lo = hi - 1
    fraction = (primary.amplitude - amplitudes[lo]) / (amplitudes[hi] - amplitudes[lo])
    mix(key) = (1-fraction)*getproperty(rows[lo].solution, key) + fraction*getproperty(rows[hi].solution, key)
    frequency, horizon, outer = mix(:frequency), mix(:horizon_amplitudes), mix(:outer_amplitudes)
    residual_at(w) = secondary_matching_residual(w, horizon, outer, primary.profiles,
        p.modes, p.spin, p.alpha, setup.secondary_modes, setup.couplings,
        setup.polar, setup.azimuthal, setup.grid; reading=:averaged)
    # Use the existing branch's frequency-ULP floor rule, recording the tolerance.
    floor = norm(residual_at(frequency + eps(real(frequency))) - residual_at(frequency))
    tolerance = max(1e-8, 2floor)
    println("SECONDARY primary_energy=", primary.energy, " tolerance=", tolerance); flush(stdout)
    solution = solve_secondary_mode(frequency, primary.profiles, p.modes, p.spin,
        p.alpha, setup.secondary_modes, setup.couplings, setup.polar, setup.azimuthal,
        setup.grid; horizon_guess=horizon, outer_tail_guess=outer[2:end],
        tolerance, jacobian_backend=:forwarddiff, reading=:averaged)
    profiles = secondary_matching_residual(solution.frequency, solution.horizon_amplitudes,
        solution.outer_amplitudes, primary.profiles, p.modes, p.spin, p.alpha,
        setup.secondary_modes, setup.couplings, setup.polar, setup.azimuthal,
        setup.grid; reading=:averaged, save_profile=true)
    @assert norm(profiles.residual) <= tolerance
    @assert all(isfinite, profiles.from_horizon.states) && all(isfinite, profiles.from_outer.states)
    energy = measured_secondary_energy(setup, primary.profiles, profiles, solution.frequency)
    @assert isfinite(energy) && energy > 0
    result = (; solution, profiles, energy, tolerance, frequency_ulp_residual=floor)
    serialize(prefix * "-secondary-unit.jls", result)
    return result
end

function measured_secondary_energy(setup, primary, secondary, frequency)
    p = setup.primary
    energy = 0.0
    for side in (:from_horizon, :from_outer)
        branch, background = getproperty(secondary, side), getproperty(primary, side)
        positions = branch.positions
        order = issorted(positions) ? eachindex(positions) : reverse(eachindex(positions))
        densities = [secondary_energy_per_tortoise_radius(view(branch.states, :, i), positions[i],
            primary_values_at(background, positions[i], p.modes), frequency, p.spin, p.alpha,
            setup.secondary_modes, setup.polar, setup.azimuthal, setup.grid) for i in order]
        energy += integrate_uniform(densities, abs(positions[2]-positions[1]))
    end
    energy
end

function fundamental_sampler(profiles)
    horizon = profiles.from_horizon
    outer = (; positions=reverse(profiles.from_outer.positions),
             states=reverse(profiles.from_outer.states; dims=2))
    position -> begin
        branch = position <= last(horizon.positions) ? horizon : outer
        x, states = branch.positions, branch.states
        first(x) <= position <= last(x) || error("Sampling outside profile domain")
        i = clamp(searchsortedlast(x, position), 1, length(x)-1)
        h = x[i+1]-x[i]
        t = (position-x[i])/h
        n = size(states, 1) ÷ 2
        r0, r1, q0, q1 = states[1,i], states[1,i+1], states[n+1,i], states[n+1,i+1]
        value = (2t^3-3t^2+1)*r0 + (t^3-2t^2+t)*h*q0 + (-2t^3+3t^2)*r1 + (t^3-t^2)*h*q1
        slope = (6t^2-6t)/h*r0 + (3t^2-4t+1)*q0 + (-6t^2+6t)/h*r1 + (3t^2-2t)*q1
        (; value, slope)
    end
end

function export_reference(setup, primary, secondary, case, prefix, provenance)
    p = setup.primary
    targets = case["target_energies"]
    factor = sqrt(targets[2] / secondary.energy)
    scaled = (; (side => (; positions=getproperty(secondary.profiles, side).positions,
        states=factor .* getproperty(secondary.profiles, side).states)
        for side in (:from_horizon, :from_outer))...,
        residual=factor .* secondary.profiles.residual)
    measured = [primary.energy, measured_secondary_energy(setup, primary.profiles, scaled,
        secondary.solution.frequency)]
    @assert maximum(abs.(measured ./ targets .- 1)) < 1e-6
    metadata = Dict("passed"=>true, "model"=>"omiya_multimode_primary_secondary",
        "normalization"=>"independent_cloud_energies", "alpha"=>p.alpha, "spin"=>p.spin,
        "inner"=>p.options.inner, "outer"=>p.options.outer, "match"=>p.options.match,
        "primary_step"=>p.options.step, "secondary_step"=>2p.options.step,
        "primary_modes"=>string.(p.modes), "secondary_modes"=>string.(setup.secondary_modes),
        "target_energies"=>targets, "measured_energies"=>measured,
        "relative_energy_errors"=>measured ./ targets .- 1,
        "cloud_amplitudes"=>[primary.amplitude, factor],
        "primary_matching_residual"=>norm(primary.profiles.residual)/primary.amplitude,
        "secondary_matching_residual"=>norm(secondary.profiles.residual),
        "secondary_matching_tolerance"=>secondary.tolerance,
        "primary_energy_rescaled"=>false, "secondary_energy_scale"=>factor,
        "energy_convention"=>"isolated cosine primary; quadratic secondary on averaged primary",
        "frequencies_real"=>real.([primary.solution.frequency,secondary.solution.frequency]),
        "frequencies_imaginary"=>imag.([primary.solution.frequency,secondary.solution.frequency]),
        "spatial_convergence_established"=>false, "provenance"=>provenance)
    serialize(prefix * ".jls", (; primary, secondary, scaled_secondary=scaled, metadata))
    samplers = (fundamental_sampler(primary.profiles), fundamental_sampler(scaled))
    open(prefix * ".csv", "w") do stream
        println(stream, "curve,ell,m,r_star_over_M,R_real,R_imag,R_abs,Q_real,Q_imag")
        for radius in p.options.inner:0.25:p.options.outer, i in 1:2
            state = samplers[i](radius)
            name = i == 1 ? "omiya_primary" : "omiya_secondary"
            @printf(stream, "%s,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                name, i, i, radius, real(state.value), imag(state.value), abs(state.value),
                real(state.slope), imag(state.slope))
        end
    end
    open(prefix * ".toml", "w") do stream
        TOML.print(stream, metadata; sorted=true)
    end
    println("ACCEPTED ", case["name"], " energies=", measured, " targets=", targets); flush(stdout)
end

function main()
    length(ARGS) == 2 || error("usage: INPUT_DIRECTORY OUTPUT_DIRECTORY")
    inputs, output = abspath.(ARGS)
    ispath(output) && error("Preserve existing output; choose a fresh directory")
    config = TOML.parsefile(joinpath(inputs, "cases.toml"))
    primary_path, secondary_path = joinpath(inputs, "primary.jls"), joinpath(inputs, "secondary.jls")
    saved_secondary = deserialize(secondary_path)
    setup = load_primary_secondary_setup(primary_path; secondary_modes=saved_secondary.secondary_modes)
    p = setup.primary
    @assert Symbol(p.potential) == :cosine && saved_secondary.reading == "averaged"
    @assert p.spin == saved_secondary.spin && p.alpha == saved_secondary.alpha
    @assert p.options == saved_secondary.primary_options
    for case in config["cases"]
        @assert case["alpha"] == p.alpha && case["spin"] == p.spin
        @assert case["inner"] == p.options.inner && case["outer"] == p.options.outer
    end
    provenance = Dict("created_utc"=>string(now(UTC)), "slurm_job_id"=>get(ENV,"SLURM_JOB_ID","local"),
        "source_snapshot"=>normpath(joinpath(@__DIR__,"../../../..")),
        "input_sha256"=>Dict(basename(path)=>bytes2hex(sha256(read(path))) for path in
            (primary_path, secondary_path, joinpath(inputs,"cases.toml"))))
    mkpath(output)
    for (index, target) in enumerate(unique(case["target_energies"][1] for case in config["cases"]))
        prefix = joinpath(output, "background$index")
        primary = primary_target(setup, target, prefix)
        secondary = secondary_target(setup, saved_secondary, primary, prefix)
        for case in config["cases"]
            case["target_energies"][1] == target || continue
            export_reference(setup, primary, secondary, case, joinpath(output, case["name"]), provenance)
        end
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
