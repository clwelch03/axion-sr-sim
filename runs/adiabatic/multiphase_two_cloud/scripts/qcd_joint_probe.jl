# Joint cosine equations in the physical, microscopic-second-cloud limit.
# Keep the nonlinear dominant cloud and the coupled minor/(2,-1) response.
# Unlike the Omiya reference, the driven component feeds back into the minor mode.
# Usage: julia --project=. qcd_joint_probe.jl INPUT_DIRECTORY OUTPUT_PREFIX
using Serialization, LinearAlgebra, Printf, SHA, TOML, Dates
include(joinpath(@__DIR__, "../../../../src/adiabatic/multiphase/load.jl"))
include(joinpath(@__DIR__, "qcd_secondary_preview.jl"))

function joint_probe_main()
    length(ARGS) == 2 || error("usage: INPUT_DIRECTORY OUTPUT_PREFIX")
    inputs, prefix = abspath.(ARGS)
    any(ispath(prefix * suffix) for suffix in (".jls", ".toml", ".csv", "-primary.jls")) &&
        error("Preserve existing output; choose a fresh prefix")
    mkpath(dirname(prefix))
    primary = deserialize(joinpath(inputs, "primary.jls"))
    preview = deserialize(joinpath(inputs, "preview.jls"))
    @assert primary.report["passed"] && preview.report["passed"]
    case = primary.report["case"]
    spin, alpha = case["spin"], case["alpha"]
    dominant_m = case["dominant_m"]
    minor_m = 3-dominant_m
    dominant_mode = case["dominant_mode"]
    minor_mode = dominant_mode == "211" ? "322" : "211"
    minor_modes = [(minor_m+2j,minor_m) for j in 0:2]
    driven_m = 2dominant_m-minor_m
    maximum_degree = max(maximum(first,primary.modes),minor_m+4,driven_m+4)
    grid = angular_grid(case["polar_nodes"], case["azimuthal_nodes"])
    polar = build_polar_values(maximum_degree, grid.cosine_theta_nodes)
    azimuthal = build_azimuthal_values(maximum_degree, grid.azimuth_nodes)
    couplings = build_cosine_squared_couplings(maximum_degree, polar, grid)
    cloud_m = [dominant_m, minor_m]
    primary_components = [fourier_component((m÷dominant_m, 0), [l for (l, mm) in primary.modes if mm == m], cloud_m)
                          for m in unique(last.(primary.modes))]
    frequencies = [primary.solution.frequency, preview.secondary_solution.frequency]
    integration = (; primary.options..., horizon_order=3, outer_order=7)
    primary_problem = TwoPhaseProblem(; spin, alpha, components=primary_components,
        cloud_m, couplings, polar_values=polar, grid, phase1_nodes=grid.azimuth_nodes,
        phase2_nodes=[0.0], potential=:cosine, coefficient_scale=primary.amplitude,
        free_outer_coefficient_indices=collect(2:length(primary.modes)),
        initial_bound=BitVector([real(sum(c.q .* frequencies)) < alpha for c in primary_components]),
        integration)
    primary_profile = TwoPhaseProfile(frequencies, primary.solution.horizon_amplitudes,
        primary.solution.outer_amplitudes)
    println("Integrating nonlinear primary with independent joint RHS"); flush(stdout)
    joint_primary = two_phase_matching_residual(primary_profile, primary_problem; save_radial_result=true)
    primary_matching = norm(joint_primary.residual) / primary.amplitude
    @assert primary_matching <= primary.report["matching_tolerance"]
    primary_difference = maximum(maximum(abs, getproperty(joint_primary, side).states -
        getproperty(primary.profiles, side).states) / maximum(abs, getproperty(primary.profiles, side).states)
        for side in (:from_horizon, :from_outer))
    joint_primary_energy = sum(integrate_energy_branch(branch, frequencies[1], spin, alpha,
        primary.modes, polar, azimuthal, grid; potential=:cosine,fundamental_m=dominant_m)
        for branch in (joint_primary.from_horizon, joint_primary.from_outer))
    @assert isapprox(joint_primary_energy, preview.report["energy_"*dominant_mode]; rtol=1e-10)
    serialize(prefix * "-primary.jls", (; primary_problem, primary_profile, joint_primary,
        primary_matching, primary_difference, joint_primary_energy))

    # In the retained joint tower, these are all components linear in the
    # small cloud. The (-1,2) mixed component starts at its amplitude squared;
    # its feedback on either plotted fundamental is higher order still.
    probe_components = [fourier_component((0, 1),first.(minor_modes),cloud_m),
                        fourier_component((2, -1),collect(driven_m:2:driven_m+4),cloud_m)]
    components = vcat(primary_components, probe_components)
    labels = [(l, c.q) for c in components for l in c.l_values]
    primary_count = length(primary.modes)
    probe_count = sum(length(c.l_values) for c in probe_components)
    count = primary_count + probe_count
    probe_indices = vcat(collect(primary_count+1:count), collect(count+primary_count+1:2count))
    # Three phase-2 nodes integrate this derivative exactly: both response
    # factors have |q2|=1, and the background has q2=0.
    basis = build_two_phase_source_basis(labels, cloud_m, polar, grid,
        grid.azimuth_nodes, 2pi .* (0:2) ./ 3)
    response_at(values, position; source_basis=basis) = begin
        geometry = tortoise_geometry(position, spin)
        coefficients = vcat(values, zeros(ComplexF64, probe_count))
        parameters = vcat(pack_complex(coefficients), [geometry.radius, spin, alpha])
        _, derivative = local_source_jacobian(parameters, source_basis; potential=:cosine)
        derivative[probe_indices, probe_indices]
    end

    # Check the response against independent full-source differences, a denser
    # phase grid, the Omiya diagonal block, and the exact zero-background limit.
    check_position = 25.0
    background = joint_primary.from_outer
    check_values_dict = primary_values_at(background, check_position, primary.modes)
    check_values = [check_values_dict[mode] for mode in primary.modes]
    response = response_at(check_values, check_position)
    dense_basis = build_two_phase_source_basis(labels, cloud_m, polar, grid,
        grid.azimuth_nodes, 2pi .* (0:6) ./ 7)
    phase_error = norm(response-response_at(check_values, check_position; source_basis=dense_basis)) / norm(response)
    @assert phase_error < 1e-12
    @assert iszero(norm(response_at(zero.(check_values), check_position)))
    direction = ComplexF64[1+.2im, -.03+.01im, .01-.02im, .3-.1im, .02+.03im, -.01+.01im]
    geometry = tortoise_geometry(check_position, spin)
    diagonal_direction = vcat(direction[1:3], zeros(ComplexF64, 3))
    diagonal_response = unpack_complex(response * pack_complex(diagonal_direction))[1:3]
    omiya_source = secondary_background_source(check_values_dict,
        Dict(mode=>direction[i] for (i,mode) in enumerate(minor_modes)), geometry.radius,
        spin, alpha, minor_modes, polar, azimuthal, grid)
    expected_diagonal = [omiya_source[mode] for mode in minor_modes]
    diagonal_error = norm(diagonal_response-expected_diagonal) / norm(expected_diagonal)
    @assert diagonal_error < 1e-12
    base = vcat(check_values, zeros(ComplexF64, probe_count))
    perturbation = vcat(zeros(ComplexF64, primary_count), direction)
    full_source_at(epsilon) = nonlinear_two_phase_source(vcat(pack_complex(base + epsilon * perturbation),
        [geometry.radius, spin, alpha]), dense_basis; potential=:cosine)[probe_indices]
    expected = response * pack_complex(direction)
    check_amplitude_scale = maximum(abs,check_values) / 0.014
    finite_errors = [norm((full_source_at(h)-full_source_at(-h))/(2h)-expected)/norm(expected)
                     for h in check_amplitude_scale .* (1e-5, 1e-6)]
    @assert finite_errors[2] < 2e-7 && finite_errors[2] < finite_errors[1] / 20
    println("CHECKS phase=", phase_error, " Omiya block=", diagonal_error,
        " full-source difference errors=", finite_errors); flush(stdout)

    println("Tabulating the joint cosine response on both radial branches"); flush(stdout)
    response_tables = Dict{Symbol,Vector{Matrix{Float64}}}()
    for side in (:from_horizon, :from_outer)
        branch = getproperty(joint_primary, side)
        response_tables[side] = [response_at(view(branch.states, 1:primary_count, i), position)
            for (i,position) in enumerate(branch.positions)]
    end
    function integrate_probe(state; save_profile=false, feedback=true)
        secondary_frequency = state[1]
        horizon = state[2:probe_count+1]
        outer = vcat(one(secondary_frequency), state[probe_count+2:end])
        trial_frequencies = [frequencies[1], secondary_frequency]
        function branch_result(side)
            branch = getproperty(joint_primary, side)
            position = first(branch.positions)
            boundary_side = side == :from_horizon ? :horizon : :outer
            coefficients = side == :from_horizon ? horizon : outer
            blocks = [component_boundary_state(boundary_side, position, spin, alpha, component,
                trial_frequencies, cloud_m, couplings, coefficients[3k-2:3k];
                horizon_order=3, outer_order=7) for (k, component) in enumerate(probe_components)]
            boundary_state = vcat((block[1:3] for block in blocks)..., (block[4:6] for block in blocks)...)
            function rhs!(derivative, radial_state, position)
                sample = round(Int, (position-first(branch.positions))/(branch.positions[2]-branch.positions[1]))+1
                @assert abs(branch.positions[sample]-position) < 1e-9
                response = response_tables[side][sample]
                packed = pack_complex(radial_state[1:probe_count])
                source = unpack_complex(response * packed)
                if !feedback
                    # Independent check: removing only m=0 -> 322 feedback
                    # must reproduce the existing Omiya eigenmode.
                    without_driven = copy(packed)
                    without_driven[vcat(4:6, 10:12)] .= 0
                    source[1:3] .= unpack_complex(response * without_driven)[1:3]
                end
                geometry = tortoise_geometry(position, spin)
                factor = geometry.delta / geometry.radial_factor^2
                derivative[1:probe_count] .= radial_state[probe_count+1:end]
                for (k, component) in enumerate(probe_components)
                    indices = 3k-2:3k
                    operator = component_tortoise_operator(position, spin, alpha, component,
                        trial_frequencies, cloud_m, couplings)
                    derivative[probe_count .+ indices] .= -operator.slope_coefficient .* radial_state[probe_count .+ indices] -
                        operator.profile_matrix * radial_state[indices] - factor .* source[indices]
                end
                return derivative
            end
            step = 2abs(branch.positions[2]-branch.positions[1]) * (1+1e-9)
            return integrate_radial(rhs!, boundary_state, position, last(branch.positions); step, save_profile)
        end
        from_horizon, from_outer = branch_result(:from_horizon), branch_result(:from_outer)
        if save_profile
            branch_mismatch = from_horizon.states[:,end]-from_outer.states[:,end]
            return (; from_horizon, from_outer, residual=branch_mismatch)
        end
        return from_horizon-from_outer
    end
    seed = preview.secondary_solution
    initial = vcat(seed.frequency, seed.horizon_amplitudes, zeros(ComplexF64,3),
        seed.outer_amplitudes[2:end], zeros(ComplexF64,3))
    control = integrate_probe(initial; save_profile=true, feedback=false)
    omiya_reference = only(filter(reference -> reference.cloud == minor_mode && reference.model == "omiya",
        preview.references))
    secondary_rows = vcat(1:3, probe_count .+ (1:3))
    control_difference = maximum(begin
        expected_states = getproperty(omiya_reference.profiles,side).states ./
            preview.report["secondary_omiya_amplitude_factor"]
        maximum(abs, getproperty(control,side).states[secondary_rows,:]-expected_states) /
            maximum(abs,expected_states)
    end for side in (:from_horizon,:from_outer))
    println("Omiya recovery with return coupling disabled: ",control_difference); flush(stdout)
    # The longer alpha=.2 domain has a measured 2e-9 discrepancy between
    # equivalent radial implementations. Retain this measured floor and flag
    # joint/Omiya differences near it instead of interpreting them as physics.
    control_tolerance = 1e-8
    @assert control_difference < control_tolerance
    @assert initial[1] == seed.frequency
    @assert length(initial) == 2probe_count
    progress = (iteration, residual) -> begin
        println("JOINT PROBE NEWTON ",iteration," residual=",residual); flush(stdout)
    end
    residual = state -> pack_complex(integrate_probe(unpack_complex(state)))
    driven_bound = real(2frequencies[1]-seed.frequency) < alpha
    admissible = state -> all(isfinite,state) && 0 < real(unpack_complex(state)[1]) < alpha &&
        (real(2frequencies[1]-unpack_complex(state)[1]) < alpha) == driven_bound
    coarse = newton_solve(residual, pack_complex(initial); tolerance=1e-5,
        max_iterations=15, jacobian_backend=:forwarddiff, admissible, progress_callback=progress)
    coarse_state = unpack_complex(coarse.state)
    center = integrate_probe(coarse_state)
    frequency_step = eps(real(coarse_state[1]))
    precision_floor = maximum(norm(integrate_probe(vcat(coarse_state[1]+sign*frequency_step,coarse_state[2:end]))-center)
        for sign in (-1,1))
    tolerance = max(1e-8, 2precision_floor)
    corrected = newton_solve(residual, coarse.state; tolerance, max_iterations=15,
        jacobian_backend=:forwarddiff, admissible, progress_callback=progress)
    solution = unpack_complex(corrected.state)
    unit_profiles = integrate_probe(solution; save_profile=true)
    @assert norm(unit_profiles.residual) <= tolerance
    secondary_profiles = (; (side => (; positions=getproperty(unit_profiles,side).positions,
        states=getproperty(unit_profiles,side).states[secondary_rows,:])
        for side in (:from_horizon,:from_outer))..., residual=unit_profiles.residual[secondary_rows])
    setup = (; spin, alpha, modes=minor_modes, polar, azimuthal, grid)
    primary_for_energy = (; profiles=joint_primary, modes=primary.modes)
    unit_energy = preview_secondary_energy(secondary_profiles, solution[1], primary_for_energy, setup)
    @assert isfinite(unit_energy) && unit_energy > 0
    target_energy = preview.report["energy_"*minor_mode]
    scale = sqrt(target_energy/unit_energy)
    joint_secondary = scaled_preview_profiles(secondary_profiles,scale)
    joint_secondary_energy = preview_secondary_energy(joint_secondary, solution[1],primary_for_energy,setup)
    @assert isapprox(joint_secondary_energy,target_energy; rtol=1e-12,atol=0)
    physical_probe_peak = maximum(maximum(abs, getproperty(unit_profiles,side).states[1:probe_count,:])*scale
        for side in (:from_horizon,:from_outer))
    @assert physical_probe_peak < 1e-30
    references = (preview.references...,
        (; cloud=dominant_mode, model="joint", profiles=joint_primary, modes=primary.modes),
        (; cloud=minor_mode, model="joint", profiles=joint_secondary, modes=setup.modes))
    report = Dict("passed"=>true, "case"=>case, "energy_"*dominant_mode=>joint_primary_energy,
        "energy_"*minor_mode=>joint_secondary_energy, "joint_computed"=>true,
        "joint_regime"=>"Nonlinear primary plus coupled infinitesimal-secondary response of the joint cosine equations; plotted at the retained finite seed normalization.",
        "cloud_m"=>cloud_m,
        "joint_components"=>[Dict("q"=>collect(c.q),"l_values"=>c.l_values) for c in components],
        "omitted_component"=>"(-1,2) starts at secondary amplitude squared; its fundamental feedback is beyond the retained infinitesimal-secondary order.",
        "primary_matching_norm"=>primary_matching, "primary_profile_difference_from_omiya"=>primary_difference,
        "secondary_unit_matching_norm"=>norm(unit_profiles.residual),
        "secondary_matching_tolerance"=>tolerance, "frequency_ulp_residual"=>precision_floor,
        "secondary_frequency_real"=>real(solution[1]),"secondary_frequency_imaginary"=>imag(solution[1]),
        "secondary_amplitude_factor"=>scale,"maximum_physical_probe_coefficient"=>physical_probe_peak,
        "secondary_unit_energy"=>unit_energy,"phase_grid_check"=>phase_error,
        "omiya_diagonal_check"=>diagonal_error,"full_source_difference_errors"=>finite_errors,
        "omiya_profile_recovery_without_feedback"=>control_difference,
        "omiya_profile_recovery_tolerance"=>control_tolerance,
        "spatial_convergence_established"=>false,"time_evolution_performed"=>false,
        "population"=>preview.report["population"],
        "provenance"=>Dict("slurm_job_id"=>get(ENV,"SLURM_JOB_ID","local"),
            "created_utc"=>string(now(UTC)), "source_snapshot"=>normpath(joinpath(@__DIR__,"../../../..")),
            "input_sha256"=>Dict(name=>bytes2hex(sha256(read(joinpath(inputs,name))))
                for name in ("primary.jls","preview.jls"))))
    serialize(prefix * ".jls", (; references, solution, unit_profiles, primary_problem, primary_profile,
        cloud_m, components, report))
    open(prefix * ".toml","w") do stream
        TOML.print(stream,report; sorted=true)
    end
    open(prefix * ".csv","w") do stream
        println(stream,"cloud,model,branch,ell,m,r_star_over_M,R_real,R_imag,R_abs,Q_real,Q_imag")
        for reference in references, side in (:from_horizon,:from_outer)
            branch = getproperty(reference.profiles,side)
            for (j,position) in enumerate(branch.positions), (i,(degree,m)) in enumerate(reference.modes)
                value,slope = branch.states[i,j],branch.states[length(reference.modes)+i,j]
                @printf(stream,"%s,%s,%s,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                    reference.cloud,reference.model,string(side),degree,m,position,real(value),imag(value),
                    abs(value),real(slope),imag(slope))
            end
        end
    end
    println("FINISHED joint comparison, matched energies=",joint_primary_energy,", ",joint_secondary_energy)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && joint_probe_main()
