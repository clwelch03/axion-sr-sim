# Omiya cosine control on the paper's [-100,100] domain at alpha=0.42.
# Reuse the established branch solves and energy/angular-momentum equations.
using TOML
include(joinpath(@__DIR__, "omiya_comparison.jl"))

function write_omiya_growth_report(prefix, report)
    open(prefix * "-report.toml.tmp", "w") do io
        TOML.print(io, report; sorted=true)
    end
    mv(prefix * "-report.toml.tmp", prefix * "-report.toml"; force=true)
end

function omiya_growth_amplitudes(initial_amplitude)
    # Include the exact initial amplitude, a lower buffer, and half an extra
    # linear e-fold beyond the common ten-e-fold time interval. Never extrapolate.
    return initial_amplitude .* exp.(collect(-2:210) ./ 20)
end

function omiya_growth_problem(template)
    # Omiya et al. 2211.01949v1, p. 11: endpoints +/-100 at alpha=0.42.
    # Match the joint run's seeds and time interval, not its numerical domain.
    integration = (; template.integration..., inner=-100.0, outer=100.0)
    return (; (name=>getproperty(template, name) for name in propertynames(template))...,
        integration)
end

# Existing reproduction prescription for the outgoing m=3 solution. The
# far-zone, source-free leg is separate from the [-100,100] cloud domain.
omiya_growth_radiative_options() = (; far_position=5000.0, far_step=0.25, far_order=16)

function run_omiya_growth(seed_path, reference_path, prefix)
    reference = TOML.parsefile(reference_path)
    mkpath(dirname(prefix))
    mkdir(prefix * ".lock")
    report = Dict{String,Any}()
    try
        template, seed = deserialize(seed_path)
        problem = omiya_growth_problem(template)
        problem.potential == :cosine && problem.cloud_m == [1, 2] || error("unsupported model")
        problem.alpha == reference["alpha"] == 0.42 &&
            problem.spin == reference["spin"] == 0.99 || error("reference background mismatch")
        [(component.q, component.l_values) for component in problem.components] ==
            [((1,0),[1]), ((0,1),[2]), ((2,-1),[0]), ((-1,2),[3])] ||
            error("require the pilot's four components")
        problem.integration.horizon_order == 3 && problem.integration.outer_order == 7 ||
            error("require boundary orders 3/7")
        initial_amplitudes = Float64.(reference["initial_amplitudes"])
        initial_amplitudes == 0.01 .* seed.solution.amplitudes || error("reference seed mismatch")
        amplitudes = omiya_growth_amplitudes(initial_amplitudes[1])
        times = collect(range(0, reference["requested_final_time"]; length=reference["steps"]+1))
        configuration = (;
            manifest=read(joinpath(@__DIR__, "../../../../provenance.json"), String),
            seed_sha256=bytes2hex(sha256(read(seed_path))),
            reference=read(reference_path, String), amplitudes, times,
            potential=:cosine, primary_step=0.1, secondary_step=0.2,
            integration=problem.integration,
            radiative_options=omiya_growth_radiative_options(),
            subdivisions=8, time_tolerance=1e-8, table_tolerance=1e-5)
        checkpoint_path = prefix * "-branch.jls"
        points = if isfile(checkpoint_path)
            saved = deserialize(checkpoint_path)
            saved.configuration == configuration || error("comparison configuration changed")
            saved.points
        else
            any(isfile(prefix * suffix) for suffix in ("-report.toml", "-trajectory.csv")) &&
                error("existing outputs lack their configuration")
            NamedTuple[]
        end
        report = Dict{String,Any}(
            "model"=>"Omiya primary/secondary cosine, fixed BH, exploratory finite domain",
            "potential"=>"cosine", "alpha"=>problem.alpha, "spin"=>problem.spin,
            "initial_amplitudes"=>initial_amplitudes,
            "normalization"=>"matched outer coefficients; total charges need not match",
            "requested_final_time"=>last(times), "time_unit"=>"t/M",
            "output_intervals"=>reference["steps"], "rk4_subdivisions"=>8,
            "primary_step"=>0.1, "secondary_step"=>0.2,
            "inner"=>problem.integration.inner, "outer"=>problem.integration.outer,
            "spatial_settings"=>Dict(string(k)=>v for (k,v) in pairs(problem.integration)),
            "reference_spatial_settings"=>reference["spatial_settings"],
            "radiative_boundary"=>Dict(string(k)=>v for (k,v) in pairs(configuration.radiative_options)),
            "branch_points_requested"=>length(amplitudes),
            "branch_points_completed"=>length(points), "termination"=>"building_branch")
        write_omiya_growth_report(prefix, report)
        (; inner, outer, match) = problem.integration
        guesses = unpack_complex(seed.solution.state)[1:2]
        linear_primary = solve_linear_mode(guesses[1], problem.spin, problem.alpha,
            [1], 1, problem.couplings; inner, outer, match, step=0.1, tolerance=1e-9)
        linear_secondary = solve_linear_mode(guesses[2], problem.spin, problem.alpha,
            [2], 2, problem.couplings; inner, outer, match, step=0.2, tolerance=1e-8)
        for amplitude in amplitudes[length(points)+1:end]
            point = omiya_branch_point(amplitude, problem, linear_primary, linear_secondary;
                radiative_options=configuration.radiative_options)
            isempty(points) || point.row.primary_energy > last(points).row.primary_energy ||
                error("primary energy branch is not monotone")
            push!(points, point)
            save_comparison_checkpoint(checkpoint_path, (; configuration, points))
            write_comparison_rows(prefix * "-branch.csv", [p.row for p in points])
            report["branch_points_completed"] = length(points)
            report["last_branch_amplitude"] = amplitude
            write_omiya_growth_report(prefix, report)
        end
        report["termination"] = "integrating"
        write_omiya_growth_report(prefix, report)
        rows = [point.row for point in points]
        trajectory = omiya_trajectory(rows, initial_amplitudes, times; subdivisions=8)
        coarser_time = omiya_trajectory(rows, initial_amplitudes, times; subdivisions=4)
        initial_index = findfirst(row -> row.amplitude == initial_amplitudes[1], rows)
        coarse_indices = sort(unique(vcat(collect(1:2:length(rows)), initial_index, length(rows))))
        coarser_table = omiya_trajectory(rows[coarse_indices], initial_amplitudes, times)
        discrepancy(other) = maximum(abs(getproperty(a, name)/getproperty(b, name)-1)
            for (a,b) in zip(trajectory, other) for name in (:A1, :A2))
        report["time_refinement_amplitude_difference"] = discrepancy(coarser_time)
        report["table_thinning_amplitude_difference"] = discrepancy(coarser_table)
        write_omiya_growth_report(prefix, report)
        report["time_refinement_amplitude_difference"] < configuration.time_tolerance ||
            error("time refinement failed")
        report["table_thinning_amplitude_difference"] < configuration.table_tolerance ||
            error("table thinning failed")
        write_comparison_rows(prefix * "-trajectory.csv.tmp", trajectory)
        mv(prefix * "-trajectory.csv.tmp", prefix * "-trajectory.csv"; force=true)
        save_comparison_checkpoint(prefix * "-summary.jls",
            (; configuration, initial_amplitudes, final=last(trajectory), report))
        report["termination"] = "final_time"
        report["last_time"] = last(trajectory).time
        report["final_amplitudes"] = [last(trajectory).A1, last(trajectory).A2]
        report["amplitude_growth_factors"] = report["final_amplitudes"] ./ initial_amplitudes
        write_omiya_growth_report(prefix, report)
        println("COMPLETED Omiya cosine growth: ", report["amplitude_growth_factors"])
        flush(stdout)
    catch exception
        open(prefix * "-failure.txt", "w") do io
            println(io, sprint(showerror, exception, catch_backtrace()))
        end
        if !isempty(report)
            report["termination"] = "failed"
            report["failure"] = sprint(showerror, exception)
            write_omiya_growth_report(prefix, report)
        end
        rethrow()
    finally
        rm(prefix * ".lock")
    end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    length(ARGS) == 3 || error("usage: omiya_growth.jl SEED_JLS REFERENCE_TOML OUTPUT_PREFIX")
    run_omiya_growth(abspath.(ARGS)...)
end
