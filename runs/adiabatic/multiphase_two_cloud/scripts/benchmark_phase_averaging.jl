# Compare analytic phase averaging with the unchanged direct tensor sums on
# dilute and late profiles, for both potentials and the current [-80,180] domain.
using Serialization, TOML, SHA, LinearAlgebra
include(joinpath(@__DIR__, "../../../../src/adiabatic/multiphase/load.jl"))

function phase_benchmark_write(path, report)
    open(path * ".tmp", "w") do io
        TOML.print(io, report; sorted=true)
    end
    mv(path * ".tmp", path; force=true)
end

phase_benchmark_relative(a,b) = norm(a-b)/max(norm(a),norm(b),floatmin(Float64))

function phase_benchmark_comparison(candidate, reference)
    relative = phase_benchmark_relative
    flux_names = (:horizon_energy_flux, :outer_energy_flux,
        :horizon_angular_momentum_flux, :outer_angular_momentum_flux)
    fluxes(sample) = [getproperty(sample.response.observables,k) for k in flux_names]
    rate_errors = abs.(candidate.balance.log_rates-reference.balance.log_rates) ./
        max.(abs.(candidate.balance.log_rates),abs.(reference.balance.log_rates),floatmin(Float64))
    derivative_errors = [relative(candidate.response.observable_derivatives[:,i],
        reference.response.observable_derivatives[:,i]) for i in 1:2]
    result = Dict{String,Any}(
        "charge_relative_error"=>relative(candidate.charges,reference.charges),
        "flux_vector_relative_error"=>relative(fluxes(candidate),fluxes(reference)),
        "individual_flux_absolute_differences"=>fluxes(candidate)-fluxes(reference),
        "charge_response_column_relative_errors"=>derivative_errors,
        "state_relative_error"=>relative(candidate.solution.state,reference.solution.state),
        "log_rate_relative_errors"=>rate_errors,
        "direct_log_rates"=>reference.balance.log_rates,
        "averaged_log_rates"=>candidate.balance.log_rates)
    result["passed"] = result["charge_relative_error"] < 1e-11 &&
        result["flux_vector_relative_error"] < 1e-8 &&
        maximum(derivative_errors) < 1e-9 && maximum(rate_errors) < 1e-6
    return result
end

function benchmark_phase_averaging(data_root, prefix)
    any(ispath(prefix*suffix) for suffix in (".toml",".jls")) && error("use a fresh benchmark prefix")
    mkpath(dirname(prefix))
    report = Dict{String,Any}("status"=>"running", "change"=>"analytic quadratic phase averages only; AD batching unchanged",
        "source_snapshot"=>normpath(joinpath(@__DIR__,"../../../..")),
        "repetitions"=>3, "cases"=>Dict{String,Any}[])
    write_report() = phase_benchmark_write(prefix*".toml",report)
    write_report()
    try
        for potential in (:chiral,:cosine)
            initial_name = "$(potential)_growth_inner80_20261102_v1"
            initial_path = joinpath(data_root,"results",initial_name*"-initial.jls")
            saved = deserialize(initial_path)
            problem, dilute = saved[1],saved[2]
            problem.potential == potential && problem.integration.inner == -80 &&
                problem.integration.outer == 180 || error("unexpected benchmark model")
            late_path=joinpath(data_root,"work/validation_20261102_v1/outputs/seeds",string(potential)*"-late.jls")
            old=deserialize(late_path)
            old_problem, late = old isa Tuple ? (old[1],old[2]) : (old.problem,old.sample)
            old_profile=unpack_two_phase_profile(late.solution.state,late.solution.amplitudes,old_problem)
            late_guess=pack_two_phase_profile(old_profile,problem)
            for (stage,amplitudes,guess) in (("dilute",dilute.solution.amplitudes,dilute.solution.state),
                    ("late",late.solution.amplitudes,late_guess))
                case=Dict{String,Any}("name"=>"$(potential)_$stage", "potential"=>string(potential),
                    "amplitudes"=>amplitudes, "integration"=>Dict(string(k)=>v for (k,v) in pairs(problem.integration)),
                    "input_sha256"=>Dict(basename(p)=>bytes2hex(sha256(read(p))) for p in (initial_path,late_path)),
                    "samples"=>Dict{String,Any}[], "steps"=>Dict{String,Any}[])
                push!(report["cases"],case);write_report()
                println("CASE ",case["name"]);flush(stdout)
                solution=correct_two_phase_profile(problem,amplitudes,guess;tolerance=1e-9,max_iterations=20)
                perturbed=copy(solution.state);perturbed[1]+=1e-9
                samples=Dict{Symbol,Any}()
                for repetition in 1:3, method in (isodd(repetition) ? (:direct,:averaged) : (:averaged,:direct))
                    GC.gc()
                    measured=@timed two_phase_evolution_sample(problem,amplitudes,copy(perturbed);
                        tolerance=1e-9,max_iterations=20,observable_method=method)
                    samples[method]=measured.value
                    push!(case["samples"],Dict("method"=>string(method),"repetition"=>repetition,
                        "seconds"=>measured.time,"bytes"=>measured.bytes,"iterations"=>measured.value.solution.iterations))
                    write_report()
                end
                comparison=phase_benchmark_comparison(samples[:averaged],samples[:direct])
                case["comparison"]=comparison;write_report()
                comparison["passed"] || error("sample equivalence failed for $(case["name"])")
                # Same accepted initial state for both methods; each step independently
                # corrects its midpoint and endpoint. Repetition 1 warms compilation.
                initial=samples[:direct]
                spacing=min(66495.25029285335,0.01/maximum(abs.(initial.balance.log_rates)))
                case["step_spacing"]=spacing
                histories=Dict{Symbol,Any}()
                for repetition in 1:3, method in (isodd(repetition) ? (:direct,:averaged) : (:averaged,:direct))
                    evaluate=(a,g)->two_phase_evolution_sample(problem,a,g;tolerance=1e-9,observable_method=method)
                    GC.gc()
                    measured=@timed evolve_two_phase(evaluate,initial;final_time=spacing,time_step=spacing)
                    histories[method]=measured.value
                    push!(case["steps"],Dict("method"=>string(method),"repetition"=>repetition,
                        "seconds"=>measured.time,"bytes"=>measured.bytes))
                    write_report()
                end
                direct,averaged=histories[:direct],histories[:averaged]
                case["step_comparison"]=phase_benchmark_comparison(averaged.samples[end],direct.samples[end])
                case["step_amplitude_relative_errors"]=abs.(averaged.samples[end].solution.amplitudes ./
                    direct.samples[end].solution.amplitudes .- 1)
                case["step_comparison"]["passed"] && maximum(case["step_amplitude_relative_errors"]) < 1e-9 ||
                    error("step equivalence failed for $(case["name"])")
                for kind in ("samples","steps")
                    best(method)=minimum(row["seconds"] for row in case[kind]
                        if row["method"]==method && row["repetition"]>1)
                    case[kind*"_speedup"]=best("direct")/best("averaged")
                end
                write_report()
                println("PASS ",case["name"]," sample speedup=",case["samples_speedup"],
                    " step speedup=",case["steps_speedup"]," rate errors=",comparison["log_rate_relative_errors"])
                flush(stdout)
            end
        end
        report["status"]="complete";report["passed"]=true;write_report()
    catch exception
        report["status"]="failed";report["failure"]=sprint(showerror,exception);write_report()
        rethrow()
    end
    return report
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS)==2 || error("usage: benchmark_phase_averaging.jl CAMPAIGN_DATA_ROOT OUTPUT_PREFIX")
    benchmark_phase_averaging(abspath.(ARGS)...)
end
