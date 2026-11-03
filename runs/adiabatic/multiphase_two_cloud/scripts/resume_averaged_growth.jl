# Continue an accepted direct-evaluator trajectory with explicit mixed-source
# provenance. Original samples and the canonical output time grid are retained.
using SHA, TOML
include(joinpath(@__DIR__, "chiral_growth.jl"))

function resume_averaged_growth(parent, prefix)
    any(ispath(prefix*s) for s in ("-configuration.jls", "-steps1000.jls", ".lock")) &&
        error("use a fresh continuation prefix")
    mkdir(prefix*".lock")
    report = Dict{String,Any}()
    try
        manifest = read(joinpath(@__DIR__, "../../../../provenance.json"), String)
        parent_configuration_bytes = read(parent*"-configuration.jls")
        parent_configuration = deserialize(IOBuffer(parent_configuration_bytes))
        report = TOML.parsefile(parent*"-report.toml")
        steps = report["steps"]
        initial_bytes = read(parent*"-initial.jls")
        problem, initial, frequencies = deserialize(IOBuffer(initial_bytes))
        checkpoint_bytes = read(parent*"-steps$(steps).jls")
        history = deserialize(IOBuffer(checkpoint_bytes))
        times = collect(range(0, report["requested_final_time"]; length=steps+1))
        n = length(history.times)
        1 <= n <= length(times) && history.times == times[1:n] || error("checkpoint time-grid mismatch")
        problem.potential == :chiral || error("this continuation requires the chiral trajectory")
        report["spatial_settings"] == Dict(string(k)=>v for (k,v) in pairs(problem.integration)) || error("domain mismatch")
        configuration = (; source_manifest=manifest, parent_prefix=parent,
            parent_source_manifest=parent_configuration.source_manifest,
            parent_configuration_sha256=bytes2hex(sha256(parent_configuration_bytes)),
            parent_checkpoint_sha256=bytes2hex(sha256(checkpoint_bytes)),
            parent_initial_sha256=bytes2hex(sha256(initial_bytes)),
            resume_step=n-1, observable_method=:averaged)
        save_checkpoint(prefix*"-configuration.jls", configuration)
        write(prefix*"-resume-source.jls", checkpoint_bytes)
        write(prefix*"-initial.jls", initial_bytes)
        cp(parent*"-linear-reference.csv", prefix*"-linear-reference.csv")
        report["termination"] = "running"
        report["parent_prefix"] = parent
        report["resume_step"] = n-1
        report["observable_method_after_resume"] = "averaged"
        report["resume_source_sha256"] = configuration.parent_checkpoint_sha256
        write_chiral_report(prefix, report)
        previous = history.samples[end]
        evaluate = (a,g) -> two_phase_evolution_sample(problem,a,g;
            tolerance=report["matching_tolerance"], observable_method=:averaged)
        current = evaluate(previous.solution.amplitudes, previous.solution.state)
        differences = abs.(current.balance.log_rates-previous.balance.log_rates) ./
            max.(abs.(current.balance.log_rates),abs.(previous.balance.log_rates),floatmin(Float64))
        maximum(differences) < 1e-6 || error("averaged/direct rate mismatch at continuation: $differences")
        report["resume_rate_relative_differences"] = differences
        write_chiral_report(prefix, report)
        println("RESUME VERIFIED at step ", n-1, "; relative rate differences=",differences)
        flush(stdout)
        save_checkpoint(prefix*"-steps$(steps).jls", history)
        write_trajectory(prefix*"-steps$(steps).csv", history)
        for index in n+1:length(times)
            spacing = times[index]-times[index-1]
            trial = evolve_two_phase(evaluate, current; final_time=spacing, time_step=spacing)
            current = trial.samples[end]
            cumulative = history.integrated_losses[end]+trial.integrated_losses[end]
            push!(history.times,times[index])
            push!(history.samples,current)
            push!(history.integrated_losses,cumulative)
            push!(history.balance_errors,current.charges-initial.charges+cumulative)
            save_checkpoint(prefix*"-steps$(steps).jls",history)
            write_trajectory(prefix*"-steps$(steps).csv",history)
            print_sample("RESUMED step=$(index-1)/$steps, t=$(times[index])",current)
        end
        report["termination"] = "final_time"
        report["last_time"] = history.times[end]
        report["final_amplitudes"] = current.solution.amplitudes
        report["amplitude_growth_factors"] = current.solution.amplitudes ./ initial.solution.amplitudes
        report["final_charges"] = current.charges
        report["relative_charge_budget"] = history.balance_errors[end] ./ max.(abs.(initial.charges),abs.(current.charges))
        write_chiral_report(prefix,report)
        println("COMPLETED averaged chiral continuation")
    catch exception
        write(prefix*"-failure.txt",sprint(showerror,exception,catch_backtrace()))
        if !isempty(report)
            report["termination"]="failed"
            report["failure"]=sprint(showerror,exception)
            write_chiral_report(prefix,report)
        end
        rethrow()
    finally
        rm(prefix*".lock")
    end
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    length(ARGS)==2 || error("usage: resume_averaged_growth.jl PARENT_PREFIX NEW_PREFIX")
    resume_averaged_growth(abspath.(ARGS)...)
end
