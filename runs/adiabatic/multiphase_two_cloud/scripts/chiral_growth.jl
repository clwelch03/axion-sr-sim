using SHA, TOML
include(joinpath(@__DIR__, "response_and_evolution.jl"))

# Exploratory growth: preserve the pilot's fixed BH and numerical settings
# except for the radial boundaries, start 100 times lower in both amplitudes,
# and follow ten initial 211 amplitude e-folding times. Time is t/M.
chiral_growth_settings() = (;
    seed_factor=0.01, initial_efolds=10.0, steps=1000,
    inner=-80.0, outer=180.0,
    matching_tolerance=1e-9, initial_rate_relative_tolerance=0.01,
)

function chiral_growth_problem(template, frequencies; potential=:chiral,
        inner=chiral_growth_settings().inner,
        outer=chiral_growth_settings().outer)
    potential in (:chiral, :cosine) || error("unsupported growth potential")
    template.cloud_m == [1, 2] || error("require cloud order 211, 322")
    template.spin == 0.99 && template.alpha == 0.42 ||
        error("this pilot requires spin=0.99 and alpha=0.42")
    [(only(component.l_values), component.q) for component in template.components] ==
        [(1,(1,0)), (2,(0,1)), (0,(2,-1)), (3,(-1,2))] ||
        error("require the pilot's two fundamentals and two driven components")
    (; free_outer_coefficient_indices) = two_phase_coefficient_indices(template.components)
    return TwoPhaseProblem(;
        template.spin, template.alpha, template.components, template.cloud_m,
        template.couplings, template.polar_values, template.grid,
        template.phase1_nodes, template.phase2_nodes,
        potential,
        # Keep the original numerical coefficient scale fixed throughout the
        # run. Changing the physical seeds does not require changing this scale.
        template.coefficient_scale, free_outer_coefficient_indices,
        initial_bound=BitVector([real(sum(component.q .* frequencies)) < template.alpha
            for component in template.components]),
        integration=(; template.integration..., inner, outer),
    )
end

function chiral_growth_initial(template, seed, settings; potential=:chiral)
    amplitudes = settings.seed_factor .* seed.solution.amplitudes
    guesses = unpack_complex(seed.solution.state)[1:2]
    template.integration.horizon_order == 3 && template.integration.outer_order == 7 ||
        error("linear reference requires boundary series orders 3/7")
    linear_grid = (; inner=settings.inner, outer=settings.outer,
        template.integration.match, template.integration.step)
    linear = map(1:2) do cloud
        solve_linear_mode(guesses[cloud], template.spin, template.alpha, [cloud], cloud,
            template.couplings; linear_grid..., tolerance=settings.matching_tolerance)
    end
    frequencies = [mode.frequency for mode in linear]
    problem = chiral_growth_problem(template, frequencies;
        potential, inner=settings.inner, outer=settings.outer)
    # These scaled linear modes are only a guess. Re-solve the joint nonlinear
    # field, including both driven components, before evaluating any rate.
    profile = TwoPhaseProfile(frequencies,
        ComplexF64[amplitudes[1]*only(linear[1].horizon_amplitudes),
            amplitudes[2]*only(linear[2].horizon_amplitudes), 0, 0],
        ComplexF64[amplitudes; 0; 0])
    initial = two_phase_evolution_sample(problem, amplitudes,
        pack_two_phase_profile(profile, problem); tolerance=settings.matching_tolerance)
    return problem, initial, frequencies
end

function write_chiral_report(prefix, report)
    open(prefix * "-report.toml.tmp", "w") do io
        TOML.print(io, report; sorted=true)
    end
    mv(prefix * "-report.toml.tmp", prefix * "-report.toml"; force=true)
end

"""
    run_joint_growth(seed_path, prefix; potential=:chiral, reference_path=nothing)

Run a fixed-BH trajectory using the old pilot only for its numerical
configuration, amplitude normalization, and linear frequency guesses, with
the radial boundaries set by `chiral_growth_settings`. Fresh linear modes on
that domain seed a corrected profile. An optional reference report fixes the
comparison's seeds, background, spatial settings, step count, and end time.
Checkpoints belong to one frozen
source/configuration and retain separate boundary fluxes and charge budgets.

The existing fixed-step midpoint method has a log-step guard, not an adaptive
time-error estimate. This finite-domain run is exploratory. A failed profile
or balance solve preserves accepted checkpoints and reports the error; it does
not relax tolerances or stop merely because a cloud starts to decline.
"""
function run_joint_growth(seed_path, prefix; potential=:chiral, reference_path=nothing)
    settings = chiral_growth_settings()
    reference = isnothing(reference_path) ? nothing : TOML.parsefile(reference_path)
    mkpath(dirname(prefix))
    mkdir(prefix * ".lock")
    report = Dict{String,Any}()
    try
        manifest_path = joinpath(@__DIR__, "../../../../provenance.json")
        isfile(manifest_path) || error("run from a frozen source snapshot")
        configuration = (; source_manifest=read(manifest_path, String),
            seed_sha256=bytes2hex(sha256(read(seed_path))), settings,
            potential, up_quark_mass=UP_QUARK_MASS,
            down_quark_mass=DOWN_QUARK_MASS, kappa=CHIRAL_SQRT_PREFACTOR,
            reference=isnothing(reference_path) ? "" : read(reference_path, String))
        configuration_path = prefix * "-configuration.jls"
        if isfile(configuration_path)
            deserialize(configuration_path) == configuration ||
                error("changed configuration; use a distinct output prefix")
        else
            any(isfile(prefix * suffix) for suffix in
                ("-initial.jls", "-steps$(settings.steps).jls", "-report.toml")) &&
                error("existing outputs lack their configuration")
            save_checkpoint(configuration_path, configuration)
        end

        initial_path = prefix * "-initial.jls"
        if isfile(initial_path)
            problem, initial, frequencies = deserialize(initial_path)
        else
            template, seed = deserialize(seed_path)
            problem, initial, frequencies = chiral_growth_initial(template, seed, settings; potential)
            save_checkpoint(initial_path, (problem, initial, frequencies))
        end
        print_sample("Fresh $potential initial state", initial)
        linear_rates = imag.(frequencies)
        all(linear_rates .> 0) || error("require two superradiant linear reference modes")
        initial.balance.log_rates[1] > 0 || error("initial 211 amplitude is not growing")
        rate_difference = initial.balance.log_rates ./ linear_rates .- 1
        final_time = if isnothing(reference)
            settings.initial_efolds / initial.balance.log_rates[1]
        else
            reference["alpha"] == problem.alpha && reference["spin"] == problem.spin ||
                error("reference background mismatch")
            reference["initial_amplitudes"] == initial.solution.amplitudes ||
                error("reference seed mismatch")
            reference["steps"] == settings.steps || error("reference step-count mismatch")
            reference["coefficient_scale"] == problem.coefficient_scale ||
                error("reference coefficient scale mismatch")
            reference["spatial_settings"] == Dict(string(k)=>v for (k,v) in pairs(problem.integration)) ||
                error("reference spatial settings mismatch")
            reference["requested_final_time"]
        end
        report = Dict{String,Any}(
            "model"=>"joint $potential, fixed BH, exploratory finite domain",
            "potential"=>string(potential), "alpha"=>problem.alpha, "spin"=>problem.spin,
            "quartic_coefficient"=>potential == :chiral ? 1-3CHIRAL_SQRT_PREFACTOR/4 : 1.0,
            "initial_amplitudes"=>initial.solution.amplitudes,
            "initial_charges"=>initial.charges,
            "seed_prescription"=>"small classical outer amplitudes; not occupation-one seeds",
            "linear_amplitude_rates"=>linear_rates,
            "initial_log_amplitude_rates"=>initial.balance.log_rates,
            "initial_rate_relative_difference"=>rate_difference,
            "requested_initial_211_efolds"=>settings.initial_efolds,
            "requested_final_time"=>final_time, "time_unit"=>"t/M",
            "time_reference"=>isnothing(reference) ? "own initial 211 rate" :
                get(reference, "time_reference", "matched chiral run"),
            "steps"=>settings.steps, "time_step"=>final_time/settings.steps,
            "time_error_control"=>"fixed midpoint steps; max log-step guard 0.05",
            "matching_tolerance"=>settings.matching_tolerance,
            "coefficient_scale"=>problem.coefficient_scale,
            "retained_components"=>["211", "322", "q=(2,-1), l=0", "q=(-1,2), l=3"],
            "spatial_settings"=>Dict(string(k)=>v for (k,v) in pairs(problem.integration)),
            "termination"=>"running",
        )
        potential == :chiral && (report["quark_mass_ratio"] = UP_QUARK_MASS/DOWN_QUARK_MASS)
        write_chiral_report(prefix, report)
        println("Initial E/J rates relative to linear modes: ", rate_difference)
        println("Requested final t/M=", final_time, "; steps=", settings.steps)
        flush(stdout)
        maximum(abs, rate_difference) <= settings.initial_rate_relative_tolerance ||
            error("dilute initial rate differs from linear reference by more than 1%")

        open(prefix * "-linear-reference.csv", "w") do io
            println(io, "cloud,initial_amplitude,linear_frequency_real,linear_growth_rate,initial_log_rate")
            for cloud in 1:2
                println(io, join((cloud, initial.solution.amplitudes[cloud],
                    real(frequencies[cloud]), linear_rates[cloud], initial.balance.log_rates[cloud]), ','))
            end
        end
        history = run_trajectory(problem, initial, prefix, final_time, settings.steps)
        report["termination"] = "final_time"
        report["last_time"] = last(history.times)
        report["final_amplitudes"] = history.samples[end].solution.amplitudes
        report["amplitude_growth_factors"] = report["final_amplitudes"] ./ initial.solution.amplitudes
        report["final_charges"] = history.samples[end].charges
        report["relative_charge_budget"] = history.balance_errors[end] ./
            max.(abs.(initial.charges), abs.(history.samples[end].charges))
        write_chiral_report(prefix, report)
        println("COMPLETED $potential trajectory: ", report["amplitude_growth_factors"])
    catch exception
        failure = sprint(showerror, exception, catch_backtrace())
        open(prefix * "-failure.txt", "w") do io
            println(io, failure)
        end
        # Do not overwrite an existing report on a rejected resume/configuration.
        if !isempty(report)
            report["termination"] = "failed"
            report["failure"] = sprint(showerror, exception)
            checkpoint = prefix * "-steps$(settings.steps).jls"
            if isfile(checkpoint)
                report["last_time"] = last(deserialize(checkpoint).times)
            end
            write_chiral_report(prefix, report)
        end
        rethrow()
    finally
        rm(prefix * ".lock")
    end
end

run_chiral_growth(seed_path, prefix) = run_joint_growth(seed_path, prefix)

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    length(ARGS) in (2, 3) ||
        error("usage: chiral_growth.jl SEED_JLS [REFERENCE_TOML] OUTPUT_PREFIX")
    run_joint_growth(abspath(ARGS[1]), abspath(ARGS[end]);
        reference_path=length(ARGS) == 3 ? abspath(ARGS[2]) : nothing)
end
