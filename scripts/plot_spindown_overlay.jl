# Overlay of two to four finished spindown runs with the same initial state (a₀, α₀, f_a, seed occupations), e.g. the
# cosine and free potentials of scripts/evolution_spindown.jl and the hydrogenic evolution of scripts/hydrogenic_spindown.jl.
# Reads results/<prefix>_diagnostics.csv for each. Rows: occupations 𝒩_i (f_a² M₀²; same f_a in all runs, so directly
# comparable), spin a = j/b² with each run's |211⟩ threshold, mass M/M₀. Line style marks the run (solid, dashed,
# dash-dot in argument order). With --markers the last run is drawn as open circles instead, for a benchmark expected
# to lie on another run (e.g. the Leaver ODE of scripts/leaver_spindown.py on the free run). Left column: linear time
# up to the end of the shortest run. Right column: all runs on a log time axis; a run that ends before the window does
# is continued by a thin dotted line at its final value. On the log axis lines pass through points interpolated between
# accepted states (plot_helpers.jl). --compact: only the log-time column, occupation and spin rows, larger type, for
# slides; output name gets _compact.
#
# Rows from evolution_spindown.jl log every rhs evaluation; only accepted states are kept, with the filter of
# scripts/plot_spindown.jl (see its header). A CSV without repeated times (one row per saved state, as from
# hydrogenic_spindown.jl) is used as is. The threshold is the CSV's `threshold` column when present (hydrogenic: the
# zero of its rate formula), otherwise a = 4Ω/(1 + 4Ω²) with Ω = Re ω₂₁₁.
#
# Run from the package root:
#   julia --project=scripts scripts/plot_spindown_overlay.jl [--markers] [--compact] <prefix 1> <label 1> <prefix 2> <label 2> [<prefix 3> <label 3> [<prefix 4> <label 4>]] [output.png]
# e.g. spindown_a0.9_alpha0.3_fa3e-3_t1e10 cosine spindown_free_a0.9_alpha0.3_fa3e-3 free \
#      hydrogenic_spindown_a0.9_alpha0.3_fa3e-3 hydrogenic
using CairoMakie, DelimitedFiles, Printf
include(joinpath(@__DIR__, "plot_helpers.jl"))

flags = filter(startswith("--"), ARGS)
issubset(flags, ["--markers", "--compact"]) || error("unknown option among $flags")
last_as_markers, compact = "--markers" in flags, "--compact" in flags
arguments = filter(!startswith("--"), ARGS)
4 ≤ length(arguments) ≤ 9 ||
    error("usage: plot_spindown_overlay.jl [--markers] <prefix 1> <label 1> <prefix 2> <label 2> [<prefix 3> <label 3> [<prefix 4> <label 4>]] [output.png]")
run_count = length(arguments) ÷ 2
prefixes = Tuple(arguments[2index - 1] for index in 1:run_count)
labels = Tuple(arguments[2index] for index in 1:run_count)
output_path = isodd(length(arguments)) ? arguments[end] :
              joinpath(@__DIR__, "..", "plots", "spindown_overlay_$(join(labels, "_"))_" *
                       replace(prefixes[1], r"^spindown_(free_)?" => "", r"_t[^_]+$" => "") * (compact ? "_compact" : "") * ".png")

"""Accepted states of one run as a named tuple of columns."""
function read_run(prefix)
    table, header = readdlm(joinpath(@__DIR__, "..", "results", prefix * "_diagnostics.csv"), ','; header=true)
    column(name) = table[:, findfirst(==(name), vec(header))]
    raw_times = column("time")
    earliest_later_time = similar(raw_times)
    running_minimum = Inf
    for index in reverse(eachindex(raw_times))
        earliest_later_time[index] = running_minimum
        running_minimum = min(running_minimum, raw_times[index])
    end
    logs_every_stage = any(raw_times[2:end] .== raw_times[1:end-1])
    accepted = [!logs_every_stage || index == 1 ||
                (raw_times[index] ≤ earliest_later_time[index] && raw_times[index] == raw_times[index - 1])
                for index in eachindex(raw_times)]
    times = raw_times[accepted]
    issorted(times) || error("$prefix: accepted states are not in time order")
    if "threshold" in vec(header)
        threshold = column("threshold")[accepted]
    else
        frequency_211 = column("Re_omega1")[accepted]
        threshold = @. 4frequency_211 / (1 + 4frequency_211^2)
    end
    (time=times, occupation_211=column("N1")[accepted], occupation_322=column("N2")[accepted],
     spin=column("a")[accepted], mass=column("b")[accepted], threshold)
end
runs = read_run.(prefixes)
for (run, label) in zip(runs, labels)
    peak = argmax(run.occupation_211)
    @printf("%-8s %4d accepted states to t = %.3e; peak 𝒩₂₁₁ = %.1f at %.3e; final 𝒩 = (%.4g, %.4g), a = %.5f (threshold %.5f), M/M₀ = %.6f, J/M₀² = %.5f\n", label, length(run.time), run.time[end],
            run.occupation_211[peak], run.time[peak], run.occupation_211[end], run.occupation_322[end], run.spin[end],
            run.threshold[end], run.mass[end], run.spin[end] * run.mass[end]^2)
end

# case parameters from the first prefix
spin_label = match(r"_a([0-9.]+)_", prefixes[1])[1]
mass_coupling_label = match(r"_alpha([0-9.]+)_", prefixes[1])[1]
decay_constant_ratio = parse(Float64, match(r"_fa([0-9.eE+-]+?)(_t|$)", prefixes[1])[1])
decay_constant_gev = decay_constant_ratio * 1.2209e19

"""`x` as "3.7 × 10¹⁶" (or "10¹⁶" when the mantissa is 1); as in scripts/plot_spindown.jl."""
function scientific(x; digits=2)
    exponent = floor(Int, log10(abs(x)))
    mantissa = round(x / 10.0^exponent; sigdigits=digits)
    mantissa ≥ 10 && (mantissa /= 10; exponent += 1)
    superscript = join(Dict('-' => '⁻', '0' => '⁰', '1' => '¹', '2' => '²', '3' => '³', '4' => '⁴', '5' => '⁵',
                            '6' => '⁶', '7' => '⁷', '8' => '⁸', '9' => '⁹')[c] for c in string(exponent))
    exponent == 0 && return string(mantissa)
    mantissa == 1 ? "10" * superscript : @sprintf("%g × 10%s", mantissa, superscript)
end

# repo palette (main: src/common/plotstyle.py): |211⟩ green, |322⟩ blue, neutral ink
mode_211_color, mode_322_color, ink = colorant"#3DAE2B", colorant"#00A0DF", colorant"#0b0b0b"
# solid, dashed, dash-dot in argument order; with --markers the last run is open circles
last_as_markers || run_count < 4 || error("four runs need --markers (there are three line styles)")
line_styles = [last_as_markers && index == run_count ? :markers : (:solid, :dash, :dashdot)[index] for index in 1:run_count]
grid = (:black, 0.06)
figure = compact ? Figure(size=(900, 850), fontsize=20, figure_padding=(16, 48, 16, 16)) : Figure(size=(1300, 1000), fontsize=15)
line_width = compact ? 2.5 : 2
title = compact ? "$(join(labels, " vs ")): a₀ = $spin_label, α₀ = $mass_coupling_label, f_a = $(scientific(decay_constant_gev)) GeV" :
        "Spin-down, $(join(labels, " / ")): a₀ = $spin_label, α₀ = $mass_coupling_label, " *
        "f_a = $(scientific(decay_constant_gev)) GeV (f_a/M_Pl = $(scientific(decay_constant_ratio))), equal initial occupations"
Label(figure[0, compact ? 1 : (1:2)], title; fontsize=compact ? 20 : 17, color=ink, halign=:left, tellwidth=false)

shared_end = minimum(run.time[end] for run in runs)
longest_end = maximum(run.time[end] for run in runs)
# the log-time column starts at 0.2 × the shortest run's end: earlier accepted states are sparse (straight segments on
# a log axis) and that span is already in the left column
windows = ((0.0, shared_end, false, "All runs, t ≤ $(scientific(shared_end)) M₀"),
           (0.2shared_end, longest_end, true, "To $(scientific(longest_end)) M₀ (log time); dotted: final value of a shorter run"))
compact && (windows = ((windows[2][1], windows[2][2], true, "Log time; dotted: final value of the shorter run"),))
for (column_index, (time_min, time_max, log_axis, window_title)) in enumerate(windows)
    first_column = column_index == 1
    time_scale = log_axis ? log10 : identity
    # log axis: ticks at 1, 2, 5 × 10ⁿ inside the window (Makie's default gives fractional exponents like 10^8.7)
    time_ticks = Makie.automatic
    if log_axis
        tick_values = [m * 10.0^n for n in floor(Int, log10(time_min)):ceil(Int, log10(time_max)) for m in (1, 2, 5)]
        tick_values = filter(t -> time_min ≤ t ≤ time_max, tick_values)
        time_ticks = (tick_values, scientific.(tick_values))
    end
    axis_occupation = Axis(figure[1, column_index]; yscale=log10, xscale=time_scale, xticks=time_ticks, title=window_title, titlealign=:left,
                           ylabel=first_column ? L"occupation $\mathcal{N}_i$  [$f_a^2 M_0^2$]" : "", xgridcolor=grid, ygridcolor=grid)
    axis_spin = Axis(figure[2, column_index]; xscale=time_scale, xticks=time_ticks, ylabel=first_column ? L"black-hole spin $a = J/M^2$" : "",
                     xlabel=compact ? L"t/M_0" : "", xgridcolor=grid, ygridcolor=grid)
    # compact: no mass row (the mass calls below go to an axis that is never placed in the figure)
    axis_mass = compact ? Axis(Figure()[1, 1]) :
                Axis(figure[3, column_index]; xscale=time_scale, xticks=time_ticks, xlabel=L"t/M_0",
                     ylabel=first_column ? L"black-hole mass $M/M_0$" : "", xgridcolor=grid, ygridcolor=grid)
    time_axes = compact ? [axis_occupation, axis_spin] : [axis_occupation, axis_spin, axis_mass]
    for run_index in eachindex(runs)
        run, style = runs[run_index], line_styles[run_index]
        shown = time_min .≤ run.time .≤ time_max
        if style == :markers   # about 40 open circles per window, spaced evenly in the plotted time coordinate
            shown_indices = findall(shown)
            coordinate = log_axis ? log10.(max.(run.time[shown_indices], time_min)) : run.time[shown_indices]
            targets = range(coordinate[1], coordinate[end]; length=40)
            picked = unique([shown_indices[argmin(abs.(coordinate .- target))] for target in targets])
            marker_style = (marker=:circle, markersize=7, strokewidth=1.2, color=:transparent)
            scatter!(axis_occupation, run.time[picked], run.occupation_211[picked]; strokecolor=mode_211_color, marker_style...)
            scatter!(axis_occupation, run.time[picked], run.occupation_322[picked]; strokecolor=mode_322_color, marker_style...)
            scatter!(axis_spin, run.time[picked], run.threshold[picked]; strokecolor=(mode_211_color, 0.8), marker_style...)
            scatter!(axis_spin, run.time[picked], run.spin[picked]; strokecolor=ink, marker_style...)
            scatter!(axis_mass, run.time[picked], run.mass[picked]; strokecolor=ink, marker_style...)
            continue
        end
        plot_times, plot_211, plot_322, plot_threshold, plot_spin, plot_mass = log_axis ?
            log_axis_resample(run.time, run.occupation_211, run.occupation_322, run.threshold, run.spin, run.mass;
                              logarithmic=(true, true, false, false, false), time_min, time_max) :
            (run.time[shown], run.occupation_211[shown], run.occupation_322[shown], run.threshold[shown], run.spin[shown],
             run.mass[shown])
        lines!(axis_occupation, plot_times, plot_211; color=mode_211_color, linewidth=line_width, linestyle=style)
        lines!(axis_occupation, plot_times, plot_322; color=mode_322_color, linewidth=line_width, linestyle=style)
        lines!(axis_spin, plot_times, plot_threshold; color=(mode_211_color, 0.8), linewidth=0.75line_width, linestyle=style)
        lines!(axis_spin, plot_times, plot_spin; color=ink, linewidth=line_width, linestyle=style)
        lines!(axis_mass, plot_times, plot_mass; color=ink, linewidth=line_width, linestyle=style)
        if run.time[end] < time_max   # continue a run that ended early at its final value
            continuation = [run.time[end], time_max]
            lines!(axis_occupation, continuation, fill(run.occupation_211[end], 2); color=mode_211_color, linestyle=:dot, linewidth=1.2)
            lines!(axis_occupation, continuation, fill(run.occupation_322[end], 2); color=mode_322_color, linestyle=:dot, linewidth=1.2)
            lines!(axis_spin, continuation, fill(run.spin[end], 2); color=ink, linestyle=:dot, linewidth=1.2)
            lines!(axis_mass, continuation, fill(run.mass[end], 2); color=ink, linestyle=:dot, linewidth=1.2)
        end
    end
    for axis in time_axes
        xlims!(axis, time_min, time_max)
    end
    linkxaxes!(time_axes...)
    for axis in time_axes[1:end-1]   # tick labels only on the bottom axis
        hidexdecorations!(axis; grid=false, ticks=false)
    end
    if first_column
        level_entries = [LineElement(color=mode_211_color, linewidth=2), LineElement(color=mode_322_color, linewidth=2)]
        run_entries = [line_styles[index] == :markers ?
                       MarkerElement(marker=:circle, markersize=7, strokewidth=1.2, color=:transparent, strokecolor=ink) :
                       LineElement(color=ink, linewidth=2, linestyle=line_styles[index]) for index in eachindex(runs)]
        axislegend(axis_occupation, [level_entries, run_entries], [[L"|211\rangle", L"|322\rangle"], collect(labels)],
                   ["level", "run"]; position=:rc, framevisible=false, orientation=:horizontal, nbanks=2)
        axislegend(axis_spin, [LineElement(color=ink, linewidth=2), LineElement(color=(mode_211_color, 0.8), linewidth=1.5)],
                   ["spin a", L"|211\rangle threshold ($\omega_{211} = \Omega_H$)"]; position=:rt, framevisible=false)
    end
end
rowsize!(figure.layout, 1, Relative(compact ? 0.55 : 0.42))
compact || colgap!(figure.layout, 40)
mkpath(dirname(output_path))
save(output_path, figure)
println("saved $output_path")
