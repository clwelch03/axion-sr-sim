# Mode growth and black-hole spin-down along a finished spindown evolution, from the diagnostics CSV written by
# scripts/evolution_spindown.jl (one row per rhs evaluation). Rows: occupations N_i of |211⟩ and |322⟩ (reference
# units f_a² M₀²); black-hole spin a = j/b² with the |211⟩ superradiance threshold; black-hole mass b = M/M₀.
# Left: the whole run; right: a window around the |211⟩ peak.
#
# The diagnostics log every rhs evaluation: trial stages, including those of rejected steps (error-estimate rejections
# and failed radial solves). Only accepted states are plotted. Tsit5 evaluates its last two stages at the same time
# t + dt, and the second of that pair is the accepted state (first-same-as-last); a rejected attempt is followed by a
# retry at earlier times. So a row is kept if it is the first row (the initial state), or if it repeats the previous
# row's time and no later row has an earlier time. Intermediate stages are not on the trajectory: near the |211⟩
# threshold they scatter by ~3e-7 in b, which showed up as spurious jumps when they were plotted.
#
# |211⟩ threshold: superradiance stops when ω/m = Ω_H, the horizon angular velocity. With Ω_H = a/(2r₊),
# r₊ = 1 + √(1 − a²) (units of the current mass), the inverse is a = 4Ω/(1 + 4Ω²). ω is Re ω₂₁₁ in current units
# (as recorded), so the threshold moves as α = α₀ b and the cloud's frequency shift change.
#
# In-progress runs: pass the Slurm log (results/<prefix>-<jobid>.out) instead of a prefix. The verbose log has t, A, Im ω,
# b, j per rhs evaluation but not 𝒩 or Re ω, so the top row shows the amplitudes A_i, and the |211⟩ threshold line is
# replaced by shading where Im ω₂₁₁ < 0. That is net |211⟩ decay, not the threshold: in the nonlinear two-mode regime
# Im ω₂₁₁ includes the transfer to |322⟩ (e.g. a₀ = 0.99, f_a/M_Pl = 3e-3: Im ω₂₁₁ < 0 at a ≈ 0.98, threshold ≈ 0.877).
#
# Run from the package root:
#   julia --project=scripts scripts/plot_spindown.jl <prefix | log.out> [output.png]
# e.g. spindown_a0.9_alpha0.3_fa3e-3  → reads results/<prefix>_diagnostics.csv, writes plots/<prefix>.png
#      results/spindown_a0.99_alpha0.3_fa3e-3-61102850.out → writes plots/spindown_a0.99_alpha0.3_fa3e-3_partial.png
using CairoMakie, DelimitedFiles, Printf

from_log = endswith(ARGS[1], ".out")
prefix = from_log ? replace(basename(ARGS[1]), r"-[0-9]+\.out$" => "") : ARGS[1]
output_path = length(ARGS) ≥ 2 ? ARGS[2] :
              joinpath(@__DIR__, "..", "plots", prefix * (from_log ? "_partial" : "") * ".png")
planck_mass_gev = 1.2209e19

if from_log
    number = raw"([-+0-9.eE]+)"
    log_line = Regex("^t = $number\\s+A = \\[$number, $number\\]\\s+Im ω = \\[$number, $number\\].*b = $number\\s+j = $number")
    log_rows = [parse.(Float64, captured.captures) for captured in filter(!isnothing, match.(log_line, eachline(ARGS[1])))]
    log_names = ["time", "A1", "A2", "Im_omega1", "Im_omega2", "b", "j"]
    column(name) = name == "a" ? column("j") ./ column("b") .^ 2 : [row[findfirst(==(name), log_names)] for row in log_rows]
else
    table, header = readdlm(joinpath(@__DIR__, "..", "results", prefix * "_diagnostics.csv"), ','; header=true)
    column(name) = table[:, findfirst(==(name), vec(header))]
end

# keep accepted states only (see the header)
raw_times = column("time")
earliest_later_time = similar(raw_times)
running_minimum = Inf
for index in reverse(eachindex(raw_times))
    earliest_later_time[index] = running_minimum
    global running_minimum = min(running_minimum, raw_times[index])
end
not_rejected = raw_times .≤ earliest_later_time
accepted_state = [index == 1 || (not_rejected[index] && raw_times[index] == raw_times[index - 1]) for index in eachindex(raw_times)]
series(name) = column(name)[accepted_state]
times = raw_times[accepted_state]
issorted(times) || error("accepted states are not in time order; the stage-pair assumption (Tsit5) does not hold")
# from a log, "occupation" holds the amplitudes A_i (labelled as such below)
occupation_211, occupation_322 = from_log ? (series("A1"), series("A2")) : (series("N1"), series("N2"))
black_hole_mass, black_hole_spin = series("b"), series("a")
if from_log
    below_threshold_211 = series("Im_omega1") .< 0
else
    frequency_211 = series("Re_omega1")                      # current-mass units, m = 1
    threshold_spin_211 = @. 4frequency_211 / (1 + 4frequency_211^2)
end

# case parameters from the file name
spin_label = match(r"_a([0-9.]+)_", prefix)[1]
mass_coupling_label = match(r"_alpha([0-9.]+)_", prefix)[1]
decay_constant_ratio = parse(Float64, match(r"_fa([0-9.eE+-]+?)(_t|$)", prefix)[1])

"""`x` as "3.7 × 10¹⁶" (or "10¹⁶" when the mantissa is 1)."""
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
grid = (:black, 0.06)
figure = Figure(size=(1300, 1000), fontsize=15)
title = "Spin-down: a₀ = $spin_label, α₀ = $mass_coupling_label, f_a = $(scientific(decay_constant_ratio * planck_mass_gev)) GeV " *
        "(f_a/M_Pl = $(scientific(decay_constant_ratio))), $(occursin("spindown_free_", prefix) ? "free" : "cosine") potential, equal initial occupations" *
        (from_log ? "\nin progress: Slurm log to t = $(scientific(times[end]; digits=3)) M₀ (amplitudes, not occupations)" : "")
Label(figure[0, 1:2], title; fontsize=17, color=ink, halign=:left, tellwidth=false)

# peak: first time within 0.1% of the maximum, so a run that saturates (no peak) centres on the start of the plateau
# instead of on a 1e-7 wiggle along it
peak_index = findfirst(≥((1 - 1e-3) * maximum(occupation_211)), occupation_211)
span = times[end] - times[1]
windows = ((times[1], times[end], "Whole run"),
           (max(times[1], times[peak_index] - 0.08span), min(times[end], times[peak_index] + 0.12span), "Around the |211⟩ peak"))
for (column_index, (time_min, time_max, window_title)) in enumerate(windows)
    shown = time_min .≤ times .≤ time_max
    first_column = column_index == 1
    axis_occupation = Axis(figure[1, column_index]; yscale=log10, title=window_title, titlealign=:left,
                           ylabel=!first_column ? "" : from_log ? L"amplitude $A_i$" : L"occupation $\mathcal{N}_i$  [$f_a^2 M_0^2$]",
                           xgridcolor=grid, ygridcolor=grid)
    lines!(axis_occupation, times[shown], occupation_211[shown]; color=mode_211_color, linewidth=2, label=L"|211\rangle")
    lines!(axis_occupation, times[shown], occupation_322[shown]; color=mode_322_color, linewidth=2, label=L"|322\rangle")

    axis_spin = Axis(figure[2, column_index]; ylabel=first_column ? L"black-hole spin $a = J/M^2$" : "",
                     xgridcolor=grid, ygridcolor=grid)
    lines!(axis_spin, times[shown], black_hole_spin[shown]; color=ink, linewidth=2, label="spin a")
    if from_log
        # shade accepted-step stages below the |211⟩ threshold (Im ω₂₁₁ < 0), each to the next stage
        shaded_times = times[shown]
        shaded = findall(below_threshold_211[shown][1:end-1])
        isempty(shaded) || vspan!(axis_spin, shaded_times[shaded], shaded_times[shaded .+ 1]; color=(mode_211_color, 0.25),
                                  label=L"|211\rangle decaying (Im $\omega_{211} < 0$)")
    else
        lines!(axis_spin, times[shown], threshold_spin_211[shown]; color=mode_211_color, linewidth=2, linestyle=:dash,
               label=L"|211\rangle threshold ($\omega_{211} = \Omega_H$)")
    end

    axis_mass = Axis(figure[3, column_index]; xlabel=L"t/M_0", ylabel=first_column ? L"black-hole mass $M/M_0$" : "",
                     xgridcolor=grid, ygridcolor=grid)
    lines!(axis_mass, times[shown], black_hole_mass[shown]; color=ink, linewidth=2)

    for axis in (axis_occupation, axis_spin, axis_mass)
        xlims!(axis, time_min, time_max)
    end
    linkxaxes!(axis_occupation, axis_spin, axis_mass)
    hidexdecorations!(axis_occupation; grid=false, ticks=false)
    hidexdecorations!(axis_spin; grid=false, ticks=false)
    if first_column
        axislegend(axis_occupation; position=:rc, framevisible=false)   # right-centre is empty in every run so far
        axislegend(axis_spin; position=:rt, framevisible=false)
    end
end
rowsize!(figure.layout, 1, Relative(0.45))
colgap!(figure.layout, 40)
mkpath(dirname(output_path))
save(output_path, figure)
println("saved $output_path ($(length(times)) accepted states of $(length(raw_times)) rhs evaluations)")
if from_log
    @printf("peak A₂₁₁ = %.4e at t = %.4e; last: t = %.4e, A = (%.4e, %.4e), a = %.5f (Im ω₂₁₁ %s 0), M/M₀ = %.6f\n",
            occupation_211[peak_index], times[peak_index], times[end], occupation_211[end], occupation_322[end],
            black_hole_spin[end], below_threshold_211[end] ? "<" : "≥", black_hole_mass[end])
else
    @printf("peak 𝒩₂₁₁ = %.4e at t = %.4e; final: 𝒩 = (%.4e, %.4e), a = %.5f (|211⟩ threshold %.5f), M/M₀ = %.6f\n",
            occupation_211[peak_index], times[peak_index], occupation_211[end], occupation_322[end], black_hole_spin[end],
            threshold_spin_211[end], black_hole_mass[end])
end
