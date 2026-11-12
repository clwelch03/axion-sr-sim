# Our plateau spin-down at α = 0.187 on the axes of Baryakhtar, Galanis, Lasenby & Simon, arXiv:2011.11646, Fig. 9 right
# panel (10 M☉, μ = 2.5e-12 eV), with their curves digitized by scripts/digitize_baryakhtar_fig9.py.
#
# Ours: the quasi-static plateau spin-down of scripts/plateau_spindown_time.jl, evaluated densely in a. The plateau is
# taken to form after the free |211⟩ growth from one quantum to the plateau occupation, t_growth = ln N/(2 ω_I) with
# N = κ 𝒩₂₁₁ (M/M_Pl)² (𝒩₂₁₁ in f_a² M² units at a = 0.9, ω_I the free rate at a = 0.9), and the plateau spin-down
# starts then:
#     t(a) = t_growth + (M/κ) ∫_a^{0.9} da′ / D(a′),   D = F_J − 2a F_E (horizon fluxes, f_a² units),
# with log D linear in a between scan spins. M and α are held fixed. Our curves stop at the last scan spin, a = 0.66,
# just above the |211⟩ threshold (D → 0 there). Only the four smallest f are drawn: larger f needs the full evolution.
# Reference lines (ours, as in their figure): |211⟩ and |322⟩ superradiance boundaries a = 4Ω/(1 + 4Ω²), Ω = Re ω/m,
# from the scan's frequencies at a = 0.66; free spin-down times ln(GM²/ħc)/(2 Im ω) with spheroidal Leaver Im ω
# (main:src/common/leaver.py), |211⟩ at a = 0.9 and |322⟩ at a = 0.654.
#
# Run from the package root:
#   julia --project=scripts scripts/plot_baryakhtar_fig9_comparison.jl
# Reads results/plateau_scan_alpha0.187.csv and results/baryakhtar_fig9_right_digitized.csv; writes
# plots/baryakhtar_fig9_comparison.png and prints t(a = 0.78) and t(a = 0.70), ours vs theirs.
using CairoMakie, DelimitedFiles, Printf

results_directory = joinpath(@__DIR__, "..", "results")
planck_mass_gev = 1.2209e19                 # 1/√G
hole_mass_solar = 10.0
mass_unit_years = hole_mass_solar * 4.925491e-6 / 3.15576e7
log_hole_mass_squared = 2log(hole_mass_solar * 1.98892e30 / 2.176434e-8)   # ln(G M²/ħc) = 179.4
decay_constants_gev = [1e11, 5e12, 5e14, 5e15]
# free spin-down times (yr) from spheroidal Leaver at α = 0.187: Im ω₂₁₁(a = 0.9) = 9.41679e-10,
# Im ω₃₂₂(a = 0.654) = 3.33120e-15 (see docs/session-2026-11-10.md, 2026-11-11 entries)
free_time_211 = log_hole_mass_squared / (2 * 9.41679e-10) * mass_unit_years
free_time_322 = log_hole_mass_squared / (2 * 3.33120e-15) * mass_unit_years

# plateau scan
table, header = readdlm(joinpath(results_directory, "plateau_scan_alpha0.187.csv"), ','; header=true)
column(name) = Float64.(table[:, findfirst(==(name), vec(header))])
converged = column("converged") .== 1
order = sortperm(column("spin")[converged]; rev=true)
pick(name) = column(name)[converged][order]
spins = pick("spin")
spin_rate = pick("inner_angular_momentum_flux") .- 2 .* spins .* pick("inner_energy_flux")   # D(a)
growth_rate_211 = pick("growth_rate_211")[1]
plateau_occupation_211 = pick("N1")[1]
boundary(frequency, m) = (Ω = frequency / m; 4Ω / (1 + 4Ω^2))
boundary_211, boundary_322 = boundary(pick("Re_omega1")[end], 1), boundary(pick("Re_omega2")[end], 2)

"""∫_{a_low}^{a_high} da/D with D = D_high exp(s (a − a_high)), s from the bracketing scan spins (as in plateau_spindown_time.jl)."""
function partial_integral(a_high, a_low, rate_high, slope)
    abs(slope * (a_low - a_high)) < 1e-10 && return (a_high - a_low) / rate_high
    (exp(-slope * (a_low - a_high)) - 1) / (slope * rate_high)
end
"""κ t / M from a = spins[1] down to `spin` (≥ spins[end])."""
function scaled_time(spin)
    total = 0.0
    for index in 1:length(spins)-1
        a_high, a_low = spins[index], spins[index + 1]
        slope = log(spin_rate[index + 1] / spin_rate[index]) / (a_low - a_high)
        if spin ≥ a_low
            return total + partial_integral(a_high, spin, spin_rate[index], slope)
        end
        total += partial_integral(a_high, a_low, spin_rate[index], slope)
    end
    error("spin $spin below the last scan spin $(spins[end])")
end
dense_spins = collect(range(spins[1], spins[end]; length=600))
dense_scaled_times = scaled_time.(dense_spins)
coupling(decay_constant) = (decay_constant / planck_mass_gev)^2
growth_time(decay_constant) =
    log(coupling(decay_constant) * plateau_occupation_211) / (2growth_rate_211) * mass_unit_years +
    log_hole_mass_squared / (2growth_rate_211) * mass_unit_years
our_times(decay_constant) = growth_time(decay_constant) .+ dense_scaled_times ./ coupling(decay_constant) .* mass_unit_years

# digitized curves
digitized, _ = readdlm(joinpath(results_directory, "baryakhtar_fig9_right_digitized.csv"), ','; header=true)
their_curve(decay_constant) = (rows = findall(isapprox.(digitized[:, 1], decay_constant; rtol=1e-6));
                               (10.0 .^ digitized[rows, 2], Float64.(digitized[rows, 3])))
"""Time at which a curve (times, spins; spins decreasing) crosses `level`, linear in log t."""
function crossing_time(times, curve_spins, level)
    index = findfirst(i -> curve_spins[i] ≥ level > curve_spins[i + 1], 1:length(curve_spins)-1)
    index === nothing && return NaN
    fraction = (curve_spins[index] - level) / (curve_spins[index] - curve_spins[index + 1])
    exp10((1 - fraction) * log10(times[index]) + fraction * log10(times[index + 1]))
end

@printf("free spin-down times: |211⟩ %.3e yr, |322⟩ %.3e yr; boundaries: |211⟩ a = %.4f, |322⟩ a = %.4f\n",
        free_time_211, free_time_322, boundary_211, boundary_322)
@printf("%10s %12s %12s %12s %8s %12s %12s %8s\n", "f [GeV]", "t_growth", "ours(0.78)", "theirs(0.78)", "ratio",
        "ours(0.70)", "theirs(0.70)", "ratio")
for decay_constant in decay_constants_gev
    their_times, their_spins = their_curve(decay_constant)
    ours = [crossing_time(our_times(decay_constant), dense_spins, level) for level in (0.78, 0.70)]
    theirs = [crossing_time(their_times, their_spins, level) for level in (0.78, 0.70)]
    @printf("%10.0e %12.3e %12.3e %12.3e %8.2f %12.3e %12.3e %8.2f\n", decay_constant, growth_time(decay_constant),
            ours[1], theirs[1], ours[1] / theirs[1], ours[2], theirs[2], ours[2] / theirs[2])
end

# figure: repo palette (main: src/common/plotstyle.py), one colour per f
colors = [colorant"#B14FC5", colorant"#ED40A9", colorant"#F38B00", colorant"#00A0DF"]
ink = colorant"#0b0b0b"
superscript_digits = Dict('-' => '⁻', '0' => '⁰', '1' => '¹', '2' => '²', '3' => '³', '4' => '⁴', '5' => '⁵', '6' => '⁶',
                          '7' => '⁷', '8' => '⁸', '9' => '⁹')
function power_label(x)
    exponent = floor(Int, log10(x)); mantissa = round(Int, x / 10.0^exponent)
    power = "10" * join(superscript_digits[c] for c in string(exponent))
    mantissa == 1 ? power : "$mantissa × $power"
end
figure = Figure(size=(1100, 760), fontsize=20, figure_padding=(16, 40, 16, 16))
Label(figure[0, 1], "Plateau spin-down vs Baryakhtar et al. Fig. 9 (right): 10 M☉, μ = 2.5 × 10⁻¹² eV (α = 0.187)";
      fontsize=20, halign=:left, tellwidth=false)
tick_values = [10.0^n for n in -2:2:10]
axis = Axis(figure[1, 1]; xscale=log10, xlabel="t (yr)", ylabel=L"black-hole spin $a$",
            xticks=(tick_values, power_label.(tick_values)), xgridcolor=(:black, 0.06), ygridcolor=(:black, 0.06))
hlines!(axis, [boundary_211, boundary_322]; color=ink, linestyle=:dash, linewidth=1.5)
vlines!(axis, [free_time_211, free_time_322]; color=(:gray50, 0.9), linestyle=:dash, linewidth=1.5)
for (index, decay_constant) in enumerate(decay_constants_gev)
    their_times, their_spins = their_curve(decay_constant)
    lines!(axis, their_times, their_spins; color=(colors[index], 0.3), linewidth=10)
    times = our_times(decay_constant)
    lines!(axis, vcat(1e-2, times), vcat(spins[1], dense_spins); color=colors[index], linewidth=2.5)
    scatter!(axis, [times[end]], [dense_spins[end]]; color=colors[index], markersize=9)
end
xlims!(axis, 1e-2, 2e10)
ylims!(axis, 0.25, 0.95)
style_entries = [LineElement(color=(ink, 0.3), linewidth=10), LineElement(color=ink, linewidth=2.5)]
f_entries = [LineElement(color=color, linewidth=3) for color in colors]
axislegend(axis, [f_entries, style_entries],
           [["f = $(power_label(f)) GeV" for f in decay_constants_gev], ["Baryakhtar et al. (digitized)", "ours (plateau)"]],
           ["", ""]; position=(0.97, 0.37), framevisible=false, nbanks=1, rowgap=2)   # empty band a ≈ 0.40-0.62, right
text!(axis, 1.3e-2, boundary_211 + 0.008; text="|211⟩ boundary", fontsize=16, color=ink)
text!(axis, 1.3e-2, boundary_322 + 0.008; text="|322⟩ boundary", fontsize=16, color=ink)
output_path = joinpath(@__DIR__, "..", "plots", "baryakhtar_fig9_comparison.png")
save(output_path, figure)
println("saved $output_path")
