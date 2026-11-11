# Plot the joint-deformation comparison written by scripts/profile_comparison_regimes.jl (layout of `main`'s
# runs/adiabatic/multiphase_two_cloud/scripts/plot_collocation_comparison.py), with the axion's physical parameters
# in the title. Upper panel: energy per unit tortoise coordinate dE/dr* of each cloud's fundamental component (free
# field, in f_a²; every curve of one cloud encloses the same E). Lower panels: change in the enclosed-energy fraction
# E(<r*)/E against the free-field cloud (V = ½μ²φ², no self-interaction; `linear` in the CSV column names), i.e. the
# share of the cloud's energy that the nonlinearity has moved inside r*.
# The profiles are in f_a-scaled units, so f_a and m_a only enter the caption:
#   M = α/(G m_a) (α = 0.3);  cloud mass E (f_a/M_Pl)² M;  spin taken Δa = J (f_a/M_Pl)²  (M_Pl = 1.2209e19 GeV).
#
# Run from the package root:
#   julia --project=scripts scripts/plot_profile_comparison_regimes.jl <label> "<axion name>" <f_a in GeV> <m_a in eV>
# e.g. qcd "QCD axion" 5.7e17 1e-11   |   stringy "String axion" 1e16 1e-12
# Reads results/profile_comparison_<label>.csv, writes plots/profile_comparison_<label>.png.
using CairoMakie, DelimitedFiles, Printf

label, axion_name = ARGS[1], ARGS[2]
decay_constant_gev, axion_mass_ev = parse(Float64, ARGS[3]), parse(Float64, ARGS[4])
planck_mass_gev = 1.2209e19

results_path = joinpath(@__DIR__, "..", "results", "profile_comparison_$(label).csv")
header_lines = filter(startswith("#"), readlines(results_path))
settings = Dict(m[1] => m[2] for m in eachmatch(r"(\w+) = ([^,\s]+)", header_lines[2]))
number(key) = parse(Float64, settings[key])
spin, mass_coupling, energy = number("spin"), number("alpha"), number("energy")
angular_momentum, threshold_spin, potential_name = number("angular_momentum"), number("threshold_spin"), settings["potential"]
table, header = readdlm(results_path, ','; header=true, comments=true)
column(name) = table[:, findfirst(==(name), vec(header))]
rstar = column("r_star_over_M")

"""`x` as "5.7 × 10¹⁷" (or "10⁻¹¹" when the mantissa is 1)."""
function scientific(x; digits=2)
    exponent = floor(Int, log10(abs(x)))
    mantissa = round(x / 10.0^exponent; sigdigits=digits)
    mantissa ≥ 10 && (mantissa /= 10; exponent += 1)
    superscript = join(Dict('-' => '⁻', '0' => '⁰', '1' => '¹', '2' => '²', '3' => '³', '4' => '⁴', '5' => '⁵',
                            '6' => '⁶', '7' => '⁷', '8' => '⁸', '9' => '⁹')[c] for c in string(exponent))
    exponent == 0 && return string(mantissa)
    mantissa == 1 ? "10" * superscript : @sprintf("%g × 10%s", mantissa, superscript)
end

ratio_squared = (decay_constant_gev / planck_mass_gev)^2
black_hole_mass = mass_coupling / (7.4803e9 * axion_mass_ev)          # M/M☉ from α = G M m_a/(ħ c³)
spin_taken = angular_momentum * ratio_squared
spin_room = spin - threshold_spin
spin_taken > spin_room && @warn "these clouds need more spin than the hole has: Δa = $spin_taken > $spin_room"
title = string(@sprintf("%s:  fₐ = %s GeV,  mₐ = %s eV,  M ≈ %.2g M☉", axion_name, scientific(decay_constant_gev),
                        scientific(axion_mass_ev), black_hole_mass), "\n",
               @sprintf("α = %g,  a/M = %g,  %s potential,  E₁ = E₂ = %g fₐ²M  (each cloud %s M)", mass_coupling, spin,
                        potential_name, energy, scientific(energy * ratio_squared)), "\n",
               @sprintf("the two clouds have taken Δa = %s of the %.3f spin above the |211⟩ threshold",
                        spin_taken ≥ 0.01 ? @sprintf("%.3f", spin_taken) : @sprintf("%.2g", spin_taken), spin_room))

# repo palette (main: src/common/plotstyle.py): brand colours in cycle order, C0 green |211⟩, C1 blue |322⟩
colors = (primary=colorant"#3DAE2B", secondary=colorant"#00A0DF")
ink = colorant"#0b0b0b"
styles = (linear=:dot, omiya=:dash, joint=:solid)
# `linear` is the free field. It is labelled "Free field" so it is not confused with Omiya's secondary, which is
# linear in its own amplitude but sits on the self-interacting primary.
model_names = (linear="Free field", omiya="Omiya", joint="Joint")
# The r² in dE/dr* leaves almost nothing near the horizon (|dE/dr*| ~ 3e-6 of the peak there), and the |322⟩ cloud
# keeps 5% of its energy beyond r* = 200M, so the window is shifted outward from the |R| plot's [-100, 200].
plot_range = (-20, 350)
domain = (rstar .>= plot_range[1]) .& (rstar .<= plot_range[2])
figure = Figure(size=(1300, 1150), fontsize=22)
Label(figure[0, 1:2], title; fontsize=23, color=ink, tellwidth=false, lineheight=1.15)
axis_options = (xlabel=L"r_*/M", xlabelsize=26, ylabelsize=26, xticklabelsize=21, yticklabelsize=21,
                xticksmirrored=true, yticksmirrored=true)
upper = Axis(figure[1, 1:2]; ylabel=L"dE/dr_* \; [f_a^2]", axis_options...)
lower = [Axis(figure[2, index]; ylabel="Change in E(<r*)/E [%]\nrelative to free field", title=name, titlesize=24, axis_options...)
         for (index, name) in enumerate(("Primary |211⟩", "Secondary |322⟩"))]

"""Enclosed-energy fraction E(<r*)/E from the left numerical endpoint (trapezoid rule on the full mesh)."""
function enclosed_fraction(density)
    cumulative = [0.0; cumsum(diff(rstar) .* (density[1:end-1] .+ density[2:end]) ./ 2)]
    return cumulative ./ cumulative[end]
end

for (index, cloud) in enumerate((:primary, :secondary))
    color = colors[cloud]
    linear_fraction = enclosed_fraction(column("energy_density_linear_$(cloud)"))
    difference = Dict{Symbol, Vector{Float64}}()
    for model in (:linear, :omiya, :joint)
        density = column("energy_density_$(model)_$(cloud)")
        lines!(upper, rstar[domain], density[domain]; color, linestyle=styles[model], linewidth=3,
               label=@sprintf("%s: %s", titlecase(string(cloud)), model_names[model]))
        difference[model] = 100 .* (enclosed_fraction(density) .- linear_fraction)[domain]
        model == :linear && continue
        lines!(lower[index], rstar[domain], difference[model]; color, linestyle=styles[model], linewidth=2.6,
               label=model_names[model])
    end
    lines!(lower[index], rstar[domain], difference[:joint] .- difference[:omiya]; color=ink, linewidth=2,
           label="Joint − Omiya")
    hlines!(lower[index], 0; color=RGBf(0.6, 0.6, 0.6), linestyle=:dot, linewidth=1)
    axislegend(lower[index]; framevisible=false, labelsize=19, patchsize=(35, 12))
end
axislegend(upper; position=:rt, nbanks=2, framevisible=false, labelsize=20, patchsize=(40, 12))
ylims!(upper; low=0)
for axis in (upper, lower...)
    xlims!(axis, plot_range...)
end
rowsize!(figure.layout, 1, Relative(0.55))
plot_path = joinpath(@__DIR__, "..", "plots", "profile_comparison_$(label).png")
save(plot_path, figure; px_per_unit=2)
println(plot_path)
