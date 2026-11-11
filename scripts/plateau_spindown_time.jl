# Black-hole spin-down time while the |211⟩ + |322⟩ cloud sits at its quasi-equilibrium (plateau), from the plateau
# scan's horizon fluxes (scripts/plateau_scan.jl, columns inner_energy_flux and inner_angular_momentum_flux).
#
# Adiabatic tracking: at small f_a the plateau relaxes on the superradiance time, while the hole changes ~(M_Pl/f_a)²
# times more slowly, so the cloud is taken to sit at the plateau of the current spin. In f_a-scaled units the plateau
# does not depend on f_a. With κ = (f_a/M_Pl)² (M_Pl = 1/√G, the convention of the evolution code and of Baryakhtar et
# al. 2011.11646), current-mass units (G = c = 1, time in M), and F_E, F_J the outward fluxes at the inner node in f_a²
# units:
#     dM/dt = −κ F_E,   dJ/dt = −κ M F_J,   a = J/M²   ⇒   da/dt = −κ (F_J − 2a F_E) ≡ −κ D(a)
# (D(a) is my notation: the spin-down rate per unit κ). The hole's mass change (≈ 1 %) and the resulting drift of
# α = μM are neglected: M stays M₀ and α stays the scan's value. Then
#     t(a) = (M₀/κ) ∫_a^{a₀} da′ / D(a′),
# so every small-f_a curve is one curve shifted along log t by (M_Pl/f_a)². D is interpolated between scan spins
# linearly in log D. D → 0 at the |211⟩ threshold, so t → ∞ there: times are reported at the scan spins only.
#
# Run from the package root:
#   julia --project=. scripts/plateau_spindown_time.jl <plateau csv> [M in M☉ = 10] [f_a in GeV ...]
using DelimitedFiles, Printf

csv_path = ARGS[1]
hole_mass_solar = length(ARGS) ≥ 2 ? parse(Float64, ARGS[2]) : 10.0
decay_constants_gev = length(ARGS) ≥ 3 ? parse.(Float64, ARGS[3:end]) : [1e11, 5e12, 5e14]

planck_mass_gev = 1.2209e19                  # 1/√G
solar_mass_seconds = 4.925491e-6             # G M☉ / c³
year_seconds = 3.15576e7
mass_unit_years = hole_mass_solar * solar_mass_seconds / year_seconds   # one M of time, in years

table, header = readdlm(csv_path, ','; header=true)
header = vec(header)
column(name) = Float64.(table[:, findfirst(==(name), header)])
converged = column("converged") .== 1
spins = column("spin")[converged]
order = sortperm(spins; rev=true)
spins = spins[order]
pick(name) = column(name)[converged][order]
energy_211, energy_322 = pick("E211"), pick("E322")
inner_energy_flux, inner_angular_momentum_flux = pick("inner_energy_flux"), pick("inner_angular_momentum_flux")
growth_rate_211 = pick("growth_rate_211")
spin_rate = inner_angular_momentum_flux .- 2 .* spins .* inner_energy_flux       # D(a)
all(spin_rate .> 0) || error("D(a) ≤ 0 at some scan spin; the plateau does not spin the hole down there")

# ∫ da / D between neighbouring scan spins with log D linear in a: for D = D₁ exp(s (a − a₁)), s = ln(D₂/D₁)/(a₂ − a₁),
# ∫_{a₂}^{a₁} da / D = (exp(−s (a₂ − a₁)) − 1) / (s D₁) (→ (a₁ − a₂)/D₁ as s → 0).
function interval_integral(a1, a2, rate1, rate2)
    slope = log(rate2 / rate1) / (a2 - a1)
    abs(slope * (a2 - a1)) < 1e-10 && return (a1 - a2) / rate1
    return (exp(-slope * (a2 - a1)) - 1) / (slope * rate1)
end
scaled_time = zeros(length(spins))           # κ t / M₀
for index in 2:length(spins)
    scaled_time[index] = scaled_time[index - 1] +
        interval_integral(spins[index - 1], spins[index], spin_rate[index - 1], spin_rate[index])
end

alpha = column("alpha")[1]
@printf("α = %.4g, M = %g M☉ (1 M = %.4e yr); plateau spin-down from a₀ = %.3f\n", alpha, hole_mass_solar,
        mass_unit_years, spins[1])
@printf("free-field |211⟩ e-fold time of N at a₀: 1/(2 ω_I) = %.3e M = %.3e yr\n",
        1 / (2growth_rate_211[1]), mass_unit_years / (2growth_rate_211[1]))
println()
@printf("%6s %10s %10s %11s %11s %11s %11s", "a", "E211", "E322", "F_E,inner", "F_J,inner", "D(a)", "κt/M₀")
for decay_constant in decay_constants_gev
    @printf(" %12s", @sprintf("t[yr] %.0e", decay_constant))
end
println()
for index in eachindex(spins)
    @printf("%6.3f %10.4g %10.4g %11.4e %11.4e %11.4e %11.4e", spins[index], energy_211[index], energy_322[index],
            inner_energy_flux[index], inner_angular_momentum_flux[index], spin_rate[index], scaled_time[index])
    for decay_constant in decay_constants_gev
        coupling = (decay_constant / planck_mass_gev)^2
        @printf(" %12.4e", scaled_time[index] / coupling * mass_unit_years)
    end
    println()
end
