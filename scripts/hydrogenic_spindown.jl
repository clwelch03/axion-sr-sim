# Free-field spin-down with leading small-α (hydrogenic) rates and frequencies: the analytic counterpart of
# `AXION_POTENTIAL=free scripts/evolution_spindown.jl`, for overlays (scripts/plot_spindown_overlay.jl).
#
# Rates: leading small-α superradiant rate of Ng et al., arXiv:1908.02312, Eq. (1), ported from
# main:src/common/hydrogenic.py (`growth_rate`, convention M ω_I = Γ/2, horizon factor m a − 2α r₊ with ω ≈ μ).
# For |211⟩ this is Detweiler's Γ = α⁹ (a − 2α r₊)[(1 − a²) + (a − 2α r₊)²]/24. Its zero is at Ω_H = μ/m, not at
# Ω_H = Re ω/m, so the hydrogenic evolution stops at a = 4(α/m)/(1 + 4(α/m)²) (current α), slightly above the
# accurate threshold; that is the threshold written to the CSV.
# Frequencies (energy per quantum): ω = μ(1 − α²/2n²), as `hydrogenic_frequency` in src/modes.jl.
#
# Units and state as in evolve_with_spindown: time in M₀, occupations 𝒩_i in f_a² M₀², b = M/M₀, j = J/M₀²,
# current α = α₀ b, a = j/b², κ = (f_a/M_Pl)²:
#   d𝒩_i/dt = 2 (M ω_I,i)/b 𝒩_i,   db/dt = −κ Σ_i (μM₀)(1 − α²/2n_i²) d𝒩_i/dt,   dj/dt = −κ Σ_i m_i d𝒩_i/dt.
# Same seed as the driver: 𝒩 = 3e-6 for |211⟩ and |322⟩. The modes are uncoupled (no self-interaction).
#
# Optional gravitational-wave (GW) annihilation, AXION_GW=1 (default off, so the earlier results reproduce).
# Two quanta of one level annihilate into a graviton that leaves to infinity. The power is
#   P_i = scale × power_coefficient_i(α) (M_cloud,i / M)²,   M_cloud,i / M₀ = κ (ω_i M₀) 𝒩_i,
# with the coefficient from the SuperRad table (scripts/gw_power.jl) at the current α, and `scale` = AXION_GW_SCALE
# (default 1; for testing the sensitivity to the coefficient). The level loses energy P_i, that is quanta at the rate
# P_i / (κ ω_i M₀):
#   d𝒩_i/dt = 2 (M ω_I,i)/b 𝒩_i − scale × power_coefficient_i(α) κ (ω_i M₀) 𝒩_i² / b².
# The hole equations keep only the first (superradiant) term: the radiated energy and angular momentum do not return
# to the hole. The radiated energy E_GW/M₀ = ∫ Σ_i P_i dt is a fifth state variable, for the output.
#
# Run from the package root:
#   julia --project=. scripts/hydrogenic_spindown.jl <f_a/M_Pl> [initial spin, default 0.9] [final t/M₀, default 1e10]
# Writes results/hydrogenic_spindown_a<spin>_alpha0.3_fa<f_a/M_Pl>_diagnostics.csv (one row per saved time);
# with AXION_GW=1 the name ends in _gw_diagnostics.csv (or _gw<scale>_diagnostics.csv for a scale other than 1).
using OrdinaryDiffEqTsit5, Printf
include(joinpath(@__DIR__, "gw_power.jl"))

length(ARGS) in (1, 2, 3) || error("usage: hydrogenic_spindown.jl <f_a/M_Pl> [initial spin] [final t/M₀]")
decay_constant_ratio = parse(Float64, ARGS[1])
spin_label = length(ARGS) ≥ 2 ? ARGS[2] : "0.9"
initial_spin = parse(Float64, spin_label)
final_time = length(ARGS) ≥ 3 ? parse(Float64, ARGS[3]) : 1e10
initial_mass_coupling = 0.3
coupling_squared = decay_constant_ratio^2
seed_occupations = [3e-6, 3e-6]
modes = ((n=2, l=1, m=1), (n=3, l=2, m=2))   # |211⟩, |322⟩ (radial overtone n − l − 1 = 0)
gw_emission = get(ENV, "AXION_GW", "0") == "1"
gw_scale_label = get(ENV, "AXION_GW_SCALE", "1")
gw_scale = gw_emission ? parse(Float64, gw_scale_label) : 0.0
gw_power_table = read_gw_power_table()

log_factorial(k) = sum(log, 1:k; init=0.0)   # all gamma-function arguments below are integers: Γ(k + 1) = k!

"""Leading small-α M ω_I (Ng et al. Eq. 1 / 2); port of main:src/common/hydrogenic.py `growth_rate`."""
function growth_rate(mass_coupling, spin; n, l, m)
    radial_overtone = n - l - 1
    horizon_radius = 1 + sqrt(1 - spin^2)
    horizon_factor = m * spin - 2mass_coupling * horizon_radius
    log_prefactor = (4l + 2) * log(2) + log_factorial(2l + 1 + radial_overtone) - (2l + 4) * log(n) -
                    log_factorial(radial_overtone) + 2 * (log_factorial(l) - log_factorial(2l) - log_factorial(2l + 1))
    angular_product = prod(k^2 * (1 - spin^2) + horizon_factor^2 for k in 1:l)
    exp(log_prefactor) * mass_coupling^(4l + 5) * horizon_factor * angular_product / 2
end

# check against the |211⟩ closed form (Detweiler, 1/24)
let α = 0.3, a = 0.9, x = a - 2α * (1 + sqrt(1 - a^2))
    detweiler = α^9 * x * ((1 - a^2) + x^2) / 24 / 2
    isapprox(growth_rate(α, a; modes[1]...), detweiler; rtol=1e-13) || error("growth_rate does not match Detweiler")
end

"""
Right-hand side. `parameters.superradiance` switches the exchange with the hole on or off (off only in the check
below); `parameters.gw_scale` multiplies the GW power (0: no GW emission).
"""
function rates!(derivative, state, parameters, _)
    occupations, mass, angular_momentum = state[1:2], state[3], state[4]
    mass_coupling, spin = initial_mass_coupling * mass, angular_momentum / mass^2
    mass_change, angular_momentum_change, radiated_power = 0.0, 0.0, 0.0
    for (index, mode) in enumerate(modes)
        quantum_energy = initial_mass_coupling * (1 - mass_coupling^2 / (2mode.n^2))   # ω M₀ (μ fixed)
        superradiant_rate = parameters.superradiance ?
                            2growth_rate(mass_coupling, spin; mode...) / mass * occupations[index] : 0.0
        cloud_mass_fraction = coupling_squared * quantum_energy * occupations[index] / mass   # M_cloud / M
        gw_power = parameters.gw_scale == 0 ? 0.0 :
                   parameters.gw_scale * gw_power_coefficient(gw_power_table, index, mass_coupling) * cloud_mass_fraction^2
        derivative[index] = superradiant_rate - gw_power / (coupling_squared * quantum_energy)
        mass_change -= coupling_squared * quantum_energy * superradiant_rate
        angular_momentum_change -= coupling_squared * mode.m * superradiant_rate
        radiated_power += gw_power
    end
    derivative[3], derivative[4], derivative[5] = mass_change, angular_momentum_change, radiated_power
    nothing
end

# Check of the GW term with an exact answer: a single |211⟩ cloud with the exchange with the hole switched off. The
# hole does not change, so the coefficient is constant and M_cloud(t) = M_cloud(0) / (1 + t/τ), with
# τ = M₀ / (power_coefficient × M_cloud(0)/M₀). Run over 3τ from a cloud of 1.85% of the hole mass.
let cloud_mass_fraction = 0.0185, quantum_energy = initial_mass_coupling * (1 - initial_mass_coupling^2 / 8)
    occupation = cloud_mass_fraction / (coupling_squared * quantum_energy)
    decay_time = 1 / (gw_power_coefficient(gw_power_table, 1, initial_mass_coupling) * cloud_mass_fraction)
    check = solve(ODEProblem(rates!, [occupation, 0.0, 1.0, initial_spin, 0.0], (0.0, 3decay_time),
                             (superradiance=false, gw_scale=1.0)), Tsit5(); reltol=1e-10, abstol=1e-20)
    final_fraction = coupling_squared * quantum_energy * check.u[end][1]
    relative_error = final_fraction / (cloud_mass_fraction / 4) - 1
    energy_error = (final_fraction + check.u[end][5]) / cloud_mass_fraction - 1   # cloud + radiated = initial cloud
    @printf("GW check: τ = %.4e M₀; M_cloud(3τ)/M_cloud(0) − 1/4: relative error %.1e; cloud + radiated energy: %.1e\n",
            decay_time, relative_error, energy_error)
    abs(relative_error) < 1e-8 && abs(energy_error) < 1e-8 || error("GW term does not reproduce M_cloud(0)/(1 + t/τ)")
end

initial_state = [seed_occupations..., 1.0, initial_spin, 0.0]
saved_times = sort(unique(vcat(range(0, min(final_time, 2e9); length=2001), exp10.(range(log10(2e9), log10(final_time); length=400)))))
problem = ODEProblem(rates!, initial_state, (0.0, final_time), (superradiance=true, gw_scale=gw_scale))
solution = solve(problem, Tsit5(); reltol=1e-10, abstol=1e-20, saveat=saved_times)
println("retcode: $(solution.retcode), $(length(solution.t)) saved states")

output_path = joinpath(@__DIR__, "..", "results", "hydrogenic_spindown_a$(spin_label)_alpha0.3_fa$(ARGS[1])" *
                       (gw_emission ? "_gw" * (gw_scale_label == "1" ? "" : gw_scale_label) : "") * "_diagnostics.csv")
open(output_path, "w") do io
    println(io, "time,N1,N2,Im_omega1,Im_omega2,b,a,threshold" * (gw_emission ? ",radiated_energy" : ""))
    previous_time = -Inf
    for (time, state) in zip(solution.t, solution.u)
        # one row per time: the solver saves t = 0 twice, and the linear and log grids meet at 2e9
        time ≤ previous_time * (1 + 1e-12) && continue
        previous_time = time
        mass, spin = state[3], state[4] / state[3]^2
        mass_coupling = initial_mass_coupling * mass
        rates = [growth_rate(mass_coupling, spin; mode...) / mass for mode in modes]   # M₀ units
        threshold = 4mass_coupling / (1 + 4mass_coupling^2)                            # zero of the |211⟩ rate, m = 1
        println(io, join((@sprintf("%.10e", value) for value in (time, state[1], state[2], rates..., mass, spin, threshold, (gw_emission ? (state[5],) : ())...)), ","))
    end
end
final = solution.u[end]
@printf("initial M ω_I: |211⟩ %.6e, |322⟩ %.6e\n", growth_rate(initial_mass_coupling, initial_spin; modes[1]...),
        growth_rate(initial_mass_coupling, initial_spin; modes[2]...))
@printf("final: t = %.3e, 𝒩 = (%.6g, %.6g), a = %.5f, M/M₀ = %.6f\nwrote %s\n", solution.t[end], final[1], final[2],
        final[4] / final[3]^2, final[3], output_path)
gw_emission && @printf("radiated in GW: E_GW/M₀ = %.6e (GW scale %s)\n", final[5], gw_scale_label)
