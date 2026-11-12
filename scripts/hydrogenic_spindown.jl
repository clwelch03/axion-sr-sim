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
# Run from the package root:
#   julia --project=. scripts/hydrogenic_spindown.jl <f_a/M_Pl> [initial spin, default 0.9] [final t/M₀, default 1e10]
# Writes results/hydrogenic_spindown_a<spin>_alpha0.3_fa<f_a/M_Pl>_diagnostics.csv (one row per saved time).
using OrdinaryDiffEqTsit5, Printf

length(ARGS) in (1, 2, 3) || error("usage: hydrogenic_spindown.jl <f_a/M_Pl> [initial spin] [final t/M₀]")
decay_constant_ratio = parse(Float64, ARGS[1])
spin_label = length(ARGS) ≥ 2 ? ARGS[2] : "0.9"
initial_spin = parse(Float64, spin_label)
final_time = length(ARGS) ≥ 3 ? parse(Float64, ARGS[3]) : 1e10
initial_mass_coupling = 0.3
coupling_squared = decay_constant_ratio^2
seed_occupations = [3e-6, 3e-6]
modes = ((n=2, l=1, m=1), (n=3, l=2, m=2))   # |211⟩, |322⟩ (radial overtone n − l − 1 = 0)

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

function rates!(derivative, state, _, _)
    occupations, mass, angular_momentum = state[1:2], state[3], state[4]
    mass_coupling, spin = initial_mass_coupling * mass, angular_momentum / mass^2
    mass_change, angular_momentum_change = 0.0, 0.0
    for (index, mode) in enumerate(modes)
        occupation_rate = 2growth_rate(mass_coupling, spin; mode...) / mass * occupations[index]
        derivative[index] = occupation_rate
        quantum_energy = initial_mass_coupling * (1 - mass_coupling^2 / (2mode.n^2))   # ω M₀ (μ fixed)
        mass_change -= coupling_squared * quantum_energy * occupation_rate
        angular_momentum_change -= coupling_squared * mode.m * occupation_rate
    end
    derivative[3], derivative[4] = mass_change, angular_momentum_change
    nothing
end

initial_state = [seed_occupations..., 1.0, initial_spin]
saved_times = sort(unique(vcat(range(0, min(final_time, 2e9); length=2001), exp10.(range(log10(2e9), log10(final_time); length=400)))))
problem = ODEProblem(rates!, initial_state, (0.0, final_time))
solution = solve(problem, Tsit5(); reltol=1e-10, abstol=1e-20, saveat=saved_times)
println("retcode: $(solution.retcode), $(length(solution.t)) saved states")

output_path = joinpath(@__DIR__, "..", "results", "hydrogenic_spindown_a$(spin_label)_alpha0.3_fa$(ARGS[1])_diagnostics.csv")
open(output_path, "w") do io
    println(io, "time,N1,N2,Im_omega1,Im_omega2,b,a,threshold")
    previous_time = -Inf
    for (time, state) in zip(solution.t, solution.u)
        # one row per time: the solver saves t = 0 twice, and the linear and log grids meet at 2e9
        time ≤ previous_time * (1 + 1e-12) && continue
        previous_time = time
        mass, spin = state[3], state[4] / state[3]^2
        mass_coupling = initial_mass_coupling * mass
        rates = [growth_rate(mass_coupling, spin; mode...) / mass for mode in modes]   # M₀ units
        threshold = 4mass_coupling / (1 + 4mass_coupling^2)                            # zero of the |211⟩ rate, m = 1
        println(io, join((@sprintf("%.10e", value) for value in (time, state[1], state[2], rates..., mass, spin, threshold)), ","))
    end
end
final = solution.u[end]
@printf("initial M ω_I: |211⟩ %.6e, |322⟩ %.6e\n", growth_rate(initial_mass_coupling, initial_spin; modes[1]...),
        growth_rate(initial_mass_coupling, initial_spin; modes[2]...))
@printf("final: t = %.3e, 𝒩 = (%.6g, %.6g), a = %.5f, M/M₀ = %.6f\nwrote %s\n", solution.t[end], final[1], final[2],
        final[4] / final[3]^2, final[3], output_path)
