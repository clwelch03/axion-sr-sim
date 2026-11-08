# Outer-boundary check for the near-threshold (−1, 2) component at a = 0.99, α = 0.3 (joint-solution
# frequencies): integrate the *linear Kerr* radial equation inward from a radius where the 1/r series is
# converged, and compare at r = 400 with (a) the truncated series the solver imposes now, (b) the optimally
# truncated series, (c) the bare Coulomb function H⁺_l(η, kr)/r normalised at the far radius.
# Run from the package root:  julia --project=. validation/outer_integration_check.jl
using AxionCloud, Printf, LinearAlgebra

background = KerrBackground(0.99)
mass_coupling = 0.3
mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)])
component_index = findfirst(c -> c.q == [-1, 2], mode_set.components)
component = mode_set.components[component_index]
frequency = -(0.295948 - 4.64e-6im) + 2(0.298350 + 1.54e-6im)      # joint solution, NOTES day 2
mixing = mixing_matrix(AngularGrid(12, mode_set), mode_set, component_index)
l, m, a = component.ls[1], component.m, background.a
k = outer_wavenumber(frequency, mass_coupling)
nu = AxionCloud.outer_exponent(frequency, mass_coupling, k)
eta = -(2frequency^2 - mass_coupling^2) / k
@printf("component q = %s, l = %d, m = %d, ω_q = %.6f%+.2ei, k = %.5f%+.1ei, ν = %.3f%+.3fi, η = %.3f%+.3fi\n",
        component.q, l, m, real(frequency), imag(frequency), real(k), imag(k), real(nu), imag(nu), real(eta), imag(eta))

series(order) = OuterSeries(background, mass_coupling, frequency, m, component.ls, mixing, [1.0 + 0im], order)
value_at(ser, r) = endpoint_values(ser, background, r - background.rplus)[1][1]
# dR/dr from the tortoise slope: Q = dR/dr* = (Δ/(r²+a²)) dR/dr
slope_at(ser, r) = endpoint_values(ser, background, r - background.rplus)[2][1] / AxionCloud.dr_drstar(background, r - background.rplus)

# --- far radius: where is the series converged? ---------------------------------------------
for r_far in (600.0, 1000.0, 1500.0)
    v25, v35, v45 = value_at(series(25), r_far), value_at(series(35), r_far), value_at(series(45), r_far)
    @printf("r_far = %6.0f (ρ = %5.1f): |order 35 − 25|/|R| = %.1e, |order 45 − 35|/|R| = %.1e\n",
            r_far, real(k) * r_far, abs(v35 - v25) / abs(v35), abs(v45 - v35) / abs(v45))
end

# --- inward integration of the linear equation in r: y = (R, Δ dR/dr) ----------------------------
function linear_rhs(r, y)
    delta = r^2 - 2r + a^2
    B = ((r^2 + a^2) * frequency - a * m)^2 / delta + 2a * m * frequency - a^2 * frequency^2 - mass_coupling^2 * r^2 - l * (l + 1)
    return [y[2] / delta, -(B + a^2 * (frequency^2 - mass_coupling^2) * mixing[1, 1]) * y[1]]
end
function integrate_inward(r_far, r_end, step)
    r = r_far
    far = series(40)
    y = [value_at(far, r), (r^2 - 2r + a^2) * slope_at(far, r)]
    n_steps = round(Int, (r_far - r_end) / step)
    h = -(r_far - r_end) / n_steps
    for _ in 1:n_steps
        k1 = linear_rhs(r, y); k2 = linear_rhs(r + h / 2, y + h / 2 * k1)
        k3 = linear_rhs(r + h / 2, y + h / 2 * k2); k4 = linear_rhs(r + h, y + h * k3)
        y += h / 6 * (k1 + 2k2 + 2k3 + k4); r += h
    end
    return y[1], y[2] / (r^2 - 2r + a^2)      # R, dR/dr at r_end
end
r_end = 400.0
R_ref, dR_ref = integrate_inward(1500.0, r_end, 0.02)
R_alt, dR_alt = integrate_inward(1000.0, r_end, 0.04)
@printf("\nreference at r = %.0f (ρ = %.2f): R = %.10e%+.10ei, R'/R = %.8f%+.8fi\n", r_end, real(k) * r_end, real(R_ref), imag(R_ref), real(dR_ref / R_ref), imag(dR_ref / R_ref))
@printf("integration self-consistency (r_far 1000 vs 1500, step 0.04 vs 0.02): |ΔR|/|R| = %.1e, |Δ(R'/R)| = %.1e\n",
        abs(R_alt - R_ref) / abs(R_ref), abs(dR_alt / R_alt - dR_ref / R_ref))
t = @elapsed integrate_inward(1000.0, r_end, 0.04)
@printf("cost of one inward integration (r 1000 → 400, step 0.04): %.1f ms\n", 1e3t)

# --- what the solver imposes now: truncated series at r = 400 ------------------------------------
println("\nseries at r = 400 versus the integrated solution:")
for order in (7, 11, 15, 20, 30)
    ser = series(order)
    R_s, dR_s = value_at(ser, r_end), slope_at(ser, r_end)
    last_term = abs(ser.coefficients[end, 1]) / r_end^order / abs(R_s / exp(im * k * r_end + nu * log(r_end)))
    @printf("  J_B = %2d: |ΔR|/|R| = %.1e   |Δ(R'/R)| = %.1e   (last term/sum = %.1e)\n", order,
            abs(R_s - R_ref) / abs(R_ref), abs(dR_s / R_s - dR_ref / R_ref), last_term)
end

# --- the bare Coulomb function, normalised at the far radius --------------------------------------
function log_gamma(z::Complex)
    c = (0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
         -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7)
    real(z) < 0.5 && return log(pi / sin(pi * z)) - log_gamma(1 - z)
    z -= 1
    acc = c[1] + sum(c[j + 1] / (z + j) for j in 1:8)
    t = z + 7.5
    return 0.5log(2pi) + (z + 0.5) * log(t) - t + log(acc)
end
function tricomi_u_outgoing(a_, b, rho; h=0.1, v_min=-14.0, v_max=6.0)
    total = zero(complex(a_) * rho)
    for v in v_min:h:v_max
        x = exp(v)
        total += exp(a_ * v - x) * (1 + im * x / (2rho))^(b - a_ - 1)
    end
    return total * h * im * (2rho)^(-a_) * exp(im * pi * (a_ - 1) / 2) * exp(-log_gamma(complex(a_)))
end
function coulomb_outgoing(l_, eta_, rho)
    a_ = l_ + 1 + im * eta_
    theta = rho - eta_ * log(2rho) - l_ * pi / 2 + imag(log_gamma(complex(a_)))
    return exp(im * theta) * (-2im * rho)^a_ * tricomi_u_outgoing(a_, 2l_ + 2, rho)
end
coulomb_profile(r) = coulomb_outgoing(l, eta, k * r) / r
for r_far in (1000.0, 1500.0, 3000.0)
    normalisation = value_at(series(40), r_far) / coulomb_profile(r_far)
    R_c = normalisation * coulomb_profile(r_end)
    h = 1e-3
    dR_c = normalisation * (coulomb_profile(r_end + h) - coulomb_profile(r_end - h)) / 2h
    @printf("bare Coulomb H⁺_%d(η, kr)/r normalised at r = %.0f: at r = 400 |ΔR|/|R| = %.1e, |Δ(R'/R)| = %.1e\n",
            l, r_far, abs(R_c - R_ref) / abs(R_ref), abs(dR_c / R_c - dR_ref / R_ref))
end
