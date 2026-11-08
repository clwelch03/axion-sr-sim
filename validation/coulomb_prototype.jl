# Prototype: outgoing Coulomb function H⁺_l(η, ρ) for complex η, ρ, as candidate exact outer data
# for propagating Fourier components (DLMF 33.2.11 with Tricomi's U from DLMF 13.4.4 on a rotated
# contour, so that the integrand decays like e^{−2ρs}). Checks: the Coulomb ODE residual, the
# large-ρ asymptotic series currently used by the solver, and the cost.
# Run from the package root:  julia --project=. validation/coulomb_prototype.jl
using FastGaussQuadrature, Printf, ForwardDiff

"""Complex log-gamma by the Lanczos approximation (g = 7, n = 9), valid for Re z > 1/2 via reflection otherwise."""
function log_gamma(z::Complex)
    coefficients = (0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
                    -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7)
    real(z) < 0.5 && return log(pi / sin(pi * z)) - log_gamma(1 - z)
    z -= 1
    accumulated = coefficients[1] + sum(coefficients[k + 1] / (z + k) for k in 1:8)
    t = z + 7.5
    return 0.5log(2pi) + (z + 0.5) * log(t) - t + log(accumulated)
end

"""
U(a, b, −2iρ) = (1/Γ(a)) ∫₀^∞ e^{−2ρs} (is)^{a−1} (1+is)^{b−a−1} i ds (Re a > 0), with x = 2ρs = e^v:

    U = i (2ρ)^{−a} e^{iπ(a−1)/2} / Γ(a) ∫_{−∞}^{∞} e^{a v − e^v} (1 + i e^v / 2ρ)^{b−a−1} dv,

whose integrand is analytic with double-exponential decay, so the trapezoidal rule converges
exponentially in the step `h` (the oscillation e^{i Im(a) v} is uniform in v).
"""
function tricomi_u_outgoing(a, b, rho; h=0.1, v_min=-14.0, v_max=6.0)
    total = zero(complex(a) * rho)
    for v in v_min:h:v_max
        x = exp(v)
        total += exp(a * v - x) * (1 + im * x / (2rho))^(b - a - 1)
    end
    return total * h * im * (2rho)^(-a) * exp(im * pi * (a - 1) / 2) * exp(-log_gamma(complex(a)))
end

"""H⁺_l(η, ρ) = e^{iθ_l} (−2iρ)^{l+1+iη} U(l+1+iη, 2l+2, −2iρ), θ_l = ρ − η ln 2ρ − lπ/2 + arg Γ(l+1+iη)."""
function coulomb_outgoing(l, eta, rho)
    a = l + 1 + im * eta
    sigma = imag(log_gamma(complex(a)))
    theta = rho - eta * log(2rho) - l * pi / 2 + sigma
    return exp(im * theta) * (-2im * rho)^a * tricomi_u_outgoing(a, 2l + 2, rho)
end

l, eta = 3, -4.28          # the near-threshold (−1, 2) component at a = 0.99, α = 0.3
println("ODE residual |u'' + (1 − 2η/ρ − l(l+1)/ρ²) u| / |u| by central differences:")
for rho in (5.0, 8.5, 21.0, 60.0, 200.0)
    h = 1e-2
    u(x) = coulomb_outgoing(l, eta, x)
    second = (u(rho + h) - 2u(rho) + u(rho - h)) / h^2
    @printf("  ρ = %6.1f: |u| = %.4e  residual = %.2e\n", rho, abs(u(rho)), abs(second + (1 - 2eta / rho - l * (l + 1) / rho^2) * u(rho)) / abs(u(rho)))
end
println("comparison with the asymptotic series U ~ z^{-a} Σ (a)_n (a−b+1)_n / (n! (−z)^n), truncated at the smallest term:")
for rho in (8.5, 21.0, 60.0, 200.0)
    a, b, z = l + 1 + im * eta, 2l + 2, -2im * rho
    term, series, smallest, n = one(z), one(z), Inf, 0
    while n < 40
        term *= (a + n) * (a - b + 1 + n) / ((n + 1) * (-z))
        abs(term) > smallest && break
        smallest = abs(term); series += term; n += 1
    end
    exact = tricomi_u_outgoing(a, b, rho) * z^a
    @printf("  ρ = %6.1f: series (%2d terms, last %.1e) = %.8f%+.8fi   exact z^a U = %.8f%+.8fi   rel. diff = %.1e\n",
            rho, n, smallest, real(series), imag(series), real(exact), imag(exact), abs(series - exact) / abs(exact))
end
println("quadrature convergence at ρ = 8.5 (step h):")
for h in (0.4, 0.2, 0.1, 0.05)
    @printf("  h = %.2f: %s\n", h, tricomi_u_outgoing(l + 1 + im * eta, 2l + 2, 8.5; h))
end
println("ForwardDiff through ρ works: ", ForwardDiff.derivative(x -> real(coulomb_outgoing(l, eta, x)), 8.5))
t = @elapsed for _ in 1:1000; coulomb_outgoing(l, eta, 8.5); end
@printf("cost: %.1f μs per evaluation\n", 1e3t)
