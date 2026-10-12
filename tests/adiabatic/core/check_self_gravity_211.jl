# julia --project=. tests/adiabatic/core/check_self_gravity_211.jl
# Standalone dilute Newtonian benchmark, not a nonlinear/Kerr gravity solver.
# Derivation: docs/adiabatic-self-gravity-review.md, "211 benchmark checkpoint".
using Test
using Printf
include(joinpath(@__DIR__, "../../../src/adiabatic/core/load.jl"))

# x = r/a_B, a_B = 1/(mu*alpha), rho = (M_cloud/a_B^3) f(x,theta).
# f = x^2 exp(-x)/24 * |Y_11|^2, with integral f d^3x = 1.
# Legendre convention: f = sum_L f_L P_L, U = (G*M_cloud/a_B) sum_L v_L P_L.
function hydrogenic_211_poisson_check(radial_step, outer_radius, angular_count)
    radius = collect(0:radial_step:outer_radius)
    cosine_theta, angular_weights = gausslegendre(angular_count)
    # Repo polar harmonics integrate to one over d(cos theta), hence 1/(2pi).
    angular_density = polar_harmonic.(1, 1, cosine_theta).^2 ./ (2pi)
    legendre(degree) = sqrt(2 / (2degree + 1)) .* polar_harmonic.(degree, 0, cosine_theta)
    angular_coefficient(degree) = (2degree + 1) / 2 *
        sum(angular_weights .* angular_density .* legendre(degree))
    radial_density = radius.^2 .* exp.(-radius) ./ 24
    density_multipoles = [angular_coefficient(degree) .* radial_density for degree in (0, 2)]
    potentials = Vector{Float64}[]
    shift_coefficients = Float64[]

    for (degree, density) in zip((0, 2), density_multipoles)
        interior = cumulative_trapezoid(radius, radius.^(degree + 2) .* density)
        exterior_integrand = zeros(length(radius))
        exterior_integrand[2:end] = radius[2:end].^(1 - degree) .* density[2:end]
        # Backward integration avoids subtracting almost equal cumulative totals.
        exterior = reverse(cumulative_trapezoid(
            -reverse(radius), reverse(exterior_integrand)
        ))
        potential = zeros(length(radius))
        potential[2:end] = -4pi / (2degree + 1) .* (
            interior[2:end] ./ radius[2:end].^(degree + 1) .+
            radius[2:end].^degree .* exterior[2:end]
        )
        potential[1] = degree == 0 ? -4pi * exterior[1] : 0.0
        push!(potentials, potential)
        # Angular integral P_L^2 = 4pi/(2L+1); positive coefficient of -G*M/a_B.
        push!(shift_coefficients, -4pi / (2degree + 1) * last(cumulative_trapezoid(
            radius, radius.^2 .* density .* potential
        )))
    end
    mass = 4pi * last(cumulative_trapezoid(radius, radius.^2 .* first(density_multipoles)))
    return (; mass, shift_coefficients, potentials, radius,
            angular_coefficients=angular_coefficient.((0, 2, 4)))
end

# Independent exact evaluation of the symmetric radial double integral:
# J_L = (2/24^2) integral dx x^(3-L)e^-x integral_0^x dy y^(4+L)e^-y.
# Expand the inner integral as a finite incomplete-gamma polynomial.
function exact_211_radial_integral(degree)
    polynomial_integral = sum(
        factorial(3 - degree + power) //
        (factorial(power) * 2^(4 - degree + power))
        for power in 0:(4 + degree)
    )
    return factorial(4 + degree) // 288 *
        (factorial(3 - degree) - polynomial_integral)
end

@testset "Dilute 211 self-gravity" begin
    exact_monopole = exact_211_radial_integral(0)
    exact_quadrupole = exact_211_radial_integral(2) / 25
    @test exact_monopole == 93 // 512
    @test exact_quadrupole == 9 // 2560
    @test exact_monopole + exact_quadrupole == 237 // 1280
    exact_total = Float64(exact_monopole + exact_quadrupole)

    println("x_max = 40; radius is in Bohr radii; coefficients multiply -alpha^2 M_cloud/M_BH")
    println("dx        mass           monopole       quadrupole     total error")
    errors = Float64[]
    for radial_step in (0.08, 0.04, 0.02, 0.01)
        result = hydrogenic_211_poisson_check(radial_step, 40.0, 8)
        monopole, quadrupole = result.shift_coefficients
        error = abs(monopole + quadrupole - exact_total)
        push!(errors, error)
        @printf("%.3f     %.10f   %.10f   %.10f   %.3e\n",
                radial_step, result.mass, monopole, quadrupole, error)
    end
    # The radial potential uses the composite trapezoidal rule: O(dx^2).
    @test all(3.8 .< errors[1:end-1] ./ errors[2:end] .< 4.2)
    result = hydrogenic_211_poisson_check(0.01, 40.0, 8)
    @test result.mass ≈ 1 rtol=1e-10
    @test result.shift_coefficients[1] ≈ Float64(exact_monopole) rtol=5e-6
    @test result.shift_coefficients[2] ≈ Float64(exact_quadrupole) rtol=2e-5
    @test result.angular_coefficients[1] ≈ 1/(4pi) rtol=1e-13
    @test result.angular_coefficients[2] ≈ -1/(4pi) rtol=1e-13
    @test abs(result.angular_coefficients[3]) < 1e-14
    @test first(result.potentials[1]) ≈ -1/4 rtol=1e-8
    @test first(result.potentials[2]) == 0
    @test last(result.radius) * last(result.potentials[1]) ≈ -1 rtol=1e-10
    @test last(result.radius)^3 * last(result.potentials[2]) ≈ 6 rtol=1e-9

    # Angular and outer-boundary checks independent of radial refinement.
    angular_refined = hydrogenic_211_poisson_check(0.01, 40.0, 16)
    @test angular_refined.shift_coefficients ≈ result.shift_coefficients rtol=1e-13
    outer_refined = hydrogenic_211_poisson_check(0.01, 60.0, 8)
    @test outer_refined.shift_coefficients ≈ result.shift_coefficients rtol=1e-11
    println("Exact total = 237/1280 = ", exact_total)
    println("Angular (8 -> 16 nodes) and outer-radius (40 -> 60) checks passed.")
end
