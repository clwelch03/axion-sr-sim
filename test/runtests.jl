using Test
using AxionCloud
using ForwardDiff, LinearAlgebra, SparseArrays, FastGaussQuadrature
const AC = AxionCloud

@testset "AxionCloud" begin

@testset "Kerr geometry" begin
    background = KerrBackground(0.9)
    @test background.rplus^2 + background.a^2 ≈ 2background.rplus
    for rstar in (-80.0, -5.0, 3.0, 180.0)
        offset = offset_from_tortoise(background, rstar)
        @test tortoise(background, offset) ≈ rstar atol=1e-10
    end
    # Schwarzschild: r* = r + 2 ln(r − 2)
    schwarzschild = KerrBackground(0.0)
    offset = offset_from_tortoise(schwarzschild, 50.0)
    @test radius(schwarzschild, offset) + 2log(offset) ≈ 50.0
    @test_throws DomainError KerrBackground(1.0)
end

@testset "Potentials" begin
    for potential in (CosinePotential(), QuarticPotential(1.0), ChiralPotential())
        small = 1e-3
        @test AC.potential_energy(potential, small) ≈ small^2 / 2 rtol=1e-5
        # h = Φ − U' must be consistent with U
        value = 0.9
        @test value - AC.nonlinearity(potential, value) ≈ ForwardDiff.derivative(x -> AC.potential_energy(potential, x), value)
    end
    @test AC.nonlinearity(QuarticPotential(1.0), 0.01) ≈ 0.01^3 / 6
    @test AC.nonlinearity(FreePotential(), 0.3) == 0
end

@testset "Mode bookkeeping" begin
    mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)])
    @test AC.n_components(mode_set) == 8
    @test all(isodd(AC.fourier_order(component)) for component in mode_set.components)
    @test all(sum(component.q) > 0 for component in mode_set.components)
    @test fundamental_channels(mode_set) == [1, 2]
    @test [component.m for component in mode_set.components[1:2]] == [1, 2]
    # parity: (2,−1) has m = 0 and even parity so l = 0
    static = findfirst(component -> component.q == [2, -1], mode_set.components)
    @test mode_set.components[static].ls == [0]
    # φ₁²φ₂* has complex frequency 2ω̃₁ − ω̃₂*: it grows with γ₂ as well as γ₁.
    mode_frequencies = [0.29 + 1e-3im, 0.295 - 2e-4im]
    @test AC.component_frequencies(mode_set, mode_frequencies)[static] ≈ 2mode_frequencies[1] - conj(mode_frequencies[2])
    @test AC.component_frequencies(mode_set, real.(mode_frequencies))[static] ≈ 2 * 0.29 - 0.295
    with_angular = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)]; n_angular=2)
    @test with_angular.components[1].ls == [1, 3]
    @test_throws ArgumentError BoundMode(1, 1, 0)
end

@testset "Polar harmonics" begin
    x, weights = gausslegendre(40)
    for (l, m) in ((1, 1), (3, 2), (4, 0), (5, 3))
        @test sum(weights .* polar_harmonic.(l, m, x) .^ 2) ≈ 1
    end
    @test abs(sum(weights .* polar_harmonic.(3, 1, x) .* polar_harmonic.(1, 1, x))) < 1e-14
    theta = 0.7
    finite_difference = (polar_harmonic(3, 2, cos(theta + 1e-6)) - polar_harmonic(3, 2, cos(theta - 1e-6))) / 2e-6
    @test AC.polar_harmonic_dtheta(3, 2, cos(theta)) ≈ finite_difference rtol=1e-7
    mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)]; n_angular=2)
    mixing = mixing_matrix(AngularGrid(12, mode_set), mode_set, 1)
    @test mixing ≈ mixing'
    @test mixing[1, 1] ≈ 1 / 5   # ∫ Ỹ₁₁² x² dx with Ỹ₁₁² = (3/4)(1 − x²)
end

@testset "Simpson weights on non-uniform grids" begin
    uniform = collect(range(0.0, 2.0; length=11))
    @test sum(AC.simpson_weights(uniform) .* uniform .^ 3) ≈ 4.0          # cubic exact on a uniform grid, odd n
    even_count = collect(range(0.0, 2.0; length=10))
    @test sum(AC.simpson_weights(even_count) .* even_count .^ 2) ≈ 8 / 3   # quadratic exact with a trailing interval
    stretched = vcat(range(0.0, 1.0; length=6), 1.0 .+ cumsum(0.2 .* 1.3 .^ (1:7)))
    @test sum(AC.simpson_weights(stretched) .* stretched .^ 2) ≈ stretched[end]^3 / 3
    @test sum(AC.simpson_weights(stretched)) ≈ stretched[end]
    background = KerrBackground(0.9)
    mesh = RadialMesh(background; rstar_min=-10.0, rstar_max=100.0, step=1.0, stretch_from=20.0, stretch_factor=1.5, max_step=10.0)
    @test mesh.rstar[1] == -10.0 && mesh.rstar[end] == 100.0
    @test all(diff(mesh.rstar)[1:30] .≈ 1.0) && maximum(diff(mesh.rstar)) ≤ 10.0 + 1e-12
end

@testset "Domain estimate" begin
    background = KerrBackground(0.99)
    mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)])
    estimate = domain_estimate(background, 0.3, mode_set; integrated_outer=false)
    @test estimate.slowest_decay ≈ 0.3^2 / 3 rtol=0.05           # κ₃₂₂ ≈ α²/n
    @test 0 < estimate.smallest_wavenumber < 0.03                   # the near-threshold (−1, 2) component
    @test estimate.capped                                           # its series would need r ≫ the usable radius
    @test estimate.normalization_radius < 1000                      # |211⟩ is negligible at the boundary beyond this
    @test estimate.asymptotic_radius ≈ estimate.normalization_radius  # so the normalization limit sets the outer endpoint
    @test offset_from_tortoise(background, estimate.rstar_min) ≈ 1e-6 * background.horizon_gap rtol=1e-6
    # With integrated outer data the mesh only has to contain the bound components.
    integrated = domain_estimate(background, 0.3, mode_set)
    @test !integrated.capped
    @test radius(background, offset_from_tortoise(background, integrated.rstar_max)) ≈ integrated.containment_radius rtol=1e-6
    @test integrated.containment_radius ≈ estimate.containment_radius
    @test integrated.far_radius > integrated.containment_radius     # the (−1, 2) series only converges far out
    @test integrated.rstar_min == estimate.rstar_min
    # containment coordinate: r^ν e^{−κr} at x = κr must be `containment` × its peak
    x = AC._containment_coordinate(1.8, 1e-4)
    @test (x / 1.8)^1.8 * exp(-(x - 1.8)) ≈ 1e-4 rtol=1e-8
    mesh = RadialMesh(background, 0.3, mode_set; max_radius=300.0)
    @test mesh.rstar[1] ≈ estimate.rstar_min
    @test radius(background, mesh.offsets[end]) ≈ max(300.0, estimate.containment_radius) rtol=1e-6
    @test minimum(diff(mesh.rstar)) ≈ 0.4 rtol=1e-2   # the uniform step is rounded to fit the interval
end

@testset "Nonlinear source" begin
    background = KerrBackground(0.9)
    mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)])
    quadrature = FieldQuadrature(mode_set)
    profiles = zeros(ComplexF64, AC.n_channels(mode_set))
    profiles[1] = 0.1
    source = nonlinear_source(quadrature, QuarticPotential(1.0), background, 0.3, 5.0, profiles)
    # Φ = 2RỸ cos ξ: h = Φ³/6 has e^{-iξ} coefficient R³Ỹ³/2 and e^{-3iξ} coefficient R³Ỹ³/6.
    x, weights = gausslegendre(40)
    harmonic = polar_harmonic.(1, 1, x)
    rho_squared = AC.rho2.(Ref(background), 5.0, x)
    @test source[1] ≈ 0.3^2 * sum(weights .* harmonic .^ 4 .* rho_squared) * 0.1^3 / 2
    @test source[3] ≈ 0.3^2 * sum(weights .* polar_harmonic.(3, 3, x) .* harmonic .^ 3 .* rho_squared) * 0.1^3 / 6
    @test all(abs(source[channel]) < 1e-18 for channel in (2, 4, 5, 6, 7, 8))
    @test all(iszero, nonlinear_source(quadrature, FreePotential(), background, 0.3, 5.0, profiles))
end

@testset "Field kernels match the complex formulas" begin
    background = KerrBackground(0.9)
    mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)]; n_angular=2)
    quadrature = FieldQuadrature(mode_set; n_polar=8, n_phase=8)
    profiles = 0.4 .* randn(ComplexF64, AC.n_channels(mode_set))
    phase_factors = quadrature.phase.phase_factors
    component_values = AC.component_profiles(quadrature, profiles)
    @test component_values ≈ reduce(hcat, quadrature.angular.harmonics[:, channels] * profiles[channels]
                                    for channels in mode_set.component_channels)
    field = reconstruct_field(quadrature, profiles)
    @test field ≈ 2 .* real.(phase_factors * transpose(component_values)) rtol=1e-14
    @test AC.project_phases(quadrature, sin.(field)) ≈ (phase_factors' * sin.(field)) ./ AC.n_points(quadrature.phase) rtol=1e-14
    potential = ChiralPotential()
    @test AC.phase_averaged_potential(quadrature, potential, profiles) ≈
          vec(sum(AC.potential_energy.(Ref(potential), field); dims=1)) ./ AC.n_points(quadrature.phase) rtol=1e-14
    # Source Jacobian against forward-mode differentiation of the source in the real and imaginary parts.
    s = 4.0
    source_reals(reals) = AC._realify(nonlinear_source(quadrature, potential, background, 0.3, s,
                                                       complex.(reals[1:2:end], reals[2:2:end])))
    autodiff = ForwardDiff.jacobian(source_reals, AC._realify(profiles))
    analytic = AC.nonlinear_source_jacobian(quadrature, potential, background, 0.3, s, profiles)
    @test AC._realify(vec(analytic.wrt_real)) ≈ vec(autodiff[:, 1:2:end]) rtol=1e-12
    @test AC._realify(vec(analytic.wrt_imag)) ≈ vec(autodiff[:, 2:2:end]) rtol=1e-12
    # Dual-number source (value from the plain kernel, partials analytic): along a profile direction it must equal
    # the source Jacobian applied to that direction, and through the offset a central difference.
    direction = randn(ComplexF64, length(profiles))
    along_profiles = ForwardDiff.derivative(e -> nonlinear_source(quadrature, potential, background, 0.3, s, profiles .+ e .* direction), 0.0)
    @test along_profiles ≈ analytic.wrt_real * real.(direction) .+ analytic.wrt_imag * imag.(direction) rtol=1e-12
    along_offset = ForwardDiff.derivative(offset -> nonlinear_source(quadrature, potential, background, 0.3, offset, profiles), s)
    step = 1e-4
    @test along_offset ≈ (nonlinear_source(quadrature, potential, background, 0.3, s + step, profiles) .-
                          nonlinear_source(quadrature, potential, background, 0.3, s - step, profiles)) ./ 2step rtol=1e-7
end

@testset "Dual-number source against finite differences" begin
    # Independent of `nonlinearity_derivative`: central differences of the plain kernel in α, spin, offset and profiles.
    mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)]; n_angular=2)
    quadrature = FieldQuadrature(mode_set; n_polar=8, n_phase=8)
    profiles = [0.3 * cis(0.7index) / index for index in 1:AC.n_channels(mode_set)]
    for potential in (CosinePotential(), ChiralPotential(), QuarticPotential(0.5))
        function source_reals(inputs)
            mass_coupling, spin, offset = inputs[1], inputs[2], inputs[3]
            AC._realify(nonlinear_source(quadrature, potential, KerrBackground(spin), mass_coupling, offset,
                                         complex.(inputs[4:2:end], inputs[5:2:end])))
        end
        inputs = vcat(0.3, 0.9, 4.0, AC._realify(profiles))
        dual = ForwardDiff.jacobian(source_reals, inputs)
        step = 1e-5
        difference = reduce(hcat, [(source_reals(inputs .+ step .* (1:length(inputs) .== k)) .-
                                    source_reals(inputs .- step .* (1:length(inputs) .== k))) ./ 2step for k in eachindex(inputs)])
        @test maximum(abs, dual .- difference) < 1e-7 * maximum(abs, difference)
    end
end

@testset "Asymptotic series satisfy the radial equation" begin
    background = KerrBackground(0.9)
    mass_coupling = 0.3
    frequency = 0.28 + 1e-6im
    mixing = fill(1 / 5, 1, 1)
    function ode_residual(series, s)
        profile(offset) = endpoint_values(series, background, offset)[1][1]
        step = 1e-4 * s
        flux(offset) = AC.delta(background, offset) * (profile(offset + step) - profile(offset - step)) / (2step)
        dflux = (flux(s + step) - flux(s - step)) / (2step)
        r = radius(background, s)
        a = background.a
        B = ((r^2 + a^2) * frequency - a)^2 / AC.delta(background, s) + 2a * frequency - a^2 * frequency^2 - mass_coupling^2 * r^2 - 2
        return abs(dflux + B * profile(s) + a^2 * (frequency^2 - mass_coupling^2) * mixing[1, 1] * profile(s)) / abs(B * profile(s))
    end
    horizon = HorizonSeries(background, mass_coupling, frequency, 1, [1], mixing, [1.0 + 0im], 5)
    @test ode_residual(horizon, 1e-2) < 1e-5
    outer = OuterSeries(background, mass_coupling, frequency, 1, [1], mixing, [1.0 + 0im], 7)
    @test ode_residual(outer, 60.0) < 1e-5
    # Slopes are consistent with the profiles.
    for (series, s) in ((horizon, 1e-2), (outer, 60.0))
        R, Q = endpoint_values(series, background, s)
        step = 1e-5 * s
        dR_ds = (endpoint_values(series, background, s + step)[1][1] - endpoint_values(series, background, s - step)[1][1]) / (2step)
        @test Q[1] ≈ dR_ds * AC.dr_drstar(background, s) rtol=1e-7
    end
end

@testset "Integrated outer data" begin
    # (a) Where the series converges at the endpoint, integrating from far out reproduces it.
    background = KerrBackground(0.9)
    mass_coupling = 0.3
    frequency = 0.8 + 1e-6im                     # k ≈ 0.74: the series is converged at r ≈ 62
    mixing = fill(1 / 5, 1, 1)
    offset = 60.0
    reference = OuterSeries(background, mass_coupling, frequency, 1, [1], mixing, [1.0 + 0im], 30)
    R_series, Q_series = endpoint_values(reference, background, offset)
    forced = integrated_outer_values(background, mass_coupling, frequency, 1, [1], mixing, [1.0 + 0im],
                                     OuterBoundary(convergence_margin=2000.0), offset)
    @test forced.n_steps > 1000                  # genuinely integrated over a dozen wavelengths
    @test forced.profiles[1] ≈ R_series[1] rtol=1e-11
    @test forced.slopes[1] ≈ Q_series[1] rtol=1e-11
    direct = integrated_outer_values(background, mass_coupling, frequency, 1, [1], mixing, [1.0 + 0im], OuterBoundary(), offset)
    @test direct.n_steps == 0                    # no integration needed: the far series is evaluated at the endpoint
    @test direct.profiles[1] ≈ R_series[1] rtol=1e-13

    # (b) The Gauss–Legendre stepper reproduces the outgoing Coulomb function (Tricomi-U representation,
    #     validation/coulomb_prototype.jl) when integrating the Coulomb equation inward.
    function log_gamma(z::Complex)
        coefficients = (0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
                        -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7)
        real(z) < 0.5 && return log(pi / sin(pi * z)) - log_gamma(1 - z)
        z -= 1
        accumulated = coefficients[1] + sum(coefficients[k + 1] / (z + k) for k in 1:8)
        t = z + 7.5
        return 0.5log(2pi) + (z + 0.5) * log(t) - t + log(accumulated)
    end
    function tricomi_u_outgoing(a, b, rho; h=0.1, v_min=-14.0, v_max=6.0)
        total = zero(complex(a) * rho)
        for v in v_min:h:v_max
            x = exp(v)
            total += exp(a * v - x) * (1 + im * x / (2rho))^(b - a - 1)
        end
        return total * h * im * (2rho)^(-a) * exp(im * pi * (a - 1) / 2) * exp(-log_gamma(complex(a)))
    end
    function coulomb_outgoing(l, eta, rho)
        a = l + 1 + im * eta
        theta = rho - eta * log(2rho) - l * pi / 2 + imag(log_gamma(complex(a)))
        return exp(im * theta) * (-2im * rho)^a * tricomi_u_outgoing(a, 2l + 2, rho)
    end
    l, eta = 3, -4.28
    rho_far, rho_end = 60.0, 8.5
    far_value = coulomb_outgoing(l, eta, rho_far)
    far_slope = complex(ForwardDiff.derivative(x -> real(coulomb_outgoing(l, eta, x)), rho_far),
                        ForwardDiff.derivative(x -> imag(coulomb_outgoing(l, eta, x)), rho_far))
    coulomb_matrix(rho) = [0 1; -(1 - 2eta / rho - l * (l + 1) / rho^2) 0]
    integrated = gauss_legendre_integrate(coulomb_matrix, rho_far, rho_end, [far_value, far_slope], 2000)
    @test integrated[1] ≈ coulomb_outgoing(l, eta, rho_end) rtol=1e-9
    # Order 6: halving the step reduces the (tiny) error by ≈ 64.
    coarse = gauss_legendre_integrate(coulomb_matrix, rho_far, rho_end, [far_value, far_slope], 250)[1]
    fine = gauss_legendre_integrate(coulomb_matrix, rho_far, rho_end, [far_value, far_slope], 500)[1]
    exact = coulomb_outgoing(l, eta, rho_end)
    @test 30 < abs(coarse - exact) / abs(fine - exact) < 130

    # (c) Far-radius and step independence for the near-threshold (−1, 2) component at a = 0.99, α = 0.3,
    #     and (d) the frequency derivative by ForwardDiff matches a finite difference.
    background = KerrBackground(0.99)
    mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)])
    component_index = findfirst(component -> component.q == [-1, 2], mode_set.components)
    component = mode_set.components[component_index]
    mixing = mixing_matrix(AngularGrid(12, mode_set), mode_set, component_index)
    frequency = 0.300752 + 7.72e-6im
    offset = 400.0 - background.rplus
    data_for(settings, frequency) = integrated_outer_values(background, mass_coupling, frequency, component.m, component.ls,
                                                            mixing, [1.0 + 0im], settings, offset)
    default = data_for(OuterBoundary(), frequency)
    @test default.n_steps > 100 && default.far_series_tail < 1e-15
    for settings in (OuterBoundary(convergence_margin=120.0), OuterBoundary(steps_per_wavelength=96),
                     OuterBoundary(convergence_margin=120.0, steps_per_wavelength=96, far_order=50))
        other = data_for(settings, frequency)
        @test other.profiles[1] ≈ default.profiles[1] rtol=1e-10
        @test other.slopes[1] / other.profiles[1] ≈ default.slopes[1] / default.profiles[1] rtol=1e-10
    end
    profile_real(w) = real(data_for(OuterBoundary(), complex(w[1], w[2])).profiles[1])
    gradient_ad = ForwardDiff.gradient(profile_real, [real(frequency), imag(frequency)])
    step = 1e-8
    finite = [(profile_real([real(frequency) + step, imag(frequency)]) - profile_real([real(frequency) - step, imag(frequency)])) / 2step,
              (profile_real([real(frequency), imag(frequency) + step]) - profile_real([real(frequency), imag(frequency) - step])) / 2step]
    @test gradient_ad ≈ finite rtol=1e-6
    # Bound components can be integrated too (decaying solution, stable inward): the |322⟩ fundamental.
    bound_index = findfirst(component -> component.q == [0, 1], mode_set.components)
    bound = mode_set.components[bound_index]
    bound_mixing = mixing_matrix(AngularGrid(12, mode_set), mode_set, bound_index)
    bound_frequency = 0.298350 + 1.54e-6im
    bound_data = integrated_outer_values(background, mass_coupling, bound_frequency, bound.m, bound.ls, bound_mixing, [1.0 + 0im],
                                         OuterBoundary(include_bound=true), offset)
    bound_other = integrated_outer_values(background, mass_coupling, bound_frequency, bound.m, bound.ls, bound_mixing, [1.0 + 0im],
                                          OuterBoundary(include_bound=true, convergence_margin=120.0, steps_per_wavelength=96), offset)
    @test bound_data.n_steps > 50
    @test bound_data.profiles[1] ≈ bound_other.profiles[1] rtol=1e-10
    bound_series = endpoint_values(OuterSeries(background, mass_coupling, bound_frequency, bound.m, bound.ls, bound_mixing, [1.0 + 0im], 7),
                                   background, offset)[1][1]
    @test abs(bound_series / bound_data.profiles[1] - 1) < 1e-3   # the κr ≈ 12 series is close but not exact
end

@testset "Analytic right-hand-side Jacobian matches automatic differentiation" begin
    background = KerrBackground(0.99)
    mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)]; n_angular=2)
    mesh = RadialMesh(background; rstar_min=-10.0, rstar_max=30.0, step=10.0)
    for potential in (ChiralPotential(), CosinePotential(), FreePotential())
        problem = RadialProblem(background, potential, 0.3, mode_set; mesh, n_phase=8, n_polar=8)
        block = 0.5 .* randn(AC.node_block_size(AC.UnknownLayout(problem)))
        frequency_reals = [0.29, 1e-3, 0.295, -2e-4]
        for s in (0.3, 20.0)
            analytic = AC.rhs_derivatives(problem, s, block, frequency_reals, 0.05)
            autodiff = AC.rhs_derivatives_autodiff(problem, s, block, frequency_reals, 0.05)
            @test maximum(abs, analytic.wrt_state .- autodiff.wrt_state) < 1e-12 * maximum(abs, autodiff.wrt_state)
            @test maximum(abs, analytic.wrt_frequencies .- autodiff.wrt_frequencies) < 1e-12 * maximum(abs, autodiff.wrt_frequencies)
        end
        for value in (0.3, 1.7)
            @test AC.nonlinearity_derivative(potential, value) ≈ ForwardDiff.derivative(x -> AC.nonlinearity(potential, x), value)
        end
    end
end

@testset "Structured Jacobian matches automatic differentiation" begin
    background = KerrBackground(0.9)
    mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)])
    mesh = RadialMesh(background; rstar_min=-10.0, rstar_max=30.0, step=10.0)
    problem = RadialProblem(background, QuarticPotential(1.0), 0.3, mode_set; mesh, n_phase=8, n_polar=8)
    layout = AC.UnknownLayout(problem)
    amplitudes, scale = [0.05, 0.02], 0.05
    unknowns = AC.initial_unknowns(problem, layout, amplitudes, scale) .+ 0.01 .* randn(AC.n_unknowns(layout))
    structured = AC.jacobian(problem, layout, unknowns, amplitudes, scale)
    naive = ForwardDiff.jacobian(u -> AC.residual(problem, layout, u, amplitudes, scale), unknowns)
    @test maximum(abs, Matrix(structured) .- naive) < 1e-10 * maximum(abs, naive)
end

@testset "Jacobian blocks, in-place assembly and boundary derivatives" begin
    background = KerrBackground(0.9)
    mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)]; n_angular=2)
    mesh = RadialMesh(background; rstar_min=-10.0, rstar_max=30.0, step=10.0)
    problem = RadialProblem(background, CosinePotential(), 0.3, mode_set; mesh, n_phase=8, n_polar=8)
    layout = AC.UnknownLayout(problem)
    amplitudes, scale = [0.05, 0.02], 0.05
    unknowns = AC.initial_unknowns(problem, layout, amplitudes, scale) .+ 0.01 .* randn(AC.n_unknowns(layout))
    workspace, pattern = AC.JacobianWorkspace(problem, layout), AC.JacobianPattern(layout)
    in_place = AC.jacobian(problem, layout, unknowns, amplitudes, scale; workspace, pattern)
    displaced = unknowns .+ 1e-3 .* randn(length(unknowns))
    AC.jacobian_blocks!(workspace, problem, layout, displaced, amplitudes, scale)
    @test AC.jacobian(problem, layout, unknowns, amplitudes, scale; workspace, pattern) == in_place   # refill is complete
    @test AC.jacobian(problem, layout, unknowns, amplitudes, scale) == in_place
    dense, blocks = Matrix(in_place), workspace.blocks
    parameter_columns = AC.parameter_range(layout)
    interval = 2
    @test dense[AC.interval_rows(layout, interval), AC.node_range(layout, interval)] == blocks.wrt_left[:, :, interval]
    @test dense[AC.interval_rows(layout, interval), AC.node_range(layout, interval + 1)] == blocks.wrt_right[:, :, interval]
    @test dense[AC.interval_rows(layout, interval), parameter_columns[AC.frequency_range(layout)]] ==
          blocks.wrt_frequencies[:, :, interval]
    @test dense[AC.outer_rows(layout), parameter_columns] == blocks.outer_wrt_parameters
    @test blocks.horizon_wrt_node == I
    n_state = AC.node_block_size(layout)
    @test AC.nnz(in_place) == (layout.n_nodes - 1) * n_state * (2n_state + 2layout.n_modes) +
                              2n_state * (n_state + AC.n_parameters(layout))
    # Threaded assembly (run the tests with `julia -t 4` to exercise it) against one sequential pass of the kernels.
    frequencies = AC.component_frequencies(mode_set, AC._complexify(unknowns[parameter_columns][AC.frequency_range(layout)]))
    node_workspace = AC.NodeWorkspace{Float64}(problem)
    node_rhs = [AC.rhs_block!(zeros(n_state), node_workspace, problem, mesh.offsets[node],
                              unknowns[AC.node_range(layout, node)], frequencies, scale) for node in 1:layout.n_nodes]
    node_derivatives = [AC.rhs_derivatives!(zeros(n_state, n_state), zeros(n_state, 2layout.n_modes), node_workspace, problem,
                                            mesh.offsets[node], unknowns[AC.node_range(layout, node)], frequencies, scale)
                        for node in 1:layout.n_nodes]
    assembled_residual = AC.residual(problem, layout, unknowns, amplitudes, scale)
    for interval in 1:layout.n_nodes-1
        left, right = unknowns[AC.node_range(layout, interval)], unknowns[AC.node_range(layout, interval + 1)]
        rows = AC.interval_residual!(zeros(n_state), node_workspace, problem, interval, left, right,
                                     node_rhs[interval], node_rhs[interval + 1], frequencies, scale)
        @test rows == assembled_residual[AC.interval_rows(layout, interval)]
        wrt_left, wrt_right, wrt_frequencies = AC.interval_jacobian!(zeros(n_state, n_state), zeros(n_state, n_state),
            zeros(n_state, 2layout.n_modes), node_workspace, problem, interval, node_derivatives[interval],
            node_derivatives[interval + 1], left, right, node_rhs[interval], node_rhs[interval + 1], frequencies, scale)
        @test (wrt_left, wrt_right, wrt_frequencies) == (blocks.wrt_left[:, :, interval], blocks.wrt_right[:, :, interval],
                                                         blocks.wrt_frequencies[:, :, interval])
    end
    # Boundary blocks against differentiation of the whole boundary residual in the parameter block.
    parameter_block = unknowns[parameter_columns]
    for (side, node, block) in ((:horizon, 1, blocks.horizon_wrt_parameters), (:outer, layout.n_nodes, blocks.outer_wrt_parameters))
        node_block = unknowns[AC.node_range(layout, node)]
        autodiff = ForwardDiff.jacobian(p -> AC.boundary_residual(problem, layout, side, node_block, p, amplitudes, scale),
                                        parameter_block)
        @test maximum(abs, block .- autodiff) < 1e-12 * maximum(abs, autodiff)
    end
end

@testset "Block factorization agrees with sparse LU" begin
    background = KerrBackground(0.9)
    mesh = RadialMesh(background; rstar_min=-10.0, rstar_max=30.0, step=5.0)
    for mode_set in (ModeSet([BoundMode(2, 1, 1)]), ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)]))
        problem = RadialProblem(background, CosinePotential(), 0.3, mode_set; mesh, n_phase=8, n_polar=8)
        layout = AC.UnknownLayout(problem)
        amplitudes, scale = fill(0.05, length(mode_set.modes)), 0.05
        unknowns = AC.initial_unknowns(problem, layout, amplitudes, scale) .+ 0.01 .* randn(AC.n_unknowns(layout))
        jacobian = AC.jacobian(problem, layout, unknowns, amplitudes, scale)
        scaling = AC.NewtonScaling(problem, layout, AC._complexify(unknowns[AC.parameter_range(layout)][AC.frequency_range(layout)]))
        block = AC.ScaledFactorization(jacobian, scaling; method=:block)
        reference = AC.ScaledFactorization(jacobian, scaling; method=:sparse_lu)
        right_hand_sides = randn(AC.n_unknowns(layout), 3)
        @test block isa AC.ScaledFactorization{AC.BlockFactorization}
        @test maximum(abs, block \ right_hand_sides .- reference \ right_hand_sides) < 1e-9 * maximum(abs, reference \ right_hand_sides)
        @test block \ right_hand_sides[:, 1] ≈ reference \ right_hand_sides[:, 1] rtol=1e-9
        @test maximum(abs, jacobian * (block \ right_hand_sides) .- right_hand_sides) < 1e-8 * maximum(abs, right_hand_sides)
        # Factorising the blocks directly (what solve_radial does) gives the same factors as the sparse matrix.
        blocks = AC.jacobian_blocks!(AC.JacobianWorkspace(problem, layout), problem, layout, unknowns, amplitudes, scale)
        @test AC.ScaledFactorization(blocks, scaling) \ right_hand_sides == block \ right_hand_sides
        @test AC.ScaledFactorization(blocks, scaling; method=:sparse_lu) \ right_hand_sides == reference \ right_hand_sides
        blocks.outer_wrt_parameters[1, 2layout.n_modes + 1] = 1.0   # an outer row depending on a horizon coefficient
        @test_throws ErrorException AC.ScaledFactorization(blocks, scaling)
    end
    # Three modes and two angular channels per component, with weights and a matrix that needs row exchanges:
    # random values in every structural slot and a zeroed diagonal entry in every fifth row.
    for mode_set in (ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2), BoundMode(4, 3, 3)]; max_order=1),
                     ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)]; n_angular=2))
        problem = RadialProblem(background, CosinePotential(), 0.3, mode_set; mesh, n_phase=6, n_polar=8)
        layout = AC.UnknownLayout(problem)
        blocks = AC.JacobianWorkspace(problem, layout).blocks
        for block in AC._block_arrays(blocks)
            block .= randn(size(block))
        end
        n_frequencies, n_horizon = 2layout.n_modes, 2layout.n_channels
        blocks.horizon_wrt_parameters[:, n_frequencies + n_horizon + 1:end] .= 0      # horizon rows: no outer coefficients
        blocks.outer_wrt_parameters[:, n_frequencies + 1:n_frequencies + n_horizon] .= 0   # outer rows: no horizon coefficients
        matrix = AC.sparse_jacobian(blocks, AC.JacobianPattern(layout))
        size_total = AC.n_unknowns(layout)
        for index in 1:5:size_total
            matrix[index, index] = 0.0
        end
        row_weights, column_weights = exp.(randn(size_total)), exp.(randn(size_total))
        factors = AC.BlockFactorization(matrix, layout, row_weights, column_weights)
        scaled = Diagonal(row_weights) * matrix * Diagonal(column_weights)
        right_hand_sides = randn(size_total, 2)
        solution = factors \ right_hand_sides
        @test maximum(abs, scaled * solution .- right_hand_sides) < 1e-10 * maximum(abs, right_hand_sides)
        @test solution ≈ lu(scaled) \ right_hand_sides rtol=1e-6   # random matrices are not well conditioned
        @test any(factors.middle_pivots .!= (1:size(factors.middle_pivots, 1)))   # rows were exchanged
    end
    # An entry outside the block structure is rejected.
    layout = AC.UnknownLayout(RadialProblem(background, CosinePotential(), 0.3, ModeSet([BoundMode(2, 1, 1)]); mesh, n_phase=8, n_polar=8))
    size_total = AC.n_unknowns(layout)
    misplaced = sparse([1], [size_total - 40], [1.0], size_total, size_total)
    @test_throws ErrorException AC.stamp_panels(misplaced, layout, ones(size_total), ones(size_total))
end

@testset "Linear single mode: eigenvalue, fluxes, sensitivities" begin
    background = KerrBackground(0.9)
    mode_set = ModeSet([BoundMode(2, 1, 1)])
    mesh = RadialMesh(background; rstar_min=-60.0, rstar_max=160.0, step=0.8)
    problem = RadialProblem(background, FreePotential(), 0.3, mode_set; mesh)
    solution = solve_radial(problem, [1e-6])
    @test solution.converged
    frequency = mode_frequencies(solution)[1]
    @test real(frequency) ≈ AC.hydrogenic_frequency(BoundMode(2, 1, 1), 0.3) rtol=2e-3
    @test imag(frequency) > 0   # superradiant
    charges = cloud_charges(solution)
    losses = boundary_losses(solution)
    @test charges.angular_momentum / charges.energy ≈ 1 / real(frequency) rtol=1e-4
    @test -losses.energy ≈ 2imag(frequency) * charges.energy rtol=1e-4
    @test -losses.angular_momentum ≈ 2imag(frequency) * charges.angular_momentum rtol=1e-4
    jacobian = charge_jacobian(solution)
    @test jacobian[1, 1] ≈ 2charges.energy / 1e-6 rtol=1e-8
    @test jacobian[2, 1] ≈ 2charges.angular_momentum / 1e-6 rtol=1e-8
    # Linear problem: the amplitude only rescales the solution.
    rescaled = solve_radial(problem, [3e-6]; guess=solution)
    @test rescaled.iterations ≤ 2
    @test mode_frequencies(rescaled)[1] ≈ frequency
end


@testset "Warm starts: same point and tangent predictor" begin
    background = KerrBackground(0.9)
    mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)])
    mesh = RadialMesh(background; rstar_min=-60.0, rstar_max=250.0, step=0.8)
    problem = RadialProblem(background, CosinePotential(), 0.3, mode_set; mesh)
    amplitudes = [1e-6, 1e-7]
    solution = solve_radial(problem, amplitudes)
    @test solution.converged

    # Same system: equal inputs give equal discretisations, a different spin does not.
    @test AC._same_discretisation(problem, with_background(problem, KerrBackground(0.9), 0.3))
    @test !AC._same_discretisation(problem, with_background(problem, KerrBackground(0.9 + 1e-9), 0.3))

    # A solve at the guess's own point takes no step and keeps the guess's factorisation.
    again = solve_radial(problem, amplitudes; guess=solution)
    @test again.converged
    @test again.iterations == 1
    @test again.unknowns == solution.unknowns
    @test again.jacobian_factorization === solution.jacobian_factorization

    # Warm start with a change: profile-like unknowns rescaled, frequencies only moved.
    layout = solution.layout
    change = 1e-3 .* sin.(1:length(solution.unknowns))
    moved = AC.warm_start_unknowns(solution, 2solution.scale; change)
    frequency_indices = AC.parameter_range(layout)[AC.frequency_range(layout)]
    @test moved[frequency_indices] == (solution.unknowns .+ change)[frequency_indices]
    @test moved[1:10] ≈ (solution.unknowns .+ change)[1:10] ./ 2

    # The tangent prediction is second order in the change of amplitudes and background.
    cache = EvolutionCache(problem)
    state = BackgroundState(1.0, 0.9)
    tangents = solution_tangents(cache, solution, state)
    function prediction_errors(relative_change)
        target_amplitudes = amplitudes .* (1 .+ relative_change .* [1.0, -2.0])
        target = BackgroundState(1 - 0.3relative_change, 0.9 * (1 - relative_change))
        target_problem = AC.current_problem(cache, target)
        reference = solve_radial(target_problem, target_amplitudes; guess=solution)
        predicted = AC.warm_start_unknowns(solution, reference.scale;
                                           change=predicted_change(tangents, target_amplitudes, target))
        plain = AC.warm_start_unknowns(solution, reference.scale)
        warm = solve_radial(target_problem, target_amplitudes; guess=solution,
                            guess_change=predicted_change(tangents, target_amplitudes, target))
        return norm(predicted .- reference.unknowns, Inf), norm(plain .- reference.unknowns, Inf), warm, reference
    end
    predicted_error, plain_error, warm, reference = prediction_errors(1e-4)
    halved_error, = prediction_errors(5e-5)
    @test predicted_error < 1e-2plain_error
    @test 3 < predicted_error / halved_error < 5
    @test warm.converged
    @test norm(warm.unknowns .- reference.unknowns, Inf) < 1e-9
    @test growth_rates(warm) ≈ growth_rates(reference) rtol=1e-6
end

@testset "Wave actions: identities and balances" begin
    background = KerrBackground(0.9)
    mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)])
    mesh = RadialMesh(background; rstar_min=-60.0, rstar_max=250.0, step=0.8)
    problem = RadialProblem(background, CosinePotential(), 0.3, mode_set; mesh)
    solution = solve_radial(problem, [1e-6, 1e-7])
    frequencies = real.(mode_frequencies(solution))
    charges, actions = cloud_charges(solution), phase_actions(solution)
    losses, action_loss = boundary_losses(solution), action_losses(solution)
    # J = Σ m_i N_i holds exactly; the fluxes are exact combinations of the action fluxes.
    @test charges.angular_momentum ≈ dot([1, 2], actions) rtol=1e-12
    @test losses.energy ≈ dot(frequencies, action_loss) rtol=1e-12
    @test losses.angular_momentum ≈ dot([1, 2], action_loss) rtol=1e-12
    # Linear regime: the action balance reproduces A' = γA mode by mode.
    rates = amplitude_rates(solution; balance=:action)
    @test rates ≈ growth_rates(solution) .* solution.amplitudes rtol=1e-4
    @test amplitude_rates(solution) == rates   # the action balance is the default
    # and dE = Σ ω dN once the cloud is (nearly) contained.
    @test charge_jacobian(solution)[1, :] ≈ vec(transpose(frequencies) * action_jacobian(solution)) rtol=1e-3
    @test_throws ArgumentError amplitude_rates(solution; balance=:momentum)
end

@testset "Balance terms: one dual pass, optional energy pieces" begin
    background = KerrBackground(0.9)
    mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)])
    mesh = RadialMesh(background; rstar_min=-60.0, rstar_max=250.0, step=0.8)
    problem = RadialProblem(background, CosinePotential(), 0.3, mode_set; mesh)
    solution = solve_radial(problem, [3e-3, 1e-4])
    @test solution.converged
    state = BackgroundState(1.0, background.a)

    # Dual-profile potential average: the analytic partials equal ForwardDiff through the field reconstruction.
    channel_count = AC.n_channels(mode_set)
    profiles = [0.3 * cis(0.7index) / index for index in 1:channel_count]
    direction = [cis(1.3index) * (1 + 0.1index) for index in 1:channel_count]
    for potential in (CosinePotential(), ChiralPotential(), QuarticPotential(0.5))
        reference = ForwardDiff.derivative(epsilon -> AC.phase_averaged_potential(problem.quadrature, potential,
                                                                                   profiles .+ epsilon .* direction), 0.0)
        tagged = ForwardDiff.derivative(epsilon -> AC.averaged_potential(problem.quadrature, potential,
                                                                          profiles .+ epsilon .* direction), 0.0)
        @test tagged ≈ reference rtol=1e-12
    end

    # The threaded node integrals reproduce the node-by-node bilinear path.
    layout, unknowns, scale = solution.layout, solution.unknowns, solution.scale
    frequencies = AC.real_component_frequencies(problem, layout, unknowns)
    weights = AC.simpson_weights(problem.mesh.rstar)
    by_hand = sum(1:layout.n_nodes) do node
        node_state = AC.physical_state(layout, unknowns, scale, node)
        bilinears = AC.field_bilinears(problem, problem.mesh.offsets[node], node_state.profiles, node_state.slopes, frequencies)
        weights[node] * collect(AC.charge_densities(problem, problem.mesh.offsets[node], bilinears))
    end
    @test collect(cloud_charges(solution)) ≈ by_hand rtol=1e-13
    # Node bilinears against the complex formulas.
    node = 40
    node_state, offset = AC.physical_state(layout, unknowns, scale, node), problem.mesh.offsets[node]
    bilinears = AC.field_bilinears(problem, offset, node_state.profiles, node_state.slopes, frequencies)
    component_values = AC.component_profiles(problem.quadrature, node_state.profiles)   # [polar node, component]
    component_dr = AC.component_profiles(problem.quadrature, node_state.slopes) ./ AC.dr_drstar(background, offset)
    component_dtheta = reduce(hcat, problem.quadrature.angular.harmonics_dtheta[:, channels] * node_state.profiles[channels]
                                    for channels in mode_set.component_channels)
    azimuthal = [component.m for component in mode_set.components]
    @test bilinears.dt_dt ≈ 2abs2.(component_values) * frequencies .^ 2 rtol=1e-13
    @test bilinears.dphi_dphi ≈ 2abs2.(component_values) * azimuthal .^ 2 rtol=1e-13
    @test bilinears.dt_dphi ≈ -(2abs2.(component_values) * (frequencies .* azimuthal)) rtol=1e-13
    @test bilinears.dr_dr ≈ vec(sum(2abs2.(component_dr); dims=2)) rtol=1e-13
    @test bilinears.dtheta_dtheta ≈ vec(sum(2abs2.(component_dtheta); dims=2)) rtol=1e-13
    @test bilinears.potential ≈ AC.phase_averaged_potential(problem.quadrature, problem.potential, node_state.profiles) rtol=1e-13
    # A dual number in the spin alone (plain unknowns) propagates through the densities.
    densities_at(spin) = AC.integrated_densities(with_background(problem, KerrBackground(spin), problem.mass_coupling), layout,
                                                 unknowns, scale)
    spin_step = 1e-4   # smaller steps meet the roundoff of the near-horizon offsets
    @test ForwardDiff.derivative(densities_at, background.a) ≈
          (densities_at(background.a + spin_step) .- densities_at(background.a - spin_step)) ./ 2spin_step rtol=1e-5

    # Default action balance leaves the energy pieces out; the amplitude rates are the same with them recorded.
    cache = EvolutionCache(problem; decay_constant_ratio=1e-2)
    recorded_cache = EvolutionCache(problem; decay_constant_ratio=1e-2, record_metric_work=true)
    terms, recorded = balance_terms(cache, solution, state), balance_terms(recorded_cache, solution, state)
    @test all(isnan, terms.charge_jacobian) && isnan(terms.metric_work) && isnan(terms.derivatives.charges_wrt_mass.energy)
    @test !any(isnan, recorded.charge_jacobian) && !isnan(recorded.metric_work)
    @test amplitude_rates(terms, cache) ≈ amplitude_rates(recorded, recorded_cache) rtol=1e-12
    @test terms.derivatives.actions_wrt_mass ≈ recorded.derivatives.actions_wrt_mass rtol=1e-12
    @test terms.actions ≈ recorded.actions rtol=1e-12
    @test terms.charges.energy ≈ recorded.charges.energy rtol=1e-12
    @test_throws ArgumentError amplitude_rates(terms; balance=:energy)
    # The spindown right-hand side records the work as NaN unless asked, and reuses the charges and actions of the terms.
    for (label, rhs_cache) in (("default", EvolutionCache(problem; decay_constant_ratio=1e-2)), ("recorded", recorded_cache))
        empty!(rhs_cache.diagnostics)
        rhs_cache.accepted_solution = solution
        rates = zeros(4)
        AC._spindown_rhs!(rates, [log.(solution.amplitudes); state.mass; state.angular_momentum], rhs_cache, 0.0)
        logged = only(rhs_cache.diagnostics)
        @test isnan(logged.metric_work) == (label == "default")
        @test logged.actions ≈ terms.actions rtol=1e-10
    end

    # Exact explicit partial against the central difference (reference) and against four re-solved problems.
    exact = recorded.derivatives
    central = background_charge_derivatives(recorded_cache, solution, state; relative_step=1e-6)
    resolved = balance_terms(recorded_cache, solution, state; resolve_background=true).derivatives
    for field in fieldnames(BackgroundDerivatives)
        exact_value, central_value, resolved_value = (collect(getfield(d, field)) for d in (exact, central, resolved))
        scale_of_field = maximum(abs, exact_value)
        @test maximum(abs, exact_value .- central_value) ≤ 1e-6 * scale_of_field
        @test maximum(abs, exact_value .- resolved_value) ≤ 1e-4 * scale_of_field
    end
    @test recorded.charge_jacobian ≈ charge_jacobian(solution) .* [state.mass, state.mass^2] rtol=1e-10
    @test recorded.action_jacobian ≈ state.mass^2 .* action_jacobian(solution) rtol=1e-10
end

@testset "Evolution: logarithmic state and trajectory" begin
    background = KerrBackground(0.9)
    mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)])
    mesh = RadialMesh(background; rstar_min=-60.0, rstar_max=250.0, step=0.8)
    problem = RadialProblem(background, CosinePotential(), 0.3, mode_set; mesh)
    amplitudes = [1e-6, 1e-7]
    solution = solve_radial(problem, amplitudes)
    final_time = 0.5 / growth_rates(solution)[1]
    result = evolve_fixed_background(problem, amplitudes, (0.0, final_time); balance=:action, saveat=final_time)
    trajectory = amplitude_trajectory(result)
    @test result.log_amplitudes
    @test trajectory.times ≈ [0.0, final_time]
    @test trajectory.amplitudes[1, :] ≈ amplitudes
    @test trajectory.amplitudes[end, :] ≈ amplitudes .* exp.(growth_rates(solution) .* final_time) rtol=1e-4
    @test length(result.diagnostics) < 25
end

@testset "Finite-difference Jacobian: steps, columns, failures, reuse" begin
    settings = FiniteDifferenceJacobian(; amplitude_step=1e-3, mass_step=1e-6, angular_momentum_step=2e-6)
    @test AC._difference_steps(settings, [0.0, 0.0, 1.0, 0.9], 2, true) == [1e-3, 1e-3, 1e-6, 2e-6]
    @test AC._difference_steps(settings, [2.0, -4.0], 2, false) == [2e-3, 4e-3]
    @test AC._difference_steps(settings, [0.0, -4.0], 2, false) == [4e-3, 4e-3]   # zero amplitude: largest |A|
    @test AC._difference_steps(settings, [0.0, 0.0, 1.0, 0.9], 2, false) == [1e-3, 1e-3, 1e-6, 2e-6]
    @test_throws ArgumentError FiniteDifferenceJacobian(; difference=:backward)

    # Quadratic rates: forward differences are off by step × curvature, central ones exact.
    linear_part = [1.0 2.0; -3.0 0.5]
    rates_at(state) = linear_part * state .+ [state[1]^2, 0.0]
    base_state, steps = [0.3, -0.2], [1e-4, 1e-4]
    evaluate(perturbation) = (; rates=rates_at(base_state .+ perturbation), solution=nothing)
    exact = linear_part .+ [2base_state[1] 0.0; 0.0 0.0]
    forward, _ = AC._difference_jacobian(evaluate, base_state, rates_at(base_state), steps, FiniteDifferenceJacobian())
    central, _ = AC._difference_jacobian(evaluate, base_state, rates_at(base_state), steps,
                                         FiniteDifferenceJacobian(; difference=:central))
    @test forward ≈ exact atol=2e-4
    @test central ≈ exact atol=1e-9

    # A solve that fails for steps above 3e-5 succeeds after one reduction (1e-4 / 4); one that always fails gives nothing.
    fragile(perturbation) = maximum(abs, perturbation) > 3e-5 ? (; rates=nothing, solution=nothing) : evaluate(perturbation)
    reduced, results = AC._difference_jacobian(fragile, base_state, rates_at(base_state), steps,
                                               FiniteDifferenceJacobian(; step_reductions=1))
    @test reduced ≈ exact atol=1e-4
    @test all(result -> result.attempts == 2 && result.step == 2.5e-5, results)
    failing(perturbation) = (; rates=nothing, solution=nothing)
    @test AC._difference_jacobian(failing, base_state, rates_at(base_state), steps, settings)[1] === nothing

    # Reuse policy: exact at the computed state, stale for at most max_reuse_steps new states, never after a rejection.
    store = AC.JacobianStore(2)
    reuse_twice = FiniteDifferenceJacobian(; max_reuse_steps=2)
    @test !AC._reuse_stored_jacobian(store, reuse_twice, [0.0, 0.0])   # nothing stored yet
    store.computed_state, store.age = [0.0, 0.0], 0
    @test AC._reuse_stored_jacobian(store, reuse_twice, [0.0, 0.0])    # retry at the computed state
    @test AC._reuse_stored_jacobian(store, reuse_twice, [1.0, 0.0])    # age 1
    @test !AC._reuse_stored_jacobian(store, reuse_twice, [1.0, 0.0])   # rejected with the stale Jacobian
    store.computed_state, store.age = [1.0, 0.0], 0
    @test AC._reuse_stored_jacobian(store, reuse_twice, [2.0, 0.0])    # age 1
    @test AC._reuse_stored_jacobian(store, reuse_twice, [3.0, 0.0])    # age 2
    @test !AC._reuse_stored_jacobian(store, reuse_twice, [4.0, 0.0])   # age 3
    @test !AC._reuse_stored_jacobian(store, FiniteDifferenceJacobian(), [5.0, 0.0])
end

@testset "Spindown diagnostics: E and J once per accepted step" begin
    background = KerrBackground(0.9)
    mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)])
    mesh = RadialMesh(background; rstar_min=-60.0, rstar_max=250.0, step=0.8)
    problem = RadialProblem(background, CosinePotential(), 0.3, mode_set; mesh)
    amplitudes = [1e-6, 1e-7]
    solution = solve_radial(problem, amplitudes)
    final_time = 0.1 / growth_rates(solution)[1]
    run(; options...) = evolve_with_spindown(problem, amplitudes, (0.0, final_time), 0.01; initial_guess=solution,
                                             dt=final_time / 2, options...)
    sparse, full = run(), run(charges_every_evaluation=true)
    energies, full_energies = [entry.charges.energy for entry in sparse.diagnostics], [entry.charges.energy for entry in full.diagnostics]
    filled = findall(!isnan, energies)
    # The rates do not use E and J, so both runs take the same evaluations.
    @test sparse.ode_solution.u == full.ode_solution.u && length(energies) == length(full_energies)
    @test !any(isnan, full_energies)
    @test first(filled) == 1 && last(filled) == length(energies)
    @test length(filled) == 1 + sparse.ode_solution.stats.naccept
    @test energies[filled] ≈ full_energies[filled] rtol=1e-12
    @test [entry.charges.angular_momentum for entry in sparse.diagnostics[filled]] ≈
          [entry.charges.angular_momentum for entry in full.diagnostics[filled]] rtol=1e-12
    # Tsit5 evaluates at the accepted state last, so the filled rows are at the accepted times.
    @test [entry.time for entry in sparse.diagnostics[filled]][2:end] ≈ sparse.ode_solution.t[2:end]
end

@testset "Evolution with stiff integrators" begin
    background = KerrBackground(0.9)
    mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)])
    mesh = RadialMesh(background; rstar_min=-60.0, rstar_max=250.0, step=0.8)
    problem = RadialProblem(background, CosinePotential(), 0.3, mode_set; mesh)
    amplitudes = [1e-6, 1e-7]
    solution = solve_radial(problem, amplitudes)
    final_time = 0.1 / growth_rates(solution)[1]
    expected = amplitudes .* exp.(growth_rates(solution) .* final_time)
    for algorithm in (Rosenbrock23(), Rodas5P(), FBDF())
        result = evolve_with_spindown(problem, amplitudes, (0.0, final_time), 0.01; algorithm, initial_guess=solution,
                                      dt=final_time / 2)
        counts, stats = result.solve_counts, result.ode_solution.stats
        final_state = result.ode_solution.u[end]
        @test amplitude_trajectory(result).amplitudes[end, :] ≈ expected rtol=1e-4
        @test counts.stage == length(result.diagnostics) == stats.nf
        @test counts.failed == 0
        @test counts.jacobian == 4counts.jacobians_computed
        # The final radial solution is at the final state, whatever the integrator's last evaluation was.
        cache = AC.EvolutionCache(problem; decay_constant_ratio=0.01)
        @test AC._solution_at_state(result.final_radial_solution, final_state, cache)
        if algorithm isa FBDF
            # FBDF asks for a Jacobian only when its Newton iteration converges slowly, at its last accepted state;
            # that state's derivative comes from the iteration, so each new Jacobian costs one accepted-state solve.
            @test counts.jacobians_computed <= stats.naccept + stats.nreject
            @test counts.accepted_state <= 2 + counts.jacobians_computed
            continue
        end
        # One Jacobian per accepted state; a retried step reuses the one at its (unchanged) starting state.
        @test counts.jacobians_computed == stats.naccept
        @test counts.jacobians_reused == stats.nreject
        # Rosenbrock23 is first-same-as-last. Rodas5P needs up to one solve per accepted step (plus the initial one
        # and the final one); its last stage sits within the matching tolerance of the accepted state when the step
        # changes little.
        algorithm isa Rosenbrock23 && @test counts.accepted_state == 1
        algorithm isa Rodas5P && @test counts.accepted_state <= 2 + stats.naccept
    end
    fixed = evolve_fixed_background(problem, amplitudes, (0.0, final_time); algorithm=Rosenbrock23(), dt=final_time / 2)
    @test amplitude_trajectory(fixed).amplitudes[end, :] ≈ expected rtol=1e-4
    @test fixed.solve_counts.jacobian == 2fixed.solve_counts.jacobians_computed
    fixed_rodas = evolve_fixed_background(problem, amplitudes, (0.0, final_time); algorithm=Rodas5P(), dt=final_time / 2)
    @test AC._solution_at_state(fixed_rodas.final_radial_solution, fixed_rodas.ode_solution.u[end],
                                AC.EvolutionCache(problem))
end

"""Stand-in for an `OrdinaryDiffEq` integrator: the fields `_spindown_rhs!` reads."""
mutable struct IntegratorStub
    iter::Int
    opts::NamedTuple{(:adaptive,), Tuple{Bool}}
end

@testset "Failed radial solves are rejected per step attempt" begin
    background = KerrBackground(0.9)
    mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)])
    mesh = RadialMesh(background; rstar_min=-60.0, rstar_max=250.0, step=0.8)
    problem = RadialProblem(background, CosinePotential(), 0.3, mode_set; mesh)
    solution = solve_radial(problem, [1e-6, 1e-7])
    state_vector = [log.(solution.amplitudes); 1.0; 0.9]
    cache = AC.EvolutionCache(problem; decay_constant_ratio=0.01)
    cache.accepted_solution = solution
    stub = IntegratorStub(3, (adaptive=true,))
    cache.integrator = stub
    rates = ones(4)
    # A tolerance below roundoff makes the solve fail: zero rates, the attempt is marked and later evaluations of it skipped.
    cache.newton_tol = 1e-30
    AC._spindown_rhs!(rates, state_vector, cache, 0.0)
    @test rates == zeros(4) && cache.failed_attempt == 3 && cache.solve_counts.failed == 1
    @test AC._reject_failed_step(state_vector, cache, 0.0)
    AC._spindown_rhs!(rates, state_vector, cache, 0.0)
    @test cache.solve_counts.stage == 1
    # The next attempt solves again even if `isoutofdomain` was never called (the integrator's own Newton gave up).
    cache.newton_tol = 1e-10
    stub.iter = 4
    @test !AC._reject_failed_step(state_vector, cache, 0.0)
    AC._spindown_rhs!(rates, state_vector, cache, 0.0)
    @test cache.solve_counts.stage == 2 && all(!iszero, rates)
    # With a fixed step the failure cannot be rejected: the evolution stops.
    stub.opts = (adaptive=false,)
    stub.iter = 5
    cache.newton_tol = 1e-30
    @test_throws ErrorException AC._spindown_rhs!(rates, state_vector, cache, 0.0)
    # A trial state with negative spin is rejected without a radial solve (and without a DomainError).
    stub.opts = (adaptive=true,)
    stub.iter = 6
    cache.newton_tol = 1e-10
    solves_before = total_solves(cache.solve_counts)
    negative_spin_state = vcat(state_vector[1:2], state_vector[3], -0.1)
    AC._spindown_rhs!(rates, negative_spin_state, cache, 0.0)
    @test rates == zeros(4) && cache.failed_attempt == 6 && total_solves(cache.solve_counts) == solves_before
    @test AC._reject_failed_step(negative_spin_state, cache, 0.0)
end

@testset "solve_radial control flow and floor_tol" begin
    # Decisions at a full step that does not reduce the residual, and after an accepted step.
    @test AC._floor_action(false, 1e-12, 1e-4, 1e-10, 1e-10) === nothing                # step above √tol: backtrack
    @test AC._floor_action(true, 1e-12, 1e-11, 1e-10, 1e-10) == :converged              # chord at its floor
    @test AC._floor_action(true, 1e-12, 1e-7, 1e-10, 1e-10) == :refactorise             # chord step above tol
    @test AC._floor_action(true, 1e-9, 1e-11, 1e-10, 1e-8) == :refactorise              # floor_tol is for Newton steps
    @test AC._floor_action(false, 1e-12, 1e-7, 1e-10, 1e-10) == :converged
    @test AC._floor_action(false, 1e-9, 1e-11, 1e-10, 1e-10) == :stuck                  # floor above tol
    @test AC._floor_action(false, 1e-9, 1e-11, 1e-10, 1e-8) == :converged               # ... below floor_tol
    @test AC._floor_action(false, 1e-9, 1e-7, 1e-10, 1e-10) === nothing
    @test AC._convergence_status(1e-11, 1e-11, 1.0, 0, false, 1e-10) == (:converged, 0)
    @test AC._convergence_status(1e-11, 1e-9, 1e-9, 0, true, 1e-10) == (:refactorise, 0)
    @test AC._convergence_status(1e-11, 1e-9, 1e-9, 1, false, 1e-10) == (:continue, 2)
    @test AC._convergence_status(1e-11, 1e-9, 1e-9, 2, false, 1e-10) == (:converged, 3)
    @test AC._convergence_status(1e-11, 1e-4, 1e-4, 2, false, 1e-10) == (:ill_conditioned, 3)
    @test AC._convergence_status(1e-9, 1e-9, 1e-9, 2, false, 1e-10) == (:continue, 0)

    background = KerrBackground(0.9)
    mode_set = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)])
    mesh = RadialMesh(background; rstar_min=-60.0, rstar_max=250.0, step=0.8)
    problem = RadialProblem(background, CosinePotential(), 0.3, mode_set; mesh)
    amplitudes = [3e-3, 1e-4]
    solution = solve_radial(problem, amplitudes)
    @test solution.converged
    layout = solution.layout
    target = amplitudes .* [1.0, 1.3]

    # A bad prediction falls back to the plain warm start: same iterations and unknowns as no prediction.
    plain = solve_radial(problem, target; guess=solution)
    bad_change = 10 .* solution.unknowns   # profiles eleven times too large, parameters unchanged
    bad_change[AC.parameter_range(layout)] .= 0
    bad = solve_radial(problem, target; guess=solution, guess_change=bad_change)
    @test bad.converged && bad.iterations == plain.iterations && bad.unknowns == plain.unknowns

    # Chord steps from the guess's factorisation reach the same solution as Newton steps; after a large change a chord
    # step is rejected and the Jacobian refactorised.
    log_path, log_io = mktemp()
    changed = amplitudes .* [1.5, 1.0]
    chord = redirect_stdout(() -> solve_radial(problem, changed; guess=solution, verbose=true), log_io)
    close(log_io)
    newton = solve_radial(problem, changed; guess=solution, reuse_factorization=false)
    @test chord.converged && newton.converged
    @test maximum(abs, chord.unknowns .- newton.unknowns) < 1e-8 * maximum(abs, newton.unknowns)
    @test occursin("rejected", read(log_path, String))

    # A tolerance below the residual floor: the solve stops unconverged, unless floor_tol accepts the floor; the
    # evolution cache passes floor_tol to its solves.
    unreachable = 1e-17
    @test !(@test_logs (:warn,) match_mode=:any solve_radial(problem, amplitudes; guess=solution, tol=unreachable, maxiter=4)).converged
    @test solve_radial(problem, amplitudes; guess=solution, tol=unreachable, floor_tol=1e-9).converged
    floor_cache = EvolutionCache(problem; newton_tol=unreachable, newton_floor_tol=1e-9, last_solution=solution)
    @test floor_cache.newton_floor_tol == 1e-9 && EvolutionCache(problem).newton_floor_tol == EvolutionCache(problem).newton_tol
    @test AC._solve_cached!(floor_cache, problem, amplitudes).converged
end

@testset "Threaded kernels and concurrent evaluations (subprocesses with 1 and 4 threads)" begin
    # The suite itself may run on one thread; these subprocesses run the chunked and spawned paths and compare.
    script = joinpath(@__DIR__, "threads_check.jl")
    outputs = map((1, 4)) do threads
        path = tempname()
        run(`$(Base.julia_cmd()) --project=$(Base.active_project()) -t $threads $script $path`)
        reinterpret(Float64, read(path))
    end
    @test outputs[1][1] == 1.0 && outputs[2][1] == 1.0   # concurrent evaluations equal sequential ones
    @test outputs[1] == outputs[2]                         # results do not depend on the thread count
end

end
