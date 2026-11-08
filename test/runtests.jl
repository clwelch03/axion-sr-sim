using Test
using AxionCloud
using ForwardDiff, LinearAlgebra, FastGaussQuadrature
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
    estimate = domain_estimate(background, 0.3, mode_set)
    @test estimate.slowest_decay ≈ 0.3^2 / 3 rtol=0.05           # κ₃₂₂ ≈ α²/n
    @test 0 < estimate.smallest_wavenumber < 0.03                   # the near-threshold (−1, 2) component
    @test estimate.capped                                           # its series would need r ≫ the usable radius
    @test estimate.normalization_radius < 1000                      # |211⟩ is negligible at the boundary beyond this
    @test estimate.asymptotic_radius ≈ estimate.normalization_radius  # so the normalization limit sets the outer endpoint
    @test offset_from_tortoise(background, estimate.rstar_min) ≈ 1e-6 * background.horizon_gap rtol=1e-6
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

end
