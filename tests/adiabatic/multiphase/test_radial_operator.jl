# Load the existing single-cloud assembly only for the comparison below.
include(joinpath(@__DIR__, "../../../src/adiabatic/single_cloud/load.jl"))

@testset "One-component radial operator" begin
    grid = angular_grid(16, 31)
    polar = build_polar_values(5, grid.cosine_theta_nodes)
    couplings = build_cosine_squared_couplings(5, polar, grid)
    frequencies, cloud_m = [0.405 + 1e-7im, 0.415 + 2e-7im], [1, 2]
    position, spin, alpha = 4.0, 0.99, 0.42
    for (q, l_values) in (((1, 0), [1, 3, 5]), ((0, 1), [2, 4]), ((-1, 2), [3, 5]))
        component = fourier_component(q, l_values, cloud_m)
        numbers = component_frequency_and_m(q, frequencies, cloud_m)
        operator = component_tortoise_operator(position, spin, alpha,
            component, frequencies, cloud_m, couplings)
        count = length(l_values)
        state = ComplexF64[0.2index + 0.03im * index^2 for index in 1:2count]
        reference = similar(state)
        linear_tortoise_rhs!(reference, state, position, spin, alpha,
            numbers.frequency, [(l, numbers.m) for l in l_values], couplings;
            fundamental_m=numbers.m)
        @test -operator.slope_coefficient * state[count+1:end] -
            operator.profile_matrix * state[1:count] ≈ reference[count+1:end]
    end

    # Independent Schwarzschild m=0 answer, including nonzero frequency.
    component = fourier_component((2, -1), [0, 2, 4], cloud_m)
    operator = component_tortoise_operator(position, 0.0, alpha,
        component, real.(frequencies), cloud_m, couplings)
    radius = tortoise_geometry(position, 0.0).radius
    lapse = 1 - 2/radius
    frequency = 0.395
    expected = [frequency^2 - lapse*(alpha^2 + l*(l+1)/radius^2)
        for l in component.l_values]
    @test operator.slope_coefficient ≈ 2lapse/radius
    @test operator.profile_matrix ≈ Diagonal(expected)

    # Under (omega,m) -> (-omega,-m), the linear coefficients are unchanged.
    positive = fourier_component((1, 0), [1, 3], cloud_m)
    negative = fourier_component((-1, 0), [1, 3], cloud_m)
    positive_operator = component_tortoise_operator(position, spin, alpha,
        positive, real.(frequencies), cloud_m, couplings)
    negative_operator = component_tortoise_operator(position, spin, alpha,
        negative, real.(frequencies), cloud_m, couplings)
    @test positive_operator.profile_matrix ≈ negative_operator.profile_matrix

    # Relabeling angular entries must permute both matrix axes.
    reordered = fourier_component((1, 0), [3, 1], cloud_m)
    reordered_operator = component_tortoise_operator(position, spin, alpha,
        reordered, real.(frequencies), cloud_m, couplings)
    @test reordered_operator.profile_matrix ≈ positive_operator.profile_matrix[[2, 1], [2, 1]]
end
