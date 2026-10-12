# Diagnose the absolute shooting-residual floor before changing acceptance.
include(joinpath(@__DIR__, "../../../../src/adiabatic/load.jl"))
using LinearAlgebra
alpha = 0.014966276123144448
spin = .99
grid = angular_grid(12, 25)
polar = build_polar_values(7, grid.cosine_theta_nodes)
couplings = build_cosine_squared_couplings(7, polar, grid)
frequencies = ComplexF64[]
for step in (.25, .125)
    solution = solve_linear_mode(complex(alpha * (1-alpha^2/8)), spin, alpha,
        [1, 3, 5], 1, couplings; step, inner=-100., outer=200904., match=40.,
        tolerance=1e-3, difference_step=1e-3 * alpha^3)
    modes = [(1,1), (3,1), (5,1)]
    rhs! = (derivative, state, position) -> linear_tortoise_rhs!(
        derivative, state, position, spin, alpha, solution.frequency, modes, couplings)
    horizon = horizon_boundary_data(-100., spin, alpha, solution.frequency,
        [1,3,5], 1, couplings, solution.horizon_amplitudes)
    state = integrate_radial(rhs!, vcat(horizon.values, horizon.tortoise_slopes),
        -100., 40.; step)
    println((; step, frequency=solution.frequency, residual=solution.residual,
             matching_state=state, relative_residual=solution.residual_norm/norm(state)))
    push!(frequencies, solution.frequency)
    flush(stdout)
end
println("REAL_RELATIVE_CHANGE=", abs(real(frequencies[2])/real(frequencies[1])-1))
println("IMAG_RELATIVE_CHANGE=", abs(imag(frequencies[2])/imag(frequencies[1])-1))
