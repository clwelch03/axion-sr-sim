# julia --project=. runs/adiabatic/multimode_secondary/scripts/prepare_alpha030_seed.jl
# Seed for the Mμ = 0.3 cosine primary (Omiya et al. 2211.01949, Figs. 6 right, 7):
# the linear l = m = 1 mode scaled to the first amplitude of the Mμ = 0.42 branch,
# with zero harmonics. continue_branch.jl re-solves it as its first point because
# the seed's potential label differs from the requested one.
using Serialization
include(joinpath(@__DIR__, "../../../../src/adiabatic/load.jl"))

spin, alpha = 0.99, 0.3
modes = [(1, 1), (3, 1), (5, 1), (3, 3), (5, 3), (5, 5)]
n_polar, n_azimuth = 16, 31
options = (step=0.1, inner=-100.0, outer=300.0, match=5.0)
amplitude = 0.0011390625

grid = angular_grid(n_polar, n_azimuth)
polar = build_polar_values(5, grid.cosine_theta_nodes)
couplings = build_cosine_squared_couplings(5, polar, grid)
linear = solve_linear_mode(alpha * (1 - alpha^2 / 8) + 0im, spin, alpha, [1, 3, 5], 1, couplings;
                           options..., tolerance=1e-8)
println("linear omega1 = ", linear.frequency)

zero_harmonics = zeros(ComplexF64, 3)
solution = (; frequency=linear.frequency,
              horizon_amplitudes=vcat(amplitude .* linear.horizon_amplitudes, zero_harmonics),
              outer_amplitudes=vcat(amplitude .* linear.outer_amplitudes, zero_harmonics))
seed = (; spin, alpha, modes, options, n_polar, n_azimuth, potential="linear_seed",
          continuation=[(; amplitude, solution)], rows=[])
path = joinpath(@__DIR__, "../seeds/alpha030_linear_seed.jls")
serialize(path, seed)
println("wrote ", path)
