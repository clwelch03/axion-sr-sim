# Run with `julia --project=. src/adiabatic/continuation.jl` from the repository root.
# Solve a small-amplitude seed, then continue through six increasing amplitudes.
# Overwrite runs/adiabatic/continuation.jls with the seed and each completed step.
using Serialization

include("angular.jl")
include("nonlinear.jl")
include("tortoise.jl")
include("radial.jl")
include("boundary_series.jl")
include("boundaries.jl")
include("integration.jl")
include("matching.jl")

output_path = joinpath(@__DIR__, "../../runs/adiabatic/continuation.jls")
mkpath(dirname(output_path))

# Build the angular grid and couplings from angular.jl and radial.jl
grid = angular_grid(16, 31)
polar_values = build_polar_values(5, grid.cosine_theta_nodes)
azimuthal_values = build_azimuthal_values(5, grid.azimuth_nodes)
couplings = build_cosine_squared_couplings(5, polar_values, grid)
spin, alpha = 0.99, 0.42
modes = [(1, 1), (3, 1), (5, 1), (3, 3), (5, 3), (5, 5)]

# Find the linear fundamental mode, normalized to outer amplitude 1.
linear = solve_linear_mode(
    0.4088 + 1e-7im, spin, alpha, [1, 3, 5], 1, couplings,
)

# Scale its m=1 amplitudes; initialize the higher harmonics to zero.
amplitude = 1e-4
horizon_seed = amplitude .* vcat(
    linear.horizon_amplitudes, zeros(ComplexF64, 3),
)
outer_seed = amplitude .* vcat(
    linear.outer_amplitudes, zeros(ComplexF64, 3),
)

# Solve the nonlinear problem at this small amplitude.
nonlinear = solve_nonlinear_mode(
    linear.frequency, amplitude, spin, alpha, modes, couplings,
    polar_values, azimuthal_values, grid;
    horizon_guess=horizon_seed,
    outer_tail_guess=outer_seed[2:end],
    tolerance=1e-11,
)

# Keep every solution, including its boundary amplitudes and diagnostics.
starting_amplitude = real(nonlinear.outer_amplitudes[1])
continuation = [(amplitude=starting_amplitude, solution=nonlinear)]
# Save the seed too, so it survives if the first continuation step fails.
serialize(output_path, (; spin, alpha, modes, linear, continuation))

for next_amplitude in starting_amplitude .* 1.5 .^ (1:6)
    previous = last(continuation)
    amplitude_ratio = next_amplitude / previous.amplitude

    # Rescale the previous solution to construct the next initial guess.
    # Newton then adjusts each harmonic independently.
    solution = solve_nonlinear_mode(
        previous.solution.frequency,
        next_amplitude, spin, alpha, modes, couplings,
        polar_values, azimuthal_values, grid;
        horizon_guess=amplitude_ratio .*
            previous.solution.horizon_amplitudes,
        outer_tail_guess=amplitude_ratio .*
            previous.solution.outer_amplitudes[2:end],
        tolerance=1e-11,
    )

    push!(continuation, (amplitude=next_amplitude, solution=solution))
    serialize(output_path, (; spin, alpha, modes, linear, continuation))

    frequency_shift = real(solution.frequency - linear.frequency)
    println(
        "A0=$next_amplitude, " *
        "shift/A0²=$(frequency_shift / next_amplitude^2), " *
        "residual=$(solution.scaled_residual_norm), " *
        "iterations=$(solution.iterations)"
    )
    flush(stdout)
end
