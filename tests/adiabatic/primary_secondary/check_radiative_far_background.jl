# julia --project=. tests/adiabatic/primary_secondary/check_radiative_far_background.jl
# Justifies leaving the primary's averaged background out of the m = 3 field's
# far outer leg (radiative_outer_boundary_state, piece 5c).
#
# Beyond its outer boundary the primary is its own outer series with the saved
# outer amplitudes (the condition its solve imposed). Evaluate it at every RK4
# stage of the far leg, add <1 - cos ϕ₁> to the m = 3 equation, and compare the
# state at the primary's outer boundary with and without it, for unit amplitude
# in each l. A change inside the span of the three outer solutions only
# relabels the outer amplitudes; the part outside the span is what can change
# the solution. Takes about a minute.
using Printf, LinearAlgebra
include(joinpath(@__DIR__, "../../../src/adiabatic/primary_secondary/load.jl"))

modes = [(3, 3), (5, 3), (7, 3)]
far_position = 5000.0
# (primary branch, point, ω₃): Mμ = 0.42 at the largest ω₃ and at the largest
# radiation tail (A1 = 13); Mμ = 0.3 at the end of its branch.
cases = (("runs/adiabatic/fig2_20260922/results/cosine_branch.jls", 60, 0.426645),
         ("runs/adiabatic/fig2_20260922/results/cosine_branch.jls", 232, 0.422331),
         ("runs/adiabatic/multimode_secondary/results/alpha030_cosine_branch.jls", 83, 0.300859))

for (path, index, frequency) in cases
    setup = load_secondary_setup(joinpath(@__DIR__, "../..", path))
    (; primary, grid) = setup
    polar = build_polar_values(9, grid.cosine_theta_nodes)
    azimuthal = build_azimuthal_values(9, grid.azimuth_nodes)
    couplings = build_cosine_squared_couplings(9, polar, grid)
    spin, alpha = primary.spin, primary.alpha
    solution = primary.continuation[index].solution
    boundary = primary.options.outer
    count = length(primary.modes)
    primary_from_series(position) = Dict(zip(primary.modes, boundary_state(
        :outer, position, spin, alpha, solution.frequency, primary.modes, couplings,
        solution.outer_amplitudes)[1:count]))

    # The series at the boundary must be the saved primary's first outer sample.
    (; profiles) = primary_profiles_at(setup, index)
    saved = profiles.from_outer.states[1:count, 1]
    series = [primary_from_series(boundary)[mode] for mode in primary.modes]
    @printf("\nMmu=%.2f A1=%.4g omega3=%.6f boundary r*=%.0f  |series - saved primary|/|saved| = %.1e\n",
            alpha, primary.continuation[index].amplitude, frequency, boundary,
            norm(series - saved) / norm(saved))
    for position in (boundary, 2boundary, 1000.0, far_position)
        deficit = averaged_mass_deficit(primary_from_series(position), polar, azimuthal, grid)
        @printf("   max <1 - cos phi1> at r*=%6.0f: %.2e\n", position, maximum(abs, deficit))
    end

    with_background(amplitudes) = integrate_radial(
        (derivative, state, position) -> secondary_tortoise_rhs!(
            derivative, state, position, primary_from_series(position), spin, alpha,
            frequency, modes, couplings, polar, azimuthal, grid),
        boundary_state(:outer, far_position, spin, alpha, frequency, modes, couplings,
                       amplitudes; fundamental_m=3, outer_order=16),
        far_position, boundary; step=0.25)

    unit_amplitudes = ([1.0 + 0im, 0, 0], [0, 1.0 + 0im, 0], [0, 0, 1.0 + 0im])
    without = hcat([radiative_outer_boundary_state(boundary, spin, alpha, frequency, modes,
                                                   couplings, amplitudes; far_position)
                    for amplitudes in unit_amplitudes]...)
    with = hcat([with_background(amplitudes) for amplitudes in unit_amplitudes]...)
    basis = Matrix(qr(without).Q)
    for column in 1:3
        change = with[:, column] - without[:, column]
        outside = change - basis * (basis' * change)
        @printf("   unit amplitude l=%d: |change|/|state| = %.2e, outside span = %.2e\n",
                modes[column][1], norm(change) / norm(without[:, column]),
                norm(outside) / norm(without[:, column]))
    end
end
