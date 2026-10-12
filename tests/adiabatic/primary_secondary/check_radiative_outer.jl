# julia --project=. tests/adiabatic/primary_secondary/check_radiative_outer.jl
# Checks of radiative_outer_boundary_state (piece 5c): the m = 3 field's far outer leg,
# at the smallest and largest ω₃ on the saved Mμ = 0.42 and 0.3 branches.
#   1. Start from the outgoing series at r* = 5000, integrate inward to 2500, and
#      compare with the series evaluated there. The difference should be RK4
#      error, falling by 16 per halving of the step.
#   2. The state at the primary's outer boundary (r* = 150 at Mμ = 0.42, 300 at
#      0.3) should not depend on the far start or the series
#      order (the amplitudes are defined at infinity), only on the step. The
#      last line per case is what the outgoing series used directly at that
#      boundary would have given.
using Printf, LinearAlgebra
include(joinpath(@__DIR__, "../../../src/adiabatic/load.jl"))

grid = angular_grid(16, 31)
polar = build_polar_values(9, grid.cosine_theta_nodes)
couplings = build_cosine_squared_couplings(9, polar, grid)
modes = [(3, 3), (5, 3), (7, 3)]
spin = 0.99
unit_amplitudes = ([1.0 + 0im, 0, 0], [0, 1.0 + 0im, 0], [0, 0, 1.0 + 0im])
series_state(position, alpha, frequency, amplitudes, order) = boundary_state(
    :outer, position, spin, alpha, frequency, modes, couplings, amplitudes;
    fundamental_m=3, outer_order=order)

# Smallest and largest ω₃ = 2 Re ω₂ - Re ω₁ on each saved branch, and the
# primary's outer boundary there.
cases = ((0.42, 0.422331, 150.0), (0.42, 0.426645, 150.0),
         (0.30, 0.300600, 300.0), (0.30, 0.300859, 300.0))

println("Check 1: series at 5000 integrated to 2500 vs series at 2500 (order 20)")
for (alpha, frequency, _) in cases
    for step in (4.0, 2.0, 1.0, 0.5)
        differences = map(unit_amplitudes) do amplitudes
            integrated = radiative_outer_boundary_state(2500.0, spin, alpha, frequency, modes, couplings,
                amplitudes; far_position=5000.0, far_step=step, outer_order=20)
            direct = series_state(2500.0, alpha, frequency, amplitudes, 20)
            norm(integrated - direct) / norm(direct)
        end
        @printf("alpha=%.2f omega3=%.6f step %.1f: %.2e %.2e %.2e\n",
                alpha, frequency, step, differences...)
    end
end

println("\nCheck 2: state at the outer boundary, relative change from the default (far 5000, order 16, step 0.25)")
for (alpha, frequency, boundary) in cases
    amplitudes = unit_amplitudes[1]
    reference = radiative_outer_boundary_state(boundary, spin, alpha, frequency, modes, couplings, amplitudes)
    for (far_position, order, step) in ((2500.0, 16, 0.25), (10000.0, 16, 0.25),
                                        (5000.0, 12, 0.25), (5000.0, 20, 0.25),
                                        (5000.0, 16, 0.5), (5000.0, 16, 0.0625))
        state = radiative_outer_boundary_state(boundary, spin, alpha, frequency, modes, couplings, amplitudes;
                                      far_position, far_step=step, outer_order=order)
        @printf("alpha=%.2f omega3=%.6f far %6.0f order %2d step %.4f: %.2e\n",
                alpha, frequency, far_position, order, step,
                norm(state - reference) / norm(reference))
    end
    naive = series_state(boundary, alpha, frequency, amplitudes, 7)
    @printf("alpha=%.2f omega3=%.6f series directly at %3.0f (order 7):        %.2e\n",
            alpha, frequency, boundary, norm(naive - reference) / norm(reference))
end
