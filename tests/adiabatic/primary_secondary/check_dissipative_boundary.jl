# julia --project=. tests/adiabatic/primary_secondary/check_dissipative_boundary.jl
# Check of dissipative_boundary_state (piece 4c, docs/multimode-piece4-guide.md).
# Start the source-free m = 0 field from its asymptotic series at a domain end
# (r_* = -100 or 150), with unit amplitude in one l; integrate 10M into the
# domain; compare with the series evaluated there. If the series satisfies the
# radial equation, the difference is RK4 error and falls by about 16 per halving
# of the step. At large A1, rows without the primary's background separate the
# part due to the primary's tail, which the linear series omits.
using Printf, Serialization, LinearAlgebra
include(joinpath(@__DIR__, "../../../src/adiabatic/primary_secondary/load.jl"))

primary_path = joinpath(@__DIR__, "../../../runs/adiabatic/fig2_20260922/results/cosine_branch.jls")
secondary_path = joinpath(@__DIR__, "../../../runs/adiabatic/multimode_secondary/results/cosine_secondary.jls")

setup = load_primary_secondary_setup(primary_path)
(; primary, secondary_modes, grid, polar, azimuthal, couplings) = setup
secondary_branch = deserialize(secondary_path)
dissipative_modes = [(0, 0), (2, 0), (4, 0)]
spin, alpha = primary.spin, primary.alpha
offset = 2length(secondary_modes)
(; inner, outer) = primary.options
zero_primary = Dict(mode => 0.0im for mode in primary.modes)
unit_amplitudes = ([1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0])

function series_difference(primary_lookup, frequency, secondary_frequency, side,
                           start, stop, step, amplitudes)
    rhs! = (derivative, state, position) -> dissipative_tortoise_rhs!(
        derivative, state, position, primary_lookup(position), spin, alpha,
        secondary_frequency, frequency, secondary_modes, dissipative_modes, couplings,
        polar, azimuthal, grid; include_source=false,
    )
    initial = dissipative_boundary_state(
        side, start, spin, alpha, frequency, dissipative_modes, couplings, amplitudes)
    integrated = integrate_radial(
        rhs!, vcat(zeros(ComplexF64, offset), initial), start, stop; step)[offset+1:end]
    series = dissipative_boundary_state(
        side, stop, spin, alpha, frequency, dissipative_modes, couplings, amplitudes)
    return norm(integrated - series) / norm(series)
end

println("|integrated - series| / |series| for unit amplitude in l = 0, 2, 4")
let secondary_frequency = real(secondary_branch.linear.frequency),
    frequency = 2real(primary.continuation[1].solution.frequency) - secondary_frequency
    for (side, start, stop) in ((:horizon, inner, inner + 10), (:outer, outer, outer - 10))
        for step in (0.4, 0.2, 0.1)
            differences = [series_difference(_ -> zero_primary, frequency, secondary_frequency,
                                             side, start, stop, step, amplitudes)
                           for amplitudes in unit_amplitudes]
            @printf("A1=0     %-8s %6.1f -> %6.1f step %.1f: %.2e %.2e %.2e\n",
                    side, start, stop, step, differences...)
        end
    end
end

for index in (100, 232)
    (; profiles) = primary_profiles_at(setup, index)
    secondary_frequency = real(secondary_branch.rows[index].solution.frequency)
    frequency = 2real(primary.continuation[index].solution.frequency) - secondary_frequency
    amplitude = primary.continuation[index].amplitude
    for (side, start, stop, branch) in ((:horizon, inner, inner + 10, profiles.from_horizon),
                                        (:outer, outer, outer - 10, profiles.from_outer))
        lookup = position -> primary_values_at(branch, position, primary.modes)
        for step in (0.4, 0.2)
            differences = [series_difference(lookup, frequency, secondary_frequency,
                                             side, start, stop, step, amplitudes)
                           for amplitudes in unit_amplitudes]
            @printf("A1=%-5.3g %-8s %6.1f -> %6.1f step %.1f: %.2e %.2e %.2e\n",
                    amplitude, side, start, stop, step, differences...)
        end
        differences = [series_difference(_ -> zero_primary, frequency, secondary_frequency,
                                         side, start, stop, 0.2, amplitudes)
                       for amplitudes in unit_amplitudes]
        @printf("         (same omega0, no background)        step 0.2: %.2e %.2e %.2e\n",
                differences...)
    end
end
