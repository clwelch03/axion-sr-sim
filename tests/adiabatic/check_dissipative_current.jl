# julia --project=. tests/adiabatic/check_dissipative_current.jl
# Check of dissipative_tortoise_rhs! (piece 4b, docs/multimode-piece4-guide.md).
# With the source off, the m = 0 field's radial current
#     (r^2 + a^2) Σ_l Im(conj(R_l) Q_l)
# is exactly constant in r_* for any solution (real ω₀; real symmetric operator).
# Integrate random data from r_* = 5 outward and inward and report the drift,
# relative to (r^2 + a^2) Σ_l |R_l||Q_l| (the size of the terms that cancel).
# The drift should fall by about 2^4 = 16 per halving of the step.
using Printf, Random, Serialization
include(joinpath(@__DIR__, "../../src/adiabatic/drivers/continue_secondary.jl"))

primary_path = joinpath(@__DIR__, "../../runs/adiabatic/fig2_20260922/results/cosine_branch.jls")
secondary_path = joinpath(@__DIR__, "../../runs/adiabatic/multimode_secondary/results/cosine_secondary.jls")

setup = load_secondary_setup(primary_path)
(; primary, secondary_modes, grid, polar, azimuthal, couplings) = setup
secondary_branch = deserialize(secondary_path)
dissipative_modes = [(0, 0), (2, 0), (4, 0)]
spin, alpha = primary.spin, primary.alpha
count = length(dissipative_modes)
offset = 2length(secondary_modes)

Random.seed!(1)
initial_dissipative = randn(ComplexF64, 2count)

function current_drift(primary_lookup, dissipative_frequency, secondary_frequency,
                       start, stop, step)
    rhs! = (derivative, state, position) -> dissipative_tortoise_rhs!(
        derivative, state, position, primary_lookup(position), spin, alpha,
        secondary_frequency, dissipative_frequency, secondary_modes, dissipative_modes,
        couplings, polar, azimuthal, grid; include_source=false,
    )
    initial = vcat(zeros(ComplexF64, offset), initial_dissipative)
    path = integrate_radial(rhs!, initial, start, stop; step, save_profile=true)

    drifts = map(eachindex(path.positions)) do sample
        state = path.states[:, sample]
        radial_factor = tortoise_geometry(path.positions[sample], spin).radial_factor
        current = radial_factor * sum(
            imag(conj(state[offset+i]) * state[offset+count+i]) for i in 1:count)
        scale = radial_factor * sum(
            abs(state[offset+i]) * abs(state[offset+count+i]) for i in 1:count)
        (current, scale)
    end
    initial_current = first(first(drifts))
    return maximum(abs(current - initial_current) / scale for (current, scale) in drifts)
end

println("A1 = 0 (no background), linear frequencies")
zero_primary = Dict(mode => 0.0im for mode in primary.modes)
secondary_frequency = real(secondary_branch.linear.frequency)
dissipative_frequency = 2real(primary.continuation[1].solution.frequency) - secondary_frequency
@printf("omega0 = %.6f\n", dissipative_frequency)
for stop in (61.0, -59.0)
    drifts = [current_drift(_ -> zero_primary, dissipative_frequency, secondary_frequency,
                            5.0, stop, step) for step in (0.4, 0.2, 0.1)]
    @printf("  5 -> %5.1f  drift at steps 0.4, 0.2, 0.1: %.3e %.3e %.3e  ratios %.1f %.1f\n",
            stop, drifts..., drifts[1] / drifts[2], drifts[2] / drifts[3])
end

println("\nOn saved primary points (steps 0.4 and 0.2 read saved primary samples)")
for index in (100, 232)
    local secondary_frequency, dissipative_frequency
    (; profiles) = primary_profiles_at(setup, index)
    secondary_frequency = real(secondary_branch.rows[index].solution.frequency)
    dissipative_frequency =
        2real(primary.continuation[index].solution.frequency) - secondary_frequency
    for (branch, stop) in ((profiles.from_outer, 61.0), (profiles.from_horizon, -59.0))
        lookup = position -> primary_values_at(branch, position, primary.modes)
        drifts = [current_drift(lookup, dissipative_frequency, secondary_frequency,
                                5.0, stop, step) for step in (0.4, 0.2)]
        @printf("  A1=%-7.4g omega0=%.6f  5 -> %5.1f  drift %.3e %.3e  ratio %.1f\n",
                primary.continuation[index].amplitude, dissipative_frequency, stop,
                drifts..., drifts[1] / drifts[2])
    end
end
