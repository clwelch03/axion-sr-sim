# julia --project=. tests/adiabatic/primary_secondary/check_radiative_rhs.jl
# Checks of radiative_source and radiative_tortoise_rhs! (piece 5b) on the saved
# Mμ = 0.42 cosine primary and its averaged secondary.
#   1. With the source off, the m = 3 field's radial current
#      (r^2 + a^2) Σ_l Im(conj(R_l) Q_l) is constant in r_* for any solution
#      (real ω₃, real symmetric operator). Drift relative to
#      (r^2 + a^2) Σ_l |R_l||Q_l| should fall by about 16 per halving of the step.
#      Inward the ratio is 4, 11 at steps 0.4 -> 0.1 and reaches 15.5 only by
#      step 0.0125: near the horizon the m = 3 wave has |ω₃ - 3Ω_H| ≈ 0.88.
#   2. The source is even in x = cos θ, so its projection onto even-l, m = 3
#      harmonics (odd in x) vanishes; report the share in each odd l.
using Printf, Random, Serialization, LinearAlgebra
include(joinpath(@__DIR__, "../../../src/adiabatic/primary_secondary/load.jl"))

primary_path = joinpath(@__DIR__, "../../../runs/adiabatic/fig2_20260922/results/cosine_branch.jls")
secondary_path = joinpath(@__DIR__, "../../../runs/adiabatic/multimode_secondary/results/cosine_secondary.jls")

setup = load_secondary_setup(primary_path)
(; primary, secondary_modes, grid) = setup
# Tables up to l = 9, for the m = 3 modes and the parity check.
maximum_l = 9
polar = build_polar_values(maximum_l, grid.cosine_theta_nodes)
azimuthal = build_azimuthal_values(maximum_l, grid.azimuth_nodes)
couplings = build_cosine_squared_couplings(maximum_l, polar, grid)

secondary_branch = deserialize(secondary_path)
radiative_modes = [(3, 3), (5, 3), (7, 3)]
spin, alpha = primary.spin, primary.alpha
count = length(radiative_modes)
offset = 2length(secondary_modes)

Random.seed!(1)
initial_radiative = randn(ComplexF64, 2count)

function current_drift(primary_lookup, radiative_frequency, secondary_frequency,
                       start, stop, step)
    rhs! = (derivative, state, position) -> radiative_tortoise_rhs!(
        derivative, state, position, primary_lookup(position), spin, alpha,
        secondary_frequency, radiative_frequency, secondary_modes, radiative_modes,
        couplings, polar, azimuthal, grid; include_source=false,
    )
    initial = vcat(zeros(ComplexF64, offset), initial_radiative)
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

println("Check 1: current drift with the source off")
println("A1 = 0 (no background), linear frequencies")
zero_primary = Dict(mode => 0.0im for mode in primary.modes)
let secondary_frequency = real(secondary_branch.linear.frequency),
    radiative_frequency = 2secondary_frequency - real(primary.continuation[1].solution.frequency)
    @printf("omega3 = %.6f\n", radiative_frequency)
    for stop in (61.0, -59.0)
        drifts = [current_drift(_ -> zero_primary, radiative_frequency, secondary_frequency,
                                5.0, stop, step) for step in (0.4, 0.2, 0.1)]
        @printf("  5 -> %5.1f  drift at steps 0.4, 0.2, 0.1: %.3e %.3e %.3e  ratios %.1f %.1f\n",
                stop, drifts..., drifts[1] / drifts[2], drifts[2] / drifts[3])
    end
end

println("\nOn saved primary points (steps 0.4 and 0.2 read saved primary samples)")
for index in (100, 232)
    (; profiles) = primary_profiles_at(setup, index)
    secondary_frequency = real(secondary_branch.rows[index].solution.frequency)
    radiative_frequency =
        2secondary_frequency - real(primary.continuation[index].solution.frequency)
    for (branch, stop) in ((profiles.from_outer, 61.0), (profiles.from_horizon, -59.0))
        lookup = position -> primary_values_at(branch, position, primary.modes)
        drifts = [current_drift(lookup, radiative_frequency, secondary_frequency,
                                5.0, stop, step) for step in (0.4, 0.2)]
        @printf("  A1=%-7.4g omega3=%.6f  5 -> %5.1f  drift %.3e %.3e  ratio %.1f\n",
                primary.continuation[index].amplitude, radiative_frequency, stop,
                drifts..., drifts[1] / drifts[2])
    end
end

println("\nCheck 2: source projections |S_l| / max over odd l, at saved secondary samples")
all_modes = [(l, 3) for l in 3:maximum_l]
for index in (1, 100, 232)
    (; profiles) = primary_profiles_at(setup, index)
    solution = secondary_branch.rows[index].solution
    secondary = secondary_matching_residual(
        solution.frequency, solution.horizon_amplitudes, solution.outer_amplitudes,
        profiles, primary.modes, spin, alpha, secondary_modes, couplings,
        polar, azimuthal, grid; save_profile=true,
    )
    for position in (-20.0, 5.0, 20.0, 50.0)
        primary_branch, secondary_path_branch = position <= primary.options.match ?
            (profiles.from_horizon, secondary.from_horizon) :
            (profiles.from_outer, secondary.from_outer)
        sample = argmin(abs.(secondary_path_branch.positions .- position))
        abs(secondary_path_branch.positions[sample] - position) < 1e-9 || error("misaligned")
        secondary_values = Dict(mode => secondary_path_branch.states[i, sample]
                                for (i, mode) in enumerate(secondary_modes))
        primary_values = primary_values_at(primary_branch, position, primary.modes)
        radius = tortoise_geometry(position, spin).radius
        source = radiative_source(primary_values, secondary_values, radius, spin, alpha,
                                  all_modes, polar, azimuthal, grid)
        scale = maximum(abs(source[(l, 3)]) for l in 3:2:maximum_l)
        @printf("A1=%-7.4g r*=%5.1f  max|S|=%.2e  odd l 3,5,7,9: %s  even l 4,6,8: %s\n",
                primary.continuation[index].amplitude, position, scale,
                join([@sprintf("%.1e", abs(source[(l, 3)]) / scale) for l in 3:2:maximum_l], " "),
                join([@sprintf("%.0e", abs(source[(l, 3)]) / scale) for l in 4:2:maximum_l], " "))
    end
end
