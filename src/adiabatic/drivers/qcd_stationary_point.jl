# Stationary inputs for the preliminary QCD spin-down campaign. This composes
# the existing Omiya solvers; the primary is nonlinear and the secondary uses
# the phase-averaged primary background. No BH time evolution occurs here.
using Serialization, SHA
include(joinpath(@__DIR__, "radiative_fluxes.jl"))

"""Reconstruct one primary solution and integrate its code energy.

Keep this helper outside the stationary solve so diagnostic evaluations cannot
capture and overwrite the caller's primary energy or profiles.
"""
function qcd_primary_profiles_and_energy(solution;
        alpha, spin, modes, couplings, polar, azimuthal, grid, options)
    profiles = nonlinear_matching_residual(
        solution.frequency, solution.horizon_amplitudes, solution.outer_amplitudes,
        spin, alpha, modes, couplings, polar, azimuthal, grid;
        options..., save_profile=true,
    )
    energy = sum(
        integrate_energy_branch(
            branch, solution.frequency, spin, alpha, modes, polar, azimuthal, grid,
        )
        for branch in (profiles.from_horizon, profiles.from_outer)
    )
    return profiles, energy
end

"""Solve the nonlinear primary, averaged secondary, and both driven fields.

Inputs alpha=M_BH*mu and spin=J_BH/M_BH^2 use the current BH mass as unit.
target_energy is the primary energy divided by G*fa^2*M_BH. Returned fluxes
are factored by their cloud-energy powers, as in two_cloud_evolution.jl.
The potential is the repository's cosine potential with QCD mass/fa mapping.
"""
function qcd_stationary_point(alpha, spin, target_energy; step=0.25, outer=2000.0,
                              difference_step=1e-6, secondary_tolerance=1e-5,
                              field_diagnostics=false)
    alpha > 0 || throw(ArgumentError("alpha must be positive"))
    0 <= spin < 1 || throw(ArgumentError("spin must be prograde and subextremal"))
    target_energy > 0 || throw(ArgumentError("target primary energy must be positive"))

    grid = angular_grid(12, 25)
    polar = build_polar_values(7, grid.cosine_theta_nodes)
    azimuthal = build_azimuthal_values(7, grid.azimuth_nodes)
    couplings = build_cosine_squared_couplings(7, polar, grid)
    modes = [(1, 1), (3, 1), (5, 1), (3, 3), (5, 3), (5, 5)]
    secondary_modes = [(2, 2), (4, 2), (6, 2)]
    dissipative_modes = [(0, 0), (2, 0), (4, 0)]
    radiative_modes = [(3, 3), (5, 3), (7, 3)]
    options = (; step, inner=-100.0, outer, match=40.0)
    profile_setup = (; alpha, spin, modes, couplings, polar, azimuthal, grid, options)
    started = time()

    linear = solve_linear_mode(
        complex(alpha * (1 - alpha^2 / 8)), spin, alpha,
        [1, 3, 5], 1, couplings; options..., tolerance=1e-7, difference_step,
    )
    println("LINEAR $(linear.frequency), residual=$(linear.residual_norm)")
    flush(stdout)

    # An arbitrarily small amplitude gives the normalization E/A^2. Then
    # correct the actual finite-amplitude solution to the requested energy.
    amplitude = 1e-8
    primary = (;
        frequency=linear.frequency,
        horizon_amplitudes=amplitude .* vcat(linear.horizon_amplitudes, zeros(ComplexF64, 3)),
        outer_amplitudes=amplitude .* vcat(linear.outer_amplitudes, zeros(ComplexF64, 3)),
    )
    profiles, energy = qcd_primary_profiles_and_energy(primary; profile_setup...)
    for iteration in 1:4
        amplitude_ratio = sqrt(target_energy / energy)
        amplitude *= amplitude_ratio
        primary = solve_nonlinear_mode(
            primary.frequency, amplitude, spin, alpha, modes, couplings, polar, azimuthal, grid;
            horizon_guess=amplitude_ratio .* primary.horizon_amplitudes,
            outer_tail_guess=amplitude_ratio .* primary.outer_amplitudes[2:end],
            options..., tolerance=1e-7, jacobian_backend=:forwarddiff,
        )
        profiles, energy = qcd_primary_profiles_and_energy(primary; profile_setup...)
        println("PRIMARY E=$energy, A=$amplitude, residual=$(primary.scaled_residual_norm)")
        flush(stdout)
        abs(energy / target_energy - 1) < 1e-6 && break
    end
    abs(energy / target_energy - 1) < 1e-6 || error("primary energy inversion failed")

    # The 322 solution is linear in its own normalization, but its radial
    # operator includes the phase-averaged nonlinear 211 background.
    secondary_linear = solve_linear_mode(
        complex(alpha * (1 - alpha^2 / 18)), spin, alpha,
        [2, 4, 6], 2, couplings; options..., step=2step,
        tolerance=secondary_tolerance, difference_step,
    )
    secondary = solve_secondary_mode(
        secondary_linear.frequency, profiles, modes, spin, alpha,
        secondary_modes, couplings, polar, azimuthal, grid;
        horizon_guess=secondary_linear.horizon_amplitudes,
        outer_tail_guess=secondary_linear.outer_amplitudes[2:end],
        tolerance=secondary_tolerance, jacobian_backend=:forwarddiff, reading=:averaged,
    )
    secondary_observables = secondary_energy_and_fluxes(
        secondary, profiles, modes, spin, alpha,
        secondary_modes, couplings, polar, azimuthal, grid,
    )
    println("SECONDARY $(secondary.frequency), residual=$(secondary.residual_norm)")
    flush(stdout)

    # The cloud angular-momentum balance needs d(omega2)/d(E1) at FIXED
    # alpha and spin. A derivative along the evolving BH trajectory would be
    # a different quantity. Make another actual primary/secondary solve here.
    # A finite energy increment resolves the tiny QCD frequency response.
    energy_increment = 0.01
    shifted_amplitude = amplitude * sqrt((energy + energy_increment) / energy)
    shifted_primary = solve_nonlinear_mode(
        primary.frequency, shifted_amplitude, spin, alpha,
        modes, couplings, polar, azimuthal, grid;
        horizon_guess=(shifted_amplitude / amplitude) .* primary.horizon_amplitudes,
        outer_tail_guess=(shifted_amplitude / amplitude) .* primary.outer_amplitudes[2:end],
        options..., tolerance=1e-7, jacobian_backend=:forwarddiff,
    )
    shifted_profiles, shifted_energy = qcd_primary_profiles_and_energy(
        shifted_primary; profile_setup...,
    )
    shifted_secondary = solve_secondary_mode(
        secondary.frequency, shifted_profiles, modes, spin, alpha,
        secondary_modes, couplings, polar, azimuthal, grid;
        horizon_guess=secondary.horizon_amplitudes,
        outer_tail_guess=secondary.outer_amplitudes[2:end],
        tolerance=secondary_tolerance, jacobian_backend=:forwarddiff, reading=:averaged,
    )
    secondary_slope = real(shifted_secondary.frequency - secondary.frequency) /
        (shifted_energy - energy)
    omega1, omega2 = real(primary.frequency), real(secondary.frequency)

    # These inhomogeneous fields mediate energy exchange and scalar emission.
    # The m=0 channel carries energy into the hole without angular momentum.
    dissipative = solve_dissipative_mode(
        2omega1 - omega2, secondary, profiles, modes, spin, alpha,
        secondary_modes, dissipative_modes, couplings, polar, azimuthal, grid,
    )
    radiative = solve_radiative_mode(
        2omega2 - omega1, secondary, profiles, modes, spin, alpha,
        secondary_modes, radiative_modes, couplings, polar, azimuthal, grid,
    )
    function primary_current(branch)
        state = branch.states[:, 1]
        geometry = tortoise_geometry(first(branch.positions), spin)
        return 4π * geometry.radial_factor * sum(
            real(im * m * primary.frequency * state[index] * conj(state[length(modes) + index]))
            for (index, (_, m)) in enumerate(modes)
        )
    end

    peak_field = maximum(
        maximum(abs, reconstruct_field(
            Dict(mode => branch.states[index, sample] for (index, mode) in enumerate(modes)),
            polar, azimuthal, grid,
        ))
        for branch in (profiles.from_horizon, profiles.from_outer)
        for sample in 1:10:length(branch.positions)
    )
    # Source-work identities test the driven solves independently of their
    # matching residuals. Positive boundary flux is energy lost by the clouds.
    dissipative_work = dissipative_source_work(
        dissipative, profiles, modes, 2omega1 - omega2, spin, alpha,
        secondary_modes, dissipative_modes, polar, azimuthal, grid,
    )
    radiative_work = radiative_source_work(
        radiative, profiles, modes, 2omega2 - omega1, spin, alpha,
        secondary_modes, radiative_modes, polar, azimuthal, grid,
    )
    row = (; alpha, spin, energy, amplitude, omega1, omega2, secondary_slope,
        primary_imag=imag(primary.frequency), secondary_imag=imag(secondary.frequency),
        primary_horizon=-primary_current(profiles.from_horizon) / energy,
        primary_outward=primary_current(profiles.from_outer) / energy,
        secondary_horizon=secondary_observables.horizon_flux / secondary_observables.energy,
        secondary_outward=secondary_observables.outer_flux / secondary_observables.energy,
        m0_horizon=dissipative.horizon_flux / (secondary_observables.energy * energy^2),
        m0_outward=dissipative.outer_flux / (secondary_observables.energy * energy^2),
        m3_horizon=radiative.horizon_flux / (secondary_observables.energy^2 * energy),
        m3_outward=radiative.infinity_flux / (secondary_observables.energy^2 * energy),
        primary_residual=primary.scaled_residual_norm, secondary_residual=secondary.residual_norm,
        m0_relative_residual=dissipative.residual_norm / dissipative.source_residual_norm,
        m3_relative_residual=radiative.residual_norm / radiative.source_residual_norm,
        m0_power_ratio=(dissipative.horizon_flux + dissipative.outer_flux) / dissipative_work,
        m3_power_ratio=(radiative.horizon_flux + radiative.infinity_flux) / radiative_work,
        peak_field, step, outer, elapsed_seconds=time() - started,
    )
    all(isfinite, values(row)) || error("nonfinite stationary observable")
    abs(row.energy / target_energy - 1) < 1e-6 || error("energy changed during diagnostics")
    max(row.m0_relative_residual, row.m3_relative_residual) < 1e-8 ||
        error("driven matching failure")
    max(abs(row.m0_power_ratio - 1), abs(row.m3_power_ratio - 1)) < 1e-3 ||
        error("driven power identity failure")
    if field_diagnostics
        # Maximize the polar profile on a denser grid including the equator.
        # Summing absolute values of distinct m blocks bounds their phase
        # interference at each sampled radius and polar angle. Summing the
        # separate cloud maxima later bounds the two independent phases too.
        diagnostic_polar = build_polar_values(7, collect(range(-1.0, 1.0; length=257)))
        function phase_envelope(branches, active_modes)
            largest = 0.0
            for branch in branches, sample in 1:10:length(branch.positions)
                envelope = zeros(257)
                for m in unique(last.(active_modes))
                    profile = sum(branch.states[k, sample] .* diagnostic_polar[mode]
                                  for (k, mode) in enumerate(active_modes) if last(mode) == m)
                    envelope .+= 2abs.(profile)
                end
                largest = max(largest, maximum(envelope))
            end
            return largest
        end
        secondary_profiles = secondary_matching_residual(
            secondary.frequency, secondary.horizon_amplitudes, secondary.outer_amplitudes,
            profiles, modes, spin, alpha, secondary_modes, couplings, polar, azimuthal, grid;
            save_profile=true, reading=:averaged,
        )
        row = (; row...,
            primary_peak_envelope=phase_envelope(
                (profiles.from_horizon, profiles.from_outer), modes),
            secondary_peak_per_sqrt_energy=phase_envelope(
                (secondary_profiles.from_horizon, secondary_profiles.from_outer), secondary_modes
            ) / sqrt(secondary_observables.energy),
        )
    end
    return (; row, primary, secondary)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    # Header: label,alpha,spin,energy,step,outer. One writer per output folder.
    length(ARGS) == 2 || error("usage: qcd_stationary_point.jl TARGETS.csv OUTPUT_DIRECTORY")
    targets_path, directory = ARGS
    mkpath(directory)
    source_hash = bytes2hex(sha256(read(@__FILE__)))
    Threads.@threads for line in readlines(targets_path)[2:end]
        label, numbers... = split(line, ',')
        alpha, spin, energy, step, outer = parse.(Float64, numbers)
        path = joinpath(directory, label * ".jls")
        if isfile(path)
            saved = deserialize(path)
            saved.source_hash == source_hash || error("source changed; use a new output folder")
            saved.row.alpha == alpha && saved.row.spin == spin &&
                saved.row.step == step && saved.row.outer == outer &&
                abs(saved.row.energy / energy - 1) < 1e-6 ||
                error("target changed for $label; use a new output folder")
            continue
        end
        println("START $label $numbers")
        flush(stdout)
        result = qcd_stationary_point(alpha, spin, energy; step, outer)

        # Each label has one writer. Publish the binary checkpoint atomically
        # and retain the source hash so incompatible resumes are rejected.
        serialize(path * ".tmp", (; source_hash, result...))
        mv(path * ".tmp", path)
        open(joinpath(directory, label * ".csv"), "w") do io
            println(io, join(keys(result.row), ','))
            println(io, join(values(result.row), ','))
        end
        println("ACCEPTED $label $(result.row)")
        flush(stdout)
    end
end
