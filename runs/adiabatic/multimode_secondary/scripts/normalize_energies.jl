include(joinpath(
    @__DIR__, "../../../../src/adiabatic/drivers/two_cloud_evolution.jl"
))

"""
    read_energy_trajectory(path)

Read the saved growth-scaled times and the two cloud energies from a CSV.
"""
function read_energy_trajectory(path)
    columns = read_csv_columns(path)
    return (
        time_growth_units = columns["time"],
        primary_energy = columns["primary_energy"],
        secondary_energy = columns["secondary_energy"],
    )
end

"""
    scalar_accounting(alpha, label)

Reconstruct and integrate all retained scalar boundary fluxes for one saved
trajectory. Return energies, angular momenta, and time in code units.

Horizon fluxes are positive into the BH; net extraction has the opposite
sign. Outward loss includes the finite-boundary currents and m = 3 infinity
flux used in the original evolution. No physical normalization is applied.
"""
function scalar_accounting(alpha, label)
    trajectory_path = joinpath(
        @__DIR__, "../results/two_cloud_evolution_$(label).csv"
    )
    trajectory = read_energy_trajectory(trajectory_path)

    primary_path, _, secondary_path, dissipative_path, radiative_path =
        TWO_CLOUD_TABLES[alpha]
    primary_branch = read_csv_columns(primary_path)
    secondary_table = read_csv_columns(secondary_path)
    dissipative_table = read_csv_columns(dissipative_path)
    radiative_table = read_csv_columns(radiative_path)

    primary_energy = trajectory.primary_energy
    secondary_energy = trajectory.secondary_energy

    # Reconstruct each channel on the same trajectory.
    primary = primary_boundary_fluxes(primary_branch, primary_energy)
    secondary = secondary_boundary_fluxes(
        primary_branch, secondary_table, primary_energy, secondary_energy
    )
    dissipative = dissipative_boundary_fluxes(
        primary_branch, dissipative_table, primary_energy, secondary_energy
    )
    radiative = radiative_boundary_fluxes(
        primary_branch, radiative_table, primary_energy, secondary_energy
    )

    coefficients = two_cloud_coefficients(read_two_cloud_table(alpha))
    torques = scalar_boundary_torques(
        coefficients, primary_energy, primary, secondary, radiative
    )

    horizon_energy_flux = primary.horizon .+ secondary.horizon .+
        dissipative.horizon .+ radiative.horizon
    outward_energy_flux = primary.outer .+ secondary.outer .+
        dissipative.outer .+ radiative.infinity

    # Saved time is tau = gamma₁ * t/M; all reconstructed fluxes use t/M.
    linear_growth_rate = first(primary_branch["omega_imag"])
    time_over_mass = trajectory.time_growth_units ./ linear_growth_rate

    extracted_energy = -cumulative_trapezoid(time_over_mass, horizon_energy_flux)
    outward_energy = cumulative_trapezoid(time_over_mass, outward_energy_flux)
    extracted_angular_momentum = -cumulative_trapezoid(time_over_mass, torques.horizon)
    outward_angular_momentum = cumulative_trapezoid(time_over_mass, torques.outward)

    # Initial cloud energy + extraction = remaining cloud energy + outward loss.
    cloud_energy = primary_energy .+ secondary_energy
    energy_residual = cloud_energy .- first(cloud_energy) .-
        extracted_energy .+ outward_energy
    relative_energy_error = maximum(abs.(energy_residual)) / maximum(cloud_energy)

    return (;
        time_over_mass, primary_energy, secondary_energy,
        horizon_energy_flux, horizon_angular_momentum_flux=torques.horizon,
        extracted_energy, outward_energy,
        extracted_angular_momentum, outward_angular_momentum,
        relative_energy_error,
    )
end

"""
    normalize_accounting(accounting, reference_units, decay_constant, initial_spin)

Convert one trajectory's scalar accounting to physical time, horizon power
and torque, and fractions of the initial BH mass and angular momentum.
Supply reference units from `black_hole_reference_units` and 0 < initial_spin <= 1.

The reservoir flag tests whether cumulative extraction exceeds the initial
Kerr rotational energy or angular momentum anywhere on the trajectory.
A false flag does not establish fixed-background validity.
"""
function normalize_accounting(accounting, reference_units, decay_constant, initial_spin)
    0 < initial_spin <= 1 ||
        throw(ArgumentError("initial spin must satisfy 0 < chi <= 1"))

    # E_phys = epsilon_f * E₀ * E_code; J_phys = epsilon_f * J₀,unit * J_code.
    energy_scale = cloud_mass_fraction(1.0, decay_constant)
    angular_momentum_scale = energy_scale / initial_spin
    power_scale = energy_scale * reference_units.energy / reference_units.time
    torque_scale = energy_scale * reference_units.angular_momentum / reference_units.time

    primary_mass_fraction = energy_scale .* accounting.primary_energy
    secondary_mass_fraction = energy_scale .* accounting.secondary_energy
    extracted_mass_fraction = energy_scale .* accounting.extracted_energy
    extracted_angular_momentum_fraction =
        angular_momentum_scale .* accounting.extracted_angular_momentum

    # The reversible Kerr rotational-energy bound, in units of M₀*c².
    rotational_energy_fraction = 1 - sqrt((1 + sqrt(1 - initial_spin^2)) / 2)
    exceeds_energy_reservoir = maximum(extracted_mass_fraction) > rotational_energy_fraction
    exceeds_angular_momentum_reservoir = maximum(extracted_angular_momentum_fraction) > 1

    time_days = uconvert.(u"d", accounting.time_over_mass .* reference_units.time)
    horizon_power = uconvert.(u"W", accounting.horizon_energy_flux .* power_scale)
    horizon_torque = uconvert.(u"N*m", accounting.horizon_angular_momentum_flux .* torque_scale)

    return (;
        time_days, primary_mass_fraction, secondary_mass_fraction,
        horizon_power, horizon_torque, extracted_mass_fraction,
        outward_mass_fraction=energy_scale .* accounting.outward_energy,
        extracted_angular_momentum_fraction,
        outward_angular_momentum_fraction=
            angular_momentum_scale .* accounting.outward_angular_momentum,
        exceeds_reservoir=exceeds_energy_reservoir || exceeds_angular_momentum_reservoir,
    )
end

"""
    print_normalization_summary(alpha, accounting, decay_constants, reports)

Print the duration, energy-accounting error, and one extraction summary per
decay constant. `reports` contains normalized records in the same order as
`decay_constants`. Ratios refer to the initial BH, not its evolved state.
"""
function print_normalization_summary(alpha, accounting, decay_constants, reports)
    time_days = first(reports).time_days
    println("\nalpha = ", alpha, ": elapsed time = ", last(time_days) - first(time_days))
    println("Maximum energy residual / peak cloud energy = ", accounting.relative_energy_error)
    println("fa [GeV]   peak cloud/M0   final extracted E/M0   final extracted J/J_BH0   exceeds reservoir")

    for (decay_constant, report) in zip(decay_constants, reports)
        peak_cloud_fraction = maximum(
            report.primary_mass_fraction .+ report.secondary_mass_fraction
        )
        final_energy_fraction = last(report.extracted_mass_fraction)
        final_angular_momentum_fraction = last(report.extracted_angular_momentum_fraction)

        @printf("%.0e     %.6g          %.6g               %.6g                 %s\n",
            ustrip(u"GeV", decay_constant), peak_cloud_fraction,
            final_energy_fraction, final_angular_momentum_fraction,
            string(report.exceeds_reservoir))
    end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    # Campaign choices: fixed-background diagnostics, with the Regge figure's axion mass.
    axion_energy = 3e-13u"eV"
    initial_spin = 0.99
    decay_constants = [1e15, 1e16, 1e17, 1e18] .* u"GeV"
    normalization_reports = Dict()

    for (alpha, label) in ((0.42, "alpha042"), (0.3, "alpha030"))
        accounting = scalar_accounting(alpha, label)
        reference_units = black_hole_reference_units(alpha, axion_energy)
        reports = [
            normalize_accounting(accounting, reference_units, decay_constant, initial_spin)
            for decay_constant in decay_constants
        ]

        print_normalization_summary(alpha, accounting, decay_constants, reports)
        for (decay_constant, report) in zip(decay_constants, reports)
            normalization_reports[(alpha, decay_constant)] = report
        end
    end
end
