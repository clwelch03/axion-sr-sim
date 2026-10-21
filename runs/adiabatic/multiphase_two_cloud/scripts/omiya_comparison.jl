#= Build a dilute Omiya-approximation branch and compare at the joint run's times.

The initial outer coefficients, Kerr parameters, domain, angular support and
boundary-series orders match the joint pilot. The existing primary/secondary
solver needs primary samples at every secondary RK4 stage, so it uses radial
steps 0.1/0.2. This is a finite-domain comparison of the two approximations,
not a reproduction of the paper's full-domain numerical settings.
=#
using Serialization, SHA, Printf, LinearAlgebra
include(joinpath(@__DIR__, "../../../../src/adiabatic/primary_secondary/load.jl"))
include(joinpath(@__DIR__, "../../../../src/adiabatic/multiphase/load.jl"))

function save_comparison_checkpoint(path, value)
    serialize(path * ".tmp", value)
    mv(path * ".tmp", path; force=true)
end

function write_comparison_rows(path, rows)
    open(path, "w") do io
        println(io, join(string.(keys(first(rows))), ','))
        for row in rows
            println(io, join(values(row), ','))
        end
    end
end

function omiya_branch_point(amplitude, problem, linear_primary, linear_secondary)
    (; spin, alpha, grid, polar_values, couplings) = problem
    azimuthal = build_azimuthal_values(3, grid.azimuth_nodes)
    primary_modes, secondary_modes = [(1, 1)], [(2, 2)]
    (; inner, outer, match, horizon_order, outer_order) = problem.integration
    options = (; inner, outer, match, step=0.1, horizon_order, outer_order)
    primary = solve_nonlinear_mode(linear_primary.frequency, amplitude, spin, alpha,
        primary_modes, couplings, polar_values, azimuthal, grid;
        horizon_guess=amplitude .* linear_primary.horizon_amplitudes,
        outer_tail_guess=ComplexF64[], options..., tolerance=1e-9,
        jacobian_backend=:forwarddiff)
    profiles = nonlinear_matching_residual(primary.frequency,
        primary.horizon_amplitudes, primary.outer_amplitudes, spin, alpha,
        primary_modes, couplings, polar_values, azimuthal, grid;
        options..., save_profile=true)
    primary_energy = sum(integrate_energy_branch(branch, primary.frequency,
        spin, alpha, primary_modes, polar_values, azimuthal, grid)
        for branch in (profiles.from_horizon, profiles.from_outer))
    secondary = solve_secondary_mode(linear_secondary.frequency, profiles,
        primary_modes, spin, alpha, secondary_modes, couplings, polar_values,
        azimuthal, grid; horizon_guess=linear_secondary.horizon_amplitudes,
        outer_tail_guess=ComplexF64[], tolerance=1e-8, jacobian_backend=:forwarddiff,
        reading=:averaged, horizon_order, outer_order)
    secondary_observables = secondary_energy_and_fluxes(secondary, profiles,
        primary_modes, spin, alpha, secondary_modes, couplings, polar_values,
        azimuthal, grid)
    omega1, omega2 = real(primary.frequency), real(secondary.frequency)
    dissipative = solve_dissipative_mode(2omega1 - omega2, secondary, profiles,
        primary_modes, spin, alpha, secondary_modes, [(0, 0)], couplings,
        polar_values, azimuthal, grid; horizon_order, outer_order)
    # Match the joint pilot's outer-series prescription at r*=80. The old
    # Fig. 6 reproduction instead propagates the outgoing solution from far out.
    radiative = solve_radiative_mode(2omega2 - omega1, secondary, profiles,
        primary_modes, spin, alpha, secondary_modes, [(3, 3)], couplings,
        polar_values, azimuthal, grid; horizon_order, outer_order,
        far_position=outer, far_step=0.2, far_order=outer_order)
    current(branch) = outward_secondary_flux(branch.states[:, 1],
        first(branch.positions), primary.frequency, spin, 1)
    primary_flux = -current(profiles.from_horizon) + current(profiles.from_outer)
    secondary_energy = secondary_observables.energy
    row = (; amplitude, primary_energy, secondary_energy,
        primary_omega=omega1, secondary_omega=omega2,
        primary_flux_per_energy=primary_flux / primary_energy,
        secondary_flux_per_energy=(secondary_observables.horizon_flux +
            secondary_observables.outer_flux) / secondary_energy,
        m0_flux_coefficient=(dissipative.horizon_flux + dissipative.outer_flux) /
            (secondary_energy * primary_energy^2),
        m3_flux_coefficient=(radiative.horizon_flux + radiative.infinity_flux) /
            (secondary_energy^2 * primary_energy),
        primary_imag=imag(primary.frequency), secondary_imag=imag(secondary.frequency),
        primary_residual=primary.scaled_residual_norm,
        secondary_residual=secondary.residual_norm,
        m0_relative_residual=dissipative.residual_norm / dissipative.source_residual_norm,
        m3_relative_residual=radiative.residual_norm / radiative.source_residual_norm)
    all(isfinite, values(row)) || error("nonfinite branch point")
    max(row.m0_relative_residual, row.m3_relative_residual) < 1e-7 ||
        error("driven matching failed")
    @printf("Omiya branch A1=%.9g E1=%.9g E2/A2^2=%.9g residuals=%.3g/%.3g\n",
        amplitude, primary_energy, secondary_energy,
        row.primary_residual, row.secondary_residual)
    flush(stdout)
    return (; row, primary, secondary)
end

function omiya_trajectory(rows, initial_amplitudes, times; subdivisions=8)
    table = NamedTuple{keys(first(rows))}(Tuple(
        [getproperty(row, name) for row in rows] for name in keys(first(rows))))
    coefficients = two_cloud_coefficients(table)
    initial_index = findfirst(==(initial_amplitudes[1]), table.amplitude)
    isnothing(initial_index) && error("initial amplitude must be a solved branch point")
    initial_energy = table.primary_energy[initial_index]
    # Smooth normalized maps avoid interpolating sqrt(E1) directly.
    primary_map = cubic_spline(table.primary_energy,
        table.amplitude.^2 ./ table.primary_energy)
    secondary_map = cubic_spline(table.primary_energy, table.secondary_energy)
    state = log.([initial_energy,
        initial_amplitudes[2]^2 * table.secondary_energy[initial_index]])
    function log_rates(state)
        energies = exp.(state)
        rates = two_cloud_rates(coefficients, energies...)
        return [rates.primary_rate, rates.secondary_rate] ./ energies
    end
    output = NamedTuple[]
    for index in eachindex(times)
        if index > 1
            spacing = (times[index] - times[index-1]) / subdivisions
            for _ in 1:subdivisions
                k1 = log_rates(state)
                k2 = log_rates(state + spacing/2 * k1)
                k3 = log_rates(state + spacing/2 * k2)
                k4 = log_rates(state + spacing * k3)
                state += spacing/6 * (k1 + 2k2 + 2k3 + k4)
            end
        end
        primary_energy, secondary_energy = exp.(state)
        A1 = sqrt(primary_energy * evaluate_spline(primary_map, primary_energy).value)
        A2 = sqrt(secondary_energy / evaluate_spline(secondary_map, primary_energy).value)
        push!(output, (; time=times[index], A1, A2, primary_energy, secondary_energy))
    end
    maximum(abs.([output[1].A1, output[1].A2] ./ initial_amplitudes .- 1)) < 1e-12 ||
        error("initial amplitude normalization mismatch")
    return output
end

function main(joint_prefix, output_prefix)
    mkdir(output_prefix * ".lock")
    try
        problem, initial = deserialize(joinpath(@__DIR__, "../../../../inputs/initial.jls"))
        problem.potential == :cosine && problem.cloud_m == [1, 2] || error("unsupported model")
        [(component.q, component.l_values) for component in problem.components] ==
            [((1, 0), [1]), ((0, 1), [2]), ((2, -1), [0]), ((-1, 2), [3])] ||
            error("this comparison requires the pilot's four components")
        problem.integration.horizon_order == 3 && problem.integration.outer_order == 7 ||
            error("observable wrapper assumes series orders 3/7")
        trajectory_path = joint_prefix * "-steps64.csv"
        times = read_csv_columns(trajectory_path)["time"]
        all(diff(times) .> 0) && first(times) == 0 || error("invalid joint times")
        initial_amplitudes = initial.solution.amplitudes
        amplitudes = sort(unique(vcat(initial_amplitudes[1], initial_amplitudes[1] .*
            sqrt.(collect(range(0.9^2, 2.2^2; length=17))))))
        configuration = (; manifest=read(joinpath(@__DIR__, "../../../../provenance.json"), String),
            trajectory_sha256=bytes2hex(sha256(read(trajectory_path))), amplitudes,
            spin=problem.spin, alpha=problem.alpha, primary_step=0.1, secondary_step=0.2)
        checkpoint_path = output_prefix * "-branch.jls"
        points = if isfile(checkpoint_path)
            saved = deserialize(checkpoint_path)
            saved.configuration == configuration || error("comparison configuration changed")
            saved.points
        else
            NamedTuple[]
        end
        (; inner, outer, match) = problem.integration
        reference_frequencies = unpack_complex(initial.solution.state)[1:2]
        linear_primary = solve_linear_mode(reference_frequencies[1], problem.spin,
            problem.alpha, [1], 1, problem.couplings;
            inner, outer, match, step=0.1, tolerance=1e-9)
        linear_secondary = solve_linear_mode(reference_frequencies[2], problem.spin,
            problem.alpha, [2], 2, problem.couplings;
            inner, outer, match, step=0.2, tolerance=1e-8)
        for amplitude in amplitudes[length(points)+1:end]
            push!(points, omiya_branch_point(amplitude, problem, linear_primary, linear_secondary))
            save_comparison_checkpoint(checkpoint_path, (; configuration, points))
            write_comparison_rows(output_prefix * "-branch.csv", [point.row for point in points])
        end
        rows = [point.row for point in points]
        trajectory = omiya_trajectory(rows, initial_amplitudes, times)
        half_steps = omiya_trajectory(rows, initial_amplitudes, times; subdivisions=4)
        initial_index = findfirst(row -> row.amplitude == initial_amplitudes[1], rows)
        coarse_indices = sort(unique(vcat(collect(1:2:length(rows)), initial_index, length(rows))))
        coarse_table = omiya_trajectory(rows[coarse_indices], initial_amplitudes, times)
        discrepancy(other) = maximum(abs(getproperty(a, name) / getproperty(b, name) - 1)
            for (a, b) in zip(trajectory, other) for name in (:A1, :A2))
        time_error, table_error = discrepancy(half_steps), discrepancy(coarse_table)
        println("Time refinement amplitude discrepancy: ", time_error)
        println("Table thinning amplitude discrepancy: ", table_error)
        time_error < 1e-9 || error("time refinement failed")
        table_error < 1e-7 || error("table thinning failed")
        write_comparison_rows(output_prefix * "-trajectory.csv", trajectory)
        save_comparison_checkpoint(output_prefix * "-summary.jls",
            (; configuration, initial_amplitudes, final=last(trajectory), time_error, table_error))
        println("COMPLETED Omiya comparison: ", last(trajectory))
    finally
        rm(output_prefix * ".lock"; recursive=true)
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS...)
