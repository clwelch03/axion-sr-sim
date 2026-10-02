# julia --project=. src/adiabatic/drivers/continue_secondary.jl primary_branch.jls output_prefix
# Solve the secondary's Eq. (61) (Omiya et al. 2211.01949) at every saved primary
# point, in amplitude order. Set AXION_MAX_POINTS to stop after that many points.
using Serialization, Printf, LinearAlgebra
include(joinpath(@__DIR__, "../load.jl"))

"""
    load_secondary_setup(primary_path; secondary_modes=[(2, 2), (4, 2), (6, 2)])

Read a primary branch checkpoint and build angular tables covering both the
primary's and the secondary's modes.
"""
function load_secondary_setup(primary_path; secondary_modes=[(2, 2), (4, 2), (6, 2)])
    primary = deserialize(primary_path)
    maximum_l = max(maximum(first, primary.modes), maximum(first, secondary_modes))
    grid = angular_grid(primary.n_polar, primary.n_azimuth)
    polar = build_polar_values(maximum_l, grid.cosine_theta_nodes)
    azimuthal = build_azimuthal_values(maximum_l, grid.azimuth_nodes)
    couplings = build_cosine_squared_couplings(maximum_l, polar, grid)
    return (; primary_path, primary, secondary_modes, grid, polar, azimuthal, couplings)
end

"""
    primary_profiles_at(setup, index)

Reintegrate saved primary point `index` with its own radial settings. Return the
profiles and the change in matching residual relative to the saved solve,
divided by A1 (should be round-off).
"""
function primary_profiles_at(setup, index)
    (; primary, grid, polar, azimuthal, couplings) = setup
    entry = primary.continuation[index]
    solution = entry.solution
    profiles = nonlinear_matching_residual(
        solution.frequency, solution.horizon_amplitudes, solution.outer_amplitudes,
        primary.spin, primary.alpha, primary.modes, couplings, polar, azimuthal, grid;
        potential=Symbol(primary.potential), save_profile=true, primary.options...,
    )
    residual_change = maximum(abs, profiles.residual - solution.residual) / entry.amplitude
    return (; profiles, residual_change)
end

# Newton unknowns: frequency, horizon amplitudes, free outer amplitudes.
secondary_unknowns(solution) =
    vcat(solution.frequency, solution.horizon_amplitudes, solution.outer_amplitudes[2:end])

function solve_secondary_point(setup, profiles, guess; tolerance, jacobian_backend, reading)
    (; primary, secondary_modes, grid, polar, azimuthal, couplings) = setup
    count = length(secondary_modes)
    return solve_secondary_mode(
        guess[1], profiles, primary.modes, primary.spin, primary.alpha,
        secondary_modes, couplings, polar, azimuthal, grid;
        horizon_guess=guess[2:count+1], outer_tail_guess=guess[count+2:end],
        tolerance, jacobian_backend, reading,
    )
end

"""
    save_secondary_branch(output_prefix, setup, settings, linear, rows)

Write the full solutions (`.jls`, for resuming) and a table (`.csv`), each
replaced only after its temporary file is complete.
"""
function save_secondary_branch(output_prefix, setup, settings, linear, rows)
    data = (;
        primary_path=setup.primary_path, secondary_modes=setup.secondary_modes,
        spin=setup.primary.spin, alpha=setup.primary.alpha,
        primary_options=setup.primary.options, settings..., linear, rows,
    )
    temporary = output_prefix * ".jls.tmp"
    serialize(temporary, data)
    mv(temporary, output_prefix * ".jls"; force=true)

    temporary = output_prefix * ".csv.tmp"
    open(temporary, "w") do io
        println(io, "primary_index,amplitude,omega1_real,omega1_imag,omega2_real,",
                    "omega2_imag,residual,iterations,primary_residual_change")
        for row in rows
            @printf(io, "%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.6g,%d,%.6g\n",
                    row.primary_index, row.amplitude, real(row.omega1), imag(row.omega1),
                    real(row.solution.frequency), imag(row.solution.frequency),
                    row.solution.residual_norm, row.solution.iterations,
                    row.primary_residual_change)
        end
    end
    mv(temporary, output_prefix * ".csv"; force=true)
end

"""
    continue_secondary(primary_path, output_prefix; tolerance=1e-8,
                       ulp_tolerance_factor=nothing, jacobian_backend=:forwarddiff,
                       reading=:averaged, maximum_points=nothing)

Solve the secondary at each primary point in amplitude order. Predict each
guess by secant extrapolation in A1 from the two previous solutions, with the
linear secondary as the A1 = 0 point. Save after every point; resume an
existing output. Stop with an error at the first failed solve.

`reading` selects the background term of [`secondary_background_source`](@ref).
Use the default `:averaged`; `:fixed_time` makes omega2 depend on the
secondary's normalization phase (see the campaign README).

With `ulp_tolerance_factor`, each point after the first uses
`max(tolerance, factor * |res(ω + ulp) - res(ω)|)`, the residual's change from one
unit in the last place of Re ω₂, measured at the previous solution. No
representable ω₂ gives a smaller residual. This sensitivity grows as the
secondary contracts away from its fixed outer amplitude, so a fixed absolute
`tolerance` (used otherwise, and always for the linear A1 = 0 solve) eventually
falls below the floor.
"""
function continue_secondary(primary_path, output_prefix; tolerance=1e-8,
                            ulp_tolerance_factor=nothing, jacobian_backend=:forwarddiff,
                            reading=:averaged, maximum_points=nothing)
    setup = load_secondary_setup(primary_path)
    (; primary, secondary_modes, couplings) = setup
    settings = (; tolerance, ulp_tolerance_factor, jacobian_backend=String(jacobian_backend),
                  reading=String(reading))

    if isfile(output_prefix * ".jls")
        previous = deserialize(output_prefix * ".jls")
        previous.primary_path == primary_path && previous.tolerance == tolerance &&
            get(previous, :ulp_tolerance_factor, nothing) == ulp_tolerance_factor &&
            get(previous, :reading, "fixed_time") == String(reading) ||
            error("existing output has a different primary, tolerance, or reading")
        linear, rows = previous.linear, previous.rows
        println("resuming after ", length(rows), " points")
    else
        # A1 = 0: the linear secondary on the primary's domain, at twice its step.
        # Guess: hydrogenic n = 3 frequency alpha * (1 - alpha^2 / 18).
        (; inner, outer, match, step) = primary.options
        linear = solve_linear_mode(
            primary.alpha * (1 - primary.alpha^2 / 18) + 0im, primary.spin, primary.alpha,
            first.(secondary_modes), last(first(secondary_modes)), couplings;
            step=2step, inner, outer, match, tolerance,
        )
        rows = []
        save_secondary_branch(output_prefix, setup, settings, linear, rows)
        @printf("  0 A1=0         omega2=%.12f %+.6ei (linear)\n",
                real(linear.frequency), imag(linear.frequency))
    end

    last_index = maximum_points === nothing ? length(primary.continuation) :
                 min(length(primary.continuation), maximum_points)

    for index in (length(rows) + 1):last_index
        amplitude = primary.continuation[index].amplitude

        # Secant prediction in A1 from the two most recent solutions. The first
        # point has only the linear solve behind it and starts from that.
        history = vcat(
            [(0.0, secondary_unknowns(linear))],
            [(row.amplitude, secondary_unknowns(row.solution)) for row in rows],
        )
        if length(history) >= 2
            (older_amplitude, older), (newer_amplitude, newer) = history[end-1], history[end]
            guess = newer .+ (amplitude - newer_amplitude) /
                             (newer_amplitude - older_amplitude) .* (newer .- older)
        else
            guess = last(history)[2]
        end

        (; profiles, residual_change) = primary_profiles_at(setup, index)
        # The residual cannot fall below its change from one ulp of Re omega2;
        # measure that at the previous solution on this primary.
        point_tolerance = tolerance
        if ulp_tolerance_factor !== nothing && !isempty(rows)
            previous = last(rows).solution
            residual_at(frequency) = secondary_matching_residual(
                frequency, previous.horizon_amplitudes, previous.outer_amplitudes,
                profiles, primary.modes, primary.spin, primary.alpha, secondary_modes,
                couplings, setup.polar, setup.azimuthal, setup.grid; reading,
            )
            ulp_sensitivity = norm(
                residual_at(previous.frequency + eps(real(previous.frequency))) -
                residual_at(previous.frequency))
            point_tolerance = max(tolerance, ulp_tolerance_factor * ulp_sensitivity)
        end
        elapsed = @elapsed solution = solve_secondary_point(
            setup, profiles, guess; tolerance=point_tolerance, jacobian_backend, reading,
        )

        push!(rows, (; primary_index=index, amplitude,
                       omega1=primary.continuation[index].solution.frequency,
                       solution, primary_residual_change=residual_change,
                       tolerance=point_tolerance))
        save_secondary_branch(output_prefix, setup, settings, linear, rows)

        @printf("%3d A1=%-9.5g omega2=%.12f %+.6ei it=%d res=%.1e t=%.1fs\n",
                index, amplitude, real(solution.frequency), imag(solution.frequency),
                solution.iterations, solution.residual_norm, elapsed)
        flush(stdout)
    end
    return rows
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    length(ARGS) == 2 || error("require primary_branch.jls output_prefix")
    maximum_points = haskey(ENV, "AXION_MAX_POINTS") ?
        parse(Int, ENV["AXION_MAX_POINTS"]) : nothing
    # Mμ = 0.3 (outer r* = 300) needs AXION_SECONDARY_ULP_TOLERANCE_FACTOR=2.
    tolerance = parse(Float64, get(ENV, "AXION_SECONDARY_TOLERANCE", "1e-8"))
    ulp_tolerance_factor = haskey(ENV, "AXION_SECONDARY_ULP_TOLERANCE_FACTOR") ?
        parse(Float64, ENV["AXION_SECONDARY_ULP_TOLERANCE_FACTOR"]) : nothing
    continue_secondary(ARGS[1], ARGS[2]; tolerance, ulp_tolerance_factor, maximum_points)
end
