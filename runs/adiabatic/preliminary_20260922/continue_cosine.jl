# Preliminary run recipe. Run from repository root with julia --project=.
using Serialization, LinearAlgebra
root = normpath(joinpath(@__DIR__, "../../.."))
for filename in ("angular", "nonlinear", "tortoise", "radial", "boundary_series", "boundaries", "integration", "matching", "energy")
    include(joinpath(root, "src/adiabatic", filename * ".jl"))
end
const output = joinpath(@__DIR__, "cosine_branch")
const seed_path = isfile(output * ".jls") ? output * ".jls" :
    joinpath(root, "runs/adiabatic/preliminary_cosine_20260922.jls")
const saved = deserialize(seed_path)
const grid = angular_grid(saved.n_polar, saved.n_azimuth)
const polar = build_polar_values(5, grid.cosine_theta_nodes)
const azimuthal = build_azimuthal_values(5, grid.azimuth_nodes)
const couplings = build_cosine_squared_couplings(5, polar, grid)
const options = saved.options
const modes = saved.modes
const spin, alpha = saved.spin, saved.alpha

# Finite-radius outward current, with positive-m and conjugate normalization.
function outward_flux(state, position, frequency)
    count = length(modes)
    geometry = tortoise_geometry(position, spin)
    return 4π * geometry.radial_factor * sum(
        real(im * m * frequency * state[index] * conj(state[count+index]))
        for (index, (l, m)) in enumerate(modes))
end

function diagnostics(amplitude, solution)
    profiles = nonlinear_matching_residual(solution.frequency,
        solution.horizon_amplitudes, solution.outer_amplitudes,
        spin, alpha, modes, couplings, polar, azimuthal, grid;
        save_profile=true, options...)
    energy = sum(integrate_energy_branch(branch, solution.frequency,
        spin, alpha, modes, polar, azimuthal, grid)
        for branch in (profiles.from_horizon, profiles.from_outer))
    horizon_flux = -outward_flux(profiles.from_horizon.states[:,1], options.inner, solution.frequency)
    outer_flux = outward_flux(profiles.from_outer.states[:,1], options.outer, solution.frequency)
    peak_field = 0.0
    for branch in (profiles.from_horizon, profiles.from_outer), index in 1:10:length(branch.positions)
        values = Dict(mode=>branch.states[k,index] for (k,mode) in enumerate(modes))
        peak_field = max(peak_field, maximum(abs, reconstruct_field(values, polar, azimuthal, grid)))
    end
    return (; amplitude, energy, normalized_energy=energy/amplitude^2,
        omega_real=real(solution.frequency), omega_imag=imag(solution.frequency),
        horizon_flux, outer_flux, total_flux=horizon_flux+outer_flux,
        balance_ratio=-(horizon_flux+outer_flux)/(2imag(solution.frequency)*energy),
        peak_field, residual=solution.scaled_residual_norm, iterations=solution.iterations)
end

function checkpoint(continuation, rows)
    serialize(output * ".jls.tmp", (; spin, alpha, modes, options, n_polar=saved.n_polar,
        n_azimuth=saved.n_azimuth, potential="cosine", tolerance=1e-9, continuation, rows))
    mv(output * ".jls.tmp", output * ".jls"; force=true)
    open(output * ".csv.tmp", "w") do io
        println(io, join(string.(keys(first(rows))), ','))
        for row in rows
            println(io, join(values(row), ','))
        end
    end
    mv(output * ".csv.tmp", output * ".csv"; force=true)
end

function main()
    started = time()
    runtime_seconds = parse(Float64, get(ENV, "AXION_RUN_SECONDS", "1000"))
    isfinite(runtime_seconds) && runtime_seconds > 0 ||
        throw(ArgumentError("AXION_RUN_SECONDS must be positive and finite"))
    continuation, rows = copy(saved.continuation), copy(saved.rows)
    amplitude, previous = last(continuation).amplitude, last(continuation).solution
    factor = amplitude < 0.03 ? 2.0 : 1.5
    checkpoint(continuation, rows)
    while time()-started < runtime_seconds && amplitude < 13
        target = min(13.0, amplitude*factor)
        ratio = target / amplitude
        println("Starting A0=$target (elapsed $(round(time()-started)) s)")
        flush(stdout)
        solution = try
            solve_nonlinear_mode(previous.frequency, target,
                spin, alpha, modes, couplings, polar, azimuthal, grid;
                horizon_guess=ratio .* previous.horizon_amplitudes,
                outer_tail_guess=ratio .* previous.outer_amplitudes[2:end],
                tolerance=1e-9, max_iterations=10, options...)
        catch exception
            println("FAILED A0=$target: ", sprint(showerror, exception))
            flush(stdout)
            factor = sqrt(factor)
            factor < 1.015 && break
            continue
        end
        row = diagnostics(target, solution)
        push!(continuation, (; amplitude=target, solution))
        push!(rows, row)
        checkpoint(continuation, rows)
        println("SAVED ", row, " elapsed=", round(time()-started))
        flush(stdout)
        amplitude, previous = target, solution
        maximum_factor = amplitude < 0.03 ? 2.0 : amplitude < 0.3 ? 1.5 : 1.25
        factor = min(maximum_factor, factor*1.1)
    end
    println("Finished at A0=", amplitude, " with ", length(rows), " saved points.")
end
main()
