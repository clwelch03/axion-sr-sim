# Prepare fresh off-grid solutions, then assess interpolation against them.
# julia --project=. tests/adiabatic/check_branch_interpolation.jl prepare|assess
# The prepared checkpoint is evaluated by branch_energy_derivatives.jl on Slurm.
using Serialization, SHA
include(joinpath(@__DIR__, "../../src/adiabatic/load.jl"))

const interpolation_campaign = normpath(joinpath(@__DIR__,
    "../../runs/adiabatic/amplitude_evolution"))
const interpolation_seed = joinpath(interpolation_campaign,
    "seeds/interpolation_checks_20260924.jls")
const interpolation_retry_seed = joinpath(interpolation_campaign,
    "seeds/interpolation_retry_20260924.jls")

function prepare_interpolation_checks()
    isfile(interpolation_seed) && error("validation seed already exists; reuse it or choose a new path")
    source_path = joinpath(interpolation_campaign, "seeds/cosine_branch_20260924.jls")
    source = deserialize(source_path)
    source_sha256 = bytes2hex(sha256(read(source_path)))
    parts = [deserialize(joinpath(interpolation_campaign,
        "results/cosine_derivatives_part$(index).jls")) for index in 0:3]
    all(part.metadata.input_sha256 == source_sha256 for part in parts) ||
        error("derivative shards do not match the frozen source")
    rows = sort(vcat([part.rows for part in parts]...); by=row -> row.branch_index)
    [row.branch_index for row in rows] == collect(eachindex(source.continuation)) ||
        error("require complete endpoint derivative data")
    entries = Dict(entry.branch_index => entry for part in parts for entry in part.solutions)
    branch = (; amplitudes=[row.amplitude for row in rows], samples=rows)
    continuation, validation_cases = NamedTuple[], NamedTuple[]

    # Three sparse regions and a finely sampled high-amplitude comparison.
    for target in (0.0820125, 2.564697905478641, 6.261469495797465, 11.8)
        left_index = argmin(abs.(branch.amplitudes .- target))
        left, right = rows[left_index], rows[left_index+1]
        left_solution = entries[left.branch_index].solution
        right_solution = entries[right.branch_index].solution
        normalized_state(solution, amplitude) = vcat(solution.frequency,
            solution.horizon_amplitudes ./ amplitude,
            solution.outer_amplitudes[2:end] ./ amplitude)
        left_state = normalized_state(left_solution, left.amplitude)
        right_state = normalized_state(right_solution, right.amplitude)
        mode_count = length(source.modes)
        for fraction in (0.25, 0.5, 0.75)
            amplitude = (1-fraction) * left.amplitude + fraction * right.amplitude
            prediction = (1-fraction) .* left_state .+ fraction .* right_state
            solution = (; frequency=prediction[1],
                horizon_amplitudes=amplitude .* prediction[2:mode_count+1],
                outer_amplitudes=amplitude .* vcat(one(prediction[1]), prediction[mode_count+2:end]))
            push!(continuation, (; amplitude, solution))
            interpolated = interpolate_branch_energy(amplitude, left, right)
            horizon_flux = (1-fraction)*left.horizon_flux + fraction*right.horizon_flux
            outer_flux = (1-fraction)*left.outer_flux + fraction*right.outer_flux
            push!(validation_cases, (; amplitude, fraction,
                left_amplitude=left.amplitude, right_amplitude=right.amplitude,
                energy=interpolated.energy, energy_derivative=interpolated.energy_derivative,
                horizon_flux, outer_flux, total_flux=horizon_flux+outer_flux,
                amplitude_rate=-(horizon_flux + outer_flux) /
                               interpolated.energy_derivative))
        end
    end
    # Guesses are corrected by the production derivative driver. They are not
    # accepted solutions and carry no inferred energy/flux diagnostics.
    serialize(interpolation_seed, merge(source, (; continuation, rows=NamedTuple[],
        validation_cases, source_sha256, relative_target=1e-3)))
    println("Prepared ", length(continuation), " off-grid checks in ", interpolation_seed)
    for row in validation_cases
        println("QUERY ", (; row.amplitude, row.fraction, row.left_amplitude, row.right_amplitude))
    end
end

function prepare_interpolation_retry()
    # Call only after the original tasks stop; preserve all accepted references.
    isfile(interpolation_retry_seed) && error("retry seed already exists")
    seed = deserialize(interpolation_seed)
    completed = Set{Int}()
    for index in 0:3
        path = joinpath(interpolation_campaign, "results/interpolation_checks_part$(index).jls")
        if isfile(path)
            union!(completed, [row.branch_index for row in deserialize(path).rows])
        end
    end
    missing = setdiff(collect(eachindex(seed.validation_cases)), collect(completed))
    isempty(missing) && error("no missing references")
    source_parts = [deserialize(joinpath(interpolation_campaign,
        "results/cosine_derivatives_part$(index).jls")) for index in 0:3]
    entries = Dict(entry.amplitude => entry for part in source_parts for entry in part.solutions)
    continuation = NamedTuple[]
    for index in missing
        query = seed.validation_cases[index]
        left, right = entries[query.left_amplitude], entries[query.right_amplitude]
        # Only the Newton guess changes; observable predictions remain frozen.
        solution = interpolate_branch_state(query.amplitude, left, right)
        push!(continuation, (; amplitude=query.amplitude, solution))
    end
    serialize(interpolation_retry_seed, merge(seed, (; continuation,
        validation_cases=seed.validation_cases[missing], original_indices=missing,
        predictor=:hermite_state)))
    println("Prepared retries for original indices ", missing)
end

function assess_interpolation_checks()
    seed = deserialize(interpolation_seed)
    seed_sha256 = bytes2hex(sha256(read(interpolation_seed)))
    parts = [deserialize(joinpath(interpolation_campaign,
        "results/interpolation_checks_part$(index).jls")) for index in 0:3]
    all(part.metadata.input_sha256 == seed_sha256 for part in parts) ||
        error("validation output differs from prepared seed")
    actual = vcat([part.rows for part in parts]...)
    retry_output = joinpath(interpolation_campaign, "results/interpolation_retry.jls")
    if isfile(retry_output)
        retry_seed = deserialize(interpolation_retry_seed)
        retry = deserialize(retry_output)
        retry.metadata.input_sha256 == bytes2hex(sha256(read(interpolation_retry_seed))) ||
            error("retry output differs from retry seed")
        retry_seed.validation_cases == seed.validation_cases[retry_seed.original_indices] ||
            error("retry changed the interpolation predictions")
        for row in retry.rows
            push!(actual, merge(row, (; branch_index=retry_seed.original_indices[row.branch_index])))
        end
    end
    sort!(actual; by=row -> row.branch_index)
    [row.branch_index for row in actual] == collect(eachindex(seed.validation_cases)) ||
        error("missing or duplicate validation solutions")
    relative_error(estimate, reference) = abs(estimate-reference) / abs(reference)
    comparisons = NamedTuple[]
    for (expected, computed) in zip(seed.validation_cases, actual)
        expected.amplitude == computed.amplitude || error("validation amplitude differs")
        computed.residual <= 1e-10 || error("reference solution did not converge")
        derivative_error = relative_error(expected.energy_derivative, computed.energy_derivative)
        growth_error = relative_error(expected.amplitude_rate, computed.amplitude_rate)
        row = (; amplitude=computed.amplitude, expected.left_amplitude,
            expected.right_amplitude, expected.fraction,
            energy_relative_error=relative_error(expected.energy, computed.energy),
            derivative_interpolated=expected.energy_derivative,
            derivative_reference=computed.energy_derivative,
            derivative_relative_error=derivative_error,
            horizon_flux_relative_error=relative_error(expected.horizon_flux, computed.horizon_flux),
            outer_flux_relative_error=relative_error(expected.outer_flux, computed.outer_flux),
            total_flux_relative_error=relative_error(expected.total_flux, computed.total_flux),
            growth_interpolated=expected.amplitude_rate,
            growth_reference=computed.amplitude_rate,
            growth_relative_error=growth_error,
            derivative_pass=derivative_error <= seed.relative_target,
            growth_pass=growth_error <= seed.relative_target,
            residual=computed.residual)
        push!(comparisons, row)
        println("COMPARISON ", row)
    end
    path = joinpath(interpolation_campaign, "results/interpolation_checks.csv")
    open(path, "w") do io
        println(io, join(string.(keys(first(comparisons))), ','))
        for row in comparisons
            println(io, join(values(row), ','))
        end
    end
    println("Derivative within 0.1%: ", count(row -> row.derivative_pass, comparisons), "/", length(comparisons))
    println("Growth within 0.1%: ", count(row -> row.growth_pass, comparisons), "/", length(comparisons))
    println("Output: ", path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("specify prepare, prepare-retry, or assess")
    if only(ARGS) == "prepare"
        prepare_interpolation_checks()
    elseif only(ARGS) == "assess"
        assess_interpolation_checks()
    elseif only(ARGS) == "prepare-retry"
        prepare_interpolation_retry()
    else
        error("specify prepare, prepare-retry, or assess")
    end
end
