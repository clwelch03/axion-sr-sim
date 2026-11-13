# Regression check of the numbers produced by AxionCloud against stored references.
#
#   julia --project=. benchmarks/regression.jl --write-reference     # (re)write benchmarks/reference.toml
#   julia --project=. benchmarks/regression.jl                       # compare against it; exit code 1 on FAIL
#
# Options: --states linear,plateau,...   restrict to named states (default: all)
#          --skip-evolution              skip the short end-to-end evolution
#          --reference PATH              reference file (default benchmarks/reference.toml)
# See benchmarks/README.md for what each check establishes and the tolerance table.
include(joinpath(@__DIR__, "common.jl"))
using .BenchmarkHarness
using .BenchmarkHarness: Quantity, collect_quantities, relative_difference, load_state, selected_specs, run_metadata
using Printf, TOML

const DEFAULT_REFERENCE_PATH = joinpath(@__DIR__, "reference.toml")
const TOLERANCE_PATH = joinpath(@__DIR__, "tolerances.toml")

struct Options
    write_reference::Bool
    state_names::Vector{String}
    skip_evolution::Bool
    reference_path::String
end

function parse_options(arguments)
    write_reference, skip_evolution = false, false
    state_names, reference_path = String[], DEFAULT_REFERENCE_PATH
    remaining = collect(arguments)
    while !isempty(remaining)
        argument = popfirst!(remaining)
        if argument == "--write-reference"
            write_reference = true
        elseif argument == "--skip-evolution"
            skip_evolution = true
        elseif argument == "--states" && !isempty(remaining)
            state_names = split(popfirst!(remaining), ',')
        elseif argument == "--reference" && !isempty(remaining)
            reference_path = popfirst!(remaining)
        else
            error("unknown argument '$argument'; see the header of benchmarks/regression.jl")
        end
    end
    return Options(write_reference, state_names, skip_evolution, reference_path)
end

"""Quantities of all selected states, in a fixed order."""
function collect_all(options::Options)
    table = Pair{String, Quantity}[]
    for spec in selected_specs(options.state_names)
        state = load_state(spec)
        evolution = spec.include_evolution && !options.skip_evolution
        elapsed = @elapsed append!(table, collect_quantities(state; evolution))
        @printf("collected '%s' in %.0f s\n", spec.name, elapsed); flush(stdout)
    end
    return table
end

# --- reference file --------------------------------------------------------------

function write_reference(path::AbstractString, table)
    metadata = run_metadata()
    document = Dict{String, Any}("metadata" => Dict{String, Any}(string(key) => value for (key, value) in pairs(metadata)))
    for (name, quantity) in table
        document[name] = Dict{String, Any}("class" => quantity.class, "anchor" => quantity.anchor, "values" => quantity.values)
    end
    open(path, "w") do io
        println(io, "# Regression reference written by benchmarks/regression.jl --write-reference. Do not edit by hand.")
        TOML.print(io, document; sorted=true)
    end
    println("wrote ", length(table), " quantities to ", path)
end

function read_reference(path::AbstractString)
    document = TOML.parsefile(path)
    return Dict(name => Quantity(Float64.(entry["values"]), entry["class"], entry["anchor"])
                for (name, entry) in document if name != "metadata"), document["metadata"]
end

read_tolerances() = Dict(name => Float64(entry["relative"]) for (name, entry) in TOML.parsefile(TOLERANCE_PATH)["class"])

# --- comparison ------------------------------------------------------------------

"""Print one line per quantity and return the number of failures."""
function compare(table, reference, tolerances, state_names, skip_evolution)
    failures = 0
    current_state = ""
    for (name, quantity) in table
        state_name = first(split(name, '/'))
        if state_name != current_state
            current_state = state_name
            println("\n== state '$state_name' ==")
            @printf("%-50s %-20s %11s %11s  %s\n", "quantity", "class", "rel. diff", "tolerance", "")
        end
        haskey(reference, name) || (println("  missing from the reference: $name"); failures += 1; continue)
        stored = reference[name]
        difference = relative_difference(quantity.values, stored.values, stored.anchor)
        class = stored.class
        haskey(tolerances, class) || error("class '$class' of $name is not in $TOLERANCE_PATH")
        tolerance = tolerances[class]
        passed = difference <= tolerance
        failures += !passed
        item = join(split(name, '/')[2:end], '/')
        note = class == "reported" ? @sprintf("(reported: %d, reference %d)", quantity.values[1], stored.values[1]) : ""
        @printf("%-50s %-20s %11.3e %11.3e  %s %s\n", item, class, difference, tolerance, passed ? "ok" : "FAIL", note)
    end
    names = Set(first(pair) for pair in table)
    for name in keys(reference)
        first(split(name, '/')) in state_names && !(name in names) && !(skip_evolution && contains(name, "/evolution.")) &&
            (println("  not produced by this run: $name"); failures += 1)
    end
    return failures
end

function main(arguments)
    options = parse_options(arguments)
    table = collect_all(options)
    if options.write_reference
        write_reference(options.reference_path, table)
        return 0
    end
    reference, metadata = read_reference(options.reference_path)
    println("\nreference written on ", get(metadata, "hostname", "?"), " (", get(metadata, "cpu_model", "?"), "), commit ",
            get(metadata, "commit", "?"))
    state_names = [spec.name for spec in selected_specs(options.state_names)]
    failures = compare(table, reference, read_tolerances(), state_names, options.skip_evolution)
    println(failures == 0 ? "\nPASS: all $(length(table)) quantities within tolerance" :
                            "\nFAIL: $failures of $(length(table)) quantities outside tolerance")
    return failures == 0 ? 0 : 1
end

exit(main(ARGS))
