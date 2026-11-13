# Markdown tables from the CSV files of validation/parameter_audit.jl.
#   julia validation/parameter_audit_tables.jl <output-dir> <study> [reference-label]
# Each quantity is given as its relative change from the reference setting (default: the last setting of the
# study that solved, which is the finest), |x − x_ref| / |x_ref|; cost columns are the minimum wall time of the
# perturbed right-hand-side call at 1 and 8 Julia threads, and its ratio to the production setting.
using Printf

const PRODUCTION_LABEL = "np16_npol12_h0.4_sf1.03_c0.0001_hr1e-06_mo3_na1"
const COLUMNS = [("Re ω (max)", ["Re_omega_211", "Re_omega_322"]), ("Im ω₂₁₁", ["Im_omega_211"]), ("Im ω₃₂₂", ["Im_omega_322"]),
                 ("d ln A₂₁₁/dt", ["dlnA211_dt"]), ("d ln A₃₂₂/dt", ["dlnA322_dt"]), ("db/dt", ["db_dt"]), ("dj/dt", ["dj_dt"])]

function read_rows(path)
    lines = readlines(path)
    header = split(lines[1], ",")
    return [Dict(zip(header, split(line, ","; limit=length(header)))) for line in lines[2:end]]
end

value(row, key) = parse(Float64, row[key])
relative(x, reference) = abs(x - reference) / abs(reference)
format_change(change) = change == 0 ? "0" : @sprintf("%.1e", change)

"""The settings tokens of a label that differ from production (`n_phase` np, `n_polar` npol, step h, stretch sf, containment c, horizon ratio hr, max_order mo, n_angular na)."""
function short_label(label)
    label == PRODUCTION_LABEL && return "production"
    parts = split(label, "_")
    tokens = [token for (token, production) in zip(parts, split(PRODUCTION_LABEL, "_")) if token != production]
    append!(tokens, parts[length(split(PRODUCTION_LABEL, "_"))+1:end])   # lattice rule, replaces the tensor grid
    return join(tokens, " ")
end

function main(directory, study, reference_label=nothing)
    rows = read_rows(joinpath(directory, "$study.csv"))
    labels = unique(row["label"] for row in rows)
    states = unique(row["state"] for row in rows)
    println("| state | setting | nodes | ", join(first.(COLUMNS), " | "), " | Δ ln A per 10¹⁰ M | t(1 thr) s | t(8 thr) s | cost vs prod (1/8 thr) |")
    println("|", repeat("---|", 7 + length(COLUMNS)))
    for state in states
        solved = Dict(row["label"] => row for row in rows if row["state"] == state && row["kind"] == "quantities")
        failed = Set(row["label"] for row in rows if row["state"] == state && row["kind"] == "failed")
        timing(label, threads) = (matches = [value(row, "seconds") for row in rows if row["state"] == state && row["kind"] == "timing" &&
                                             row["label"] == label && row["threads"] == string(threads)];
                                  isempty(matches) ? NaN : minimum(matches))
        reference = something(reference_label, last(filter(label -> haskey(solved, label), labels)))
        for label in labels
            if label in failed && !haskey(solved, label)
                println("| $state | $(short_label(label)) | | failed |", repeat(" |", length(COLUMNS) + 3))
                continue
            end
            haskey(solved, label) || continue
            changes = [maximum(relative(value(solved[label], key), value(solved[reference], key)) for key in keys) for (_, keys) in COLUMNS]
            one, eight = timing(label, 1), timing(label, 8)
            ratio = @sprintf("%.2f / %.2f", one / timing(PRODUCTION_LABEL, 1), eight / timing(PRODUCTION_LABEL, 8))
            marker = label == reference ? " (ref)" : ""
            drift = 1e10 * maximum(abs(value(solved[label], key) - value(solved[reference], key)) for key in ("dlnA211_dt", "dlnA322_dt"))
            println("| $state | $(short_label(label))$marker | $(solved[label]["nodes"]) | ", join(format_change.(changes), " | "),
                    " | ", format_change(drift), @sprintf(" | %.3f | %.3f | %s |", one, eight, ratio))
        end
    end
    segments = [row for row in rows if row["kind"] == "segment"]
    isempty(segments) && return
    println("\nSegment 1.5e9 → 1.9e9 from the plateau state (relative change of the final state from the reference setting):\n")
    println("| setting | rhs calls | ln A₂₁₁ (abs) | ln A₃₂₂ (abs) | b (rel) | j (rel) | wall s (8 thr) |")
    println("|---|---|---|---|---|---|---|")
    reference_segment = only(filter(row -> row["label"] == something(reference_label, segments[end]["label"]), segments))
    final(row) = [value(row, key) for key in ("Im_omega_322", "dlnA211_dt", "dlnA322_dt", "db_dt", "dj_dt")]
    for row in segments
        state, base = final(row)[2:5], final(reference_segment)[2:5]
        @printf("| %s | %s | %.1e | %.1e | %.1e | %.1e | %.2f |\n", short_label(row["label"]), row["Re_omega_211"], abs(state[1] - base[1]),
                abs(state[2] - base[2]), relative(state[3], base[3]), relative(state[4], base[4]), value(row, "seconds"))
    end
end

main(ARGS...)
