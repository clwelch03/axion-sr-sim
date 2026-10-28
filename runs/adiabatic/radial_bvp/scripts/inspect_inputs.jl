using Serialization, LinearAlgebra
root = normpath(joinpath(@__DIR__, "../../../.."))
campaign = joinpath(root, "runs/adiabatic/multiphase_two_cloud")
snapshot = joinpath(campaign, "work/amplitude_continuation_20261028_v1")
include(joinpath(snapshot, "src/adiabatic/primary_secondary/load.jl"))
include(joinpath(snapshot, "src/adiabatic/multiphase/load.jl"))
for name in ("direct_fig7_20261026_v5-input", "ad_seed_20261027_v1-seed",
        "driven_field_seed_20261028_v1-forced")
    record = deserialize(joinpath(campaign, "results", name * ".jls"))
    println(name, ": ", keys(record))
    for key in keys(record)
        value = getproperty(record, key)
        value isa NamedTuple && println("  ", key, ": ", keys(value))
    end
    if hasproperty(record, :problem)
        @show record.problem.integration record.problem.source_basis.mode_labels
        @show record.problem.coefficient_scale
    end
    if hasproperty(record, :profile)
        @show record.profile.frequencies record.profile.horizon_coefficients record.profile.outer_coefficients
    end
    if hasproperty(record, :radial_result)
        @show maximum(abs.(record.radial_result.from_horizon.states); dims=2)
        @show maximum(abs.(record.radial_result.from_outer.states); dims=2)
    end
    if hasproperty(record, :dissipative)
        for key in (:primary_profiles, :dissipative, :radiative_near)
            value = getproperty(record, key)
            radial = key == :primary_profiles ? value : value.profiles
            @show key size(radial.from_horizon.states) size(radial.from_outer.states)
            @show extrema(radial.from_horizon.positions) extrema(radial.from_outer.positions)
        end
    end
end
