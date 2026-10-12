# Load primary_secondary once into the calling module.
if !isdefined(@__MODULE__, :ADIABATIC_PRIMARY_SECONDARY_LOADED)
    include(joinpath(@__DIR__, "../single_cloud/load.jl"))
    include(joinpath(@__DIR__, "background.jl"))
    include(joinpath(@__DIR__, "secondary_mode.jl"))
    include(joinpath(@__DIR__, "secondary_energy.jl"))
    include(joinpath(@__DIR__, "dissipative_mode.jl"))
    include(joinpath(@__DIR__, "radiative_mode.jl"))
    include(joinpath(@__DIR__, "evolution.jl"))
    include(joinpath(@__DIR__, "coupled_balance.jl"))
    include(joinpath(@__DIR__, "scalar_energy_accounting.jl"))
    include(joinpath(@__DIR__, "continuation.jl"))
    include(joinpath(@__DIR__, "secondary_tables.jl"))
    include(joinpath(@__DIR__, "dissipative_tables.jl"))
    include(joinpath(@__DIR__, "radiative_tables.jl"))
    include(joinpath(@__DIR__, "stationary_point.jl"))
    include(joinpath(@__DIR__, "reference_tables.jl"))
    const ADIABATIC_PRIMARY_SECONDARY_LOADED = true
end
