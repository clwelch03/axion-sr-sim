# Load the optional MIRK backend into the calling module once.
if !isdefined(@__MODULE__, :ADIABATIC_MULTIPHASE_MIRK_LOADED)
    include(joinpath(@__DIR__, "load.jl"))
    using BoundaryValueDiffEqMIRK, NonlinearSolve, ADTypes, SciMLBase
    include(joinpath(@__DIR__, "mirk_corrections.jl"))
    apply_two_phase_mirk_corrections!()
    include(joinpath(@__DIR__, "collocation_mirk.jl"))
    const ADIABATIC_MULTIPHASE_MIRK_LOADED = true
end
