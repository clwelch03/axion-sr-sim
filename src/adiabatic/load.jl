for filename in (
    "angular",
    "stationary/nonlinear",
    "multimode/background",
    "multimode/secondary_mode",
    "multimode/secondary_energy",
    "multimode/dissipative_mode",
    "multimode/radiative_mode",
    "multimode/two_cloud_evolution",
    "tortoise",
    "stationary/radial",
    "stationary/boundary_series",
    "stationary/boundaries",
    "stationary/integration",
    "stationary/matching",
    "stationary/continuation_predictor",
    "observables/energy",
    "amplitude_evolution/branch_sensitivity",
    "amplitude_evolution/interpolation",
    "amplitude_evolution/time_integration",
)
    include(joinpath(@__DIR__, filename * ".jl"))
end
