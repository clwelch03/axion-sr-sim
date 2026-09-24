for filename in (
    "angular",
    "stationary/nonlinear",
    "tortoise",
    "stationary/radial",
    "stationary/boundary_series",
    "stationary/boundaries",
    "stationary/integration",
    "stationary/matching",
    "stationary/continuation_predictor",
    "observables/energy",
)
    include(joinpath(@__DIR__, filename * ".jl"))
end