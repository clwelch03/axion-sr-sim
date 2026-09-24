for filename in (
    "angular",
    "nonlinear",
    "tortoise",
    "radial",
    "boundary_series",
    "boundaries",
    "integration",
    "matching",
    "continuation_predictor",
    "energy",
)
    include(joinpath(@__DIR__, filename * ".jl"))
end