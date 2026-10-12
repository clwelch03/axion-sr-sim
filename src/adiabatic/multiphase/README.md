# Mutually deforming independent clouds — planned

This is the home for the proposed general multiphase model. No multiphase
profile solver, action closure, or evolution is implemented yet.

`load.jl` currently loads only `core/`, keeping the dependency boundary usable
without importing the primary/secondary approximation. Add tested components
here as the [implementation plan](../../../docs/multiphase-implementation-plan.md)
progresses. The first physical target is two mutually deforming clouds;
three-cloud evolution will test the additional closure equations.
