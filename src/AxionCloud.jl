"""
    AxionCloud

Multiphase adiabatic evolution of a superradiant axion cloud around a Kerr
black hole, treating self-interactions non-perturbatively.

The pipeline is: choose bound modes and a potential (`ModeSet`,
`AxionPotential`), solve the coupled radial boundary-value problem at fixed
amplitudes (`solve_radial`), extract charges and fluxes (`cloud_charges`,
`boundary_losses`), and advance the amplitudes in time (`evolve`).
"""
module AxionCloud

using LinearAlgebra
using SparseArrays
using ForwardDiff
using FastGaussQuadrature: gausslegendre
using OrdinaryDiffEqTsit5

include("kerr.jl")
include("potential.jl")
include("modes.jl")
include("angular.jl")
include("field.jl")
include("series.jl")
include("radial.jl")
include("bvp.jl")
include("initial_guess.jl")
include("charges.jl")
include("evolution.jl")
include("calibration.jl")

export KerrBackground, tortoise, offset_from_tortoise, radius
export AxionPotential, FreePotential, CosinePotential, QuarticPotential, ChiralPotential
export BoundMode, ModeSet, FourierComponent, fundamental_channels
export AngularGrid, polar_harmonic, mixing_matrix
export FieldQuadrature, nonlinear_source, reconstruct_field
export HorizonSeries, OuterSeries, endpoint_values, outer_wavenumber
export OuterBoundary, IntegratedOuterData, integrated_outer_values, is_propagating, uses_integrated_data, gauss_legendre_integrate
export RadialMesh, DomainEstimate, domain_estimate, RadialProblem, solve_radial, RadialSolution, profiles, slopes, mode_frequencies, growth_rates
export Charges, cloud_charges, component_charges, outward_flux, boundary_losses, charge_jacobian
export phase_actions, action_flux, action_losses, action_jacobian, amplitude_response
export with_background, amplitude_rates, evolve_fixed_background, evolve_with_spindown
export BackgroundState, BackgroundDerivatives, balance_terms, background_charge_derivatives, metric_work, boyer_lindquist_metric_work,
       EvolutionCache, EvolutionResult, EvolutionDiagnostics,
       AmplitudeTrajectory, amplitude_trajectory
export amplitudes_for_energy, fundamental_energies, with_potential

end
