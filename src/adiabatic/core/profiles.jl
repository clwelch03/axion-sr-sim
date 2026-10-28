"""
    TwoPhaseProfile(frequencies, horizon_coefficients, outer_coefficients)

Physical boundary data for a stationary trial or converged solution.

This container can hold one or more cloud frequencies; individual solvers
determine which cloud counts they support. The corresponding branch settings
or prepared problem specify the cloud and radial-mode ordering.

Construction stores the supplied arrays without copying them, integrating
the radial equations, or checking consistency or convergence.
"""
struct TwoPhaseProfile{T<:Real}
    # One complex base frequency per cloud represented by this profile.
    frequencies::Vector{Complex{T}}

    # Physical boundary coefficients, in the same radial-mode order.
    horizon_coefficients::Vector{Complex{T}}
    outer_coefficients::Vector{Complex{T}}
end