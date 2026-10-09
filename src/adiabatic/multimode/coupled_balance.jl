# First piece of coupled BH/cloud evolution: instantaneous balance equations.
# No stationary solves, interpolation, time stepping, or file I/O live here.
# Derivation: docs/coupled-bh-cloud-balance.md, sections 1–4.
#
# We use G = c = hbar = 1 and keep the INITIAL BH mass M0 as the reference:
#   tau = t/M0, b = M_BH/M0, j = J_BH/M0^2,
#   epsilon = (fa/M_Planck)^2, q_i = E_i_physical/(epsilon*M0).
# M_Planck is unreduced. All returned rates are derivatives with respect to tau.
# Stationary solvers instead use the CURRENT BH mass: e_i = q_i/b,
# omega_hat = M_BH*omega_physical, chi = j/b^2, alpha = alpha0*b.
"""
    stationary_state(; initial_alpha, mass_ratio, angular_momentum,
                       primary_reference_energy, secondary_reference_energy)

Convert the evolving state from initial-BH units to current-BH units.
Inputs are alpha0 = M0*mu, b = M_BH/M0, j = J_BH/M0^2, and
q_i = E_i_physical/(epsilon*M0), where epsilon = (fa/M_Planck)^2.
Outputs are alpha = M_BH*mu, spin = J_BH/M_BH^2, and
the code cloud energies e_i = E_i_physical/(epsilon*M_BH).

This helper uses our prograde, subextremal domain: 0 <= spin < 1.
Zero cloud energy is allowed here; no stationary solve is performed.
"""
function stationary_state(; initial_alpha, mass_ratio, angular_momentum,
                          primary_reference_energy, secondary_reference_energy)
    inputs = (initial_alpha, mass_ratio, angular_momentum,
              primary_reference_energy, secondary_reference_energy)
    all(isfinite, inputs) || throw(ArgumentError("state inputs must be finite"))
    initial_alpha > 0 || throw(DomainError(initial_alpha, "require alpha0 > 0"))
    mass_ratio > 0 || throw(DomainError(mass_ratio, "require positive BH mass"))
    primary_reference_energy >= 0 ||
        throw(DomainError(primary_reference_energy, "require nonnegative primary energy"))
    secondary_reference_energy >= 0 ||
        throw(DomainError(secondary_reference_energy, "require nonnegative secondary energy"))

    # The particle mass mu stays fixed as the BH mass changes.
    alpha = initial_alpha * mass_ratio
    spin = angular_momentum / mass_ratio^2

    # Change the energy unit from epsilon*M0 to epsilon*M_BH.
    primary_energy = primary_reference_energy / mass_ratio
    secondary_energy = secondary_reference_energy / mass_ratio

    converted = (; alpha, spin, primary_energy, secondary_energy)
    all(isfinite, values(converted)) ||
        throw(DomainError(converted, "conversion produced nonfinite values"))
    0 <= spin < 1 || throw(DomainError(spin, "require 0 <= spin < 1"))
    return converted
end

"""
    cloud_balance_matrix(primary_frequency, secondary_frequency,
                         secondary_slope, secondary_energy)

Return the 2-by-2 matrix multiplying the primary and secondary energy rates.
Frequencies are real frequencies in instantaneous BH units; `secondary_energy`
is e₂ = q₂/b, and `secondary_slope` is ∂omega_hat₂/∂e₁ at fixed alpha and chi.

The first row expresses energy balance. The second expresses angular-momentum
balance for m₁ = 1 and m₂ = 2. Its first entry includes the change of the
secondary angular momentum 2e₂/omega_hat₂ as the primary energy changes.
"""
function cloud_balance_matrix(primary_frequency, secondary_frequency,
                              secondary_slope, secondary_energy)
    primary_angular_response = 1 / primary_frequency -
        2 * secondary_energy * secondary_slope / secondary_frequency^2
    secondary_angular_response = 2 / secondary_frequency
    return [1.0 1.0; primary_angular_response secondary_angular_response]
end

"""
    sum_boundary_channels(channels)

Sum the energy and angular-momentum fluxes at each boundary separately.
Each channel supplies `m`, real instantaneous `frequency`, and signed energy
fluxes `horizon` and `outward`, already scaled to the current cloud energies.
Positive horizon flux enters the BH; positive outward flux leaves the cloud.
Thus a superradiant horizon flux is negative. Each torque is (m/frequency)
times its energy flux. The primary's locked harmonics share that same ratio.

Returned fluxes retain instantaneous code units: physical energy flux is
epsilon times code flux, and physical torque is epsilon*M_BH times code torque.
"""
function sum_boundary_channels(channels)
    return (;
        horizon_energy=sum(channel.horizon for channel in channels),
        outward_energy=sum(channel.outward for channel in channels),
        horizon_torque=sum(channel.m / channel.frequency * channel.horizon for channel in channels),
        outward_torque=sum(channel.m / channel.frequency * channel.outward for channel in channels),
    )
end

"""
    projected_cloud_rates(balance, boundary_fluxes, mass_ratio;
                          work_per_epsilon, background_angular_rate)

Solve the two balance equations for [dq₁/dtau, dq₂/dtau]. Here `mass_ratio` is b.
For the leading approximation both keyword inputs are zero.

The extra inputs let the algebra checkpoint test prescribed metric response:
`work_per_epsilon` is physical metric-work power divided by epsilon;
`background_angular_rate` is K_b*b′ + K_j*j′, where epsilon*M0²*K is cloud
angular momentum and the partial derivatives hold q₁,q₂ fixed. The angular
balance uses this rate divided by b to match instantaneous torque units.
This function does not calculate either response; those calculations remain
unimplemented for our stationary clouds.
"""
function projected_cloud_rates(balance, boundary_fluxes, mass_ratio;
                               work_per_epsilon, background_angular_rate)
    energy_input = -boundary_fluxes.horizon_energy - boundary_fluxes.outward_energy +
        work_per_epsilon
    angular_input = -boundary_fluxes.horizon_torque - boundary_fluxes.outward_torque -
        background_angular_rate / mass_ratio

    # Solve the simultaneous energy and angular-momentum equations.
    return balance \ [energy_input, angular_input]
end

"""
    leading_coupled_rates(mass_ratio, angular_momentum, epsilon, balance, boundary_fluxes)

Return instantaneous BH and cloud rates in fixed initial-BH units.
Inputs `mass_ratio` and `angular_momentum` are b and j (not physical M and J).
Outputs `mass_rate`, `angular_rate`, `spin_rate`, `primary_rate`, and
`secondary_rate` are respectively b′, j′, chi′, q₁′, and q₂′.

This leading small-cloud-mass model omits metric work and background angular
response together. It includes only the channels supplied in `boundary_fluxes`;
it does not add GW losses. In this approximation the energy identity is

    b′ + epsilon*(q₁′ + q₂′) = -epsilon*outward_energy.

This is an instantaneous rate evaluator, not a time integrator or a claim of
exact conservation on a changing Kerr metric.
"""
function leading_coupled_rates(mass_ratio, angular_momentum, epsilon, balance, boundary_fluxes)
    mass_ratio > 0 && epsilon >= 0 || error("require b > 0 and epsilon >= 0")
    spin = angular_momentum / mass_ratio^2

    # Only horizon flux changes the BH directly. Torque has one extra M_BH
    # factor, which becomes b when expressed in fixed M0 units.
    mass_rate = epsilon * boundary_fluxes.horizon_energy
    angular_rate = epsilon * mass_ratio * boundary_fluxes.horizon_torque

    # q_i already removes epsilon from the physical cloud energy. There is
    # no extra epsilon or b multiplying these reference-energy rates.
    cloud_rates = projected_cloud_rates(balance, boundary_fluxes, mass_ratio;
        work_per_epsilon=0.0, background_angular_rate=0.0)

    # Differentiate chi = j/b^2, allowing both j and b to change.
    spin_rate = angular_rate / mass_ratio^2 - 2 * spin * mass_rate / mass_ratio
    return (; mass_rate, angular_rate, spin_rate,
        primary_rate=cloud_rates[1], secondary_rate=cloud_rates[2])
end
