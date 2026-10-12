using PhysicalConstants.CODATA2022: NewtonianConstantOfGravitation,
    SpeedOfLightInVacuum, ReducedPlanckConstant
using Unitful

# Unreduced Planck energy: sqrt(hbar * c^5 / G).
const PLANCK_ENERGY = uconvert(u"GeV", sqrt(
    float(ReducedPlanckConstant) * float(SpeedOfLightInVacuum)^5 /
    float(NewtonianConstantOfGravitation)
))

function cloud_mass_fraction(energy_code, decay_constant)
    energy_scale = ustrip(u"GeV", decay_constant) / ustrip(u"GeV", PLANCK_ENERGY)
    return energy_code * energy_scale^2
end

"""
    black_hole_reference_units(alpha, axion_energy)

Return reference `(; time, energy, angular_momentum)` in s, J, and J*s.
Supply positive alpha = G*M₀*mu/(hbar*c), and the positive axion rest energy
mu*c² as a Unitful energy (for example, 3e-13u"eV").

The units are t₀ = G*M₀/c³ = alpha*hbar/(mu*c²), E₀ = M₀*c²,
and J₀,unit = G*M₀²/c = E₀*t₀. The initial BH angular momentum is
chi₀ * J₀,unit. These reference units stay fixed throughout a trajectory.
"""
function black_hole_reference_units(alpha, axion_energy)
    isfinite(alpha) && alpha > 0 ||
        throw(ArgumentError("alpha must be finite and positive"))
    energy_ev = ustrip(u"eV", axion_energy)
    isfinite(energy_ev) && energy_ev > 0 ||
        throw(ArgumentError("axion rest energy must be finite and positive"))

    reference_time = uconvert(u"s", alpha * float(ReducedPlanckConstant) / axion_energy)
    reference_energy = uconvert(u"J", reference_time *
        float(SpeedOfLightInVacuum)^5 / float(NewtonianConstantOfGravitation))
    reference_angular_momentum = uconvert(u"J*s", reference_energy * reference_time)
    return (time=reference_time, energy=reference_energy,
            angular_momentum=reference_angular_momentum)
end
