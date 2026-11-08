"""
Axion self-interaction potentials.

With `Φ = a_phys / f_a` and `U = V/(μ² f_a²)` the field equation is
`□Φ − α²Φ + α² h(Φ) = 0`, so a potential is fully specified by its
nonlinearity `h(Φ) = Φ − U'(Φ)` and (for the energy) `U(Φ)` itself. Every
potential is normalised so that `U ≈ Φ²/2` for small `Φ`.
"""

abstract type AxionPotential end

"""Dimensionless potential `U(Φ) = V/(μ² f_a²)`, used when computing the cloud energy."""
function potential_energy end

"""Nonlinearity `h(Φ) = Φ − U'(Φ)`; vanishes for a free field."""
function nonlinearity end

"""
    nonlinearity_derivative(potential, Φ)

`h'(Φ) = 1 − U''(Φ)`, needed for the Jacobian of the nonlinear source. The
generic method differentiates `nonlinearity` automatically; the shipped
potentials override it with closed forms.
"""
nonlinearity_derivative(potential::AxionPotential, field_value) =
    ForwardDiff.derivative(value -> nonlinearity(potential, value), field_value)

"""Free (Klein–Gordon) field, `U = Φ²/2`. Useful for linear-regime checks."""
struct FreePotential <: AxionPotential end
potential_energy(::FreePotential, field_value) = field_value^2 / 2
nonlinearity(::FreePotential, field_value) = zero(field_value)
nonlinearity_derivative(::FreePotential, field_value) = zero(field_value)

"""Standard cosine potential `U = 1 − cos Φ`, giving `h = Φ − sin Φ`."""
struct CosinePotential <: AxionPotential end
potential_energy(::CosinePotential, field_value) = 1 - cos(field_value)
nonlinearity(::CosinePotential, field_value) = field_value - sin(field_value)
nonlinearity_derivative(::CosinePotential, field_value) = 1 - cos(field_value)

"""
    QuarticPotential(λ)

Truncated potential `U = Φ²/2 − λΦ⁴/24`, giving `h = λΦ³/6`. With `λ = 1` this
is the leading self-interaction of the cosine potential.
"""
struct QuarticPotential{T<:Real} <: AxionPotential
    quartic_coupling::T
end
potential_energy(potential::QuarticPotential, field_value) = field_value^2 / 2 - potential.quartic_coupling * field_value^4 / 24
nonlinearity(potential::QuarticPotential, field_value) = potential.quartic_coupling * field_value^3 / 6
nonlinearity_derivative(potential::QuarticPotential, field_value) = potential.quartic_coupling * field_value^2 / 2

"""
    ChiralPotential(; quark_mass_ratio = 0.48)

Chiral QCD axion potential obtained by integrating out the pions,

    V(φ) = m_π² f_π² [1 − √(1 − κ sin²(φ/2f_a))],   κ = 4z/(1+z)²,   z = m_u/m_d.

Normalising by `μ² f_a² = m_π² f_π² z/(1+z)²` gives `U = (4/κ)[1 − √(1 − κ sin²(Φ/2))]`
and `U'(Φ) = sin Φ / √(1 − κ sin²(Φ/2))`. The default `z = 0.48` is the
lattice light-quark mass ratio.
"""
struct ChiralPotential{T<:Real} <: AxionPotential
    quark_mass_ratio::T
    pion_coupling::T   # κ = 4z/(1+z)²
end

function ChiralPotential(; quark_mass_ratio=0.48)
    ratio = quark_mass_ratio
    ChiralPotential(ratio, 4ratio / (1 + ratio)^2)
end

function potential_energy(potential::ChiralPotential, field_value)
    pion_coupling = potential.pion_coupling
    return (4 / pion_coupling) * (1 - sqrt(1 - pion_coupling * sin(field_value / 2)^2))
end

function nonlinearity(potential::ChiralPotential, field_value)
    pion_coupling = potential.pion_coupling
    return field_value - sin(field_value) / sqrt(1 - pion_coupling * sin(field_value / 2)^2)
end

"""`h' = 1 − cos Φ / D − κ sin²Φ / (4D³)` with `D = √(1 − κ sin²(Φ/2))`."""
function nonlinearity_derivative(potential::ChiralPotential, field_value)
    pion_coupling = potential.pion_coupling
    root = sqrt(1 - pion_coupling * sin(field_value / 2)^2)
    return 1 - cos(field_value) / root - pion_coupling * sin(field_value)^2 / (4root^3)
end

"""
    is_even(potential)

Whether `U(−Φ) = U(Φ)`. Even potentials only generate odd Fourier orders in the
phase expansion, which halves the number of components we need to retain.
All potentials shipped here are even.
"""
is_even(::AxionPotential) = true
