"""
Hydrogenic starting point for the Newton iteration: each fundamental channel is
seeded with the non-relativistic bound-state profile at the prescribed
amplitude, every driven channel with zero.
"""

"""Generalised Laguerre polynomial `L_k^{(a)}(x)` by three-term recurrence."""
function laguerre(k::Int, a, x)
    k == 0 && return one(x)
    previous, current = one(x), 1 + a - x
    for j in 1:k-1
        previous, current = current, ((2j + 1 + a - x) * current - (j + a) * previous) / (j + 1)
    end
    return current
end

"""
    hydrogenic_profile(mode, mass_coupling, r)

Small-`α` radial profile of `|nlm⟩` normalised to the outer convention
`R → r^{n−1} e^{−κr}`, `κ = α²/n`, so that its leading outer coefficient is one.
"""
function hydrogenic_profile(mode::BoundMode, mass_coupling, r)
    n, l = mode.n, mode.l
    k = n - l - 1
    scaled_r = 2mass_coupling^2 * r / n
    leading_coefficient = (-1)^k / factorial(k) * (2mass_coupling^2 / n)^(n - 1)
    return scaled_r^l * exp(-scaled_r / 2) * laguerre(k, 2l + 1, scaled_r) / leading_coefficient
end

"""
    initial_unknowns(problem, layout, amplitudes, scale)

Unknown vector seeded with hydrogenic profiles and frequencies, zero driven
components, and zero leading boundary coefficients (which the first Newton
step then fills in).
"""
function initial_unknowns(problem::RadialProblem, layout::UnknownLayout, amplitudes::AbstractVector, scale)
    mode_set, background, mass_coupling = problem.mode_set, problem.background, problem.mass_coupling
    unknowns = zeros(n_unknowns(layout))
    n_chan = layout.n_channels
    for (mode_index, mode) in enumerate(mode_set.modes)
        channel = layout.fundamental_channels[mode_index]
        profile(r) = (amplitudes[mode_index] / scale) * hydrogenic_profile(mode, mass_coupling, r)
        for node in 1:layout.n_nodes
            offset = problem.mesh.offsets[node]
            r = radius(background, offset)
            block = node_range(layout, node)
            unknowns[block[4channel - 3]] = profile(r)
            unknowns[block[4channel - 1]] = dr_drstar(background, offset) * ForwardDiff.derivative(profile, r)
        end
    end
    parameter_block = view(unknowns, parameter_range(layout))
    parameter_block[frequency_range(layout)] = _realify(complex.(hydrogenic_frequency.(mode_set.modes, mass_coupling)))
    return unknowns
end
