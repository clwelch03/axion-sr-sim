"""
The first-order radial system in the tortoise coordinate, the mesh it lives on,
and the problem definition that ties geometry, potential, modes and quadrature
together.
"""

"""
    RadialMesh(background; rstar_min=-80, rstar_max=180, step=0.4,
               stretch_from=nothing, stretch_factor=1.03, max_step=8step)

Mesh in the tortoise coordinate, uniform with spacing `step` up to
`stretch_from` and geometrically stretched beyond it (each interval
`stretch_factor` times the previous, capped at `max_step`) until `rstar_max`,
which is always the last node. With `stretch_from = nothing` the whole mesh is
uniform. Stretching lets the outer boundary sit far out, where the profiles are
smooth on scales of tens of `M`, at little cost. Stores the horizon offsets
`s = r − r₊` at the nodes and at the interval midpoints (needed by the MIRK4 stage).
"""
struct RadialMesh{T<:Real}
    rstar::Vector{T}
    offsets::Vector{T}
    midpoint_offsets::Vector{T}
end

function RadialMesh(background::KerrBackground, rstar::AbstractVector)
    midpoints = (rstar[1:end-1] .+ rstar[2:end]) ./ 2
    offsets = offset_from_tortoise.(Ref(background), rstar)
    # The element type follows the background (dual numbers when differentiating through it).
    RadialMesh(convert(Vector{eltype(offsets)}, collect(rstar)), offsets, offset_from_tortoise.(Ref(background), midpoints))
end

function RadialMesh(background::KerrBackground; rstar_min=-80.0, rstar_max=180.0, step=0.4,
                    stretch_from=nothing, stretch_factor=1.03, max_step=8step)
    uniform_end = stretch_from === nothing ? rstar_max : min(stretch_from, rstar_max)
    nodes = collect(range(rstar_min, uniform_end; length=round(Int, (uniform_end - rstar_min) / step) + 1))
    width = step
    while nodes[end] < rstar_max - 1e-9 * abs(rstar_max)
        width = min(width * stretch_factor, max_step)
        push!(nodes, min(nodes[end] + width, rstar_max))
    end
    nodes[end] = rstar_max
    RadialMesh(background, nodes)
end

n_nodes(mesh::RadialMesh) = length(mesh.rstar)
n_intervals(mesh::RadialMesh) = n_nodes(mesh) - 1

"""
    DomainEstimate

Physically motivated tortoise endpoints for a parameter point `(a, α, modes)`,
with the quantities that set them (all from hydrogenic frequency estimates):

* `rstar_min`: where the horizon offset is `horizon_offset_ratio` times the
  gap `r₊ − r₋`, the radius of convergence of the Frobenius series (the linear
  rates are converged to 1e-10 already at a ratio of 1e-4, see `validation/`);
* `containment_radius`: Boyer–Lindquist radius beyond which every bound
  component `r^ν e^{−κr}` has fallen below `containment` times its peak;
* `asymptotic_radius`: radius beyond which the outer series of every
  propagating component converges comfortably, `2 k r ≥ asymptotic_margin |ν|²`
  (capped at `max_radius`, in which case `capped` is true);
* `normalization_radius`: radius beyond which the fastest-decaying *fundamental*
  component is below `normalization_floor` of its peak. The amplitudes are
  imposed through the leading outer coefficients, so past this radius the
  boundary condition no longer constrains the amplitude at the Newton
  tolerance and the system turns singular (the eigenvalue comes out wrong);
* `far_radius`: radius from which the integrated outer data of the
  propagating components start (`2 k r ≥ |η|² + convergence_margin`); it is not
  part of the mesh but sets the cost of the inward integration, which grows
  like `1/k²` towards threshold;
* `rstar_max`: with `integrated_outer` (the default, matching
  `OuterBoundary()`), the tortoise coordinate of the containment radius; with
  the series at the endpoint, that of the larger of the containment and
  asymptotic radii, limited by the normalization radius. `capped` is true when
  the series would have needed a larger radius than that limit (or
  `max_radius`) allowed;
* `slowest_decay`: smallest bound decay constant `κ`, `smallest_wavenumber`:
  smallest propagating `k` (`Inf` if no component propagates), and
  `peak_radius`: largest `ν/κ` among bound components (where to start stretching).
"""
struct DomainEstimate{T<:Real}
    rstar_min::T
    rstar_max::T
    containment_radius::T
    asymptotic_radius::T
    normalization_radius::T
    far_radius::T
    capped::Bool
    slowest_decay::T
    smallest_wavenumber::T
    peak_radius::T
end

"""Solve `ν ln(x/ν) − (x − ν) = ln ε` for the `x = κ r > ν` at which `r^ν e^{−κr}` is `ε` times its peak."""
function _containment_coordinate(long_range_power, containment)
    long_range_power ≤ 0 && return -log(containment)
    target = log(containment)
    excess(x) = long_range_power * log(x / long_range_power) - (x - long_range_power) - target
    low, high = long_range_power, long_range_power - 2target + 2long_range_power * log(2 - 2target)
    while excess(high) > 0
        high *= 2
    end
    for _ in 1:100
        middle = (low + high) / 2
        excess(middle) > 0 ? (low = middle) : (high = middle)
    end
    return high
end

"""
    domain_estimate(background, mass_coupling, mode_set; containment=1e-4, horizon_offset_ratio=1e-6,
                    asymptotic_margin=10, max_radius=5000, normalization_floor=1e-8,
                    integrated_outer=true, convergence_margin=40)

Estimate the numerical domain needed at this parameter point; see `DomainEstimate`.
Component frequencies are the hydrogenic `q·ω` with `ω_i = α(1 − α²/2n_i²)`, so
the estimate is only indicative at large `α` or for strongly deformed clouds.
With the series at the endpoint (`integrated_outer=false`) a component within
`|ω_q − α| ≲ 1e-3 α` of threshold hits the cap; with integrated outer data the
mesh only has to contain the bound components, and a warning is issued when
the far radius of a near-threshold component exceeds `20 max_radius`.
"""
function domain_estimate(background::KerrBackground, mass_coupling, mode_set::ModeSet; containment=1e-4,
                         horizon_offset_ratio=1e-6, asymptotic_margin=10, max_radius=5000.0, normalization_floor=1e-8,
                         integrated_outer::Bool=true, convergence_margin=40.0)
    frequencies = component_frequencies(mode_set, hydrogenic_frequency.(mode_set.modes, mass_coupling))
    containment_radius, asymptotic_radius, normalization_radius, far_radius = 0.0, 0.0, Inf, 0.0
    slowest_decay, smallest_wavenumber, peak_radius = Inf, Inf, 0.0
    for (index, frequency) in enumerate(frequencies)
        frequency > 0 || continue
        if frequency < mass_coupling
            decay = sqrt(mass_coupling^2 - frequency^2)
            long_range_power = -1 + (2frequency^2 - mass_coupling^2) / decay
            containment_radius = max(containment_radius, _containment_coordinate(long_range_power, containment) / decay)
            slowest_decay = min(slowest_decay, decay)
            peak_radius = max(peak_radius, long_range_power / decay)
            if is_fundamental(mode_set.components[index])
                normalization_radius = min(normalization_radius,
                                           _containment_coordinate(long_range_power, normalization_floor) / decay)
            end
        else
            wavenumber = sqrt(frequency^2 - mass_coupling^2)
            long_range_power = (2frequency^2 - mass_coupling^2) / wavenumber
            asymptotic_radius = max(asymptotic_radius, asymptotic_margin * long_range_power^2 / (2wavenumber))
            far_radius = max(far_radius, (long_range_power^2 + convergence_margin) / (2wavenumber))   # |η| = (2ω²−α²)/k
            smallest_wavenumber = min(smallest_wavenumber, wavenumber)
        end
    end
    limit = min(max_radius, normalization_radius)
    if integrated_outer
        capped = false
        outer_radius = containment_radius          # as before, containment is never cut by the limit
        far_radius > 20max_radius && @warn "domain_estimate: a component is close to threshold (k = $(round(smallest_wavenumber, sigdigits=2))); " *
            "its integrated outer data start at r ≈ $(round(Int, far_radius)) and the integration cost grows like 1/k² (k → 0 needs a different expansion)" maxlog=1
    else
        capped = asymptotic_radius > limit
        asymptotic_radius = min(asymptotic_radius, limit)
        outer_radius = max(containment_radius, asymptotic_radius)
        capped && @warn "domain_estimate: a propagating component needs r > $(round(Int, limit)) for its outer series to be asymptotic, " *
                        "beyond the usable radius set by the amplitude normalization or max_radius; its rates will depend on the endpoint" maxlog=1
    end
    rstar_max = tortoise(background, outer_radius - background.rplus)
    rstar_min = tortoise(background, horizon_offset_ratio * background.horizon_gap)
    DomainEstimate(rstar_min, rstar_max, containment_radius, asymptotic_radius, normalization_radius, far_radius, capped,
                   slowest_decay, smallest_wavenumber, peak_radius)
end

"""
    RadialMesh(background, mass_coupling, mode_set; step=0.4, stretch_factor=1.03, max_step=nothing, kwargs...)

Mesh on the domain of `domain_estimate(background, mass_coupling, mode_set; kwargs...)`:
uniform with spacing `step` out to twice the largest bound-profile peak radius,
then geometrically stretched with intervals capped at `max_step` (default: the
smaller of `8 step` and a quarter wavelength of the fastest propagating component).
"""
function RadialMesh(background::KerrBackground, mass_coupling, mode_set::ModeSet; step=0.4, stretch_factor=1.03,
                    max_step=nothing, kwargs...)
    estimate = domain_estimate(background, mass_coupling, mode_set; kwargs...)
    frequencies = component_frequencies(mode_set, hydrogenic_frequency.(mode_set.modes, mass_coupling))
    largest_wavenumber = maximum(sqrt(max(frequency^2 - mass_coupling^2, 0.0)) for frequency in frequencies)
    quarter_wavelength = largest_wavenumber > 0 ? pi / (2largest_wavenumber) : Inf
    cap = max_step === nothing ? min(8step, quarter_wavelength) : max_step
    stretch_from = tortoise(background, 2estimate.peak_radius - background.rplus)
    RadialMesh(background; rstar_min=estimate.rstar_min, rstar_max=estimate.rstar_max, step,
               stretch_from=max(stretch_from, estimate.rstar_min + step), stretch_factor, max_step=cap)
end

"""
    RadialProblem(background, potential, α, mode_set; kwargs...)

Fixed data of the coupled radial boundary-value problem at given spin, mass
coupling `α = Mμ`, retained modes, and numerical settings. The amplitudes are
deliberately *not* part of the problem: they are passed to `solve_radial`, so
one problem serves a whole time evolution.

Keyword arguments: `mesh` (default `RadialMesh(background)`), `n_polar`,
`n_phase` (quadrature sizes), `horizon_order` (default 3) and `outer_order`
(default 7) for the asymptotic expansions, `outer_boundary` (an `OuterBoundary`;
how propagating components get their endpoint data) and the shorthand
`outer_data = :integrated | :series` overriding its `data` field.
"""
struct RadialProblem{T<:Real, Q<:Real, P<:AxionPotential}
    background::KerrBackground{T}
    potential::P
    mass_coupling::T
    mode_set::ModeSet
    quadrature::FieldQuadrature{Q}
    mesh::RadialMesh{T}
    mixing::Vector{Matrix{Q}}   # angular coupling matrix of each component
    horizon_order::Int
    outer_order::Int
    outer_boundary::OuterBoundary
end

function RadialProblem(background::KerrBackground{T}, potential::AxionPotential, mass_coupling, mode_set::ModeSet;
                       mesh::RadialMesh=RadialMesh(background), n_polar::Int=12, n_phase::Int=16,
                       horizon_order::Int=3, outer_order::Int=7, outer_boundary::OuterBoundary=OuterBoundary(),
                       outer_data::Union{Nothing, Symbol}=nothing) where {T}
    quadrature = FieldQuadrature(mode_set; n_polar, n_phase)
    mixing = [mixing_matrix(quadrature.angular, mode_set, index) for index in 1:n_components(mode_set)]
    if outer_data !== nothing
        outer_boundary = OuterBoundary(; data=outer_data, include_bound=outer_boundary.include_bound,
                                       convergence_margin=outer_boundary.convergence_margin, far_order=outer_boundary.far_order,
                                       steps_per_wavelength=outer_boundary.steps_per_wavelength)
    end
    RadialProblem(background, potential, T(mass_coupling), mode_set, quadrature, mesh, mixing, horizon_order, outer_order, outer_boundary)
end

"""
    with_background(problem, background, mass_coupling)

The same problem (modes, potential, quadrature, tortoise mesh in units of the
current mass) on a different Kerr background. Used when the black hole evolves.
The background and mass coupling may be dual numbers, so that the residual can
be differentiated through the background with ForwardDiff.
"""
function with_background(problem::RadialProblem, background::KerrBackground{B}, mass_coupling) where {B}
    mesh = RadialMesh(background, problem.mesh.rstar)
    RadialProblem(background, problem.potential, convert(B, mass_coupling), problem.mode_set,
                  problem.quadrature, mesh, problem.mixing, problem.horizon_order, problem.outer_order, problem.outer_boundary)
end

"""
    with_potential(problem, potential)

The same problem with a different self-interaction potential (same modes,
mesh and quadrature). Used to seed nonlinear solves from free-field ones.
"""
with_potential(problem::RadialProblem, potential::AxionPotential) =
    RadialProblem(problem.background, potential, problem.mass_coupling, problem.mode_set, problem.quadrature,
                  problem.mesh, problem.mixing, problem.horizon_order, problem.outer_order, problem.outer_boundary)

n_channels(problem::RadialProblem) = n_channels(problem.mode_set)
n_modes(problem::RadialProblem) = n_modes(problem.mode_set)

"""
    radial_rhs(problem, s, profiles, slopes, component_frequencies, scale)

Right-hand side of the first-order system at offset `s`:

    dR/dr* = Q
    dQ/dr* = −2rΔ/(r²+a²)² Q − [ΔB R + Δ(a²(ω̃_q²−α²) C R + N)] / (r²+a²)²

`ΔB` is formed as a polynomial so that nothing diverges at the horizon. The
profiles are stored in units of `scale` (see `RadialSolution`), and the
nonlinear source is evaluated on the physical field `scale · R`.
Returns the pair `(dR, dQ)` of vectors over channels.
"""
function radial_rhs(problem::RadialProblem, s, profiles::AbstractVector, slopes::AbstractVector,
                    component_frequencies::AbstractVector, scale)
    background, mode_set, mass_coupling = problem.background, problem.mode_set, problem.mass_coupling
    a = background.a
    r = radius(background, s)
    horizon_function, sum_r2a2 = delta(background, s), r2a2(background, s)
    source = nonlinear_source(problem.quadrature, problem.potential, background, mass_coupling, s, scale .* profiles) ./ scale

    forcing = similar(source)
    for (component_index, component) in enumerate(mode_set.components)
        channels = mode_set.component_channels[component_index]
        frequency = component_frequencies[component_index]
        m = component.m
        coupled = problem.mixing[component_index] * profiles[channels]
        for (local_index, channel) in enumerate(channels)
            l = component.ls[local_index]
            delta_B = (sum_r2a2 * frequency - a * m)^2 + horizon_function * (2a * m * frequency - a^2 * frequency^2 - mass_coupling^2 * r^2 - l * (l + 1))
            forcing[channel] = delta_B * profiles[channel] +
                               horizon_function * (a^2 * (frequency^2 - mass_coupling^2) * coupled[local_index] + source[channel])
        end
    end
    dR = slopes
    dQ = (-2r * horizon_function / sum_r2a2^2) .* slopes .- forcing ./ sum_r2a2^2
    return dR, dQ
end
