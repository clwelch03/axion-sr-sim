using Serialization, LinearAlgebra, ForwardDiff, SHA

const EXPERIMENT = get(ENV,"RADIAL_BVP_ROOT",normpath(joinpath(@__DIR__, "..")))
const REPOSITORY = normpath(joinpath(EXPERIMENT, "../../.."))
const REFERENCE_CAMPAIGN = joinpath(REPOSITORY, "runs/adiabatic/multiphase_two_cloud")
const SOURCE_SNAPSHOT = joinpath(REFERENCE_CAMPAIGN, "work/amplitude_continuation_20261028_v1")
include(joinpath(SOURCE_SNAPSHOT, "src/adiabatic/primary_secondary/load.jl"))
include(joinpath(SOURCE_SNAPSHOT, "src/adiabatic/multiphase/load.jl"))

reference_record(name) = deserialize(joinpath(REFERENCE_CAMPAIGN, "results", name * ".jls"))
file_hash(path) = bytes2hex(sha256(read(path)))

"""Interpolate one [R; Q] block, using Q as the derivative of R.

This is only an initial-guess interpolant. It does not smooth across the
shooting join or make a mismatched reference into a solved profile.
"""
function sample_block(radial, position, offset, count)
    branch = position <= last(radial.from_horizon.positions) ?
        radial.from_horizon : radial.from_outer
    forward = first(branch.positions) < last(branch.positions)
    positions = forward ? branch.positions : reverse(branch.positions)
    index = clamp(searchsortedlast(positions, position), 1, length(positions)-1)
    left = forward ? index : length(positions)+1-index
    right = forward ? index+1 : length(positions)-index
    x0, x1 = positions[index], positions[index+1]
    x0-1e-10 <= position <= x1+1e-10 || error("Initial profile does not cover $position")
    h, t = x1-x0, (position-x0)/(x1-x0)
    r0 = branch.states[offset+1:offset+count, left]
    r1 = branch.states[offset+1:offset+count, right]
    q0 = branch.states[offset+count+1:offset+2count, left]
    q1 = branch.states[offset+count+1:offset+2count, right]
    r = (2t^3-3t^2+1).*r0 + (t^3-2t^2+t)*h.*q0 +
        (-2t^3+3t^2).*r1 + (t^3-t^2)*h.*q1
    q = (6t^2-6t)/h.*r0 + (3t^2-4t+1).*q0 +
        (-6t^2+6t)/h.*r1 + (3t^2-2t).*q1
    return r, q
end

"""Keep the current physical boundary series, mode ordering, and amplitudes."""
function endpoint_state(profile, problem, side)
    position = side == :horizon ? problem.integration.inner : problem.integration.outer
    coefficients = side == :horizon ? profile.horizon_coefficients : profile.outer_coefficients
    values, slopes = Vector{eltype(coefficients)}(), Vector{eltype(coefficients)}()
    offset = 0
    for component in problem.components
        count = length(component.l_values)
        state = component_boundary_state(side, position, problem.spin, problem.alpha,
            component, profile.frequencies, problem.cloud_m, problem.couplings,
            coefficients[offset+1:offset+count];
            horizon_order=problem.integration.horizon_order,
            outer_order=problem.integration.outer_order)
        append!(values, state[1:count]); append!(slopes, state[count+1:end])
        offset += count
    end
    return pack_complex(vcat(values, slopes))
end

"""Load a known joint solution, or full reference profiles at the target."""
function initial_case(name)
    input = reference_record("direct_fig7_20261026_v5-input")
    weak = reference_record("ad_seed_20261027_v1-seed")
    problem = input.problem
    count = length(problem.source_basis.mode_labels)
    @assert count == 14
    if name == "weak" || name == "scaled_weak"
        amplitudes = name == "weak" ? weak.cloud_amplitudes : input.cloud_amplitudes
        ratio = amplitudes[2]/weak.cloud_amplitudes[2]
        factors = [ratio^abs(q[2]) for (_,q) in problem.source_basis.mode_labels]
        profile = TwoPhaseProfile(copy(weak.profile.frequencies),
            factors .* weak.profile.horizon_coefficients,
            factors .* weak.profile.outer_coefficients)
        guess = position -> begin
            r,q = sample_block(weak.radial_result, position, 0, count)
            return vcat(factors.*r, factors.*q)
        end
        sources = ["direct_fig7_20261026_v5-input", "ad_seed_20261027_v1-seed"]
    elseif name in ("target_zero", "target_driven")
        amplitudes = input.cloud_amplitudes
        forced = reference_record("driven_field_seed_20261028_v1-forced")
        boundary_seed = reference_record("driven_field_seed_20261028_v1-" *
            (name == "target_driven" ? "both_near" : "zero"))
        profile = boundary_seed.profile
        guess = position -> begin
            r1,q1 = sample_block(forced.primary_profiles, position, 0, 6)
            r2,q2 = sample_block(forced.dissipative.profiles, position, 0, 3)
            r0,q0 = sample_block(forced.dissipative.profiles, position, 6, 3)
            r3,q3 = sample_block(forced.radiative_near.profiles, position, 6, 2)
            driven = name == "target_driven"
            r = vcat(r1, amplitudes[2].*r2,
                (driven ? amplitudes[2] : 0).*r0,
                (driven ? amplitudes[2]^2 : 0).*r3)
            q = vcat(q1, amplitudes[2].*q2,
                (driven ? amplitudes[2] : 0).*q0,
                (driven ? amplitudes[2]^2 : 0).*q3)
            return vcat(r,q)
        end
        sources = ["direct_fig7_20261026_v5-input", "driven_field_seed_20261028_v1-forced",
            "driven_field_seed_20261028_v1-" * (name == "target_driven" ? "both_near" : "zero")]
    else
        error("Unknown case $name")
    end
    # A fixed floor prevents a zero initial driven field from acquiring zero
    # error weight. Record these scales; they are not matching tolerances.
    scale_guess = name == "target_zero" ? initial_case("target_driven").guess : guess
    samples = hcat([scale_guess(x) for x in range(-100,300; length=801)]...)
    peaks = vec(maximum(abs.(samples); dims=2))
    ratio = amplitudes[2]/input.cloud_amplitudes[2]
    floors = [1e-8 * (q[2] == 0 ? 1.0 : ratio) for (_,q) in problem.source_basis.mode_labels]
    field_scales = max.(peaks[1:count], floors)
    slope_scales = max.(peaks[count+1:end], field_scales./100)
    state_scales = repeat(vcat(field_scales, slope_scales), 2)
    packed_reference = pack_two_phase_profile(profile, problem)
    # Normalize each free coefficient by its scaled boundary influence.
    # Frequencies instead use an explicit 0.001 M^-1 displacement scale.
    boundaries = packed -> begin
        trial = unpack_two_phase_profile(packed, amplitudes, problem)
        vcat(endpoint_state(trial, problem, :horizon)./state_scales,
            endpoint_state(trial, problem, :outer)./state_scales)
    end
    boundary_jacobian = ForwardDiff.jacobian(boundaries, packed_reference)
    parameter_scales = 1 ./ max.(vec(sqrt.(sum(abs2,boundary_jacobian; dims=1))),1e-30)
    parameter_scales[[1,2,29,30]] .= 1e-3
    provenance = (; snapshot=SOURCE_SNAPSHOT,
        snapshot_manifest_sha256=file_hash(joinpath(SOURCE_SNAPSHOT,"provenance.json")),
        inputs=Dict(s => file_hash(joinpath(REFERENCE_CAMPAIGN,"results",s*".jls")) for s in sources))
    return (; name, problem, amplitudes, profile, guess, state_scales,
        parameter_scales, packed_reference, provenance)
end

"""Expose the same equations in real scaled coordinates for collocation.

Packing follows pack_complex([R; Q]): [Re R; Re Q; Im R; Im Q].
Parameter zero is the supplied frequency/boundary guess. All scales are fixed
throughout a solve, including mesh refinement.
"""
function collocation_callbacks(case)
    (; problem, amplitudes, state_scales, parameter_scales, packed_reference) = case
    profile_at(p) = unpack_two_phase_profile(packed_reference .+ parameter_scales.*p,
        amplitudes, problem)
    function rhs!(du,u,p,t)
        physical = unpack_complex(state_scales.*view(u,1:56))
        frequencies = profile_at(p).frequencies
        derivative = similar(physical)
        nonlinear_two_phase_rhs!(derivative, physical, t, problem.spin, problem.alpha,
            problem.components, frequencies, problem.cloud_m, problem.couplings,
            problem.source_basis; potential=problem.potential)
        du[1:56] .= pack_complex(derivative)./state_scales
        return nothing
    end
    # MIRK passes its augmented state to BC callbacks, but the original p.
    # Read the current parameter unknowns from that state when present.
    current_parameters(u,p) = length(u) == 112 ? view(u,57:112) : p
    bca! = (res,u,p) -> (res .= view(u,1:56) .- endpoint_state(
        profile_at(current_parameters(u,p)),problem,:horizon)./state_scales; nothing)
    bcb! = (res,u,p) -> (res .= view(u,1:56) .- endpoint_state(
        profile_at(current_parameters(u,p)),problem,:outer)./state_scales; nothing)
    guess = (p,t) -> pack_complex(case.guess(t))./state_scales
    return (; rhs!, bca!, bcb!, guess, profile_at)
end
