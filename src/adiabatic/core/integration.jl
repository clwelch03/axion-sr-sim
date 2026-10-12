"""
    integrate_radial(rhs!, initial_state, start, stop; step=0.05, save_profile=false)

Integrate from `start` to `stop` using fixed-step fourth-order Runge–Kutta.

The callback `rhs!(derivative, state, position)` must fill `derivative` without
changing `state`. Return only the final state as a complex floating-point array;
leave `initial_state` unchanged. Integrate in either direction and adjust the
spacing to land on `stop`.

For radial evolution, positions are tortoise coordinates and `state = [R; Q]`
with `Q = dR/dr_*`. Require positive finite `step`, which sets the maximum
spacing, not an error tolerance. Check accuracy by reducing it.

With `save_profile=true`, return `(positions, states)` including both endpoints.
Column `j` of `states` is the state at `positions[j]`; positions follow the
integration direction. The default stores only the final state, as needed
during Newton iterations.
"""
function integrate_radial(rhs!, initial_state, start, stop; step=0.05, save_profile=false)
    isfinite(step) && step > 0 ||
        throw(ArgumentError("step must be positive and finite"))

    # Adjust the spacing to land exactly on the destination
    # It becomes negative when integrating inward
    count = max(1, ceil(Int, abs(stop - start) / step))
    spacing = (stop - start) / count
    state = complex.(float.(initial_state))
    k1, k2, k3, k4 = (similar(state) for _ in 1:4)
    positions = save_profile ? collect(range(start, stop; length=count + 1)) : nothing
    states = save_profile ? Matrix{eltype(state)}(undef, length(state), count + 1) : nothing
    save_profile && (states[:, 1] = state)

    for index in 0:(count - 1)
        position = start + index * spacing

        # Four estimates of the derivative across this step
        rhs!(k1, state, position)
        rhs!(k2, state .+ (spacing / 2) .* k1, position + spacing / 2)
        rhs!(k3, state .+ (spacing / 2) .* k2, position + spacing / 2)
        rhs!(k4, state .+ spacing .* k3, position + spacing)

        state .+= (spacing / 6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
        save_profile && (states[:, index + 2] = state)
    end

    return save_profile ? (; positions, states) : state
end
