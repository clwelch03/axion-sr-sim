"""
Choosing amplitudes by physical content rather than by the asymptotic
coefficient, e.g. "both clouds carry 650 f_a² M of energy".
"""

"""Free-field energy of each bound mode's fundamental component (see `component_charges`)."""
function fundamental_energies(solution::RadialSolution)
    mode_set = solution.problem.mode_set
    fundamentals = [mode_set.channel_component[channel] for channel in fundamental_channels(mode_set)]
    return [component_charges(solution, component).energy for component in fundamentals]
end

"""
    amplitudes_for_energy(problem, target_energies; initial=nothing, guess=nothing, tol=1e-6,
                          maxiter=40, max_ratio=3.0, verbose=false)

Find amplitudes such that the free-field energy of each fundamental component
equals `target_energies[i]`, returning the converged `RadialSolution`.

Unless `initial` amplitudes are given, the starting point is a linear
estimate: one free-field solve at small amplitudes, rescaled mode by mode with
`A_i ← A_i √(target_i / E_i)` (exact for a linear cloud). The nonlinear
problem is then solved from that estimate and the amplitudes are iterated in
`ln A_i` with a per-mode secant estimate of the logarithmic slope
`d ln E_i / d ln A_i` (which is 2 for a linear cloud but noticeably smaller for
a deformed one). Each step is warm-started from the previous solution and the
per-step amplitude change is limited to a factor `max_ratio` so that the Newton
iteration stays in its basin of convergence.
"""
function amplitudes_for_energy(problem::RadialProblem, target_energies::AbstractVector; initial=nothing,
                               guess=nothing, tol=1e-6, maxiter=40, max_ratio=3.0, verbose=false)
    if initial === nothing
        linear_problem = with_potential(problem, FreePotential())
        linear_solution = solve_radial(linear_problem, fill(1e-6, n_modes(problem)); guess)
        amplitudes = linear_solution.amplitudes .* sqrt.(target_energies ./ fundamental_energies(linear_solution))
        verbose && println("linear estimate: A = $amplitudes")
        guess = linear_solution
    else
        amplitudes = collect(float.(initial))
    end
    solution = solve_radial(problem, amplitudes; guess)
    energies = fundamental_energies(solution)
    log_slopes = fill(2.0, n_modes(problem))   # d ln E_i / d ln A_i, linear value until measured
    for iteration in 1:maxiter
        verbose && println("calibration $iteration: A = $(solution.amplitudes)  E/target = $(energies ./ target_energies)  converged = $(solution.converged)")
        all(abs.(energies ./ target_energies .- 1) .< tol) && solution.converged && return solution
        factors = clamp.(exp.(log.(target_energies ./ energies) ./ log_slopes), 1 / max_ratio, max_ratio)
        previous_solution = solution
        new_amplitudes = amplitudes .* factors
        solution = solve_radial(problem, new_amplitudes; guess=previous_solution)
        # If Newton leaves its basin, halve the step in log-amplitude until it converges.
        backoff = 0
        while !solution.converged && backoff < 4
            backoff += 1
            factors = sqrt.(factors)
            new_amplitudes = amplitudes .* factors
            verbose && println("  step did not converge; backing off to A = $new_amplitudes")
            solution = solve_radial(problem, new_amplitudes; guess=previous_solution)
        end
        new_energies = fundamental_energies(solution)
        # Secant update of the logarithmic slopes, kept within a physically sensible range.
        for mode_index in eachindex(log_slopes)
            amplitude_change = log(new_amplitudes[mode_index] / amplitudes[mode_index])
            energy_change = log(new_energies[mode_index] / energies[mode_index])
            abs(amplitude_change) > 1e-12 && (log_slopes[mode_index] = clamp(energy_change / amplitude_change, 0.3, 3.0))
        end
        amplitudes, energies = new_amplitudes, new_energies
    end
    @warn "amplitudes_for_energy did not reach the target energies to tolerance $tol"
    return solution
end
