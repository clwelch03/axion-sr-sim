# The chord iteration (reused factorisation) must reproduce full Newton along an amplitude sequence.
# Run from the package root:  julia --project=. validation/chord_vs_newton.jl
using AxionCloud
background = KerrBackground(0.9)
mode_set = ModeSet([BoundMode(2,1,1), BoundMode(3,2,2)])
problem = RadialProblem(background, CosinePotential(), 0.3, mode_set)
amplitudes = [2e-2, 1e-4]
reference = solve_radial(problem, amplitudes)
previous_chord, previous_full = reference, reference
for factor in (1.001, 1.005, 1.03)
    amps = amplitudes .* [factor, 1.0]
    t_chord = @elapsed chord = solve_radial(problem, amps; guess=previous_chord)
    t_full = @elapsed full = solve_radial(problem, amps; guess=previous_full, reuse_factorization=false)
    println("factor $factor: chord $(round(t_chord, digits=2))s ($(chord.iterations) it)  full $(round(t_full, digits=2))s ($(full.iterations) it)",
            "  |Δω| = ", maximum(abs, mode_frequencies(chord) .- mode_frequencies(full)),
            "  rel Δ(dA/dt) = ", maximum(abs, amplitude_rates(chord) ./ amplitude_rates(full) .- 1))
    flush(stdout)
    global previous_chord, previous_full = chord, full
end
