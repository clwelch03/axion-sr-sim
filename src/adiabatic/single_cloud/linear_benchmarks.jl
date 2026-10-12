

"""
    write_growth_rate_points(seed_path, output_path;
                             spin=0.99, step=0.05, decay_lengths=15.0)

Solve seeded linear bound modes and write their growth rates and diagnostics to CSV.

Read headerless rows `(radial_overtone, alpha, omega_real, omega_imag)` from
`seed_path`; use the real frequency as the branch seed and initialize growth
at zero. Overwrite `output_path` with a header and one row per seed, flushing
rows and printing progress as each solve completes. Return `nothing`.

Use `spin` in units `M = 1`, positive maximum RK4 spacing `step`, and outer
endpoint `max(150.0, decay_lengths/sqrt(alpha^2-omega_real^2))`. Fix the inner
endpoint at `-120.0` and use `l_values = [1, 3, 5]`, `m = 1`. Throw an error
if the real frequency differs from its seed by `1e-4` or more.
"""
function write_growth_rate_points(seed_path, output_path;
                                  spin=0.99, step=0.05, decay_lengths=15.0)
    grid = angular_grid(16, 31)
    couplings = build_cosine_squared_couplings(
        5, build_polar_values(5, grid.cosine_theta_nodes), grid)
    open(output_path, "w") do output
        println(output, "radial_overtone,alpha,omega_real,omega_imag,residual_norm,iterations,step,inner,outer,spin,l_max")
        for line in eachline(seed_path)
            overtone, alpha, reference_real, _ = parse.(Float64, split(line, ','))
            # Leaver selects the radial branch through the real frequency.
            # Start the imaginary part at zero: the matching solve finds growth.
            guess = complex(reference_real)
            decay_rate = sqrt(alpha^2 - reference_real^2)
            # The overtone at low alpha is more extended than the peak cloud.
            # Place the outer endpoint many exponential decay lengths away.
            outer = max(150.0, decay_lengths / decay_rate)
            solution = solve_linear_mode(guess, spin, alpha, [1, 3, 5], 1, couplings;
                step, inner=-120.0, outer)
            abs(real(solution.frequency) - reference_real) < 1e-4 ||
                error("Adiabatic solve left the seeded branch at alpha=$alpha, n=$overtone")
            row = (Int(overtone), alpha, real(solution.frequency), imag(solution.frequency),
                   solution.residual_norm, solution.iterations, step, -120.0, outer, spin, 5)
            println(output, join(row, ','))
            flush(output)
            println("  adiabatic n=$(Int(overtone)) alpha=$alpha: omega=$(solution.frequency)")
            flush(stdout)
        end
    end
end
