# Batch entry point for the Python figure driver. One Julia process handles
# both radial overtones, so compilation is paid only once.
for file in ("angular", "tortoise", "radial", "boundary_series", "boundaries", "integration", "matching")
    include(joinpath(@__DIR__, file * ".jl"))
end

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

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 5 || error("usage: fig13_points.jl SEEDS.csv OUTPUT.csv SPIN STEP DECAY_LENGTHS")
    write_growth_rate_points(ARGS[1], ARGS[2]; spin=parse(Float64, ARGS[3]),
        step=parse(Float64, ARGS[4]), decay_lengths=parse(Float64, ARGS[5]))
end
