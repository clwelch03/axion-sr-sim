# Slower, standalone check: julia --project=. tests/adiabatic/single_cloud/check_nonlinear_convergence.jl
println("Julia started; loading adiabatic source files.")
flush(stdout)
using Test
using LinearAlgebra: norm

include(joinpath(@__DIR__, "../../../src/adiabatic/load.jl"))
println("Source files loaded; starting convergence check.")
flush(stdout)

@testset "Nonlinear radial step and endpoint convergence" begin
    grid = angular_grid(16, 31)
    polar_values = build_polar_values(5, grid.cosine_theta_nodes)
    azimuthal_values = build_azimuthal_values(5, grid.azimuth_nodes)
    couplings = build_cosine_squared_couplings(5, polar_values, grid)
    modes = [(1, 1), (3, 1), (5, 1), (3, 3), (5, 3), (5, 5)]
    spin, alpha, amplitude = 0.99, 0.42, 1e-4
    settings = (
        (label="baseline", step=0.05, inner=-100.0, outer=100.0),
        (label="half_step", step=0.025, inner=-100.0, outer=100.0),
        (label="inner_farther", step=0.05, inner=-140.0, outer=100.0),
        (label="outer_farther", step=0.05, inner=-100.0, outer=150.0),
        (label="outer_farther_again", step=0.05, inner=-100.0, outer=200.0),
    )
    results = NamedTuple[]
    output_path = joinpath(@__DIR__, "../../../runs/adiabatic/small_amplitude/results/nonlinear_convergence.csv")
    mkpath(dirname(output_path))
    open(output_path, "w") do output
        println(output, "case,amplitude,spin,alpha,step,inner,outer,match,linear_real,linear_imag,nonlinear_real,nonlinear_imag,shift_real,shift_imag,scaled_residual,iterations")
        for setting in settings
            println("Starting $(setting.label).")
            flush(stdout)
            options = (; setting.step, setting.inner, setting.outer, match=5.0)
            # Recompute the linear reference at identical numerical settings.
            linear = solve_linear_mode(0.4088 + 1e-7im, spin, alpha, [1, 3, 5], 1,
                couplings; options...)
            horizon_seed = amplitude .* vcat(linear.horizon_amplitudes, zeros(ComplexF64, 3))
            outer_seed = amplitude .* vcat(linear.outer_amplitudes, zeros(ComplexF64, 3))
            nonlinear = solve_nonlinear_mode(linear.frequency, amplitude, spin, alpha,
                modes, couplings, polar_values, azimuthal_values, grid;
                horizon_guess=horizon_seed, outer_tail_guess=outer_seed[2:end],
                tolerance=1e-11, options...)
            shift = nonlinear.frequency - linear.frequency
            push!(results, (; setting.label, frequency=nonlinear.frequency, shift))
            @test nonlinear.scaled_residual_norm < 1e-11
            row = (setting.label, amplitude, spin, alpha, setting.step, setting.inner,
                setting.outer, 5.0, real(linear.frequency), imag(linear.frequency),
                real(nonlinear.frequency), imag(nonlinear.frequency), real(shift), imag(shift),
                nonlinear.scaled_residual_norm, nonlinear.iterations)
            println(output, join(row, ','))
            flush(output)
            println("$(setting.label): omega=$(nonlinear.frequency), shift=$shift, " *
                    "scaled residual=$(nonlinear.scaled_residual_norm)")
            flush(stdout)
        end
    end

    for result in results[2:end]
        # Resolve the nonlinear real shift to 0.1%, and the total growth rate
        # to one part per million. This does not constrain the tiny growth shift.
        @test isapprox(real(result.shift), real(results[1].shift); rtol=1e-3, atol=0)
        @test isapprox(imag(result.frequency), imag(results[1].frequency); rtol=1e-6, atol=0)
    end
    println("Saved numerical settings and results to $output_path")
end
