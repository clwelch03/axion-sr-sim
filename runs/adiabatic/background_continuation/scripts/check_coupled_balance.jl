# Algebra and saved-data checkpoint, not a time integrator or production library.
using Test, Serialization, LinearAlgebra, SHA
include(joinpath(@__DIR__, "../../../../src/adiabatic/primary_secondary/load.jl"))

# Balance functions are loaded from src/adiabatic/primary_secondary/coupled_balance.jl.
# This driver only assembles saved inputs and checks the equations.

function check_coupled_balance(; output=joinpath(@__DIR__, "../results/balance"))
    mkpath(output)
    primary_path = joinpath(@__DIR__, "../results/amplitude2052/nearby_backgrounds.jls")
    response_path = joinpath(@__DIR__, "../results/amplitude2052/secondary/secondary_backgrounds.jls")
    primary = only(filter(entry -> entry.row.label == "baseline", deserialize(primary_path).entries)).row
    response = only(filter(entry -> entry.row.label == "baseline", deserialize(response_path).entries)).row
    table = read_two_cloud_table(0.42)
    coefficients = two_cloud_coefficients(table)
    primary_energy, secondary_energy = primary.energy, 0.1
    secondary_slope = evaluate_spline(coefficients.secondary_omega, primary_energy).derivative
    balance = cloud_balance_matrix(primary.omega_real, response.omega2_real, secondary_slope, secondary_energy)
    channels = (
        (m=1, frequency=primary.omega_real, horizon=primary.horizon_flux, outward=primary.outer_flux),
        (m=2, frequency=response.omega2_real,
            horizon=secondary_energy*response.secondary_horizon_per_energy,
            outward=secondary_energy*response.secondary_outer_flux/response.secondary_energy),
        (m=0, frequency=response.omega0,
            horizon=secondary_energy*response.m0_horizon_per_energy,
            outward=secondary_energy*response.m0_outer_flux/response.secondary_energy),
        (m=3, frequency=response.omega3,
            horizon=secondary_energy^2*response.m3_horizon_per_energy_squared,
            outward=secondary_energy^2*response.m3_infinity_per_energy_squared),
    )
    fluxes = sum_boundary_channels(channels)
    epsilon = cloud_mass_fraction(1.0, 1e16u"GeV")
    rates = leading_coupled_rates(1.0, 0.99, epsilon, balance, fluxes)
    zero_fluxes = (horizon_energy=0.0, outward_energy=0.0, horizon_torque=0.0, outward_torque=0.0)

    tests = @testset "Coupled balance checkpoint" begin
        @testset "Saved nonlinear fixed-BH limit" begin
            original = two_cloud_rates(coefficients, primary_energy, secondary_energy)
            @test isapprox(rates.primary_rate, original.primary_rate; rtol=1e-11)
            @test isapprox(rates.secondary_rate, original.secondary_rate; rtol=1e-11)
            @test isapprox(fluxes.horizon_energy+fluxes.outward_energy, original.energy_flux; rtol=1e-11)
            @test isapprox(fluxes.horizon_torque+fluxes.outward_torque, original.angular_momentum_flux; rtol=1e-11)
            @test isapprox(rates.mass_rate+epsilon*(rates.primary_rate+rates.secondary_rate),
                           -epsilon*fluxes.outward_energy; atol=1e-22)
            frozen = leading_coupled_rates(1.0, 0.99, 0.0, balance, fluxes)
            @test frozen.mass_rate == frozen.angular_rate == frozen.spin_rate == 0.0
            @test frozen.primary_rate == rates.primary_rate
            @test frozen.secondary_rate == rates.secondary_rate
        end
        @testset "Dilute channel counting and BH signs" begin
            frequency1, frequency2 = 0.4, 0.415
            dilute = cloud_balance_matrix(frequency1, frequency2, 0.0, 0.1)
            # Signed occupations lost/gained per emitted or absorbed scalar.
            for (m, frequency, horizon, outward, occupation_changes) in (
                (1, frequency1, -0.03, 0.0, (-1.0, 0.0)),
                (2, frequency2, -0.03, 0.0, (0.0, -1.0)),
                (0, 2frequency1-frequency2, 0.03, 0.0, (-2.0, 1.0)),
                (3, 2frequency2-frequency1, 0.0, 0.03, (1.0, -2.0)),
            )
                boundary = sum_boundary_channels(((; m, frequency, horizon, outward),))
                result = leading_coupled_rates(1.0, 0.99, epsilon, dilute, boundary)
                number_flux = (horizon+outward)/frequency
                @test isapprox(result.primary_rate, frequency1*occupation_changes[1]*number_flux; atol=1e-16)
                @test isapprox(result.secondary_rate, frequency2*occupation_changes[2]*number_flux; atol=1e-16)
                @test isapprox(result.mass_rate, epsilon*horizon; atol=1e-22)
                @test isapprox(result.angular_rate, epsilon*m/frequency*horizon; atol=1e-22)
                if m == 0
                    @test result.mass_rate > 0 && result.angular_rate == 0 && result.spin_rate < 0
                elseif m == 3
                    @test result.mass_rate == result.angular_rate == result.spin_rate == 0
                else
                    @test result.mass_rate < 0 && result.angular_rate < 0 && result.spin_rate < 0
                end
            end
        end
        @testset "Reference units at nonunit BH mass" begin
            mass_ratio, spin = 0.83, 0.91
            angular_momentum = spin*mass_ratio^2
            result = leading_coupled_rates(mass_ratio, angular_momentum, epsilon, balance, fluxes)
            @test result.primary_rate == rates.primary_rate
            @test result.secondary_rate == rates.secondary_rate
            @test result.mass_rate == rates.mass_rate
            @test isapprox(result.angular_rate, mass_ratio*rates.angular_rate; rtol=1e-14)
            @test isapprox(result.spin_rate, epsilon/mass_ratio*(fluxes.horizon_torque-2spin*fluxes.horizon_energy); rtol=1e-14)
            reference_energies = mass_ratio.*[primary_energy,secondary_energy]
            @test reference_energies./mass_ratio ≈ [primary_energy,secondary_energy]
            code_rates = ([result.primary_rate,result.secondary_rate] .-
                          [primary_energy,secondary_energy].*result.mass_rate)./mass_ratio
            reconstructed = mass_ratio.*code_rates .+ [primary_energy,secondary_energy].*result.mass_rate
            @test reconstructed ≈ [result.primary_rate,result.secondary_rate]
            local_angular_rate = mass_ratio*dot(balance[2,:],[result.primary_rate,result.secondary_rate])
            @test isapprox(result.angular_rate+epsilon*local_angular_rate,
                           -epsilon*mass_ratio*fluxes.outward_torque; atol=1e-22)
        end
        @testset "Moving-background work and angular response" begin
            # Independent, exactly solvable dilute two-oscillator family.
            # These are externally prescribed background rates, not a BH trajectory.
            frequencies(background) = begin
                mass_ratio, angular_momentum = background
                alpha = 0.42mass_ratio
                spin = angular_momentum/mass_ratio^2
                [0.42*(1-alpha^2/8+1e-3spin), 0.42*(1-alpha^2/18+2e-4spin)]
            end
            background = [0.9,0.7]
            background_rates = [-0.002,-0.007]
            occupations = [2.0,3.0]
            frequencies_now = frequencies(background)
            reference_energies = occupations.*frequencies_now
            frequency_rates = ForwardDiff.jacobian(frequencies,background)*background_rates
            expected_rates = occupations.*frequency_rates
            work = sum(expected_rates)
            angular_charge(background_state) =
                sum([1.0,2.0].*reference_energies./frequencies(background_state))
            angular_response = dot(ForwardDiff.gradient(angular_charge,background),background_rates)
            mass_ratio = background[1]
            local_matrix = cloud_balance_matrix((mass_ratio.*frequencies_now)...,0.0,0.1)
            computed = projected_cloud_rates(local_matrix,zero_fluxes,mass_ratio;
                work_per_epsilon=work,background_angular_rate=angular_response)
            @test computed ≈ expected_rates
            @test isapprox(sum(computed),work; atol=1e-17)
            @test isapprox(mass_ratio*dot(local_matrix[2,:],computed)+angular_response,0.0; atol=1e-17)
            wrong = projected_cloud_rates(local_matrix,zero_fluxes,mass_ratio;
                work_per_epsilon=0.0,background_angular_rate=angular_response)
            @test norm(wrong-expected_rates) > 1e-5
            @test norm(wrong./frequencies_now .- occupations.*frequency_rates./frequencies_now) > 1e-5
            leading = projected_cloud_rates(local_matrix,zero_fluxes,mass_ratio;
                work_per_epsilon=0.0,background_angular_rate=0.0)
            @test leading == [0.0,0.0]
            @test mass_ratio*dot(local_matrix[2,:],leading)+angular_response == angular_response
            @test isapprox(epsilon*sum(computed),epsilon*work; atol=1e-22)
        end
    end

    row = (; alpha=0.42, spin=0.99, mass_ratio=1.0, primary_energy, secondary_energy,
        decay_constant_GeV=1e16, epsilon, cloud_mass_fraction=epsilon*(primary_energy+secondary_energy),
        secondary_frequency_energy_slope=secondary_slope,
        balance_determinant=det(balance), balance_condition=cond(balance),
        fluxes..., rates..., alpha_rate=0.42rates.mass_rate)
    open(joinpath(output,"leading_rhs_example.csv"),"w") do io
        println(io,join(keys(row),',')); println(io,join(values(row),','))
    end
    metadata = (; approximation="leading instantaneous background; work and background angular response omitted",
        primary_sha256=bytes2hex(sha256(read(primary_path))),
        response_sha256=bytes2hex(sha256(read(response_path))),
        script_sha256=bytes2hex(sha256(read(@__FILE__))),
        balance_sha256=bytes2hex(sha256(read(joinpath(@__DIR__, "../../../../src/adiabatic/primary_secondary/coupled_balance.jl")))),
        revision=strip(read(`git rev-parse HEAD`,String)),
        table_sha256=[bytes2hex(sha256(read(path))) for path in TWO_CLOUD_TABLES[0.42]])
    serialize(joinpath(output,"balance_checkpoint.jls"),(;metadata,row))
    println("EXAMPLE ",row)
    println("CHECK COUNTS ",Test.get_test_counts(tests))
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    check_coupled_balance()
end
