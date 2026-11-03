using Test, ForwardDiff
include(joinpath(@__DIR__, "../../../src/adiabatic/primary_secondary/load.jl"))
include(joinpath(@__DIR__, "../../../runs/adiabatic/multiphase_two_cloud/scripts/comparison_secondary_energy.jl"))

@testset "Potential-consistent comparison energy" begin
    grid = angular_grid(8, 17)
    polar = build_polar_values(3, grid.cosine_theta_nodes)
    azimuthal = build_azimuthal_values(3, grid.azimuth_nodes)
    primary = Dict((1,1)=>0.4+0.07im, (3,1)=>0.04-0.03im, (3,3)=>0.06+0.01im)
    modes = [(2,2), (3,2)]
    state = ComplexF64[0.23+0.08im, 0.05-0.02im, 0.003im, 0.001]
    position, spin, alpha, frequency = 30.0, 0.99, 0.3, 0.298+1e-8im
    args = (state,position,primary,frequency,spin,alpha,modes,polar,azimuthal,grid)
    original = secondary_energy_per_tortoise_radius(args...)
    @test comparison_secondary_energy_density(args...) == original
    @test comparison_secondary_energy_density(args...;potential=:cosine) == original
    free = energy_per_tortoise_radius(state,position,frequency,spin,alpha,modes,
        polar,azimuthal,grid;fundamental_m=2,potential=:quadratic)
    @test comparison_secondary_energy_density(args...;potential=:quadratic) ≈ free rtol=3e-15
    zero_primary = Dict(mode=>zero(value) for (mode,value) in primary)
    @test comparison_secondary_energy_density(state,position,zero_primary,frequency,
        spin,alpha,modes,polar,azimuthal,grid;potential=:chiral) ≈ free rtol=3e-15
    # Independent check: directly average U(phi1+epsilon*phi2) over both phases,
    # then extract its quadratic coefficient using centered finite differences.
    phi1 = reconstruct_field(primary,polar,azimuthal,grid)
    secondary_profile = sum(state[i].*polar[mode] for (i,mode) in enumerate(modes))
    geometry = tortoise_geometry(position,spin)
    rho_squared = geometry.radius^2 .+ spin^2 .* grid.cosine_theta_nodes.^2
    kinetic = energy_per_tortoise_radius(state,position,frequency,spin,0.0,modes,
        polar,azimuthal,grid;fundamental_m=2)
    for potential in (:cosine,:chiral,:quartic)
        expected = comparison_secondary_energy_density(args...;potential)-kinetic
        errors = Float64[]
        for epsilon in (0.05,0.025,0.0125)
            measured = 0.0
            for i in eachindex(grid.polar_weights), j in axes(phi1,2), k in 0:18
                second = 2real(secondary_profile[i]*cis(2pi*k/19))
                first = phi1[i,j]
                coefficient = (potential_value(first+epsilon*second,potential)+
                    potential_value(first-epsilon*second,potential)-2potential_value(first,potential))/(2epsilon^2)
                measured += 2pi*geometry.dr_drstar*alpha^2*grid.polar_weights[i]*
                    rho_squared[i]*grid.azimuth_weight/19*coefficient
            end
            push!(errors,abs(measured-expected)/abs(expected))
        end
        @test errors[end] < 1e-5
        @test 3.5 < errors[1]/errors[2] < 4.5
        @test 3.3 < errors[2]/errors[3] < 4.7
    end
    @test_throws ArgumentError comparison_secondary_energy_density(args...;potential=:invalid)
end
