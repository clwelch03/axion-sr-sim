using ForwardDiff

@testset "Analytic phase averages preserve joint tensor and AD" begin
    grid = angular_grid(8, 9)
    polar = build_polar_values(4, grid.cosine_theta_nodes)
    phases = (grid.azimuth_nodes .+ 0.31, 2pi .* (0:10) ./ 11 .- 0.29)
    cloud_m, frequencies = [1, 2], [0.4, 0.415]
    components = [fourier_component(q, ls, cloud_m) for (q, ls) in
        [((1,0),[1,3]), ((0,1),[2,4]), ((2,-1),[0,2]), ((-1,2),[3]), ((3,-2),[1])]]
    labels = [(l,c.q) for c in components for l in c.l_values]
    basis = build_two_phase_source_basis(labels, cloud_m, polar, grid, phases...)
    count = length(labels)
    radial = [0.13cis(0.7i)/i for i in 1:count]
    slopes = [(0.01+0.03im*i)*radial[i] for i in 1:count]
    parameters = vcat(real.(radial), imag.(radial), real.(slopes), imag.(slopes), frequencies)
    for spin in (0.0, 0.7, 0.99), potential in (:quadratic, :cosine, :chiral, :quartic, :sextic)
        cache = phase_averaged_observable_cache(spin, 0.42, components, cloud_m, polar, grid, phases, basis)
        function evaluate(p, direct)
            r = complex.(p[1:count], p[count+1:2count])
            q = complex.(p[2count+1:3count], p[3count+1:4count])
            values = direct ? two_phase_observables_at_radius(Dict(zip(labels,r)), Dict(zip(labels,q)),
                8.0, spin, 0.42, p[end-1:end], cloud_m, polar, grid, phases...; potential) :
                phase_averaged_observables_at_radius(r,q,8.0,p[end-1:end],cache; potential)
            collect(Base.values(values))
        end
        @test evaluate(parameters, false) ≈ evaluate(parameters, true) rtol=2e-12
        @test ForwardDiff.jacobian(p -> evaluate(p,false), parameters) ≈
            ForwardDiff.jacobian(p -> evaluate(p,true), parameters) rtol=3e-12 atol=1e-13
        # Repeated evaluations with different coefficients must not reuse field data.
        changed = copy(parameters)
        changed[1:4count] .*= 0.017
        @test evaluate(changed,false) ≈ evaluate(changed,true) rtol=2e-12
        @test length(cache.radial_weights) == 1
    end

    cache = phase_averaged_observable_cache(0.7,0.42,components,cloud_m,polar,grid,phases,basis)
    weights = phase_averaged_radius_weights(8.0, cache).potential_weights
    coefficients = vcat(real.(radial),imag.(radial))
    for potential in (:cosine,:chiral)
        explicit = c -> sum(weights .* potential_value.(2 .* (basis.basis*c), Ref(potential)))
        averaged = c -> phase_averaged_potential_energy(c,basis,weights,potential)
        @test ForwardDiff.gradient(averaged,coefficients) ≈ ForwardDiff.gradient(explicit,coefficients) rtol=2e-12
        @test ForwardDiff.hessian(averaged,coefficients) ≈ ForwardDiff.hessian(explicit,coefficients) rtol=2e-12
    end
    @test_throws DimensionMismatch phase_averaged_observables_at_radius(radial,slopes,8.0,[0.4],cache)
    @test_throws ArgumentError phase_averaged_observables_at_radius(radial,slopes,8.0,complex.(frequencies),cache)
    @test_throws DimensionMismatch phase_averaged_observables_at_radius(radial[2:end],slopes,8.0,frequencies,cache)

    # The two driven labels differ by (3,-3), so three-point grids alias them.
    aliased = (2pi .* (0:2)./3, 2pi .* (0:2)./3)
    @test_throws ArgumentError phase_averaged_observable_cache(0.7,0.42,components,cloud_m,polar,grid,aliased,basis)
    @test_throws ArgumentError phase_averaged_observable_cache(0.7,0.42,
        [components;components[1]],cloud_m,polar,grid,phases,basis)
    @test_throws ArgumentError phase_averaged_observable_cache(0.7,0.42,
        [fourier_component((0,0),[0],cloud_m)],cloud_m,polar,grid,phases,basis)
    @test_throws ArgumentError phase_averaged_observable_cache(0.7,0.42,
        [components;fourier_component((-1,0),[1],cloud_m)],cloud_m,polar,grid,phases,basis)
    @test_throws ArgumentError phase_averaged_observable_cache(0.7,0.42,components,
        cloud_m,polar,grid,(Float64[],phases[2]),basis)
    @test_throws DimensionMismatch phase_averaged_observable_cache(0.7,0.42,components,
        cloud_m,polar,grid,(phases[1],),basis)
end

@testset "Averaged integration preserves quadrature, orientation and signs" begin
    grid = angular_grid(8,9)
    polar = build_polar_values(1,grid.cosine_theta_nodes)
    components = [fourier_component((1,0),[1],[1,2])]
    phases = (grid.azimuth_nodes,[0.0])
    basis = build_two_phase_source_basis([(1,(1,0))],[1,2],polar,grid,phases...)
    cache = phase_averaged_observable_cache(0.0,0.42,components,[1,2],polar,grid,phases,basis)
    function branch(start,stop)
        positions=collect(range(start,stop; length=7))
        r=[sqrt((1+x^3)/(4pi*0.4*tortoise_geometry(x,0.0).radius^2)) for x in positions]
        (; positions, states=transpose(hcat(r,0.07im .* r)))
    end
    radial=(; from_horizon=branch(3,5),from_outer=branch(8,5))
    direct = two_phase_integrated_observables(radial,0.0,0.42,[0.4,0.415],components,
        [1,2],polar,grid,phases...; potential=:chiral)
    result = phase_averaged_integrated_observables(radial,[0.4,0.415],cache;potential=:chiral)
    @test collect(values(result)) ≈ collect(values(direct)) rtol=2e-12
    @test result.angular_momentum ≈ (8+8^4/4)-(3+3^4/4) rtol=2e-12
    @test result.horizon_energy_flux ≈ -0.07*(1+3^3) rtol=2e-12
    @test result.outer_energy_flux ≈ 0.07*(1+8^3) rtol=2e-12
    reverse_branch(b)=(; positions=reverse(b.positions),states=reverse(b.states;dims=2))
    reversed=(;from_horizon=reverse_branch(radial.from_horizon),from_outer=reverse_branch(radial.from_outer))
    @test collect(values(phase_averaged_integrated_observables(reversed,[0.4,0.415],cache;potential=:chiral))) ≈ collect(values(result))
    bad=(;from_horizon=branch(3,4),from_outer=radial.from_outer)
    @test_throws ArgumentError phase_averaged_integrated_observables(bad,[0.4,0.415],cache)
    bad_states=copy(radial.from_horizon.states);bad_states[1,1]=NaN
    @test_throws ArgumentError phase_averaged_integrated_observables(
        (;from_horizon=(;radial.from_horizon.positions,states=bad_states),radial.from_outer),[0.4,0.415],cache)
end

@testset "Three independent phase labels use the same averaging kernel" begin
    grid=angular_grid(6,9)
    polar=build_polar_values(3,grid.cosine_theta_nodes)
    cloud_m,frequencies=[1,2,3],[0.4,0.41,0.42]
    components=[fourier_component(q,[l],cloud_m) for (q,l) in
        [((1,0,0),1),((0,1,0),2),((0,0,1),3),((1,1,-1),0)]]
    phases=ntuple(i -> 2pi .* (0:4)./5 .+ 0.13i,3)
    labels=[(only(c.l_values),c.q) for c in components]
    phase_tuples=collect(Iterators.product(phases...))
    samples=[(i,phase) for phase in phase_tuples for i in eachindex(grid.polar_weights)]
    # Independent test construction; production two-phase source builder is unchanged.
    basis_matrix=[polar[(l,abs(sum(q .* cloud_m)))][i] *
        (part==1 ? cos(sum(q .* phase)) : sin(sum(q .* phase)))
        for (i,phase) in samples, (part,(l,q)) in [(part,label) for part in 1:2 for label in labels]]
    weights=[grid.polar_weights[i]/length(phase_tuples) for (i,_) in samples]
    basis=TwoPhaseSourceBasis(labels,basis_matrix,weights,[grid.cosine_theta_nodes[i]^2 for (i,_) in samples])
    cache=phase_averaged_observable_cache(0.7,0.42,components,cloud_m,polar,grid,phases,basis)
    r=ComplexF64[0.2+0.1im,0.1-0.07im,-0.03im,0.02+0.01im]
    slopes=(0.01+0.04im).*r
    position=8.0
    geometry=tortoise_geometry(position,0.7)
    theta_table=build_polar_theta_derivatives(polar,grid.cosine_theta_nodes)
    for potential in (:cosine,:chiral)
        reference=zeros(4)
        for (index,(i,phase)) in enumerate(samples)
            field=0.0;derivatives=zeros(4)
            for (j,(l,q)) in enumerate(labels)
                m=sum(q .* cloud_m);omega=sum(q .* frequencies)
                wave=cis(-sum(q .* phase))
                y=polar[(l,abs(m))][i]
                field+=2real(r[j]*wave)*y
                derivatives .+= 2 .* real.([-im*omega*r[j]*y,
                    slopes[j]/geometry.dr_drstar*y,r[j]*theta_table[(l,abs(m))][i],im*m*r[j]*y] .* wave)
            end
            metric=kerr_metric(position,0.7,grid.cosine_theta_nodes[i])
            mixed=metric.inverse_metric*scalar_stress_tensor(field,derivatives,metric.metric,metric.inverse_metric,0.42;potential)
            weight=2pi*weights[index]*(geometry.radius^2+0.7^2*grid.cosine_theta_nodes[i]^2)
            reference .+= weight .* [-geometry.dr_drstar*mixed[1,1],geometry.dr_drstar*mixed[1,4],-mixed[2,1],mixed[2,4]]
        end
        result=phase_averaged_observables_at_radius(r,slopes,position,frequencies,cache;potential)
        @test collect(values(result)) ≈ reference rtol=2e-12
    end
end
