# PAUSED DRAFT: not a validated trajectory recipe. See ../README.md.
# Rebuild in small, explained pieces in src before resuming this campaign.
# Limited, fixed-reference-unit time integrator for the leading coupled balance.
include(joinpath(@__DIR__,"../../background_continuation/scripts/check_coupled_balance.jl"))

const LOCAL_FIELDS=(:energy,:omega1,:omega2,:primary_horizon,:primary_outward,
    :secondary_horizon,:secondary_outward,:m0_horizon,:m0_outward,:m3_horizon,:m3_outward)

function local_stationary_model(path)
    saved=deserialize(path)
    rows=[sample.row for sample in saved.samples]
    row_for(label)=only(filter(row -> row.label==label,rows))
    amplitude_rows=sort(filter(row -> row.alpha==0.42 && row.spin==0.99,rows);by=row->row.amplitude)
    length(amplitude_rows)==5 && length(rows)==9 || error("require complete nine-point stencil")
    splines=NamedTuple{LOCAL_FIELDS}(Tuple(cubic_spline([row.amplitude for row in amplitude_rows],
        [getproperty(row,key) for row in amplitude_rows]) for key in LOCAL_FIELDS))
    slopes(direction,radius)=NamedTuple{LOCAL_FIELDS}(Tuple(
        (getproperty(row_for(direction*"_plus"),key)-getproperty(row_for(direction*"_minus"),key))/(2radius)
        for key in LOCAL_FIELDS))
    return (;metadata=saved.metadata,center=row_for("center"),splines,
        alpha_slopes=slopes("alpha",saved.metadata.alpha_radius),
        spin_slopes=slopes("spin",saved.metadata.spin_radius),
        amplitude_min=first(amplitude_rows).amplitude,amplitude_max=last(amplitude_rows).amplitude)
end

function stationary_values(model,amplitude,alpha,spin)
    bounds=model.metadata
    model.amplitude_min<=amplitude<=model.amplitude_max || throw(DomainError(amplitude,"outside amplitude stencil"))
    abs(alpha-bounds.alpha0)<=bounds.alpha_radius || throw(DomainError(alpha,"outside alpha stencil"))
    abs(spin-bounds.spin0)<=bounds.spin_radius || throw(DomainError(spin,"outside spin stencil"))
    values=NamedTuple{LOCAL_FIELDS}(Tuple(evaluate_spline(getproperty(model.splines,key),amplitude).value+
        (alpha-bounds.alpha0)*getproperty(model.alpha_slopes,key)+
        (spin-bounds.spin0)*getproperty(model.spin_slopes,key) for key in LOCAL_FIELDS))
    energy_slope=evaluate_spline(model.splines.energy,amplitude).derivative
    energy_slope>0 || throw(DomainError(energy_slope,"nonmonotonic primary energy"))
    secondary_slope=evaluate_spline(model.splines.omega2,amplitude).derivative/energy_slope
    return (;values,secondary_slope)
end

function amplitude_for_energy(model,energy,alpha,spin)
    low,high=model.amplitude_min,model.amplitude_max
    minimum_energy=stationary_values(model,low,alpha,spin).values.energy
    maximum_energy=stationary_values(model,high,alpha,spin).values.energy
    minimum_energy<=energy<=maximum_energy || throw(DomainError(energy,"outside local energy range"))
    for _ in 1:48
        middle=(low+high)/2
        if stationary_values(model,middle,alpha,spin).values.energy<energy
            low=middle
        else
            high=middle
        end
    end
    return (low+high)/2
end

function instantaneous_rates(values,secondary_slope,mass_ratio,angular_momentum,secondary_energy,epsilon)
    omega0,omega3=2values.omega1-values.omega2,2values.omega2-values.omega1
    channels=(
        (m=1,frequency=values.omega1,horizon=values.primary_horizon,outward=values.primary_outward),
        (m=2,frequency=values.omega2,horizon=secondary_energy*values.secondary_horizon,outward=secondary_energy*values.secondary_outward),
        (m=0,frequency=omega0,horizon=secondary_energy*values.m0_horizon,outward=secondary_energy*values.m0_outward),
        (m=3,frequency=omega3,horizon=secondary_energy^2*values.m3_horizon,outward=secondary_energy^2*values.m3_outward),
    )
    fluxes=sum_boundary_channels(channels)
    balance=cloud_balance_matrix(values.omega1,values.omega2,secondary_slope,secondary_energy)
    cond(balance)<1e5 || throw(DomainError(cond(balance),"ill-conditioned cloud balance"))
    rates=leading_coupled_rates(mass_ratio,angular_momentum,epsilon,balance,fluxes)
    return (;rates,fluxes)
end

function local_rhs(model,state,epsilon;initial_secondary=0.1)
    # Evolve increments to resolve small BH changes without subtracting O(1) numbers.
    mass_ratio=1+state[1]
    angular_momentum=model.metadata.spin0+state[2]
    mass_ratio>0 || throw(DomainError(mass_ratio,"nonpositive BH mass"))
    spin=angular_momentum/mass_ratio^2
    0<spin<1 || throw(DomainError(spin,"outside Kerr spin domain"))
    alpha=model.metadata.alpha0*mass_ratio
    primary_reference=model.center.energy+state[3]
    secondary_reference=initial_secondary+state[4]
    primary_energy,secondary_energy=primary_reference/mass_ratio,secondary_reference/mass_ratio
    primary_energy>0 && secondary_energy>0 || throw(DomainError((primary_energy,secondary_energy),"nonpositive cloud energy"))
    epsilon*(primary_energy+secondary_energy)<0.01 || throw(DomainError(state,"cloud fraction exceeds local-model limit"))
    secondary_energy/primary_energy<0.01 || throw(DomainError(state,"secondary fraction exceeds demonstration limit"))
    amplitude=amplitude_for_energy(model,primary_energy,alpha,spin)
    (;values,secondary_slope)=stationary_values(model,amplitude,alpha,spin)
    0<2values.omega1-values.omega2<alpha && 2values.omega2-values.omega1>alpha ||
        throw(DomainError(state,"changed driven-channel boundary type"))
    (;rates,fluxes)=instantaneous_rates(values,secondary_slope,mass_ratio,angular_momentum,secondary_energy,epsilon)
    derivative=[rates.mass_rate,rates.angular_rate,rates.primary_rate,rates.secondary_rate,
        fluxes.outward_energy,mass_ratio*fluxes.outward_torque]
    all(isfinite,derivative) || error("nonfinite coupled RHS")
    return (;derivative,amplitude,alpha,spin,mass_ratio,angular_momentum,
        primary_reference,secondary_reference,primary_energy,secondary_energy,values,secondary_slope,rates,fluxes)
end

function rk4_coupled_step(rhs,state,step)
    step>0 || throw(ArgumentError("step must be positive"))
    k1=rhs(state)
    k2=rhs(state+step*k1/2)
    k3=rhs(state+step*k2/2)
    k4=rhs(state+step*k3)
    return state+step*(k1+2k2+2k3+k4)/6
end

function trajectory_row(model,state,time,epsilon)
    point=local_rhs(model,state,epsilon)
    return (;time,point.mass_ratio,point.angular_momentum,point.alpha,point.spin,point.amplitude,
        point.primary_reference,point.secondary_reference,point.primary_energy,point.secondary_energy,
        mass_change=state[1],angular_momentum_change=state[2],primary_change=state[3],secondary_change=state[4],
        outward_energy_integral=state[5],outward_torque_integral=state[6],
        cloud_mass_fraction=epsilon*(point.primary_energy+point.secondary_energy),
        energy_audit=state[1]+epsilon*(state[3]+state[4]+state[5]),
        horizon_energy=point.fluxes.horizon_energy,horizon_torque=point.fluxes.horizon_torque,
        outward_energy=point.fluxes.outward_energy,spin_rate=point.rates.spin_rate)
end

function write_local_csv(path,rows)
    open(path*".tmp","w") do io
        println(io,join(keys(first(rows)),','))
        for row in rows
            println(io,join(values(row),','))
        end
    end
    mv(path*".tmp",path;force=true)
end

function integrate_local(model,stencil_path,output_prefix;time_step,final_time=500.0)
    epsilon=cloud_mass_fraction(1.0,1e16u"GeV")
    metadata=(;time_step,final_time,epsilon,decay_constant_GeV=1e16,initial_secondary=0.1,
        initial_primary=model.center.energy,alpha0=model.metadata.alpha0,spin0=model.metadata.spin0,
        approximation="leading light-cloud balance, local cubic-amplitude plus linear-background surrogate",
        stencil_sha256=bytes2hex(sha256(read(stencil_path))),
        script_sha256=bytes2hex(sha256(read(@__FILE__))),
        balance_sha256=bytes2hex(sha256(read(joinpath(@__DIR__,"../../background_continuation/scripts/check_coupled_balance.jl")))),
        revision=strip(read(`git rev-parse HEAD`,String)))
    checkpoint=output_prefix*".jls"
    if isfile(checkpoint)
        saved=deserialize(checkpoint)
        saved.metadata==metadata || error("trajectory resume metadata differs")
        times,states,rows=saved.times,saved.states,saved.rows
    else
        times=[0.0];states=[zeros(6)];rows=[trajectory_row(model,states[1],0.0,epsilon)]
    end
    rhs=state -> local_rhs(model,state,epsilon).derivative
    while last(times)<final_time
        next_time=min(last(times)+time_step,final_time)
        next_state=rk4_coupled_step(rhs,last(states),next_time-last(times))
        row=trajectory_row(model,next_state,next_time,epsilon)
        # An invalid RK stage or endpoint throws before replacing the last checkpoint.
        push!(times,next_time);push!(states,next_state);push!(rows,row)
        serialize(checkpoint*".tmp",(;metadata,times,states,rows))
        mv(checkpoint*".tmp",checkpoint;force=true)
        write_local_csv(output_prefix*".csv",rows)
    end
    println("TRAJECTORY step=",time_step," final=",last(rows));flush(stdout)
    return (;metadata,times,states,rows)
end

function run_local_checkpoint()
    directory=joinpath(@__DIR__,"../results")
    stencil_path=joinpath(directory,"stencil/stationary.jls")
    model=local_stationary_model(stencil_path)
    mkpath(joinpath(directory,"trajectories"))
    runs=[integrate_local(model,stencil_path,joinpath(directory,"trajectories/step$(step)");time_step=step)
          for step in (10.0,5.0,2.5)]
    checks=NamedTuple[]
    function record(name,error,limit)
        push!(checks,(;name,error,limit,passed=isfinite(error)&&error<=limit))
    end
    tests=@testset "Local coupled time integration" begin
        exact_errors=Float64[]
        for steps in (8,16,32)
            state=[1.0]
            for _ in 1:steps
                state=rk4_coupled_step(identity,state,1/steps)
            end
            push!(exact_errors,abs(only(state)-exp(1)))
        end
        @test 14<exact_errors[1]/exact_errors[2]<18
        @test 14<exact_errors[2]/exact_errors[3]<18
        @test rk4_coupled_step(state->zeros(length(state)),[1.0,2.0],0.3)==[1.0,2.0]
        @test_throws DomainError stationary_values(model,model.amplitude_max+0.01,0.42,0.99)
        @test_throws DomainError stationary_values(model,model.metadata.amplitude0,0.421,0.99)
        @test_throws DomainError local_rhs(model,[-1.0,0,0,0,0,0],runs[1].metadata.epsilon)
        for amplitude in (model.amplitude_min,model.metadata.amplitude0,model.amplitude_max)
            energy=stationary_values(model,amplitude,0.42,0.99).values.energy
            @test isapprox(amplitude_for_energy(model,energy,0.42,0.99),amplitude;atol=1e-12)
        end
        for run in runs
            @test last(run.times)==500.0
            error=maximum(abs(row.energy_audit) for row in run.rows)
            record("energy_audit_step$(run.metadata.time_step)",error,1e-17)
            @test error<1e-17
        end
        for index in (1,2), (coordinate,key) in enumerate((:mass_change,:angular_momentum_change,:primary_change,
                :secondary_change,:outward_energy_integral,:outward_torque_integral))
            value=last(runs[index].states)[coordinate]
            reference=last(runs[3].states)[coordinate]
            error=abs((value-reference)/reference)
            record("step$(runs[index].metadata.time_step)_$(key)",error,2e-6)
            @test error<2e-6
        end
    end
    write_local_csv(joinpath(directory,"trajectories/time_checks.csv"),checks)
    fine=last(runs)
    midpoint=only(filter(row -> row.time==250.0,fine.rows))
    endpoint=last(fine.rows)
    targets=[(label="midpoint",amplitude=midpoint.amplitude,alpha=midpoint.alpha,spin=midpoint.spin),
        (label="endpoint",amplitude=endpoint.amplitude,alpha=endpoint.alpha,spin=endpoint.spin),
        (label="endpoint_amplitude_minus",amplitude=endpoint.amplitude-0.001,alpha=endpoint.alpha,spin=endpoint.spin),
        (label="endpoint_amplitude_plus",amplitude=endpoint.amplitude+0.001,alpha=endpoint.alpha,spin=endpoint.spin)]
    serialize(joinpath(directory,"validation_targets.jls"),(;targets,
        trajectory_sha256=bytes2hex(sha256(read(joinpath(directory,"trajectories/step2.5.jls"))))))
    println("PASS: time integration; ",Test.get_test_counts(tests))
    println("Fresh stationary targets written to results/validation_targets.jls")
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    run_local_checkpoint()
end
