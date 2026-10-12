# PAUSED DRAFT: the stencil job was cancelled for incremental development.
# See ../README.md before using this recipe or its partial output.
# Build or validate the local stationary surrogate; no time integration here.
using Serialization, SHA, Dates
include(joinpath(@__DIR__, "../../../../src/adiabatic/primary_secondary/load.jl"))

function sample_stationary(output_directory; targets_path=nothing)
    repository = normpath(joinpath(@__DIR__, "../../../.."))
    campaign = joinpath(repository,"runs/adiabatic/background_continuation/results/amplitude2052")
    primary_path = joinpath(campaign,"nearby_backgrounds.jls")
    secondary_path = joinpath(campaign,"secondary/secondary_backgrounds.jls")
    primary_saved, secondary_saved = deserialize(primary_path), deserialize(secondary_path)
    initial_primary = only(filter(entry -> entry.row.label == "baseline_half_step",primary_saved.entries))
    initial_secondary = only(filter(entry -> entry.row.label == "baseline_half_step",secondary_saved.entries))
    source_path = joinpath(repository,primary_saved.metadata.seed_path)
    setup = load_secondary_setup(source_path)
    amplitude0 = initial_primary.row.amplitude
    amplitude_radius, alpha_radius, spin_radius = 0.02, 1e-4, 1e-4
    targets = isnothing(targets_path) ? [
        (label="center",amplitude=amplitude0,alpha=0.42,spin=0.99),
        [(label="amplitude_$(offset)",amplitude=amplitude0+offset*amplitude_radius/2,alpha=0.42,spin=0.99)
            for offset in (-2,-1,1,2)]...,
        (label="alpha_minus",amplitude=amplitude0,alpha=0.42-alpha_radius,spin=0.99),
        (label="alpha_plus",amplitude=amplitude0,alpha=0.42+alpha_radius,spin=0.99),
        (label="spin_minus",amplitude=amplitude0,alpha=0.42,spin=0.99-spin_radius),
        (label="spin_plus",amplitude=amplitude0,alpha=0.42,spin=0.99+spin_radius),
    ] : deserialize(targets_path).targets
    modes, secondary_modes = setup.primary.modes, setup.secondary_modes
    dissipative_modes, radiative_modes = [(0,0),(2,0),(4,0)],[(3,3),(5,3),(7,3)]
    grid=setup.grid
    polar=build_polar_values(7,grid.cosine_theta_nodes)
    azimuthal=build_azimuthal_values(7,grid.azimuth_nodes)
    couplings=build_cosine_squared_couplings(7,polar,grid)
    options=merge(setup.primary.options,(;step=0.05))
    metadata=(;amplitude0,amplitude_radius,alpha0=0.42,alpha_radius,spin0=0.99,spin_radius,
        options,primary_tolerance=1e-11,secondary_tolerance=1e-8,targets,
        modes,secondary_modes,dissipative_modes,radiative_modes,
        n_polar=setup.primary.n_polar,n_azimuth=setup.primary.n_azimuth,
        primary_sha256=bytes2hex(sha256(read(primary_path))),
        secondary_sha256=bytes2hex(sha256(read(secondary_path))),
        source_sha256=bytes2hex(sha256(read(source_path))),
        script_sha256=bytes2hex(sha256(read(@__FILE__))),
        revision=strip(read(`git -C $repository rev-parse HEAD`,String)))
    mkpath(output_directory)
    checkpoint=joinpath(output_directory,"stationary.jls")
    samples=if isfile(checkpoint)
        saved=deserialize(checkpoint)
        saved.metadata==metadata || error("stationary resume metadata differs")
        saved.samples
    else
        NamedTuple[]
    end
    println("METADATA ",metadata);flush(stdout)
    function primary_current(branch,frequency,spin)
        count=length(modes)
        state=branch.states[:,1]
        geometry=tortoise_geometry(first(branch.positions),spin)
        return 4π*geometry.radial_factor*sum(real(im*m*frequency*state[i]*conj(state[count+i]))
            for (i,(_,m)) in enumerate(modes))
    end
    for target in targets
        any(sample -> sample.row.label==target.label,samples) && continue
        started=time();println("START ",target," at ",now(UTC));flush(stdout)
        seed=isempty(samples) ? (;row=(amplitude=amplitude0,alpha=0.42,spin=0.99),
            primary=initial_primary.solution,secondary=initial_secondary.secondary) :
            samples[argmin([((sample.row.amplitude-target.amplitude)/amplitude_radius)^2+
                ((sample.row.alpha-target.alpha)/alpha_radius)^2+
                ((sample.row.spin-target.spin)/spin_radius)^2 for sample in samples])]
        ratio=target.amplitude/seed.row.amplitude
        primary=solve_nonlinear_mode(
            complex(real(seed.primary.frequency)*target.alpha/seed.row.alpha,imag(seed.primary.frequency)),
            target.amplitude,target.spin,target.alpha,modes,couplings,polar,azimuthal,grid;
            horizon_guess=ratio.*seed.primary.horizon_amplitudes,
            outer_tail_guess=ratio.*seed.primary.outer_amplitudes[2:end],
            potential=:cosine,jacobian_backend=:forwarddiff,tolerance=metadata.primary_tolerance,options...)
        profiles=nonlinear_matching_residual(primary.frequency,primary.horizon_amplitudes,primary.outer_amplitudes,
            target.spin,target.alpha,modes,couplings,polar,azimuthal,grid;
            potential=:cosine,save_profile=true,options...)
        energy1=sum(integrate_energy_branch(branch,primary.frequency,target.spin,target.alpha,
            modes,polar,azimuthal,grid;potential=:cosine) for branch in (profiles.from_horizon,profiles.from_outer))
        secondary=solve_secondary_mode(
            complex(real(seed.secondary.frequency)*target.alpha/seed.row.alpha,imag(seed.secondary.frequency)),
            profiles,modes,target.spin,target.alpha,secondary_modes,couplings,polar,azimuthal,grid;
            horizon_guess=seed.secondary.horizon_amplitudes,outer_tail_guess=seed.secondary.outer_amplitudes[2:end],
            tolerance=metadata.secondary_tolerance,jacobian_backend=:forwarddiff,reading=:averaged)
        energy2=secondary_energy_and_fluxes(secondary,profiles,modes,target.spin,target.alpha,
            secondary_modes,couplings,polar,azimuthal,grid)
        omega1,omega2=real(primary.frequency),real(secondary.frequency)
        omega0,omega3=2omega1-omega2,2omega2-omega1
        m0=solve_dissipative_mode(omega0,secondary,profiles,modes,target.spin,target.alpha,
            secondary_modes,dissipative_modes,couplings,polar,azimuthal,grid)
        m3=solve_radiative_mode(omega3,secondary,profiles,modes,target.spin,target.alpha,
            secondary_modes,radiative_modes,couplings,polar,azimuthal,grid)
        work0=dissipative_source_work(m0,profiles,modes,omega0,target.spin,target.alpha,
            secondary_modes,dissipative_modes,polar,azimuthal,grid)
        work3=radiative_source_work(m3,profiles,modes,omega3,target.spin,target.alpha,
            secondary_modes,radiative_modes,polar,azimuthal,grid)
        row=(;target...,energy=energy1,omega1,omega2,
            primary_horizon=-primary_current(profiles.from_horizon,primary.frequency,target.spin),
            primary_outward=primary_current(profiles.from_outer,primary.frequency,target.spin),
            secondary_horizon=energy2.horizon_flux/energy2.energy,
            secondary_outward=energy2.outer_flux/energy2.energy,
            m0_horizon=m0.horizon_flux/energy2.energy,m0_outward=m0.outer_flux/energy2.energy,
            m3_horizon=m3.horizon_flux/energy2.energy^2,m3_outward=m3.infinity_flux/energy2.energy^2,
            primary_residual=primary.scaled_residual_norm,secondary_residual=secondary.residual_norm,
            secondary_balance=-(energy2.horizon_flux+energy2.outer_flux)/(2imag(secondary.frequency)*energy2.energy),
            m0_power_ratio=(m0.horizon_flux+m0.outer_flux)/work0,
            m3_power_ratio=(m3.horizon_flux+m3.infinity_flux)/work3,
            m0_relative_residual=m0.residual_norm/m0.source_residual_norm,
            m3_relative_residual=m3.residual_norm/m3.source_residual_norm,
            elapsed_seconds=time()-started)
        row.primary_residual<=metadata.primary_tolerance && row.secondary_residual<=metadata.secondary_tolerance || error("matching failure")
        abs(row.secondary_balance-1)<1e-4 && abs(row.m0_power_ratio-1)<1e-5 && abs(row.m3_power_ratio-1)<1e-5 || error("power balance failure")
        max(row.m0_relative_residual,row.m3_relative_residual)<1e-9 || error("driven matching failure")
        omega3>target.alpha && 0<omega0<target.alpha || error("changed driven boundary type")
        push!(samples,(;row,primary,secondary))
        serialize(checkpoint*".tmp",(;metadata,samples));mv(checkpoint*".tmp",checkpoint;force=true)
        open(joinpath(output_directory,"stationary.csv.tmp"),"w") do io
            println(io,join(keys(row),','))
            for sample in samples
                println(io,join(values(sample.row),','))
            end
        end
        mv(joinpath(output_directory,"stationary.csv.tmp"),joinpath(output_directory,"stationary.csv");force=true)
        println("ACCEPTED ",row);flush(stdout)
    end
    println("PASS: ",length(samples)," stationary samples")
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    1<=length(ARGS)<=2 || error("usage: sample_stationary.jl OUTPUT_DIRECTORY [TARGETS.jls]")
    sample_stationary(ARGS[1];targets_path=length(ARGS)==2 ? ARGS[2] : nothing)
end
