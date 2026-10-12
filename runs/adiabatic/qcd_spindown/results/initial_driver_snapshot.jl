# Stationary inputs for the preliminary QCD spin-down campaign. This composes
# the existing Omiya solvers; the primary is nonlinear and the secondary uses
# the phase-averaged primary background. No BH time evolution occurs here.
using Serialization, SHA
include(joinpath(@__DIR__, "../load.jl"))

function qcd_stationary_point(alpha, spin, target_energy; step=0.25, outer=2000.0)
    grid = angular_grid(12, 25)
    polar = build_polar_values(7, grid.cosine_theta_nodes)
    azimuthal = build_azimuthal_values(7, grid.azimuth_nodes)
    couplings = build_cosine_squared_couplings(7, polar, grid)
    modes = [(1,1),(3,1),(5,1),(3,3),(5,3),(5,5)]
    secondary_modes = [(2,2),(4,2),(6,2)]
    options = (; step, inner=-100.0, outer, match=40.0)
    started = time()
    linear = solve_linear_mode(complex(alpha*(1-alpha^2/8)), spin, alpha,
        [1,3,5], 1, couplings; options..., tolerance=1e-7)
    println("LINEAR ",linear.frequency," residual=",linear.residual_norm); flush(stdout)

    # An arbitrarily small amplitude gives the normalization E/A^2. Then
    # correct the actual finite-amplitude solution to the requested energy.
    amplitude = 1e-8
    primary = (; frequency=linear.frequency,
        horizon_amplitudes=amplitude.*vcat(linear.horizon_amplitudes,zeros(ComplexF64,3)),
        outer_amplitudes=amplitude.*vcat(linear.outer_amplitudes,zeros(ComplexF64,3)))
    function profiles_and_energy(solution)
        profiles = nonlinear_matching_residual(solution.frequency,
            solution.horizon_amplitudes,solution.outer_amplitudes,spin,alpha,
            modes,couplings,polar,azimuthal,grid; options...,save_profile=true)
        energy = sum(integrate_energy_branch(branch,solution.frequency,spin,alpha,
            modes,polar,azimuthal,grid) for branch in (profiles.from_horizon,profiles.from_outer))
        return profiles, energy
    end
    profiles, energy = profiles_and_energy(primary)
    for iteration in 1:4
        ratio = sqrt(target_energy/energy)
        amplitude *= ratio
        primary = solve_nonlinear_mode(primary.frequency,amplitude,spin,alpha,
            modes,couplings,polar,azimuthal,grid;
            horizon_guess=ratio.*primary.horizon_amplitudes,
            outer_tail_guess=ratio.*primary.outer_amplitudes[2:end],
            options...,tolerance=1e-7,jacobian_backend=:forwarddiff)
        profiles, energy = profiles_and_energy(primary)
        println("PRIMARY E=",energy," A=",amplitude," residual=",primary.scaled_residual_norm);flush(stdout)
        abs(energy/target_energy-1)<1e-6 && break
    end
    abs(energy/target_energy-1)<1e-6 || error("primary energy inversion failed")

    secondary_linear = solve_linear_mode(complex(alpha*(1-alpha^2/18)),spin,alpha,
        [2,4,6],2,couplings;options...,step=2step,tolerance=1e-5)
    secondary = solve_secondary_mode(secondary_linear.frequency,profiles,modes,
        spin,alpha,secondary_modes,couplings,polar,azimuthal,grid;
        horizon_guess=secondary_linear.horizon_amplitudes,
        outer_tail_guess=secondary_linear.outer_amplitudes[2:end],
        tolerance=1e-5,jacobian_backend=:forwarddiff,reading=:averaged)
    secondary_observables = secondary_energy_and_fluxes(secondary,profiles,modes,
        spin,alpha,secondary_modes,couplings,polar,azimuthal,grid)
    println("SECONDARY ",secondary.frequency," residual=",secondary.residual_norm);flush(stdout)
    omega1, omega2 = real(primary.frequency),real(secondary.frequency)
    dissipative = solve_dissipative_mode(2omega1-omega2,secondary,profiles,modes,
        spin,alpha,secondary_modes,[(0,0),(2,0),(4,0)],couplings,polar,azimuthal,grid)
    radiative = solve_radiative_mode(2omega2-omega1,secondary,profiles,modes,
        spin,alpha,secondary_modes,[(3,3),(5,3),(7,3)],couplings,polar,azimuthal,grid)
    function primary_current(branch)
        state = branch.states[:,1]
        geometry = tortoise_geometry(first(branch.positions),spin)
        return 4π*geometry.radial_factor*sum(real(im*m*primary.frequency*
            state[i]*conj(state[length(modes)+i])) for (i,(_,m)) in enumerate(modes))
    end
    peak_field = maximum(maximum(abs,reconstruct_field(
        Dict(mode=>branch.states[i,k] for (i,mode) in enumerate(modes)),polar,azimuthal,grid))
        for branch in (profiles.from_horizon,profiles.from_outer)
        for k in 1:10:length(branch.positions))
    row = (; alpha,spin,energy,amplitude,omega1,omega2,
        primary_imag=imag(primary.frequency),secondary_imag=imag(secondary.frequency),
        primary_horizon= -primary_current(profiles.from_horizon)/energy,
        primary_outward= primary_current(profiles.from_outer)/energy,
        secondary_horizon=secondary_observables.horizon_flux/secondary_observables.energy,
        secondary_outward=secondary_observables.outer_flux/secondary_observables.energy,
        m0_horizon=dissipative.horizon_flux/(secondary_observables.energy*energy^2),
        m0_outward=dissipative.outer_flux/(secondary_observables.energy*energy^2),
        m3_horizon=radiative.horizon_flux/(secondary_observables.energy^2*energy),
        m3_outward=radiative.infinity_flux/(secondary_observables.energy^2*energy),
        primary_residual=primary.scaled_residual_norm,secondary_residual=secondary.residual_norm,
        m0_relative_residual=dissipative.residual_norm/dissipative.source_residual_norm,
        m3_relative_residual=radiative.residual_norm/radiative.source_residual_norm,
        peak_field,step,outer,elapsed_seconds=time()-started)
    return (;row,primary,secondary)
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    # Header: label,alpha,spin,energy,step,outer. One writer per output folder.
    targets_path, directory = ARGS
    mkpath(directory)
    source_hash = bytes2hex(sha256(read(@__FILE__)))
    for line in readlines(targets_path)[2:end]
        label, numbers... = split(line,',')
        alpha,spin,energy,step,outer = parse.(Float64,numbers)
        path = joinpath(directory,label*".jls")
        if isfile(path)
            saved=deserialize(path)
            saved.source_hash==source_hash || error("source changed; use a new output folder")
            continue
        end
        println("START ",label," ",numbers);flush(stdout)
        result=qcd_stationary_point(alpha,spin,energy;step,outer)
        serialize(path*".tmp",(;source_hash,result...));mv(path*".tmp",path)
        open(joinpath(directory,label*".csv"),"w") do io
            println(io,join(keys(result.row),','));println(io,join(values(result.row),','))
        end
        println("ACCEPTED ",label," ",result.row);flush(stdout)
    end
end
