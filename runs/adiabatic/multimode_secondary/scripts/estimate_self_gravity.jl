# julia --project=. runs/adiabatic/multimode_secondary/scripts/estimate_self_gravity.jl
# Published hydrogenic shifts evaluated on saved cloud energies, with unchanged
# profiles, Kerr parameters, and scalar trajectories. No Poisson/profile solve.
include(joinpath(@__DIR__, "normalize_energies.jl"))

"""
    estimate_self_gravity(alpha, accounting, coefficients, decay_constant; spin=0.99)

Apply Baryakhtar et al., arXiv:2011.11646, Appendix G's rounded coefficients:
delta_omega1/mu = -alpha^2 (0.19 q1 + 0.11 q2),
delta_omega2/mu = -alpha^2 (0.11 q1 + 0.09 q2), q_i = epsilon_f E_i.
The conversion uses M_cloud = mu N at leading nonrelativistic order. The 0.19
coefficient is the rounded 237/1280; retain the published precision throughout.

All returned frequencies and detunings are multiplied by GM_BH (c=hbar=1).
Compare the dilute estimates with the *saved nonlinear* frequencies. This is
an extrapolated sensitivity diagnostic, not a corrected nonlinear solution.
"""
function estimate_self_gravity(alpha, accounting, coefficients, decay_constant; spin=0.99)
    primary_energy = accounting.primary_energy
    primary_omega = [evaluate_spline(coefficients.primary_omega, energy).value
                     for energy in primary_energy]
    secondary_omega = [evaluate_spline(coefficients.secondary_omega, energy).value
                       for energy in primary_energy]
    primary_mass = cloud_mass_fraction.(primary_energy, decay_constant)
    secondary_mass = cloud_mass_fraction.(accounting.secondary_energy, decay_constant)
    primary_self = -alpha^3 .* 0.19 .* primary_mass
    primary_cross = -alpha^3 .* 0.11 .* secondary_mass
    secondary_cross = -alpha^3 .* 0.11 .* primary_mass
    secondary_self = -alpha^3 .* 0.09 .* secondary_mass
    primary_shift = primary_self .+ primary_cross
    secondary_shift = secondary_self .+ secondary_cross
    horizon_omega = spin / (2 * (1 + sqrt(1 - spin^2)))
    return (;
        primary_mass, secondary_mass, primary_omega, secondary_omega,
        primary_self, primary_cross, secondary_self, secondary_cross,
        primary_shift, secondary_shift,
        primary_sr_margin=horizon_omega .- primary_omega,
        secondary_sr_margin=2horizon_omega .- secondary_omega,
        splitting=secondary_omega .- primary_omega,
        splitting_shift=secondary_shift .- primary_shift,
        m3_margin=2secondary_omega .- primary_omega .- alpha,
        m3_shift=2secondary_shift .- primary_shift,
        # m=0 driven field lies below the positive mass threshold on these runs.
        m0_margin=alpha .- (2primary_omega .- secondary_omega),
        m0_shift=2primary_shift .- secondary_shift,
    )
end

function estimate_saved_self_gravity()
    output_directory = joinpath(@__DIR__, "../results")
    summary_rows = []
    println("Published dilute self-gravity shifts on unchanged scalar trajectories.")
    println("No relativistic/nonlinear profile correction; no BH evolution.")
    for (alpha, label) in ((0.42, "alpha042"), (0.3, "alpha030"))
        accounting = scalar_accounting(alpha, label)
        table = read_two_cloud_table(alpha)
        coefficients = two_cloud_coefficients(table)
        units = black_hole_reference_units(alpha, 3e-13u"eV")
        println("\nalpha = ", alpha)
        println("fa [GeV]  peak_cloud/M0  max|dw1|/mu  max|dw2|/mu  max|d_split/split|  max|d_m3/m3_margin|  max|dw1/SR1|  max|dw2/SR2|  exceeds_reservoir")
        for decay_constant in [1e15, 1e16, 1e17, 1e18] .* u"GeV"
            estimate = estimate_self_gravity(alpha, accounting, coefficients, decay_constant)
            (; primary_shift, secondary_shift) = estimate
            # A threshold crossing must be reported explicitly rather than masked
            # by an absolute-value denominator or interpreted as a reliable new state.
            @assert all(estimate.splitting .> 0)
            @assert all(estimate.m3_margin .> 0)
            @assert all(estimate.m0_margin .> 0)
            @assert all(estimate.primary_sr_margin .> 0)
            @assert all(estimate.secondary_sr_margin .> 0)
            report = normalize_accounting(accounting, units, decay_constant, 0.99)
            row = (;
                alpha, decay_constant_GeV=ustrip(u"GeV", decay_constant),
                peak_cloud_fraction=maximum(estimate.primary_mass .+ estimate.secondary_mass),
                max_primary_shift_over_mu=maximum(abs.(primary_shift)) / alpha,
                max_secondary_shift_over_mu=maximum(abs.(secondary_shift)) / alpha,
                max_splitting_fraction=maximum(abs.(estimate.splitting_shift ./ estimate.splitting)),
                max_m3_threshold_fraction=maximum(abs.(estimate.m3_shift ./ estimate.m3_margin)),
                max_primary_sr_fraction=maximum(abs.(primary_shift ./ estimate.primary_sr_margin)),
                max_secondary_sr_fraction=maximum(abs.(secondary_shift ./ estimate.secondary_sr_margin)),
                max_m0_threshold_fraction=maximum(abs.(estimate.m0_shift ./ estimate.m0_margin)),
                splitting_crossed=any(estimate.splitting .+ estimate.splitting_shift .<= 0),
                m3_threshold_crossed=any(estimate.m3_margin .+ estimate.m3_shift .<= 0),
                m0_threshold_crossed=any(estimate.m0_margin .- estimate.m0_shift .<= 0),
                exceeds_reservoir=report.exceeds_reservoir,
            )
            push!(summary_rows, row)
            @printf("%.0e  %.4g  %.4g  %.4g  %.4g  %.4g  %.4g  %.4g  %s\n",
                row.decay_constant_GeV, row.peak_cloud_fraction,
                row.max_primary_shift_over_mu, row.max_secondary_shift_over_mu,
                row.max_splitting_fraction, row.max_m3_threshold_fraction,
                row.max_primary_sr_fraction, row.max_secondary_sr_fraction,
                string(row.exceeds_reservoir))

            # Save the baseline trajectory once; every shift and q scales as fa^2.
            if row.decay_constant_GeV == 1e16
                columns = merge((;
                    time_growth_units=accounting.time_over_mass .* table.linear_growth_rate,
                    time_days=ustrip.(u"d", report.time_days),
                ), estimate)
                path = joinpath(output_directory, "self_gravity_estimate_$(label)_fa1e16.csv")
                write_estimate_columns(path, columns)
                for (name, values) in (("splitting", estimate.splitting),
                                       ("m3 margin", estimate.m3_margin),
                                       ("primary SR margin", estimate.primary_sr_margin))
                    @printf("  %s / mu range: %.5g to %.5g\n", name,
                        minimum(values) / alpha, maximum(values) / alpha)
                end
            end
        end
    end
    summary_path = joinpath(output_directory, "self_gravity_estimate_summary.csv")
    open(summary_path, "w") do io
        println(io, join(string.(keys(first(summary_rows))), ","))
        for row in summary_rows
            println(io, join(values(row), ","))
        end
    end
    return summary_rows
end

function write_estimate_columns(path, columns)
    open(path, "w") do io
        println(io, join(string.(keys(columns)), ","))
        for index in eachindex(first(values(columns)))
            println(io, join((column[index] for column in values(columns)), ","))
        end
    end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    estimate_saved_self_gravity()
end
