# Piece 6 (docs/multimode-adiabatic-design.md): evolution of the primary and
# secondary energies by energy and angular-momentum balance (Omiya et al.
# 2211.01949 Eqs. (62)-(63)), for the paper's Fig. 6.
using Printf
include(joinpath(@__DIR__, "../load.jl"))

const MULTIMODE_RESULTS = joinpath(@__DIR__, "../../../runs/adiabatic/multimode_secondary/results")

# Saved tables per branch: primary branch, secondary, secondary energy, m = 0, m = 3.
const TWO_CLOUD_TABLES = Dict(
    0.42 => (joinpath(@__DIR__, "../../../runs/adiabatic/fig2_20260922/results/cosine_branch.csv"),
             joinpath(MULTIMODE_RESULTS, "cosine_secondary.csv"),
             joinpath(MULTIMODE_RESULTS, "cosine_secondary_energy.csv"),
             joinpath(MULTIMODE_RESULTS, "cosine_dissipative_m0.csv"),
             joinpath(MULTIMODE_RESULTS, "cosine_radiative_m3.csv")),
    0.3 => (joinpath(MULTIMODE_RESULTS, "alpha030_cosine_branch.csv"),
            joinpath(MULTIMODE_RESULTS, "alpha030_cosine_secondary.csv"),
            joinpath(MULTIMODE_RESULTS, "alpha030_cosine_secondary_energy.csv"),
            joinpath(MULTIMODE_RESULTS, "alpha030_cosine_dissipative_m0.csv"),
            joinpath(MULTIMODE_RESULTS, "alpha030_cosine_radiative_m3.csv")),
)

"""
    read_csv_columns(path)

Read a numeric CSV with a header line into a dictionary of column vectors.
"""
function read_csv_columns(path)
    lines = readlines(path)
    names = split(first(lines), ",")
    rows = [parse.(Float64, split(line, ",")) for line in lines[2:end]]
    return Dict(name => [row[k] for row in rows] for (k, name) in enumerate(names))
end

"""
    read_two_cloud_table(alpha)

Join the saved piece 1-5 tables for one branch point by point and return the
coefficients of the balance equations at each saved point, as vectors:

- `primary_energy` E₁ and `primary_omega` Re ω₁;
- `linear_growth_rate` Im ω₁ at the first (smallest-amplitude) point, the
  time unit ω_{1,I} of Fig. 6;
- `secondary_omega` Re ω₂;
- `primary_flux_per_energy` F₁/E₁, the primary's own total flux (horizon and
  radiation) per unit E₁; tends to -2 Im ω₁ at small E₁;
- `secondary_flux_per_energy` (F_H + F_out)/Ê₂ of the secondary itself;
  tends to -2 Im ω₂;
- `m0_flux_coefficient` F⁽⁰⁾/(Ê₂ E₁²), the m = 0 flux per unit E₂ and E₁²;
  tends to the paper's F₀ (Eq. (36));
- `m3_flux_coefficient` (F∞ + F_H)⁽³⁾/(Ê₂² E₁), the m = 3 flux per unit E₂²
  and E₁; tends to the paper's F₃.

Dividing each flux by its small-field power makes it tend to a constant as
E₁ → 0, which keeps the splines accurate there. Positive flux means energy
lost. The secondary-dependent columns are independent of its normalization.
"""
function read_two_cloud_table(alpha)
    branch_path, secondary_path, energy_path, m0_path, m3_path = TWO_CLOUD_TABLES[alpha]
    branch = read_csv_columns(branch_path)
    secondary = read_csv_columns(secondary_path)
    energy = read_csv_columns(energy_path)
    m0 = read_csv_columns(m0_path)
    m3 = read_csv_columns(m3_path)

    amplitude = branch["amplitude"]
    for table in (secondary, energy, m0, m3)
        maximum(abs.(table["amplitude"] ./ amplitude .- 1)) < 1e-12 ||
            error("tables for alpha = $alpha are not aligned point by point")
    end

    primary_energy = branch["energy"]
    secondary_energy = energy["energy"]
    return (;
        amplitude,
        primary_energy,
        primary_omega=branch["omega_real"],
        linear_growth_rate=first(branch["omega_imag"]),
        secondary_omega=secondary["omega2_real"],
        primary_flux_per_energy=branch["total_flux"] ./ primary_energy,
        secondary_flux_per_energy=(energy["horizon_flux"] .+ energy["outer_flux"]) ./
                                  secondary_energy,
        m0_flux_coefficient=(m0["horizon_flux"] .+ m0["outer_flux"]) ./
                            (secondary_energy .* primary_energy .^ 2),
        m3_flux_coefficient=(m3["infinity_flux"] .+ m3["horizon_flux"]) ./
                            (secondary_energy .^ 2 .* primary_energy),
    )
end

# Fig. 6 initial data: E1 ~ 0.27 (Mμ = 0.42) or 12 (0.3), and E2 = 1e-8, in Fa² M.
const FIG6_INITIAL_PRIMARY_ENERGY = Dict(0.42 => 0.27, 0.3 => 12.0)
const FIG6_INITIAL_SECONDARY_ENERGY = 1e-8

"""
    fig6_evolution(alpha; time_step=1e-3, final_time=10.0, points=nothing)

Evolve the two clouds from Fig. 6's initial data on the saved branch for
`alpha`, with time in units of the linear growth rate ω_{1,I}. Pass `points`
to build the coefficient splines from a subset of the saved points.
"""
function fig6_evolution(alpha; time_step=1e-3, final_time=10.0, points=nothing)
    table = read_two_cloud_table(alpha)
    coefficients = points === nothing ? two_cloud_coefficients(table) :
                   two_cloud_coefficients(table; points)
    return evolve_two_clouds(
        coefficients, FIG6_INITIAL_PRIMARY_ENERGY[alpha], FIG6_INITIAL_SECONDARY_ENERGY;
        time_unit=table.linear_growth_rate, final_time, time_step,
    )
end

"""
    evolution_features(evolution)

Return the peak of E₁ and of E₂ (value and time) and the final energies.
"""
function evolution_features(evolution)
    (; times, primary_energies, secondary_energies) = evolution
    primary_peak = argmax(primary_energies)
    secondary_peak = argmax(secondary_energies)
    return (; primary_peak=primary_energies[primary_peak], primary_peak_time=times[primary_peak],
              secondary_peak=secondary_energies[secondary_peak],
              secondary_peak_time=times[secondary_peak],
              primary_final=last(primary_energies), secondary_final=last(secondary_energies))
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    # julia --project=. src/adiabatic/drivers/two_cloud_evolution.jl
    # Writes results/two_cloud_evolution_alpha042.csv and ..._alpha030.csv.
    for (alpha, label) in ((0.42, "alpha042"), (0.3, "alpha030"))
        evolution = fig6_evolution(alpha)
        path = joinpath(MULTIMODE_RESULTS, "two_cloud_evolution_$(label).csv")
        open(path, "w") do io
            println(io, "time,primary_energy,secondary_energy")
            for (time, primary, secondary) in zip(evolution.times[1:10:end],
                                                  evolution.primary_energies[1:10:end],
                                                  evolution.secondary_energies[1:10:end])
                @printf(io, "%.4f,%.10g,%.10g\n", time, primary, secondary)
            end
        end
        features = evolution_features(evolution)
        @printf("Mμ = %.2f: E1 peak %.1f at %.3f, E2 peak %.1f at %.3f, final E1 %.1f, E2 %.1f -> %s\n",
                alpha, features.primary_peak, features.primary_peak_time,
                features.secondary_peak, features.secondary_peak_time,
                features.primary_final, features.secondary_final, path)
    end
end
