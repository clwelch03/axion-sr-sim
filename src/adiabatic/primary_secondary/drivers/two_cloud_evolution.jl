# Piece 6 (docs/multimode-adiabatic-design.md): evolution of the primary and
# secondary energies by energy and angular-momentum balance (Omiya et al.
# 2211.01949 Eqs. (62)-(63)), for the paper's Fig. 6.
using Printf
include(joinpath(@__DIR__, "../load.jl"))

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    # julia --project=. src/adiabatic/primary_secondary/drivers/two_cloud_evolution.jl
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
