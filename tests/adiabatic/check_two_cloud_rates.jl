# julia --project=. tests/adiabatic/check_two_cloud_rates.jl
# Step 6b checks of two_cloud_rates (piece 6, docs/multimode-adiabatic-design.md):
# exact limits of the energy and angular-momentum balance, one channel at a time,
# at saved points of both branches.
#   1. Primary alone, E2 = 0: dE1/dt = -F1, and F1 = -2 Im ω1 E1 to the branch's
#      own balance accuracy (the initial straight line of Fig. 6).
#   2. Secondary alone, any E2: dE1/dt = 0 exactly, dE2/dt = -F2.
#   3. m = 0 alone: F/ω0 quanta per unit time, each taking two primary quanta
#      and giving one secondary quantum: dE1/dt = -(2ω1/ω0) F, dE2/dt = (ω2/ω0) F.
#   4. m = 3 alone: F/ω3 quanta per unit time, each taking two secondary quanta
#      and giving one primary quantum: dE1/dt = (ω1/ω3) F, dE2/dt = -(2ω2/ω3) F.
# 3 and 4 neglect the dω2/dE1 term, which is O(E2); at E2 = 1e-8 they should
# hold to about that size relative.
using Printf
include(joinpath(@__DIR__, "../../src/adiabatic/drivers/two_cloud_evolution.jl"))

for alpha in (0.42, 0.3)
    table = read_two_cloud_table(alpha)
    coefficients = two_cloud_coefficients(table)
    @printf("Mμ = %.2f\n", alpha)
    for index in (1, length(table.primary_energy) ÷ 3, length(table.primary_energy))
        primary_energy = table.primary_energy[index]
        omega1 = table.primary_omega[index]
        omega2 = table.secondary_omega[index]
        omega0, omega3 = 2omega1 - omega2, 2omega2 - omega1
        @printf("  E1 = %.4g\n", primary_energy)

        rates = two_cloud_rates(coefficients, primary_energy, 0.0; channels=(:primary,))
        @printf("    1. primary, E2 = 0:  |dE1/dt + F_E| / F_E = %.1e, dE2/dt = %.1e\n",
                abs(rates.primary_rate + rates.energy_flux) / abs(rates.energy_flux),
                rates.secondary_rate)

        rates = two_cloud_rates(coefficients, primary_energy, 50.0; channels=(:secondary,))
        @printf("    2. secondary, E2 = 50:  dE1/dt / F_E = %.1e, |dE2/dt + F_E| / F_E = %.1e\n",
                rates.primary_rate / rates.energy_flux,
                abs(rates.secondary_rate + rates.energy_flux) / abs(rates.energy_flux))

        rates = two_cloud_rates(coefficients, primary_energy, 1e-8; channels=(:m0,))
        quanta = rates.energy_flux / omega0
        @printf("    3. m = 0, E2 = 1e-8:  relative errors %.1e (dE1/dt), %.1e (dE2/dt)\n",
                abs(rates.primary_rate / (-2omega1 * quanta) - 1),
                abs(rates.secondary_rate / (omega2 * quanta) - 1))

        rates = two_cloud_rates(coefficients, primary_energy, 1e-8; channels=(:m3,))
        quanta = rates.energy_flux / omega3
        @printf("    4. m = 3, E2 = 1e-8:  relative errors %.1e (dE1/dt), %.1e (dE2/dt)\n",
                abs(rates.primary_rate / (omega1 * quanta) - 1),
                abs(rates.secondary_rate / (-2omega2 * quanta) - 1))
    end
    # Item 1, second half: F1/E1 at the first point against -2 Im ω1 (the time unit).
    @printf("  F1/E1 / (-2 Im ω1) at the first point: %.8f\n",
            table.primary_flux_per_energy[1] / (-2table.linear_growth_rate))
end
