# Driver: Witte & Mummery's solve_system (arXiv 2412.03655, github.com/SamWitte/Axion_SR_TDE, commit
# bac01b20d2d14230ae4ea065aac1ac16329c0292) restricted to the |211> and |322> levels, started from our seeds.
# Modelled on src/Nlevels_Runs/run_Nlevels.jl of their repository; the fixed numerical parameters are theirs.
#
# Usage:
#   julia --project=<witte env> run_witte_two_level.jl <witte src dir> <f_a/M_Pl> <end time in M0> <channel list file> <output csv>
# Output csv columns (their units): time_years, u_211, u_322, spin, mass_Msun, with u = N/(G M_BH^2) (M_BH the
# initial mass, 10 Msun here).
using DelimitedFiles, Printf
witte_src_dir, ratio_text, end_time_text, channel_list, output_csv = ARGS
cd(witte_src_dir)
# Julia 1.13 no longer makes `OrdinaryDiffEq.LinearSolve` resolvable from `using LinearSolve` inside OrdinaryDiffEq, but
# solve_system_unified.jl (radau branch) writes OrdinaryDiffEq.LinearSolve.LUFactorization(). Bind that name from outside,
# to the same LinearSolve module; their source is not edited and the solver is the same.
using OrdinaryDiffEq
let linear_solve_module = Base.loaded_modules[Base.PkgId(Base.UUID("7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"), "LinearSolve")]
    isdefined(OrdinaryDiffEq, :LinearSolve) || @eval OrdinaryDiffEq const LinearSolve = $linear_solve_module
end
include(joinpath(witte_src_dir, "super_rad.jl"))

planck_ratio = parse(Float64, ratio_text)           # f_a / M_Pl, our label
end_time_in_initial_mass = parse(Float64, end_time_text)
alpha_initial, spin_initial, mass_bh = 0.3, 0.9, 10.0
decay_constant_gev = planck_ratio * M_pl            # their M_pl = 1.22e19 GeV
mu_ev = alpha_initial / (GNew * mass_bh)            # driver's formula
kappa = planck_ratio^2
seed_in_our_units = 3e-6                            # N_211 = N_322 = 3e-6 f_a^2 M0^2
initial_occupation = kappa * seed_in_our_units                    # u = kappa * N (initial mass is the normalisation mass)

# their time unit: years of 3.15e7 s (Core/constants.jl, YEAR_IN_SECONDS), M0 = G M_BH hbar seconds
seconds_per_initial_mass = GNew * mass_bh * hbar
years_per_initial_mass = seconds_per_initial_mass / 3.15e7
end_time_years = end_time_in_initial_mass * years_per_initial_mass
println("f_a = ", decay_constant_gev, " GeV, kappa = ", kappa, ", initial_occupation u = ", initial_occupation, ", years per M0 = ", years_per_initial_mass,
        ", end time = ", end_time_years, " yr, bosenova cap u(211) = ", 1024 * pi * kappa * (2 / 2)^4 / (9 * alpha_initial^3))

keep_two_levels = (n, l, m) -> (n, l, m) in ((2, 1, 1), (3, 2, 2))
initial_log_state = log.([initial_occupation, initial_occupation, spin_initial, mass_bh])

time_years, states, modes, spin, mass = solve_system(
    mu_ev, decay_constant_gev, spin_initial, mass_bh, end_time_years;
    impose_low_cut = 1e-100, return_all_info = true, n_times = 1000000, eq_threshold = 1e-100, abstol = 1e-30,
    non_rel = false, debug = true, high_p = true, N_pts_interp = 100, N_pts_interpL = 100, Nmax = 3, cheby = true,
    lm_only = false, u0_override = initial_log_state, t_start = 0.0, gw_model = :off, f_edd = 0.0, track_alpha = true,
    save_factor = 5.0, save_dlna = 1e-3, self_grav = false, rate_input = channel_list, track_alpha_refresh = 0.0,
    solver = :radau, analytic_jac = false, keep_levels = keep_two_levels, maxiters = 5e6, timescale_saves = true,
    regime_split = false, t_save_start = 1e-9, emax2_mode = nothing)

println("modes evolved: ", modes)
println("final t = ", time_years[end], " of ", end_time_years, " yr; spin = ", spin[end], "  mass = ", mass[end])
open(output_csv, "w") do io
    println(io, "time_years,u_211,u_322,spin,mass_Msun")
    for i in eachindex(time_years)
        @printf(io, "%.15e,%.15e,%.15e,%.15e,%.15e\n", time_years[i], states[1][i], states[2][i], spin[i], mass[i])
    end
end
println("wrote ", output_csv, " (", length(time_years), " rows)")
if time_years[end] < end_time_years * (1 - 1e-12)
    println("INCOMPLETE: stopped before the end time")
end
