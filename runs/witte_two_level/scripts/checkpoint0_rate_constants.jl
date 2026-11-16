# Checkpoint 0: compare Witte & Mummery's load_rate_coeffs for the two-channel list with our CSV
# (data/witte_mummery_quartic_rates.csv).  Usage (repo root):
#   julia --project=<witte env> runs/witte_two_level/scripts/checkpoint0_rate_constants.jl <witte repo src dir>
using DelimitedFiles, Interpolations
witte_src_dir = ARGS[1]
repository_root = abspath(joinpath(@__DIR__, "..", "..", ".."))
cd(witte_src_dir)
include(joinpath(witte_src_dir, "super_rad.jl"))

two_channel_list = abspath(joinpath(@__DIR__, "two_channel_list.txt"))
alpha_initial, spin_initial = 0.3, 0.9
mass_bh = 10.0
mu = alpha_initial / (GNew * mass_bh)
coupling_ratio_to_planck = 3e-3                       # f_a / M_Pl of case A
decay_constant = coupling_ratio_to_planck * M_pl      # GeV, with their M_pl = 1.22e19
rates = load_rate_coeffs(mu, mass_bh, spin_initial, decay_constant, 3, Float64[]; non_rel=false,
                         rate_input=two_channel_list, regime_split=false)
println("their keys: ", collect(keys(rates)))

csv_lines = [ln for ln in eachline(joinpath(repository_root, "data", "witte_mummery_quartic_rates.csv"))
             if !startswith(ln, "#") && !startswith(ln, "alpha")]
table = reduce(vcat, [parse.(Float64, split(ln, ","))' for ln in csv_lines])
log_interpolate(column) = 10^LinearInterpolation(log10.(table[:, 1]), log10.(table[:, column]))(log10(alpha_initial))
planck_factor = (M_pl / decay_constant)^4
horizon_ratio = (1 + sqrt(1 - spin_initial^2)) / (1 + sqrt(1 - 0.95^2))
expected = Dict("211_211^322^BH" => log_interpolate(2) * planck_factor * horizon_ratio,
                "322_322^211^Inf" => log_interpolate(3) * planck_factor)
worst = 0.0
for (key, value) in expected
    theirs = rates[key]
    relative = theirs / value - 1
    global worst = max(worst, abs(relative))
    println(key, ": their load_rate_coeffs = ", theirs, "  our CSV (log-log interp, x(Mpl/fa)^4 [x r_P ratio]) = ", value, "  relative difference = ", relative)
end
println("alpha used = ", mu * GNew * mass_bh, "  f_a = ", decay_constant, " GeV  largest |relative difference| = ", worst)
