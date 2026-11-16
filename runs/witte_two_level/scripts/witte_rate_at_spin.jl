# One relativistic quartic scattering rate from Witte & Mummery's code (arXiv 2412.03655,
# github.com/SamWitte/Axion_SR_TDE, commit bac01b20d2d14230ae4ea065aac1ac16329c0292) at a chosen (alpha, spin).
# Their src/Compute_all_rates.jl computes the same rate on an alpha grid with the spin fixed at a = 0.95 (its line 97);
# this driver repeats the call of its main() for a single alpha and any spin, with the same numerical settings.
# Their source is not edited.
#
# Usage:
#   julia --project=<witte env> witte_rate_at_spin.jl <witte src dir> <S1> <S2> <S3> <BH|Inf> <alpha> <spin> <NptsCh>
# NptsCh: their Chebyshev point count, which main() spaces linearly from 40 to 70 over the alpha grid (65 at the table
# point alpha = 0.30066 of the 19-point grid of 211_211_322_BH_LvrHc_.dat). Prints the rate in the units of their table.
using Printf
witte_src_dir, state_1, state_2, state_3, destination, alpha_text, spin_text, chebyshev_points_text = ARGS
cd(witte_src_dir)
include(joinpath(witte_src_dir, "solve_sr_rates.jl"))
include(joinpath(witte_src_dir, "Core/constants.jl"))
include(joinpath(witte_src_dir, "state_utils.jl"))

coupling = parse(Float64, alpha_text)
spin = parse(Float64, spin_text)
chebyshev_points = parse(Int, chebyshev_points_text)

# settings of main() in Compute_all_rates.jl
kpts, rpts, rmaxT, Nang, Npts_Bnd = 14, 1500, 100, 5000000, 2000
M = 1.0
Ntot_safe = 30000
NON_REL = false
prec = 200
cvg_acc = 1e-20
iterC = 10
Lcheb = 4
der_acc = 1e-20
run_leaver, use_heunc, BHlmax, ang_method = true, true, 0, :mc   # their command-line defaults

n1, l1, m1 = parse_state_from_args(state_1)
n2, l2, m2 = parse_state_from_args(state_2)
n3, l3, m3 = parse_state_from_args(state_3)
to_inf = destination != "BH"
mu = coupling / (M * GNew)
maxN = maximum([n1 n2 n3])
rmax_1 = (2.0 .^ (2.0 .* 3 .- 2 .* (1 .+ 3)) .* gamma(2 .+ 2 .* 3) ./ (0.03 .^ 2) ./ factorial(2 .* 3 - 1) .* 7.0)
rmax_ratio = (2.0 .^ (2.0 .* maxN .- 2 .* (1 .+ maxN)) .* gamma(2 .+ 2 .* maxN) ./ coupling .^ 2 ./ factorial(big(2 .* maxN - 1)) .* 7.0) ./ rmax_1
h_mve = 0.2 ./ rmax_ratio

start = time()
rate = gf_radial(mu, M, spin, n1, l1, m1, n2, l2, m2, n3, l3, m3; rpts=rpts, Npts_Bnd=Npts_Bnd, eps_fac=1e-3,
                 Ntot_safe=Ntot_safe, Nang=Nang, NON_REL=NON_REL, h_mve=h_mve, to_inf=to_inf, rmaxT=rmaxT,
                 run_leaver=run_leaver, NptsCh=chebyshev_points, cvg_acc=cvg_acc, prec=prec, iterC=iterC, debug=true,
                 der_acc=der_acc, Lcheb=Lcheb, use_heunc=use_heunc, BHlmax=BHlmax, self_grav=false, ang_method=ang_method)
@printf("RESULT %s_%s_%s_%s alpha = %.10g spin = %.6f NptsCh = %d rate = %.10e (%.0f s)\n",
        state_1, state_2, state_3, destination, coupling, spin, chebyshev_points, Float64(rate), time() - start)
