#!/usr/bin/env python
"""Coupled (l, m=1) angular tower for the 1+1D nonlinear axion QBS (README Sec 10.5 [B]).

Generalizes evolve.py from the single l=m=1 mode to the odd-l tower
    Phi = (1/r) sum_{l in [1,3,5,...]} psi_l(t,x) Y_{l1}(theta) e^{i phi}.
The RWA quartic (3/4)|Psi|^2 Psi is NOT diagonal in l: cubing an l=1 cloud
generates l=3 at the SAME frequency (angular overlap Q^3_{111} != 0), l=3 feeds
l=5, and so on. THIS multi-mode angular cascade -- not single-mode detuning --
is the physically dominant saturation channel (README 10.1/10.6). This module
lets us watch it and answer "how many angular modes does saturation need?" by
converging in l_max.

What is and isn't coupled here:
  * Each channel l carries its OWN linear radial operator: the single-mode
    operator of evolve.py with that l (l(l+1) and the diagonal metric weight
    c_ll = <Y_l1|cos^2|Y_l1>). Higher l are deeper-bound and grow much more
    slowly on their own; they are populated mainly by the cubic from l=1.
  * The NONLINEAR quartic couples all channels via angular.coupling_tensors.
    This is the new physics and the deliverable.
  * The OFF-DIAGONAL LINEAR (cos^2 theta) metric coupling l -> l+-2 is O(c^2) tiny
    for a QBS (README Sec 2, ~1.6e-3) and would couple psi_tt across channels
    (a channel mass-matrix); it is DROPPED here, matching the README's
    single-mode-truncation stance. `angular.cos2_coupling` provides the
    off-diagonals if we ever need to restore them.

The horizon PML (Dolan near-horizon reduced form) is l-INDEPENDENT -- deep in the
tortoise region V->0 for every l and the reduced equation
psi_tt = psi_xx + 2im Om(psi_x - psi_t) has no l -- so a single shared absorber
serves all channels.

REGRESSION GUARANTEE: with l_max=1 (one channel) every array here is identical to
evolve.run(...) for the same arguments, because Qbar[0,0,0,0]=1 and the operator
is the same single-mode operator. tower_test.py checks this to machine precision.
"""
import numpy as np, sys, argparse, os

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))           # radial1d/  (radial1d.py, evolve.py)
sys.path.insert(0, "/global/scratch/projects/pc_heptheory/clwelch/software/GRChombo-Axion-SR/Examples/AxionSuperradiance/linear_rates")
import leaver
from radial1d import build_grid, coeffs, d_dx, d2_dx2, horizons, tortoise_x
import angular


def run(alpha=0.50, a=0.90, l_max=3, xmin=-70., xmax=220., N=4000, T=3000., cfl=0.4,
        x0=20., sigx=6., eps_match=0.02, gmax=1.0, probe_r=None, nrec=6000,
        lambda_quartic=0.0, init_amplitude=1.0, seed_channels=None, seed="gaussian"):
    """Evolve the coupled (l, m=1) tower up to l_max.

    l_max          : highest angular index; channels are [1,3,...,l_max] (odd sector).
    lambda_quartic : the SAME signed knob as evolve.run -- normalized so the
                     one-channel tower matches single-mode exactly. The channel-j
                     interior momentum RHS gains the Hermitian cubic
                        -Delta * lambda_quartic * sum_{abc} Qbar^j_{abc} psi_a psi_b* psi_c.
    init_amplitude : Gaussian seed amplitude placed on the l=1 channel.
    seed_channels  : optional dict {l: amplitude} to seed higher channels too
                     (default: only l=1 seeded at init_amplitude; the tower is
                     meant to GROW the higher channels from the cubic).

    Returns a dict with per-channel time series (probe psi_l(t), interior L2 norm
    cloud_norm_l(t)) plus the shared metadata, so one can watch l=3 rise out of l=1.
    """
    w_ref = complex(leaver.solve(alpha, a))
    rp, _ = horizons(a); OmH = a / (2 * rp); m = 1
    x, r = build_grid(a, xmin, xmax, N); dx = x[1] - x[0]
    ls, _c_matrix, Qbar = angular.coupling_tensors(l_max, m)
    nl = len(ls)

    # per-channel linear coefficients (only D and Vc depend on l; A2,A1,At are m-only)
    D = np.empty((nl, N)); Vc = np.empty((nl, N), complex)
    for i, l in enumerate(ls):
        c_ll = angular.cos2_coupling(l, l, m)
        Di, A2, A1, At, Vci, Delta, P = coeffs(r, a, alpha, l=l, m=m, c_ll=c_ll)
        D[i] = Di; Vc[i] = Vci
    invD = 1.0 / D                                    # A2,A1,At,Delta,P are shared (l-independent)

    # PML region and attenuation (shared across channels -- reduced eqn has no l)
    x_match = tortoise_x(rp + eps_match, a)
    pml = x < x_match
    g = np.zeros_like(x)
    g[pml] = gmax * ((x_match - x[pml]) / (x_match - xmin)) ** 2
    xo = xmax - 40; out = x > xo
    gout = np.zeros_like(x); gout[out] = 0.8 * ((x[out] - xo) / 40) ** 2
    imO = 1j * m * OmH
    interior = ~pml
    phys = interior & ~out
    lam_interior = (interior * Delta * lambda_quartic)   # shared radial prefactor of the cubic

    if probe_r is None: probe_r = min(2.0 / alpha**2, 0.6 * r[-1])
    ip = int(np.argmin(np.abs(r - probe_r)))
    dt = cfl * dx; nsteps = int(T / dt)

    def rhs(psi, pi, chi):
        # psi,pi,chi are (nl, N) complex. Derivatives per channel (nl is small).
        dpsi = np.empty_like(psi); dpi = np.empty_like(pi); dchi = np.empty_like(chi)
        # nonlinear cubic force per channel (interior): F_j = sum_abc Qbar[j,a,b,c] psi_a psi_b* psi_c
        if lambda_quartic != 0.0:
            force = np.einsum('jabc,ax,bx,cx->jx', Qbar, psi, np.conj(psi), psi)
        for i in range(nl):
            px = d_dx(psi[i], dx); pxx = d2_dx2(psi[i], dx); cx = d_dx(chi[i], dx)
            dpi_int = invD[i] * (A2 * pxx + A1 * px + At * pi[i] - Vc[i] * psi[i])
            if lambda_quartic != 0.0:
                dpi_int = dpi_int - invD[i] * lam_interior * force[i]
            dpi_pml = pxx + 2 * imO * (px - pi[i]) - cx - imO * chi[i] - g * (pi[i] + imO * psi[i])
            dpi[i] = np.where(pml, dpi_pml, dpi_int) - gout * pi[i]
            dpsi[i] = pi[i] - gout * psi[i]
            dchi[i] = np.where(pml, g * (px - chi[i] + imO * psi[i]) - imO * chi[i], 0.0)
        return dpsi, dpi, dchi

    # seed: l=1 (Gaussian or clean eigenmode) plus any requested higher channels;
    # higher channels are meant to GROW from the cubic, so default them to zero.
    seed_amps = dict(seed_channels) if seed_channels else {}
    seed_amps.setdefault(ls[0], init_amplitude)
    psi = np.zeros((nl, N), complex)
    pi = np.zeros_like(psi); chi = np.zeros_like(psi)
    if seed == "eigenmode":
        from radial1d import eigenmode_seed
        # clean l=1 eigenmode (psi and pi=-i w psi) on channel 0; the l=1
        # eigenfrequency w_ref matches the ls[0]=1 channel operator exactly.
        psi[0], pi[0] = eigenmode_seed(x, r, alpha, a, w_ref, l=ls[0], m=m,
                                       amplitude=seed_amps[ls[0]])
        higher = [(i, l) for i, l in enumerate(ls) if i != 0 and l in seed_amps]
    else:
        higher = list(enumerate(ls))
    for i, l in higher:
        if l in seed_amps:
            psi[i] = seed_amps[l] * np.exp(-((x - x0) ** 2) / (2 * sigx**2))

    trec = np.linspace(0, T, nrec); rec_i = (trec / dt).astype(int); rp_ = 0
    tsig = []; psig = []; cloud = []       # psig[k] = psi[:,ip] at record k ; cloud[k] = per-channel L2
    for n in range(nsteps):
        if rp_ < nrec and n >= rec_i[rp_]:
            tsig.append(n * dt); psig.append(psi[:, ip].copy())
            cloud.append(np.sqrt(np.sum(np.abs(psi[:, phys]) ** 2, axis=1) * dx))
            rp_ += 1
        k1 = rhs(psi, pi, chi)
        k2 = rhs(psi + 0.5 * dt * k1[0], pi + 0.5 * dt * k1[1], chi + 0.5 * dt * k1[2])
        k3 = rhs(psi + 0.5 * dt * k2[0], pi + 0.5 * dt * k2[1], chi + 0.5 * dt * k2[2])
        k4 = rhs(psi + dt * k3[0], pi + dt * k3[1], chi + dt * k3[2])
        psi = psi + (dt / 6) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        pi = pi + (dt / 6) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        chi = chi + (dt / 6) * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2])
        psi[:, 0] = 0; pi[:, 0] = 0; psi[:, -1] = 0; pi[:, -1] = 0
        if not np.all(np.isfinite(psi)): print(f"  NaN at t={n*dt:.1f}"); break

    return dict(w_ref=w_ref, OmH=OmH, ls=np.array(ls), l_max=l_max,
                t=np.array(tsig), psi=np.array(psig).T,            # psi[channel, time]
                cloud_norm=np.array(cloud).T,                      # cloud_norm[channel, time]
                probe_r=r[ip], dx=dx, dt=dt, x_match=x_match,
                lambda_quartic=lambda_quartic, init_amplitude=init_amplitude)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--alpha", type=float, default=0.50); ap.add_argument("--a", type=float, default=0.90)
    ap.add_argument("--l_max", type=int, default=3, help="highest angular index (channels 1,3,...,l_max)")
    ap.add_argument("--T", type=float, default=3000.); ap.add_argument("--N", type=int, default=4000)
    ap.add_argument("--xmin", type=float, default=-70.); ap.add_argument("--xmax", type=float, default=220.)
    ap.add_argument("--gmax", type=float, default=1.0); ap.add_argument("--eps_match", type=float, default=0.02)
    ap.add_argument("--lambda", type=float, default=0.0, dest="lambda_quartic",
                    help="signed quartic self-coupling (<0 attractive/axion, >0 repulsive)")
    ap.add_argument("--amp", type=float, default=1.0, dest="init_amplitude", help="l=1 seed amplitude")
    ap.add_argument("--save", type=str, default=None)
    q = ap.parse_args()
    import evolve  # for the shared inst_gamma extractor
    o = run(q.alpha, q.a, q.l_max, q.xmin, q.xmax, q.N, q.T, gmax=q.gmax,
            eps_match=q.eps_match, lambda_quartic=q.lambda_quartic, init_amplitude=q.init_amplitude)
    ls = o["ls"]; w = o["w_ref"]
    print(f"alpha={q.alpha} a={q.a} l_max={q.l_max} channels={list(ls)} N={q.N} T={q.T} "
          f"lambda={q.lambda_quartic:+.3g} amp={q.init_amplitude:g}  (probe r={o['probe_r']:.1f})")
    print(f"  Leaver (l=1): E={w.real:.6f}  Gamma={w.imag:+.4e}")
    # final channel populations (fraction of total interior norm) and per-channel Gamma
    final_norm = o["cloud_norm"][:, -1]
    frac = final_norm**2 / np.sum(final_norm**2)
    for i, l in enumerate(ls):
        t_g, g_t = evolve.inst_gamma(o["t"], np.maximum(o["cloud_norm"][i], 1e-300))
        late = t_g > 0.6 * t_g[-1]
        print(f"  l={l:2d}: ||psi||_final={final_norm[i]:.4e}  norm-frac={frac[i]:.3e}  "
              f"Gamma(t) late-mean={np.mean(g_t[late]):+.4e}")
    if q.save is not None:
        np.savez(q.save, **o); print("saved", q.save)


if __name__ == "__main__": main()
