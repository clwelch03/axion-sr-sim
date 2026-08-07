"""Stage D: the actual deliverable -- superradiant growth peak, a=0.99,
alpha=0.42 (Leaver: E~0.4088, Gamma~+1.52e-7, e-fold ~6.6e6 M). logfit is
infeasible here (would need ~1e7 M to separate growth from transient by
passive decay), so this leans on the now-fixed anchor-selecting prony().
Near superradiant threshold E~m*OmH, so the near-horizon wavelength
lambda=2pi/|E-m*OmH| is large (~251 here, vs ~38 for the damped a=0.9
benchmark) -- PML needs to be wider in absolute x than before.
"""
import numpy as np, sys
sys.path.insert(0, "/global/scratch/projects/pc_heptheory/clwelch/software/GRChombo-Axion-SR/Examples/AxionSuperradiance/linear_rates")
import leaver
import evolve as ev

alpha, a = 0.42, 0.99
w_leaver = complex(leaver.solve(alpha, a))
print(f"Leaver target: E={w_leaver.real:.6f}  Gamma={w_leaver.imag:+.4e}", flush=True)

o = ev.run(alpha=alpha, a=a, N=7700, xmin=-670., xmax=250., T=10000.,
           eps_match=0.002, gmax=1.5, nrec=8000)
print(f"probe r={o['probe_r']:.2f}  dx={o['dx']:.4f}  x_match={o['x_match']:.2f}", flush=True)

t, psi = o["t"], o["psi"]
E1, G1 = ev.logfit(t, psi)
print(f"logfit : E={E1:.6f}  Gamma={G1:+.4e}  (Gamma {G1/w_leaver.imag:+.3f}x, expect garbage -- too short for logfit)", flush=True)

for tmin_frac in [0.1, 0.35, 0.6]:
    try:
        E2, G2, om_all, c_all = ev.prony(t, psi, p=12, tmin_frac=tmin_frac, anchor=w_leaver, full=True)
        print(f"\nprony (tmin_frac={tmin_frac}): E={E2:.6f} Gamma={G2:+.4e}  (Gamma {G2/w_leaver.imag:+.3f}x)", flush=True)
        order = np.argsort(np.abs(om_all - w_leaver))
        print("  candidates nearest Leaver anchor (E, Gamma, |c|, dist):", flush=True)
        for k in order[:6]:
            print(f"    E={om_all[k].real:+.6f}  Gamma={om_all[k].imag:+.4e}  |c|={abs(c_all[k]):.3e}"
                  f"  dist={abs(om_all[k]-w_leaver):.4e}", flush=True)
    except Exception as e:
        print(f"prony (tmin_frac={tmin_frac}) failed: {e}", flush=True)

np.savez("stage_d_growth.npz", t=t, psi=psi, w_leaver=w_leaver,
         phi_flux=o["phi_flux"], phidot_flux=o["phidot_flux"], phix_flux=o["phix_flux"],
         flux_r=o["flux_r"])
print("saved stage_d_growth.npz", flush=True)
