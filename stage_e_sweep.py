"""Stage E: Gamma(alpha) sweep at a=0.9 (damped regime) vs Leaver, using the
validated PML config (eps_match=0.002, xmin=-120, N=4000, gmax=1.5) and the
fixed anchor-selecting prony() for extraction. Overlays the 3D DMD point
(a=0.5 damped, 0.92x Leaver, from the GRChombo-3D kerr_schild_reslice work)
to show where the 1+1D tool reaches that 3D could not.
"""
import numpy as np, sys, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, "/global/scratch/projects/pc_heptheory/clwelch/software/GRChombo-Axion-SR/Examples/AxionSuperradiance/linear_rates")
import leaver
import evolve as ev

a = 0.90
alphas = [0.30, 0.40, 0.50, 0.60]   # damped regime at a=0.9 (E > m*OmH there)

results = []
for alpha in alphas:
    w = complex(leaver.solve(alpha, a))
    o = ev.run(alpha=alpha, a=a, N=4000, xmin=-120., xmax=220., T=3000.,
               eps_match=0.002, gmax=1.5)
    t, psi = o["t"], o["psi"]
    E1, G1 = ev.logfit(t, psi)
    try:
        E2, G2 = ev.prony(t, psi, p=10, tmin_frac=0.35, anchor=w)
    except Exception as e:
        E2, G2 = float("nan"), float("nan")
        print(f"alpha={alpha}: prony failed: {e}", flush=True)
    E3, G3 = ev.freqfilter(t, psi, E_peak=w.real, kbins=5.0)   # validated best
    results.append((alpha, w.real, w.imag, E1, G1, E2, G2, E3, G3))
    print(f"alpha={alpha}: Leaver Gamma={w.imag:+.4e}  logfit {G1/w.imag:+.3f}x"
          f"  prony {G2/w.imag:+.3f}x  freqfilter {G3/w.imag:+.3f}x", flush=True)

results = np.array(results)
np.savez("stage_e_sweep.npz", results=results,
         cols="alpha,E_leaver,Gamma_leaver,E_logfit,Gamma_logfit,E_prony,Gamma_prony,E_freqfilter,Gamma_freqfilter")

fig, ax = plt.subplots(figsize=(7, 5))
al, Gl = results[:, 0], results[:, 2]
Gf, Gp, Gff = results[:, 4], results[:, 6], results[:, 8]
ax.plot(al, np.abs(Gl), "k--", marker="s", label="Leaver (truth)")
ax.plot(al, np.abs(Gf), "C0", marker="o", label="logfit")
ax.plot(al, np.abs(Gp), "C1", marker="^", label="prony (anchored)")
ax.plot(al, np.abs(Gff), "C2", marker="D", label="freqfilter (validated)")
ax.scatter([0.50], [5.2050e-04 * 0.92], color="red", marker="*", s=150,
           label="3D DMD (a=0.9,alpha=0.5 damped, 0.92x)", zorder=5)
ax.set_yscale("log")
ax.set_xlabel("alpha"); ax.set_ylabel("|Gamma|")
ax.set_title(f"Gamma(alpha) at a={a}: 1+1D vs Leaver vs 3D DMD")
ax.legend()
plt.tight_layout()
plt.savefig("stage_e_sweep.png", dpi=130)
print("saved stage_e_sweep.png", flush=True)
