"""Test the fixed prony() against the a=0.9 alpha=0.5 production config
(the same beating data floor_check.py examined). Prints the full candidate
list so we can verify Prony genuinely finds a mode near Leaver, not just
defaulting to the anchor with nothing nearby."""
import numpy as np
import evolve as ev

o = ev.run(alpha=0.50, a=0.90, N=4000, xmin=-120., xmax=220., T=3000.,
           eps_match=0.002, gmax=1.5)
w = o["w_ref"]; t = o["t"]; psi = o["psi"]

E1, G1 = ev.logfit(t, psi)
print(f"Leaver : E={w.real:.6f}  Gamma={w.imag:+.4e}")
print(f"logfit : E={E1:.6f}  Gamma={G1:+.4e}   ({G1/w.imag:+.3f}x)")

for tmin_frac in [0.1, 0.35, 0.6]:
    E2, G2, om_all, c_all = ev.prony(t, psi, p=10, tmin_frac=tmin_frac, anchor=w, full=True)
    print(f"\nprony (tmin_frac={tmin_frac}): E={E2:.6f} Gamma={G2:+.4e} ({G2/w.imag:+.3f}x)")
    order = np.argsort(np.abs(om_all - w))
    print("  candidates nearest Leaver anchor (E, Gamma, |c|):")
    for k in order[:6]:
        print(f"    E={om_all[k].real:+.6f}  Gamma={om_all[k].imag:+.4e}  |c|={abs(c_all[k]):.3e}"
              f"  dist_to_Leaver={abs(om_all[k]-w):.4e}")
