"""Validate freqfilter() on the a=0.9, alpha=0.5 DAMPED benchmark where the
answer is known (Leaver Gamma=-5.205e-4, prony gives 0.98-1.00x). If the
frequency-filter also recovers ~1.0x here, we trust it for Stage D (a=0.99
growth peak) where logfit beats and prony is Gamma-resolution-limited.
"""
import numpy as np
import evolve as ev

o = ev.run(alpha=0.50, a=0.90, N=4000, xmin=-120., xmax=220., T=3000.,
           eps_match=0.002, gmax=1.5)
w = o["w_ref"]; t, psi = o["t"], o["psi"]
print(f"Leaver     : E={w.real:.6f}  Gamma={w.imag:+.4e}")
E1, G1 = ev.logfit(t, psi)
print(f"logfit     : E={E1:.6f}  Gamma={G1:+.4e}  ({G1/w.imag:+.3f}x)")
E2, G2 = ev.prony(t, psi, anchor=w)
print(f"prony      : E={E2:.6f}  Gamma={G2:+.4e}  ({G2/w.imag:+.3f}x)")
for kb in [3.0, 5.0, 8.0]:
    E3, G3 = ev.freqfilter(t, psi, E_peak=w.real, kbins=kb)
    print(f"freqfilter (kbins={kb}): E={E3:.6f}  Gamma={G3:+.4e}  ({G3/w.imag:+.3f}x)")
