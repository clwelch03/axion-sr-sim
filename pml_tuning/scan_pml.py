"""Tune the PML (gmax, eps_match, xmin, resolution) so Gamma -> Leaver at a=0.5."""
import numpy as np, evolve as ev

# Dolan requires g(x) to vary smoothly over "a few wavelengths" of the local
# near-horizon advective wave (lambda = 2pi/|E - m*OmH| ~ 38 at a=0.9,
# alpha=0.5). Prior configs were only 1.5-2.65 lambda wide -- too abrupt.
# Target 4-5 lambda, holding dx~0.073 fixed (already shown resolution-
# converged) by scaling N with the larger domain.
configs=[
  dict(N=5600,T=3000,xmin=-190,gmax=1.5,eps_match=0.001),  # width~163 = 4.3 lambda
  dict(N=6000,T=3000,xmin=-220,gmax=1.5,eps_match=0.001),  # width~190 = 5.0 lambda
]
for c in configs:
    o=ev.run(alpha=0.5,a=0.9,xmax=220.,**c)
    E,G=ev.logfit(o["t"],o["psi"]); w=o["w_ref"]
    print(f"N={c['N']} xmin={c['xmin']} gmax={c['gmax']} eps={c['eps_match']}: "
          f"E {E/w.real:.4f}x  Gamma {G/w.imag:+.3f}x  (G={G:+.3e} vs {w.imag:+.3e})",flush=True)
