"""Tune the PML (gmax, eps_match, xmin, resolution) so Gamma -> Leaver at a=0.5."""
import numpy as np, evolve as ev
configs=[
  dict(N=4000,T=3000,xmin=-90,gmax=1.5,eps_match=0.005),
  dict(N=8000,T=3000,xmin=-90,gmax=1.5,eps_match=0.005),
  dict(N=4000,T=3000,xmin=-120,gmax=1.5,eps_match=0.002),
]
for c in configs:
    o=ev.run(alpha=0.5,a=0.9,xmax=220.,**c)
    E,G=ev.logfit(o["t"],o["psi"]); w=o["w_ref"]
    print(f"N={c['N']} xmin={c['xmin']} gmax={c['gmax']} eps={c['eps_match']}: "
          f"E {E/w.real:.4f}x  Gamma {G/w.imag:+.3f}x  (G={G:+.3e} vs {w.imag:+.3e})",flush=True)
