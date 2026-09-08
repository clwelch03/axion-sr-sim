"""Eigenvalue-only convergence of Gamma vs outer boundary x_right and resolution n.

Ran to test whether the eigenvalue Gamma's offset from the (flat-Lambda) Leaver
value was outer-boundary truncation (we use x_right=250; Dolan uses ~500). Result:
Im(omega) is FLAT across x_right and n -- not the boundary, not dx. The offset was
the benchmark, not the solver: the flat Lambda=l(l+1) Leaver value omits the
spheroidal angular correction that our operator carries. Against the full spheroidal
Leaver value (leaver.py, spheroidal=True, evaluated at runtime below) the eigenvalue
matches to ~0.01%. Kept as the record of that check; pure sparse eigen-solves, no
evolution, so it is fast.

a=0.9, alpha=0.5, l=m=1, deep PML (bump at -100, width 10).
"""
import numpy as np
from solver import Tower
from qbs_eigenmode import solve_qbs
from measure_gamma import leaver_benchmark


def eigen_gamma(x_right, n):
    tower = Tower(spin=0.9, alpha=0.5, m=1, l_max=1,
                  x_left=-200.0, x_right=x_right, n=n,
                  x_pml_center=-100.0, pml_width=10.0)
    omega, _ = solve_qbs(tower)
    return omega


if __name__ == "__main__":
    # runtime benchmark: full spheroidal Leaver at this case's parameters (no constant)
    ref_tower = Tower(spin=0.9, alpha=0.5, m=1, l_max=1, x_left=-200.0,
                      x_right=250.0, n=8, x_pml_center=-100.0, pml_width=10.0)
    w_leaver = leaver_benchmark(ref_tower)
    print(f"Leaver (spheroidal, runtime): omega_R={w_leaver.real:.6f}, "
          f"Gamma={w_leaver.imag:+.4e}\n")

    for n in (8, 16):
        for x_right in (250.0, 400.0, 600.0):
            omega = eigen_gamma(x_right, n)
            print(f"n={n:3d}  x_right={x_right:6.0f}  "
                  f"omega_R={omega.real:.6f} ({omega.real/w_leaver.real:.5f}x)  "
                  f"Gamma={omega.imag:+.4e} ({omega.imag/w_leaver.imag:.4f}x)")
        print()
