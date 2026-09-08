"""Grid convergence test (Dolan Sec. III.G, check 1: the ratio test).

Holds the deep-PML domain fixed (x_left=-200, bump at -100 width 10, x_right=250)
and refines the grid n over a factor-of-two sequence n = 4, 8, 16, i.e. dx =
M/4, M/8, M/16.  For a 4th-order-accurate scheme the error should fall by 2^4=16
per halving, so the Richardson ratio

    ratio = (f_coarse - f_medium) / (f_medium - f_fine)

should approach 16 for both the decay rate Gamma (time-domain slope, Dolan's
method) and the real frequency omega_R.  The Richardson extrapolation

    f_extrap = f_fine + (f_fine - f_medium) / (2^p - 1)     (p = 4)

removes the leading finite-dx error.  NOTE: omega_R converges cleanly (ratio
~16); the *fitted* Gamma does NOT extrapolate onto Leaver, because its residual
is not discretization -- it is PML-absorption contamination (the time-domain fit
depends on PML width; see diagnose_gamma.py).  The eigenvalue Gamma is the clean,
PML-independent quantity and matches spheroidal Leaver to ~0.01%.  The benchmark
is evaluated from leaver.py at runtime (no hardcoded constant).

Single channel (l_max=1): at a=0.9, alpha=0.5 the tower is negligible (verified
in measure_gamma), so we isolate the pure dx dependence here.
"""
import numpy as np
from measure_gamma import measure_gamma, leaver_benchmark

RICHARDSON_ORDER = 4                      # expected 4th-order convergence
resolutions = [4, 8, 16]                  # dx = M/n, factor-of-two refinement


def richardson_ratio(values):
    """(coarse - medium) / (medium - fine); -> 2^p for order p on halving dx."""
    coarse, medium, fine = values
    return (coarse - medium) / (medium - fine)


def richardson_extrapolate(medium, fine, order=RICHARDSON_ORDER):
    """Leading-error-removed estimate from the two finest grids."""
    return fine + (fine - medium) / (2**order - 1)


if __name__ == "__main__":
    w_leaver = None
    gamma_of_n = {}
    omega_R_of_n = {}
    for n in resolutions:
        tower, times, norms, gamma, omega = measure_gamma(l_max=1, n=n)
        if w_leaver is None:                      # params fixed across n; solve once
            w_leaver = leaver_benchmark(tower)
        gamma_of_n[n] = gamma
        omega_R_of_n[n] = omega.real
        print(f"n={n:3d}  dx={tower.dx:.4f}  "
              f"omega_R={omega.real:.6f} ({omega.real/w_leaver.real:.5f}x)  "
              f"Gamma={gamma:+.4e} ({gamma/w_leaver.imag:.4f}x)")

    gamma_values = [gamma_of_n[n] for n in resolutions]
    omega_values = [omega_R_of_n[n] for n in resolutions]

    print(f"\nRichardson ratio (expect ~{2**RICHARDSON_ORDER} for 4th order):")
    print(f"  omega_R:  {richardson_ratio(omega_values):.2f}")
    print(f"  Gamma:    {richardson_ratio(gamma_values):.2f}")

    omega_R_extrap = richardson_extrapolate(omega_values[1], omega_values[2])
    gamma_extrap = richardson_extrapolate(gamma_values[1], gamma_values[2])
    print(f"\nRichardson extrapolation (dx -> 0) vs Leaver:")
    print(f"  omega_R: {omega_R_extrap:.6f}  (Leaver {w_leaver.real:.6f}, "
          f"{omega_R_extrap/w_leaver.real:.5f}x)")
    print(f"  Gamma:   {gamma_extrap:+.4e}  (Leaver {w_leaver.imag:+.4e}, "
          f"{gamma_extrap/w_leaver.imag:.4f}x)")
