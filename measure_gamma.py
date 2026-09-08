"""Measure the QBS decay rate Gamma by TIME-DOMAIN evolution (Dolan's method).

Steps:
  1. build the Tower,
  2. get a transient-free eigenmode seed from qbs_eigenmode.solve_qbs,
  3. evolve it with Tower.rk4_step (Dolan Sec. III),
  4. fit the slope of log(interior norm) vs t  ->  Gamma.

Because the seed is a (near-)eigenmode, the interior norm decays as a clean
straight line in log from t=0 -- no shedding transient to wait out.

Benchmark: we DON'T hardcode a Leaver number -- we evaluate leaver.py at the
tower's own (alpha, spin, l, m) at runtime and compare.  Use spheroidal=True:
our operator carries the full a^2 mu^2 cos^2(theta) angular coupling, so the
apples-to-apples Leaver value is the spheroidal one, not the flat Lambda=l(l+1)
approximation (the flat one biases Gamma ~1% -- that gap is the missing angular
correction, not the solver).
"""
import numpy as np
from solver import Tower
from qbs_eigenmode import solve_qbs
import leaver


def leaver_benchmark(tower, n_overtone=2, spheroidal=True):
    """Full-spheroidal Leaver QBS frequency for this tower's parameters.

    Evaluated at runtime, so there is no stale constant to drift.  n_overtone is
    leaver.py's principal number (2 = nodeless fundamental = Fig.-13 'n=0').  The
    base angular channel is l=|m|.  Returns complex omega = omega_R + i Gamma.
    """
    l = tower.l_values[0]
    return leaver.solve(tower.alpha, tower.spin, l=l, m=tower.m,
                        n=n_overtone, spheroidal=spheroidal)


def interior_norm(psi, physical):
    """L2 norm of psi over the physical region (excludes the PML strip)."""
    return np.sqrt(np.sum(np.abs(psi[:, physical])**2))


def measure_gamma(l_max=1, n=8, t_final=400.0, fit_from=50.0):
    tower = Tower(spin=0.9, alpha=0.5, m=1, l_max=l_max,
                  x_left=-200.0, x_right=250.0, n=n,
                  x_pml_center=-100.0, pml_width=10.0)

    # transient-free seed: psi = eigenmode profile, pi = -i omega psi, chi = 0
    omega, psi_profile = solve_qbs(tower)
    psi = psi_profile.copy()
    pi = -1j * omega * psi_profile
    chi = np.zeros_like(psi)

    dt = 0.8 * tower.dx
    nsteps = int(t_final / dt)
    # measurement region: everything to the right of the PML bump
    physical = tower.x > tower.x_pml_center + tower.pml_width / 2.0

    times = np.empty(nsteps)
    norms = np.empty(nsteps)
    state = (psi, pi, chi)
    for step in range(nsteps):
        times[step] = step * dt
        norms[step] = interior_norm(state[0], physical)
        state = tower.rk4_step(state, dt)

    # Gamma = slope of log(norm) over the fit window
    window = times > fit_from
    gamma, _ = np.polyfit(times[window], np.log(norms[window]), 1)
    return tower, times, norms, gamma, omega


if __name__ == "__main__":
    w_leaver = None
    for l_max in (1, 3, 5):
        tower, times, norms, gamma, omega = measure_gamma(l_max=l_max)
        if w_leaver is None:                      # l_max-independent; solve once
            w_leaver = leaver_benchmark(tower)
        print(f"l_max={l_max} ({tower.n_channels} channels {tower.l_values}):  "
              f"fitted Gamma={gamma:+.4e} ({gamma/w_leaver.imag:.4f}x)  "
              f"eigenvalue Gamma={omega.imag:+.4e} ({omega.imag/w_leaver.imag:.4f}x)")
    print(f"Leaver (spheroidal): omega_R={w_leaver.real:.6f}  Gamma={w_leaver.imag:+.4e}")
