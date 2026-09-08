"""Localize the 1.26% Gamma overshoot (a=0.9, alpha=0.5, l=m=1, deep PML).

The convergence test showed Gamma is converged in dx but sits 1.26% above Leaver,
while omega_R is dead-on. So the systematic lives in the imaginary part. Two cheap
decisive checks, both at n=8:

  (1) eigenvalue-Gamma vs fitted-Gamma.
      - fitted Gamma      = time-domain slope of log(interior L2 norm)  [Dolan]
      - eigenvalue Gamma  = Im(omega) from solve_qbs  [the operator's own QBS]
      Agreement at 1.0126x  => bias is in the discretized operator (PML/BC).
      Eigenvalue clean, fit high => bias is the evolution/fit (window/transient).

  (2) PML-parameter sensitivity (Dolan Sec. III.G, check 6).
      Vary pml_width and target_reflection. If Gamma moves, the layer is
      over/under-absorbing flux; if Gamma is flat, the PML is innocent.
"""
import numpy as np
from solver import Tower, pml_profile
from qbs_eigenmode import solve_qbs
from measure_gamma import leaver_benchmark


def run_one(n=8, pml_width=10.0, target_reflection=1e-6, x_pml_center=-100.0,
            t_final=400.0, fit_from=50.0):
    """Return (omega, fitted_gamma) for one (pml_width, target_reflection)."""
    tower = Tower(spin=0.9, alpha=0.5, m=1, l_max=1,
                  x_left=-200.0, x_right=250.0, n=n,
                  x_pml_center=x_pml_center, pml_width=pml_width)
    # override the damping profile so target_reflection can be swept
    tower.pml_gamma = pml_profile(tower.x, x_pml_center, pml_width,
                                  target_reflection=target_reflection)

    omega, psi_profile = solve_qbs(tower)
    psi = psi_profile.copy()
    pi = -1j * omega * psi_profile
    chi = np.zeros_like(psi)

    dt = 0.8 * tower.dx
    nsteps = int(t_final / dt)
    physical = tower.x > x_pml_center + pml_width / 2.0

    times = np.empty(nsteps)
    norms = np.empty(nsteps)
    state = (psi, pi, chi)
    for step in range(nsteps):
        times[step] = step * dt
        norms[step] = np.sqrt(np.sum(np.abs(state[0][:, physical])**2))
        state = tower.rk4_step(state, dt)

    window = times > fit_from
    fitted_gamma, _ = np.polyfit(times[window], np.log(norms[window]), 1)
    return omega, fitted_gamma


if __name__ == "__main__":
    # runtime benchmark: evaluate leaver.py at this case's parameters (no constant)
    ref_tower = Tower(spin=0.9, alpha=0.5, m=1, l_max=1, x_left=-200.0,
                      x_right=250.0, n=8, x_pml_center=-100.0, pml_width=10.0)
    w_leaver = leaver_benchmark(ref_tower)
    print(f"Leaver (spheroidal, runtime): omega_R={w_leaver.real:.6f}  "
          f"Gamma={w_leaver.imag:+.4e}\n")

    print("=== (1) eigenvalue-Gamma vs fitted-Gamma  (n=8, default PML) ===")
    omega, fitted = run_one()
    print(f"  omega_R           = {omega.real:.6f}  ({omega.real/w_leaver.real:.5f}x)")
    print(f"  eigenvalue Gamma  = {omega.imag:+.4e}  ({omega.imag/w_leaver.imag:.4f}x)")
    print(f"  fitted     Gamma  = {fitted:+.4e}  ({fitted/w_leaver.imag:.4f}x)")

    print("\n=== (2a) PML width sensitivity  (target_reflection=1e-6) ===")
    for w in (6.0, 10.0, 16.0):
        omega, fitted = run_one(pml_width=w)
        print(f"  width={w:5.1f}:  eig={omega.imag:+.4e} ({omega.imag/w_leaver.imag:.4f}x)"
              f"   fit={fitted:+.4e} ({fitted/w_leaver.imag:.4f}x)")

    print("\n=== (2b) target_reflection sensitivity  (width=10) ===")
    for reflection in (1e-4, 1e-6, 1e-8):
        omega, fitted = run_one(target_reflection=reflection)
        print(f"  R0={reflection:.0e}:  eig={omega.imag:+.4e} ({omega.imag/w_leaver.imag:.4f}x)"
              f"   fit={fitted:+.4e} ({fitted/w_leaver.imag:.4f}x)")
