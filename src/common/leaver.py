#!/usr/bin/env python
"""
Leaver continued-fraction solver for massive-scalar quasi-bound states on Kerr.

Reference: S. R. Dolan, "Instability of the massive Klein-Gordon field on the
Kerr spacetime", Phys. Rev. D 76, 084001 (2007), arXiv:0705.2880 — eqs. (33)-(46).

This is the Project 1.1 *linear input*: it returns the complex frequency
omega = E_nlm + i Gamma_nlm of the fastest-growing superradiant bound state
(default n=2, l=m=1) around a Kerr BH of mass M=1 and spin a, for a scalar of
mass mu. In these units the gravitational fine-structure constant is
    alpha = M mu = mu.
A bound state has Re(omega) < mu and decays at infinity (Re(q) < 0); it is
*superradiant/growing* when Gamma_nlm > 0, which requires E_nlm < m Omega_H.

The NR sweep (Project 1.2-A, see ../PROJECT_1.2A_LINEAR_RATES.md) benchmarks its
3+1D E_nlm, Gamma_nlm against the values returned here.

Angular eigenvalue: c^2 = a^2(omega^2 - mu^2) is tiny for bound states
(omega ~ mu), so Lambda_lm = l(l+1) to leading order; the O(c^2) spheroidal
correction is ~1e-3 here and below the NR comparison precision. Set
`spheroidal=True` to include the leading Seidel term as a check.

Requires mpmath (arbitrary precision — Gamma_nlm ranges over ~1e-7..1e-13).
"""
import argparse
import mpmath as mp

mp.mp.dps = 30  # decimal digits of working precision


def _coeffs(omega, alpha, a, m, l, Lambda):
    """Dolan (2007) intermediate constants c0..c4, M=1, mu=alpha."""
    mu = alpha
    I = mp.mpc(0, 1)
    b = mp.sqrt(1 - a * a)          # (r+ - r-)/2 ;  Dolan eq (45)
    # q: bound state -> Re(q) < 0 so R ~ e^{qr} decays at infinity
    q = mp.sqrt(mu * mu - omega * omega)
    if mp.re(q) > 0:
        q = -q
    # The Dolan eq (34) exponents sigma = 2 rp (omega-wc)/(rp-rm) and
    # chi = (mu^2 - 2 omega^2)/q set the R ansatz prefactor, not the c_n below,
    # so they are not needed here.

    w = omega
    wam = w - a * m / 2                          # recurring (omega - a m / 2)
    wiq = w - I * q                              # recurring (omega - i q)

    c0 = 1 - 2 * I * w - (2 * I / b) * wam
    c1 = (-4 + 4 * I * (w - I * q * (1 + b)) + (4 * I / b) * wam
          - 2 * (w * w + q * q) / q)
    c2 = 3 - 2 * I * w - 2 * (q * q - w * w) / q - (2 * I / b) * wam
    c3 = (2 * I * wiq**3 / q + 2 * wiq**2 * b + q * q * a * a + 2 * I * q * a * m
          - Lambda - 1 - wiq**2 / q + 2 * q * b
          + (2 * I / b) * (wiq**2 / q + 1) * wam)
    c4 = (wiq**4 / q**2 + 2 * I * w * wiq**2 / q
          - (2 * I / b) * (wiq**2 / q) * wam)
    return c0, c1, c2, c3, c4


def _cf_residual(omega, alpha, a, m, l, Lambda, N):
    """f(omega) = beta_0 - alpha_0 * CF; roots are the bound-state frequencies.

    alpha_n = n^2 + (c0+1)n + c0     (Dolan eq 37)
    beta_n  = -2n^2 + (c1+2)n + c3   (eq 38)
    gamma_n = n^2 + (c2-3)n + c4     (eq 39)
    CF = gamma_1/(beta_1 - alpha_1 gamma_2/(beta_2 - ...)), by backward recursion.
    """
    c0, c1, c2, c3, c4 = _coeffs(omega, alpha, a, m, l, Lambda)

    def alpha_n(n):
        return n * n + (c0 + 1) * n + c0

    def beta_n(n):
        return -2 * n * n + (c1 + 2) * n + c3

    def gamma_n(n):
        return n * n + (c2 - 3) * n + c4

    # backward recurrence for the tail continued fraction u_1
    u = mp.mpf(0)
    for n in range(N, 0, -1):
        u = gamma_n(n) / (beta_n(n) - alpha_n(n) * u)
    return beta_n(0) - alpha_n(0) * u


def _spheroidal_leading(alpha, a, l, m, omega):
    """Leading Seidel O(c^2) correction to Lambda_lm (s=0). Small for bound
    states; used only when spheroidal=True."""
    c2 = a * a * (omega * omega - alpha * alpha)
    # s=0 first-order coefficient (Abramowitz-Stegun / Seidel 1989):
    f1 = (2 * m * m - 2 * l * (l + 1) + 1) / ((2 * l - 1) * (2 * l + 3))
    return l * (l + 1) + f1 * c2


def solve(alpha, a, l=1, m=1, n=2, N=500, spheroidal=False, guess=None):
    """Return complex omega for the (n,l,m) massive-scalar bound state on Kerr.

    alpha = M mu (=mu, M=1); a = spin (0<=a<1); n = overtone (n=2 is the
    fastest-growing l=m=1 mode). N = continued-fraction depth.
    """
    mu = alpha
    if guess is None:
        # hydrogenic seed: E_nlm ~ mu(1 - alpha^2/2n^2)
        guess = mu * (1 - alpha * alpha / (2 * n * n))
    guess = mp.mpf(guess) if not isinstance(guess, mp.mpc) else guess

    Lambda = l * (l + 1)

    def f(omega):
        L = (_spheroidal_leading(alpha, a, l, m, omega)
             if spheroidal else Lambda)
        return _cf_residual(omega, alpha, a, m, l, L, N)

    # Muller needs 3 points; cluster them TIGHTLY around the hydrogenic guess so
    # it converges to the physical bound state and does not wander to a spurious
    # root elsewhere in the complex plane. Muller then resolves the tiny Im(omega).
    seeds = (mp.mpc(guess, 0), mp.mpc(guess + 1e-5, 0), mp.mpc(guess, 1e-6))
    omega = mp.findroot(f, seeds, solver='muller')
    return complex(omega)


# ------------------------------- validation --------------------------------

def _validate():
    print("=== Leaver massive-scalar Kerr bound states: validation ===\n")

    # 1) hydrogenic E_nlm at small alpha (a=0.9, l=m=1, n=2)
    print("[1] E_nlm -> hydrogenic mu(1-alpha^2/8) at small alpha (a=0.9):")
    for alpha in (0.05, 0.1, 0.2):
        w = solve(alpha, 0.9)
        hyd = alpha * (1 - alpha**2 / 8)
        print(f"    alpha={alpha:<5} E_nlm={w.real:.8f}  hydrogenic={hyd:.8f}"
              f"  (rel {abs(w.real-hyd)/hyd*100:.2f}%)  Gamma_nlm={w.imag:.3e}")

    # 2) sign: Schwarzschild (a=0) must DECAY (Gamma_nlm<0); spinning GROWS
    print("\n[2] sign of Gamma_nlm (l=m=1, alpha=0.3):")
    for a in (0.0, 0.5, 0.9, 0.99):
        w = solve(0.3, a)
        Omega_H = a / (2 * (1 + (1 - a*a)**0.5))
        sr = "superradiant (wR<mOmega_H)" if w.real < Omega_H else "not superradiant"
        print(f"    a={a:<5} Gamma_nlm={w.imag:+.3e}  Omega_H={Omega_H:.4f}  {sr}")

    # 3) small-alpha scaling: Detweiler Gamma_nlm ~ (m Omega_H - E_nlm)(mu M)^(4l+4);
    #    the (mu M)^(4l+4) factor is alpha^8 for l=1, so the local log-slope ~8
    #    (mildly above 8 from the varying threshold factor at finite alpha).
    print("\n[3] growth-rate scaling: d ln(Gamma_nlm)/d ln(alpha) -> ~4l+4 = 8 (l=1):")
    import math
    als = [0.08, 0.10, 0.12]
    ws = [solve(al, 0.9) for al in als]
    for i in range(len(als) - 1):
        slope = (math.log(ws[i+1].imag) - math.log(ws[i].imag)) / \
                (math.log(als[i+1]) - math.log(als[i]))
        print(f"    alpha {als[i]}->{als[i+1]}: slope={slope:.2f}")

    # 4) Dolan's peak: max Gamma_nlm ~ 1.5e-7 near alpha~0.42 at a=0.99
    print("\n[4] peak growth rate at a=0.99 (Dolan: ~1.5e-7 near alpha~0.42):")
    best = (0, 0)
    for i in range(30, 55, 3):
        alpha = i / 100
        w = solve(alpha, 0.99)
        if w.imag > best[1]:
            best = (alpha, w.imag)
        print(f"    alpha={alpha:.2f}  E_nlm={w.real:.5f}  Gamma_nlm={w.imag:.4e}")
    print(f"    -> peak Gamma_nlm={best[1]:.4e} at alpha={best[0]:.2f}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--validate", action="store_true", help="run validation suite")
    ap.add_argument("--alpha", type=float, help="M mu")
    ap.add_argument("--a", type=float, default=0.9, help="spin (default 0.9)")
    ap.add_argument("--l", type=int, default=1)
    ap.add_argument("--m", type=int, default=1)
    ap.add_argument("--n", type=int, default=2)
    ap.add_argument("--N", type=int, default=500, help="continued-fraction depth")
    ap.add_argument("--spheroidal", action="store_true")
    args = ap.parse_args()

    if args.validate or args.alpha is None:
        _validate()
    else:
        w = solve(args.alpha, args.a, args.l, args.m, args.n, args.N,
                  args.spheroidal)
        print(f"alpha={args.alpha} a={args.a} (l,m,n)=({args.l},{args.m},{args.n}): "
              f"E_nlm={w.real:.10f}  Gamma_nlm={w.imag:.6e}")
