#!/usr/bin/env python
"""Angular coupling coefficients for the coupled (l, m=1) tower (README Sec 10.5 [B]).

The tower ansatz is  Phi = (1/r) sum_l psi_l(t,r) Y_{l1}(theta) e^{i phi}, m=1 fixed.
Two kinds of angular coupling connect the channels l = 1, 3, 5, ... (odd sector):

  1. LINEAR (metric) coupling: Dolan's cos^2(theta) factor couples l -> l, l+-2.
     Matrix element c_{jl} = <Y_{j1}|cos^2 theta|Y_{l1}>. The DIAGONAL c_{ll} is
     the number the single-mode code already carries (c_11 = 1/5); the
     off-diagonals c_{l,l+-2} are the O(c^2) metric coupling the README notes as
     tiny for a QBS (omega ~ mu). Provided here for completeness / future use.

  2. NONLINEAR (quartic) coupling: the RWA-resonant cubic (3/4)|Psi|^2 Psi,
     projected onto channel j, couples channels through the four-harmonic overlap
        Q^j_{abc} = INT Y*_{j1} Y_{a1} Y*_{b1} Y_{c1} dOmega.
     This is the coupling that lets a pure l=1 cloud SOURCE l=3 at the SAME
     frequency (Q^3_{111} != 0) -- the multi-mode saturation channel that is the
     real deliverable of step [B].

Both are computed by exact Gauss-Legendre quadrature in u = cos(theta): the
theta-parts Lambda_{lm}(u) = N_{lm} P_l^m(u) are polynomials (times sin^m), so a
modest fixed node count integrates them exactly. The phi integral is trivial
(the net azimuthal exponent is 0 for every term kept: m conservation
-1+1-1+1 = 0), giving an overall 2*pi that cancels in the normalized tensor.

Validated (see module __main__): orthonormality = 1, c_11 = 1/5 (matches the
single-mode c_ll), and the normalized tensor Qbar^1_{111} = 1 so a one-channel
tower reproduces evolve.py exactly. The tensor is Hermitian (norm-conserving):
sum_j Im(psi_j* F_j) = 0 to machine precision for the cubic force
F_j = sum_{abc} Qbar^j_{abc} psi_a psi_b* psi_c.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import lpmv, factorial

# 64 Gauss-Legendre nodes: exact for the (low-degree) polynomial integrands here,
# with generous headroom up to l ~ 15.
_U, _W = leggauss(64)


def _Lambda(l, m, u):
    """Real theta-part of the spherical harmonic Y_{lm} (m >= 0):
    Lambda_{lm}(u) = N_{lm} P_l^m(u), with u = cos(theta) and the Condon-Shortley
    phase carried by scipy's lpmv. Normalized so INT |Lambda|^2 2*pi = 1."""
    norm = np.sqrt((2 * l + 1) / (4 * np.pi) * factorial(l - m) / factorial(l + m))
    return norm * lpmv(m, l, u)


def cos2_coupling(j, l, m=1):
    """Linear metric coupling c_{jl} = <Y_{jm}|cos^2 theta|Y_{lm}>.
    Diagonal j=l reproduces the single-mode c_ll (c_11 = 1/5); off-diagonal
    j = l+-2 is the (tiny, O(c^2)) metric channel coupling."""
    return 2 * np.pi * np.sum(_W * _Lambda(j, m, _U) * _Lambda(l, m, _U) * _U**2)


def quartic_overlap(j, a, b, c, m=1):
    """Un-normalized four-harmonic overlap
        Q^j_{abc} = INT Y*_{jm} Y_{am} Y*_{bm} Y_{cm} dOmega
    (all azimuthal indices = m, so the phi integral is 2*pi). This is the
    resonant-cubic angular coupling: channel j is sourced by the product
    psi_a psi_b* psi_c."""
    return 2 * np.pi * np.sum(
        _W * _Lambda(j, m, _U) * _Lambda(a, m, _U) * _Lambda(b, m, _U) * _Lambda(c, m, _U)
    )


def channels(l_max, m=1):
    """The odd-sector channel list [1, 3, 5, ...] up to l_max (the l=m=1 mode
    couples only to l = 3, 5, ... via both cos^2 theta and the cubic)."""
    return list(range(m, l_max + 1, 2))


def coupling_tensors(l_max, m=1):
    """Assemble the coupling data for channels ls = channels(l_max).

    Returns (ls, c_matrix, Qbar) where:
      ls        : list of l values, e.g. [1, 3] for l_max=3.
      c_matrix  : (nl, nl) linear cos^2 theta coupling, c_matrix[i,k] = c_{ls[i],ls[k]}.
      Qbar      : (nl, nl, nl, nl) quartic tensor NORMALIZED by Q^{m}_{mmm}
                  (the single-mode overlap the scalar knob lambda_quartic absorbs),
                  so Qbar[0,0,0,0] = 1 and the one-channel tower reduces exactly to
                  the single-mode evolve.py with the same lambda_quartic. The cubic
                  force on channel j is F_j = sum_{abc} Qbar[j,a,b,c] psi_a psi_b* psi_c.
    """
    ls = channels(l_max, m)
    nl = len(ls)
    c_matrix = np.array([[cos2_coupling(j, k, m) for k in ls] for j in ls])
    Q0 = quartic_overlap(m, m, m, m, m)  # single-mode overlap <|Y_m1|^2 Y_m1|Y_m1>
    Qbar = np.empty((nl, nl, nl, nl))
    for ij, j in enumerate(ls):
        for ia, a in enumerate(ls):
            for ib, b in enumerate(ls):
                for ic, c in enumerate(ls):
                    Qbar[ij, ia, ib, ic] = quartic_overlap(j, a, b, c, m) / Q0
    return ls, c_matrix, Qbar


def _selftest():
    print("orthonormality <Y_11|Y_11>        =",
          2 * np.pi * np.sum(_W * _Lambda(1, 1, _U) ** 2), "(want 1)")
    print("diagonal metric c_11              =", cos2_coupling(1, 1), "(want 0.2)")
    print("diagonal metric c_33              =", cos2_coupling(3, 3))
    print("off-diagonal metric c_13          =", cos2_coupling(1, 3), "(tiny-coupling weight)")
    ls, c_matrix, Qbar = coupling_tensors(l_max=7)
    print("channels up to l_max=7            =", ls)
    print("Qbar^1_111 (normalization)        =", Qbar[0, 0, 0, 0], "(want 1)")
    print("Qbar^3_111 (l=1 SOURCES l=3)      =", Qbar[1, 0, 0, 0], "(want != 0)")
    print("Qbar^5_111 (selection rule)       =", Qbar[2, 0, 0, 0], "(want ~0: 1x1x1 -> l<=3)")
    # Hermiticity / norm conservation of the cubic force, on a random state.
    rng = np.random.default_rng(0)
    psi = rng.standard_normal(len(ls)) + 1j * rng.standard_normal(len(ls))
    force = np.einsum('jabc,a,b,c->j', Qbar, psi, np.conj(psi), psi)
    dnorm = np.sum(np.imag(np.conj(psi) * force))
    print("Hermiticity sum_j Im(psi_j* F_j)  =", dnorm, "(want ~0)")


if __name__ == "__main__":
    _selftest()
