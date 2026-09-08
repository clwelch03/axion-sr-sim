"""Seed generator: the quasi-bound-state eigenmode profile for a Tower.

This is a UTILITY, not the physics.  measure_gamma.py uses Dolan's time-domain
evolution; it just needs a transient-free starting profile, and the cleanest
such profile is the QBS eigenmode.

How it works, in one idea: the linearised system is  d/dt u = G u  with
u = [psi; pi; chi].  A single mode evolves as u ~ e^{lambda t}, so the QBS is an
eigenvector of G with eigenvalue lambda = Gamma - i omega_R.  We assemble the
matrix form of Tower.rhs and ask a sparse eigensolver for the eigenvector whose
frequency sits nearest the expected QBS.  The tower's A couples channels, so we
solve the generalised problem  K v = lambda M v  (M carries A) instead of
inverting A.  The public entry point is solve_qbs(tower).
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


def _stencil_matrices(num_points, dx):
    """Sparse D1, D2 reproducing Tower.first/second_derivative exactly.

    4th-order 5-point centred; the two boundary rows each side stay zero, matching
    the frozen Dirichlet layers in Tower.rhs.
    """
    interior = np.arange(2, num_points - 2)

    def band(offsets_and_coeffs, scale):
        rows, cols, data = [], [], []
        for offset, coeff in offsets_and_coeffs:
            rows.append(interior)
            cols.append(interior + offset)
            data.append(np.full(interior.size, coeff / scale))
        return sp.csr_matrix((np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))),
                             shape=(num_points, num_points))

    first = band([(-2, 1), (-1, -8), (1, 8), (2, -1)], 12 * dx)
    second = band([(-2, -1), (-1, 16), (0, -30), (1, 16), (2, -1)], 12 * dx**2)
    return first, second


def _channel_block_matrix(n_channels, diagonal_for_block):
    """Assemble an (n_channels x n_channels) block matrix whose (i, j) block is a
    diagonal N x N matrix.  diagonal_for_block(i, j) returns the diagonal (length
    N) for that block, or None if the block is empty."""
    blocks = []
    for i in range(n_channels):
        row = []
        for j in range(n_channels):
            diagonal = diagonal_for_block(i, j)
            row.append(None if diagonal is None else sp.diags(diagonal))
        blocks.append(row)
    return sp.bmat(blocks, format="csr")


def _generalized_matrices(tower):
    """Return (K, M) such that the QBS satisfies K v = lambda M v.

    Blocks act on u = [Psi; Pi; Chi], each stacked channel-major (all nodes of
    channel 0, then channel 1, ...).  Built directly from Tower.rhs term by term,
    with the pi-equation multiplied through by A so we never invert A:

        Psi_dot = Pi
        A Pi_dot = (spatial + mu-coupling) Psi + coeff_pi Pi
                   + A ( -Chi' - i m Omega Chi - gamma (Pi + i m Omega Psi) )
        Chi_dot = gamma (Psi' + i m Omega Psi) - (gamma + i m Omega) Chi
    """
    N, n_ch = tower.num_points, tower.n_channels
    D1, D2 = _stencil_matrices(N, tower.dx)
    D1_full = sp.kron(sp.identity(n_ch), D1)
    D2_full = sp.kron(sp.identity(n_ch), D2)
    im_Omega = 1j * tower.m * tower.Omega_H
    I_full = sp.identity(n_ch * N, format="csr")
    zero = sp.csr_matrix((n_ch * N, n_ch * N))

    # a weight that is the same array for every channel -> block-diagonal
    def per_node(weight):
        return sp.diags(np.tile(weight, n_ch))

    gamma = per_node(tower.pml_gamma)

    # A = Sigma2_0 I + a2Delta * c   (the LHS coupling matrix of Eq. 25)
    A = _channel_block_matrix(n_ch, lambda i, j:
        (tower.Sigma2_0 if i == j else 0.0 * tower.Sigma2_0) + tower.a2Delta * tower.c_matrix[i, j]
        if (i == j or tower.c_matrix[i, j] != 0) else None)

    # spatial operator (coeffs are the same per channel; potential is per channel)
    spatial = (per_node(tower.coeff_psi_xx) @ D2_full
               + per_node(tower.coeff_psi_x) @ D1_full
               - sp.diags(tower.potential.ravel()))

    # neighbour-psi coupling  -alpha^2 a2Delta c_{l,l+-2}  (absent for a single channel)
    off_diagonal_c = tower.c_matrix - np.diag(np.diag(tower.c_matrix))
    if off_diagonal_c.any():
        mu_coupling = -tower.alpha**2 * _channel_block_matrix(n_ch, lambda i, j:
            tower.a2Delta * off_diagonal_c[i, j] if off_diagonal_c[i, j] != 0 else None)
    else:
        mu_coupling = zero

    K = sp.bmat([
        [zero,                            I_full,                     zero],
        [spatial + mu_coupling - A @ (im_Omega * gamma),
                                          per_node(tower.coeff_pi) - A @ gamma,
                                                                      -A @ (D1_full + im_Omega * I_full)],
        [gamma @ D1_full + im_Omega * gamma,
                                          zero,                       -(gamma + im_Omega * I_full)],
    ], format="csc")
    M = sp.block_diag([I_full, A, I_full], format="csc")
    return K, M


def solve_qbs(tower, omega_R_guess=0.4782, omega_guess=None, maxiter=None):
    """Return (omega, psi_profile) for the quasi-bound state of `tower`.

    omega is complex (omega_R + i Gamma); psi_profile is shaped
    (n_channels, num_points), normalised to unit peak, ready to seed as
    psi with pi = -i omega psi.

    omega_guess (complex) targets the full frequency omega_R + i Gamma, which is
    needed to pick a specific overtone near threshold -- the fundamental and the
    first overtone share almost the same omega_R and separate only in Gamma, so a
    real-only target cannot tell them apart. Falls back to omega_R_guess (real)
    when omega_guess is None.

    maxiter caps the ARPACK iterations. The scipy default is 10*N (~1e5 here), so
    a point where shift-invert cannot isolate the mode spins for minutes; pass a
    modest maxiter (e.g. 2000) to make such a point raise ArpackNoConvergence fast
    instead of hanging, and let the caller record it and move on.
    """
    K, M = _generalized_matrices(tower)
    # lambda = Gamma - i omega_R = -i * omega, so target lambda = -i * (guess omega)
    guess = omega_guess if omega_guess is not None else omega_R_guess
    target = -1j * guess
    # We need only the eigenpair nearest the accurate complex target.  Asking
    # ARPACK for a wider cluster makes near-threshold solves pathological: the
    # QBS converges quickly, but unnecessary neighbouring modes do not.
    eigenvalues, eigenvectors = spla.eigs(K, M=M, k=1, sigma=target, which="LM",
                                          maxiter=maxiter)

    # pick the eigenvalue closest to the expected QBS frequency
    best = np.argmin(np.abs(eigenvalues - target))
    lam = eigenvalues[best]
    omega = 1j * lam                                   # omega = omega_R + i Gamma

    # first block of the eigenvector is Psi (channel-major) -> (n_channels, N)
    N, n_ch = tower.num_points, tower.n_channels
    psi_profile = eigenvectors[:n_ch * N, best].reshape(n_ch, N)
    psi_profile = psi_profile / np.max(np.abs(psi_profile))
    return omega, psi_profile
