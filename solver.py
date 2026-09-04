"""Time-domain solver for Dolan 1212.1477 Sec. III, built up in steps.

Step 1 (this file, current): a single decoupled angular channel (l = m = 1),
Eq. (25) with the l -> l +/- 2 coupling switched off.  Method of Lines:
4th-order 5-point spatial stencils, RK4 in time, hard Dirichlet (psi = pi = 0)
at both edges.  No PML yet -- the left edge reflects, so step-1 runs must stop
before reflections return.

Coordinates: r is Boyer-Lindquist areal radius, x is the Kerr tortoise
coordinate with dx/dr = (r^2 + a^2) / Delta.  Units: M = 1.

Later steps add: the PML near-horizon region (Eq. 44), then the coupled-l
tower (full Eq. 25) via a per-node Thomas solve.
"""

import numpy as np

def cos2_theta_helper(l, m):
    """
    We use the fact that
    cosθ Y_lm  =  a(l,m) Y_{l+1,m}  +  a(l-1,m) Y_{l-1,m}
    where this function encodes a(l,m).

    :param l: angular momentum quantum number
    :param m: azimuthal quantum number
    """
    if l < abs(m): return 0
    return np.sqrt( ((l+1)**2 - m**2) / ((2*l+1)*(2*l+3)) )

def cos2_theta_matrix(m, l_max):
    """
    The matrix element c_{jl}^m from Dolan's equation 23 is zero unless j ∈ {l-2, l, l+2}.
    We can thus precompute all the values we may need, saving us time in the long run.

    In c_{jl}^m:
        * l is the angular momentum quantum number; we fix it when we decide what mode we're investigating
        * m is the azimuthal quantum number; it only serves to put a lower bound on j
        * j is the index for the source channel; we sum over all j to get the full wavefunction

    May want to vectorize this in the future?
    """
    # construct the tower of possible nonzero l-values for a set m
    l_values = list(range(abs(m), l_max + 1, 2))
    num_l_values = len(l_values)

    matrix = np.zeros((num_l_values, num_l_values))

    for i, l in enumerate(l_values):
        matrix[i, i] = cos2_theta_helper(l, m)**2 + cos2_theta_helper(l-1, m)**2
        if i+1 < num_l_values:
            matrix[i, i+1] = matrix[i+1, i] = cos2_theta_helper(l, m) * cos2_theta_helper(l+1, m)

    return matrix, l_values

class KerrTortoise:
    """Analytic Kerr tortoise map x(r) and its numerical inverse r(x).

    dx/dr = (r^2 + a^2) / Delta,  Delta = (r - r_plus)(r - r_minus).
    Partial fractions give
        x(r) = r + coeff_plus  * ln((r - r_plus ) / M)
                 + coeff_minus * ln((r - r_minus) / M)
    with coeff_plus = 2 M r_plus / (r_plus - r_minus), and coeff_minus the
    same with r_plus <-> r_minus (so coeff_minus < 0).  The integration
    constant is irrelevant here -- we only ever invert x(r).
    """

    def __init__(self, spin, mass=1.0):
        self.mass = mass
        self.spin = spin
        root = np.sqrt(mass * mass - spin * spin)
        self.r_plus = mass + root
        self.r_minus = mass - root
        gap = self.r_plus - self.r_minus
        self.coeff_plus = 2.0 * mass * self.r_plus / gap
        self.coeff_minus = -2.0 * mass * self.r_minus / gap

    def delta(self, r):
        return r**2 - 2 * self.mass * r + self.spin**2 

    def x_r(self, r):
        log_plus_term = self.r_plus * np.log(np.abs(r - self.r_plus) / (2 * self.mass))
        log_minus_term = self.r_minus * np.log(np.abs(r - self.r_minus) / (2 * self.mass))

        return r + (2 * self.mass) / (self.r_plus - self.r_minus) * (log_plus_term - log_minus_term)

    def dx_dr(self, r):
        return (r**2 + self.spin**2) / self.delta(r)

    def r_x(self, x_target, tol=1e-10, max_iter=100):
        """Invert x(r) by Newton's method, vectorized over x_target.

        Grid is uniform in x; we need r at each node.  Seed from the two
        analytic limits -- far out x ~ r, near-horizon x ~ coeff_plus *
        ln((r - r_plus)/M) -- then Newton with r' = 1/dx_dr.  Clamp r above
        the horizon each step so the log arguments stay positive.

        Convergence is on the *relative* x-residual |res| / (1 + |x|).  Near
        the horizon dx/dr blows up, so a relative error in r at the float64
        floor leaves an absolute x-residual ~1e-9 no matter what; the relative
        criterion is scale-free and reaches the float64 floor everywhere.
        """
        x_target = np.asarray(x_target, dtype=float)
        near_horizon_seed = self.r_plus + self.mass * np.exp(x_target / self.coeff_plus)
        far_seed = np.maximum(x_target, self.r_plus + 1.0)
        r = np.where(x_target > 0.0, far_seed, near_horizon_seed)

        for _ in range(max_iter):
            residual = self.x_r(r) - x_target
            r = r - residual / self.dx_dr(r)
            r = np.maximum(r, self.r_plus * (1.0 + 1e-14))
            if np.max(np.abs(residual) / (1.0 + np.abs(x_target))) < tol:
                break
        else:
            raise RuntimeError(
                f"r_x Newton failed to converge: max relative residual "
                f"{np.max(np.abs(residual) / (1.0 + np.abs(x_target))):.3e}"
            )
        return r


class Tower:
    """The coupled-l tower of Dolan Eq. (25) on a uniform tortoise grid.

    One object holds all same-parity channels l in {|m|, |m|+2, ..., l_max}.
    l_max = |m| gives a single rung; the code path is identical, so the
    single-channel validation and the full tower run through exactly the same
    RHS -- there is no separate single-channel version to replace later.

    Channel i is angular momentum l_values[i]; state arrays are shaped
    (n_channels, n_points).  Units: M = 1, and alpha = M*mu is Dolan's mass
    coupling (so mu = alpha here).
    """

    def __init__(self, spin, alpha, m, l_max, x_left, x_right, n, x_pml_start):
        self.mass = 1.0
        self.spin = spin
        self.alpha = alpha
        self.m = m
        self.n = n
        self.dx = self.mass / n
        self.c_matrix, self.l_values = cos2_theta_matrix(m, l_max)
        self.n_channels = len(self.l_values)
        
        self.tortoise = KerrTortoise(spin, self.mass)

        # uniform grid in the tortoise coordinate x, then r at each node
        self.num_points = round((x_right - x_left) / self.dx) + 1
        self.x = x_left + self.dx * np.arange(self.num_points)
        self.r = self.tortoise.r_x(self.x)

        # --- l-independent node arrays (functions of r only) ---
        a = self.spin
        r = self.r
        Delta = self.tortoise.delta(r)
        self.Delta = Delta
        self.a2Delta = a**2 * Delta                       # scales c_matrix inside A
        self.Sigma2_0 = (r**2 + a**2)**2 - a**2 * Delta

        self.coeff_psi_xx = (r**2 + a**2)**2                                    # * psi''
        self.coeff_psi_x = 2j * a * m * (r**2 + a**2) - 2.0 * a**2 * Delta / r  # * psi'
        self.coeff_pi = -4j * a * m * self.mass * r                            # * pi

        # --- l-dependent potential V_l, Eq. (26), one row per channel ---
        # V_l = Delta * [ l(l+1) + (2M/r)(1 - a^2/(Mr)) + 2iam/r
        #                 + mu^2 (r^2 + a^2 c_ll) ]
        self.potential = np.empty((self.n_channels, self.num_points), dtype=complex)
        for i, l in enumerate(self.l_values):
            c_ll = self.c_matrix[i, i]
            self.potential[i] = Delta * (
                l * (l + 1)
                + (2.0 * self.mass / r) * (1.0 - a**2 / (self.mass * r))
                + 2j * a * m / r
                + self.alpha**2 * (r**2 + a**2 * c_ll)
            )

        # PML: horizon angular velocity Omega_H and the damping profile gamma(x).
        # Eq. (44) uses the constant horizon Omega since the PML lives at r -> r_+.
        self.Omega_H = self.spin / (self.tortoise.r_plus**2 + self.spin**2)
        self.pml_gamma = pml_profile(self.x, x_pml_start)

    def first_derivative(self, f):
        """4th-order centered 5-point d/dx, Dolan Eq. (37).

        Acts along the node axis (axis -1), so it works row-wise on the
        (n_channels, n_points) state.  The two boundary layers each side are
        left at zero -- the hard Dirichlet edges never need their derivative.
        """
        d = np.zeros_like(f)
        d[..., 2:-2] = (
            -f[..., 4:] + 8 * f[..., 3:-1] - 8 * f[..., 1:-3] + f[..., :-4]
        ) / (12.0 * self.dx)
        return d

    def second_derivative(self, f):
        """4th-order centred 5-point d^2/dx^2, Dolan Eq. (37)."""
        d = np.zeros_like(f)
        d[..., 2:-2] = (
            -f[..., 4:] + 16 * f[..., 3:-1] - 30 * f[..., 2:-2]
            + 16 * f[..., 1:-3] - f[..., :-4]
        ) / (12.0 * self.dx**2)
        return d

    def build_A_diagonals(self):
        """The three diagonals of A = Sigma2_0 I + a2Delta c, Eq. (25) LHS.

        A is tridiagonal in channel and static in time (geometry only), so we
        rebuild it here from the stored node arrays rather than caching a copy.
        Each diagonal is shaped (n_channels, n_points); the unused ends of sub
        (row 0) and sup (last row) stay zero.
        """
        sub = np.zeros((self.n_channels, self.num_points))
        diag = np.empty((self.n_channels, self.num_points))
        sup = np.zeros((self.n_channels, self.num_points))
        for i in range(self.n_channels):
            diag[i] = self.Sigma2_0 + self.a2Delta * self.c_matrix[i, i]
            if i + 1 < self.n_channels:
                sup[i] = self.a2Delta * self.c_matrix[i, i + 1]
            if i - 1 >= 0:
                sub[i] = self.a2Delta * self.c_matrix[i, i - 1]
        return sub, diag, sup

    def rhs(self, state):
        """Full Eq. (25) interior + Eq. (44) PML source.  state = (psi, pi, chi),
        each complex (n_channels, n_points).  Returns (psi_dot, pi_dot, chi_dot).

        pi_dot = A^{-1} b   (the full coupled Eq. 25, computed everywhere)
                 + [ -chi' - imOmega chi - gamma (pi + imOmega psi) ]   (PML source)
        chi_dot = gamma (psi' - chi + imOmega psi) - imOmega chi.

        Only the three genuinely-new PML terms (those carrying chi or gamma) are
        added; the base 2imOmega(psi' - pi) already falls out of A^{-1} b in the
        near-horizon limit, so adding it again would double-count.  chi and gamma
        are both zero in the interior, so there the equations reduce exactly to
        the two-field Eq. 25 solver.
        """
        psi, pi, chi = state
        psi_x = self.first_derivative(psi)
        psi_xx = self.second_derivative(psi)
        chi_x = self.first_derivative(chi)

        # explicit RHS b_l (everything in Eq. 25 that is not a pi_dot).
        # coeff_* are (n_points,), broadcasting over the channel axis.
        b = (
            self.coeff_psi_xx * psi_xx
            + self.coeff_psi_x * psi_x
            + self.coeff_pi * pi
            - self.potential * psi
        )
        # neighbor-psi coupling: -a^2 mu^2 Delta (c_{l,l+2} psi_{l+2} + c_{l,l-2} psi_{l-2}).
        # the channel-space contraction picks out the l +/- 2 neighbors via
        # the off-diagonal of c; a2Delta already carries a^2 Delta.
        if self.n_channels > 1:
            offdiag_c = self.c_matrix - np.diag(np.diag(self.c_matrix))
            coupling = np.einsum("ij,jn->in", offdiag_c, psi)
            b = b - self.alpha**2 * self.a2Delta * coupling

        # solve A pi_dot = b (tridiagonal in channel, every node at once)
        sub, diag, sup = self.build_A_diagonals()
        pi_dot = thomas_solve(sub, diag, sup, b)

        # PML source (Eq. 44).  imOmega uses the constant horizon Omega; gamma is
        # (n_points,) and broadcasts over channels.  Both gamma and chi vanish in
        # the interior, so this is identically zero there.
        imOmega = 1j * self.m * self.Omega_H
        pi_dot = pi_dot - chi_x - imOmega * chi - self.pml_gamma * (pi + imOmega * psi)
        chi_dot = self.pml_gamma * (psi_x - chi + imOmega * psi) - imOmega * chi

        psi_dot = pi.copy()

        # hard Dirichlet: freeze the two boundary layers each side
        for arr in (psi_dot, pi_dot, chi_dot):
            arr[..., :2] = 0.0
            arr[..., -2:] = 0.0
        return (psi_dot, pi_dot, chi_dot)

    def rk4_step(self, state, dt):
        """One classical RK4 step, generic over the number of fields in state."""
        def add(s, k, factor):
            return tuple(si + factor * ki for si, ki in zip(s, k))

        k1 = self.rhs(state)
        k2 = self.rhs(add(state, k1, dt / 2))
        k3 = self.rhs(add(state, k2, dt / 2))
        k4 = self.rhs(add(state, k3, dt))
        return tuple(
            si + dt / 6 * (a + 2 * b + 2 * c + d)
            for si, a, b, c, d in zip(state, k1, k2, k3, k4)
        )

def pml_profile(x, x_pml_start, target_reflection=1e-6, power=3):
    """Canonical polynomial-graded PML damping profile gamma(x).

    Standard Berenger / Collino-Monk grading: gamma ramps from exactly 0 at
    the interior edge x_pml_start up to gamma_max at the left grid edge,

        gamma(x) = gamma_max * ((x_pml_start - x) / L_pml)^power   (x < x_pml_start)
                    = 0                                              (x >= x_pml_start)

    with L_pml = x_pml_start - x_left.  gamma_max is fixed by the theoretical
    reflection coefficient of a graded layer,

        gamma_max = -(power + 1) ln(R_0) / (2 L_pml),

    valid because the near-horizon advective equation (Eq. 43) has unit
    characteristic speed in the tortoise coordinate.  Starting at exactly 0
    is essential -- a jump in gamma reflects.  Returns an array over nodes.

    :param x_pml_start: interior edge of the PML layer; must sit left of the
        entire physical region so gamma is zero everywhere we measure.
    :param target_reflection: R_0, target residual reflection (default 1e-6).
    :param power: grading exponent p; 3 (cubic) is the usual default.
    """
    pml_thickness = x_pml_start - x[0]
    depth = (x_pml_start - x) / pml_thickness      # 0 at x_pml_start and 1 at left grid edge
    gamma_max = -(power+1) * np.log(target_reflection) / (2 * pml_thickness)

    # compute the profile of the PML
    gamma = np.where(x < x_pml_start, gamma_max * depth**power, 0.0)
    return gamma

def thomas_solve(sub, diag, sup, rhs):
    """Solve a tridiagonal system A x = rhs by the Thomas algorithm.

    A is tridiagonal in the FIRST axis (channel); every operation is
    elementwise over the remaining node axis, so all nodes are solved at once.
    sub/diag/sup/rhs are shaped (n_channels, n_points); sub[0] and sup[-1] are
    unused.  n_channels = 1 degenerates to a single divide.
    """
    n_channels = diag.shape[0]
    c_prime = np.empty_like(rhs)
    d_prime = np.empty_like(rhs)
    x = np.empty_like(rhs)

    c_prime[0] = sup[0] / diag[0]
    d_prime[0] = rhs[0] / diag[0]
    for i in range(1, n_channels):
        denom = diag[i] - sub[i] * c_prime[i - 1]
        c_prime[i] = sup[i] / denom
        d_prime[i] = (rhs[i] - sub[i] * d_prime[i - 1]) / denom

    x[n_channels - 1] = d_prime[n_channels - 1]
    for i in range(n_channels - 2, -1, -1):
        x[i] = d_prime[i] - c_prime[i] * x[i + 1]
    return x
