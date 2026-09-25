## 1. Purpose and scope

This project reproduces the adiabatic treatment of a self-interacting axion cloud around a rotating black hole described by [Omiya, Takahashi & Tanaka](https://arxiv.org/pdf/2201.04382v3). The aim is to calculate how the cloud’s spatial profile changes as its amplitude grows, and then determine its slow evolution from the balance between energy extracted from the black hole and energy radiated outward.

The calculation has two parts. First, we solve for a family of approximately stationary field configurations, parameterized by the fundamental mode’s outer amplitude $A_0$. Each configuration has its own radial profiles and frequency. Second, we calculate the cloud energy and boundary fluxes along this family to determine $A_0(t)$.

We build the implementation in stages:

1. Define the mode expansion and angular basis.
2. Construct the projected radial equations and boundary conditions.
3. Validate the linear limit against independent frequency calculations.
4. Add self-interactions and continue the solutions in amplitude.
5. Calculate energy and fluxes, then evolve the amplitude.

The current Julia implementation includes linear matching and its comparison with Leaver’s continued-fraction method, nonlinear continuation for three potentials, retained profiles, and energy/flux diagnostics. Slow amplitude evolution remains to be implemented. See [the reproduction status](adiabatic-repro.md) for validation limits and [the campaign README](../runs/adiabatic/fig2_20260922/README.md) for saved results. The initial nonlinear target is $a/M=0.99$, $M\mu=0.42$, starting from the fundamental $l_0=m_0=1$ cloud on a fixed Kerr background.

## 2. Field expansion and conventions

We begin with the paper’s mode expansion, Eq. (16):

$$ \phi(t,r,\theta,\varphi;A_0) = \sum_{n\geq1}\sum_{l\geq nm_0} e^{-in[\omega_0(A_0)t-m_0\varphi]} \widetilde R_{nl}(r;A_0) Y_{l nm_0}(\cos\theta) +\mathrm{c.c.} $$

Here, $A_0$ labels a configuration, $\omega_0$ is its fundamental frequency, and $\widetilde R_{nl}$ are its complex radial profiles. Adding the complex conjugate makes the field real.

The integer $n$ labels harmonics of the fundamental pattern, that is, *not* the radial overtone convention of Dolan. Each harmonic has

$$ m_n=nm_0,\qquad \omega_n=n\omega_0. $$

Throughout the radial equations and boundary conditions, $m_0$ denotes the fundamental azimuthal number and $m_n$ denotes the azimuthal number of its $n$th harmonic. For example, with $m_0=1$, the third harmonic has $m_3=3$ and $\omega_3=3\omega_0$. The code's mode index `m` corresponds to $m_n$, while `fundamental_m` corresponds to $m_0$.

We use natural units and set $M=1$. Consequently, the code parameters satisfy

$$ \texttt{spin}=a/M,\qquad \texttt{alpha}=M\mu. $$

### Angular normalization

The polar harmonics take on their standard values

$$ \sqrt{\frac{(2l+1)(l-m)!}{2(l+m)!}}\, P_l^m(x) $$

with canonical normalization. Here $m$ is a generic angular-basis index; for the $n$th active harmonic, substitute $m=m_n$. The azimuthal dependence is stored separately for computational reasons. We precompute harmonics for all $0\leq m\leq l\leq l_{\max}$, and select the active modes separately. Our present linear calculation retains $(l,m_n)=(1,1),(3,1),(5,1)$, all with $n=1$ and $m_0=1$.

### Radial and amplitude conventions

The code’s radial variable $R$ is the paper’s $\widetilde R$, with no additional $1/r$ rescaling. Following Eq. (20), $A_0$ is the coefficient of the fundamental mode’s decaying outer solution. For $M=1$,

$$ \widetilde R_{1l_0}(r;A_0) \sim A_0 e^{-qr} r^{-1+(2\omega_0^2-\mu^2)/q}, \qquad q=\sqrt{\mu^2-\omega_0^2},\quad \mathrm{Re}(q) >0. $$

In the linear problem, the overall amplitude is arbitrary, so we fix this coefficient to one. In the nonlinear problem, $A_0$ is the continuation parameter that determines the configuration.


## 3. Projected radial equations

Substituting the mode expansion into the field equation and projecting onto the angular basis gives a coupled system of radial equations. We write the system corresponding to PDF Eq. (21) as
$$ \Delta \frac{d^2R_l}{dr^2} +\Delta'\frac{dR_l}{dr} +B_l R_l +a^2(\omega_n^2-\mu^2) \sum_{l'}C^{m_n}_{ll'}R_{l'} +N_l=0. \qquad \qquad (*) $$

Here we suppress the harmonic index on $R_l$ while considering a fixed block with $m_n=nm_0$ and $\omega_n=n\omega_0$. Derivatives in this equation are with respect to the Boyer–Lindquist radius $r$.

### Local radial coefficients

In units $M=1$, define

$$ \Delta=r^2-2r+a^2, \qquad \Delta'=2r-2, \qquad K=(r^2+a^2)\omega_n-am_n. $$

The coefficient multiplying each radial profile locally is

$$ B_l = \frac{K^2}{\Delta} -\mu^2r^2 +2am_n\omega_n -a^2\omega_n^2 -l(l+1). $$

This coefficient collects the frequency, mass, rotation, and angular-momentum terms that act on $R_l$ without mixing it with another angular mode.

### Angular coupling

The remaining linear angular dependence enters through

$$ C^{m_n}_{ll'} = \int_{-1}^{1} Y_{l m_n}(x)\,x^2\,Y_{l'm_n}(x)\,dx. $$

Multiplication by $x^2=\cos^2\theta$ couples a mode to itself and to its $l\pm2$ neighbors. It preserves $m_n$, so different azimuthal blocks remain independent in the linear problem.

The term $N_l$ denotes the projected nonlinear remainder after the mass term has been included in the linear operator. The linear solver sets $N_l=0$. The nonlinear solver reconstructs the real angular field and projects $\mu^2(r^2+a^2\cos^2\theta)[\phi-U'(\phi)]$ onto the retained modes. The same potential selector is used in the energy integral.

### Conversion to the tortoise coordinate

Near the horizon, evaluating the radial equation directly introduces divisions by the small quantity $\Delta$. We instead integrate in the tortoise coordinate $r_*$, defined by

$$ \frac{dr_*}{dr}=\frac{r^2+a^2}{\Delta}. $$

Let $A=r^2+a^2$ and $Q_l=\frac{dR_l}{dr_*}.$ The chain rule then gives

$$ \Delta\frac{dR_l}{dr}=AQ_l, $$

and therefore

$$ \frac{d}{dr} \left(\Delta\frac{dR_l}{dr}\right) = 2rQ_l+\frac{A^2}{\Delta}\frac{dQ_l}{dr_*}. $$

Substituting this expression into the radial equation produces two first-order equations:

$$ \frac{dR_l}{dr_*}=Q_l, $$

$$ \frac{dQ_l}{dr_*} = -\frac{2r\Delta}{A^2}Q_l -\frac{\Delta}{A^2} \left[ B_l R_l +a^2(\omega_n^2-\mu^2) \sum_{l'}C^{m_n}_{ll'}R_{l'} +N_l \right]. $$

Before evaluating these coefficients numerically, we simplify the apparently singular product:

$$ \frac{\Delta}{A^2}B_l = \left(\frac{K}{A}\right)^2 + \frac{\Delta}{A^2} \left[ -\mu^2r^2+2am_n\omega_n-a^2\omega_n^2-l(l+1) \right]. $$

This expression remains finite at the horizon and avoids computing a large $K^2/\Delta$ only to multiply it by $\Delta$ afterward.

### Integration state

For $N_{\mathrm{modes}}$ retained angular modes, we pack the variables as

$$ \mathbf{u} = (R_1,\ldots,R_{N_{\mathrm{modes}}}, Q_1,\ldots,Q_{N_{\mathrm{modes}}})^\mathsf{T}. $$

Both halves follow the same explicit mode ordering. The radial right-hand-side function evaluates $d\mathbf{u}/dr_*$, which the integrator uses to advance the profiles through space.

## 4. Boundary conditions and radial integration

We seek radial profiles that are ingoing at the event horizon and decaying or outgoing at infinity, following PDF Eqs. (22)–(23). These conditions select the physically allowed solutions of the radial equations. The horizon radii and angular velocity are
$$ r_\pm=1\pm\sqrt{1-a^2}, \qquad \Omega_H=\frac{a}{r_+^2+a^2}. $$

Because numerical integration begins at finite radii, we use asymptotic series to supply both the field values $R_l$ and their tortoise-coordinate slopes $Q_l$ at each endpoint. The nonlinear terms in $(*)$ are subleading as $r \rightarrow r_+$ and $r \rightarrow \infty$; we tackle each of these cases separately.

### Ingoing horizon boundary: physics
Define $s=r-r_+$ and $g=r_+-r_-$. In Eq. (22), the printed $m$ inside terms of the form $n(\omega_0-m\Omega_H)$ must be read as the fundamental $m_0$ to agree with the expansion in Eq. (16). We write it explicitly as $m_0$ below and retain $M$ to match the paper; the implementation sets $M=1$.

As $r \rightarrow r_+$, equation (22) says
$$ \tilde{R}_{nl} \rightarrow H_l(r) s^{p_H}, $$
where
$$ H_l(r) = A_{nl}^{\text{in}} \left(\frac{r-r_-}{M}\right)^{in\frac{2Mr_-}{g}(\omega_0-m_0\frac{a}{2Mr_-}) - (1-2nM\omega_0i)-iM (\mu^2-2n^2\omega_0^2)/\sqrt{n^2\omega_0^2-\mu^2}} \times e^{i\sqrt{n^2\omega_0^2 - \mu^2}(r-r_-)} \times M^{-p_H} $$
collects all the factors that are smooth and finite at $r=r_+$ and
$$ p_H = -in\frac{2Mr_+}{g}(\omega_0 - m_0\Omega_H). $$

Using $2Mr_+=r_+^2+a^2$, $\omega_n=n\omega_0$, and $m_n=nm_0$ gives

$$ \begin{aligned} p_H &=-i\frac{r_+^2+a^2}{g}(n\omega_0-nm_0\Omega_H)\\ &=-i\frac{r_+^2+a^2}{g}(\omega_n-m_n\Omega_H)\\ &=-i\frac{(r_+^2+a^2)\omega_n-am_n}{g}. \end{aligned} $$

Thus the factor $n$ multiplies both the frequency and the fundamental azimuthal number; it is absorbed into $\omega_n$ and $m_n$ together. To express this behavior in the tortoise coordinate, start from
$$ \frac{dr_*}{dr} = \frac{r^2+a^2}{\Delta}, \qquad \Delta=s(s+g). $$

Near the horizon, $r^2+a^2\to r_+^2+a^2$ and $s+g\to g$, giving
$$ \frac{dr_*}{ds} \sim \frac{r_+^2+a^2}{gs}. $$

We integrate this and find
$$ r_* = \frac{r_+^2+a^2}{g}\ln s +\text{constant}+O(s). $$

Because of this, and the fact that $H_l$ approaches a constant as we near the horizon, we have
$$ R_l\propto e^{-i(\omega_n-m_n\Omega_H)r_*}. $$

Differentiating with respect to $r_*$ then gives
$$ Q_l=\frac{dR_l}{dr_*} \sim-i(\omega_n-m_n\Omega_H)R_l. $$

### Ingoing horizon boundary: integration
The numerical integration begins at a finite coordinate $r_{*,\mathrm{inner}}$, corresponding to a small but nonzero $s=r-r_+$. We therefore include corrections to the asymptotic ingoing behavior by writing
$$ R_l(s)=s^{p_H}P_l(s), \qquad P_l(s)=\sum_{j=0}^{J_H}v_{lj}s^j. $$

Here $J_H$ is the highest retained power in the correction polynomial. Our implementation uses $J_H=3$. The leading coefficients $v_{l0}$ are the free horizon amplitudes. The higher coefficients are determined by substituting this expansion into the coupled linear radial equations and equating powers of $s$. To obtain the initial slope, differentiate:
$$ \frac{dR_l}{ds} = s^{p_H-1}\left[p_HP_l(s)+sP_l'(s)\right], $$

where $P_l'(s)$ means differentiation with respect to $s$. Since $ds/dr=1$ and $\Delta=s(s+g)$, the tortoise-coordinate slope is
$$ Q_l(s) = \frac{\Delta}{r^2+a^2}\frac{dR_l}{ds} = \frac{s+g}{r^2+a^2}s^{p_H} \left[p_HP_l(s)+sP_l'(s)\right]. $$

Evaluating these expressions at the inner endpoint supplies the initial $R_l$ and $Q_l$ for outward integration. The final expression cancels the potentially troublesome factor $1/s$ before numerical
evaluation.

### Outgoing outer boundary: physics

At large radius, define the radial wave number by
$$ k_n^2=\omega_n^2-\mu^2. $$

Combining the radial powers in Eq. (23), its asymptotic form is
$$ \widetilde R_{nl} \sim A_{nl}^{\mathrm{out}} e^{ik_n r} \left(\frac{r}{M}\right)^{p_\infty} \left[ 1+\frac{a_1}{r/M} +\frac{a_2}{(r/M)^2} +\cdots \right], $$

where
$$ p_\infty = -1-iM\frac{\mu^2-2\omega_n^2}{k_n} = -1+iM\frac{2\omega_n^2-\mu^2}{k_n}. $$

Here $\omega_n=n\omega_0$. We have retained $M$ to match the paper; the implementation sets $M=1$.

The choice of square-root branch determines the physical boundary condition. For a real positive frequency above the mass threshold, $\omega_n>\mu$, we choose $k_n>0$. With time dependence $e^{-i\omega_n t}$, the factor
$$ e^{i(k_n r-\omega_n t)} $$

then describes an outward-moving wave. For complex frequencies, we
continue this outgoing branch.

For a bound mode, we instead choose

$$ k_n=iq_n, \qquad q_n=\sqrt{\mu^2-\omega_n^2}, \qquad \mathrm{Re}(q_n)>0. $$

The exponential becomes

$$ e^{ik_n r}=e^{-q_n r}, $$

which decays at infinity. Thus the same asymptotic form supplies
either an outgoing radiative boundary or a decaying bound-state
boundary, depending on the wave-number branch. The present linear
cloud calculation uses the decaying branch.

The exact threshold $k_n=0$ requires a different expansion and is
excluded from the current implementation.

### Outgoing outer boundary: integration

The numerical outer boundary is at a finite coordinate
$r_{*,\mathrm{outer}}$, corresponding to a large Boyer–Lindquist
radius $r$. Returning to $M=1$, we write the truncated outer series as

$$ R_l(r)=E(r)W_l(z), \qquad E(r)=e^{ik_n r}r^{p_\infty}, \qquad z=\frac{1}{r}, $$

with

$$ W_l(z)=\sum_{j=0}^{J_\infty}w_{lj}z^j. $$

The leading coefficients $w_{l0}$ are the free outer amplitudes.
Higher coefficients are determined by substituting this expansion
into the coupled linear radial equations and matching powers of
$1/r$. Our implementation uses $J_\infty=7$, following the expansion
order stated below Eq. (23).

To calculate the slope, note that

$$ \frac{dE}{dr} = \left(ik_n+\frac{p_\infty}{r}\right)E, \qquad \frac{dz}{dr}=-\frac{1}{r^2}. $$

The product and chain rules therefore give

$$ \frac{dR_l}{dr} = E(r) \left[ \left(ik_n+\frac{p_\infty}{r}\right)W_l(z) -\frac{W_l'(z)}{r^2} \right], $$

where $W_l'(z)$ means differentiation with respect to $z$, not $r$.
Converting to the tortoise-coordinate slope,

$$ Q_l(r) = \frac{\Delta}{r^2+a^2}E(r) \left[ \left(ik_n+\frac{p_\infty}{r}\right)W_l(1/r) -\frac{W_l'(1/r)}{r^2} \right]. $$

Evaluating $R_l$ and $Q_l$ at the outer endpoint supplies the initial
state for inward integration. The outer radius and series order
must be checked for convergence, particularly for weakly bound
modes with a long decay length.

### Why integrate from both boundaries?

The radial problem imposes conditions at both ends of the domain. For a trial frequency and trial boundary amplitudes, the horizon series supplies an initial state satisfying the ingoing condition, while the outer series supplies a separate initial state satisfying the decaying or outgoing condition. Neither initial state alone guarantees that the solution satisfies the condition at the opposite boundary.

For a bound mode, a general large-radius solution contains both decaying and growing exponentials. If we integrate outward across this region, a small numerical admixture of the growing solution can overwhelm the desired decaying solution. Starting with the decaying boundary data at the outer endpoint and integrating inward makes the desired exponential grow in the direction of integration, reducing this source of contamination.

We therefore integrate outward from $r_{*,\mathrm{inner}}$ and inward from $r_{*,\mathrm{outer}}$ to a common matching coordinate, chosen so that

$$ r_{*,\mathrm{inner}}<r_{*,\mathrm{match}}<r_{*,\mathrm{outer}}. $$

The two integrations represent one solution only when both the profiles and their slopes agree there:

$$ R_l^{\mathrm{inner}}(r_{*,\mathrm{match}})=R_l^{\mathrm{outer}}(r_{*,\mathrm{match}}),\qquad Q_l^{\mathrm{inner}}(r_{*,\mathrm{match}})=Q_l^{\mathrm{outer}}(r_{*,\mathrm{match}}). $$

The labels here identify the boundary from which each solution was integrated. Adjusting the frequency and amplitudes to satisfy these equalities is the matching problem described in Section 5.

### What RK4 advances

During each radial integration, the trial frequency and boundary amplitudes are fixed. The independent variable is $r_*$, and the evolving numerical state is the complex vector of all retained profiles and their slopes:

$$ \mathbf{u}=(R_1,\ldots,R_{N_{\mathrm{modes}}},Q_1,\ldots,Q_{N_{\mathrm{modes}}})^{\mathsf T},\qquad \frac{d\mathbf{u}}{dr_*}=\mathbf{f}(r_*,\mathbf{u}). $$

The first half of $\mathbf{f}$ contains the $Q_l$ values; the second half contains the expressions for $dQ_l/dr_*$ derived in Section 3, including the angular couplings. Thus RK4 advances every retained mode together through space. Physical time dependence is already represented by the mode expansion, and slow evolution of $A_0(t)$ will be a separate calculation.

For a signed radial step $h$, classical RK4 forms four derivative estimates:

$$ \mathbf{k}_1=\mathbf{f}(r_{*,j},\mathbf{u}_j), $$
$$ \mathbf{k}_2=\mathbf{f}(r_{*,j}+h/2,\mathbf{u}_j+(h/2)\mathbf{k}_1), $$
$$ \mathbf{k}_3=\mathbf{f}(r_{*,j}+h/2,\mathbf{u}_j+(h/2)\mathbf{k}_2), $$
$$ \mathbf{k}_4=\mathbf{f}(r_{*,j}+h,\mathbf{u}_j+h\mathbf{k}_3). $$

These give the update

$$ \mathbf{u}_{j+1}=\mathbf{u}_j+\frac{h}{6}(\mathbf{k}_1+2\mathbf{k}_2+2\mathbf{k}_3+\mathbf{k}_4),\qquad r_{*,j+1}=r_{*,j}+h. $$

Here $j$ counts radial steps, and the bold $\mathbf{k}_i$ are derivative vectors, distinct from the outer wave number $k_n$. Outward integration uses $h>0$; inward integration uses $h<0$. The implementation chooses a uniform spacing whose magnitude is no larger than the requested `step` and which lands exactly on the matching coordinate. It returns the final state by default; `save_profile=true` also retains every position and state for observables and inspection.

The standalone matching calculation defaults to $r_{*,\mathrm{inner}}=-100$, $r_{*,\mathrm{outer}}=100$, $r_{*,\mathrm{match}}=5$, and maximum spacing $0.05$, in units $M=1$. The Fig. 13 driver extends the outer domain for the more weakly bound clouds. These choices require convergence checks: the requested spacing is not an error tolerance. In a regime dominated by RK4 integration error, halving the spacing should reduce the accumulated error by approximately $2^4=16$. Boundary truncation and matching errors must be checked separately.


## 5. Matching and amplitude continuation

`solve_linear_mode` adjusts frequency and boundary amplitudes to match both
`R` and `Q`. `solve_nonlinear_mode` fixes the first outer amplitude to the real
value `A0`, packs the remaining complex unknowns into real coordinates, and
solves the mismatch scaled by `A0`. Real coordinates are necessary because
reconstructing a real field introduces complex conjugation.

The Newton solver supports centered finite differences (`:finite`) and
ForwardDiff (`:forwarddiff`), with descent and stagnation checks and a correction
limit. These are solve controls, separate from radial discretization. Campaigns
can use secant or spaced-history cubic prediction, then correct each proposed
amplitude with Newton. Failed trials are not accepted branch points. See the
campaign README for current defaults; saved runs may use older settings.

The available dimensionless potentials and remainders are:

| Selector | U(phi) | phi − U′(phi) |
|---|---|---|
| `:cosine` | `1 - cos(phi)` | `phi - sin(phi)` |
| `:quartic` | `phi^2/2 - phi^4/24` | `phi^3/6` |
| `:sextic` | `phi^2/2 - phi^4/24 + phi^6/720` | `phi^3/6 - phi^5/120` |

The code evaluates the cosine potential and remainder with small-field
cancellation safeguards. Source and energy functions supply the mass factor
`alpha^2` separately.

## 6. Observables and code structure

After accepting a solution, reconstruct profiles and integrate energy on both
radial branches with composite Simpson quadrature. The campaign drivers also
evaluate finite-radius boundary currents. Their signed total flux enters the
planned slow evolution equation `E′(A0) dA0/dt = -F_tot`. That time evolution
is not implemented by radial RK4 or by amplitude continuation.

Load the library with `include("src/adiabatic/load.jl")` from the repository root.
Paths below are relative to `src/adiabatic/`:

| Location | Responsibility |
|---|---|
| `angular.jl`, `tortoise.jl` | Angular basis/quadrature and Kerr geometry |
| `stationary/nonlinear.jl`, `stationary/radial.jl` | Potential projection and radial equations |
| `stationary/boundary_series.jl`, `stationary/boundaries.jl` | Finite-radius boundary states |
| `stationary/integration.jl`, `stationary/matching.jl` | Radial RK4, residuals, and Newton solves |
| `stationary/continuation_predictor.jl` | Predictions from accepted branch history |
| `observables/energy.jl` | Energy density and branch quadrature |
| `drivers/` | Small- and finite-amplitude continuation, branch derivatives, profile inspection, and Fig. 13 points |
| `fig13_points.py`, `plot_profiles.py` | Python adapter and profile plotting |

The finite-amplitude continuation drivers are `drivers/continue_cosine.jl` and
`drivers/continue_polynomial.jl`. Campaign submit scripts, seeds, results, and
diagnostics live under `runs/adiabatic/<campaign>/`. Tests live in
`tests/adiabatic/`.
