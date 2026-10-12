# Piece 4 guide: the driven $m = 0$ field

Written 2026-10-02 as a companion to the piece table in
[multimode-adiabatic-design.md](multimode-adiabatic-design.md). It explains
what piece 4 computes, why, and how, assuming only the physics of the two
clouds and none of the code's internals. Equation numbers refer to the PDF
[2211.01949v1](papers/2211.01949v1.pdf). Units are $M = 1$ and $\alpha = M\mu$,
so in code $\mu$ and $\alpha$ are the same number.

## 1. What piece 4 is for

There are two bound clouds on the Kerr black hole:

- the **primary** $\phi_1$, with $l = m = 1$ and frequency $\omega_1$. It is
  nonlinear, and it is the stationary branch solved by `continue_branch.jl`;
- the **secondary** $\phi_2$, with $l = m = 2$ and frequency $\omega_2$. It is
  linear in its own amplitude $A_2$ but feels the nonlinear primary (pieces 1–3).

Both are superradiant: $\omega - m\Omega_H < 0$, where $\Omega_H$ is the
horizon's angular velocity. Each therefore *gains* energy from the black hole
and grows. Something has to stop the growth. In the paper the main brake is the
left diagram of Fig. 3. Two quanta of the primary convert into one quantum of
the secondary plus one quantum of a new field, the **$m = 0$ field**. That field
falls into the black hole and carries energy away for good.

Piece 4 computes this $m = 0$ field and the energy per unit time it carries
into the horizon, $F_H^{(0)}$. In piece 6 that flux is one of the loss terms
$F_E$ in the balance equations (62)–(63). In the perturbative theory it is the
$F_0 E_1^2 E_2$ term of Eq. (36), with $F_0$ defined in Eq. (38) and (A11).

## 2. Where the $m = 0$ field comes from: bookkeeping of phases

Write each cloud's time and azimuth dependence as a single phase (my notation):

$$
\psi_1 \equiv \omega_1 t - \varphi, \qquad \psi_2 \equiv \omega_2 t - 2\varphi ,
$$

where $\varphi$ is the Boyer–Lindquist azimuth (not a field). Up to profiles,
$\phi_1 \sim e^{-i\psi_1} + \text{c.c.}$ and $\phi_2 \sim e^{-i\psi_2} + \text{c.c.}$
A product of fields has the phase that is the sum of its factors' phases. For
the product $\phi_1\phi_1\overline{\phi_2}$ (two primaries and one conjugate
secondary):

$$
-2i\psi_1 + i\psi_2 = -i(2\omega_1 - \omega_2)t + i(2-2)\varphi .
$$

So this product oscillates at

$$
\omega_0 \equiv 2\omega_1 - \omega_2
$$

with azimuthal number $m = 0$. This is what "projecting onto $m = 0$
automatically sets $\omega$ to $2\omega_1 - \omega_2$" means on p. 11. Following
footnote 2, we use real parts: $\omega_0 = 2\,\mathrm{Re}\,\omega_1 - \mathrm{Re}\,\omega_2$.

At $M\mu = 0.42$ in the linear limit, $\omega_1 \approx 0.4088$ and
$\omega_2 \approx 0.4156$, so $\omega_0 \approx 0.402$. This has two
consequences:

- $\omega_0 < \mu = 0.42$: the $m = 0$ field cannot propagate to infinity. Far
  away it decays like $e^{-\kappa_0 r}$, with $\kappa \equiv \sqrt{\mu^2-\omega^2}$
  (my notation) and $\kappa_0 \approx 0.12$.
- $m = 0$ means $\omega_0 - m\Omega_H = \omega_0 > 0$: the field is **not**
  superradiant. Waves that reach the horizon carry positive energy into it.
  This is the loss.

## 3. The nonlinear source: a harmonic of $1 - \cos\phi_1$

The field equation is $\Box\phi = \mu^2\sin\phi$. Put $\phi = \phi_1 + \phi_2 + \phi_0$,
with $\phi_0$ the $m = 0$ field, and expand $\sin\phi$ about $\phi_1$. The terms
linear in $\phi_2$ and $\phi_0$ are

$$
\mu^2\cos\phi_1\,(\phi_2 + \phi_0) = \mu^2(\phi_2+\phi_0) - \mu^2(1-\cos\phi_1)(\phi_2+\phi_0).
$$

The first term is the ordinary mass term. The second is what is new.

At a fixed point in space, $1 - \cos\phi_1$ depends on time and azimuth only
through $\psi_1$, periodically. So it has a Fourier series in $\psi_1$ (my
notation for the coefficients):

$$
1 - \cos\phi_1 = \sum_k [1-\cos\phi_1]_k(r, x)\, e^{-ik\psi_1}, \qquad x = \cos\theta .
$$

- $k = 0$ is the period average $\langle 1-\cos\phi_1\rangle$, which pieces 1–3
  use. It lowers the effective mass seen by the secondary.
- $k = 2$ multiplies $\overline{\phi_2} \sim e^{+i\psi_2}$ to give phase
  $-2i\psi_1 + i\psi_2$, which is exactly the $m = 0$ phase of Section 2.

The source of the $m = 0$ field is therefore

$$
[1-\cos\phi_1]_2 \;\overline{P_2}, \qquad P_2 \equiv \sum_l R^2_l(r)\, Y_{l2}(x),
$$

where $P_2$ (my notation) is the secondary's angular profile at radius $r$. In
the small-field limit $1 - \cos\phi_1 \approx \phi_1^2/2$, so
$[1-\cos\phi_1]_2 \to \tfrac12 P_1^2$ with $P_1 \equiv \sum_l R^1_{l1}Y_{l1}$.
That is the perturbative source of (A5). The nonperturbative calculation simply
keeps the full harmonic.

**Where this is in the paper.** The paper never writes the projected
nonperturbative source out. It is assembled from three places:

- the source $S_r$ of Eq. (52): write $\cos\phi_1 = 1 - (1-\cos\phi_1)$; the
  part of $\mu^2\cos\phi_1\,\phi_2$ at $(\omega_0, m=0)$ is
  $-\mu^2[1-\cos\phi_1]_2\,\overline{P_2}$;
- the projection described on p. 11, done as in Eq. (61);
- its small-field form, written out in (A5) and (A9).

Remove the Green's function and the radial integral from (A9). What remains
matches the code term by term:

| (A9) | code (`dissipative_source`) |
|---|---|
| $\mu^2$ | `alpha^2` |
| $r'^2 + a^2\cos^2\theta'$ | `rho_squared` |
| $\tfrac12 R_1^2 S_1^2$ | `harmonic` $= [1-\cos\phi_1]_2$, whose small-field limit is $\tfrac12 P_1^2$ |
| $R_2^* S_2^*$ | `conj.(secondary_profile)` |
| $\int d\cos\theta'\, S_{l0\omega_0}$ | `sum(weighted_term .* polar_values[mode])` |

There are two differences. The full harmonic replaces $\tfrac12 P_1^2$; that is
the nonperturbative step. The projection is onto spherical $Y_{l0}$ with $l$
coupled through $C_{ll'}$, instead of onto spheroidal $S_{l0\omega_0}$; the
paper makes the same choice in (53), (60), and (61). The overall sign does not
affect the flux.

The other terms do not contribute. The $\tfrac12\sin\phi_1\,\phi_2^2$ term in
Eq. (52) has only odd harmonics of $\psi_1$, so it cannot give $-2i\psi_1 + i\psi_2$.

**How it is computed (step 4a, done).** At $t = 0$, $e^{-ik\psi_1} = e^{+ik\varphi}$,
so $[1-\cos\phi_1]_2$ is the ordinary azimuthal Fourier coefficient. Rebuild
$\phi_1$ on the angular grid, compute $1 - \cos\phi_1$, multiply by
$e^{-2i\varphi}$, and average over $\varphi$. This is `mass_deficit_harmonic(..., 2)`.
The checks (Section 8) confirm the $\tfrac12 P_1^2$ limit and that the angular
grid is fine enough.

## 4. The equation for the $m = 0$ field

Write the field as a sum over angular components:

$$
\phi_0 = \sum_l e^{-i\omega_0 t}\, R^0_l(r)\, Y_{l0}(x) + \text{c.c.}
$$

Project the field equation onto each $Y_{l0}$. In the repo's convention
(`core/radial_operator.jl`), every cloud obeys a radial equation of the same shape:

$$
\frac{d}{dr}\!\left(\Delta \frac{dR_l}{dr}\right) + B_l(r)\, R_l
+ a^2(\omega^2-\mu^2)\sum_{l'} C_{ll'} R_{l'} + N_l = 0,
\qquad \Delta = r^2 - 2r + a^2 .
$$

- $B_l$ is the single-$l$ Kerr radial potential. For $m = 0$ it is
  $B_l = \omega_0^2 (r^2+a^2)^2/\Delta - \mu^2 r^2 - a^2\omega_0^2 - l(l+1)$.
- $C_{ll'}$ is the matrix element of $\cos^2\theta$ between $Y_{l0}$ and
  $Y_{l'0}$. On Kerr the angular dependence is spheroidal, not a single
  spherical harmonic, so each $l$ couples to $l \pm 2$. This is why there are
  several components.
- $N_l$ holds everything beyond the linear Kerr operator. Here it is

$$
N_l = \alpha^2 \!\int_{-1}^{1}\! dx\; Y_{l0}\,\rho^2 \Big(
\langle 1-\cos\phi_1\rangle P_0 \;+\; [1-\cos\phi_1]_2\,\overline{P_2}\Big),
\qquad \rho^2 = r^2 + a^2x^2, \quad P_0 \equiv \sum_l R^0_l Y_{l0}.
$$

  The first term is the primary's averaged background, which the $m = 0$ field
  feels just as the secondary does (the left side of Eq. (51) has $\cos\phi_1$).
  The second is the source from Section 3.

We keep only the $k = 0$ (averaged) part of $1 - \cos\phi_1$ in the operator.
Other harmonics would move $\phi_0$ to other frequencies. For example, $k = 2$
sends it back to $(\omega_2, m = 2)$, a correction to the secondary that Sec. III
neglects. This is the same averaged reading chosen for Eq. (61) in pieces 1–3.

**Which $l$.** The source is even under $x \to -x$, and $Y_{l0}$ has parity
$(-1)^l$, so only even $l$ appear. We keep $l = 0, 2, 4$, the analogue of the
primary's $1, 3, 5$ and the secondary's $2, 4, 6$. The paper says $l = 0$
dominates; we will measure its share.

**Boundary conditions.** These are the same physical conditions as for the
clouds:

- **Horizon:** purely ingoing, $R \sim e^{-i\omega_0 r_*}$ as $r_* \to -\infty$.
  Nothing comes out of a black hole. ($r_*$ is the tortoise coordinate, which
  pushes the horizon to $r_* = -\infty$.)
- **Outer:** decaying, $R \sim e^{-\kappa_0 r}\, r^{p}$, with $p$ a
  frequency-dependent power. Nothing comes in from infinity, and below $\mu$
  nothing can go out.

These are conditions on the *homogeneous* (source-free) field, so they are
valid only where the source is negligible. That holds at both ends:

- Far out, the source falls like $e^{-(2\kappa_1+\kappa_2) r}$ with
  $2\kappa_1 + \kappa_2 \approx 0.25$, faster than the $m = 0$ field's own
  $e^{-0.12\, r}$.
- Near the horizon, the source carries a factor $\Delta$, which vanishes
  exponentially in $r_*$.

## 5. How it is solved numerically, from the ground up

Where each step's code is explained:

| Step | Function (`primary_secondary/dissipative_mode.jl`, 4a in `background.jl`) | Explained in |
|---|---|---|
| 4a | `mass_deficit_harmonic` | Section 3 |
| 4b | `dissipative_source`, `dissipative_tortoise_rhs!` | Section 4; 5.3; end of 5.8 |
| 4c | `dissipative_boundary_state` | 5.5, 5.6 |
| 4d | `dissipative_matching_residual` | 5.7 |
| 4d | `solve_dissipative_mode` | 5.8; Section 6 (flux, power identity) |
| 4e | driver along the branch | Section 7 |

### 5.1 The unknowns and their labels

The $m = 0$ field is described completely by three complex radial functions,
$R^0_0(r)$, $R^0_2(r)$, $R^0_4(r)$: one for each angular component $Y_{l0}$ in
the sum of Section 4. The code calls each pair $(l, m)$ a **mode**, so
`dissipative_modes = [(0, 0), (2, 0), (4, 0)]` means "the components with
$l = 0, 2, 4$, all with $m = 0$". The first number is always $l$, the degree of
the polar function $Y_{lm}(x)$ of Eq. (54); the second is $m$.

There is **no $n$** for this field. Two different $n$'s appear elsewhere, and
neither applies:

- Bound eigenmodes carry a hydrogenic label $n = l + 1 + (\text{number of radial nodes})$.
  The primary and secondary are the nodeless states for their $l$. The label
  exists because a bound state can only sit at discrete frequencies. The
  $m = 0$ field is not an eigenmode. Its frequency $\omega_0$ is imposed by the
  source, and its radial shape is whatever the source drives. There is no
  quantization and no node count.
- In the primary's code, `n = m ÷ fundamental_m` is the harmonic number $k$
  (frequency $k\omega_1$). The $m = 0$ field has a single frequency, so there
  is no harmonic number either.

### 5.2 The radial coordinate

The code integrates in the tortoise coordinate $r_*$, defined by
$dr_*/dr = (r^2+a^2)/\Delta$, with the additive constant fixed in `tortoise.jl`:

$$
r_* = r + \frac{2r_+}{r_+-r_-}\ln(r-r_+) - \frac{2r_-}{r_+-r_-}\ln(r-r_-) .
$$

Near the horizon $r_* \approx 8.09\,\ln(r-r_+)$ (at $a = 0.99$), so the horizon
is at $r_* = -\infty$. The domain and its landmarks:

| $r_*$ | $r$ | comment |
|---|---|---|
| $-100$ | $r_+ + 1.4\times10^{-6}$ | inner end of the domain |
| $5$ | $3.73$ | matching point (Section 5.7) |
| $150$ | $140.1$ | outer end of the domain |

Near the horizon a wave of frequency $\omega$ is simply $e^{\pm i\omega r_*}$,
with constant wavelength in $r_*$. A uniform step in $r_*$ resolves it
everywhere, whereas in $r$ the oscillations pile up against $r_+$.

### 5.3 From a second-order equation to a state vector

The radial equation is second order. RK4 (below) solves first-order systems,
so we introduce the slopes $Q_l \equiv dR_l/dr_*$ and write

$$
\frac{dR_l}{dr_*} = Q_l, \qquad
\frac{dQ_l}{dr_*} = \big(\text{the radial equation solved for } d^2R_l/dr_*^2\big).
$$

The second line is the radial equation of Section 4 rewritten in $r_*$
(`tortoise_operator_coefficients` in `core/radial_operator.jl` holds its
coefficients). At any $r_*$, the six complex numbers

$$
\text{state} = (R_0, R_2, R_4,\; Q_0, Q_2, Q_4)
$$

determine the solution from there on. The **right-hand side** function
(step 4b, `dissipative_tortoise_rhs!`) takes the state at one $r_*$ and returns
its derivative $d(\text{state})/dr_*$. In the code the secondary's six numbers
are placed in front (Section 5.8), making 12.

### 5.4 RK4 and its error

RK4 (fourth-order Runge–Kutta) advances the state from $r_*$ to $r_* + h$. It
evaluates the right-hand side four times: once at the start, twice at the
midpoint $r_* + h/2$, and once at the end. It then combines them with weights
$1, 2, 2, 1$ (divided by 6). Over a fixed distance its error is proportional to
$h^4$, so **halving the step divides the error by $2^4 = 16$**. Every check
below uses this. A discrepancy that shrinks by about 16 per halving is RK4
error. One that does not shrink is a mistake in the equations.

The right-hand side needs the primary at the start, the midpoint, and the end
of each step. The primary is saved only at spacing 0.1, so every field built on
it must use step 0.2. Then each midpoint lands on a saved sample.

### 5.5 How many solutions, and what the boundary conditions select

Without the source, each $l$ has two independent solutions (second order). With
three coupled $l$'s there are **six** independent source-free solutions, and
every source-free solution is a combination of them. With the source, the
general solution is one particular solution plus any combination of the six.

At each end of the domain, the six split into a physical and an unphysical
half:

| End | Physical (kept) | Unphysical (excluded) |
|---|---|---|
| horizon, $r_* \to -\infty$ | 3 ingoing, $\propto e^{-i\omega_0 r_*}$ | 3 outgoing, $\propto e^{+i\omega_0 r_*}$ (radiation from the horizon) |
| outer, $r \to \infty$ | 3 decaying, $\propto e^{-\kappa_0 r}$ | 3 growing, $\propto e^{+\kappa_0 r}$ |

Each end therefore leaves three free complex numbers (one per $l$): the
**amplitudes** $c_H = (c_{H,0}, c_{H,2}, c_{H,4})$ at the horizon and $c_O$ at
the outer end.

### 5.6 "The series": the boundary state

To start RK4 at $r_* = -100$ we need the actual numbers $(R_l, Q_l)$ of the
physical solutions there. These come from two series, one at each end. Each is
an asymptotic solution of the source-free equation, with the leading physical
behavior times correction terms.

**Horizon series** (Frobenius form), with $s \equiv r - r_+$:

$$
R_l = s^{p_H} \sum_{j=0}^{3} v_{j,l}\, s^j, \qquad
p_H = -i\omega_0\,\frac{r_+^2+a^2}{r_+-r_-} .
$$

Because $r_* \approx \frac{2r_+}{r_+-r_-}\ln s$ and $r_+^2 + a^2 = 2r_+$, the
factor $s^{p_H}$ equals $e^{-i\omega_0 r_*}$ up to a constant: ingoing. This is
Eq. (59) with $(\omega_1, m_1) \to (\omega_0, 0)$. The leading coefficients
$v_{0,l}$ are the amplitudes $c_{H,l}$, chosen freely. The corrections
$v_{1,l}, v_{2,l}, v_{3,l}$ follow by substituting the series into the
equation and setting each power of $s$ to zero. They mix the $l$'s through
$C_{ll'}$. At $r_* = -100$, $s \approx 1.4\times10^{-6}$, so the corrections are
tiny and order 3 is plenty.

**Outer series**:

$$
R_l = e^{-\kappa_0 r}\, r^{p_O} \sum_{j=0}^{7} \frac{w_{j,l}}{r^j}, \qquad
p_O = -1 + \frac{2\omega_0^2-\mu^2}{\kappa_0} .
$$

This is Eq. (58) with $\omega \to \omega_0$, $m = 0$, on the decaying branch.
The $-1$ in $p_O$ is the geometric $1/r$; the rest comes from the $1/r$
(Newtonian) part of the potential. Again $w_{0,l} = c_{O,l}$ is free, and the
rest follow power by power in $1/r$. At $r = 140$, $1/r \approx 0.007$, so terms
up to $1/r^7$ reach the $10^{-15}$ level.

The **boundary state** (step 4c, `dissipative_boundary_state`) evaluates one of
these series and its $r_*$-derivative at a given $r_*$ and amplitudes. It
returns the six numbers $(R_l, Q_l)$. The series themselves are the existing
`horizon_series` and `outer_series` (`core/boundary_series.jl`), shared
with the primary and the secondary. They solve the equation **without $N$**:
their recursion needs coefficients that are polynomials in $r$ (or $s$), and
$N$ is built from the numerical primary. That is valid only where $N$ is
negligible. Section 4 explains why the source is negligible at both ends.
Step 4c measures how much the background part matters.

### 5.7 Shooting from both ends and matching

1. Start at $r_* = -100$ from the horizon series with amplitudes $c_H$, and
   integrate outward to the **matching point** $r_* = 5$.
2. Start at $r_* = 150$ from the outer series with amplitudes $c_O$, and
   integrate inward to $r_* = 5$.
3. A true solution is smooth at $r_* = 5$. Both its value and its slope must
   agree there, for each $l$. The difference
   $(R^\text{horizon} - R^\text{outer},\; Q^\text{horizon} - Q^\text{outer})$ is
   the **matching residual**: six complex numbers, which must all vanish.

Six conditions for six unknown amplitudes. The integration runs from each end
toward the middle because that is the stable direction. Integrating the
decaying solution *outward* would let round-off seed the growing solution
$e^{+\kappa_0 r}$, which would soon swamp it. Integrated inward, the decaying
solution grows and any growing contamination shrinks. The matching point
$r_* = 5$ is the primary's choice, reused so that all three fields share one
domain.

### 5.8 Solving for the amplitudes: superposition, not Newton

For the secondary, $\omega_2$ was unknown and appears inside $B$ and $C$. The
residual depends on it nonlinearly, so we iterated: guess, integrate,
linearize, update. For the $m = 0$ field, $\omega_0$ is fixed by the two clouds,
and the only unknowns are the six amplitudes. The equation is linear in $R^0$
with a fixed source, so by superposition the residual is exactly

$$
\text{residual}(c) = r_0 + \mathsf{M}\, c, \qquad c = (c_H, c_O),
$$

where $r_0$ is the residual with all amplitudes zero (source alone), and
column $i$ of $\mathsf{M}$ is the extra residual when amplitude $i$ is set to 1.
RK4 is itself linear in the state, so this holds exactly for the discretized
equations too. Seven integrations give $r_0$ and $\mathsf{M}$. One $6\times 6$
complex solve, $\mathsf{M}c = -r_0$, gives the amplitudes, and a final
integration gives the profiles. Newton would get the same answer in one
iteration: its Jacobian is $\mathsf{M}$.

There is no normalization choice: the source sets the scale, so $R^0 \propto A_2$
and fluxes scale as $|A_2|^2$.

**What could go wrong: near-resonance.** $\mathsf{M}$ is singular only if a
source-free solution satisfies both boundary conditions at the real frequency
$\omega_0$, that is, if $\omega_0$ were exactly an $m = 0$ bound-state
frequency. Those states are damped (complex frequency), so exact singularity
cannot happen. A nearby level still enhances the response. Rough hydrogenic
estimates, $\omega \approx \mu(1 - \alpha^2/2n^2)$ with $n$ the hydrogenic label
of Section 5.1, put even-$l$, $m = 0$ levels near $0.383$ ($n = 1$) and $0.411$
($n = 2$). $\omega_0$ is $0.402$ at small $A_1$ and falls to $0.348$ at
$A_1 \approx 10$, so it crosses that region along the branch. The paper's
Eq. (38) has the same physics in its factor $1/W_{l0}(\omega_0)$. Section 8
(step 4e) reports what the runs found: no resonance at $M\mu = 0.42$; at
$M\mu = 0.3$ a resonance near $0.296$ that $\omega_0$ sits just below.

**Why the secondary is integrated alongside.** The source needs the
secondary's $R^2_l$ at every RK4 stage (Section 5.4). The secondary was saved
only at its own step points (spacing 0.2), and the midpoints fall between them.
Doubling the $m = 0$ step to 0.4 would put the midpoints on saved secondary
samples, but $105/0.4$ is not an integer, so the steps would not land on the
matching point. Instead, the secondary and the $m = 0$ field are integrated
together as one 12-number state $[R^2, Q^2;\, R^0, Q^0]$ at step 0.2. The
secondary is then computed at every stage, and the primary is read at saved
samples.

## 6. What is measured: the horizon flux

The outward energy flux through a sphere at tortoise radius $r_*$ is, for a
single-frequency field written as above (`outward_secondary_flux` in
`primary_secondary/secondary_energy.jl`),

$$
F(r_*) = 4\pi\,\omega_0\,(r^2+a^2)\sum_l \mathrm{Im}\!\left(\overline{R_l}\,Q_l\right).
$$

The factor $4\pi$ comes from the time average of the $e^{\pm i\omega_0 t}$
cross terms (a factor 2) and the azimuthal integral ($2\pi$). The energy lost
into the horizon is $F_H^{(0)} = -F$ at the inner end. The flux at the outer
end is tiny, because the field decays there.

**An exact identity used as a check (the power identity).** Multiply the radial
equation by $\overline{R_l}$, sum over $l$, and take the imaginary part. Then:

- The $B$, $C$, and background terms drop out: they are real (with $\omega_0$
  real) and symmetric in $l, l'$.
- The derivative term becomes a total derivative,
  $\tfrac{d}{dr}\,\mathrm{Im}\sum_l \Delta\,\overline{R_l}\,R_l'$.
- With $\Delta\, dR/dr = (r^2+a^2)\,Q$, this is $F/(4\pi\omega_0)$.

Integrating from the horizon outward gives

$$
F_H^{(0)} + F_\text{out}^{(0)} = -4\pi\,\omega_0 \int_{r_+}^{\infty} dr\;
\mathrm{Im}\sum_l \overline{R^0_l}\, S_l ,
$$

where $S_l$ is the source part of $N_l$. In words: the energy leaving through
the boundaries equals the work done per unit time by the source, which is the
energy taken from the two clouds. The numerical solution must satisfy this. It
is the analogue of the $2\Gamma_2 E_2 = -F$ check of piece 3.

## 7. Steps

| Step | Content | Check |
|---|---|---|
| 4a | `mass_deficit_harmonic`: $[1-\cos\phi_1]_k$ | $k = 0$ equals $\langle 1-\cos\phi_1\rangle$; $k = 2 \to \tfrac12 P_1^2$ with error $\propto A_1^2$; angular grid doubled |
| 4b | right-hand side for the combined state $[R^2,Q^2;\,R^0,Q^0]$ | with the source off, $(r^2+a^2)\sum_l\mathrm{Im}(\overline{R^0_l}Q^0_l)$ is constant in $r_*$ (Section 6 with $S = 0$) |
| 4c | $m = 0$ boundary data (the generic helper rejects $m = 0$) | homogeneous solution integrated from the series at $r_* = -100$ matches the series at $r_* = -90$; same at the outer end |
| 4d | the linear solve for the amplitudes | residual at round-off; power identity |
| 4e | driver along the saved branches ($M\mu = 0.42$ and $0.3$) | $F_H^{(0)} \propto A_1^4$ at small $A_1$; step convergence at $A_1 = 2.05$; adding $l = 6$; $l = 0$ share |

Optional, after 4e: compute the perturbative $F_0$ of Eq. (38) independently
with a Green's function, as a second method.

## 8. Results so far (2026-10-02)

All three checks run on the saved $M\mu = 0.42$ cosine branch and its averaged
secondary. Scripts: `tests/adiabatic/primary_secondary/check_mass_deficit_harmonic.jl`,
`check_dissipative_current.jl`, `check_dissipative_boundary.jl`.

### Step 4a: the mass-deficit harmonic

Evaluated at $r_* = -20, 5, 20, 50$:

| Check | Result |
|---|---|
| $k = 0$ against `averaged_mass_deficit` ($A_1 = 0.0011, 9.9, 13$) | identical |
| $k = 2$ against $\tfrac12 P_1^2$ at $r_* = 20$ | relative difference $/A_1^2 = 0.441$, constant as $A_1$ doubles from $0.0011$ to $0.009$, then $0.437, 0.425$ |
| azimuthal grid 31 → 62 points | relative change $\le 6\times10^{-16}$, even at $A_1 = 13$ |

The "relative difference" in the second row is the gap between the exact harmonic
and its small-field formula, not a numerical error: it is the nonlinear
correction the calculation is meant to keep. At $A_1 \approx 10$, where
$|\phi_1|$ reaches about 3 near $r_* = 5$, the small-field formula is too large
by about a factor of 2. The constant ratio at small $A_1$ shows the leading term
is exactly the perturbative source, with the right phase convention. The grid
test shows no aliasing, even where $|[1-\cos\phi_1]_2| \approx 0.47$. Neither
check tests the saved primary itself: its $l \le 5$ truncation and its lack of
step convergence above $A_1 \approx 9$ remain open.

### Step 4b: the right-hand side conserves the current

**Idea.** With the source off, the radial current
$J \equiv (r^2+a^2)\sum_l \mathrm{Im}(\overline{R_l}\,Q_l)$ is exactly constant
in $r_*$ for *any* solution (Section 6 with $S = 0$). The identity needs no
boundary conditions, so it tests the right-hand side alone: the Kerr
coefficients at $m = 0$, the $l \leftrightarrow l\pm2$ coupling, and the
background term.

**Procedure.**
1. Choose $\omega_0$ for the case at hand.
2. Choose six random complex numbers as the state at $r_* = 5$. Set the
   secondary to zero and switch the source off.
3. Integrate outward to $r_* = 61$ and inward to $r_* = -59$ (intervals that
   are multiples of 0.4) at steps 0.4, 0.2 and 0.1.
4. At every step compute $J$. Report the largest drift $|J - J(5)|$, divided by
   $(r^2+a^2)\sum_l |R_l||Q_l|$ at the same point.

That denominator is the size of the terms whose imaginary parts make up $J$.
Dividing by $J(5)$ instead is misleading outward: the growing solution makes
$|R|$ huge while $J$ stays small.

| Case | Direction | Step 0.4 | Step 0.2 | Step 0.1 | Ratio per halving |
|---|---|---|---|---|---|
| $A_1 = 0$ (no background) | out | $2.3\times10^{-5}$ | $1.2\times10^{-6}$ | $7.0\times10^{-8}$ | 18.8, 17.5 |
| $A_1 = 0$ | in | $9.3\times10^{-6}$ | $7.6\times10^{-7}$ | $5.3\times10^{-8}$ | 12.3, 14.3 |
| $A_1 = 9.9$ | out / in | $2.3\times10^{-5}$ / $9.3\times10^{-6}$ | $1.2\times10^{-6}$ / $7.5\times10^{-7}$ | — | 19.0 / 12.3 |
| $A_1 = 13$ | out / in | $2.2\times10^{-5}$ / $9.2\times10^{-6}$ | $1.2\times10^{-6}$ / $7.5\times10^{-7}$ | — | 19.0 / 12.4 |

The drift is RK4 error (ratios near 16). A sign or symmetry mistake in the
operator would leave a drift that does not shrink. The saved primary cannot be
read below step 0.2, so the large-$A_1$ rows stop there.

### Step 4c: the boundary series agree with the equation

**Idea.** The series (Section 5.6) claims to be a solution of the source-free
equation near each end. Start a numerical solution from the series at the
domain end, integrate it a short way in, and compare with the series evaluated
directly at the new point. If the series is right, the two differ only by RK4
error.

**Procedure** at the horizon end (the outer end is the same, from $r_* = 150$
inward to $140$):
1. Choose $\omega_0$ and amplitudes $c_H = (1, 0, 0)$: unit amplitude in
   $l = 0$ only.
2. Evaluate the horizon series at $r_* = -100$. This gives the starting state.
3. Integrate the source-free equation (secondary zero) to $r_* = -90$.
4. Evaluate the horizon series directly at $r_* = -90$ with the same amplitudes.
5. Report $|\text{integrated} - \text{series}| / |\text{series}|$.
6. Repeat with $c_H = (0,1,0)$ and $(0,0,1)$. These are the three numbers in
   each cell, below given as a range.

At large $A_1$, a second run at the same $\omega_0$ but with the primary set to
zero separates the background's effect from the rest.

| Case | Step 0.4 | Step 0.2 | Step 0.1 |
|---|---|---|---|
| $A_1 = 0$, horizon | $2.2\times10^{-5}$ | $1.4\times10^{-6}$ | $8.8\times10^{-8}$ |
| $A_1 = 0$, outer | $5\times10^{-8}$ | $3\times10^{-9}$ | $2$–$6\times10^{-10}$ |
| $A_1 = 9.9$, horizon | $1.1\times10^{-5}$ | $6.7\times10^{-7}$ (no background: $6.8\times10^{-7}$) | — |
| $A_1 = 9.9$, outer | $1.5\times10^{-6}$ | $1.1$–$1.3\times10^{-7}$ (no background: $1.0\times10^{-7}$) | — |
| $A_1 = 13$, horizon | $9.8\times10^{-6}$ | $6.0\times10^{-7}$ (no background: $6.1\times10^{-7}$) | — |
| $A_1 = 13$, outer | $2.0\times10^{-6}$ | $2.0$–$2.7\times10^{-7}$ (no background: $1.2\times10^{-7}$) | — |

- At $A_1 = 0$ the difference falls by 16 per halving. So both series satisfy
  the equation, including the $l$ coupling, at $m = 0$.
- At the horizon the background makes no difference: it carries the factor
  $\Delta \approx 4\times10^{-7}$ at $r_* = -100$.
- At the outer end the background adds about $1.5\times10^{-7}$ at $A_1 = 13$,
  which is also why the step ratio there is 7–12 rather than 16. The source is
  **not** the primary's bound part, which has decayed to $|R_{11}| \sim 10^{-12}$
  there. It is the primary's **third harmonic**: the $(3,3)$ and $(5,3)$
  components at frequency $3\omega_1 \approx 1.10 > \mu$. These are travelling
  waves, the radiation the nonlinear cloud emits to infinity. Their amplitudes
  ($|R_{33}| \approx 1.8\times10^{-4}$, $|R_{53}| \approx 2.4\times10^{-4}$)
  barely change between $r_* = 140$ and $150$. They keep
  $\langle 1-\cos\phi_1\rangle \approx 10^{-7}$ at the outer boundary, which
  the source-free series omits. That lowers the effective $\mu^2$ by
  $\alpha^2\times10^{-7}$, shifts $\kappa_0$ by about $4\times10^{-8}$, and over
  $10M$ gives a mismatch of a few $\times10^{-7}$, as seen.
- This is not worth correcting. It amounts to starting the outer solution with
  a decay rate off by parts in $10^7$, where $R^0$ has already decayed by
  $e^{-\kappa_0 \cdot 140}$. The horizon flux is set near the cloud. The
  primary's and secondary's outer boundaries carry the same approximation.
  Outer-radius convergence in step 4e would expose it if it mattered.

### Step 4d: the solve and the power identity

**Procedure** at each saved point:
1. Take $\omega_0 = 2\,\mathrm{Re}\,\omega_1 - \mathrm{Re}\,\omega_2$.
2. Run `solve_dissipative_mode` (Section 5.8).
3. Compare $F_H + F_\text{out}$ with the work done by the source,
   $-4\pi\omega_0\int dr\,\mathrm{Im}\sum_l \overline{R^0_l}S_l$ (Section 6),
   integrated over the saved profiles with Simpson's rule.

Script: `tests/adiabatic/primary_secondary/check_dissipative_solve.jl`.

| $A_1$ | $\omega_0$ | cond (scaled) | residual after / before | $F_H/\hat E_2$ | power ratio | $l = 2$ share | $l = 4$ share |
|---|---|---|---|---|---|---|---|
| 0.0011 | 0.40198 | 4.4 | $2.0\times10^{-15}$ | $5.95\times10^{-14}$ | 1.00000026 | $2.2\times10^{-6}$ | $9\times10^{-14}$ |
| 0.344 | 0.39761 | 126 | $1.2\times10^{-15}$ | $9.41\times10^{-5}$ | 1.00000021 | $1.4\times10^{-6}$ | $5\times10^{-13}$ |
| 9.14 | 0.35033 | 81 | $1.5\times10^{-15}$ | $3.29\times10^{-2}$ | 0.99999921 | $2.3\times10^{-5}$ | $2\times10^{-10}$ |
| 9.94 | 0.34807 | 78 | $3.4\times10^{-15}$ | $3.47\times10^{-2}$ | 0.99999916 | $2.5\times10^{-5}$ | $3\times10^{-10}$ |
| 10.2 | 0.34729 | 77 | $3.1\times10^{-15}$ | $3.53\times10^{-2}$ | 0.99999915 | $2.6\times10^{-5}$ | $3\times10^{-10}$ |
| 13.0 | 0.34064 | 71 | $2.5\times10^{-15}$ | $3.99\times10^{-2}$ | 0.99999903 | $3.4\times10^{-5}$ | $4\times10^{-10}$ |

- The residual after the solve is at round-off, as superposition requires.
- The power identity holds to $10^{-6}$. This is the first check that involves
  the source: the boundary fluxes equal the energy the source supplies.
- The raw condition number of $\mathsf{M}$ is about $5\times10^5$. That is
  column scale only: a unit outer amplitude grows enormously on its way in to
  $r_* = 5$. After scaling each column to unit norm it is 4–126, so there is no
  near-resonance at these points.
- $l = 0$ carries the horizon flux up to parts in $10^5$, as the paper says.
  The higher $l$ are held back by the angular-momentum barrier at
  $\omega_0 r_+ \approx 0.4$.
- $F_\text{out} \lesssim 10^{-25}$: all the loss goes into the horizon.
- $F_H/\hat E_2$ (energy into the horizon per unit time, per unit secondary
  energy) rises from $6\times10^{-14}$ to $0.04$ along the branch.

Not yet established: convergence in step, outer radius, and $l$-truncation,
and agreement with the perturbative $F_0$ (step 4e and the optional check).

### Step 4e: along both branches

The driver `src/adiabatic/primary_secondary/drivers/dissipative_fluxes.jl` runs steps 4a–4d at
every saved point. Commands, the full results table, and file names are in
the [campaign README](../runs/adiabatic/multimode_secondary/README.md), section
"Piece 4". In brief:

- **Every point passes the step 4d checks.** The residual after the solve is
  $\le 8\times10^{-15}$, and the power identity holds to $10^{-6}$
  ($M\mu = 0.42$) and $8\times10^{-8}$ ($M\mu = 0.3$).
- **$F_H/\hat E_2$ rises monotonically** along both branches. At small $A_1$ its
  log-slope against $A_1$ is $3.9997$: the $E_1^2$ law of the perturbative
  $F_0E_1^2E_2$. It bends below 4 as the nonlinearity sets in.
- **Convergence.**
  - Radial step at $A_1 = 2.05$: changes of $4\times10^{-8}$ relative, at the
    Newton-tolerance floor.
  - $l$-truncation: adding $l = 6$ changes $F_H/\hat E_2$ by up to
    $4.5\times10^{-4}$.

**The condition number at $M\mu = 0.3$.** The scaled condition number of
$\mathsf{M}$ grows from 12 to 4100 along the $M\mu = 0.3$ branch, while
$\omega_0$ barely moves. Two cheap tests settled what it means
(`tests/adiabatic/primary_secondary/check_dissipative_resonance.jl`):

1. **Scan the drive frequency** at fixed primary, secondary, and source. The
   $m = 0$ response has a real resonance near $\omega \approx 0.295$–$0.296$,
   with a sharp dip just above. That fits the hydrogenic $m = 0$ level at
   $0.2966$. $\omega_0$ sits on its lower flank at both $A_1 = 0.001$ and
   $0.06$, so the enhancement exists already in the linear limit. But
   $\sigma_\text{min}/\sigma_\text{max}$ does not dip at the resonance: it
   falls smoothly across the whole band.
2. **Find which amplitude combination is nearly in the null space** of
   $\mathsf{M}$ (the right singular vector of $\sigma_\text{min}$). At small
   $A_1$ it is horizon $l = 0$ plus outer $l = 0$: the mild resonance. At
   $A_1 = 0.06$ it lies entirely on the outer $l = 2$ and $l = 4$ amplitudes.
   A true near-resonance needs horizon and outer parts together.

The large condition number is therefore a shooting artifact. Integrated
inward through the cloud, where the background mixes even $l$, the outer
$l = 2$ and $l = 4$ solutions line up with each other. At $4\times10^3$ it
costs under 4 of 16 digits. The residual and the power identity confirm the
solution. The driver warns above $10^8$; the remedy would be re-orthogonalizing
the outer solutions during the inward integration.

## 9. How this compares with what the paper did

The paper says only that Eq. (51) is solved "in a method similar to the one
used in determining $\phi_2(A_1)$" (p. 11), with real frequencies (footnote 2)
and the flux of Eq. (66).

- **Same:** projected coupled radial equations, boundary conditions like
  (58)–(59) at $(\omega_0, 0)$, the primary's background in the operator, real
  $\omega_0$, and the flux of Eq. (66).
- **Our choices where the paper is silent:** the averaged reading of the
  background; $l = 0, 2, 4$; and the secondary normalized to unit outer
  amplitude instead of unit energy (a change of variables, converted with
  $\hat E_2$).
- **Implementation only:** the superposition solve and the combined
  integration. They give the same solution as the paper's shooting.
- **Ours only:** all the checks.

## Symbols introduced here

| Symbol | Meaning |
|---|---|
| $\psi_1, \psi_2$ | cloud phases $\omega_1 t - \varphi$, $\omega_2 t - 2\varphi$ |
| $\omega_0$ | $2\,\mathrm{Re}\,\omega_1 - \mathrm{Re}\,\omega_2$, frequency of the $m = 0$ field |
| $[1-\cos\phi_1]_k$ | coefficient of $e^{-ik\psi_1}$ in $1 - \cos\phi_1$ (code: `mass_deficit_harmonic`) |
| $P_0, P_1, P_2$ | angular profiles $\sum_l R_l Y_{lm}$ of the $m = 0$ field, primary ($m = 1$ part), secondary |
| $\kappa$ | $\sqrt{\mu^2-\omega^2}$, decay rate at large $r$ |
| $Q_l$ | $dR_l/dr_*$, the slope carried in the state vector |
| $s$ | $r - r_+$, distance from the horizon |
| $p_H, p_O$ | powers in the horizon and outer series (Section 5.6) |
| $c_H, c_O$ | free leading amplitudes of the horizon and outer series, one per $l$ |
| $J$ | radial current $(r^2+a^2)\sum_l \mathrm{Im}(\overline{R_l}Q_l)$, equal to $F/(4\pi\omega_0)$ |
| $\mathsf{M}, r_0$ | superposition matrix and source-only matching residual |
| $F_H^{(0)}$ | energy per unit time into the horizon carried by the $m = 0$ field |
