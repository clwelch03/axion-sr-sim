# Strategy checkpoint: separate fast oscillations from slow cloud evolution

2026-09-16. Research and feasibility checkpoint; no replacement solver or
orders-of-magnitude nonlinear speedup has been demonstrated.

The recommendation is to build a relativistic nonlinear profile continuation
method for long, slowly evolving cloud histories, retaining the time-domain
towers for transients and validation. This changes the time representation,
while preserving radial profiles and angular harmonics. A cubic-only source
optimization is not the recommended main route.

## Evidence from this repository

A fresh single-thread run of `src/benchmark_nonlinear.py --repeats 3`, at
a=0.99, alpha=0.42, dx=M/6, 4,267 radial nodes and six channels, gave:

| Operation | Median | Range |
|---|---:|---:|
| NumPy nonlinear source | 0.741381 s | 0.711357–0.922874 s |
| Linear coupled RK4 | 0.011893 s | 0.011562–0.012221 s |
| NumPy nonlinear coupled RK4 | 2.758012 s | 2.681483–2.806181 s |

At dt=0.8 dx, 1e5 M needs 750,000 steps. These timings imply 2.48 h
for the linear calculation and 574.59 h for the nonlinear calculation.
They exclude I/O and are repeated steps from a fixed Gaussian state.
Eliminating the ENTIRE nonlinear source would remove much more work than
eliminating just Phi-sin(Phi), but the fine radial timestep remains.
For 1e6 M at dx=M/8 there are ten million steps. A 12-hour target would
require 4.32 ms per full step, including output amortized over the run.
No existing benchmark establishes that performance.

`solve_qbs` on the same radial grid and m=1, l_max=5, with
omega_R_guess=0.409 and maxiter=500 returned

    M omega = 0.408805564684 + 1.50613981e-7 i
    solve time = 0.440 s
    oscillation period = 15.370 M
    amplitude e-folding time = 6.63949e6 M

This is a linear feasibility measurement, not a nonlinear timing prediction
or a new convergence study. It demonstrates that our existing operator can
extract slow growth without evolving millions of periods. In its linear
limit a 1e6 M run grows the amplitude only by exp(0.150614), about 16%.

## Established alternative and its scope

[Omiya, Takahashi & Tanaka (2022)](https://arxiv.org/html/2201.04382)
solve relativistic nonlinear radial profiles along an amplitude sequence,
retaining the cosine potential, then determine slow evolution from energy
fluxes. Appendix A uses the same six scalar channels as our benchmark.
The method assumes a single helical pattern; it does not resolve collapse
or arbitrary independent frequencies. Appendix B finds a roughly 30% flux
change at large amplitude on adding modes: six channels are a starting
point, not a precision guarantee.

[Omiya et al. (2023)](https://arxiv.org/html/2211.01949) extends the treatment
to a secondary cloud and finds that mode interaction can determine
saturation. Independent clouds cannot be represented merely by adding the
phase-locked m=3 harmonic of an m=1 cloud. Our current positive-odd-m
restriction also excludes a separately seeded m=2 cloud.

## Proposed mapping to our variables

For one rotating pattern, use the existing horizon-adapted coordinates and
write a candidate family as

\[
\psi_{\ell m}(t,x)\simeq q_{\ell m}(x;s(t))e^{-im\vartheta(t)},
\qquad \dot\vartheta=\Omega(s).
\]

Here s labels the profile family, not necessarily the peak amplitude. Record
max|Phi| as an observable so the requested 0.01, 0.1, 0.3, 0.7 cases retain
their meaning. Profiles may change shape with s. This is different from
multiplying one fixed eigenvector by a changing number.

At frozen s, time derivatives become -im Omega and -m^2 Omega^2. Define
L_m q as the spatial and mass terms of Eq. 25, and C_m pi as its term
linear in pi. The leading profile residual in the physical interior is

\[
\mathcal R_m[q,\Omega]
=-m^2\Omega^2 A_m q_m-L_mq_m
 +im\Omega C_mq_m-b_{{\rm NL},m}[q].
\]

The source is the existing full sine source evaluated at phase zero.
Its phase law is exact in the continuum:

\[
b_{{\rm NL},m}[\{q_{\ell k}e^{-ik\vartheta}\}]
=e^{-im\vartheta}b_{{\rm NL},m}[\{q_{\ell k}\}].
\]

This follows by changing the azimuthal integration variable. A numerical
check with all six complex channels, sampled peak |Phi|=3, 64x128 angular
points, r=7, and phases 0.37, 1.1, 2.4 gave relative source errors
3.59e-16, 3.91e-16, 3.83e-16. This verifies the nonlinear phase law;
it does not establish a solution family or its stability.

There is an essential missing part before implementing a nonlinear solve:
an exactly periodic real-frequency solution with physical radiative
boundaries cannot have unbalanced net energy flux. Slow growth terms, a
solvability condition, or the published complex-frequency prescription
must complete the frozen residual. Simply setting this residual to zero
with real Omega and arbitrary amplitude would generally overconstrain it.
Nor does assigning a common exponential growth to every harmonic give an
exact nonlinear solution. Phase and amplitude conditions also need to be
specified to remove normalization/phase degeneracies.

For a fixed background, the slow balance would be

\[
E_s\dot s=-F_E^H-F_E^\infty,
\]

with positive horizon flux INTO the hole and positive outer flux OUT of
the cloud. Define E on a consistent finite cloud region or subtract the
outgoing radiation contribution; an eternal outgoing wavetrain is not a
finite-energy cloud. Converge extraction radius and exterior matching.
The growth observable for comparisons can be
Gamma_E=(1/2)d ln(E)/dt, reducing to omega_I in the linear limit. At finite
amplitude it need not equal d ln(max|Phi|)/dt. An alpha–growth-rate plot
therefore needs amplitude and rate-convention labels.

## Keeping a route to self-gravity and backreaction

These are proposed interfaces, not implemented gravity:

- Keep radial profiles, angular mode labels, physical normalization f_a,
  and derivatives accessible. A table containing only Gamma(alpha,A)
  would discard information needed for stress-energy and self-gravity.
- Give the profile residual an explicit geometry input. A weak gravitational
  field changes derivative terms and the acceleration matrix as well as
  the algebraic source. The current `additional_b` hook is insufficient
  for a general evolving metric.
- Add cloud-sourced gravitational multipoles as separate unknowns, beginning
  with a justified weak-field or perturbative approximation. They need even
  azimuthal numbers including zero, unlike our present scalar towers.
  Solve for them together with the scalar profile, retaining a future
  time-dependent metric path. Poisson gravity is an approximation whose
  validity is not guaranteed solely by small cloud mass at alpha=0.42.
- Evolve BH mass and angular momentum from physical horizon fluxes. Fix
  reference units M0; keep mu constant as M_BH changes. Do not mutate spin
  in the cached Tower or silently change its tortoise grid.
- With evolving geometry, derive consistent total energy/angular momentum
  accounting including binding energy and work terms. Even a tentative
  family energy E(s,M,J) gives E_s sdot + E_M Mdot + E_J Jdot, not merely
  E_s sdot. A fixed-background Killing-energy balance cannot automatically
  be reused in a time-dependent metric.

Self-gravitating profile methods are independently supported by
[May, East & Siemonsen](https://arxiv.org/abs/2410.21442), who construct
general-relativistic cloud solutions under an axisymmetry/oscillation-average
assumption. This supports compatibility of the approach with gravity, not
equivalence to a dynamical real sine-Gordon Einstein system.
[Takahashi](https://arxiv.org/abs/2604.21239) provides relativistic perturbative
frequency-shift benchmarks for self-interactions and self-gravity.
See [the gravity review](gravitational-backreaction-review.md) for the
coordinate, normalization, boundary and flux issues already identified.

## Decision gates before a replacement or speed claim

1. Derive the slow-growth closure, normalization and physical boundary
   conditions in our coordinates. Recover the linear QBS and the relation
   -(F_H+F_infinity)/(2E)=omega_I. The current outer frozen boundary is
   unsuitable for nonlinear harmonics with m Omega > mu.
2. Continue to one small finite-amplitude profile with the full sine source.
   Demonstrate residual convergence, radial/angular convergence, and the
   expected leading frequency correction proportional to amplitude squared.
3. Compare the profile with a short time-domain evolution: drift, generated
   harmonics and fluxes. Add an outgoing outer treatment before claiming
   long-term radiative balance. Stability against off-family perturbations
   is a separate check from convergence of a nonlinear root solve.
4. Measure the total cost of constructing a converged amplitude sequence,
   including flux extraction, then compare with the corresponding evolution.
   A reduced number of time steps alone is not a measured speedup.
5. Use time-domain windows when the family changes rapidly, loses stability,
   or independent modes matter. If unrestricted dynamics are required
   throughout, this approximation cannot replace the time-domain solver;
   a compiled full RHS and radial-grid/time-integration redesign remain
   candidates requiring their own benchmarks.

A direct cubic contraction remains a useful controlled weak-field comparison.
But at Phi=0.7 the cubic force differs from Phi-sin(Phi) by 2.48%, and at
Phi=3 by 57.4%. These pointwise errors are not projected-source error bounds.
Replacing the sine potential with the cubic term would change the question
being answered unless its error is demonstrably acceptable throughout the run.

## Slow-growth closure checkpoint

The following is a derivation from the repository's first-order evolution,
not a claim to have implemented the published nonlinear boundary solver.
It gives an explicit consistency residual for a proposed profile family.

Let F(U) denote `coupled_tower_rhs`, with U including psi, pi and the PML
auxiliary chi. Define the azimuthal rotation R(theta) on every field of
tower m by multiplication by exp(-i m theta). For a differentiable family
U(s), substitute u(t)=R(theta(t))U(s(t)) into the full first-order system.
Rotational covariance gives

\[
 f(s)\,\partial_s U_m=F_m(U)+im\Omega(s)U_m,
 \qquad \dot s=f(s),\quad\dot\theta=\Omega(s).
\]

This is an exact invariance equation for such a family in the semidiscrete
system (up to angular quadrature error). Existence of a useful one-parameter
family is an additional assumption: it does not represent arbitrary initial
data, independent frequencies or radiation history. Also, the repository's
outer boundary remains reflective. Algebraic consistency with this F alone
does not establish physical radiative balance.

For a concrete normalization, choose a fixed complex linear functional L
supported on the m=1 psi components and require L(psi_1)=s, real and positive.
For example, use the discrete overlap with a reference eigenprofile, divided
by its squared norm. Applying L to the psi row gives

\[
 f=\operatorname{Re}L(\pi_1),\qquad
 \Omega=-\operatorname{Im}L(\pi_1)/s.
\]

The phase condition becomes singular at s=0; initialize from a nonzero
linear seed and change the normalization if its overlap vanishes. s is not
automatically peak |Phi|. Crucially, momentum satisfies

\[
\pi_m=f\,\partial_s\psi_m-im\Omega\psi_m.
\]

Setting pi=-im Omega psi would discard the growth term. Keeping all three
state fields in the family avoids introducing second derivatives with respect
to s when substituting directly into the second-order wave equation.

For the linear m=1 eigenstate V with eigenvalue lambda=Gamma-i omega_R,
U_1=sV, f=Gamma s, Omega=omega_R solves this equation. Reconstruct the PML
component of the eigenstate from the returned psi using

\[
\pi=\lambda\psi,\qquad
\chi=\frac{\gamma(\partial_x\psi+i\Omega_H\psi)}
                 {\lambda+\gamma+i\Omega_H}.
\]

A fresh check at the parameters and grid above compared this reconstructed
state against the actual `Tower.rhs` (including boundary freezing). It gave
omega=0.408805564684+1.50613981e-7 i and

    ||F_linear(V)-lambda V|| / ||V|| = 6.46e-11
    ||F_linear(V)-lambda V|| / ||Gamma V|| = 4.29e-4

The second normalization is important: a small error relative to the carrier
frequency alone need not resolve the much smaller growth rate. This is an
operator consistency check, not an energy-flux validation.

Using that same eigenprofile scaled to sampled peak A=max|Phi|, retaining
the full sine force on 32x64 angular points, gives the following residual
if the profile shape and generated modes are held fixed:

| A | Closure residual / state norm | Norm of missing m=3 acceleration |
|---|---:|---:|
| 0.01 | 7.6703e-7 | 5.9879e-7 |
| 0.1 | 7.6681e-5 | 5.9855e-4 |
| 0.3 | 6.8852e-4 | 1.6107e-2 |
| 0.7 | 3.7051e-3 | 2.0123e-1 |

These are Euclidean norms of the stored arrays in M=1 units, not energy
errors or estimates of growth-rate shifts. The small-amplitude residual
divided by state norm scales as A^2; the generated m=3 acceleration scales
as A^3. A fixed-template amplitude model therefore fails this consistency
check. A deformed profile family with additional harmonics is needed.

Solving the invariance equation by implicit continuation/collocation in s
is a candidate alternative to the published frozen complex-frequency
construction. It is not yet a validated algorithm. Multiplying existing
time-domain variables by a rotating phase alone does NOT remove the radial
CFL restriction or demonstrate a speedup. Boundary conditions, continuation
conditioning, family stability and physical flux checks remain open.

The same derivation permits further slow coordinates z=(s,M_BH,J_BH,...):
replace f U_s by sum_a zdot_a U_{z_a}. Metric multipoles transform with their
own azimuthal indices (including zero and even indices). This supplies an
explicit mathematical extension for backreaction without assuming it is
captured by one scalar forcing function. It still requires a derived matter–
gravity evolution system and consistent constraints; rapid gravity dynamics
cannot be forced into slow coordinates by this change of variables.
