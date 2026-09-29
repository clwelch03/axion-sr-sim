# Multimode adiabatic evolution: initial design

2026-09-28. Proposed framework and first checkpoint; no multimode solver is
implemented by this note. Start on a fixed Kerr background.

See the [September 28 literature review](multimode-literature-review-2026-09-28.md)
for prior work and the narrower novelty claim: mutual nonlinear profile
response and its effect on evolution. General multimode rate equations and
perturbative cross-mode frequency shifts already exist.

## Reproduction target

[Omiya et al., Section III](https://arxiv.org/html/2211.01949v1#S3)
uses a nonlinear primary cloud, a secondary cloud linearized on that background,
and driven dissipative fields. Equations (50)–(54) give

\[
\phi_1=\phi_1(A_1),\qquad
\phi_2=A_2\psi_2(A_1),\qquad
(\square_g-\mu^2\cos\phi_1)\psi_2=0.
\]

The retained dissipative channels have
\((\omega,m)=(2\omega_{1,R}-\omega_{2,R},0)\) and
\((2\omega_{2,R}-\omega_{1,R},3)\).
Equations (69)–(70) determine two amplitude rates through total energy and
angular momentum balance. This approximation neglects deformation of the
primary by the secondary. It is the initial benchmark, not a solution for
arbitrarily many comparably nonlinear clouds.

## What can be reused

The existing stationary branch supplies the primary profiles. Low-level radial
coefficients already accept an explicit frequency and azimuthal number.
Matching, asymptotic series, quadrature, and implicit observable derivatives
provide useful components.

However, the assembled stationary and observable routines assign frequencies
as `frequency * m / fundamental_m`. Their `(l,m)` entries are harmonics of one
pattern, not independent bound states. Adding an entry to that list cannot
create a second independently oscillating cloud. The current scalar amplitude
integrator likewise does not constitute a multimode evolution solver.

Keep three identities separate:

- A bound cloud is identified by its linear-limit `(n,l,m)` and continued
  eigenbranch; `n` here is the principal quantum number.
- A Fourier component is identified by an integer vector
  \(\boldsymbol q\), with
  \(\Omega_{\boldsymbol q}=\boldsymbol q\cdot\boldsymbol\omega_R\) and
  \(m_{\boldsymbol q}=\boldsymbol q\cdot\boldsymbol m\).
- Each component has its own angular expansion in \(l\). These angular entries
  do not acquire independent cloud amplitudes.

Store one representative of each real-field conjugate pair. Allow zero and
negative azimuthal numbers with explicit conjugation conventions. In
particular, the absorbing channel requires `m=0`. Distinct frequencies can
share the same `m`, especially when radial overtones are included.

## Why more clouds need additional equations

For any proposed family parameterized by \(K\) real amplitudes,

\[
\begin{pmatrix}
\partial_{A_1}E&\cdots&\partial_{A_K}E\\
\partial_{A_1}J&\cdots&\partial_{A_K}J
\end{pmatrix}\dot{\boldsymbol A}
=-\begin{pmatrix}F_E\\F_J\end{pmatrix}.
\]

This matrix has rank at most two. For five amplitudes, at least three directions
of motion are undetermined by these balances alone. Even for two amplitudes,
check invertibility and conditioning. Conservation is essential but cannot
supply the missing mode-transfer dynamics.

A useful weak-interaction starting point is to enumerate processes
\(i+j\to k+\text{emitted field}\). Associate each with

\[
\nu_{ac}=\delta_{ak}-\delta_{ai}-\delta_{aj},\qquad
\Omega_c=\omega_{i,R}+\omega_{j,R}-\omega_{k,R},\qquad
m_c=m_i+m_j-m_k.
\]

With consistently normalized occupations and event rates, assemble

\[
\dot N_a=2\omega_{a,I}N_a+\sum_c\nu_{ac}\mathcal R_c(\boldsymbol N).
\]

This automatically handles the factor of two when `i=j`. Do not separately
count `(i,j)` and `(j,i)`; put source multiplicities and normalization factors
in the rate calculation. At leading classical quartic order, isolated
stimulated channels have cubic occupation dependence. This is not an exact
closure for the full cosine potential.

At fixed linear frequencies, each column obeys
\(\sum_a\omega_{a,R}\nu_{ac}=-\Omega_c\) and
\(\sum_a m_a\nu_{ac}=-m_c\). These are cheap algebraic checks before any radial
solve. For a positive-frequency escaping channel,
\(F_{E,c}=\Omega_c\mathcal R_c\) and
\(F_{J,c}=m_c\mathcal R_c\); horizon fluxes require their signed convention.
Combine coherent sources before squaring whenever they drive the same output
field and their interference survives averaging.

[Witte and Mummery, Sections II.3–II.4](https://arxiv.org/html/2412.03655v3#S2.SS4)
provides a comparison using quartic scattering rates, including relativistic
rate calculations, and evolution with higher bound levels. Its rate-network
organization is useful; it does not replace the nonlinear background problem.

## Extension that still needs a derivation

For finite amplitudes, replacing constant rates by amplitude-dependent rates
is a proposal to test, not an established closure. Interaction energy and
profile deformation matter; in general one cannot set
\(E=\sum_a\omega_a(A)N_a\) and reuse the linear equations unchanged.

Derive mode-resolved slow equations by projection or an averaged-action
calculation, with normalization appropriate to the open Kerr boundary problem.
Require their predicted total energy and angular momentum changes to match
the boundary fluxes. A real scalar has no exact particle-number conservation
law that would automatically supply extra balances.

An amplitude-only model also requires averaging over relative phases. Test
this against resonant combinations: if a detuning is comparable to a slow
transfer rate, retain the corresponding phase and coherent coupling. Additional
Fourier components may be slaved driven fields or independently evolving
clouds; a near-resonant bound response needs special treatment.

## Incremental work

1. **Secondary eigenproblem at a frozen primary amplitude.** Construct the
   projected background correction from
   \(\mu^2(r^2+a^2\cos^2\theta)(1-\cos\phi_1)\psi_2\), with the retained
   frequency components specified explicitly. First reproduce the carrier
   truncation of Section III. As \(A_1\to0\), the correction must vanish and
   the frequency must recover the independent linear `322` benchmark. Keep
   the secondary normalization separate from the primary outer amplitude.
2. **Driven responses.** Solve for the `m=0` and `m=3` channels with ingoing
   horizon and decaying/outgoing outer conditions as appropriate. Check the
   weak-amplitude scaling, flux signs, and angular-momentum flux relation.
   Account for each channel once when combining eigenvalue growth and fluxes.
3. **Two-cloud evolution.** Add angular momentum and its derivatives, construct
   the balance matrix, and reproduce Figure 6 at `alpha=0.42`, `spin=0.99`.
   Check the single-cloud limit, time-step convergence, and spatial truncation.
4. **General mode transfers.** Validate an automatically assembled weak-field
   channel network, then derive and test its finite-amplitude extension against
   the two-cloud calculation before evolving four or five clouds.

Several weak secondary clouds could initially share a one-dimensional table
in primary amplitude. If their backreaction becomes appreciable, that
hierarchy fails and the profiles must depend on additional amplitudes.
Avoid committing to a dense five-dimensional table: local continuation along
trajectories is a candidate once the closure is established.

Five slow amplitude equations are manageable. The costly work is the radial
response calculation, channel enumeration, and a controlled choice of
closure and truncation. The existing single-cloud interpolation checks do not
establish spatial convergence for these new calculations.
