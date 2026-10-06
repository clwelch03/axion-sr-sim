# Further adiabatic work — October 6, 2026

This note evaluates six proposed extensions: additional bound modes, coupled
black-hole evolution, dependence on the axion decay constant, a handoff to
direct field evolution near a possible bosenova, and an exploration of
self-gravity, and gravitational-wave (GW) losses. It records a research proposal, not an implementation plan
already under way or a claim that new effects have been demonstrated. The
literature search was targeted rather than exhaustive.

**Recommendation:** make black-hole evolution and decay-constant dependence the
next main project. Begin with a cheap physical normalization of the existing
solutions, including diagnostics of cloud self-gravity and GW emission. Include
GW energy and angular-momentum losses in the coupled evolution; quantify any
channels omitted as negligible. Develop conservative gravitational feedback if
the self-gravity diagnostic indicates a consequential correction.
Then investigate the omitted nonlinear response and selected higher
modes where they could change the result. Use short direct PDE evolutions as
early validation; develop the full collapse-and-restart scheme if physically
accessible trajectories approach an instability.

| Stage | Work | Scientific question |
|---|---|---|
| First | Normalize existing results across decay constants | What physical cloud masses and torques do our solutions represent? |
| Alongside normalization | Estimate self-gravity on saved profiles | Could gravitational frequency shifts change a threshold or transfer channel? |
| Alongside normalization | Estimate GW power and torque on saved profiles | How much energy and angular momentum escape, and when does depletion matter? |
| Next | Evolve BH mass and spin together with the clouds, including GW losses | Does nonlinear evolution change the amount or timing of spin extraction? |
| Conditional, before precision conclusions | Add gravitational response to the nonlinear adiabatic family | Does cloud gravity change fluxes, saturation, or the spin history? |
| Then | Mutual cloud deformation and selected additional modes | Where does the two-cloud prediction become unreliable? |
| Early validation, later extension | Short PDE checks, followed conditionally by collapse and restart | Does the adiabatic trajectory survive direct evolution, and what happens when it fails? |

The numbering below preserves the original four ideas and adds self-gravity
as item 5 and GW losses as item 6; it does not specify execution order.

## Starting point and comparison with prior work

The [October 5 handoff](../HANDOFF.md) records completion of the two-cloud
reproduction at $M\mu=0.3$ and $0.42$, with dimensionless BH spin $0.99$.
The primary cloud is nonlinear; the secondary responds to that primary, while
its feedback on the primary is neglected. The evolution follows energy and
angular-momentum balance, including the driven scalar horizon and infinity
channels. Neither the existing outgoing $m=3$ scalar field nor the saved
energy trajectories include GW radiation reaction.
See the [multimode design](multimode-adiabatic-design.md) for the approximation
and the [campaign README](../runs/adiabatic/multimode_secondary/README.md) for
the numerical evidence.

Two especially close comparisons constrain what would constitute a new result:

- [Omiya et al., *Deci-Hz gravitational waves from the self-interacting axion
  cloud around the rotating stellar mass black hole* (2024)](https://arxiv.org/html/2404.16265v1)
  already combine four bound modes, BH mass and spin evolution, and varying
  decay constant. Their evolution omits nonlinear profile deformation and
  associated frequency shifts. They also neglect GW losses in the evolution
  (Sec. III.1), while computing GW signals afterward (Sec. V.1). Sections
  III–IV benchmark the scalar/BH sector, not a GW-inclusive evolution.
- [Witte & Mummery, *Stepping Up Superradiance Constraints on Axions*
  (2025)](https://arxiv.org/html/2412.03655v3) treat relativistic quartic
  scattering and relevant bound levels through principal number $n=5$, coupled
  to BH evolution. Their higher-level comparisons show that an added state can
  suppress previously growing clouds, change spin-down, and remove an apparent
  bosenova. They also identify limitations of the level truncation for higher
  azimuthal numbers.

Thus, neither adding modes nor evolving the BH is by itself a novelty claim.
The candidate contribution is a controlled comparison of nonlinear adiabatic
evolution with existing perturbative transfer models, establishing where
profile response changes an observable and where simpler models remain
accurate. The [September 28 review](multimode-literature-review-2026-09-28.md)
discusses the broader literature and related candidate gaps.

## 3. Dependence on the decay constant

Let $\phi$ be the canonically normalized real axion field, $f_a$ its decay
constant, and $\mu$ its mass. For the cosine potential,

$$
V(\phi)=\mu^2 f_a^2\left[1-\cos(\phi/f_a)\right].
$$

Define $\Phi=\phi/f_a$, the dimensionless field used by the scalar solver.
On a prescribed Kerr background its equation is

$$
\Box_g\Phi-\mu^2\sin\Phi=0,
$$

where $\Box_g$ is the wave operator for that metric. There is no explicit
$f_a$ in this equation. The stress tensor, however, has an overall $f_a^2$:

$$
T_{ab}=f_a^2\left[
\partial_a\Phi\partial_b\Phi
-g_{ab}\left(\frac12\nabla_c\Phi\nabla^c\Phi
+\mu^2(1-\cos\Phi)\right)\right].
$$

Consequently, at fixed $\mu$, geometry, and initial dimensionless field,
changing $f_a$ leaves the normalized scalar evolution unchanged. Physical
cloud energies, angular momenta, and scalar fluxes scale as $f_a^2$. This
statement applies to the fixed-background scalar problem; it does not assume
that gravitational radiation or metric response shares the same scaling.

This explains why a fixed-background scan alone would mostly re-express the
existing solutions. The dependence becomes dynamically interesting when:

- The BH responds to the physical horizon fluxes.
- Initial conditions are specified at fixed physical energy or occupation,
  rather than fixed $\Phi$. Even the time to grow from a specified seed then
  needs consistent normalization.
- Cloud self-gravity or gravitational radiation becomes relevant.
- Changing $f_a$ also changes $\mu$, as along the QCD-axion relation.

For an ALP study, first vary $f_a$ at fixed $\mu$ to isolate the interaction
dependence. Treat the QCD relation as a separate path through parameter space.

**First deliverable:** convert the existing two-cloud trajectories to physical
cloud mass fractions, horizon energy fluxes, and horizon torques for several
$f_a$. Compare the accumulated extraction with the BH's finite mass and spin
reservoir. This identifies when a frozen-BH trajectory becomes inconsistent;
it is not itself a prediction of the subsequent spin history. Record the
Planck-mass convention explicitly when converting code units.
Use the same normalized profiles for the self-gravity diagnostic in item 5.
Estimate GW power, torque, and accumulated losses as described in item 6;
post-processing diagnoses their importance but does not produce a trajectory
with radiation reaction.

## 2. Coupled black-hole and cloud evolution

This is the most direct connection to
[Ning, Safdi & Welch, *No Evidence for Superradiant Axions in LIGO-Virgo-KAGRA
GWTC-5 Binary Black Hole Spins*](https://arxiv.org/html/2607.01317v1).
Appendix F reports that its self-interaction checks at the QCD decay constant
generally leave the conclusions intact, while explicitly identifying other
$f_a$ values as future work. The natural extension is a broader ALP study with
the QCD relation as a benchmark. There is no demonstrated reason yet to expect
a substantial revision of the QCD conclusion.

In geometrized units $G=c=1$, let $M_{\rm BH}$ and $J_{\rm BH}$ be the BH mass
and angular momentum, and define its dimensionless spin by
$\chi=J_{\rm BH}/M_{\rm BH}^2$. Let $F_E^H$ and $F_J^H$ denote total energy
and angular-momentum fluxes through the horizon, positive into the BH.
At leading adiabatic order, with no external accretion,

$$
\dot M_{\rm BH}=F_E^H,\qquad
\dot J_{\rm BH}=F_J^H,
$$

$$
\dot\chi=\frac{F_J^H}{M_{\rm BH}^2}
-\frac{2\chi F_E^H}{M_{\rm BH}}.
$$

The horizon fluxes must include every retained channel, including the driven
scalar fields and any retained GW horizon flux. Infinity emission enters the cloud's balance; it is not directly a
flux into the BH. A cloud with approximately constant energy can still
continuously extract BH spin if superradiant input balances emission. Maximum
cloud energy therefore does not determine cumulative spin extraction.

The new BH equations are simple. The substantial work is making the cloud
family follow the changing background. With fixed particle mass $\mu$,

$$
\alpha(t)=M_{\rm BH}(t)\mu
$$

changes along with $\chi$. Profiles, frequencies, and channel fluxes must be
updated accordingly. Our two existing values of $\alpha$ at one spin are
insufficient for that evolution. Continuation along selected trajectories
could establish feasibility before constructing a broad interpolation table.

Energy accounting also needs a derivation at the chosen approximation order.
For example, if $E_{\rm cl}$ denotes cloud energy and $A_i$ its amplitude
coordinates, differentiating a family
$E_{\rm cl}(\{A_i\},M_{\rm BH},J_{\rm BH})$ introduces derivatives with
respect to the BH parameters as well as the amplitudes. A changing metric can
also exchange energy with the field. These effects must be treated consistently;
simply inserting a new mass into the fixed-background balance law is not a
complete derivation. Use fixed reference units so changes in the time or
energy normalization do not masquerade as evolution.

This proposal initially uses an instantaneous Kerr background with an adiabatic
GW loss prescription. Conservative cloud gravity is a separate correction. Carry the
self-gravity diagnostic in item 5 alongside this calculation, and use it to
decide whether gravitational feedback is needed before interpreting the spin
history quantitatively.

**First scientific comparison:** at matched physical initial conditions,
compare perturbative transfer equations with the present nonlinear adiabatic
approximation, using the same retained modes and loss processes. Measure
$\chi(t)$, cloud energies, cumulative horizon angular-momentum flux, and the
spin reached at specified ages. Recover the dilute limit first. Independently
checking the perturbative flux prefactors, still listed as open in the handoff,
would make that comparison interpretable.

A useful result could be a significant correction to the spin history, a
shifted boundary for effective spin-down, or an established domain in which
the simpler model is accurate. A population-level inference update would
follow a validated single-BH forward model.

## 6. Gravitational-wave losses and the combined budget

**Account for GW energy and angular momentum in the evolution.** A Poisson
self-gravity solve does not supply these radiative losses. The following is a
proposed leading adiabatic balance prescription; no GW flux calculator or
radiation reaction has been implemented in the current two-cloud solver.

**Physical balance and signs.** In geometrized units, define
$P_{\rm GW}^{\infty}$ and $Q_{\rm GW}^{\infty}$ as outgoing GW energy and
angular-momentum fluxes. The power is nonnegative; the torque is signed.
Define $P_{\rm GW}^{H},Q_{\rm GW}^{H}$ as the corresponding signed fluxes
into the BH. The horizon energy flux can be negative for superradiant waves.
Write $F_{E,\phi}^{H,\infty}$ and $F_{J,\phi}^{H,\infty}$ for the scalar
fluxes, with the same boundary sign conventions. On a fixed background,
the proposed secular cloud balances are

$$
\dot E_{\rm cl}=-F_{E,\phi}^{H}-F_{E,\phi}^{\infty}
                 -P_{\rm GW}^{H}-P_{\rm GW}^{\infty},
\qquad
\dot J_{\rm cl}=-F_{J,\phi}^{H}-F_{J,\phi}^{\infty}
                 -Q_{\rm GW}^{H}-Q_{\rm GW}^{\infty}.
$$

For an evolving BH, use $\dot M_{\rm BH}=F_{E,\phi}^{H}+P_{\rm GW}^{H}$
and $\dot J_{\rm BH}=F_{J,\phi}^{H}+Q_{\rm GW}^{H}$ in item 2. GW power
escaping to infinity is removed from the cloud, not subtracted again from
the BH. At the retained adiabatic order the combined accounting must give

$$
\frac{d}{dt}(M_{\rm BH}+E_{\rm cl})
 =-F_{E,\phi}^{\infty}-P_{\rm GW}^{\infty},\qquad
\frac{d}{dt}(J_{\rm BH}+J_{\rm cl})
 =-F_{J,\phi}^{\infty}-Q_{\rm GW}^{\infty}.
$$

These are target budget identities, not an exact decomposition of the energy
of a dynamical spacetime. Changing-background work and conservative binding
energy must be included consistently at the chosen order, as discussed in
item 2. Neglecting the GW horizon flux is a separate approximation to check.

**Which channels.** For a radiative component proportional to
$e^{-i\Omega t+im_{\rm GW}\varphi}$, define $\Omega>0$ as its angular
frequency and $m_{\rm GW}$ as its signed azimuthal number. Its fluxes obey
$Q_{\rm GW}^{b}=(m_{\rm GW}/\Omega)P_{\rm GW}^{b}$ at each boundary
$b=H,\infty$. Both members of a real-field conjugate pair represent one
physical channel and must not be counted twice. In the dilute two-cloud
limit, with real bound frequencies $\omega_1,\omega_2$ and $m_1=1,m_2=2$,
the candidate channels are

| Process | $\Omega$ | $m_{\rm GW}$ | Leading dependence on code cloud energies |
|---|---|---|---|
| Primary pair annihilation | $2\omega_1$ | $2$ | $E_1^2$ |
| Secondary pair annihilation | $2\omega_2$ | $4$ | $E_2^2$ |
| Mixed pair annihilation | $\omega_1+\omega_2$ | $3$ | $E_1E_2$ |
| Transition $2\to1$, when $\omega_2>\omega_1$ | $\omega_2-\omega_1$ | $1$ | $E_1E_2$ |

These label source frequencies, not assumed equal or nonzero coefficients.
Angular selection rules and radial overlaps determine the powers. Transitions
and annihilations require distinct coefficients; nonlinear harmonics can add
channels. The sourced spin-2 Teukolsky calculation is a route to relativistic
fluxes. [Omiya et al., Sec. V.1](https://arxiv.org/html/2404.16265v1#S5.SS1)
describe this construction. For these $211$/$322$ clouds, annihilation and
transition channels are also discussed by
[Baryakhtar et al., Sec. VII](https://arxiv.org/html/2011.11646#S7).
Do not assume a long-wavelength quadrupole estimate supplies a controlled
annihilation coefficient at $\alpha=0.3$ or $0.42$.

**Conversion to the existing balance equations.** Define our normalization
$\epsilon_f=(f_a/M_{\rm Pl})^2$ with unreduced
$M_{\rm Pl}=\sqrt{\hbar c/G}$, and fix reference mass $M_0$. In geometrized
units let $\hat t=t/M_0$, $E_{\rm cl}=\epsilon_f M_0(E_1+E_2)$, and
$J_{\rm cl}=\epsilon_f M_0^2 J_{\rm code}$. Thus the extra code fluxes are

$$
g_E=\frac{P_{\rm GW}^{H}+P_{\rm GW}^{\infty}}{\epsilon_f},\qquad
g_J=\frac{Q_{\rm GW}^{H}+Q_{\rm GW}^{\infty}}{\epsilon_f M_0}.
$$

At fixed normalized scalar profile and background, leading GW power scales
as $\epsilon_f^2$, whereas physical scalar energy and scalar power scale as
$\epsilon_f$. Therefore $g_E$ and $g_J$ scale as $\epsilon_f$. For example,
define a dimensionless dilute-channel coefficient $C_{ij}^{b}$ by
$P_{ij}^{b}=\epsilon_f^2 C_{ij}^{b}E_iE_j$; then
$g_{E,ij}^{b}=\epsilon_f C_{ij}^{b}E_iE_j$ and
$g_{J,ij}^{b}=(m_{\rm GW}/\hat\Omega)g_{E,ij}^{b}$, where
$\hat\Omega=M_0\Omega$. The $f_a^4$ physical-power scaling is consistent
with [Omiya et al., Sec. V.1](https://arxiv.org/html/2404.16265v1#S5.SS1).
The coefficients must specify their energy normalization and symmetry factors;
they are not present in our scalar flux tables.

Using the existing two-cloud approximation for $J_{\rm code}$, let primes
denote $d/d\hat t$ and take all frequencies below in code units. The extension
of the current linear balance solve is

$$
\begin{pmatrix}
1&1\\
\displaystyle\frac{1}{\omega_1}
 -\frac{2E_2}{\omega_2^2}\frac{d\omega_2}{dE_1}
 &\displaystyle\frac{2}{\omega_2}
\end{pmatrix}
\begin{pmatrix}E_1'\\E_2'\end{pmatrix}
=-\begin{pmatrix}F_E+g_E\\F_J+g_J\end{pmatrix},
$$

where $F_E,F_J$ retain their current scalar-only meanings. This preserves the
secondary-frequency derivative already in `two_cloud_rates`. For saved time
$\tau=\gamma_1\hat t$, with $\gamma_1$ the dimensionless linear primary
growth rate, divide both rates by $\gamma_1$. A GW power must not be subtracted
independently from both cloud energies. In the dilute constant-frequency
limit a pure transition of outgoing power $P$ instead gives
$\dot E_1=+\omega_1P/(\omega_2-\omega_1)$ and
$\dot E_2=-\omega_2P/(\omega_2-\omega_1)$ in consistent physical units:
the lower cloud gains energy while the total loses $P$.

**First calculation and decisive checks.** Estimate channel powers and torques
on the saved trajectories, then integrate
$E_{\rm GW}^{\infty}(t)=\int_0^t P_{\rm GW}^{\infty}dt'$ and
$J_{\rm GW}^{\infty}(t)=\int_0^t Q_{\rm GW}^{\infty}dt'$ with the correct
growth-time conversion. Compare these with cloud budgets and accumulated
horizon extraction, and compare instantaneous losses with individual scalar
inputs and outputs. Near scalar balance, the net scalar flux can vanish;
a ratio to that net flux is not a useful smallness test. Check late depletion
after superradiant input shuts off as well as the growth phase.

Validate coefficients in a matched dilute benchmark before using them for
nonlinear profiles. The first integration checks are recovery of the saved
scalar-only trajectory when GW terms vanish; pure single-cloud annihilation
with constant coefficient, $E(t)=E(0)/[1+kE(0)t]$ for $\dot E=-kE^2$;
the transition balance above; and conservation of the combined BH, cloud,
and accumulated outgoing energy and angular momentum at the retained order.
Report GW-inclusive and GW-free results at matched initial conditions.
Post-processing alone cannot justify simply subtracting accumulated GW energy
from a saved trajectory once the loss changes its amplitudes or fluxes.

## 5. An exploration of self-gravity

**Include the exploration now; make dynamical feedback conditional.** Start
with saved profiles during physical normalization, alongside the first
coupled-BH calculation. If gravity is consequential, incorporate its response
before a broad parameter scan or precision spin-down conclusions. It need not
wait for extra modes or a bosenova. The
[gravitational backreaction review](gravitational-backreaction-review.md)
describes the gravity approximations and time-domain implementation issues;
the [adiabatic self-gravity review](adiabatic-self-gravity-review.md) gives the
more direct route for the present branch solver.

BH evolution and cloud self-gravity are separate effects. Updating Kerr mass
and spin accounts for horizon exchange, while the extended cloud also sources
a gravitational field. Adding the cloud mass to the central Kerr mass does
not reproduce that field. Likewise, strong cosine self-interaction does not
by itself imply strong gravitational deformation. At fixed dimensionless
profile the additional coupling is

$$
\kappa_g=8\pi Gf_a^2=(f_a/\bar M_{\rm Pl})^2,
\qquad \bar M_{\rm Pl}=(8\pi G)^{-1/2}.
$$

This makes the decay-constant scan a natural place to assess gravity, although
the physically accessible profiles must ultimately come from the coupled
evolution.

1. **Diagnose gravity on representative saved profiles.** Record
   $M_c/M_{\rm BH}$, the estimated potential magnitude, and frequency shifts
   across the proposed $f_a$ values. Begin with cloud-region Poisson multipoles:
   monopole, then quadrupole, with independent radial and angular convergence.
   Recover the dilute hydrogenic $211$ benchmark from
   [Baryakhtar et al., Appendix G](https://arxiv.org/html/2011.11646#A7)
   before applying the diagnostic to distorted profiles. Define the averaging
   of the real-field source explicitly; averaging fast axion oscillations
   need not remove slow, nonaxisymmetric beats between clouds. Poisson gravity
   is a diagnostic whose validity requires a slowly varying, nonrelativistic
   source, not merely a small cloud mass. At $\alpha=0.3$ and $0.42$, assess
   relativistic uncertainty using matched dilute benchmarks, including
   [May, East & Siemonsen](https://arxiv.org/abs/2410.21442), before treating
   the estimated shifts as quantitative predictions for our nonlinear clouds.
2. **Decide by observable sensitivity.** Compare the gravitational shifts
   with the superradiant detuning $m\Omega_H-\omega$, relevant bound-resonance
   detunings, and the distance of each combination frequency from its emission
   threshold. Use shifts of the actual frequency combinations, since common
   shifts can cancel. Near a sensitive channel, test the resulting change in
   transfer rates and horizon torque; a small $\delta\omega/\omega$ alone
   does not establish a small effect on $\chi(t)$. State a target accuracy
   and include numerical and gravity-model uncertainty. If the estimated
   corrections could exceed it, promote gravitational feedback to the next
   calculation. If they remain below it, document the tested domain and defer
   feedback. A future phase-sensitive GW study needs a separate accumulated
   phase criterion, even when spin-down is insensitive.
3. **Perturb the nonlinear adiabatic family first.** Solve for the weak
   cloud-sourced metric, then differentiate the scalar matching problem to
   obtain frequency, profile, and flux corrections. This retains the full
   cosine response and can reuse the branch Jacobian; the new metric source,
   operator, and boundary corrections still require derivation. A frequency
   shift inserted into unchanged flux tables is insufficient. Validate the
   dilute limit and small-coupling response, then compare corrected and
   uncorrected spin histories at matched initial conditions. Include binding
   energy consistently and separate cloud gravity from the evolving Kerr
   mass/spin pieces to avoid double counting. Check the scalar response itself:
   near a poorly conditioned branch, a small metric perturbation need not
   produce a small profile change.

The first deliverable is a self-gravity sensitivity assessment accompanying
the normalized trajectories. The conditional scientific result is the change
in torque, saturation, and $\chi(t)$ caused by gravitational response of the
nonlinear clouds. If Poisson accuracy is inadequate, advance to sourced
linearized gravity on Kerr, assessing omitted stresses, currents, and
oscillating sectors. Full dynamical Einstein-scalar evolution is a separate,
later project if gravitational deformation requires it. A self-interaction
bosenova alone does not establish that need; short fixed-background PDE
validation can proceed where the gravity diagnostic supports that approximation.

## 1. Additional modes and the missing evolution equations

Here, an additional mode means an independent bound cloud labeled by its
linear-limit quantum numbers. It is distinct from adding angular harmonics or
combination-frequency radiation to represent an existing cloud more accurately.

For $K$ independent amplitudes $A_i$, energy and angular-momentum balances
provide only two equations. At fixed BH parameters they have the form

$$
\sum_{i=1}^{K}\frac{\partial E}{\partial A_i}\dot A_i=-F_E,
\qquad
\sum_{i=1}^{K}\frac{\partial J}{\partial A_i}\dot A_i=-F_J,
$$

where $E,J$ describe the retained cloud configuration and $F_E,F_J$ are its
net boundary losses. For $K>2$, these equations do not determine all amplitude
rates. The missing information is mode-resolved dynamics. Another exact global
conservation law need not exist; in particular, particle number is not an exact
conserved charge of the full real cosine field.

Whitham's averaged-Lagrangian method is a plausible way to derive those
dynamics. The idea is to average over rapid oscillations while retaining
slowly varying phases and amplitudes. If $\overline L$ denotes the spatially
integrated, phase-averaged Lagrangian and $\omega_i$ the real frequency of the
$i$th independent phase, its conjugate wave action is schematically

$$
I_i=\frac{\partial\overline L}{\partial\omega_i}.
$$

Variation of the phases supplies action evolution equations. This is the
structural reason the method can supply more information than total energy
and angular momentum. General multiphase formulations exist; see
[*Nonlinear Theory for Coalescing Characteristics in Multiphase Whitham
Modulation Theory*](https://d-nb.info/1226144489/34).

For our problem, the horizon and outgoing radiation require action source
terms. These have to be derived consistently with the slowly growing or
decaying modes. Near-resonant combinations of phases may evolve slowly and
must then remain dynamical variables rather than being averaged away.
Independent conserved actions cannot simply be imposed on the dissipative
system. This is an unresolved derivation, not a ready-to-use closure.

Before undertaking a general construction, use three milestones:

1. **Measure feedback from the secondary onto the primary.** The present
   approximation omits it. A leading response calculation tests whether
   mutual deformation changes frequencies, fluxes, or evolution along an
   accessible trajectory. Two mutually deforming clouds are already a useful
   extension before introducing the extra closure problem.
2. **Identify consequential additional modes.** Test candidate perturbations
   on the occupied background, including required frequency sidebands.
   Determine whether growth or forced population can become important within
   the available evolution time. A positive growth rate alone does not show
   that a mode changes the observable.
3. **Derive and validate a three-cloud evolution.** Require recovery of the
   known perturbative transfer equations in the dilute limit and the existing
   two-cloud balance equations in the corresponding truncation. Check energy
   and angular-momentum accounting for each process.

Choose the retained set by its effect on the desired observable, rather than
by the next angular number alone. A small population can still mediate an
important depletion channel. Conversely, a larger mode set is not automatically
converged. Near an emission threshold or a bound resonance, even a modest
frequency shift might affect a transfer rate; that is a mechanism to test,
not an effect already established for our trajectories.

## 4. Handoff to direct PDE evolution

Direct nonlinear bosenova and multimode simulations already exist, including
[Yoshino & Kodama, *Bosenova and Axiverse*](https://arxiv.org/html/1505.00714).
The candidate opportunity is a controlled connection between long adiabatic
growth and a short episode of direct field evolution, eventually followed by
return to an adiabatic description. Establishing priority for that complete
workflow would require a more focused literature review.

The first deliverable should be a validation handoff well before a suspected
collapse. Reconstruct a saved adiabatic field and its time derivative, convert
them to the PDE solver's field and momentum variables, and evolve for several
internal dynamical times. Preserve the phase choices, retained harmonics,
coordinate conventions, and appropriate radiation content. Compare the field,
boundary fluxes, and response to a small perturbation with the adiabatic
prediction. Repeat with numerical refinement to distinguish physical departure
from initialization or discretization error.

There is a concrete implementation gap: the current
[`nonlinear_tower.py`](../src/dolan/nonlinear_tower.py) wrapper only supports
positive odd-$m$ towers, with negative partners reconstructed by reality.
A two-cloud $m=1,2$ handoff requires even-$m$ and generated $m=0$ support,
together with the appropriate angular sectors. The
[nonlinear progress record](nonlinear-progress.md) documents existing source
and short-integration checks; it does not establish accurate collapse
evolution.

For an eventual automatic switch, seek evidence of a physical instability or
loss of timescale separation. Useful diagnostics include resolved growth of
perturbations, a soft internal mode, and rapid change of the configuration
relative to its internal response time. A continuation failure, branch turning
point, or large field amplitude alone is not a complete stability diagnosis.
The axion oscillation period alone is also an insufficient clock: an internal
response time can become long near an instability.

The switch should occur while the adiabatic configuration still provides
credible initial data. Compare results from several earlier switch points;
agreement after initial transients would test sensitivity to the handoff.
Freezing the BH during the short PDE episode is an approximation to assess
from its duration and accumulated horizon fluxes.

Returning to adiabatic evolution is a separate milestone. After radiation has
departed and the remnant has relaxed, determine its bound-mode content, energy,
and angular momentum, and establish whether the retained adiabatic family can
represent it. A disturbed remnant need not fit the original two-cloud family.

The decision to build the full bosenova cycle should follow evidence that
physical trajectories reach an instability after relevant depletion channels
are included. An artificially overfilled cloud can be a numerical benchmark,
but does not establish that superradiant growth reaches that state.

## Immediate next result

Start with physical normalization, GW loss estimates, and a frozen-profile
self-gravity diagnostic for the two completed trajectories, then construct a
limited coupled-BH calculation including GW losses and its matched perturbative
comparison. The first research figure
should show **$\chi(t)$ for several
$f_a$**, accompanied by cloud energies and accumulated horizon
angular-momentum flux, plus cumulative GW energy and angular momentum emitted
to infinity. Report the self-gravity sensitivity alongside it;
include gravitational feedback before precision interpretation wherever the
diagnostic indicates it could change the conclusion.

Use that comparison to choose the next calculation: gravitational response if
self-gravity shifts a consequential threshold or torque, mutual profile
response if scalar deformation changes the torque, additional modes if the
truncation controls the outcome, or direct PDE evolution if the adiabatic
solution shows evidence of instability. Short PDE validation remains useful
even if all physically accessible trajectories avoid a bosenova.
