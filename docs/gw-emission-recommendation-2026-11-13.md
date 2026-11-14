# Adding gravitational-wave losses to the cloud evolution: recommendation

Date: 2026-11-13. Status: recommendation only; nothing has been implemented.

## Recommendation

Add gravitational-wave (GW) emission as a **tabulated loss term**: the power is
a known coefficient times the product of two cloud masses, with the
$|211\rangle$ and $|322\rangle$ annihilation coefficients taken from the
relativistic (Teukolsky) fits published with SuperRad (Siemonsen, May & East,
arXiv 2211.03845), and the $|322\rangle \to |211\rangle$ transition coefficient
taken from the non-relativistic multipole formula (Baryakhtar et al.
2011.11646; Cannizzaro, Spieksma & Witte 2611.15865 Appendix C.2). This is a
few days of work and it captures the part of GW emission that matters for our
runs.

Do **not** build our own sourced spin-2 Teukolsky solver now. It is a project
of several weeks, and the cases where GW emission matters in our runs are the
cases where the cloud profile is closest to the linear one, which is exactly
what the SuperRad coefficient already assumes.

Do **not** use the flat-space annihilation formula of 2611.15865 Appendix C.1,
or the commonly quoted Schwarzschild coefficient, for $|211\rangle$
annihilation. They differ from the relativistic result by an order of magnitude
or more at $\alpha = 0.3$, and that factor decides whether GW emission is a
20% or a 70% effect in our long runs.

The reasoning and the numbers follow.

## Two corrections to how the question was posed

**arXiv 2404.16265 is not a recipe for the relativistic calculation.** That
paper (Omiya, Takahashi, Tanaka & Yoshino) explicitly neglects GW emission in
its evolution ("we neglect the effect of gravitational waves on the
evolution", Sec. III.A). It computes the GW strain afterwards, as an
observable, from the *linear* Kerr bound states. Its description of the method
is one paragraph (its Eqs. 29–34): the amplitude at infinity is a radial
integral of the ingoing spin-2 Teukolsky solution against a source built from
the scalar stress-energy, but the source is not written out ("can be found in"
Teukolsky 1973) and the solver is not described. No table or fit of the power
is given, only plots (its Fig. 10). The actual method is in Yoshino & Kodama
(1312.2326), Siemonsen, May & East (2211.03845), and Guo, Bao & Zhang
(2212.07186). So "the full relativistic road" means implementing those papers.

**Appendix C of 2611.15865 contains two different approximations of very
different quality.**

- *Transitions (C.2):* a multipole expansion (mass and current moments of the
  overlap of two hydrogenic levels). Their Fig. 17 compares it with the
  relativistic transition fluxes of 2404.16265 and reports a "moderate"
  discrepancy even at large $\alpha$. This one is usable.
- *Annihilation (C.1):* the flat-space formula for a harmonic source applied to
  hydrogenic wavefunctions. No relativistic check is shown. Arvanitaki,
  Baryakhtar & Huang (1411.2263) state that this same flat-space hydrogenic
  calculation has a cancellation of the leading term for the $2p$ level.
  Evaluating their closed form gives a power roughly 250 times below the
  relativistic value at $\alpha = 0.3$. **Caveat:** we have not checked the
  numbers actually used in 2611.15865 against this; it is an inference from the
  method they describe, and should be checked before quoting it to anyone.

Equation numbers for 2611.15865 in this document come from the arXiv HTML,
because the PDF did not parse. They need rechecking against the PDF before
being cited elsewhere. Neither paper is in `docs/papers/` yet.

## How large is the effect in our runs?

Write the annihilation power of a single level as
$P = \tilde P(\alpha)\,(M_c/M)^2$ in units $G = c = 1$, where $M_c$ is the cloud
mass, $M$ the black-hole mass, and $\tilde P$ (notation borrowed from the
SuperRad paper's normalization, not used elsewhere in this repo) is a
dimensionless coefficient. With no resupply, the cloud mass then decays as
$M_c(t) = M_c(0)/(1 + t/\tau)$ with $\tau = M / (\tilde P \, M_c(0)/M)$.

Coefficients for $|211\rangle$ at $\alpha = 0.3$:

| Source | $\tilde P_{211}$ | Ratio to SuperRad |
|---|---|---|
| SuperRad (interpolated Teukolsky data, evaluated at the saturation spin) | $8.4\times10^{-11}$ | 1 |
| Yoshino & Kodama flat-space, $\alpha^{14}/640$ | $7.5\times10^{-11}$ | 0.9 (partly coincidence) |
| Brito, Cardoso & Pani Schwarzschild, $0.0249\,\alpha^{14}$ | $1.2\times10^{-9}$ | 14 |
| Baryakhtar et al. 2021, $\gamma \simeq 10^{-2}\alpha^{14}\mu$ (converted to power) | $\approx 1\times10^{-9}$ | 11 |
| Flat-space hydrogenic with the $2p$ cancellation (Arvanitaki et al.) | $\approx 3\times10^{-13}$ | 1/250 |

Consequences for the existing spindown runs ($\alpha = 0.3$, initial spin 0.9):

| Run | Peak $M_c/M_0$ | $\tau$ (SuperRad) | Run length | Fraction of cloud radiated |
|---|---|---|---|---|
| $f_a = 3.7\times10^{16}$ GeV | 0.30% | $4\times10^{12}\,M_0$ | $\le 10^{10}\,M_0$ | below 0.3%: negligible |
| $f_a = 5.7\times10^{17}$ GeV and larger, and the free case | 1.85% | $6.5\times10^{11}\,M_0$ | $1.4\times10^{11}$ to $3\times10^{11}\,M_0$ | about 18% to 32% |

For the large-$f_a$ runs, the same estimate with the Schwarzschild or
Baryakhtar coefficient gives about 70% radiated by $1.4\times10^{11}\,M_0$.
In the present model the $|211\rangle$ cloud in the $5.7\times10^{17}$ GeV run
is instead reabsorbed by the hole by that time. So GW emission changes the
answer to "how much of the extracted energy returns to the hole", and the size
of the change depends on a coefficient that the simple formulas get wrong by a
factor of ten. The note in `docs/session-2026-11-10.md` that estimated
$\tau \approx 10^{11}\,M$ for this cloud is about six times too short by the
SuperRad coefficient.

The other channels, at $\alpha = 0.3$:

- $|322\rangle$ annihilation: $\tilde P_{322} \approx 1.2\times10^{-16}$,
  timescale above $10^{17}\,M$. Negligible, but free to include.
- $|322\rangle \to |211\rangle$ transition: moves quanta from $|322\rangle$ to
  $|211\rangle$ on a timescale of about $2\times10^{12}\,M$ for a 1.85%
  $|211\rangle$ cloud. Slower than annihilation and beyond current run
  lengths, but it carries little energy and one unit of angular momentum per
  quantum, so it changes the composition, not the energy budget. Cheap to
  include.
- Cross-annihilation (one quantum from each level, frequency
  $\omega_1 + \omega_2$): starts at octupole order; 2404.16265 calls it "much
  smaller" without numbers. Leave out and say so.
- GW flux into the horizon: not quantified by any of the sources read here.
  Leave out and say so.

These estimates use published coefficients and peak cloud masses read from our
diagnostics files. They are not a substitute for running the evolution with
the term included, because the cloud is also exchanging quanta with the hole
during this time.

## Why the tabulated coefficient is adequate here

The SuperRad coefficient assumes the linear Kerr bound state. Our code solves
for the nonlinearly deformed profile. The question is whether the deformation
matters where GW emission matters. It appears not to, for a scaling reason:

- In the code's units, the field is measured in units of $f_a$ and the cloud
  mass fraction is $M_c/M = \kappa E$, with $\kappa \equiv (f_a/M_{\rm Pl})^2$
  (code: `coupling_to_hole`) and $E$ the cloud energy in code units. The GW
  loss in code units is $\kappa\,\tilde P\,E^2$: it is proportional to
  $\kappa$, while the self-interaction losses do not carry that factor.
- So GW emission matters at large $f_a$, where a given mass fraction
  corresponds to a small field in units of $f_a$. In the $5.7\times10^{17}$ GeV
  run the peak field is roughly $|\Phi| \approx 0.02$ (scaled from the peak
  field scan), and the leading nonlinear correction to the potential is of
  relative size $\Phi^2/12 \approx 4\times10^{-5}$.
- At $3.7\times10^{16}$ GeV the field is larger ($|\Phi| \approx 0.16$), but
  there GW emission is negligible anyway.

This argument is an estimate, not a test. The decisive check is listed in the
last section.

One limitation to record: SuperRad evaluates the coefficient at the spin where
the level is exactly at the superradiant threshold and ignores the spin
argument. Its authors quote an error of at most about 24% from this for the
scalar $m = 1$ level. Our clouds sit near threshold when GW emission matters
(threshold spin 0.87 at $\alpha = 0.3$), so this is acceptable.

## What the tabulated-coefficient option involves

Physics. For a GW line built from $q_1$ quanta of level 1 and $q_2$ quanta of
level 2 (annihilation of $|211\rangle$: $q = (2, 0)$; transition:
$q = (-1, +1)$), with GW frequency $\Omega = q_1\omega_1 + q_2\omega_2$ and
power $P$:

- level $i$ loses quanta at the rate $q_i P/\Omega$;
- the cloud loses energy $P$ and angular momentum
  $(q_1 m_1 + q_2 m_2)\,P/\Omega$;
- nothing is added to the black-hole equations (the radiation goes to
  infinity). Brito et al. (1706.06311) and 2611.15865 do the same.

The evolution on this branch already balances the number of quanta in each
level against loss rates, and re-solves the cloud profile at every step, so the
GW term is one more contribution to the existing loss vector, with $\alpha$
taken from the current black-hole mass. The earlier proposal in
`docs/further_work_1006.md` (section 6) has the same structure.

Steps, each with its check:

1. **Coefficient table.** Export $\tilde P(\alpha)$ for $m = 1$ and $m = 2$
   from the SuperRad Python package to a small data file, so the Julia code
   does not depend on Python. Check: reproduce $8.4\times10^{-11}$ at
   $\alpha = 0.3$, and the small-$\alpha$ fit
   $0.01093\,\alpha^{14} - 0.02901\,\alpha^{15}$.
2. **Profile-free test first.** Add the term to the two existing scripts that
   evolve occupation numbers with linear rates and no profile solve
   (`scripts/leaver_spindown.py`, `scripts/hydrogenic_spindown.jl`). Check
   with an exact answer: a single $|211\rangle$ cloud with the superradiant
   rate switched off must follow $M_c(0)/(1 + t/\tau)$. Then compare the free
   case with the second drop in Baryakhtar et al. Fig. 9, which
   `docs/session-2026-11-12.md` already identified as the feature our model is
   missing; expect agreement in shape, not in timing, because their
   coefficient is about 11 times larger.
3. **Main evolution.** Add the same term to the loss vector in
   `src/evolution.jl`, and add the GW power and torque to the diagnostics
   output. Checks: total energy and angular momentum (cloud + hole + scalar
   radiation + GW) are conserved to the existing tolerance; setting the
   coefficient to zero reproduces the existing runs exactly; the result at
   small $f_a$ is unchanged.
4. **Transition term**, with the multipole coefficient
   ($P \simeq 3.6\times10^{-7}\,\alpha^{12} M_{c,211} M_{c,322}/M^2$,
   $\Omega = \tfrac{5}{72}\alpha^2\mu$, from Baryakhtar et al.). Check: quanta
   lost by $|322\rangle$ equal quanta gained by $|211\rangle$.

Estimated effort: a few days including the checks. This adds no measurable
cost per step.

## What the full relativistic option would involve

For comparison, computing the GW flux from our own profiles requires:

1. Spin-weight $-2$ spheroidal harmonics at the GW frequency. Our Leaver code
   handles spin weight 0 only.
2. The homogeneous spin-2 radial solution that is ingoing at the horizon, and
   its amplitude at infinity. Direct integration of the Teukolsky equation is
   badly conditioned at large radius, so the standard practice is to transform
   to the Sasaki–Nakamura equation first, or to use the series method of Mano,
   Suzuki & Takasugi. The series machinery in `src/series.jl` is the right
   kind of tool but is written for the scalar operator.
3. The source: the parts of the scalar stress-energy that oscillate at
   $2\omega_1$, $\omega_1 \pm \omega_2$, and so on, projected on the Kinnersley
   tetrad and acted on by Teukolsky's derivative operators. The code currently
   computes only time-averaged bilinears of the field. This is long but
   mechanical algebra.
4. The radial integral of (2) against (3). The cloud spans about four GW
   wavelengths and the integrand oscillates with strong cancellation, which is
   the known difficulty of this calculation.

Validation targets exist: the SuperRad value at $\alpha = 0.3$ with a linear
profile, and the Yoshino–Kodama flat-space limit at small $\alpha$.

Estimated effort: four to eight weeks at this project's incremental pace. That
is an estimate, not based on a prototype. What it would buy: the effect of
profile deformation on the annihilation power, relativistic transition and
cross-annihilation rates, the spin dependence away from threshold, and the
horizon GW flux.

## When to revisit the relativistic option

Revisit it if any of these turns out to be true after the tabulated term is in:

- The spin history or final cloud composition changes noticeably when the
  $|211\rangle$ coefficient is varied by its stated uncertainty (about 25%).
  This is the cheap decisive test: one extra run with the coefficient scaled
  by 1.25.
- A regime appears where GW emission matters and the field is not small in
  units of $f_a$ (for example much larger $\alpha$, where the coefficient
  grows as $\alpha^{14}$).
- Levels with $m \ge 3$ are added. SuperRad has no relativistic scalar
  coefficients there, and 2611.15865 reports that GW transitions can then
  decide which higher levels grow.
- The GW signal itself, not only its backreaction, becomes a result we want
  to report.

## Sources

- Omiya, Takahashi, Tanaka & Yoshino, arXiv 2404.16265 (full text read).
- Cannizzaro, Spieksma & Witte, arXiv 2611.15865 (HTML version read).
- Siemonsen, May & East, arXiv 2211.03845, and the `superrad` package v0.3.0
  (run to obtain the coefficients quoted here).
- Yoshino & Kodama, arXiv 1312.2326; Brito, Cardoso & Pani, arXiv 1411.0686;
  Brito et al., arXiv 1706.06311; Arvanitaki, Baryakhtar & Huang,
  arXiv 1411.2263; Baryakhtar et al., arXiv 2011.11646.
- Repo: `src/evolution.jl`, `src/charges.jl`, `docs/further_work_1006.md`
  section 6, `docs/session-2026-11-10.md`, `docs/session-2026-11-12.md`, and
  the spindown diagnostics under `results/`.
