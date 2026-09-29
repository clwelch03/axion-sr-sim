# Multimode axion evolution: literature and candidate research gaps

Search date: 2026-09-28. This is a targeted review of published papers and
preprints, not a complete citation-database export. Searches covered the two
proposed starting papers, related references and later work, and combinations
of axion/Kerr clouds with multimode evolution, nonlinear profiles, phase
dynamics, and adiabatic evolution. Publication/version dates were checked where
recent work affects the assessment; search-engine crawl dates are not paper
dates. Negative findings below mean no matching calculation was identified
in this review, not proof of priority over unpublished or unindexed work.

**Assessment:** the strongest candidate is a controlled adiabatic treatment of
several mutually deforming clouds of one real axion field on Kerr, together
with a quantitative test of when existing occupation-rate descriptions fail.
Four or five clouds is a computational target, not by itself a novelty claim.

## Closest prior work

| Paper | Scope checked and consequence for our project |
|---|---|
| [Yoshino & Kodama, *Bosenova and Axiverse* (2015), 1505.00714](https://arxiv.org/html/1505.00714#S4) | Section 4 includes direct nonlinear sine-Gordon evolution of superposed angular modes and radial overtones. Nonlinear multimode field simulations already exist. These are useful dynamical benchmarks, not the proposed general slow-evolution construction. |
| [Baryakhtar, Galanis, Lasenby & Simon, *Black hole superradiance of self-interacting scalar fields* (2021), 2011.11646](https://arxiv.org/html/2011.11646) | Quartic mode transfer, scalar emission, saturation, higher levels, and observational consequences; Appendix B derives perturbative frequency/rate corrections. Neither a channel network nor demonstrating that interactions populate other modes is new. |
| [Omiya et al., *Impact of multiple modes…* (2022 preprint, 2023 journal), 2211.01949](https://arxiv.org/html/2211.01949v1#S3) | Section III treats one nonlinear primary and a secondary on that background, neglecting the secondary's deformation of the primary. Section IV already examines additional modes perturbatively. The extension is mutual nonlinear response and controlled evolution, not the first test of a third cloud. |
| [Witte & Mummery, *Stepping Up Superradiance Constraints on Axions* (2024 preprint, 2025 journal), 2412.03655v3](https://arxiv.org/html/2412.03655v3#S2.SS4) | Relativistic quartic scattering and coupled occupations for relevant states through principal number n=5, with BH evolution. Section II.4.4 discusses higher-level truncation. Relativistic radial calculations do not make the scalar interaction treatment nonperturbative in cloud amplitude. |
| [Takahashi, *Relativistic frequency shifts in gravitational waves from axion clouds* (2026), 2604.21239v2](https://arxiv.org/html/2604.21239v2) | Sections III–IV compute perturbative self-interaction and self-gravity frequency shifts, including cross-mode contributions. Published August 2026. Computing the first relativistic frequency shift due to another occupied cloud is therefore not an available novelty claim. |

The two-amplitude conservation matrix in the design note is an implementation
observation, not a newly discovered research gap. Existing rate equations
already provide individual mode-transfer laws. Our unresolved task is deriving
their consistent finite-amplitude replacement when the basis profiles depend
on several amplitudes, and validating the approximation.

## Other work that narrows possible claims

| Paper | Relevance and boundary of overlap |
|---|---|
| [Chia et al., *Self-interacting gravitational atoms in the strong-gravity regime* (2023), 2212.11948](https://arxiv.org/html/2212.11948) | Nonperturbative quasi-equilibrium complex-scalar configurations with quartic interactions and metric backreaction, neglecting BH spin. Important prior work for nonlinear profiles; different from evolving several real-field frequencies on Kerr. |
| [Collaviti et al., *Observational prospects…* (2024 preprint, CQG 2025), 2407.04304](https://arxiv.org/abs/2407.04304) | Self-interacting cloud evolution applied to GW signals, frequency drift, and detector reach. Adding observational forecasts alone does not establish a new dynamical treatment. Abstract-level screening here. |
| [Takahashi, Omiya & Tanaka, *Self-interacting axion clouds…in binary systems* (2024), 2408.08349](https://arxiv.org/html/2408.08349) | Self-interaction and tidal dissipation, transitions, and binary evolution. A binary extension would need a separate novelty argument. |
| [Xie & Huang, *Self-interaction effects on the Kerr black hole superradiance…* (2025), 2503.10347](https://arxiv.org/html/2503.10347#S2.SS3) | Sections II.3–II.4 use a perturbative effective-mass prescription, assume dominance by one level, and explicitly exclude effects from populated higher levels. Relevant comparator for frequency and growth-rate corrections, not a completed mutual-profile evolution. |
| [Caputo, Franciolini & Witte, *Superradiance Constraints from GW231123* (2025 preprint, PRD 2026), 2507.21788](https://arxiv.org/html/2507.21788#S2) | Applies multimode occupation and BH evolution to spin constraints. Such evolution is already used in current phenomenology. |
| [Zhu, Piao & Zhang, *Black Hole Superradiance of Interacting Multi-Field* (2025), 2509.10892](https://arxiv.org/abs/2509.10892) | Two different scalar species with a cubic coupling. Multiple fields and multiple bound modes of one field are distinct problems. Abstract-level screening. |
| [Lyu, Cai & Liu, *Self-stimulated avalanche in superradiant axion clouds* (2026), 2606.10776v3](https://arxiv.org/html/2606.10776v3) | Coherent two-level transitions seeded by external GWs with radiation feedback. Version 3, dated September 16, substantially revises the formulation and title. Broad claims that coherent cloud transitions are unexplored would be wrong; our candidate concerns intrinsic scalar self-interaction. |
| [Mennenga, *Quartic Scalar Clouds on Fixed Kerr Backgrounds* (July 2026), 2607.28044](https://arxiv.org/abs/2607.28044) | Stationary synchronized complex-field Q-clouds, including an axion-inspired potential comparison. This does not establish dissipative real-field multimode evolution. Abstract-level screening. |
| [Li et al., *Dynamics and Frequency Conversion of Accreting Axion Clouds* (August 31, 2026), 2609.00128](https://arxiv.org/html/2609.00128) | Nonlinear accretion-grown clouds and potential-dependent harmonic emission. Strong-field potential dependence and collapse/saturation are active topics, but the driving mechanism differs from Kerr superradiance. |

Also screened [*Crowdsourcing Gravitational Waves from Superradiant Axions*
(March 2026), 2603.15734](https://arxiv.org/html/2603.15734): its evolution
calculation explicitly neglects axion self-interactions. It does not fill the
present gap.

## Candidate questions, ordered by fit to this repository

### 1. Does mutual nonlinear deformation change the evolution?

**Best fit; apparent gap, with moderate confidence in the novelty assessment.**
Construct a family in which the profiles and frequencies of multiple clouds
respond to their joint amplitudes. Retain the cosine interaction within an
explicit, converged Fourier/angular truncation. Derive the slow evolution
from that family, including loss channels, rather than assuming unchanged
linear basis functions.

Compare three descriptions at matching physical initial conditions: a
perturbative rate model, the dominant-primary approximation, and the mutually
deformed calculation. Useful outputs are maximum cloud energy, equilibrium
composition, dominant outgoing flux, and instability rates of additional
clouds. BH spin-down is a later output requiring explicit BH evolution.

This can yield a useful result even if corrections are small: a demonstrated
domain of accuracy for commonly used rate models. Strong microscopic coupling
does not guarantee large local field amplitude; early saturation can keep a
cloud dilute. First establish that an accessible trajectory reaches a regime
where the extra nonlinear response matters.

### 2. Do amplitude shifts change a channel's kinematics or resonance?

**A concrete mechanism to test within the first question.** For an emitted
combination, evaluate

\[
\Omega_c(\boldsymbol A)=\omega_{i,R}(\boldsymbol A)
 +\omega_{j,R}(\boldsymbol A)-\omega_{k,R}(\boldsymbol A),\qquad
\delta\Omega_c=\delta\omega_i+\delta\omega_j-\delta\omega_k.
\]

For positive frequency, infinity propagation changes character at
\(\Omega_c=\mu\). A useful diagnostic is
\(|\delta\Omega_c|/|\Omega_c^{(0)}-\mu|\), away from a zero denominator.
Near threshold, evaluate the radial response directly; this ratio is a flag,
not an error estimate. Loss into the horizon can remain when infinity emission
closes. Separately compare detuning from bound resonances with their widths
and the evolution timescale.

The possible new result is feedback of changing occupations on thresholds,
overlaps, and subsequent evolution, beyond already studied relativistic
kinematics and perturbative shifts. No such effect is asserted to occur in
our parameter range without calculation.

### 3. Which apparent equilibria survive additional cloud perturbations?

**Good intermediate project if computed on a nonlinear background.** Linearize
the scalar dynamics around a continued background and calculate growth or
decay of candidate extra modes, retaining the needed frequency sidebands.
Compare this with the perturbative stability result before undertaking a
large evolution campaign. Stability of a multi-frequency solution may require
phase variables; a scalar occupation growth rate is not automatically enough.

The target would be a shifted instability boundary or an error bound on the
existing boundary, not the already established possibility of higher-mode
growth. A small eigenvalue is a useful place to test sensitivity, provided
its sign is resolved against numerical error.

### 4. How many modes are sufficient for a specified observable?

**Explicitly open in part of the literature, but not exclusive to our method.**
Witte–Mummery identify missing channels into m=6 states with n>=7 when
discussing m>=3 spin-down. Thus four or five chosen modes need not form a
closed or converged model. [Section II.4.4](https://arxiv.org/html/2412.03655v3#S2.SS4.SSS4)

Build the retained set from allowed sources and measured growth/transfer
rates. Test convergence of the desired flux, instability boundary, or spin
history as states are added. A stable low-mode population alone does not
bound effects from an omitted state that can grow later. This project can
start with existing perturbative rates independently of nonlinear profiles.

### 5. When is phase averaging invalid?

**More speculative; do not make this the central claim yet.** For a candidate
four-wave interaction check both angular selection and
\(\Delta\omega=\omega_i+\omega_j-\omega_k-\omega_l\). If the detuning is large
compared with relevant slow rates, the phase-sensitive contribution averages
away. If it is comparable, derive and evolve the resonant phase. Near-degenerate
hydrogenic energies alone do not establish a surviving Kerr resonance.

A phase-dependent outcome differing from a rate model would be interesting,
but neither phase locking nor oscillatory/chaotic exchange should be promised.
The neighboring GW/tidal literature already includes coherent dynamics.

## Recommended first research result

Keep the Section III reproduction as the validation milestone. Then measure
the leading omitted response of the primary to a finite secondary, compare
its effect on fluxes and frequencies, and locate a regime where that response
is significant or demonstrate that it stays small. Next compute the growth
of a candidate third cloud on the corrected background. Use this evidence to
choose the fourth and fifth clouds.

A possible project question is:

> When do mutually deformed axion clouds depart from the usual scattering-rate
> description, and does that departure change saturation or mode selection?

This requires a defined adiabatic closure, a weak-field recovery check, and
radial/angular/Fourier convergence. Selected time-domain comparisons would
help distinguish an approximation failure from a genuine instability.
Self-gravity, binary tides, and waveform forecasts are additional projects;
none is required for the first fixed-Kerr result.
