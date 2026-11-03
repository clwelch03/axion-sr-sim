# Project handoff

Updated 2026-11-02. This is the authoritative project status: what works,
open caveats, and next steps. Campaign READMEs own exact commands, saved-result
paths, and detailed validation records; design documents own equations and
implementation decisions.

Active direction: the 36-point QCD linear/Omiya/joint static-profile campaign
is complete. The next scientific comparison is time evolution: whether mutual
deformation changes peak cloud populations, saturation timing, and two-cloud
quasi-equilibrium, including evolution of BH mass and spin. This next campaign
has not been implemented or launched. New joint evolution uses the leading-order
chiral QCD potential; retain cosine for the faithful Omiya reproduction.
The Omiya, Takahashi, Tanaka & Yoshino reproduction
([2211.01949v1](docs/papers/2211.01949v1.pdf), Sec. III) and existing limited
evolution prototypes supply starting points, not a validated three-model QCD
evolution comparison. See the
[multiphase campaign](runs/adiabatic/multiphase_two_cloud/README.md) for current
outputs and [adiabatic-repro.md](docs/adiabatic-repro.md) for the primary input.

**November 3: fresh-checkout validation before the Perlmutter move.**
All **1,399 Julia unit checks** passed in Slurm job **27537972**: 17 loader,
508 core, 257 single-cloud, 62 primary/secondary, and 555 multiphase checks.
The stale test for the standalone AD diagnostic retired in `5698692` was
removed; current source/response derivative regressions remain in the suite.
All 72 Python tools tests and 23 Dolan tests also passed. The environment
export, NERSC wrapper and transfer procedure are in the
[Perlmutter handoff](docs/perlmutter-handoff.md); remote execution remains
subject to the status recorded with the transfer artifacts.

**November 2: accelerated chiral continuation and four-model plot.**
The user clarified that they intended the faster run for chiral, and requested
switching its remaining steps. Job 27508586 resumes after step 970 under
`chiral_growth_averaged_resume_20261102_v1`. The unchanged inherited history,
original time grid and new rates were checked before cancelling old job
27501595. The continuation completed successfully in 7m03s, preserving the
original time grid at all 1,001 samples. The completed four-model plot uses
this continuation, the accelerated cosine run, the faithful Omiya control,
and the user's chosen independent exponential linear control. Chiral ends
0.41356% (211) and 0.38214% (322) above linear in outer amplitude. See the
[continuation record](runs/adiabatic/multiphase_two_cloud/README.md#november-2-chiral-continuation-and-four-model-comparison)
for source provenance, handoff checks, PNG/PDF plots and exported linear/combined
trajectories (`growth_four_models_inner80_20261102_v1`). The old chiral report
may still say `running` despite cancellation; use the completed new prefix.

**November 2: cosine speedup run 27507601 completed in 49m01s.**
The user requested a complete trajectory with the averaged evaluator. It
reuses the previous inner80 seeds, numerical settings and time grid, with
new snapshot/output prefix `cosine_growth_averaged_inner80_20261102_v1`.
The snapshot verifies all 264 files. The previous direct cosine job 27501596
completed in 2h51m47s; other trajectories were left intact. See the
[run record](runs/adiabatic/multiphase_two_cloud/README.md#november-2-cosine-evolution-with-analytic-phase-averaging)
for the final comparison and plots. Slurm confirms COMPLETED, exit 0:0;
all 1,000 steps reached the same final time. The whole-run speedup is 3.50x,
including startup, checkpoints and plotting. Maximum relative amplitude
differences are 6.67e-15 (211) and 1.88e-10 (322); maximum rate differences
are 1.17e-13 and 1.44e-7. The final E/J budget errors are about 1.735e-5.

**November 2: analytic phase averaging implemented for joint responses.**
`multiphase/averaged_observables.jl` is now the default observable path for
fresh response/evolution samples, with `observable_method=:direct` preserving
the original tensor/phase sums. The nonlinear potential retains its quadrature;
matching, AD chunk sizes, balance equations and Omiya physics are unchanged.
No running trajectory was restarted for this change. All 555 available
multiphase checks (82 new) and 17 loader checks passed; the full runner has a
then-pre-existing missing include, documented with the excluded test and saved logs
in the [benchmark record](runs/adiabatic/multiphase_two_cloud/README.md#november-2-analytic-phase-averaging).
Both benchmark jobs completed, with all eight domain/potential/profile cases
passing. On [-80,180], complete steps are 1.62–1.88x faster for chiral and
3.56–4.63x for cosine; the largest relative rate difference is 3.65e-7 in the
weak late 322 cloud. The record owns detailed timings and frozen-source
provenance. The kernel accepts longer phase labels, but the profile solver
and E/J closure still describe two independent clouds.

**November 2: outer boundary changed to 180 for fresh growth runs.**
The user selected `r_*^out/M=180` after the isolated boundary checks. The
chiral/joint-cosine growth driver now overrides the old seed's outer boundary
in both linear initialization and joint matching. Following the user's
clarification, the Omiya control retains cosine and uses its own [-100,100]
domain (paper p. 11 at alpha=0.42), with the existing far-zone construction of
the outgoing m=3 solution. It matches the joint seeds/background/end time;
the reports distinguish the spatial settings. Other numerical defaults remain
as before. Fresh source snapshots/output prefixes and a new chiral reference
report are required for the matched joint-cosine control. See the [rerun recipe](runs/adiabatic/multiphase_two_cloud/README.md#november-2-outer-boundary-default-for-fresh-growth-runs).

**Fresh three-run submission:** chiral **27498870**, joint cosine **27498913**,
and Omiya cosine **27498914**, following that recipe at the user's request.
New output prefixes in the multiphase campaign are respectively
`results/chiral_growth_outer180_20261102_v1`,
`results/cosine_growth_domains_20261102_v1`, and
`results/omiya_growth_domains_20261102_v1`. The new common end time is
t/M=6.649525029285335e7, with unchanged seeds and 1,000 output intervals.
Snapshots `work/chiral_growth_outer180_20261102_v1` and
`work/growth_controls_domains_20261102_v1` have 260/261 verified files, and
their source/dependencies/seed match. Chiral passed its initial linear-rate
check; use the [new run record](runs/adiabatic/multiphase_two_cloud/README.md#november-2-fresh-three-run-evolution-submission)
and live reports for progress. The earlier joint runs both stopped at step
794 on the log-step guard; their apparent 322 depletion is not a validated
physical result. The separate validation jobs were left as they were.

**November 2: exploratory chiral evolution submitted as Slurm 27495927.**
The user selected the existing fixed-BH alpha=0.42/spin=0.99 pilot settings,
100-fold smaller outer amplitudes (A211=1e-6, A322=1e-7), and ten initial 211
amplitude e-folding times. This uses 1,000 corrected midpoint steps and fresh
chiral initial profiles, with the old finite domain/support. It is separate
from the eventual physical QCD/BH-spin-down comparison. Snapshot:
`runs/adiabatic/multiphase_two_cloud/work/chiral_growth_20261102_v1` (257 verified
files); output prefix: `results/chiral_growth_20261102_v1` in that campaign.
The job requests one CPU, 8 GiB, and eight hours; accepted states are checkpointed.
Check the [campaign record](runs/adiabatic/multiphase_two_cloud/README.md#november-2-chiral-growth-pilot)
and Slurm/logs for its outcome before claiming a completed trajectory.

**November 2: matched cosine controls submitted at the user's request.**
Jobs **27496847** (joint cosine) and **27496848** (Omiya cosine) share the chiral
seeds and exact end time t/M=6.646859783614447e7. Their snapshot is
`work/growth_controls_20261102_v1` (261 verified files); output prefixes in the
multiphase campaign are `results/cosine_growth_20261102_v1` and
`results/omiya_growth_20261102_v1`. Joint cosine uses the same 1,000 midpoint
steps; Omiya first builds a 213-point profile/flux table, then integrates the
unchanged energy/angular-momentum equations with interpolation/time checks.
This Omiya control uses the pilot's finite domain and harmonic support, as in
the earlier October 19 comparison; it is not a full-domain paper reproduction.
No Omiya physics routines or running chiral snapshot were changed. All 85
focused configuration/evolution checks passed. See the
[control record](runs/adiabatic/multiphase_two_cloud/README.md#november-2-matched-cosine-evolution-controls)
and job reports for outcomes; do not assume completion from submission.

The user added `:chiral` to the shared potential value and remainder with
mu/md=2.2/4.79 (c4=0.352968986964824). All 159 focused potential/source checks
pass, including high-precision comparisons down to field 1e-60. Restoring the
computed `amplitude_rates` return field reconciles the existing balance tests;
all 47 evolution checks pass, and replaying saved October 15 balance inputs
reproduces their rates exactly. No Omiya physics routines were changed.
Keep the work collaborative: explain numerical experiments and choices before
starting them and bring outcomes back before expanding the scope.

Read first: [multimode-adiabatic-design.md](docs/multimode-adiabatic-design.md)
(equations, decisions, component map) and the
[secondary campaign README](runs/adiabatic/multimode_secondary/README.md)
(commands, results, checks). Cite equation numbers from the PDF in
`docs/papers/`, not the arXiv HTML.

## Where things stand

### October 30: all 36 static comparisons complete; evolution is the next priority

Verified from the saved collection: **36/36 cases are plotted**, comprising
24 initially-211 and 12 initially-322 cases. The grid is
alpha = 0.2, 0.3, 0.4; initial spin = 0.8, 0.9, 0.95, 0.99; and
BH mass = 5, 20, 100 solar masses. All 36 PDF pages are present, and all 108
source-file hashes in the collection metadata match the saved files.

- [All 36 profile comparisons](runs/adiabatic/multiphase_two_cloud/results/qcd_campaign_plots_20261030_v1/all_36_profiles.pdf).
- [Per-point manifest and measured differences](runs/adiabatic/multiphase_two_cloud/results/qcd_campaign_plots_20261030_v1/comparison_manifest.csv).
- [Collection metadata and source hashes](runs/adiabatic/multiphase_two_cloud/results/qcd_campaign_plots_20261030_v1/metadata.json).

Each figure contains both clouds and all three models, with the **same energy
for each corresponding cloud across models**. The retained population targets
are the original first-saturation estimates; the nonlinear accepted energies
set the exact comparison energies. No arbitrary macroscopic second population
was substituted. The full cosine potential remains the user's choice.

The joint curves use the **nonlinear dominant cloud plus coupled
small-secondary limit**, not a finite-population two-cloud nonlinear solution.
This is appropriate to the retained microscopic secondary seeds, but does not
test mutual deformation when both clouds carry substantial populations.
Maximum Omiya/linear magnitude differences, divided by each cloud's linear
peak, are 0.102611% for 211 and 0.229732% for 322. Maximum joint/Omiya
differences are 0.0000157748% and 0.0000756882%, respectively. These are saved
grid measurements, not spatially converged physical effects. The plotted
recovery-discrepancy bands measure implementation consistency, not spatial
error. Small or unresolved shape differences do not establish negligible
changes in fluxes or accumulated evolution.

**Interpretation agreed with the user:** these static results are inconclusive
about dynamical importance. Further searches for visibly deformed snapshots
or debates over prescribed populations are no longer the main priority.
Compare evolution from the same physical initial conditions and let the
populations be outcomes. Measure peak occupations/energies, saturation times,
whether and when both clouds coexist, quasi-equilibrium, and BH spin/mass
histories. Equal-energy matching remains a requirement for diagnostic profile
snapshots; do not force equal populations between independently evolved models.

**Omiya/QCD clarification:** Omiya's cosine framework can also approximate the
QCD axion, but this paper does not impose the QCD mass–decay-constant relation.
The dimensionless figure scales out the decay constant; it does not itself
choose a physical coupling. The string-axion example fa/Mpl ~ 1e-3 appears in
the justification for neglecting BH spin-down (PDF p.26, Eqs.95–96).
Our stellar-mass QCD choices force much larger fa, so those large dimensionless
energies are inaccessible through the available spin budget. This does not
make their evolution wrong or exclude large deformation for all QCD axions.
For example, imposing our QCD relation at fa = 1.22e16 GeV gives
mu = 4.67e-10 eV and, at alpha = 0.3, BH mass ~0.086 solar masses; this is our
physical conversion, not a mass assigned in their figure. The saved high-energy
reproduction remains a useful solver benchmark.

No new three-model QCD evolution has been run. The older
[QCD spin-down pilot](runs/adiabatic/qcd_spindown/README.md) already evolves BH
mass/spin with an Omiya-like scalar model, but uses one trajectory-specific
table and a nonrelativistic, noninteracting comparator. Its timing difference
primarily reflects different linear growth rates and does not isolate mutual
deformation. Existing joint short-evolution checkpoints also retain the
slow-growth closure and convergence caveats below. Use these as starting
material, not evidence that the proposed comparison is complete. The active
implementation sequence is in [Next steps](#next-steps); older dated next-step
records below are historical.

The new [QCD two-cloud evolution design](docs/qcd-two-cloud-evolution-design.md)
maps the reusable Omiya/joint code, defines the three-model comparison and
fixed-reference units, and separates implemented fixed-BH evolution from the
joint BH-response/slow-growth work still needed. It also specifies the
small-secondary transition, cost strategy, outputs, and incremental checks.
This is a design document; no new evolution code or job accompanies it.

### October 30: first linear/Omiya/joint QCD comparison completed

Slurm **27438309** completed the missing joint calculation for
alpha=0.3/spin=0.99/M=5 solar masses. The three-model, two-cloud figure is
`results/qcd_joint_probe_20261030_v2/both_clouds.png` in the multiphase campaign.
It uses the **small-secondary limit of the joint cosine equations**, explicitly
labeled on the figure: an independently integrated nonlinear primary plus
the coupled 322/(2,-1) response, including m=0 feedback into 322. This is not a
finite-population two-cloud nonlinear solve. At the retained microscopic seed,
higher powers of the secondary amplitude are far below floating precision.
The (-1,2) component starts at second order in that amplitude and is omitted
from this limit. Existing full finite-population joint benchmarks remain separate.

Both clouds are plotted at the same respective energies across all three
models. Primary energy agrees to printed precision; secondary energy agrees
within 1e-12 relative. The independently integrated joint primary coincides
with Omiya to roundoff. For 322, joint versus linear differs by 0.2298074% of
the linear peak, and joint versus Omiya by 0.0000756882%. Spatial convergence
of that small latter difference is not established. Joint unit-probe matching
is 1.99e-8, accepted using the measured frequency-precision floor.

Checks passed against full-source finite differences (errors 2.90e-7 then
2.90e-9), a denser second-phase grid (5.26e-16), the Omiya diagonal source
(3.76e-16), and Omiya profile recovery after disabling the return coupling
(8.29e-11). The first attempt, 27438292, passed those source checks but failed
before Newton due to a driver variable-scoping error; v2 fixes it.

The user asked why Omiya Fig. 7 has much larger deviations. The physical unit
conversion is decisive: E/M_BH = [E/(Fa^2 M)]*(fa/Mpl)^2, with unreduced Mpl.
Our fa=7.109e17 GeV gives a factor 0.00339048; Ecode=660 would be **2.24 BH
masses per cloud**, whereas the accepted QCD pilot is Ecode=9.32745 and
E/M_BH=0.0316246. Omiya fixes the Kerr background and neglects BH spin-down
and GW emission, explicitly justifying the former for a string-axion example
fa/Mpl~1e-3 (PDF p.26, Eqs.95–96). There, Ecode=660 corresponds to only
0.00066 of the BH mass. This is a different coupling regime, not evidence
that Omiya's evolution is wrong. Our saved high-energy reference at
Ecode=650.5419 already has 8.16%/19.05% linear–Omiya profile deviations for
211/322, so the profile solver does reproduce large deformation at those inputs.

All 23 additional initially-211 cases have produced passing three-model plots
in array **27438516_1-23%3** (24 completed comparisons including the first pilot),
using the frozen complete pipeline `work/qcd_211_comparisons_20261030_v2/`:
primary, Omiya probe, coupled joint probe, energy checks, and plots. Outputs
go to the matching `results/` directory, one named subdirectory per case.
The first v1 check (27438377_1) completed the primary/Omiya pair but failed an
overstrict 1e-10 cross-implementation recovery gate. Diagnostic 27438465
measured the discrepancy as 2.01827e-9 on the longer alpha=.2 domain; the
primary implementations agree to 1.00e-12 and the source identities to ~1e-15.
The v2 gate is 1e-8 in profile scale, with the measured discrepancy saved;
joint/Omiya differences below ten times that discrepancy are labeled unresolved
and the measured discrepancy band is plotted. Spatial convergence remains
unestablished independently of this consistency check.

Focused retry **27438529 succeeded** for the dominant-322 alpha=.4/spin=.95/M=5
pilot after changing the shooting match from 5 to 25 and rematching its free
reference. Three Newton corrections give normalized residual 6.79e-12 and
energy 26.59159 (1.68% below the approximate target); the dominant profile
differs from its equal-energy linear curve by 0.146626% of the linear peak.
Snapshot/output name: `qcd322_match025_20261030_v1`.

The reversed joint-probe pipeline is now implemented in campaign drivers:
cloud order [322,211], with the `(2,-1)` response now an m=3 outgoing field,
including its return coupling into 211. The 12 cases are frozen in
`work/qcd_322_comparisons_20261030_v1/`. Their production match point is 24,
adjacent to the successful 25 test and compatible with the primary step .2
and probe step .4. Both free references are rematched on those grids.
First full reversed comparison **27438650_1** and the other 11 cases in
**27438706_2-12%3** have all produced passing three-model plots. The latter
array was dependent on the first case's success. The same source, matching,
energy, and profile consistency checks apply in both dominance orders.

Collector **27438707**, submitted after the comparison arrays, has written the
complete 36-point manifest and 36-page `all_36_profiles.pdf` under
`results/qcd_campaign_plots_20261030_v1/`. Earlier partial bundles remain
historical outputs. This static campaign introduced no time evolution or
population/potential change.

### October 30: both clouds restored to the QCD preview

The user explicitly requires **both clouds even when one is tiny**. Replacing
the minor population by zero and plotting only the 211 pilot omitted part of
the requested comparison. Slurm **27438012** completed the missing Omiya 322
probe on the accepted nonlinear 211 background. The corrected two-panel
preview is `results/qcd_both_clouds_20261030_v1/both_clouds.png` in the
multiphase campaign; its CSV/JLS/TOML retain both clouds and both models.
The joint model remains outstanding and is explicitly absent from this preview.

Each cloud uses its own physical vertical scale and has identical energy
across the linear/Omiya curves. E211=9.3274498169; E322=4.2693206198e-76,
recovered from the original nonzero seed diagnostic N322=1.0120597141.
That microscopic occupation is an illustrative probe normalization, not a
precise classical-cloud prediction. No arbitrary macroscopic population was
introduced. The older population diagnostic used c4=0.3; it is preserved,
while profile solves retain the chosen full cosine. The correction restores
the omitted mode, rather than revising the population prescription.

Maximum magnitude differences divided by each cloud's linear peak are
0.10261% (211) and 0.22973% (322). The secondary unit-amplitude matching norm
is 2.92e-8, below the measured-frequency-precision tolerance 1.24e-7; both
secondary energies agree with the target within 1e-12 relative. Spatial
convergence of these small differences is still unestablished. A tiny second
population suppresses its backreaction but does not eliminate its response
to the dominant background. Do not infer the second cloud's Omiya/joint
comparison from the dominant cloud's single-cloud limit.

### October 30: visibility of QCD profile differences assessed

The accepted 5-solar-mass alpha=0.3/spin=0.99 211 pilot has maximum magnitude
difference from the equal-energy linear profile **0.10261% of the linear peak**.
Its profile panel nearly overlaps; a percent-difference panel resolves the
shape change. The preview and 36-point ranking are in multiphase campaign
`results/qcd_visibility_20261030_v1/`; the figure contains linear/Omiya curves
only, because the independent joint limit is not yet computed.

Ranking the linear-field quartic potential energy against binding energy puts
the strongest candidates at 5 solar masses: alpha=0.4/spin=0.95 (322),
alpha=0.4/spin=0.9 (322), and alpha=0.3/spin=0.99 (211). The 322 points are
estimates, not accepted nonlinear solutions; the spin-0.95 direct solve failed.
The ranking is a weak-field diagnostic, not a forecast of exact profile
differences. At fixed alpha/spin it decreases by factors 16 and 400 for 20 and
100 solar masses relative to 5. None is expected to yield three visibly
separated dominant-cloud curves at the current negligible-other-cloud
populations: Omiya and joint should coincide in that limit. Angular-average,
mass-scaling, and radial-quadrature checks passed; spatial convergence of the
small measured nonlinear difference remains separate.

At the user's explicit request, a subagent investigated a regime with a
non-negligible subdominant cloud. Its nearest extension is initial
alpha=0.2/spin=0.99/M=3 solar masses, mu=8.91e-12 eV and fa=6.40e17 GeV,
using the current cosine coupling and QCD relation. A separate read-only
evaluation of the screen confirms optimistic catch-up 1.83 years versus a
4.92-year GW half-life, with analytic declining-pump gain 1.35 times the
required gain. This is outside the 36-point grid and is only a candidate for
later coexistence; no new populations or runs were selected. More accurate
rates and competing modes need checking, and coexistence would not guarantee
large profile deformation.

### October 30: all 36 points screened for later subdominant-cloud growth

At the user's request, `scripts/check_qcd_catchup.py` evaluates the growth-versus-
depletion condition without integrating coupled population or radial evolution.
Results are in multiphase campaign `results/qcd_catchup_20261030_v1/`, including
the 36-row screen, 576 assumption variants, source/input hashes, and preserved
drivers. **0/36 pass the nominal optimistic timescale test.** For 24 initially
211-dominated points, 322 growth is too slow; for the other 12, 211 remains net
damped after 322 saturation. The latter are alpha=0.3/spin=0.8 and alpha=0.4
at spin=0.8,0.9,0.95, at every mass.

The calculation uses the chosen cosine quartic coefficient c4=1 (f=fa), a
mass-changing first-saturation endpoint, and occupations/rates consistently
in units of the final BH mass. Compare the logarithmic build-up time, including
signed subdominant horizon growth/decay and interaction pumping, with the
dominant cloud's GW-annihilation half-life. The pump is held at its peak for
this optimistic screen, following the logic of Baryakhtar et al. Eq. (42).
An analytic check allowing the pump to decrease during GW depletion also has
no nominal passes. No BH-age cut is assumed.

The five closest cases are all 5 solar masses: alpha=0.2/spin=0.9,0.95,0.99,
and alpha=0.3/spin=0.95,0.99. Their optimistic build-up/depletion ratios are
4.72,3.54,2.89,6.52,3.66, respectively. A tenfold increase in growth relative
to GW depletion moves these five across the simple timescale threshold (four
also cross the analytic half-life check). These are sensitivity flags, not
coexistence detections or confidence intervals. None passes when varying only
the endpoint approximation, reference seed 1/1e6, target fraction 1/0.1, and
binding correction. Small-alpha rates at this grid's alpha remain preliminary.
Four tests passed, including analytic gain versus independent quadrature,
post-extraction occupation/angular-momentum units, signed-rate conventions,
and quartic-coupling scaling. Existing profile population targets are unchanged.

### October 30: three-model comparison and equal-energy requirement reaffirmed

The user requires **linear, Omiya, and joint fully nonlinear profiles at each
QCD point**, with the same cloud energy across the three curves on each plot.
The 10% population-target window is not an allowed mismatch between models.
For two populated clouds, match the same pair of per-cloud comparison energies
across models; the two clouds need not have equal energies to each other. Retain
the established comparison-energy convention and report total joint-field
energy separately. Nonlinear profiles cannot be rescaled after solving.

The current pilot driver calls the single-cloud nonlinear solver: it supplies
the Omiya dominant-cloud profile and a free linear reference, not the joint
calculation or the Omiya other-mode probe. Its historical CSV label `nonlinear`
denotes that single-cloud profile and must not be plotted as the joint model.
To avoid repeated expensive joint energy targeting, the proposed comparison
recipe is a direct joint solve near the population estimate, followed by Omiya
reference solves matched to its measured energies and exact scaling of the
free linear references. Equal energy must be verified before a three-way plot
is labeled as such.

The current first-saturation population estimate leaves the other cloud
classically negligible. At zero other-cloud population the dominant Omiya and
joint equations should reduce to the same single-cloud problem at matched
settings. The existing joint solve wrapper requires two positive amplitudes;
do not insert an arbitrary second population to bypass this or silently label
an Omiya result as an independently computed joint curve. The user chose to
retain the current populations for now, and subsequently clarified that this
requires retaining and plotting the minor mode as well. The zero target in
the original population CSV was an implementation shortcut, not permission
to omit that mode. Use the retained finite seed diagnostic (with its stated
limitations) or a clearly labeled normalized probe representation. Independent
joint profiles including that response are not yet implemented or submitted.

The user also asked whether the subdominant cloud is always negligible and
what investigating a later two-cloud stage would entail. The first-peak result
does not establish negligibility at later times. A later-stage estimate would
first compare growth/transfer and depletion times and check Baryakhtar et al.'s
regime assumptions, then use quasi-equilibrium formulas where applicable.
Outside those regimes, populations require an age/history prescription or a
population-only rate evolution; such a later-stage calculation has not been
performed. Retain the chosen cosine coupling consistently if applying those
rates to the profile campaign. No later-stage population change was requested.

### October 30: two direct QCD nonlinear pilot jobs submitted

Slurm array **27437158**, tasks 1–2, attempts the two 5-solar-mass pilots:
211 at alpha=0.3/spin=0.99, target code energy 9.4075951158; and 322 at
alpha=0.4/spin=0.95, target 27.0455406277. The 211 pilot exported a passing
result after three Newton corrections (46.6 seconds in the solve): measured
energy 9.3274498169, 0.852% below target, with matching residual 2.64e-12.
Its linear comparison energy is identical at printed precision. Slurm reports
211 COMPLETED in 1m07s. The 322 pilot FAILED in 5m49s after 15 Newton corrections,
with residual 1.79793e6; no accepted 322 profile was produced. The cause needs a
focused diagnostic, and no retry has been submitted. The user questioned
repeating the slow energy-targeting loop;
these jobs instead make one fixed-amplitude nonlinear solve each, with
A=sqrt(target_energy/unit-amplitude_linear_energy). There is no amplitude ramp
or repeated energy targeting. A 10% energy window is sufficient for the
approximate population prescription; the actual energy and an occupation
estimate are saved. The comparison linear profile is normalized exactly to
the measured nonlinear energy. The nonlinear profile itself is never rescaled.

Source/config/seeds are frozen in multiphase campaign
`work/qcd_profile_pilots_20261030_v1`; results go to the matching `results/`
directory, with logs `logs/qcd-pilot-27437158_TASK.out`. Each task uses one CPU,
16 GiB, and a two-hour limit. Syntax, driver loading, and the frozen manifest
passed. The jobs check the saved linear energy normalization, nonlinear
matching at `1e-8`, and the energy of the exported linear comparison. No library
code was changed. Pilot outcomes and spatial/angular convergence remain to be
checked; this step does not compute the other-mode probe or time evolution.

### October 30: retain the cosine potential for the QCD profile campaign

The user chose to keep the current full, untruncated cosine potential and its
existing normalization for the exploratory radial profiles. The QCD mass–decay
constant relation still sets the physical parameters; the cosine is the chosen
potential approximation. A chiral-potential comparison can be a later rerun and
is not a prerequisite for these profiles. The proposed next step is two nonlinear
dominant-cloud pilots at the estimated populations, one 211 and one 322, followed
by the 36-point campaign. Time evolution remains outside this campaign.

### October 29: QCD population estimates for all 36 points

The user accepted order-of-magnitude population accuracy and requested the
literature-based estimates. `scripts/estimate_qcd_populations.py` in the
multiphase campaign produced `results/qcd_populations_20261029_v1/` with a
36-row CSV, source/input hashes, assumptions, and a preserved driver copy.
No nonlinear profiles or time trajectories were submitted in this step.

The working stage is first dominant saturation, with occupation-one seeds and
no pre-existing classical cloud. Leading spin extraction sets N; saved linear
real frequencies set energy, and the QCD relation sets fa and code units.
Eight backgrounds (24 physical points) select 211; four (12 points) select 322:
alpha=0.3/spin=0.8 and alpha=0.4/spin=0.8,0.9,0.95. Dominant mass fractions
are 0.00523–0.0598, occupations 3.00e75–2.51e79, and code energies
0.00386–27.0. An algebraic mass-changing endpoint gives occupation corrections
of 1.30–7.87; the original backgrounds and leading populations remain the
working targets, explicitly treated as approximate saturation scales.

The linearized growing-phase transfer diagnostic leaves the other cloud
microscopic at all 36 points. This remains true using saved versus hydrogenic
linear rates, leading versus corrected populations, transfer multipliers
0.1/1/10, and seeds 1/1e6. This is a qualitative first-peak result, not a
factor-ten prediction of microscopic occupations or later two-cloud coexistence.
Zero minor-cloud target energy denotes the negligible classical limit.
The recommendation is one nonlinear dominant cloud plus a normalized linear
probe for the other mode. See the campaign README for the 12-background table,
equations, approximations, and next-stage interpretation. Six focused tests
passed, including independent energy-integral/angular-conservation checks;
all 36 rows, physical/code conversions, and input hashes were verified.

### October 29: five missing linear secondaries retried

At the user's request, array **27425079** submits the five remaining secondary
references: alpha=0.3/spin=0.9 and alpha=0.2 at all four spins. The separately
recovered alpha=0.3/spin=0.95 reference remains in
`results/omiya_spin095_linear_20261029_v2/` in the multiphase campaign.
The retry snapshot is `work/linear_grid_retry_20261029_v1`, with outputs in
`results/linear_grid_retry_20261029_v1` and logs `logs/linear-grid-27425079_TASK.out`.
Only secondaries are requested; previous results are preserved.

The opt-in driver first locates the root at absolute tolerance `1e-5`, then
measures the residual change from one representable step in real frequency.
Final matching and fresh reconstruction must pass `max(1e-8, 2*measured_change)`,
the rule used in the successful spin-0.95 recovery. The coarse tolerance is
not an acceptance criterion. Requested and actual tolerances, both frequency
probes, profiles, energy normalization, fluxes, and provenance are saved.
Python/Julia/shell syntax and the frozen SHA-256 manifest passed before submission.
All five tasks **COMPLETED** in 23–26 seconds and passed their measured-precision
criteria. The alpha=0.3 residual is `1.47906e-8` at tolerance `1.06769e-7`;
alpha=0.2 residuals range from `3.61701e-7` to `1.49161e-6`, at tolerances
`1.40253e-5`–`1.54520e-5`. Fresh reconstruction, positive energies, and all
five profile/solution exports were checked. Together with the original grid and
the separate spin-0.95 recovery, all **12 linear reference pairs are available**.
These retries do not meet the original fixed `1e-8` criterion; the measured
frequency-precision allowance is explicit in each report. Spatial convergence
and tiny growth-rate accuracy remain separate checks. The user will choose
finite QCD populations together with the assistant; no nonlinear QCD run is
submitted.

### October 29: twelve-background linear reference grid attempted

At the user's request, Slurm array **27418189** attempts both linear modes at
all 12 alpha/spin backgrounds: alpha `[0.2,0.3,0.4]` and spin
`[0.8,0.9,0.95,0.99]`. All 12 tasks have finished (23–37 seconds each): all
12 m=1 profiles and six m=2 profiles passed. Six tasks COMPLETED and six FAILED
because their m=2 solve could not reduce the residual below `1e-8`. Failures
are all four alpha=0.2 cases (residuals `1.34e-7`–`1.08e-6`) and alpha=0.3
at spin=0.95/0.9 (`1.17e-8`/`1.48e-8`). All alpha=0.4 modes passed. The
alpha=0.3, spin=0.99 pair exactly reproduces the saved reference frequencies.
No tolerance was relaxed; residual scaling/precision needs a focused check
with the user before accepting the missing secondaries.
Each task calls the existing `solve_linear_mode` for m=1 with
ell `[1,3,5]` and m=2 with ell `[2,4,6]`; no model-library changes were made.
This attempts 24 dimensionless reference profiles for the proposed 36-point
QCD catalogue; 18 accepted profiles are saved. Physical cloud populations and nonlinear QCD runs remain the
next collaborative step with the user.

Frozen source/config: campaign `work/linear_grid_20261029_v1`; outputs:
`results/linear_grid_20261029_v1`; logs: `logs/linear-grid-27418189_TASK.out`.
Each mode retains unit first outer coefficient and records free-field energy
and the factor needed for unit-energy normalization, complex frequency,
profiles/slopes, signed boundary fluxes, and matching residuals. The absolute
matching tolerance is `1e-8`; radial steps are 0.1/0.2 for m=1/2. Domains are
`[-100,675]` at alpha=0.2 and `[-100,300]` otherwise. These initial settings
do not establish spatial/endpoint convergence or accurate tiny growth rates.
The driver attempts the second mode even if the first fails and saves both
outcomes. Shell/Python/Julia syntax checks passed before submission.

### October 29: alpha, spin, and independent cloud-energy sweep submitted

The user requested alpha **0.2, 0.3, 0.4**, spin **0.8, 0.9, 0.95, 0.99**,
and five independent energy pairs `(E1/E0,E2/E0)`:
**(1,1), (0.5,1), (1,0.5), (1.5,1), (1,1.5)**, with
`E0=650.5419114886128 F_a^2 M`. Slurm array **27409247** submits all **60**
configurations with at most four running, one CPU/24 GiB/four hours each.
The user explicitly requested submitting without waiting for the new driver's
small physical validation; validation job **27409219** runs alongside, with
no dependency blocking the array. At submission the array was PENDING and
validation RUNNING. Shell/Python/Julia syntax checks passed; physical validation
and sweep outcomes remain to be checked.

First status check: baseline task 1 COMPLETED; tasks 2–5 RUNNING and have each
passed their first radial solve; tasks 6–60 queued. No production failures.
Validation **27409219 FAILED** after 2m54s: all 153 grid/input checks and the
small initial radial solve passed, then energy measurement raised `KeyError: 2`.
The helper sized its shared azimuthal table from primary harmonics alone;
the small test has only primary m=1 but secondary m=2. The production primary
includes m=5, so its table already includes m=2. The live helper now sizes from
both clouds. Validation **27409399** reruns the check using frozen copy
`work/parameter_energy_20261029_validation_v2`, with an added exact-production-table
comparison. The production snapshot/array remain unchanged. Validation outputs
are `work/parameter_energy_20261029_validation_v2_results` in the campaign.

Next status check (~25 minutes into the array): the same four energy variants
have passed 3–4 radial solves each, with no production failures; 55 cases remain
queued. Validation **27409399 FAILED** after 4m26s. The angular fix passed,
and energy targeting on the coarse mesh reached relative error below `1e-6`,
but the final radial refinement shifted the primary energy by `4.43e-4`, above
the unchanged `1e-4` acceptance threshold. Its rejection is correct. The live
driver now retargets after refinement at differential tolerance `1e-6` and
energy quadrature 0.5. Validation **27410452** tests this from immutable
`work/parameter_energy_20261029_validation_v3`; outputs append `_results`.
The original production array continues on its original snapshot and will
reject any analogous final energy mismatch; preserve those checkpoints for
correction with the new driver if needed. No production tolerances were relaxed.

Validation **27410452 COMPLETED** in 4m47s; its `validation.toml` reports passed.
The refined unequal-energy test reached both targets within the unchanged
criteria after one extra refined correction. At the subsequent ~50-minute
status check, production tasks 2–5 had completed 7–8 accepted radial solves
each, with no failures; only the reused baseline was final. Per-solve times
are mostly 5–7 minutes. The driver splits half-energy changes into three
targets and 1.5x changes into two, restarting the Broyden estimate each time.
It also unnecessarily re-solved identical background/amplitudes at stage
boundaries (about 4 minutes each). The live driver now skips those unchanged
solves; the running immutable snapshot remains unchanged. This small skip
was syntax-checked; validation 27410452 predates it. Last night's single
refinement took 12m00s, while full equal-energy targeting took 41m53s.

Frozen source/config/seed: `runs/adiabatic/multiphase_two_cloud/work/parameter_energy_20261029_v1`.
Outputs: `runs/adiabatic/multiphase_two_cloud/results/parameter_energy_20261029_v1/`.
Logs: campaign `logs/parameters-27409247_TASK.out` and
`logs/parameters-27409219_4294967294.out`. See the campaign README for commands.
The campaign driver varies both energies independently, saves every radial solve
and accepted continuation stage, and reports full-field E/J and boundary fluxes.
It does not modify the model library or perform time evolution. The alpha=0.2
outer endpoint is provisionally enlarged to 675; other endpoints remain 300.
Spatial/angular/endpoint convergence and physical branch reachability across
this wider parameter range are not established.

The October 28 equal-energy job **27401326** actually **COMPLETED** in 41m53s.
Its accepted energies are `[650.5436220051258,650.5431992413581]`, and its
independent defect is `3.383702626475915e-7` at tolerance `1e-6`.
That accepted result is the sweep seed; the original equal-energy grid point
is reused after rechecking its energy, while the other 59 points are attempted.

### October 28: equal-energy joint comparison requested

The user questioned why the joint primary looks closer to the linear reference
than to Omiya's, suggested low cloud occupation, and requested matching cloud
energies instead of outer coefficients. The original joint field has total
finite-domain energy `1109.6799475`, whereas the two reference labels each have
energy `650.5419114886128`. Equal outer coefficients do not fix these energies.

Job **27401326** is running the bounded energy-targeting driver
`runs/adiabatic/multiphase_two_cloud/scripts/equal_energy_target.jl`, using frozen
snapshot `work/equal_energy_20261028_v1` and result prefix
`results/equal_energy_20261028_v1` in that campaign. It adjusts both amplitudes
and re-solves the full joint equations, starting from `src_target_20261028_v2`.
The per-cloud comparison convention is exactly the reference convention:
isolated cosine energy for the primary harmonics, and quadratic secondary
energy on the averaged primary. These labels exclude mixed driven modes;
the full interacting-field total energy is measured separately.
The driver checks the reference energy calculation before solving, records
every amplitude/energy iteration, and saves each accepted radial solution.
Energy targeting uses differential tolerance `1e-5`, followed by `1e-6`
refinement and a radial energy-quadrature check. Limits: eight solves, one CPU,
24 GiB, two hours. Log: `logs/equal-energy-27401326.out`.
Final equal-energy curves have not yet been produced.

### October 28: four reference curves corrected to equal energies

The common primary/secondary branch access is now used by
`runs/adiabatic/multiphase_two_cloud/scripts/reference_profiles.jl`. It reconstructs
both Omiya profiles and the two linear eigenmodes, then exports their radial
samples and provenance. The first export (`v1`, job **27397837**) reproduced the
saved radial arrays exactly, but used fixed outer coefficients for the linears.
The user correctly identified that this normalization was wrong for Fig. 7.

The current [PNG](runs/adiabatic/multiphase_two_cloud/figures/reference_curves_20261028_v2.png),
[PDF](runs/adiabatic/multiphase_two_cloud/figures/reference_curves_20261028_v2.pdf), and
[CSV](runs/adiabatic/multiphase_two_cloud/results/reference_curves_20261028_v2.csv)
use **equal cloud energies**. Job **27398471** completed in **42 seconds**.
Each linear mode is rescaled using its free-field energy with U(phi)=phi^2/2;
all four measured energies agree with **650.5419114886128 F_a^2 M** to roundoff.
This is the saved nonlinear point nearest the paper's approximate 660.
The nonlinear samples are unchanged; linear peaks are now **0.1151975** and
**0.0300686**. The correction passed 63 focused energy/branch/linear checks,
including analytic energy and amplitude-squared scaling tests. The earlier
interface/loader/continuation checks passed 98 assertions.

The campaign README records the frozen source, normalization, settings,
residuals, and reproduction commands. The v1 data remain as a historical
fixed-amplitude comparison. Coupled multiphase target profiles and spatial
refinement remain separate work.

### October 28: collocation backend integration in progress

The user is incorporating the radial BVP experiment into `src` through tested
code blocks that they paste themselves. Shared endpoint assembly, scaled
field/parameter coordinates, the augmented RHS, endpoint residuals, and the
callback factory are now present. The user has also pasted the optional MIRK
backend, its two version-specific library corrections, and the solve wrapper.
**The MIRK loader/corrections/solve-wrapper stage is checked off.** All three
new files match the tested blocks. Seventeen focused numerical checks passed:
the corrected coupled Jacobian agrees exactly with uncompressed AD, exact
mesh transfer is accurate to `1.11e-16`, and a two-fundamental quartic physical
solve converges on 29 nodes with boundary residual `1.83e-15`, zero parameter
drift, and independent off-mesh defect `3.55e-5` at requested tolerance `1e-4`.
Fresh root-project checks also passed for the base loader, optional backend,
repeated loading, correction idempotence, and pinned MIRK version `1.19.2`.

Review found a real include cycle: `load.jl` included `load_collocation.jl`,
which included `load.jl` before the base loader's completion flag was set.
Removed the optional-backend include from `load.jl`. Load the backend explicitly
through `src/adiabatic/multiphase/load_collocation.jl`; it loads the base library.
The root Julia project and resolved manifest now contain the MIRK dependencies
with an exact `1.19.2` compatibility pin. Validation used the experiment's
writable depot followed by the existing user depot.

Tested source blocks and `validation.json` are preserved under
`runs/adiabatic/radial_bvp/work/src_port_20261028__3pqfftq/`.
Initial-guess/scaling preparation and result assessment follow. The small
physical check establishes the wrapper's operation on that finite-domain
problem; a full target solve through the new backend remains outstanding.

The user has pasted `two_phase_collocation_scales` into `collocation.jl`.
The reviewed source and its focused checks are preserved in
`runs/adiabatic/radial_bvp/work/scales_port_20261028_44ut6s/scales_snippet.jl`.
All 21 focused checks passed, reproducing both state and parameter scales
exactly for all four original experiment cases, including per-mode weak-case
floors and the full-reference weights for a zero-driven initial guess.

**Current checkpoint:** the user requested one end-to-end test, then the joint
deformation curves, and explicitly asked to touch base between those stages.
The single four-component quartic test completed successfully: all seven
checks passed, using automatic scales, both fundamentals, and both driven modes.
It finished on 38 mesh nodes with boundary residual `2.3248e-13`, zero parameter
drift, and independent off-mesh defect `4.1154e-5 < 1e-4`. The source, test,
result, and TOML report are preserved under
`runs/adiabatic/radial_bvp/work/assessment_port_20261028_4qlEfV/`.
The staged `assessment_snippet.jl` defines `assess_two_phase_collocation` for
the user to append to `collocation.jl`. This completes the scale/solve/assessment
test on the small finite-domain problem. The user subsequently authorized the
joint-deformation run after this touch-base checkpoint.
The current equal-energy linear/Omiya curves are exported as
`results/reference_curves_20261028_v2` in the multiphase campaign; the earlier
fixed-amplitude `v1` export is preserved separately.

**Target run completed:** job **27398811 COMPLETED** in 12m00s, using frozen snapshot
`work/src_target_20261028_v2` in the multiphase campaign and output prefix
`results/src_target_20261028_v2`. It restarts the refined 14-mode target with
the new optional backend at differential tolerance `1e-6`, initial spacing
0.5, nonlinear tolerance `1e-8`, and AD chunks of 28. Limits: one CPU, 24 GiB,
55 minutes, and a 2400-second callback solve limit. The source and all inputs
are hash-checked. Its driver packages the tested assessment helper until the
user pastes that helper into the model library. Log: `logs/src-target-27398811.out`.
The 881-point result passes assessment: boundary residual `1.68835e-13`,
parameter drift `9.69e-21`, and independent defect `3.79897e-7 < 1e-6`.
Its maximum frequency change from the older refined checkpoint is
`3.21794e-13`; maximum per-mode global relative profile change, including
slopes, is `5.49e-7`. Timed solve work was 664.7 s including 127.7 s of
compilation; peak batch memory was about 8.7 GiB. This is a current-source
reproduction of the existing finite-domain target.
Earlier job 27398712 failed before solving because the old refined checkpoint
does not contain a `passed_residual_checks` summary field; v2 verifies its
saved success flag, boundary residual, defect, and parameter drift directly.

**Comparison normalization clarified by the user:** match the energies of the
linear and Omiya reference curves. Use the `reference_curves_20261028_v2`
metadata: reference E is `650.5419114886128`, and the linear amplitude factors
relative to the old fixed-amplitude references are `0.7255014992542931` for
the primary and `0.7541792626953592` for the secondary. The joint solve retains
the Omiya outer coefficients; its cloud energies have not been independently
retuned. The target driver exports the preserved fixed-amplitude baseline;
`plot_collocation_comparison.py --energy-references ...v2.toml` applies the
verified linear factors to both R and Q and labels the joint normalization.
Keep this distinction in the plot caption and exported metadata.

**Comparison artifacts generated and inspected:**
`plots/multiphase_deformation_20261028_v1.{png,pdf}` contains the six-profile
overlay and signed changes relative to each energy-normalized linear peak,
including the extra joint-minus-Omiya curve. The same prefix's `-curves.csv`
contains all six complex profiles and slopes after the linear rescaling;
`-plot-data.csv` contains the plotted magnitudes and fractional differences;
`-plot.json` records normalizations, linear factors, run diagnostics, and
input/source hashes. The original target CSV remains at its original outer
normalizations. The energy-reference conversion was checked against all
12,008 original reference rows, agreeing within `2.78e-17` absolute.
Joint cloud energies remain untuned; angular, boundary, and retained-support
systematics still require separate checks. The new curves are preliminary.
The refined target checkpoint `radial_bvp/results/target_driven_mirk4_refined_v1.jls`
is available as a warm initial guess for a current-source target reproduction.

**Required readability follow-up:** overhaul
[`src/adiabatic/multiphase/collocation_mirk.jl`](src/adiabatic/multiphase/collocation_mirk.jl).
The user finds the current wrapper difficult to understand and accepts it for
this stage. Revisit the validation, coordinate setup, initial-guess callback,
MIRK configuration, solve, and result decoding with them; make the flow easier
to follow while preserving the tested numerical behavior. Do not treat this
stage's numerical check-off as completion of that readability work.

**Learning-map update completed October 28:** updated
[`docs/multiphase-learning-map.md`](docs/multiphase-learning-map.md) after this
integration work, as explicitly requested by the user. It now includes the
collocation helpers, optional loader and corrections, scales, assessment,
function relationships, and the outstanding MIRK readability work. Every
existing understanding-checkbox state was preserved, and new entries start
unchecked. Keep the map aligned when the assessment helper is pasted into
the library or the wrapper is simplified.

### October 28: investigate a direct solve at the plot amplitudes

**Twelve additional driven-seed corrections authorized and submitted:** the user
requested twelve more corrections from the last driven-field state. Job
**27388374** resumes correction 6, residual `1692.4986686897735`, and allows
corrections **7–18**, stopping earlier on convergence or a solver failure.
It uses the same verified frozen library, full target amplitudes, join 60,
boundary prescription, residual weights, full-width local AD, and `1e-5`
tolerance. New standalone restart scripts preserve the prior run and use output
prefix `results/target_match60_driven_20261028_v2`; log:
`logs/target-match-resume-27388374.out`. Input checkpoint SHA-256:
`f680a44bb983ede293e36575f27300ffe4472d4a063ba2061611f0278fa4476f`.
Setup, syntax, source hashes, checkpoint identity, and total correction-count
checks passed. The first residual evaluation also checks fresh agreement with
the stored mismatch. Every accepted checkpoint retains the cumulative count.
Resources: one CPU, 8 GB, one-hour limit. Check this job before another extension.

**Driven-field seeds added:** at the user's request, job **27384564 COMPLETED**
in 1m33s using standalone `work/driven_field_seed_20261028_v1.jl`, the verified
prior library snapshot, and frozen reference inputs. All forced-field matching
and power checks passed; the five joint initial guesses are saved under
`results/driven_field_seed_20261028_v1`. The full-target residual is 56.5614
with zero driven fields, 3529.9916 with m=0 only, 57.0283 with m=3 only, and
3529.2624 with both fields using the joint boundary prescription. The increase
is mostly in the primary mismatch. Using far-leg m=3 reference coefficients
instead gives 3529.0943, nearly the same overall result. The near/far m=3
outer states differ by 28.3% for identical leading coefficients; this is not
a boundary-subspace comparison after refitting. See the investigation record.

Job **27384946** tested six Newton corrections with both driven fields, at the
same target amplitudes, join 60, boundary prescription, and `1e-5` tolerance.
It **FAILED at the six-correction cap after 8m36s**, ending 13:07:02 PDT:
residual `3529.2624 → 1692.4987`, about a 52% decrease. Trial evaluations per
correction were 4,3,4,3,3,1; the last correction reduced the norm by 10.7% and
accepted its first evaluated trial. This progresses more than the zero-driven
test, but remains far above tolerance and does not establish fast convergence.
Input: `results/driven_field_seed_20261028_v1-both_near.jls`; output prefix:
`results/target_match60_driven_20261028_v1`; log:
`logs/target-match-newton-27384946.out`. The last accepted state and failure
report are saved; no matched target profile exists and no further target job
was submitted at this review. A further bounded attempt can resume this state.
The production library and reference driver were left unchanged; source files
in the shared workspace have concurrent edits. At the same review continuation
**27378594** was RUNNING after 2h01m, with 13 new matched points. Latest
A2=`5.394443404258194e-6` (4.95082% of target), norm `6.00678166743477e-6`;
attempt 14 is underway with the increment now 0.78125%.

The user requested reevaluating the workflow and explicitly authorized outlining
and beginning a fixed-target shooting investigation. See the
[investigation plan](docs/multiphase-target-solve-investigation-2026-10-28.md).
Job **27380795** COMPLETED in 5m19s, comparing joins at r*=5, 10, 20, 40, 60
using the original target-amplitude guess and verified frozen source. The old
baseline reproduces exactly. At r*=60 the horizon/outer field norms are
0.06573/0.06340, the relative field/slope jumps are 0.09335/0.29657, and the
scaled residual is 56.5614 versus 1.90648e9 at r*=5. Full Newton steps are
admissible at 40 and 60. Native Jacobian conditioning remains severe; at 60
the native/equilibrated condition estimates are 1.95e16/2.99e3. The plan records
all five comparisons and their limits.

Selected r*=60 for a bounded **six-correction direct solve at full target A2**,
job **27381429**, with unchanged `1e-5` tolerance and the same frozen library.
It **FAILED at the correction cap after 13m32s**, ending at 12:27:05 PDT:
residual `56.5614 → 53.8721`, only a 4.75% decrease. Six trial integrations
were needed for the first correction and nine for each later correction. The
last accepted state and failure report are preserved under
`results/target_match60_20261028_v1`; no matched target profile exists. Log:
`logs/target-match-newton-27381429.out`. The changed join improved the initial
mismatch but did not yield useful Newton convergence for this guess.

**Saved-matrix follow-up:** explicit equilibration and a 256-bit solve of the
same saved initial Jacobian agree with the original Newton correction to
`9.46e-15` relative norm. Inaccurate linear-system solution does not explain the
first failed full step. That step predicts cancellation in the linear model but
fresh integration gives residual 6118.8 versus 56.56 initially; the first accepted
step is only 1/32 of it. The largest boundary-coefficient changes are in the
mixed driven components initialized to zero. See the investigation record and
`logs/target-match60-linear-check-20261028-v1.out`. Next candidates are testing
driven-field seeds or isolating those correction blocks, and checking the local
derivative along the actual target Newton direction. These remain hypotheses
and proposed tests. Merely raising the correction cap is not recommended.
Existing continuation **27378594** remains RUNNING; at the same review it had
saved seven new points, reaching A2=`5.069450807458185e-6` (4.65255% of target)
with weighted norm `5.828310199789505e-6`; attempt 8 was underway. These
diagnostics change no production library or driver.

### October 28: amplitude attempt 25 recovered

The goal remains jointly solved primary/secondary profiles at the plot's fixed
outer amplitudes: A1=`0.036390480353373555`, A2=`1.0896058341192524e-4`.
Continuation job **27337026** reached its four-hour limit on October 27 at
21:54:43 PDT after saving 24 matched points. Attempt 25's accepted checkpoint
had already reached tolerance at correction 8, but its radial profiles were
not saved before timeout. The queue was empty at the October 28 review.

The user explicitly requested recovering attempt 25. Recovery job **27377949**
COMPLETED in 33 s on n0003.lr7 using the verified original 201-file snapshot.
Fresh integration reproduces the checkpoint's weighted residual vector exactly;
its norm is `2.153421086390111e-6 < 1e-5`. Both radial branches are finite,
span the original endpoints, and meet at r*=5. The complete matched seed is now
`runs/adiabatic/multiphase_two_cloud/results/amplitude_step_20261027_v1-attempt25.jls`.
It preserves A1, A2=`4.54807098617202e-6` (4.17405% of target), the saved
Newton state, residual weights, and correction count. The original accepted
checkpoint is unchanged. See the campaign README for recovery provenance.

**Continuation restarted on October 28:** the user then explicitly requested
further continuation. Changed only the live driver's initial fractional step
to `0.015625` and output prefix to `amplitude_step_20261028_v1`. Snapshot
`work/amplitude_continuation_20261028_v1` verifies all 201 files; core/multiphase
source matches the previous snapshot and its frozen input is recovered attempt
25. The frozen driver setup and saved-seed checks passed. Job **27378594** is
RUNNING on **n0064.lr7**, with one CPU, 8 GB, and a four-hour limit. Its first
new amplitude is `4.6191345953309576e-6`; the initial weighted residual is
`606644.1247245938`. Check `logs/amplitude-27378594.out` and Slurm before
submitting another continuation. Results use the new prefix with attempt
numbering restarted at 1. No target profile is established by this launch.

The user raised concern about time-evolution cost. The current evolution code
corrects midpoint and endpoint profiles and evaluates charge responses at both;
it already uses an implicit amplitude-response predictor, whereas this
continuation uses only coefficient scaling with unchanged predicted frequencies.
Thus the current correction cost is a real performance warning, but is not a
benchmark of a complete evolution step. At a constant 1.5625% increment about
205 continuation steps remain. A focused next experiment is comparing a secant
or implicit-response predictor at the same amplitude increment and tolerance,
then timing a full evolution sample/step including response calculations before
committing to a long trajectory. This experiment has not been submitted.
The target-amplitude joint profiles, spatial refinement, and comparison plot
remain outstanding. This checkpoint supersedes the October 27 launch-only
status below; preserve the user's small-snippet workflow for further changes.

### October 27 run review and workplan

The [October 27 workplan](docs/multiphase-workplan-2026-10-27.md) records today's
priorities and verified outcomes. The user's end goal for today, potentially
a stretch, is the target-amplitude linear/Omiya/joint deformation plot for
both clouds, including difference panels and PDF/PNG/CSV exports. The seed,
derivative checks, continuation, and reference curves are milestones toward
that figure; prepare the plotting layout while solves run.
At approximately 13:08 PDT, Slurm showed an
empty queue for `clwelch`; accounting and campaign logs confirm all five
October 26 jobs ended that afternoon, with no additional overnight run.
At that morning review, the weak seed remained unmatched at
`2.274196817352566e-5` versus `1e-5`, and its final accepted state was not saved.
The October 26 independent diagnostic preserves
the initial guess and one directional comparison, not that final iterate.
Later October 27, the user authorized completing derivative verification
autonomously. The [verification record](docs/multiphase-ad-verification-2026-10-27.md)
tracks the fixed input, criteria, frozen snapshot, and results. All 55 focused
local checks passed. Job **27333738** COMPLETED: all 56 columns at the initial
weak guess and both finite-difference directions pass. Job **27334006** COMPLETED after
22m45s: correction 13 and fresh reintegration give `5.306800040310494e-6`, below
the unchanged `1e-5` tolerance. The seed and both radial branches are saved in
`results/ad_seed_20261027_v1-seed.jls` under the multiphase campaign. It retains
A2 at 0.1% of its target. Job **27334124** COMPLETED at 15:59:48 PDT:
all 56 matched-state columns pass (worst relative error `9.587805e-9` versus
`1e-8`) and both finite-difference directions pass. The agreed derivative
verification is complete for the tested configuration/states; target-amplitude
convergence and spatial accuracy remain separate checks. All three ran on
**n0060.lr7** using verified snapshot `ad_verification_20261027_v1`; their
artifacts survived the n0003 SSH outage on shared storage. Saved comparisons,
seed state, radial branches, and the matched-stage seed hash were checked
during recovery. The queue was empty at that recovery checkpoint.

**Amplitude continuation launched later October 27:** after building
`continue_amplitude.jl` together in pasted, commented blocks, the user authorized
creating the Slurm wrapper and submitting it. Job **27337026** started at
**17:54:21 PDT** on **n0060.lr7**, with one CPU, 8 GB, and a four-hour limit
(21:54:21 PDT). Snapshot `amplitude_continuation_20261027_v1` verifies all 201
files; its core/multiphase library is byte-for-byte identical to the verified
AD snapshot, and its frozen input is the matched weak seed. The initial A2
step is `1.1985664175311778e-7`; its weighted residual fell from
`2354.7637801303267` to `0.009104341308015403` after correction 1. At this launch
check, its accepted-state file exists but no new matched point is saved yet.
The target is `1.0896058341192524e-4`, with A1 and spatial settings fixed.
See the campaign README for the run command, output names, and adaptive rules.
The plot remains the end goal. Check the active job before submitting anything else.

**Workflow reminder:** the user rejected an autonomous driver restructuring
even though it stayed below their roughly 160-line cap. They now explicitly
prefer code to paste themselves. Keep the existing script structure, names,
and comments; show small tested edits. The exception for creating and running
the submission shell does not authorize another driver rewrite.

### Important: revisit multiphase architecture

**Explicit user priority, October 26: revisit this architecture and prevent
bloat.** Acceptance of the current ForwardDiff specialization lets us proceed
toward another convergence test; it does not close the architecture discussion.
Keep that convergence attempt the immediate priority, then return to this review.

Review the responsibilities of source evaluation, the local Jacobian, and the
ForwardDiff adapter together. Look for duplicated source mathematics, repeated
packing/unpacking, unnecessary wrappers, and confusing function/method names.
Keep the physical calculation easy to follow and justify any additional
abstraction with a concrete need. Preserve numerical agreement and the measured
speedup while simplifying. Continue the user's tested-snippet workflow for
library changes unless they explicitly delegate implementation.

### Resume here: direct mutual-deformation profiles and source optimization (October 26)

The requested [October 26 workplan](docs/multiphase-workplan-2026-10-26.md) organizes:
projection walkthrough, a small reviewed change for another convergence test,
and the fixed-amplitude linear/Omiya/joint deformation comparison. The source
walkthrough and reviewed batching edit are complete; verified v5 was submitted
as job **27304701**, which exhausted all 12 corrections at scaled residual
`1.4947468124727777e8`. No joint solution or profiles were saved. See the
checkpoint below.

The current deliverable is a Fig. 7-style primary/secondary profile overlay
with two additional jointly nonlinear curves. A direct solve fixes both outer
coefficients selected from the equal-energy Omiya reference; joint energies
are outputs. The commented driver, frozen inputs, logging, and run settings
are recorded in the [October 23 work summary](docs/multiphase-work-2026-10-23.md)
and [campaign record](runs/adiabatic/multiphase_two_cloud/README.md#direct-fig-7-profile-attempt-october-23).

Checked on October 26: jobs **27241421** (v2) and **27241725** (v3, detailed
logging) both ended with **TIMEOUT** at their four-hour limits on October 23,
at 19:13:04 and 19:29:40 PDT respectively. Both accepted two corrections,
with scaled residuals 1.90648e9 -> 6.83206e8 -> 4.65316e8, still far above
the 1e-5 tolerance. The v3 log stops during correction 3's third Jacobian
pass; its first two Jacobians took about 70.6 and 71.3 minutes. Only the
initial-input files exist: accepted Newton states were not checkpointed,
and no target joint solution or profile CSV was produced. Preserve these
attempts. Accepted-state checkpointing remains a possible follow-up.

A standalone cached-basis/local-source-Jacobian prototype measured 55.8x
faster Jacobian construction through short endpoint integrations. Full-domain
derivative agreement and Newton convergence remain unresolved; see the v4
timing and failed check below. Guided implementation
has added the basis builder, ordinary matrix source, and `local_source_jacobian`
to `projection.jl`. The builder returns an immutable `TwoPhaseSourceBasis`
containing `mode_labels`, `basis`, `weights`, and `cosine_theta_squared`.
At the user's request, the source now requires this precomputed basis data:
only `nonlinear_two_phase_source(parameters, source_basis)` remains. The
dictionary method and uncached fallback were removed. `TwoPhaseProblem`
constructs the basis once by default, and radial/matching/response calculations
reuse it. The radial RHS packs its coefficient vector and checks that its mode
order agrees with the basis. Direct quadrature remains a test reference.

The multiphase suite passed 388 checks through source, radial, matching,
observable, and profile-coordinate tests. It then stopped on two existing
evolution-test errors: `test_evolution.jl` expects `amplitude_rates`, while
`two_cloud_balance_rates` returns `log_rates`. The same return fields were
verified from the committed function; neither evolution file was changed.
After renaming the type, builder, and problem field to source-basis names,
147 focused source/RHS/matching/profile-coordinate checks passed again,
including repeated loading.
The user has now added the derivative-specific source method and Newton's
chunk-size/accepted-state callback options. The assistant's driver additions
passed 112 focused Newton/source/driver checks and were frozen in v4. User-submitted
job **27296945** matched the old initial residual in 14.63 s and built the
full-domain Jacobian in 92.51 s. Its first directional finite-difference check
reached only 2.38e-4 relative agreement against a 1e-5 threshold, so the preflight
stopped the job before Newton. The discrepancy is unresolved; the timing does
not establish correctness or convergence. Preserve the v4 snapshot and log.

**User-directed reset, October 26:** the assistant added too much driver code
without walking through it. Restored `direct_profile.jl` and
`submit_direct_profile.sh` to their prior versions and deleted the live
`preflight_direct_profile.jl`. The library changes remained; at that reset the driver used
default ForwardDiff batching and did not save accepted iterates. Resume the
request to expedite another convergence test, one reviewed block at a time,
including driver changes. The smallest proposed driver change is to pass the
existing `jacobian_chunk_size` option; the reviewed edit is recorded below. Discuss
the derivative discrepancy before choosing the next run. Preserve the user's
tested-snippet workflow and learning-map checkboxes.

**Reviewed v5 preparation, October 26:** the user added only
`jacobian_chunk_size=length(initial)` to `logged_profile_newton`'s existing
Newton call. The actual driver passed five known-answer checks, including
one derivative evaluation carrying all 56 directions. The setup-only check
passed with unchanged physical settings; 38 focused source checks also passed.
Snapshot `work/direct_fig7_20261026_v5` verifies all 194 files, and its three
reference inputs have the same hashes as v4. The user submitted job **27304701**;
Slurm confirms the frozen v5 wrapper, snapshot, and output prefix. It was
PENDING at the first check and subsequently RUNNING. Its initial scaled
residual is `1.906482699496033e9` in 12.44 s, identical to v4. At the subsequent
18-minute Slurm check it remained RUNNING; the inspected log had accepted
correction 7 at `2.805706872654804e8` and begun correction 8. The first seven
Jacobians took 51–70 s each. The residual has fallen by a factor of 6.80 but
remained `2.81e13` times the `1e-5` target, with slow recent progress.
Final status: FAILED (exit 1) after 29m04s, with all 12 corrections accepted
and final residual `1.4947468124727777e8` at matching-timer 1721.05 s.
All 12 Jacobians took 51.41–69.71 s. The input and failure report are preserved;
there is no solution, profile CSV, or saved accepted iterate to resume.
Review the failure and unresolved derivative comparison before choosing
another Newton attempt; no automatic resubmission or continuation.
The scale is the fixed secondary outer coefficient `1.0896058341192524e-4`,
not a local field norm. A fresh integration of the saved initial guess with
the frozen v5 source gives matching-point field norms 0.05055 (horizon side)
and 129715.5 (outer side), dominated by the secondary l=6 mismatch. This is
a severe initial joining failure; intermediate states are not saved, so the
same component breakdown is unavailable for the current iterate.
Follow-up controls with frozen v5 source reproduce the original primary's
unscaled mismatch `2.0333e-12` and secondary's unit-amplitude mismatch
`3.1578e-6` (`3.4408e-10` at the target A2). In the full joint system,
zeroing the secondary boundary coefficients gives `2.1184e-12`; multiplying
them by 0.1 gives `900.405`, versus `207731.467` at the target amplitude.
Halving radial spacing changes the latter by only `1.49e-6` relative.
This points to the direct finite-secondary guess and amplification during
inward shooting, rather than the reference solutions or coarse radial spacing,
as the immediate failure. It does not yet identify which added interaction
dominates or prove a remedy. See the campaign diagnostic record.
The full-domain derivative discrepancy remains unresolved; there
is no preflight gate or accepted-iterate saving. See the campaign README for
the exact submission command and run record.

**Reference-driver workflow:** the user explicitly requested building the
reference comparison together, including driver code, one small tested snippet
at a time. They have inserted the input-loading and Omiya-selection blocks in
`scripts/reference_profiles.jl`. The primary-reconstruction block is now also
inserted, together with the user's explanatory `joint_input` field tree.
It passed Slurm check 27305885: scaled primary mismatch `5.587385657966086e-11`,
with change from the saved residual `8.671068890112971e-29`.
The user resumed this work and explicitly retained manual driver edits.
Job **27337260** tested the original secondary snippet: matching passed, but
separate integration at A2 disagreed with scaled unit samples by `8.24126e-10`
on the inward branch, above the proposed `1e-10` scaling threshold. The revised
snippet reconstructs the unit solution and scales its saved samples afterward,
valid for this linear equation on a fixed primary. Job **27337856** passed:
normalized mismatch `3.157814155718417e-6`, sample-scaling errors below `2e-16`,
and independent RHS-scaling checks at 14 radial samples below `1.6e-16`.
Use `work/check_reference_secondary_20261027_v2.jl` for the tested block.
Job **27337858** also passed both linear eigenmode/profile checks; candidate
blocks are in `work/check_reference_linears_20261027_v1.jl`. Scaled matching
norms are `3.73863e-13` and `6.78907e-9` for primary and secondary. The four
reference profiles now exist in `results/reference_secondary_20261027_v2.jls`
(both Omiya profiles) and `results/reference_linears_20261027_v1.jls`.
The live driver still ends at primary reconstruction; give the user the tested
secondary block first, then walk through the linear blocks. Combined CSV/plot
export and spatial refinement remain unfinished. Preserve the tested-snippet
workflow and the learning-map boxes.

**User-authorized weak-secondary continuation, October 26:** after the direct
failure and control calculations, the user agreed to obtain a weak-secondary
seed, verify its matching Jacobian, and then increase the secondary amplitude
in small steps. **Stopped at the user's departure request; no further job or
amplitude ramp is queued by this work.** Job **27305905** ran from frozen snapshot
`work/weak_secondary_20261026_v1`. The first stage fixes the primary and all
spatial settings, using A2=`1.0896058341192524e-7` (0.1% of the target).
Its residual uses the old scale for primary-only rows and the smaller A2 for
all secondary-dependent rows, with norm tolerance `1e-5`. It saves a matched
seed before comparing the optimized Jacobian against ordinary-source AD and
centered finite differences in two directions. Source AD checks and all 38
direct-quadrature checks passed before submission. No shared library was edited.
Job 27305905 FAILED at the 12-correction cap after 14m00s: weighted residuals
`1.3255351e7 -> ... -> 124.01057 -> 0.001374344 -> 2.274196817352566e-5`.
The last value is only 2.274 times the `1e-5` tolerance, but it did not pass.
No seed or accepted iterate was saved, because this driver saves only after
convergence. Its post-convergence Jacobian check was never reached.
During this work the user explicitly requested more frequent, understandable
updates and decision points. Explain each numerical experiment before starting
it; report progress roughly every minute with its implication. The current
attempt only covered the weak seed and its derivative checks. Bring that outcome
back to the user before starting an amplitude ramp; no ramp is submitted.
After slow seed progress, the user explicitly selected an independent Jacobian
check at the original weak-secondary guess while the seed attempt finishes.
Submitted job **27306004**, snapshot `work/weak_jacobian_20261026_v2`, using the
same diagnostic helper without Newton; log `logs/weak-jacobian-27306004.out`.
This tests derivatives at the initial guess, not at a converged solution.
Job 27306004 FAILED after 2m31s at the first finite-difference check, preserving
`results/weak_jacobian_20261026_v2.jls` (initial state, full optimized Jacobian,
first direction and errors). Ordinary-source and optimized AD agreed to
`7.47785e-15` relative in that direction. Finite-difference errors were
0.5103, 0.3044, and 0.1636 at steps `1e-6`, `1e-7`, `1e-8`; direction 2 was
not reached. The predicted derivative norm is `5.8345761e17`; even step
`1e-8` predicts a change 440 times the initial residual. These probes are not
demonstrably local, so their failure does not isolate a derivative bug.
Slurm confirms the direct run and both diagnostics have finished.

Next decision on return: review a small driver change to save accepted states
using the existing callback before retrying the nearly matched weak seed with
a modestly larger iteration cap. Keep the `1e-5` tolerance. Use the saved
Jacobian to investigate smaller, properly scaled finite-difference probes;
do not automatically relax the derivative threshold or start continuation.
No retry or further diagnostic was submitted after the departure request.

### Physical-profile organization checkpoint (October 21)

The physical-profile and solver-coordinate migration is complete for this pass.
The user authorized the solver migration and the remaining annotation, test,
driver, naming, and documentation cleanup. The broader docstring overhaul is
unfinished. Ordinarily the user prefers to make library edits themselves with
small guided steps; this task-specific authorization does not replace that
preference. Preserve the user's learning-map checkboxes.

Current interfaces and responsibilities:

- `TwoPhaseProblem` stores fixed settings; `TwoPhaseProfile` stores complex
  frequencies and physical boundary coefficients. Neither stores sampled
  radial functions. Both are defined in `profiles.jl`.
- `pack_two_phase_profile` and `unpack_two_phase_profile` convert between a
  physical profile and 4N real coordinates (2N complex unknowns). Frequencies
  are unscaled; boundary coefficients use the fixed problem scale. The two
  prescribed amplitudes are separate inputs to unpacking. Numeric promotion
  preserves derivatives from either the packed state or amplitudes.
- `two_phase_matching_residual(profile, problem; save_radial_result=false)`
  evaluates physical matching, including reduced single-cloud component sets.
  Saving returns `from_horizon`, `from_outer`, and `residual`. The integration
  calls translate the flag to the core keyword `save_profile`.
- `solve_two_phase_mode` constructs a problem and packs its complex Newton
  candidates before unpacking profiles. Integration defaults, the Newton
  backend, scaling, admissibility behavior, and return fields are preserved.
- Response residuals and observables use the same unpacker and matching
  interface. `two_phase_trial_matching` and `unpack_two_phase_trial` are removed.
  Library and live-driver definitions accepting a problem now annotate it as
  `TwoPhaseProblem`; physical profile arguments use `TwoPhaseProfile`.
- `two_phase_integrated_observables` calls its integration bundle
  `radial_matching_result`; this bundle is not a `TwoPhaseProfile`.
- The unpacker, standalone solver, and live response/evolution driver use
  `two_phase_coefficient_indices`. Preserve the descriptive field names
  `fixed_outer_coefficient_indices` and `free_outer_coefficient_indices`.
  The driver's admissibility check now receives its already-packed guess
  directly.
- `radial.jl` contains `component_tortoise_operator`, `component_boundary_state`,
  and `nonlinear_two_phase_rhs!`. Keep these helpers independent of the profile
  types. The opt-in loader includes this file before `profiles.jl`.

Validation against the edited working tree:

- All **97 focused checks** passed in a fresh Julia process: 24 physical
  matching checks, 49 integrated-observables checks, and 24 coordinate checks.
  They cover single-cloud limits for three potentials, ordering with unequal
  angular blocks, radial RK4 refinement, saved versus plain residuals,
  manufactured charges, boundary signs, packing/scaling, exact Jacobians with
  respect to state and amplitudes, and invalid inputs. Command:

  ```bash
  julia --startup-file=no --project=. -e 'using Test; include("src/adiabatic/multiphase/load.jl"); include("tests/adiabatic/multiphase/test_matching_residual.jl"); include("tests/adiabatic/multiphase/test_integrated_observables.jl"); include("tests/adiabatic/multiphase/test_response_coordinates.jl")'
  ```

- The preceding solver migration passed all **20 nonlinear matching checks**,
  including a reordered-component solve. Newton took three iterations to a
  scaled residual norm of approximately 1.22e-10 on the test's modest grid.
- The live response/evolution driver loaded without running its campaign;
  two setup checks confirmed the prescribed and free coefficient indices.
- These are algorithm/interface checks, not physical spatial convergence.
  The full suite and campaign were not rerun; saved campaign results below
  remain historical evidence for their frozen source.

Deferred work and known caveats:

- Continue the user's docstring walkthrough using the
  [learning map](docs/multiphase-learning-map.md). The
  [organization proposal](docs/multiphase-profile-organization.md) includes
  later work beyond this completed interface migration.
- Boundary-assembly cleanup and making the corrector return
  `TwoPhaseMatchedSolution` remain deferred. The corrector still returns a
  named tuple with packed state, amplitudes, and diagnostics.
- Resolve admissibility semantics before sharing predicates: the standalone
  solver rejects nonpositive component frequencies and exact thresholds;
  the response predicate checks fundamental bounds and branch classification.
- `two_cloud_balance_rates` returns `log_rates`, `scaled_condition`, and
  `balance_residual`, but `test_evolution.jl` still expects `amplitude_rates`.
  That pre-existing mismatch is outside this pass and remains unresolved.
- Preserve existing working-tree changes, frozen source snapshots, and saved
  results. Saved objects require explicit compatibility review before reuse
  with new interfaces; do not resume a historical campaign under edited source.

### Current long-evolution checkpoint (October 19)

The "Review adiabatic workplan" chat ended while preparing a primary-amplitude
doubling run, after the user narrowed the scope to one trajectory and a linear
growth comparison. No job had been submitted. This has now been continued as
Slurm job **27100648**, submitted to `lr7`, with frozen snapshot
`multiphase_two_cloud/work/doubling_20261019_v2` and result prefix
`cosine_doubling_20261019_v2`. It completed successfully in 17m52s, reaching
the target at `t/M=4607986.0447`. The secondary grew 0.029340%. The primary
endpoint is 0.016464% above the linear exponential and 0.010753% below the
exponential using its initial E/J rate. Relative energy/angular-momentum
budget errors are 1.958e-7/9.954e-8 of their initial charges.

The user subsequently requested joint/linear/Omiya comparison curves with
matching initial amplitudes. The saved Figure 6 coefficient table does not
cover this dilute amplitude range. Job **27103035** computes a small branch
using the existing nonlinear-primary/linear-secondary solvers, on the joint
pilot's domain and retained support, with primary/secondary radial steps
0.1/0.2. It completed in 1m01s using snapshot `omiya_comparison_20261019_v2`;
the driver passes 18 analytic normalization/time checks. Time-step halving
and branch-table thinning change amplitudes by at most 4.66e-15/4.11e-14.
The joint endpoint is 0.0167171% higher in A1 and 0.0000539919% lower in A2
than the Omiya reference. The comparison plot now contains joint, linear,
and Omiya curves, with actual ratio axes near one (no offsets/subtractions).
The matched Omiya primary energy remains below 0.001 Fa^2 M; this pilot has
not reached the paper's strongly nonlinear cloud energies. See the
[comparison record](runs/adiabatic/multiphase_two_cloud/README.md#omiya-approximation-comparison-october-19).

The driver reuses the October 15 pilot's initial joint profile and stops when
the measured primary amplitude reaches `2e-4` from `1e-4`; the secondary starts
at `1e-5`. The initial-rate doubling estimate is `t/M=4.607271e6`, with 64
nominal steps and corrected, checkpointed endpoints. It runs one trajectory,
then plots both clouds against independently solved linear eigenmode growth
on the same grid and the initial E/J growth rates held constant. The preceding attempt passed all 47 balance/midpoint/target checks, then
failed before evolution because the linear-reference solver received unsupported
boundary-order keywords. The corrected driver passes only supported radial
settings and checks that the boundary orders agree. Both attempts are preserved. See the
[campaign record](runs/adiabatic/multiphase_two_cloud/README.md#october-19-primary-amplitude-doubling)
for the exact recipe, outputs, and restart constraints.

The immediate plan is to measure the departure from exponential growth, then
test selected profiles for spatial/domain/support sensitivity at that effect's
scale. A longer run does not by itself establish appreciable mutual cloud
deformation. The next physics milestone remains a converged effect of a populated
secondary on the primary, with a justified slow-growth formulation. Actions
for three clouds and BH/GW/self-gravity extensions remain subsequent work.

### Current multiphase checkpoint and resume point (October 15)

The joint profile solver, finite-domain E/J observables, two-amplitude implicit
responses, and a short-evolution integrator are implemented. On October 15 the
user authorized direct implementation of all three response/evolution steps,
with readable comments and docstring skeletons left for them to complete.
The matched-cloud derivative validation and short cosine evolution are complete;
results, files, and commands are in the
[joint two-cloud campaign](runs/adiabatic/multiphase_two_cloud/README.md).

- `multiphase/load.jl` now includes observables, `response.jl`, and `evolution.jl`.
  It remains opt-in; the other models do not import the joint model.
- The two positive amplitudes fix the first retained outer coefficient in
  q=(1,0) and q=(0,1). E and J remain outputs, not amplitude coordinates.
  `unpack_two_phase_trial` is shared by the original solver and response code.
- `response.jl` supplies a real-packed matching residual, a Newton corrector,
  and total E/J derivatives. The coefficient/residual scale is fixed at the
  reference amplitude rather than differentiated. `core/sensitivity.jl`
  factors the state Jacobian once for both amplitude directions and applies
  the full observable chain rule. Existing scalar-amplitude APIs remain intact.
- The observable adapter explicitly takes real parts of the complex matching
  frequencies as carriers, omitting slow profile time derivatives. This
  specifies the numerical family; it does not validate the physical
  slow-growth approximation. Charges are integrated over the supplied finite
  radial domain without tail corrections. Horizon losses negate inner currents.
- `evolution.jl` solves the 2x2 E/J balance system, reports a dimensionless
  condition number after row/column scaling, and rejects ill-conditioned
  matrices. Explicit midpoint steps evolve log amplitudes. Every stage is
  independently Newton-corrected from an implicit-response predictor.
  Independently accumulated boundary losses supply numerical charge budgets.
- All 985 Julia unit checks passed in Slurm: 13 loader, 350 core, 243
  single-cloud, 16 primary/secondary, and 363 multiphase. These include 14
  vector-response, 13 unpacking, 8 balance/scaling, and 19 midpoint/budget
  checks. The campaign then validates both E/J response columns against
  independently corrected neighbors at 1% and 0.5% offsets, checks fixed-scale
  invariance, and compares 4-step and 8-step cosine trajectories.
- Completed pilot job 26978242 uses snapshot `response_20261015_v2`, matching
  tolerance 1e-9, and prefix `cosine_pilot_20261015_v2`. Both derivative columns
  agree with independently corrected neighbors to at worst 9.50e-9 relative
  in E/J and 2.07e-9 in the state response. Changing the fixed coefficient scale
  by seven preserves the E/J derivatives.
- Four-step and eight-step cosine trajectories reach t/M=199406.6212.
  A1 grows 3.04542%, A2 grows 0.00126903%. Endpoint amplitudes differ by at most
  7.06e-11 relative when halving dt. Eight-step budget errors are 4.84e-10 of
  initial E and 2.46e-10 of initial J (2.08e-6 and 1.88e-6 of their respective
  accumulated changes), improving by a factor of four from the coarse run.
  The scaled balance condition stays between 1140 and 1211; saved endpoint
  matching residuals are below 1.97e-10. These are finite-grid numerical checks.
- The preceding job 26977745 passed all 985 unit checks, then stopped in a
  linear seed solve at residual 1.26e-10 above its requested 1e-11 tolerance.
  The corrected pilot uses the established 1e-9 matching tolerance. Both source
  snapshots and logs are preserved; see the campaign README for exact recipes.

The next scientific checkpoints remain spatial and retained-support
convergence, finite-boundary/tail sensitivity, and a consistent slow-growth
formulation. An accurate numerical E/J budget alone does not establish these.
Phase actions remain deferred to the three-independent-cloud extension.

The detailed implementation evidence and other project tracks follow.

- **Multiphase foundations (October 15):** independent-cloud labels, Fourier
  component construction, single-component linear radial assembly, and
  two-phase field reconstruction, projection, nonlinear-source evaluation, and
  full radial RHS assembly are present in the opt-in multiphase library.
  All 159 multiphase checks pass: 12 frequency/azimuthal-label
  checks, 7 radial-operator checks, 8 reconstruction/reference-quadrature checks,
  5 library-projector checks, 5 quartic-source checks, 12 single-cloud-limit checks,
  14 quartic-table checks, 18 phase-origin checks, 13 cosine-convergence checks,
  21 full-RHS checks, 24 matching/integration checks, and 20 Newton matching checks.
  The existing 10 loader checks also passed at
  the linear-assembly checkpoint. The radial checks compare to the single-cloud
  assembly and an analytic Schwarzschild m=0 result, and check signed-frequency/
  azimuthal symmetry and angular ordering. The phase benchmark separates two
  coefficients sharing l=m=3 to 3.11e-16 absolute error; an azimuth-only
  projection returns their sum. The field now comes from the library's
  `reconstruct_two_phase_field_at_radii`, which evaluates one radius per call;
  projection is checked both with an explicit reference sum and with the library's
  `project_two_phase_field`. Library checks also cover an absent component,
  conjugate recovery, and an analytic sine field on unequal phase grids.
  The source checks compare the two driven components to exact integrated
  quartic formulas and verify angular-parity, absent-phase, and zero-field zeros.
  With only q=(1,0) populated, the fundamental and generated third-harmonic
  sources agree with the existing single-cloud implementation for quartic,
  sextic, and cosine potentials; a q_2=2 output vanishes. The full quartic table,
  independent and combined phase shifts, and phase-averaged potential invariance
  now pass. A separate two-cloud cosine example converges under independent
  polar/phase refinement: 32 nodes in each refined direction agree with a 64^3
  reference to at most 1.24e-13 absolute source error. A warmed 32^3 evaluation
  including table setup took about 0.020 seconds and allocated 4.98 MB on this
  host (allocation total, not peak RSS); one real field occupies 256 KiB.
  This establishes quadrature convergence for the tested fixed support, not
  convergence of solved profiles or a universal grid choice.
  The RHS is implemented in `profiles.jl`; all 21 checks in
  `test_nonlinear_rhs.jl` are registered and pass against the actual library.
  These cover the single-cloud limit for all potentials, ordering, state
  preservation, dimensions, driven-source signs/scaling, and zero-field slopes.
  The state stores all R values then all Q values,
  ordered by component and then l; the source is computed once from all R.
  The implemented `two_phase_matching_residual` now connects this RHS to the shared
  integrator, assembles boundary blocks into global [R; Q] order, and returns
  both profile and slope mismatches. All 24 registered checks in
  `test_matching_residual.jl` pass against the library: full single-cloud integration histories for
  all potentials, saved/plain output agreement, joint ordering, input checks,
  and radial refinement in both directions (error-reduction ratios 17.60 and
  14.28 for step sizes 0.4/0.2/0.1, consistent with fourth-order convergence).
  The inserted `solve_two_phase_mode` passes 20 registered checks in
  `test_nonlinear_matching.jl` against the actual library. It fixes the first retained
  l coefficient in each fundamental q block, scales coefficients/residuals by
  the smaller positive outer normalization, and solves for two complex trial
  frequencies plus free boundary coefficients. Frequencies remain on their
  initial positive-frequency bound/outgoing branches.
  The four-component quartic trial at spin=0.99, alpha=0.42, outer normalizations
  [1e-4,1e-5], l=[1,2,0,3], inner=-30, outer=80, match=5, radial step=0.4,
  and quadrature 8x9x9 converges in 3 corrections (scaled residual 13.3874 ->
  1.2248e-10). Both driven outer coefficients are nonzero. This is a matched
  numerical prototype, not an endpoint/angular/radial convergence result or
  a validated slow-growth prescription. The October 15 library run reproduces
  the candidate result; all 159 foundation checks pass (now part of the larger suite).
  A positive-frequency
  boundary wrapper is present; its prototype checks still need adding to the
  repository suite. A converged physical joint-profile family and the action closure
  remain pending. The October 15 sequencing decision keeps total joint-field
  E/J balance for the first two-cloud evolution at fixed BH parameters, with
  full amplitude derivatives and matrix-conditioning checks. Phase actions
  are deferred until the extension to three independent clouds; the profile
  slow-growth approximation still needs to be established. Keep the profile
  solver in `multiphase/`; next validate boundaries and matched-profile
  convergence before physical continuation. Run
  `julia --project=. tests/adiabatic/multiphase/runtests.jl`.
  The first observables increment is
  `two_phase_potential_energy_density` for the joint potential contribution
  per unit tortoise radius. Its 10 checks in `test_potential_energy.jl` pass
  against the inserted helper, loaded explicitly: single-cloud agreement for all potentials,
  phase-origin invariance, an analytic quartic interaction term, and zero/input
  checks. The helper is in `multiphase/observables.jl`; loader inclusion and test
  registration are now complete. These checks are additional to the original 159.
  The complete local tensor-based observables candidate is described below.
  This potential-only check requires no matched profile or
  choice of time-derivative convention for complex trial frequencies.
  The user selected a shared stress-tensor architecture instead of implementing
  each energy contribution separately. Its first candidate helpers,
  `build_polar_theta_derivatives` and `reconstruct_two_phase_field_derivatives`,
  pass 32 checks in `test_field_derivatives.jl`. They reconstruct the real field
  and its Boyer-Lindquist coordinate derivatives (t,r,theta,varphi), with
  Q=dR/dr_* converted to dR/dr. Frozen-field time derivatives require explicitly
  supplied finite real carrier frequencies; complex trial frequencies are
  rejected, and slow amplitude/profile time derivatives are not included.
  Angular derivatives agree with automatic differentiation; all four field
  derivatives agree with an independently differentiated field to <4.5e-16.
  Signed/zero m, repeated angular patterns across distinct q, multiple l in
  one q, input preservation, and zero fields are checked. The helpers have
  been inserted into `multiphase/projection.jl`, including the user-requested
  interface that builds the polar derivative table inside reconstruction from
  `grid`. All 32 checks pass and are registered in the suite.
  The derivative checks do not validate a slow-growth model.
  Two further candidates, `kerr_metric` (for `core/geometry.jl`) and
  `scalar_stress_tensor` (for `multiphase/observables.jl`), pass 115 checks:
  86 metric checks in `core/test_kerr_metric.jl` and 29 tensor/contraction checks
  in `multiphase/test_stress_tensor.jl`. The metric accepts r_* as its position
  input but returns Boyer-Lindquist components in (t,r,theta,varphi) order,
  with signature (-+++). The pointwise tensor returns T_mu_nu/Fa^2.
  Raising its first index gives energy density -T^t_t, angular-momentum density
  T^t_phi, and outward currents -T^r_t and T^r_phi. Test quadrature includes the
  rho^2 spatial measure, 2pi azimuthal factor, phase average, and dr/dr_* only
  for charge densities. Energy agrees with the single-cloud routine for real
  carriers and all three potentials to 3.1e-15 relative error. Joint-field
  fluxes and Schwarzschild angular momentum agree with analytic expressions.
  `scalar_stress_tensor` and `kerr_metric` are now inserted; the 115 checks
  pass against the actual library and are registered. Observables load through
  the multiphase loader. Local checks
  do not establish radial-integrated charges, endpoint convergence, or evolution.
  The inserted `two_phase_observables_at_radius` now passes all 44
  checks in `test_observables_at_radius.jl` against the actual library. It returns energy and
  angular momentum per unit tortoise radius plus both outward currents from
  one quadrature loop. The single-cloud energy agrees to 2.0e-15 relative error
  for three spins and all potentials. Checks also cover joint analytic fluxes,
  Schwarzschild angular momentum, flux reversal, phase-origin invariance,
  exact quartic cross energy, zero fields, input preservation/validation,
  and automatic derivatives of a prescribed local field scaling against finite
  differences. That last check is not a matched-profile amplitude response.
  Observables loader inclusion and test registration are now complete.
  All observables use explicit real carrier frequencies.
  The `two_phase_integrated_observables` candidate passes 49 checks in
  `test_integrated_observables.jl`. It accepts saved matching histories,
  orients each branch toward increasing r_*, requires uniform spacing and a
  common match position, reuses `integrate_uniform` for total E and J, and
  returns horizon/outer energy and angular-momentum losses. Horizon losses
  reverse the inward endpoint's outward-current sign. Checks recover exact
  cubic E/J integrals, including odd interval counts and reversed storage,
  the existing single-cloud integrated energy for all potentials, endpoint
  signs, component reordering, and saved-history compatibility. Nonfinite
  profile samples or local observables are rejected. The first short-domain
  trial fixture overflowed the observables (finite profiles grew to 4.9e247);
  reducing its boundary coefficients by 1e-8 gave a finite dilute fixture.
  Manufactured snapshots and trial histories test bookkeeping, not physical
  matching or convergence. The radial wrapper is now inserted into
  `multiphase/observables.jl`; the local flux fields are named
  `radial_energy_flux` and `radial_angular_momentum_flux` (positive toward
  increasing radius). At the user's request, a tested refactor moves the
  nested helpers to top-level `ordered_radial_branch` and
  `integrate_two_phase_branch_observables`, passing their inputs explicitly.
  All 49 checks passed with that refactor, which is now inserted; comments
  and the radial-branch docstring are also present. The tests are registered.
  Outputs cover only the supplied
  finite radial domain; tails and radial/boundary/truncation convergence remain
  separate. Matched-profile amplitude responses are implemented as described
  in the current checkpoint above.
- **Adiabatic layout refactor (October 12):** branch
  `refactor/adiabatic-layout` groups source and tests by `core/`,
  `single_cloud/`, and `primary_secondary/`. At that checkpoint `multiphase/`
  had only a scope README and core-only loader; the October 15 work above
  supersedes that historical status. Shared operators,
  potentials, Newton/integration routines, quadrature, and sensitivities were
  separated from single-cloud assembly. Reusable workflow functions now load
  from model libraries; model drivers contain command execution. The root
  Julia loader remains available. See the [source guide](src/adiabatic/README.md)
  for current paths and imports, and the
  [implementation plan](docs/multiphase-implementation-plan.md) for the next
  physics milestones. Validation: 509 grouped Julia unit checks plus 10 loader
  checks, 5 Python checks, and the 44-check coupled-balance benchmark pass.
  A nonlinear profile, energy, frequencies, and both complete saved-table
  Fig. 6 trajectories agree exactly before/after the refactor. The QCD grid's
  independent reference-trajectory check retains its 1.71e-8 maximum spin
  difference. The Python wheel builds and contains all 51 Julia source/driver
  files. Numerical campaign artifacts and frozen manifests retain their
  locations and contents; source hashes intentionally reject moved/changed
  implementations. This validates the refactor, not new numerical convergence
  or multiphase physics.
- **QCD grid stopped before disconnect (October 12):** the requested 10x10 log
  grid spans 2–100 solar masses and 1e-14–1e-11 eV, with the final exclusion
  `alpha0 > 0.45`, leaving 83 cells (maximum 0.437617). The user selected a
  10-million-year age cap. Inputs, shared-source hashes, scripts, and restart
  instructions are in the [grid campaign](runs/adiabatic/qcd_grid_20261012/README.md).
  Initial stationary checks at alpha=0.149663 and 0.437617 passed; the
  alpha=0.0149663 check failed at an absolute shooting-residual floor. Array
  26875559 was automatically canceled before starting any production cases.
  The remaining diagnostic 26875867 was canceled at the user's disconnect
  request. No survey jobs or automatic retries remain active. Do not resubmit
  unchanged: the low-alpha conditioning and an end-to-end numerical cell need
  checking first. Forty planned cells are explicitly labeled linear age
  screens, not nonlinear simulations; only two screens have been run locally.
  The new time driver reproduces the saved reference spin curve to 1.71e-8.
- **First coupled QCD spin-down trajectory:** at the user's explicit request,
  a preliminary run now evolves a 20-solar-mass BH from chi=0.99 with
  `m_a=1e-12 eV`, using a nonlinear 211 primary, a 322 secondary on its averaged
  background, and driven m=0,3 scalar channels. The QCD relation fixes
  `f_a=5.7e18 GeV`; small classical seeds are specified in the
  [campaign README](runs/adiabatic/qcd_spindown/README.md).
  A fresh 21-point stationary track at alpha approximately 0.150–0.137
  supplies the leading coupled balance. It reaches chi=0.75 in 33.245 days,
  versus 29.789 days for two linear nonrelativistic clouds with the same seeds;
  the one-year spins are 0.508953 and 0.509951. PNG/PDF plots, trajectory CSVs,
  radial checks, and provenance are saved in that campaign. All 80 stationary
  comparisons pass; time refinement changes spin by at most 1.15e-8.
  Initial secondary horizon-flux radial sensitivity is 0.47%, while the
  primary's is below 5.5e-7 relative. This is a
  trajectory-specific interpolation with refreshed horizon factors, not a
  general background table. GW losses, self-gravity, and the subleading metric
  work/background angular response remain absent. The cloud reaches about
  9.2% of the remaining BH mass; timing is preliminary under that approximation.
  An October 12 field-amplitude audit restores `Phi=fa*phi` along this saved
  trajectory: sampled primary `|Phi|/fa` grows from `5.38e-6` to `9.21e-4`
  (171 times), or `3.06e13` to `5.25e15 GeV`. The cubic/mass ratio stays below
  `1.42e-7`; comparison with the binding scale gives below `2.96e-5`.
  Eight post-processing checks pass; the campaign README records the plot,
  CSV, normalization, and limits of this primary-only amplitude diagnostic.
- **Multimode adiabatic evolution:** all six Sec. III components are implemented
  at $M\mu = 0.42$ and $0.3$, $a/M = 0.99$: averaged primary background,
  secondary eigenproblem and energy, driven $m = 0$ and $m = 3$ fluxes, and
  evolution from energy and angular-momentum balance. The nonperturbative
  Fig. 6 curves agree with the paper to reading precision at $M\mu = 0.42$;
  at $0.3$ the primary peak is about 10% lower and plateaus 5–10% higher.
  Fig. 7 profile peaks at $M\mu = 0.3$ agree to reading precision. See the
  [secondary campaign](runs/adiabatic/multimode_secondary/README.md) for the
  saved tables, comparisons, and checks. Pieces 5–6 are committed in `5b1f322`.
- **Single-cloud adiabatic work:** cosine, quartic, and sextic stationary
  branches are preserved in the [Figure 2 campaign](runs/adiabatic/fig2_20260922/README.md).
  The [amplitude campaign](runs/adiabatic/amplitude_evolution/README.md) has a
  merged 426-point derivative table and a saved cosine trajectory through
  $t/M = 4\times10^7$. Its approach to $A_0 \simeq 12.96268$ is a result on
  the saved spatial grid, subject to the high-amplitude caveat below.
- **Stationary background continuation:** ten primary-only solves at fixed
  $A_0=0.0011390625$ around $(\alpha,\chi)=(0.42,0.99)$ pass twenty observable
  comparisons. Offsets of $\pm0.001$ in alpha/spin converge in two or three
  Newton corrections; two routes to the joint downward offset agree in growth
  rate to $4.01\times10^{-10}$ relative. Halving the radial step changes growth
  by at most $1.91\times10^{-8}$ relative. See the
  [background-continuation campaign](runs/adiabatic/background_continuation/README.md).
  The same ten solves and twenty comparisons also pass at
  $A_0=2.051758324382913$, where the primary peak field is about 1.446 and
  $M\omega_I$ is 6.137 times its dilute value. Half-step changes in imaginary
  frequency and horizon flux are below $1.1\times10^{-8}$ relative; the small
  outward flux changes by $6.8\times10^{-5}$ relative. These are local primary
  stationary checkpoints, not a coupled BH trajectory; angular/endpoint checks
  remain separate work.
- **Secondary and driven background response:** the averaged secondary and
  $m=0,3$ fields now solve on all ten nonlinear-primary configurations above.
  The initial 133/135 passing checks include route agreement, secondary energy
  balance, and driven power identities. Two coarse-grid $m=3$ horizon-flux
  comparisons miss the $10^{-4}$ relative tolerance. Re-solving baseline and
  the joint offset at primary step 0.025 passes all 28 follow-up checks, with
  changes in that flux reduced from about $1.6\times10^{-4}$ to at most
  $5.10\times10^{-6}$ relative. Original failures are preserved. See the
  [response results and refinement](runs/adiabatic/background_continuation/README.md#secondary-and-driven-fields-on-the-finite-amplitude-backgrounds).
  This establishes local scalar response at fixed $A_0$, not BH–cloud evolution.
- **Coupled-balance prototype:** the leading small-cloud-mass scalar-only
  balance is derived in [the balance note](docs/coupled-bh-cloud-balance.md),
  with BH mass/angular momentum and both cloud energies in fixed initial-BH
  units. The campaign's `check_coupled_balance.jl` passes 44 checks against
  the existing nonlinear RHS, dilute channel counting, unit conversion, and
  a changing-background occupation-conservation example. The latter exposes
  the metric-work term required beyond the leading approximation. Actual
  nonlinear metric work and cloud angular-momentum background derivatives
  remain unevaluated; that checkpoint checked instantaneous rates only.
  The separate QCD campaign above now integrates the leading law. See the
  [checkpoint recipe and illustrative RHS](runs/adiabatic/background_continuation/README.md#coupled-balance-checkpoint--october-8-2026).
- **Dolan time-domain work:** the coupled-tower solver, eigenvalue/Leaver
  comparison, and frequency-filter analysis have a
  [linear validation record](docs/dolan-repro.md) and
  [Fig. 13 campaign](runs/dolan/fig13/README.md). The nonlinear source and
  synchronized RK4 have [separate checks](docs/nonlinear-progress.md);
  those do not establish accurate collapse evolution.
- **Self-gravity benchmark:** a standalone Newtonian Poisson check on an
  analytic hydrogenic $211$ density recovers the monopole/quadrupole shift
  $\delta\omega/\mu=-(237/1280)\alpha^2 M_c/M_{\rm BH}$. All 16 checks pass,
  including radial refinement, angular/outer-radius checks, and an independent
  exact radial integral. See the [derivation and command](docs/adiabatic-self-gravity-review.md#211-benchmark-checkpoint--october-7-2026).
  Saved nonlinear profiles and scalar evolution still have no gravity correction.
- **Dilute self-gravity screening:** published hydrogenic self/cross shifts
  have been evaluated on both saved two-cloud energy trajectories. At
  $f_a=10^{16}$ GeV, peak $|\delta\omega_1|/\mu$ is $7.70\times10^{-6}$
  (alpha 0.42) and $5.80\times10^{-6}$ (0.30); the outward scalar threshold
  margin changes by at most 0.0219% and 0.0483%. At $10^{17}$ GeV these
  changes become 2.19% and 4.83%. See the
  [campaign recipe and limits](runs/adiabatic/multimode_secondary/README.md#dilute-self-gravity-estimates--october-7-2026).
  This is an analytic-coefficient estimate, not the nonlinear frequency
  response or a bound on scalar flux errors. BH evolution remains separate.
- **Extensions:** beyond the preliminary coupled trajectory above, GW losses, self-gravity feedback,
  mutual nonlinear cloud response, additional independent clouds, and a PDE
  handoff remain proposals. The [October 6 research note](docs/further_work_1006.md)
  discusses their ordering and scientific motivation.

## Implementation decisions

The [multimode design](docs/multimode-adiabatic-design.md#decisions-and-their-reasons)
records the averaged reading of Eq. (61), restored $\mu^2$, secondary domain
and normalization, and the energy variables used for the balance equations.
The campaign records the $M\mu = 0.3$ tolerance setting, phase test, resonance
scan, and convergence evidence. Use its averaged secondary results;
`cosine_secondary_fixed_time.*` is retained only as a superseded experiment.

## Known caveats

- The Mμ = 0.42 primary branch is **not step-converged above A₁ ≈ 9**: its
  step-0.1 solutions re-evaluated at step 0.05 have scaled residual 4e-5 at
  A₁ = 9.1 rising to 4e-2 at 13. This bears on the cosine flux zero near
  A₁ ≈ 12.96. It does not affect Fig. 6 at Mμ = 0.42, which stays below
  E₁ ≈ 1000 (A₁ ≲ 0.8, step-converged to ~1e-7).
- Unchecked: secondary l-truncation (l ≤ 6) and outer boundary at large A₁.
  The paper reports neither.
- Full high-amplitude angular-quadrature and mode-truncation convergence of
  the primary remains unestablished. Interpolation and time-step checks on
  saved tables do not supply those spatial checks.
- Pieces 4–5: the outer-radius convergence of the $m = 0$ field is unchecked, and
  so is each coefficient against an independent evaluation of the perturbative
  $F_0$ and $F_3$ (Eqs. (38)–(39), (A11), (A14), e.g. by a Green's function).
  The power laws match; the prefactors have not been compared.
- The new static QCD comparisons test the coupled small-secondary limit.
  The dynamical effects of the Sec. III approximations, including neglect of
  secondary self-interaction and its deformation of the primary when both
  populations are substantial, remain unchecked by this campaign.
- The saved two-cloud trajectories include scalar fluxes only. GW power,
  torque, and radiation reaction have not been implemented; physical energy
  normalization alone does not add them.

## Next steps

The October 30 priority is a controlled **linear-profile / Omiya / joint
evolution comparison**, beginning with a small pilot rather than immediately
launching another 36-point campaign. The completed static figures supply
initial checks and diagnostic snapshots; evolution should determine later
populations. No new evolution jobs were submitted with this handoff update.
The [evolution design](docs/qcd-two-cloud-evolution-design.md) supplies the code
reuse map and proposed numerical/physical contract for the steps below.

1. **Audit the existing evolution pieces and define the matched comparison.**
   Reuse the primary/secondary coupled balance and joint response/evolution
   work where validated. Specify common initial BH mass/spin, axion mass/decay
   constant, seed occupations, time origin, and physical units. Keep the QCD
   relation and full cosine. The original population diagnostic used c4=0.3;
   new transfer rates must consistently use the retained cosine (c4=1), not
   silently inherit that diagnostic's coupling. Define a linear-profile
   baseline with perturbative interactions to isolate deformation; any purely
   noninteracting comparator should be separately labeled. Match the dilute
   relativistic growth rates and included channels across approximations so
   growth-rate conventions are not mistaken for mutual-deformation effects.
2. **Validate instantaneous rates and balance before a long trajectory.**
   Check the dilute and one-cloud limits, scalar transfer channels, horizon
   and outgoing fluxes, and energy/angular-momentum accounting. Confirm the
   physical slow-growth closure and amplitude/occupation rate convention;
   tiny matching residuals or differing complex eigenfrequencies alone do not
   establish accurate physical transfer rates. Use focused radial/domain and
   angular/component checks on the rates that control the pilot. The current
   small-secondary solution must connect to a validated finite-population
   joint solve if the second cloud becomes substantial.
3. **Couple cloud populations to BH mass and angular momentum.** Update alpha,
   spin, frequencies, and superradiance thresholds as the BH evolves, with
   clear fixed-reference units and consistent cloud-energy accounting. Do not
   reuse a fixed-Kerr trajectory as a physical stellar-mass QCD prediction.
   Establish whether GW depletion competes on the simulated timescale; include
   it consistently if relevant, or explicitly delimit a scalar-only pilot.
   The earlier catch-up screen makes this especially relevant for late growth
   of the other cloud. Keep omitted self-gravity and subleading metric-work
   effects explicit when interpreting percent-level differences.
4. **Run and benchmark a short matched pilot, then extend through saturation.**
   A natural initial control is alpha=0.3/spin=0.99/M=5 solar masses, with the
   existing static comparison; add a reversed-dominance control once it works.
   Set and document the seed prescription and elapsed-time range for the new
   trajectory. Check a complete step's runtime, balance residuals, time-step
   sensitivity, and fresh stationary samples before extending. Reuse nearby
   solutions or validated interpolation to avoid repeated amplitude ramps;
   do not extrapolate an old trajectory-specific table into a new regime.
5. **Compare dynamical outcomes and then choose expansion points.** Record
   N211(t), N322(t), both energies, BH mass/spin, peak populations, saturation
   times, coexistence, and any quasi-equilibrium. Do not presume large changes
   from mutual deformation. Compare profiles at equal per-cloud energies for
   selected snapshots while allowing the three histories to diverge. Expand
   the parameter scan only after the pilot establishes accuracy and cost.

Supporting work and older unresolved items remain relevant where they affect
that sequence:

- **Architecture and joint response:**
  [revisit the source/derivative architecture](#important-revisit-multiphase-architecture)
  with the user's concern about bloat in view. The older v4 derivative
  discrepancy and physical slow-growth closure require review before reusing
  that evolution path. Static comparison completion does not resolve them.
- **Coupled-BH preparation:** local primary continuation in alpha and spin
  now works at both dilute amplitude and $A_0\simeq2.052$; the secondary and
  driven channels are checked on the nonlinear backgrounds, with the targeted
  refinement described above. Use primary step 0.05 or finer for the tested
  small $m=3$ horizon-flux accuracy requirement.
  The leading balance is now derived and its algebra checked. Its functions
  now live in `src/adiabatic/primary_secondary/coupled_balance.jl`, with units and
  equations documented in the source. The user stopped the time-integrator
  draft to require small, explained implementation steps. Stationary-table
  job 26757461 was cancelled; `runs/adiabatic/coupled_local/` contains paused
  drafts, not a validated trajectory. Review the balance piece together before
  resuming that incremental teaching track. Later on October 8 the user
  explicitly authorized completing the separate QCD spin-down campaign directly;
  its commented driver and Python integrator now provide a first limited
  trajectory with refreshed stationary data, fixed reference units, and explicit
  seeds/$f_a$. That exception does not remove the usual incremental workflow.
  The fixed outer amplitude used in the earlier stationary check is not
  an evolution prescription. Beyond the leading approximation, evaluate
  metric work, background angular response, and consistent gravitational
  energy/boundary accounting as described in the balance note.
- **Perturbative coefficients.** The small-$E_1$ limits of the spline
  coefficients `m0_flux_coefficient` and `m3_flux_coefficient` are the paper's
  $F_0$ and $F_3$ (Eqs. (38)–(39)). An independent evaluation of those (for
  example by the Green's-function formula) would check the prefactors of
  pieces 4–5, which so far are checked only through power laws and the power
  identity. The same coefficients would also give the perturbative curves of
  Fig. 6.
- **$M\mu = 0.3$ details.** Our $E_1$ peak looks ~10% below the figure and the
  plateaus 5–10% above; this is at the limit of reading a log axis. Not
  investigated.
- **Single-cloud saturation:** refine the high-amplitude spatial solve before
  interpreting the saved cosine flux zero or saturation trajectory physically.
- **Historical October 6 proposal:** physical normalization of the saved
  two-cloud trajectories, with GW power/torque estimates and a frozen-profile
  self-gravity diagnostic, followed by a limited coupled-BH calculation
  including GW losses and a matched perturbative
  comparison. See the [October 6 proposal](docs/further_work_1006.md) for the
  rationale and conditional extensions; these are not completed capabilities
  and do not supersede the October 30 sequence above.

## Working with the user

Follow [CLAUDE.md](CLAUDE.md) for incremental development, notation,
validation, and Slurm conventions. The
[piece 4 guide](docs/multimode-piece4-guide.md) is the model for explanations:
define each object, describe the check procedure, and state what it establishes.
Use LaTeX for math in Markdown documentation.

The September 24 handoff is in Git history. The
[September 28 branch-work record](docs/adiabatic-branch-work-2026-09-28.md)
preserves the following stage's experiments; its next-step list is historical.
