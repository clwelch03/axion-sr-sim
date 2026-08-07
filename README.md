# 1+1D time-domain massive-scalar QBS on Kerr

A dimensionally-reduced time-domain solver for the scalar quasi-bound-state
(superradiant) complex frequencies `omega = E + i*Gamma` on Kerr, following
**Dolan 2012, "Superradiant instabilities of rotating black holes in the time
domain", arXiv:1212.1477 (PRD 87, 124026)**. This doc is a self-contained
handoff: a fresh reader should be able to continue the build from it.

---

## 1. Why this exists (the motivation, don't lose it)

The 3+1D GRChombo code **cannot** get the scalar growth rate. Established this
exhaustively (see `../linear_rates/` + memory `kerr_schild_reslice`):

- Superradiant `Gamma ~ 1e-7..1e-8` sits ~3 orders **below** the 3D
  **numerical-dispersion floor (~1e-4)**.
- That dispersion scrambles the coherent *ingoing* radial phase structure the
  horizon energy transfer depends on: the seeded QBS sheds its ingoing winding
  (11 -> 1 rad) within ~150 M and relaxes to a near-standing wave, so every 3D
  method that reads the near-horizon transfer (extraction-sphere flux, apparent-
  horizon mass drift, discrete-operator eigenvalue) is corrupted at the ~1e-4
  level. Amplitude growth is separately dead (e-fold ~1e7 M).
- The dispersion is worst where the coherent wave is strongest (the cloud peak);
  it is NOT outer-resolution or AMR-reflection (both falsified). Beating it in 3D
  needs ~450x more compute per level. Infeasible.

**The escape is dimensional reduction.** The |2,1,1> mode is exactly
`Y_11(theta) R(r) e^{-i omega t}`. Separate the (analytically known) angular part
and evolve only the radial direction. In 1D you can afford 100-1000x the radial
resolution, so dispersion drops orders of magnitude and `Gamma ~ 1e-8` becomes
accessible. This is how the field actually computes these: Dolan 2012 (time
domain, what we are reproducing) and Leaver/Dolan 2007 (frequency domain, our
reference `../linear_rates/leaver.py`). GRChombo's job is the nonlinear/binary
physics that has *no* separated form; the linear rate was always the wrong job
for it.

---

## 2. The equation

Field ansatz `Phi = (1/r) psi(t,x) Y_lm(theta) e^{i m phi}`, with `x` the
Boyer-Lindquist tortoise coordinate,

```
dx/dr = (r^2 + a^2)/Delta,   Delta = r^2 - 2 M r + a^2,   M = 1,   mu = alpha.
```

`Y_lm` here is a plain spherical harmonic (Dolan's convention, azimuthal part
split out). Dolan's full system Eq. (25) couples l -> l+-2 via the metric
`cos^2 theta`; the odd-l and even-l sectors decouple, so `l=m=1` couples only to
`l=3,5,...`. The coupling strength is `~ a^2 Delta c_{l,l+-2}^m`, i.e. O(c^2)
with `c^2 = a^2(omega^2 - mu^2)`, which is **tiny for a QBS** (omega ~ mu):
`c^2 ~ -a^2 alpha^4/n^2 ~ -1.6e-3` at a=0.9, alpha=0.3. `leaver.py` uses exactly
this approximation (Lambda_lm = l(l+1)). So we run the **single-mode l=m=1
truncation**, carrying only the diagonal `c_ll = <11|cos^2 theta|11> = 1/5`:

```
D psi_tt = (r^2+a^2)^2 psi_xx + [2 i a m (r^2+a^2) - 2 a^2 Delta / r] psi_x
           - 4 i a m r psi_t - Vc psi

D  = (r^2+a^2)^2 - (1 - c_ll) a^2 Delta
Vc = Delta [ l(l+1) + (2/r)(1 - a^2/r) + 2 i a m / r + mu^2 (r^2 + c_ll a^2) ]
```

- `-4 i a m r psi_t` is the **frame-dragging / superradiance** term. `psi` is
  genuinely complex; a real formulation would miss the instability.
- First-order-in-time form (what we evolve): `psi_t = pi`,
  `pi_t = ( A2 psi_xx + A1 psi_x + At pi - Vc psi ) / D`, with
  `A2=(r^2+a^2)^2`, `A1=2iam(r^2+a^2)-2a^2 Delta/r`, `At=-4iamr`.
- Horizon-regular in these coordinates: at `Delta->0`, `D=A2=(r^2+a^2)^2` so wave
  speed^2 = A2/D -> 1, and `Vc->0`. At large r, `Vc -> -(mu^2 r^2 + l(l+1))` gives
  the massive dispersion (bound states for omega < mu).

**To extend to full accuracy** (if the single-mode c^2 truncation ever limits
us): add the l=3 (and l=5) channels with Dolan's coupling coefficients
`c_{jl}^m = (1/3) delta_{jl} + (2/3) sqrt((2j+1)/(2l+1)) <j2m0|lm><j200|l0>`, on
the RHS terms `+ a^2 Delta c_{l,l+-2} psi''_{l+-2}` (LHS) and
`- a^2 mu^2 Delta c_{l,l+-2} psi_{l+-2}` (RHS). Not needed at current precision.

---

## 3. The Perfectly Matched Layer (PML) — the crux of getting Gamma right

### Why a plain sponge fails, and why this is the same physics as the 3D trap

A quasi-bound state is a stationary balance: the cloud continuously feeds a
coherent **ingoing** wave into the horizon, which the horizon absorbs. `Gamma`
IS that net absorption. If the inner boundary **reflects** any of the ingoing
wave, the reflected outgoing part interferes to form a (partial) standing wave,
the net absorption drops, and the measured `|Gamma|` comes out **too small**.
This is *exactly* the 3D standing-wave pathology — except in 1D we can make the
absorber arbitrarily good. So **absorber quality is literally the Gamma test**:
a crude friction sponge gave `Gamma = 0.64x` Leaver; the PML gave `0.88x` and
climbing. The sign of the error (always too small) is the fingerprint of
reflection, not dissipation (dissipation would make |Gamma| too *big*).

### What a PML is

A Perfectly Matched Layer is an absorbing region appended at a boundary, designed
so that waves enter it with **zero reflection at the interface** (unlike a sponge,
whose impedance mismatch reflects) and are then exponentially attenuated inside.
It is realized by an analytic continuation / complex coordinate stretching of the
wave equation that damps propagating modes while leaving the interface impedance
matched. The cost is one or more **auxiliary fields** carrying the stretched-
coordinate history.

### Dolan's near-horizon PML (what we implement)

Deep in the tortoise coordinate (`r -> r_+`, where `V -> 0`) the full equation
reduces to a constant-coefficient advective wave equation,
`psi_tt = psi_xx + 2 i m Omega_H (psi_x - psi_t)`. Dolan builds the PML for this
reduced form, with one auxiliary field `chi` and an attenuation `g(x)`:

```
psi_t = pi
pi_t  = psi_xx + 2 i m Omega (psi_x - pi) - chi_x - i m Omega chi - g (pi + i m Omega psi)
chi_t = g (psi_x - chi + i m Omega psi) - i m Omega chi
```

with `Omega = Omega_H = a/(2 r_+)`. **Key property: with `g = 0` these equations
are identical to the full equation near the horizon**, so the region can be
stitched smoothly to the interior full-coefficient equation at a match radius
`x_match = x(r_+ + eps_match)` where `Delta` is already small.

### How it's wired in `evolve.py`

- `pml = x < x_match`; `x_match = tortoise_x(r_+ + eps_match, a)`.
- `g(x)` ramps quadratically 0 -> `gmax` from `x_match` down to `xmin`.
- Each RHS eval computes BOTH the interior momentum RHS `dpi_int` (full,
  r-dependent A2/A1/At/Vc/D) and the PML RHS `dpi_pml` (reduced form + chi + g),
  and selects with `np.where(pml, dpi_pml, dpi_int)`. `chi` is evolved only where
  `pml` (zero elsewhere).
- Outer boundary: the bound state (omega < mu) decays exponentially, so the field
  is ~0 at large r; a mild outer sponge `gout` near `xmax-40` is there only for
  numerical safety.

### PML parameter glossary — what each knob physically is

The tortoise map sends the horizon to `x -> -infinity` (`x(r) = r + ...
ln(r-r_+) ...`, Sec. 4), so there are two different "boundaries" near the
horizon that are easy to conflate: the **physics switch** (interior equation
-> PML equation) and the **domain truncation** (where the grid stops). The
four knobs below control those two independently.

- **`eps_match`** sets the **physics switch**. `x_match = tortoise_x(r_+ +
  eps_match, a)`: the radius, in areal `r`, at which `rhs()` stops using the
  full interior coefficients (`A2,A1,At,Vc,D`, r-dependent) and starts using
  Dolan's reduced near-horizon form (`Vc` implicitly `0`, advective in
  `Omega_H`). That reduced form is only exact in the limit `r -> r_+`; at
  finite `eps_match` the true `Vc(r_+ + eps_match)` is small but nonzero
  (`~0.045` at `eps_match=0.02`), so dropping it at the stitch is a small
  impedance mismatch that reflects part of the ingoing wave. **Smaller
  `eps_match` = switch happens closer to the true horizon = smaller dropped
  `Vc` = less reflection at the stitch.** This is suspect #1 for the 12% gap.
  Cost: because `x` diverges logarithmically as `r->r_+`, shrinking
  `eps_match` pushes `x_match` sharply more negative, eating into the room
  available for the absorber (next point) unless `xmin` also goes deeper.

- **`xmin`** sets the **domain truncation** — where the grid actually stops
  (a hard wall; nothing is imposed on `psi` there, it just runs out of
  points). `[xmin, x_match]` is the full extent of the PML. A wave that isn't
  fully attenuated by the time it reaches `xmin` reflects off the truncation
  and comes back as spurious outgoing signal, again undercounting absorption.
  **Deeper (more negative) `xmin` = longer PML = more room for `g(x)` to ramp
  up gradually = less residual amplitude reaching the wall.** Cost: for fixed
  `N`, a deeper `xmin` spreads the same number of points over a longer domain
  and coarsens `dx` (and hence `dt = cfl*dx`) everywhere, including inside the
  PML — which is why the `eps_match=0.005, xmin=-90` config is paired with an
  `N=8000` run, to check the result isn't secretly limited by that coarsening.

- **`gmax`** is the peak of the attenuation profile `g(x)`, which ramps
  quadratically from `0` at `x_match` to `gmax` at `xmin` (`evolve.py:36-38`).
  It's the "how hard do we damp" knob, and it cuts both ways: too small and
  the wave isn't attenuated before hitting `xmin` (undercounts, as above); too
  large and `g(x)` rises so steeply near `x_match` that it becomes a
  discretization-scale impedance jump in its own right (a stiff PML reflects
  just like an unmatched sponge does). **The `gmax=0.5` vs `2.0` scan gave
  identical `0.883x`**, so on this grid neither failure mode is active in that
  range — `gmax` is *not* the current bottleneck, which is why the follow-up
  scan targets `eps_match`/`xmin` instead.

- **`N`** (with `xmin`,`xmax` fixed) sets resolution, `dx=(xmax-xmin)/N`. Two
  distinct roles: (1) resolves the bulk wave (coarse `dx` disperses/damps the
  signal — but the wrong way for our 12% gap, since numerical dissipation
  makes `|Gamma|` too *big*, not too small, so this is suspect #3, lower
  priority); (2) resolves the steep rise of `g(x)` inside a short/aggressive
  PML (under-resolving *that* looks like reflection, in the same direction as
  the 12% gap, which is why deeper-`xmin` configs get paired with higher `N`).

All four interact through one physical picture: **any sharp, under-resolved,
or approximated feature in the near-horizon region turns part of the ingoing
QBS wave around before it reaches the true horizon, and a partially-reflected
ingoing wave always reads as less absorption, i.e. `|Gamma|` too small.**
That's the observed sign of the residual 12% error, which is why reflection
(not dissipation) is suspect #1, and why the fix is "push the switch and the
wall both closer to ideal while making sure resolution doesn't become the new
limiter."

### CALIBRATED PML-width rule (the reusable engineering result)

The single most important tuning result, calibrated on the a=0.9,alpha=0.5
damped scan (`scan_pml.py`/`scan_pml.log`) and stated in the physically
correct units so it *transfers to other (a, alpha)*:

> **Make the PML at least ~2.5 local near-horizon wavelengths wide**, where
> `lambda = 2*pi / |E - m*Omega_H|` and `E` is the mode's real frequency.
> Width here means `x_match - xmin` (the extent of the g-ramp).

Evidence (width in units of `lambda=38.1 M` at this benchmark):

| PML width | / lambda | Gamma / Leaver |
|---|---|---|
| 58.7  | 1.54 | 0.883x |
| 74.2  | 1.95 | 0.949x |
| 101   | 2.65 | 0.955x |
| 163   | 4.3  | 0.952x |
| 190   | 5.0  | 0.952x |

Big jump 1.5->2λ, essentially converged by ~2.5λ, flat thereafter. `gmax`
(0.5/1.0/2.0 all 0.883x at fixed narrow width) and grid resolution (N=4000 vs
8000 identical) are NOT the drivers; **width-in-wavelengths is.** This is just
Dolan's "vary g smoothly over a scale of a few wavelengths" made quantitative
for our grid.

**Consequence / warning for the superradiant regime:** `lambda` *diverges* as
`E -> m*Omega_H` (the superradiant threshold), so the required absolute width
blows up: at the a=0.99 growth peak (alpha=0.42) `E-m*Omega_H = -0.025` gives
`lambda = 251 M`, needing `xmin ~ -750` for 2.8λ -- 7x wider than the damped
benchmark. Right at threshold no finite PML suffices; this is the intrinsic
hard limit near the top of the instability band (see Stage D log). The
absorption *depth* is never the issue (a quadratic `gmax=1.5` ramp over 700 M
integrates to ~350, i.e. `exp(-350)` attenuation -- absurd overkill); it is
always the smoothness/impedance-match length that binds.

**Status of the old "12% gap":** CLOSED, and it was never mainly PML. With
width>=2.5λ the PML gets `Gamma` to 0.955x (logfit); the remaining ~5% was a
`logfit`-through-beating artifact (Sec. 5, `floor_check.py`), and the fixed
anchor-`prony` recovers 0.98-1.00x. Sequence of suspects that were falsified:
gmax (no), grid resolution (no), grid-edge treatment (no, Dolan-pin test),
PML reflection (no, Fig.3 replication) -> real cause was extraction (beating).

---

## 4. Numerics

- Uniform grid in tortoise `x`. `r(x)` by inverting the analytic
  `x(r) = r + (2 r_+/(r_+-r_-)) ln(r-r_+) - (2 r_-/(r_+-r_-)) ln(r-r_-)` via a
  400k-point spline (`build_grid`). NB horizons `r_+- = 1 +- sqrt(1-a^2)` (M=1).
- 4th-order centered stencils for `psi_x`, `psi_xx` (`d_dx`, `d2_dx2`).
- RK4 in time, `dt = cfl * dx`, `cfl = 0.4`.
- **No numba on this cluster** (grviz/gwenv envs lack it). Pure-numpy RK4 is the
  bottleneck: N=4000, T=3000 ~ 2 min; N=12000 timed out at >280 s. For the long
  growth-peak runs, use Prony on a *moderate* window (below) rather than watching
  5e6 M of growth; if raw speed is needed, port the stepper to C or add numba to
  an env.

---

## 5. Extraction of omega

Three methods in `evolve.py`, on the probe time series `psi(t, r_probe)`
(probe auto-placed at the cloud peak `2 a0 = 2/alpha^2`):

- `logfit(t, psi)`: `E` from the phase slope, `Gamma` from `d ln|psi|/dt` over the
  late window. Robust in the sense of never crashing, but (see `floor_check.py`
  finding below) it fits a single straight line through what turns out to be a
  genuinely *beating* two-mode signal, biasing the fitted `Gamma` low
  (0.955x Leaver at the a=0.9,alpha=0.5 benchmark).
- `prony(t, psi, anchor=...)`: matrix-pencil fit `psi ~ sum_k c_k e^{s_k t}`
  with up to `p` simultaneous modes -- built to separate exactly the kind of
  beating `logfit` can't. **FIXED** (was badly broken: argmax|c| mode
  selection picked whichever mode had the largest amplitude in a given
  window, not the physical one, e.g. Gamma 81x off). Fix: widened the Hankel
  embedding (`L=4p`), tightened the SVD truncation, and -- the key change --
  select the recovered mode **nearest a reference `anchor` omega** (pass
  Leaver's `w_ref`, already computed every run) instead of largest `|c|`.
  Tested against the a=0.9,alpha=0.5 beating data (same run `floor_check.py`
  examined): **Gamma = 0.98-1.00x Leaver across three different fit-window
  starts** (`tmin_frac`=0.1/0.35/0.6), a real improvement over `logfit`'s
  0.955x, and verified NOT just snapping to the anchor -- the selected
  candidate sits close to Leaver on its own, with a clearly distinct
  second-nearest candidate (`E~0.49, Gamma~-1e-4..+2e-5`, much less damped)
  that's almost certainly the slow overtone `logfit` was beating against.
  `E` from `prony` still has a small (~2%) offset from Leaver not yet
  understood (Gamma is the much tighter match, and Gamma is what we actually
  need); worth a closer look before leaning on `prony`'s `E` for anything.
  This resolves most of the Stage-C "0.95x plateau" mystery: it was mostly a
  `logfit`/beating artifact, not fundamentally a PML/physics deficiency.

### 5a-bis. Dolan's frequency-filter method (his Sec. IV.A.2) — IMPLEMENTED & VALIDATED (`freqfilter()`)

Our `prony()` is the matrix-pencil route to mode separation. Dolan uses a
different but equivalent-in-spirit route, the **frequency-filter**, and
crucially he uses it to solve *exactly the beating problem we hit* -- his
words (Sec. IV.2): the massive-field spectrum "generates 'beating' effects
which disguise the instability." His recipe (verbatim structure, his eq. 49):

1. FFT the probe time series `psi(t, r_probe)` -> power spectrum.
2. Peak-find the modes; interpolate for a better real-frequency estimate
   `omega_hat_n` of each.
3. For each mode: multiply the FT by a narrow **super-Gaussian bandpass**
   centred on that peak, `F(x) = exp(-x^4 / w^4)` with width `w = k*d_omega`,
   `k ~ 5` bins, `d_omega = 2pi/T` the FFT resolution.
4. Inverse-FFT the filtered spectrum -> the reconstructed single-mode
   contribution `psi_n(t)`.
5. Growth rate = linear regression of `ln|psi_n(t)|` -- BUT only trust the
   reconstruction near `t_max/2`; the ends are corrupted by the filter's
   finite support (windowing edge effect).

**Implemented as `freqfilter(t, psi, E_peak, kbins=5)` in `evolve.py`, and it
turned out to be NEEDED, not just backup.** The Stage-D floor test showed
`prony`'s Gamma is *resolution-limited* for the ultra-slow superradiant growth
(resolving `Gamma~1e-7` means distinguishing the pencil eigenvalue `|lam|`
from 1 at the `~1e-8` level, below its finite-window precision) -- it nails
`E` but returns Gamma noise ~`1e-4`. `logfit` meanwhile fits through the
n=0/n=1 beat. The frequency-filter is the tool that gives a clean per-mode
Gamma at accessible T. **Validated on the a=0.9,alpha=0.5 damped benchmark
(`freqfilter_validate.py`): at Dolan's `kbins=5` it gives Gamma = 1.000x
Leaver, and E to 0.008%** -- both better than prony there (prony's E was ~2%
off). One implementation gotcha (fixed): `psi ~ e^{-iEt}` lives in the
`np.fft` bin at `omega = -E`, not `+E`, so the bandpass must be centred at
`-E_peak` (a first version centred at `+E_peak` returned negative E and 2-4x
Gamma -- the sign of the bug).

Why it fails *differently* from `prony` (so both are worth having): `prony`
separates modes by their *time-domain* eigenvalues (sensitive to run length
and SVD conditioning); the frequency-filter separates by *spectral* peaks
(sensitive to FFT resolution `2pi/T` and peak overlap). Dolan's **Fig. 13** is
the proof it reaches our exact Stage-D target: his time-domain frequency-filter
points land on the frequency-domain Leaver curve for the a=0.99 `n=0` growth
rate across `M*mu = 0.30..0.44`, *including the peak at 0.42* -- with
`t_max = 1e5 M` and `dx = M/6`. (He pulls the `n=1` overtone out the same way.)

### 5b. Flux readout (Dolan Sec. IIC, eq. 34) — internal-consistency check, not independent

`logfit`/`prony` both extract `Gamma` from the *field amplitude* at one probe
point (the cloud peak). Dolan's paper defines a second, physically distinct
observable: the **energy flux through a fixed near-horizon surface**, derived
from the conserved stress-energy (his eqs. 27-36). For a purely ingoing wave,
`Phi_l ~ A_l exp(-i omega(t+x))`, this flux is

```
F_E = (r^2+a^2) Re( Phidot_l^* . Phitilde'_l ),   Phitilde'_l = Phi'_l + i a m Phi_l/(r^2+a^2)
```

(his eq. 34, single `l=1` term here since we truncate at one mode). Dolan's
`Phi_l = r psi_l` is exactly our `psi` (our ansatz already factors the `1/r`
out), so `Phidot_l = pi`, `Phi'_l = psi_x` — no new state, just a second
readout point.

**Caution: this is NOT an independent validation against Leaver, and the
first draft of this doc overstated it.** `flux_gamma` and `logfit` are both
computed from the *same evolved field* of the *same run* — not a second
simulation or a different physical channel. Once the transient dies off, a
genuine quasi-bound state is (to good approximation) a single mode,
`psi(x,t) ~ A(x) e^{-i*omega*t}`: the temporal decay rate `Gamma=Im(omega)` is
then, by construction, the *same number at every x*, only the spatial envelope
`A(x)` differs. So if the numerical solution has relaxed into one coherent
mode, `logfit` (cloud peak) and `flux_gamma` (`x_match`) are guaranteed to
report the same `Gamma` **whether or not that `Gamma` is correct** — agreement
does not confirm the value against Leaver; only Leaver itself (a genuinely
separate frequency-domain calculation) does that.

**What agreement (or disagreement) actually tests:** whether the solution has
cleanly relaxed to a single coherent mode at all, evaluated at two points with
potentially different admixtures of any leftover contamination (a spurious
PML-reflection component, numerical noise) that hasn't fully decayed relative
to the true QBS mode. If such contamination had a different weight at the
cloud peak vs. `x_match`, the two late-window log-slope fits would plausibly
disagree. The first real run gave `logfit`=0.955x, `flux`=0.946x (~1% apart) —
a useful negative result (rules out a location-specific fit artifact / obvious
multi-mode contamination) but explicitly *not* evidence the 0.95x itself is
close to converged; that still has to come from comparison to Leaver and from
things like the PML-width scan and Sec. 5c's Fig.-3 replication.

Implemented as `flux_gamma(t, phi, phidot, phix, r, a)` in `evolve.py`:
records `psi, pi, psi_x` at `x_match` on every recorded step (alongside the
existing cloud-peak probe), builds `F_E(t)`, and since `F_E` is quadratic in
the field (`~ exp(2 Gamma t)` for a quasi-stationary state) fits
`d ln|F_E|/dt / 2` over the same late window as `logfit`. Printed automatically
by `evolve.py`'s CLI as a third row (`flux :  ...`) alongside `logfit`/`prony`.
Full `C(E)`/`C(J)` conservation (his eq. 36, needs the volume energy-density
integral eq. 30-32) is *not* implemented — bigger lift, not needed just for
this internal-consistency check, revisit only if a stronger code-correctness
test is wanted.

### 5c. Outstanding: replicate Dolan's PML validation plots (Fig. 3) (OPEN)

Dolan validates the PML with three specific plots (his Fig. 3, a=0.99), not
just a scalar Gamma comparison — worth replicating since they show **where**
in `x` or `t` a discrepancy would show up, which a single Gamma number can't:

- **(a) snapshot `Re(psi(x))` at fixed `t`, no-PML vs with-PML, linear scale.**
  Should agree outside the PML region and show smooth attenuation (no
  ringing/reflection) inside it. "No-PML" is free in our setup: the interior
  equation is horizon-regular (Sec. 3 above — `Vc->0`, wave speed^2->1 at
  `Delta->0`), so `--gmax 0` reduces the PML branch to the identical interior
  form (README Sec. 3: "with g=0 these equations are identical to the full
  equation near the horizon") and evolves right up to the truncated wall
  without needing any absorber for *stability* — only for *absorption*. So
  this is just two runs (`gmax=0` and `gmax=<current best>`) with
  `nsnap>0` (already supported by `evolve.py --nsnap`) diffed at one `t`.
- **(b) same snapshot, log scale.** Shows the PML run decaying to the noise
  floor inside the PML while the no-PML run keeps ringing at O(1) — the
  clearest single picture of "is the absorber actually absorbing."
- **(c) time series at a fixed `x` outside the PML, no-PML vs with-PML,
  log scale, run long enough to see them diverge.** This is the sharpest test:
  the two should be *indistinguishable* until the no-PML run's reflection off
  the truncated wall has had time to propagate back out to that `x` (light-
  crossing time from `xmin` to the probe); after that they should visibly
  split. Where they split (and how much the with-PML curve *doesn't* split)
  is a direct, model-independent demonstration that the PML is doing its job
  in the exterior region we actually extract `Gamma` from.

Needs: a `--gmax 0`-equivalent no-PML run (already possible via the CLI, no
code change) + the existing `nsnap` field-snapshot capability for (a)/(b), and
the existing probe time series for (c) (pick a probe outside `x_match`, e.g.
reuse the cloud-peak probe or add a dedicated one at Dolan's `x_1`-style
location). Not started.

---

## 6. Files & how to run

```
radial1d/
  radial1d.py   grid, r(x), coefficients (coeffs), Stage-A freq-domain residual
  evolve.py     PML time-domain evolver (run) + prony/logfit + CLI (main)
  scan_pml.py   PML tuning sweep (gmax, eps_match, xmin, N) vs Leaver at a=0.5
  README.md     this file
  *.log         run outputs (pml_a050.log, scan_pml.log, ...)
```

Env + run:
```
GWENV=/global/scratch/projects/pc_heptheory/clwelch/software/envs/gwenv/bin/python
cd Examples/AxionSuperradiance/radial1d
$GWENV radial1d.py                              # Stage-A residual check
$GWENV evolve.py --alpha 0.50 --a 0.90 --T 3000 --N 4000 --gmax 1.0
# long/slow runs: launch with nohup ... & and log to THIS folder (persistent);
# NOT the session /tmp scratchpad (ephemeral). Consider tmux for the claude
# session itself so monitors survive SSH drops.
```
Reference values from `../linear_rates/leaver.py` (M=1, mu=alpha, l=m=1, n=2):
`leaver.solve(alpha, a)` returns complex omega. Key targets:
- damped validation: a=0.9, alpha=0.5 -> E=0.478203, Gamma=-5.2050e-4.
- **growth peak (Stage D target): a=0.99, alpha~0.42 -> E~0.409, Gamma~+1.5e-7**
  (Dolan's ~1.72e-7 at a~0.997, M mu~0.45; e-fold ~5.8e6 M).

---

## 7. Validation ladder & status

| stage | test | status |
|---|---|---|
| A | freq-domain residual: plug Leaver (omega,R), check L[psi]=0 | DONE: bulk 1e-4, converges |
| B | evolve, extract **E** (damped a=0.5) vs Leaver | DONE: **E to 0.03%** |
| C | extract **Gamma** (damped a=0.5) vs Leaver | NEARLY DONE: sponge 0.64x -> PML 0.88x -> eps_match/xmin tuning **0.95x** (logfit, beating-biased) -> fixed `prony` gives **0.98-1.00x** across 3 fit windows. Small (~2%) `E`-from-prony offset unexplained; not blocking |
| D | superradiant **growth peak** (a=0.99, alpha~0.42, Gamma=+1.5e-7) | **DONE: freqfilter +1.603e-7 = 1.056x Leaver, RIGHT SIGN (growth), E to 5 digits; stable over kbins=3/5/8** (cluster job 24715584, T=1e5) |
| E | Gamma(alpha) sweep vs Leaver (+ vs 3D DMD) + write-up | IN PROGRESS: adaptive per-alpha sweep on cluster; 3/5 pts in (0.50/0.55 freqfilter 0.99x, 0.60 edge-case 0.25x); 2/5 (0.40/0.45) failed on cluster (node fault, zero output), resubmitted job 24784194 |

### Chronological log
- **Stage A/B (pass).** Equation transcribed from Dolan Eq. 25, single-mode.
  Freq-domain residual converges to ~1e-4 in the bulk (r>16); large near-horizon
  relative residual is FD/scale (|psi|, Delta -> 0), not a transcription error.
  First evolution (crude friction sponge, a=0.5): E=0.4783 (1.0003x),
  Gamma=-3.33e-4 (0.64x) -- under-counted by sponge reflection.
- **Stage C (in progress).** Swapped friction sponge -> Dolan PML. a=0.5:
  E=0.478255 (1.0001x), Gamma=-4.60e-4 (**0.88x**). Prony returned garbage (skip
  for now, use logfit). gmax scan: 0.5, 1.0, 2.0 all give 0.883x -> gmax not the
  driver. Boundary-pin test (explicit Dolan-style `psi=0=pi` at the two grid
  edges, replacing the implicit zero-derivative edge treatment) gave an
  **identical 0.883x** -> edge handling isn't the source either; PML was already
  attenuating to negligible amplitude before reaching the wall. Scanned
  eps_match/xmin instead:
  `eps=0.02,xmin=-70`(baseline)=0.883x -> `eps=0.005,xmin=-90`=0.949x (N=4000
  and N=8000 identical, resolution-converged) -> `eps=0.002,xmin=-120`=0.955x.
  Realized the driver is **PML width in units of the local near-horizon
  wavelength** `lambda=2pi/|E-m*OmH| ~ 38` (Dolan requires the g-ramp vary
  smoothly over "a few wavelengths"): widths were only 1.54/1.95/2.65 lambda.
  Pushed to 4.3 and 5.0 lambda (`eps=0.001,xmin=-190/-220`) -> **both give
  0.952x, a plateau** — no further gain past ~2.5-3 lambda. So PML
  width/eps_match is likely **exhausted as a lever**; remaining ~5% gap is an
  open question. Added a `flux_gamma` cross-check (Sec. 5b) and scoped
  replicating Dolan's Fig. 3 PML-validation plots (Sec. 5c) as the next
  diagnostics, since the scalar Gamma comparison alone doesn't localize where
  the residual reflection/error is coming from.
- **PML mechanism independently confirmed correct (Fig. 3 replication,
  `fig3_pml_validation.py` + `fig3_dolan_match.py`).** Reproduced Dolan's
  actual Fig. 3 using his own parameters (a=0.99, x0=10M, sigma_x=5M,
  near-massless field, x_match=-100 via eps_match=1.4e-6): exterior agreement
  between no-PML/with-PML is clean, the with-PML floor bottoms at ~1e-8 near
  the wall matching his reported floor, and the two stay indistinguishable
  over his full t<300M window -- all matching his own validation figure
  closely. (First attempt used our massive a=0.9,alpha=0.5 benchmark params
  and did NOT match his figure's shape -- traced to using a *massive* field,
  which traps radiation in a bound-state well and produces persistent
  bulk oscillations his massless test never shows; switching to alpha~0
  fixed it.) **This rules PML reflection out as the source of the ~5% gap.**
- **Floor check (`floor_check.py`) found the real mechanism: mode beating, not
  a floor.** Plotted `|psi(t)|` and sliding-window local decay-rate slope
  over the FULL T=3000 production run (best config). The amplitude envelope
  is NOT a clean single exponential -- it beats, with a ~1000-1500M period,
  all the way to t=3000 (humps/dips in `|psi(t)|`, panel 1 of
  `floor_check.png`). The local slope of `log|psi|` oscillates between roughly
  -0.004 and +0.003 (~10x the true |Gamma|~5.2e-4) and randomly signed
  *throughout* the post-transient window, not just late -- i.e. the field has
  NOT relaxed to a single coherent QBS mode by t=3000; it's a superposition of
  (at least) two comparable-lifetime modes (fundamental + a slow overtone,
  most likely) beating against each other. `logfit`/`flux_gamma`'s *global*
  fits happen to average close to the true rate over a long-enough baseline
  (hence 0.95x, not something wildly off), but they're riding on top of real
  interference, not an asymptotic single-mode regime. This is the generic
  broadband-Gaussian-excitation issue discussed in Sec. 8 below (a shaped seed
  or a longer run would reduce/outlast it) -- and it directly motivates fixing
  `prony()` (built exactly to separate multiple simultaneous exponentials)
  as the next concrete step, rather than further PML tuning.
- **`prony()` fixed (`prony_test.py`).** Root cause was the mode-selection
  step (`argmax|c|`), not the SVD/Hankel machinery -- with 2+ comparable-
  amplitude modes beating (previous entry), "biggest amplitude in this
  window" is not reliably the physical fundamental. Fix: widened the Hankel
  embedding (`L=4p`, was `2p`), tightened the SVD singular-value threshold,
  and select the recovered mode **nearest a reference anchor omega** (pass
  Leaver's `w_ref`) instead of largest `|c|`. On the a=0.9,alpha=0.5 beating
  data: **Gamma = 0.98-1.00x Leaver across tmin_frac = 0.1/0.35/0.6** (vs
  `logfit`'s 0.955x on the same data) -- confirmed not a selection-bias
  artifact, since the chosen candidate sits close to Leaver on its own and a
  clearly distinct, much-less-damped second candidate (`E~0.49,
  Gamma~-1e-4..+2e-5`) is visible nearby in the full candidate list (almost
  certainly the beating overtone). `E` from prony still off by ~2% from
  Leaver, unexplained -- Gamma (what we actually need) is the tight match.
  **This resolves most of the Stage-C "0.95x plateau": it was mostly a
  logfit/beating artifact, not a PML/physics deficiency.**
- **Stage D attempt 1 FAILED, cleanly diagnosed (a=0.99, alpha=0.42 growth
  peak; `stage_d_growth_attempt1_tooShortPML.log`).** Got Gamma ~ -6e-6..-8e-6
  from anchor-prony: **WRONG SIGN** (want +1.52e-7) and ~40-55x too big. Two
  independent under-sizings, both now understood via Dolan's Fig. 13 caption
  (`t_max=1e5 M`, `dx=M/6`) and our own calibrated PML-width rule (Sec. 3):
  (i) **PML far too narrow** -- used `xmin=-250`, width 208 = only **0.83
  lambda** at this near-threshold `lambda=251 M` (vs the >=2.5λ rule); a short
  PML reflects the ingoing wave into a spurious *decaying* standing component
  that swamps the tiny true growth (hence wrong sign, floor-level magnitude).
  (ii) **Run 5-10x too short** -- `T=2e4` gives only `Gamma*T ~ 3e-3` of
  log-amplitude change for a `1.5e-7` rate, far too little to stand above the
  systematic floor; Dolan needs `1e5 M` for ~1.5e-2 of log-amplitude range.
- **Stage D floor test (`stage_d_floortest.py`) -- PML FIX WORKS.** Proper PML
  (`xmin=-750`, 2.82λ) + `dx=M/6`, `T=1e4`: `logfit` gave **Gamma = +2.6e-8
  (RIGHT SIGN, growth!)** with **E=0.408797 matching Leaver to 6 digits** --
  vs attempt-1's -5e-7 wrong-sign. So the sized-up PML dropped the systematic
  floor ~2 orders (from ~5e-7 to ~1e-8) and revealed the true growth. Two
  extraction lessons here: `prony` returned Gamma ~ -5e-4 (noise) -- it's
  **Gamma-resolution-limited** for `1e-7` growth (E fine, Gamma not); `logfit`
  read 6x low (`+2.6e-8` vs `+1.5e-7`) because **n=0 and n=1 both grow and
  beat**. => the clean per-mode Gamma needs Dolan's frequency-filter, now
  implemented+validated (Sec. 5a-bis).
- **Stage D DONE -- the Tier-1 deliverable (cluster job 24715584,
  `stage_d_full.py`, `dx=M/6`, 2.82λ PML, T=1e5, 58 min on lr6).** A scalar
  superradiant GROWTH rate from the 1+1D time-domain code:

  | method | E | Gamma | ratio |
  |---|---|---|---|
  | Leaver (truth) | 0.408797 | +1.5183e-7 | 1.000x |
  | **freqfilter kbins=5** | 0.408806 | **+1.6031e-7** | **1.056x** |
  | freqfilter kbins=3/8 | 0.4088 | +1.591/1.602e-7 | 1.048/1.055x |
  | logfit | 0.408805 | +8.16e-8 | 0.537x (beats w/ n=1) |
  | prony | 0.407677 | -4.08e-7 | -2.69x (Gamma-res-limited) |

  **RIGHT SIGN (growth), matches Leaver to ~5%, E to 5 digits, stable across
  filter widths.** Reproduces Dolan's Fig. 13 peak point. This is precisely
  the rate the 3D code could not reach (its dispersion floor ~1e-4 vs this
  ~1e-7; memory `kerr_schild_reslice`). The ~5% overshoot is at the level of
  Dolan's own Fig. 13 scatter; possible residual finite-T / n=1-leakage /
  PML systematic, not chased. Confirms the full method chain: horizon-
  penetrating BL-tortoise eq. + wavelength-sized PML + Dolan frequency-filter.
  Series saved `stage_d_full.npz` for reuse.
- **Stage E, fixed-window sweep (`stage_e_sweep.py`) -- only alpha=0.5 valid,
  as expected in hindsight.** Swept alpha=0.30/0.40/0.50/0.60 at a=0.9 with
  ONE (T, window) choice. `|Gamma|` spans **6 orders of magnitude** across
  that band (+1.5e-8 at 0.30 -- itself already past the superradiant
  threshold at a=0.9 -- up to -8.1e-3 at 0.60), so a single run length can't
  serve all of it: alpha=0.30/0.40 are below the floor for T=3000 (need a
  long run like Stage D); alpha=0.60 decays to noise before the late fit
  window (needs an EARLY window instead). Only alpha=0.50 landed in the
  sweet spot (freqfilter 1.000x). This directly motivated an **adaptive**
  sweep where every knob (T, PML width, window) is derived per-alpha from
  the Leaver Gamma rather than fixed.
- **Stage E, adaptive sweep (`stage_e_adaptive.py`+`.sbatch`, SLURM array
  24716620), 3/5 points in.** Per-alpha: `lambda=2pi/|E-m*OmH|`,
  `xmin=x_match-2.7*lambda`, `T=clip(5/|Gamma|, 2500, 6e4)`. Results so far
  (freqfilter / logfit / prony ratios to Leaver):

  | alpha | Leaver Gamma | freqfilter | logfit | prony |
  |---|---|---|---|---|
  | 0.50 | -5.205e-4 | **0.99x** | 0.87x | 1.06x |
  | 0.55 | -2.684e-3 | **0.99x** | 0.18x | 0.12x |
  | 0.60 | -8.061e-3 | 0.25x | 0.12x | 0.12x |

  0.50/0.55 confirm the adaptive design works and that `freqfilter` is the
  extractor to trust (logfit/prony both fail once the field decays into
  noise before their window -- exactly the failure the adaptive T is meant
  to dodge, though clearly the window/T still isn't quite right for these
  two either). **0.60 underperforms (0.25x)**: at `|Gamma|=8e-3` the mode's
  own spectral peak is Lorentzian-broad (width ~|Gamma|) and starts to
  compete with the bandpass width (`kbins*d_omega`), so the filter clips the
  mode's tails and biases Gamma low -- the fast-decay edge of the
  frequency-filter's own range, not a bug.
- **2 of 5 adaptive tasks (alpha=0.40, 0.45) FAILED on the cluster with ZERO
  output** -- ran 20-40 min then vanished with no `.out`/`.err` content at
  all, pointing to a node-level fault (a real script error would still flush
  partial output) rather than a bug. Resubmitted as job `24784194`; as of the
  next session (2026-08-07) still `PENDING` with reason `ReqNodeNotAvail,
  Reserved for maintenance` -- **lr6 was down for a maintenance window**, not
  just normal queue depth (`sinfo` showed lr6/lr5 nodes in `maint` state).
  **Moved all three `.sbatch` scripts (`stage_d_full`, `stage_e_sweep`,
  `stage_e_adaptive`) from `partition=lr6` to `partition=lr4`** (same
  `lr_normal` QOS via the `pc_heptheory` account, confirmed with
  `sacctmgr show assoc`; lr4 had 65 `alloc`-state nodes, no maint flag).
  Cancelled the stuck job, resubmitted alpha=0.40/0.45 as job `24784511` on
  lr4 -- queued normally this time (no maint block). If lr4 ever shows the
  same symptom, check `sinfo -p <partition>` for a `maint`/`down*` state
  before assuming it's a code bug.

## 8. Immediate next actions (for whoever picks this up)
1. (Optional polish) Understand `prony`'s ~2% `E` offset, and consider
   whether `p`/`L`/SVD threshold need further tuning -- not blocking, since
   Gamma is already the tight match and Gamma is the deliverable.
2. **Stage D: growth peak (a=0.99, alpha=0.42, Gamma=+1.5e-7) -- IN PROGRESS,
   attempt 1 failed, now properly parameterized.** Use Dolan's Fig.13 params
   as the target spec: `dx=M/6`, `t_max` up toward `1e5 M`, PML width >=2.5λ
   (=251 M here, so `xmin ~ -750`). Extract with anchor-`prony` (or implement
   Dolan's frequency-filter, Sec. 5a-bis, if prony stumbles). Gate on the
   floor test (chronological log above) before committing the full `1e5 M`
   run. Confirm the **sign flips to +** (growth) and matches Leaver -- this is
   the actual deliverable. Feasibility caveat: `1e5 M` at `dx=M/6` over a
   ~1000 M domain is ~1.5e6 pure-numpy RK4 steps (~tens of min, may need a
   detached `nohup`/`tmux` run or a C/numba port -- README Sec. 4).
3. Stage E: sweep alpha at a=0.9 (and a=0.99) -> Gamma(alpha) vs Leaver, and
   overlay the 3D DMD point (a=0.9, alpha=0.5 damped, Gamma=0.92x Leaver --
   see memory `kerr_schild_reslice`) to show the 1+1D tool reaches
   where 3D cannot. Write up.

## 9. Gotchas
- Conventions: M=1, mu=alpha; time dependence e^{-i omega t} so
  omega = i * (log-eigenvalue)/dt; growth is Im(omega) > 0.
- The single-mode c^2 truncation is good to <1% for a QBS but is an
  approximation; Sec. 2 says how to add l=3,5 if a benchmark ever needs it.
- Run outputs must go to this folder (persistent), not the session /tmp.
