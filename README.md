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

### The current 12% gap and how to close it (OPEN)

At a=0.5 the PML gives `Gamma = 0.88x` Leaver (E is perfect). The residual ~12%
under-count is still reflection/mismatch. Leading suspects and the knobs:
1. **The stitch drops a small potential.** At `x_match` the reduced form sets
   `Vc = 0`, but `Vc(r_+ + 0.02) ~ 0.045` is not quite zero -> a small impedance
   step -> reflection. Fix: smaller `eps_match` (push match closer to horizon;
   needs deeper `xmin` so the PML still has room). *Scan in progress:
   `eps_match=0.005, xmin=-90`.*
2. **PML too short / g-ramp too steep.** A short or steep ramp reflects. Fix:
   deeper `xmin` (longer PML), gentler ramp. `gmax` itself was NOT the driver
   (0.5 and 1.0 both gave 0.88x).
3. **Resolution / numerical dissipation of the mode.** Would push |Gamma| the
   *wrong* way (too big), so probably not the 12%; but check `N=8000`.
4. **Extraction window** (`logfit tmin_frac`, run length T). Cross-check by
   varying the fit window.

`scan_pml.py` sweeps (1)-(3); read `scan_pml.log`.

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

Two methods in `evolve.py`, on the probe time series `psi(t, r_probe)`
(probe auto-placed at the cloud peak `2 a0 = 2/alpha^2`):

- `logfit(t, psi)`: `E` from the phase slope, `Gamma` from `d ln|psi|/dt` over the
  late window. **Robust**; this is the number we trust right now.
- `prony(t, psi)`: matrix-pencil fit `psi ~ sum_k c_k e^{s_k t}`, dominant mode ->
  omega. **Currently BROKEN** (returns a spurious mode, e.g. Gamma=-4e-2 at a=0.5;
  81x off). Needs debugging: likely SVD rank/threshold or the Hankel window; also
  should down-select the physical root by nearness to a Leaver anchor
  (`E ~ mu`), not by largest |c| alone. Prony/matrix-pencil is what will get the
  tiny growth-regime Gamma from a moderate window, so **fixing it is on the
  critical path for Stage D.**

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
| C | extract **Gamma** (damped a=0.5) vs Leaver | IN PROGRESS: sponge 0.64x -> PML 0.88x; tuning to close 12% |
| D | superradiant **growth peak** (a=0.99, alpha~0.42, Gamma=+1.5e-7) | BLOCKED on C + Prony fix |
| E | Gamma(alpha) sweep vs Leaver (+ vs 3D DMD) + write-up | not started |

### Chronological log
- **Stage A/B (pass).** Equation transcribed from Dolan Eq. 25, single-mode.
  Freq-domain residual converges to ~1e-4 in the bulk (r>16); large near-horizon
  relative residual is FD/scale (|psi|, Delta -> 0), not a transcription error.
  First evolution (crude friction sponge, a=0.5): E=0.4783 (1.0003x),
  Gamma=-3.33e-4 (0.64x) -- under-counted by sponge reflection.
- **Stage C (in progress).** Swapped friction sponge -> Dolan PML. a=0.5:
  E=0.478255 (1.0001x), Gamma=-4.60e-4 (**0.88x**). Prony returned garbage (skip
  for now, use logfit). gmax scan: 0.5 and 1.0 both give 0.88x -> gmax not the
  driver; testing smaller eps_match + deeper xmin + finer N (`scan_pml.py`).

## 8. Immediate next actions (for whoever picks this up)
1. Finish the PML tuning (`scan_pml.py`): get Gamma -> ~1.00x at a=0.5 by
   reducing `eps_match` (+ deeper `xmin`) and/or lengthening the PML. Update the
   table above.
2. Fix `prony()` (physical-root selection + SVD threshold) and cross-check it
   against `logfit` on the a=0.5 case.
3. Stage D: run the growth peak (a=0.99, alpha=0.42). Gamma is +1.5e-7, so
   `logfit` needs either a very long run or (better) a working Prony on a moderate
   (~1e4-1e5 M) window. Confirm the **sign flips to +** (growth) and matches
   Leaver -- this is the actual deliverable: a scalar superradiant growth rate.
4. Stage E: sweep alpha at a=0.9 (and a=0.99) -> Gamma(alpha) vs Leaver, and
   overlay the 3D DMD point (a=0.5 damped, 0.92x) to show the 1+1D tool reaches
   where 3D cannot. Write up.

## 9. Gotchas
- Conventions: M=1, mu=alpha; time dependence e^{-i omega t} so
  omega = i * (log-eigenvalue)/dt; growth is Im(omega) > 0.
- The single-mode c^2 truncation is good to <1% for a QBS but is an
  approximation; Sec. 2 says how to add l=3,5 if a benchmark ever needs it.
- Run outputs must go to this folder (persistent), not the session /tmp.
