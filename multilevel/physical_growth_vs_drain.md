# Physically calibrated growth vs. drain

This note reports the first *calibrated* test of whether the self-interaction
level-transfer channel (`211+100 -> 322`, see `transition_demo.py`) can compete
with real superradiant growth of `|211>`. It replaces the earlier hand-tuned
demo (`plot_transition_sink_demo.py`, `Gamma_322=1.5` chosen purely to make a
drain visible) with each mode's own Leaver `Im(omega)`, with no tuning.

## Step 0: finding a parameter window with the intended sign pattern

At the existing demo defaults (`alpha=0.30`, `chi=0.99`), `MODE_322` is not
actually a linear sink: `Im(omega_211)=+2.69e-8` (growing, as expected) but
`Im(omega_322)=+1.79e-12` is also positive (growing, not damped). The intended
picture -- `211` superradiant, `322` a genuine horizon sink -- does not hold
there.

A scan at `chi=0.99` over `alpha in [0.05, 0.80]` found: `211` is superradiant
(`Im(omega)>0`) for `alpha <~ 0.44` and damped above that; `322` is damped
(`Im(omega)<0`) only for `alpha <~ 0.158`, turning superradiant above that. The
only overlap window with both signs as intended is `alpha <~ 0.158`. Used
`alpha=0.15`, `chi=0.99`:

| mode | Im(omega) | sign |
|---|---:|---|
| 211 | +1.836409e-10 | growing (superradiant) |
| 322 | -4.335580e-16 | damped (sink), but ~6 orders of magnitude smaller than 211's rate |
| 100 | -2.788465e-04 | damped -- dominates the native linear-rate landscape by ~6 orders of magnitude over the other two |

## Step 1: implementation and gate tests

Added an opt-in `use_leaver_linear_rates=False` keyword to
`transition_demo.run()`. When `True`, each channel gets its own
`+Im(omega_j)*A_j` growth/decay term, composing additively with any explicit
`loss_rates`. Default (`False`) behavior is unchanged.

Full gate suite re-run after the change, all pass: `angular_selection_test.py`,
`reduction_test.py`, `evolution_reduction_test.py`,
`sign_and_hermiticity_test.py`, `minimal_transition_test.py`, `gate_test.py`
(7/7).

## Step 2: the comparison (alpha=0.15, T=3000, native interaction_strength=1000)

- Pure-linear control (`lambda=0`): fitted `d ln|A_211|/dt = 1.8364e-10`,
  matching `Im(omega_211)` essentially exactly (sanity check).
- Full native-coupling run: fitted net `d ln|A_211|/dt = +2.841e-08` -- about
  155x larger in magnitude than the bare linear rate, but **still positive**
  (net growth, not suppression), and negligible in absolute terms: e-folding
  time ~3.5e7 M, versus `|100>`'s own e-folding time of ~1800 M.
- The cubic coupling drives visible Rabi-like sloshing (`|322>` peaks at 0.12%
  population, `|100>` at 0.21%) but it is oscillatory, not a secular drain --
  confirmed by both a direct log-slope fit and an independent envelope fit
  (period ~41M, from 73 resolved `|322>` peaks), agreeing to 3 sig figs.

## Conclusion

At native coupling strength, in the only parameter window where the intended
sign pattern holds, the self-interaction channel does **not** produce a
measurable secular drain competing with superradiant growth -- both the bare
growth rate and the coupling-driven net rate sit many orders of magnitude
below any physically accessible timescale. The dominant native linear effect
in this three-mode system is actually the `|100>` mediator's own decay, not
the intended `|322>` sink -- `Im(omega_322)` is ~6 orders of magnitude smaller
than `Im(omega_211)` at this alpha, so `322` is a very weak sink to begin
with. This is a real, useful negative result: it means the earlier hand-tuned
demo's dramatic drain required a `Gamma_322` roughly 15 orders of magnitude
larger than the physical value, and gives no support (at these parameters) for
this channel measurably correcting the naive single-mode chi_sat/Regge
picture.

**Caveats:** only one parameter window was checked (the only alpha range with
the required sign pattern at `chi=0.99`); the radial overlaps still use peak
(not energy) normalization per the outstanding item in `README.md`; and
`interaction_strength=1000` is carried over from the earlier demos, not
independently derived.

Figure: `physical_growth_vs_drain.png` (left: `|A_211|^2` fractional deviation,
linear-only vs. linear+native-coupling; right: log-scale bar comparison of all
native rates, colored by sign, including the extracted net rate).
