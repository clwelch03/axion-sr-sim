# Multilevel self-interaction scaffold

This directory defines the next Tier-1.5 extension: self-interaction-induced
transfer between distinct quasi-bound levels.  The current implementation has
the mode catalog, angular coupling tensor, and a cheap coherent amplitude
demonstration; calibrated radial transition rates are not implemented yet.

The important distinction from the existing `(l,m=1)` tower is that physical
level transitions require different frequencies and, in general, different
azimuthal sectors.  Adding more angular `l` channels at fixed `m=1` produces
same-frequency sidebands, not the full level-transfer problem.

## Mode labels

Every bound channel must carry the full hydrogenic/Leaver label

```text
(principal_n, ell, azimuthal_m)
```

Here `principal_n` is the principal-quantum-number label accepted by
`linear_rates/leaver.py` (`n=2` for `211`), while a radial overtone in a
time-domain signal is not automatically a fixed-principal-`n` state.  A mode
catalog must store the label, complex Leaver frequency, and radial eigenfunction
together; the frequency alone is not enough to identify a channel.

The first target sector is

```text
211  = (n=2, ell=1, m=1)
322  = (n=3, ell=2, m=2)

The first allowed angular demonstration also includes `100 = (1,0,0)` as an
`m=0` mediator.  A pure `211` condensate cannot directly source `322`, since
`1 - 1 + 1 = 1`, whereas `211 + 100 -> 322` satisfies the selection rule.
```

with explicit open channels for horizon absorption and continuum emission.
Additional bound levels should be added only after this sector passes its
selection-rule and rate tests.

## Quartic selection rule

For the number-conserving rotating-wave cubic, write

```text
Psi = sum_a psi_a(t, x) Y_(ell_a,m_a)(theta,phi) / r.
```

Projection of `Psi Psi* Psi` onto an output channel `j` uses

```text
Q[j,a,b,c] = integral Y*_j Y_a Y*_b Y_c dOmega,
```

and vanishes unless

```text
m_j = m_a - m_b + m_c.
```

The angular tensor is necessary but not sufficient: each allowed term also
needs a radial overlap involving the correctly normalized Leaver radial
functions.  The radial convention must be derived and tested against the
`radial1d`/Dolan field normalization before rates are interpreted physically.

## Bound modes versus open sinks

The bound-only quartic system is Hermitian.  It can exchange occupation and
phase coherently between bound modes, but it cannot by itself produce an
irreversible saturation rate.  The implementation must therefore distinguish:

1. **Bound-bound mixing:** coherent four-wave coupling between retained
   quasi-bound modes, with detunings built from their distinct real
   frequencies.
2. **Horizon/continuum loss:** non-Hermitian depletion rates or fluxes for
   channels that leave the retained bound-state basis.

The schematic target process is `211 + 211 -> 322 + horizon(m=0)`, together
with the reverse/population and continuum channels required by the chosen
effective model.  The horizon/continuum terms must not be hidden inside an
angular overlap tensor or inferred from a closed Hermitian tower.

## Required validation gates

Before a physical sweep, the implementation should pass these cheap tests:

1. **Label gate:** each channel's Leaver frequency and radial profile agree
   with the requested `(n,ell,m)` label.
2. **Selection gate:** forbidden `m` combinations remain exactly zero; allowed
   angular overlaps have the expected symmetry and normalization.
3. **Bound-only gate:** with open rates disabled, the retained-mode system
   conserves its Hermitian norm/energy to integrator accuracy.
4. **Sign gate:** adding a positive horizon/continuum loss rate depletes the
   corresponding bound occupation; reversing that sign fails the test.
5. **Reduction gate:** disabling all interlevel couplings reproduces the
   existing single-level `211` equations and growth rate.
6. **Minimal transition gate:** the `211`/`322` sector exhibits the predicted
   allowed transfer while a forbidden-sector control remains unpopulated.

Only after these gates should we compare a multilevel saturation rate with the
single-level detuning and adiabatic-spin-down baselines.  The physical
motivation is the self-interaction-induced level exchange and scalar emission
identified in [Baryakhtar et al., arXiv:2011.11646](https://arxiv.org/abs/2011.11646),
not a claim that the current RWA tower has already measured those processes.

**All six gates pass (`gate_test.py`, 7/7 including the split reduction
check):**

1. **Label:** `mode_catalog.leaver_frequency` returns finite Leaver
   frequencies for the `211`/`322` labels (`gate_test.run_label_gate`).
2. **Selection:** `angular_selection_test.py` -- signed-`m` overlaps match
   the historical `m`-only API, the `211+100->322` channel is nonzero, the
   `211+211->322` (no `100`) channel is exactly zero, and the signed-`m`
   theta factors preserve spherical-harmonic normalization.
3. **Bound-only:** `sign_and_hermiticity_test.run_bound_only_test` -- with
   every explicit loss rate zero, the three-mode `211/100/322` system
   conserves total norm to `8e-13` over `T=300`.
4. **Sign:** `sign_and_hermiticity_test.run_sign_test` -- an explicit
   `-loss_rate * A` sink (added to `transition_demo.run`, independent of the
   Hermitian coupling tensor) depletes the sourced channel for positive
   `loss_rate` and grows it for negative `loss_rate`, matching
   `exp(-2*loss_rate*t)` to `1e-6` relative.
5. **Reduction:** `reduction_test.py` (profile: reduces to
   `radial1d.eigenmode_seed` to machine precision) and
   `evolution_reduction_test.py` (full RK4 evolution: `linear_driver.run`
   with the active set `[211]` matches `evolve_tower.run(l_max=1)` to
   `~1e-15`).
6. **Minimal transition:** `minimal_transition_test.py` -- the
   `211+100->322` sector populates `322` to `~3e-3`, while the
   forbidden-sector control (`211` alone, no `100` mediator) leaves `322`
   at exactly zero population.

Run the full suite with `$GWENV multilevel/gate_test.py`.

## Current files and cheap transition check

* `mode_catalog.py` stores explicit `(n,ell,m)` labels and Leaver frequencies.
* `couplings.py` builds the mode-indexed angular tensor `Q[j,a,b,c]` and its
  radial-weighted form.  The latter uses the current reduced-field convention
  `I = integral dx u_j* u_a u_b* u_c / r^2` with peak-normalized Leaver
  profiles.
* `transition_demo.py` integrates the bound-only amplitude equations using
  real Leaver frequencies and those radial overlaps, with an optional
  explicit `loss_rates` sink (used only by the sign gate above; disabled by
  default so the demo itself stays Hermitian).  With its default small-`100`
  seed at `alpha=0.30`, `chi=0.99`, the initially empty `322` channel reaches
  population about `3e-3` in the cheap demonstration run.
* `sign_and_hermiticity_test.py`, `minimal_transition_test.py` implement
  gates [3]/[4]/[6] above.

This last result is a radial-overlap/coherent-transfer test only.  The radial
profiles still use peak normalization rather than an energy/adjoint
normalization, and the model still omits metric backreaction; the sink used
by the sign gate is an artificial diagnostic knob, not yet the physical
horizon/continuum rate (that requires deriving the open-channel flux from the
Leaver `Im(omega)` per retained mode, or an explicit horizon-absorption
boundary condition, and normalizing it consistently with the bound-bound
tensor).
