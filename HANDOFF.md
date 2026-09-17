# Handoff — axion-sr-sim setup + design doc plan

Written 2026-09-09, end of the session that extracted this repo from the
GRChombo fork. Delete this file once the setup tasks below are done and the
design doc is underway.

## What this repo is

The `radial1d/` Dolan 1212.1477 Sec. III reproduction, extracted (with its
~34 commits of history, via `git filter-repo`) from
`clwelch03/GRChombo-Axion-SR` where it lived at
`Examples/AxionSuperradiance/radial1d/`. Files are hoisted to repo root.

GRChombo is no longer used for this line of work. The fork stays parked for
possible future 3+1D validation but is not the active codebase. Even for an
eventual 2D/3D solver the plan is to write our own, not go back to GRChombo —
though that is a big future decision, not committed.

`waveform/` (old semi-analytic BBH chirp + bilby PE) was deliberately NOT
brought over — judged stale. It stays in the fork.

## Setup tasks still outstanding

1. **Push to GitHub.** Repo `clwelch03/axion-sr-sim` must be created empty on
   github.com (no `gh` CLI on Lawrencium). Then:
   ```
   git remote set-url origin https://github.com/clwelch03/axion-sr-sim.git
   git push -u origin main
   ```
   Auth is HTTPS + personal access token, same as the GRChombo fork.
2. **Tombstone the fork.** In `GRChombo-Axion-SR`:
   `git rm -r Examples/AxionSuperradiance/radial1d`, add a `MOVED.md` pointing
   here, commit. (Do NOT push the rewritten history back to the fork — a stray
   `git push` to the old local origin already got correctly rejected once.)
3. **Write this repo's `CLAUDE.md`.** Carry over from the fork's CLAUDE.md:
   the "Working style on this project" section (cheap decisive test, falsify
   before mechanism), the "Naming style" section (long informative names), and
   the "Build from-scratch reproductions incrementally / hold the user's hand"
   section — that rule governs this whole repo. Additionally - say what you mean.
   Avoid flowery language for the sake of flowery language, and metaphor for the sake of metaphor.
   DROP: Lawrencium module loads, Chombo build config, the GRChombo directory layout, the upstream-merge
   section. Add: how to run the solver + tests, the `runs/` layout.
4. **`runs/` data.** ~105 MB of `.npz` outputs were `cp`-ed into `runs/` and
   `runs/fig13/` and are gitignored (regenerable). The SLURM `.sh` scripts and
   `README.md` ARE committed (commit c86c6b3) — they define reproduction.

## The design doc (docs/DESIGN.md)

**The user is writing this themselves** — they have never written a design doc
and want the ownership/understanding that comes from doing it. Your job is
pointers and review, NOT drafting it. They have already started `DESIGN.md`
(it's tracked, at repo root for now — move to `docs/` when ready).

### What a design doc is (the guidance already given)

- Audience: the user in 6 months, or a cold collaborator who knows the physics
  but not the code.
- It answers *why*, not *what*. Captures decisions invisible in any single
  file. Test: if a sentence stays true after someone rewrites the
  implementation, it belongs here; if it'd need updating, it's a code comment.
- Allowed to drift on details and still be useful. No line numbers, no exact
  constants in it.
- 2–4 pages. One box-and-arrow diagram of the RHS pipeline is worth a page of
  prose.

### Outline (agreed) and where the raw material lives

1. **Scope & non-goals** — README.md "Why this branch exists". Non-goals: no
   3+1D, no backreaction/nonlinear yet, tower is the faithful equation not a
   Leaver gap-closer (README "What we learned").
2. **The physical problem** — solver.py module docstring (coords, units);
   `Tower.__init__` V_l block (Eq. 26); README Dolan section-by-section
   (Eq. 38 coupled system, PML Eq. 43/44, boundaries). One paragraph on why
   the horizon-avoiding slice + tortoise + PML approach at all (the
   no-finite-exact-Pi-seed problem — see the fork's memory
   `pi_seed_inconsistency` / `kerr_schild_reslice`).
3. **Module map** — one line each: solver.py (grid + operator + MoL),
   leaver.py (benchmark, `spheroidal=True` is the apples-to-apples one),
   qbs_eigenmode.py (eigenvalue Γ), measure_gamma.py + frequency_filter.py
   (time-domain Γ, the path that extends to nonlinear), fig13_growth_rates.py,
   analyze_timedomain_run.py, runs/.
4. **Core data model** — the `Tower` object. README "How the PML fits the
   pipeline" is a first draft. State `(ψ,π,χ)` shape `(n_channels, n_points)`,
   channel i ↔ `l_values[i]`, roles of each field (README bullets). What's
   precomputed in `__init__` (all `self.coeff_*`, `self.potential`,
   `self.pml_gamma`) vs rebuilt every `rhs` call (A diagonals, derivatives).
5. **Control flow** — `Tower.rhs` then `rk4_step`, step by step (stencils → b →
   Thomas → +PML source). README numbered list. Then the two extraction paths
   branch: eigenvalue vs time-domain.
6. **Key design decisions** (the heart of the doc) — for each, state the
   rejected alternative:
   - PML as an *additive source* on the full Eq. 25 solver, not a two-region
     solver with seam matching (README "Two design decisions that make this
     clean").
   - Tower-native from day one — `l_max=|m|` is just the one-rung case through
     the identical RHS, so the ℓ-tower can't get orphaned (it did last time).
   - `s = r − r₊` log-space tortoise inversion — catastrophic cancellation at
     deep PML nodes otherwise (solver.py `KerrTortoise.x_of_s` / `r_x`).
7. **Invariants & gotchas** — README "What we learned": ℓ-tower negligible at
   a=0.9–0.99 / these α (ℓ=3 ~4e-4 of ℓ=1); benchmark against SPHEROIDAL
   Leaver not flat Λ=l(l+1); the time-domain *slope fit* is PML-contaminated
   (use eigenvalue for rates, or Dolan's frequency filter). The
   `evolve_metric=0 ⇒ time-independent flux` falsification principle from the
   fork's CLAUDE.md "Working style".
8. **Known limitations / future path** — n=1 overtone 13% deficit
   (frequency-filter bleed, README next steps); hyperboloidal slicing as the
   eventual replacement for PML+wall (README "Future upgrade path").

### Suggested writing order

Sections 3, 4, 5 first (mechanical — they force a code read-through). Then 2
and 6 (the vocabulary is in place). Then 1, 7, 8 (harvest from README + the
fork's memory files).

## Repo layout (updated 2026-09-17)

The code is packaged under `src/dolan/`, `src/adiabatic/`, and `src/common/`.
Dolan contains the existing time-domain and nonlinear work; common contains
Leaver, the hydrogenic benchmark, and plotting. Adiabatic is the new reproduction
of Omiya, Takahashi & Tanaka (2201.04382), starting with the linear radial problem.
See `README.md` for commands and `docs/adiabatic-repro.md` for the next checkpoint.
Run outputs remain in place under `runs/`.

## Working conventions (carry into CLAUDE.md)

- Cheap decisive test over thorough survey. Ask what single measurement
  falsifies a mechanism before proposing it. Prefer tests with an exact
  analytic answer.
- Long informative variable names (`cloud_norm` not `anorm`). Legacy short
  names in old files can stay.
- From-scratch reproductions built incrementally, one piece at a time, user
  following and confirming each step. The user is a physicist and a relative
  novice to numerical methods, using the reproduction to learn the method —
  don't hand back finished modules, don't silently lift code, don't skip the
  hard part for a fast result.
