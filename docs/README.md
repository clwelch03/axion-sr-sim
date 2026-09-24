# Documentation guide

Current entry points, reviewed 2026-09-24:

| Document | Purpose |
|---|---|
| [Project README](../README.md) | Scope, environment setup, and quickstart |
| [Working conventions](../CLAUDE.md) | Development workflow, commands, and module map |
| [Project handoff](../HANDOFF.md) | Where to resume work |
| [Adiabatic design](adiabatic-design.md) | Equations, conventions, and implementation structure |
| [Adiabatic reproduction](adiabatic-repro.md) | Validated checkpoints and remaining work |
| [Dolan reproduction](dolan-repro.md) | Current entry points and dated linear validation record |
| [Dolan design](DESIGN.md) | Time-domain operator and state; unfinished design draft |
| [Campaign index](../runs/README.md) | Submission recipes, saved data, and preservation policy |

## Research and historical records

These documents retain the evidence and decisions from their dates. Their
benchmark timings, job states, and proposed next steps are not current defaults.

- [Nonlinear tower progress](nonlinear-progress.md): Dolan source and RK4 checks.
- [Nonlinear performance](nonlinear-performance.md): NumPy optimizations and Julia bridge experiments.
- [Nonlinear strategy](nonlinear-strategy.md): motivation for the adiabatic approach.
- [Gravitational backreaction review](gravitational-backreaction-review.md): research proposals, not implemented gravity.
- [Figure 2 history](../runs/adiabatic/fig2_20260922/HISTORY.md): continuation experiments and numerical diagnostics.
- [Figure 2 handoff, September 23](fig2-data-collection-handoff.md): earlier diagnostic state, superseded by the campaign README.
- [Setup handoff, September 9](setup-handoff-20260909.md): original repository extraction notes.

## Keeping this current

Update the relevant status guide when an implementation checkpoint changes.
Keep exact commands, defaults, result counts, and restart behavior in campaign
READMEs; record dated measurements in campaign histories. Use saved CSV/JLS
artifacts as evidence for numerical claims and label preliminary results.
When moving code, update command examples and links in the same change.
