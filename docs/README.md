# Documentation guide

Start with [the project handoff](../HANDOFF.md) for current capabilities,
open caveats, and next steps. Entry points reviewed 2026-10-12:

| Document | Purpose |
|---|---|
| [Project README](../README.md) | Scope, environment setup, and quickstart |
| [Working conventions](../CLAUDE.md) | Development workflow, commands, and module map |
| [Codex chat recovery](codex-chat-recovery.md) | Pull missing chats and messages into the current node at the start of a session |
| [Project handoff](../HANDOFF.md) | Authoritative project status and next steps |
| [Adiabatic design](adiabatic-design.md) | Equations, conventions, and implementation structure |
| [Adiabatic source guide](../src/adiabatic/README.md) | Physical-model folders, dependency boundaries, loaders, and current commands |
| [Adiabatic reproduction](adiabatic-repro.md) | Validation records and observable conventions |
| [Multimode adiabatic design](multimode-adiabatic-design.md) | Sec. III reproduction (2211.01949): equations, decisions, component map, extension ideas |
| [Multiphase implementation plan](multiphase-implementation-plan.md) | Proposed mutual-deformation and three-cloud milestones, closure checks, and repository/branch organization |
| [Multiphase radial RK4 design](multiphase-rk4-design.md) | Flowcharts for one radial RK4 step and its joint nonlinear RHS, state ordering, and function map |
| [Coupled BH–cloud balance](coupled-bh-cloud-balance.md) | Fixed reference units, leading scalar-only coupled rates, metric work and angular response, and algebraic checks |
| [QCD spin-down campaign](../runs/adiabatic/qcd_spindown/README.md) | First preliminary coupled BH/211/322 trajectory and comparison with linear nonrelativistic evolution |
| [Piece 4 guide](multimode-piece4-guide.md) | Explainer for the driven $m = 0$ field: physics, equation, solve method, checks |
| [Secondary campaign](../runs/adiabatic/multimode_secondary/README.md) | Secondary-cloud runs, results, and checks |
| [Dolan reproduction](dolan-repro.md) | Current entry points and dated linear validation record |
| [Dolan design](DESIGN.md) | Time-domain operator and state; unfinished design draft |
| [Campaign index](../runs/README.md) | Submission recipes, saved data, and preservation policy |

## Research and historical records

These documents retain the evidence and decisions from their dates. Their
benchmark timings, job states, and proposed next steps are not current defaults.

- [Multiphase work, October 23](multiphase-work-2026-10-23.md): direct Fig. 7 profile setup, source/Jacobian optimization experiments, partial library implementation, and running-job results.
- [Nonlinear tower progress](nonlinear-progress.md): Dolan source and RK4 checks.
- [Nonlinear performance](nonlinear-performance.md): NumPy optimizations and Julia bridge experiments.
- [Nonlinear strategy](nonlinear-strategy.md): motivation for the adiabatic approach.
- [Gravitational backreaction review](gravitational-backreaction-review.md): research proposals, not implemented gravity.
- [Adiabatic self-gravity review](adiabatic-self-gravity-review.md): literature and proposed metric response of the nonlinear stationary branch.
- [Multimode literature review, September 28](multimode-literature-review-2026-09-28.md): existing results, recent papers, and candidate gaps in nonlinear multimode evolution.
- [Further adiabatic work, October 6](further_work_1006.md): priorities for BH spin evolution, decay-constant dependence, additional modes, and a PDE handoff.
- [Figure 2 history](../runs/adiabatic/fig2_20260922/HISTORY.md): continuation experiments and numerical diagnostics.
- [Branch work, September 28](adiabatic-branch-work-2026-09-28.md): driver consolidation, dense interpolation checks, and the high-amplitude extension; its next-step list is historical.
- [Figure 2 handoff, September 23](fig2-data-collection-handoff.md): earlier diagnostic state, superseded by the campaign README.
- [Setup handoff, September 9](setup-handoff-20260909.md): original repository extraction notes.

## Keeping this current

Update `HANDOFF.md` when an implementation checkpoint or project priority changes.
Keep project status there; other entry points should link to it. Keep working
conventions in `CLAUDE.md` and equations and implementation decisions in design docs.
Keep exact commands, defaults, result counts, and restart behavior in campaign
READMEs; record dated measurements in campaign histories. Use saved CSV/JLS
artifacts as evidence for numerical claims and label preliminary results.
When moving code, update command examples and links in the same change.
