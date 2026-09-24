# Project handoff

Updated 2026-09-24 against the working tree and saved campaign tables.

The active adiabatic work is the preliminary three-potential Figure 2
reproduction. Linear matching, nonlinear continuation, polynomial potentials,
energy integrals, and finite-radius flux diagnostics are implemented. Slow
amplitude evolution and high-amplitude radial/angular convergence remain open.
The Dolan time-domain solver remains available for linear validation and
nonlinear transients; gravitational backreaction is not implemented.

Start with:

- [README](README.md): setup, commands, and project scope.
- [Documentation index](docs/README.md): current guides versus historical notes.
- [Adiabatic status](docs/adiabatic-repro.md): implemented checks and remaining work.
- [Figure 2 campaign](runs/adiabatic/fig2_20260922/README.md): saved tables, restart settings, and plotting.
- [Working conventions](CLAUDE.md): incremental development and numerical checks.

Load Julia library code through `src/adiabatic/load.jl`. The stationary solver
is in `src/adiabatic/stationary/`, energy integration in `observables/`, and
reusable drivers in `drivers/`. Campaign recipes stay under `runs/`.
Check `git status` before editing: this is a shared workspace.

Remaining work should distinguish branch collection from convergence. A small
matching residual does not validate spatial discretization or mode truncation.
Keep source and energy potential selections consistent, preserve amplitude
order when plotting flux versus energy, and use only one writer per checkpoint.
Heavy calculations go through Slurm; a documentation refresh requires no new
production run.

The [original setup handoff](docs/setup-handoff-20260909.md) and
[September 23 Figure 2 handoff](docs/fig2-data-collection-handoff.md) are historical
records. Their setup tasks, solver defaults, and next-step lists are superseded
by the current guides and campaign recipes.
