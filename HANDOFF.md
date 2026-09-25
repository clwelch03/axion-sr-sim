# Project handoff

Updated 2026-09-24 against the working tree and saved campaign tables.

The active adiabatic work is preparing amplitude evolution from stationary
energy and flux tables. Implicit energy derivatives are saved at all 178 cosine
branch points, and a single-point finite-difference check supports the derivative
calculation on the current spatial grid. Time integration remains unimplemented.
The preliminary three-potential Figure 2 data remain available separately.
The Dolan time-domain solver remains available for linear validation and
nonlinear transients; gravitational backreaction is not implemented.

At the September 24 closeout, all 12 independent interpolation reference points
have completed. Cubic Hermite energy derivatives meet the 0.1% target at 8/12
points; amplitude rates using linear flux interpolation meet that threshold at
4/12. Maximum discrepancies are 0.157% and 1.635%, respectively. These sampled
checks do not establish accuracy throughout the branch or spatial convergence.
All jobs for this check have finished, and no further run is needed tonight.

Suggested next step: use the saved state derivatives to compute flux amplitude
derivatives, then compare Hermite flux interpolation against the same 12 saved
references. Guide the user through this incrementally. Energy interpolation
also needs attention in sparse intervals; changing flux interpolation alone
will not resolve that error. See the
[amplitude evolution campaign](runs/adiabatic/amplitude_evolution/README.md)
and its [comparison table](runs/adiabatic/amplitude_evolution/results/interpolation_checks.csv).

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
