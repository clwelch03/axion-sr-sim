# Brief: can a stiff integrator run the spin-down evolution past the |211⟩ plateau?

For an agent starting fresh in this worktree (`/global/u1/c/clwelch/axion-sr-sim-rewrite`, branch
`rewrite/axioncloud`, Julia package `AxionCloud`). Written 2026-11-12. Read `CLAUDE.md` first for working conventions.
The user is a physicist learning numerical methods; they type library code themselves.

## The question

`evolve_with_spindown` (`src/evolution.jl`, around line 650) evolves the state `[ln A₂₁₁, ln A₃₂₂, b, j]`. Here
`A` are the mode amplitudes, `b = M/M₀`, and `j = J/M₀²`. It uses Tsit5, an explicit Runge–Kutta method. Each
right-hand-side (rhs) evaluation is a nonlinear radial solve (`solve_radial`) of about 10 s.

Once |211⟩ saturates, the hole parks at the |211⟩ superradiance threshold. The cloud and spin then sit on a plateau
that relaxes quickly, while the |322⟩ growth (Im ω₃₂₂ ≈ 1.1e-12) and the hole evolve slowly. Tsit5's step stays near
1e8 M₀ the whole time. Reaching the |322⟩ stage (≳ 1e13 M₀ at α = 0.3; 1e15–1e17 M₀ at α = 0.187) therefore takes
weeks to years of CPU.

**Hypothesis (not yet tested):** the step limit is explicit-method stability on a stiff fast mode, not accuracy. If
so, an implicit or Rosenbrock method would take steps set by the slow evolution, perhaps ≫ 1e9 M₀.

Your job is one cheap decisive test of this hypothesis, with **no edits to `src/`**.

## Evidence already in hand

- Job 61102851 (cosine potential, a₀ = 0.9, α = 0.3, f_a/M_Pl = 0.466873, i.e. f_a = 5.7e18 GeV, target t = 1e13)
  timed out after 17 h at t = 2.10e11.
  - Log: `results/spindown_a0.9_alpha0.3_fa0.466873_t1e13-61102851.out`. Each line starting `t = ` is one rhs
    evaluation, with its t, A, Im ω, Newton iterations, b, j, and w. The lines include Tsit5 stages and rejected
    trials, not only accepted states.
  - Between t = 1.5e9 and 1.16e10 there are 564 rhs lines, about 94 Tsit5 steps (about 1.6 h).
  - The plateau holds A₂₁₁ = 2.7968e-4, a = j/b² = 0.86941, and b = 0.981478 throughout.
  - A₃₂₂ grows freely: 2.723e-9 at 1.5e9 and 3.444e-9 at 2.10e11, a log slope of 1.12e-12 against Im ω₃₂₂ = 1.10e-12.
  - Steps cycle with a period of about 1.04e8 (about 6 rhs evaluations). Im ω₂₁₁ swings between +1.5e-14 and
    −2.6e-13 around zero, which looks like an explicit method at its stability limit.
- The restart state is the last row of `results/spindown_a0.9_alpha0.3_fa0.466873_solution.csv` (t = 1.5e9).
  - That CSV predates the 16-digit format: A has 9 digits, b and j have 16.
  - The restart path is `scripts/evolution_spindown.jl` with 4 arguments. It rebuilds the radial solution by 30-step
    amplitude continuation, because a cold solve can fail in the nonlinear regime.
  - Restart correctness was checked earlier in `docs/session-2026-11-10.md`.
- Context for the numbers: `docs/session-2026-11-10.md`. Search for "61102851" and "restart".

## Steps

Work in a scratch directory outside the repo, e.g. `$PSCRATCH/stiff-integrator-test/`. Do not modify `Project.toml`,
`Manifest.toml`, `src/`, or tracked files: this is a shared workspace.

1. **Environment.** Make a scratch Julia project:
   - `Pkg.develop(path="<worktree>")` for AxionCloud.
   - `Pkg.add(["OrdinaryDiffEqRosenbrock", "ADTypes"])`. Both are already in the depot; try `Pkg.offline(true)` if the
     network is slow. Let Pkg resolve against AxionCloud's `OrdinaryDiffEqCore` version.

   On Perlmutter: `module load julia/1.11.7`; on Lawrencium: `~/.juliaup/bin/julia +1.11.7`.

2. **Jacobian eigenvalues at the plateau, the decisive check, before any integration.**
   1. Rebuild the restart state as the driver does: copy the restart block of `scripts/evolution_spindown.jl`.
   2. Build an `AxionCloud.EvolutionCache` the way `evolve_with_spindown` does.
   3. Compute a 4×4 finite-difference Jacobian of `AxionCloud._spindown_rhs!` at that state. Warm-start every
      perturbed solve from the unperturbed solution.
   4. Print the eigenvalues.
   - Expected if the hypothesis holds: one or two eigenvalues with real part around −1e-8 to −1e-7. Tsit5's
     stability limit |hλ| ≈ 3 at h ≈ 1e8 fits there. The rest should be of order 1e-12 or smaller.
   - If no eigenvalue is that fast, the stiffness explanation is wrong. **Stop and report.** Something else limits
     Tsit5, for example solve noise or a real oscillation.
   - **Noise floor first.** Evaluate the rhs twice at the same state from two different warm starts: the
     continuation solution, and a solve continued from a slightly different amplitude. The difference between the
     two is the noise from `newton_tol = 1e-10`. Choose each finite-difference step so the rate change it produces
     is at least about 100× that noise. Report the noise, the steps, and the Jacobian. If no step gives a clean
     Jacobian, try `newton_tol = 1e-12` and say so.

3. **Stiff integrator in a scratch copy of `evolve_with_spindown`.**
   - Define the copy in your test script, calling AxionCloud internals by qualified name: `AxionCloud.EvolutionCache`,
     `AxionCloud._spindown_rhs!`, `AxionCloud._reject_failed_step`, `AxionCloud.current_problem`, and so on.
   - Use `Rosenbrock23(autodiff=AutoFiniteDiff())` first. It is a W-method, so it tolerates an approximate Jacobian.
     Try `Rodas5P` as well if Rosenbrock23 works.
   - Do not use ForwardDiff through `solve_radial`: the cache mutates and the warm start is stateful.
   - The system is autonomous, so pass `tgrad` returning zeros through an `ODEFunction` to save one rhs per Jacobian.
   - Things that must change from the Tsit5 version:
     - **Warm-start callback.** `_record_accepted_solution!` errors unless the last rhs evaluation was at the accepted
       state. That holds for Tsit5 (first-same-as-last) but not for a Rosenbrock method. In your copy, have the
       callback solve once at the accepted state, warm-started from `cache.last_solution`, and store that as
       `cache.accepted_solution`. This costs one extra solve per step.
     - **Diagnostics.** `cache.diagnostics` and the verbose log also record the Jacobian's perturbed evaluations. Do
       not use them as the trajectory. Use `saveat` and the ODE solution, or record diagnostics from the
       accepted-step callback.
     - **Failed solves.** A failed solve during a Jacobian evaluation sets `cache.solve_failed`, and
       `_reject_failed_step` then rejects the step. Check that this is what happens, rather than a step that silently
       uses zero rates.

4. **Short run.** Restart from t = 1.5e9 to 1.15e10 with `reltol = abstol = 1e-7`, as Tsit5 used.
   - Log every accepted step: t, step size, rhs count so far, A, b, j, Im ω₂₁₁.
   - Watch the first few steps. If the step size is not growing past about 1e9 within roughly 20 steps, stop and
     report. Do not tune for hours.
   - This is short enough for a login node only if steps are long. Otherwise use Slurm:
     `sbatch --account=m3166 --qos=shared --constraint=cpu --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=16G --time=04:00:00 --wrap="module load julia/1.11.7; export OPENBLAS_NUM_THREADS=4 JULIA_NUM_THREADS=1; cd <dir>; julia --project=. <script>"`.
     Write output to a log file so progress can be followed.

5. **Compare with Tsit5** at matched times, over the same span, using the 61102851 log.
   - Pass: steps ≫ 1e8 M₀ (say ≥ 1e9 on the plateau), and fewer total rhs evaluations than Tsit5's 564.
   - Pass: A₂₁₁, b, and a agree with the log to ~1e-6 relative.
   - Pass: ln A₃₂₂ grows at 1.10e-12 ± 1%.
   - Pass: Im ω₂₁₁ stays near zero. A smaller swing than Tsit5's is expected and fine.
   - If the short run passes, propose, but do not submit, a run from 1.5e9 to 2.1e11, so the full 61102851 record can
     be compared, and an estimate of the cost to reach 1e13.

## What to report back

Keep the report short and give numbers:
- The noise floor, the finite-difference steps, and the Jacobian eigenvalues.
- The step sizes and rhs count against Tsit5's.
- The comparison table.
- Pass or fail on each criterion.
- The scratch paths of the script and the log.

Describe the library change that would be needed, but do not make it. The user types library code, so give the
change as a diff for them to apply. The change would be:
- Let `evolve_with_spindown` take an `algorithm` keyword instead of hard-coding Tsit5.
- Generalize the warm-start callback.
- Add the stiff-solver package to `Project.toml`.

State what the test establishes and what it does not. A 1e10-M₀ window on one plateau says nothing about the
transition when |322⟩ takes over or about α = 0.187. Label the result preliminary.

## Constraints

- No `git push`. No commits unless the user asks. Check `git status` before touching anything in the worktree.
- Do not cancel or modify running Slurm jobs. Job 61154173 belongs to the user.
- Use long, descriptive variable names in new code (`CLAUDE.md`, "Naming style").
- Report f_a in GeV (f_a/M_Pl = 0.466873 ↔ 5.7e18 GeV, with the non-reduced M_Pl = 1.2209e19 GeV).
- Define any notation you introduce.
