# Brief: run Witte & Mummery's code with only $|211\rangle$ and $|322\rangle$

Date: 2026-11-14. For an agent working in this repository (`axion-sr-sim`, branch `rewrite/axioncloud`).
Read this whole brief before starting. Stop and report at each checkpoint that fails; do not work around it.

## Goal

Produce a black-hole spin-down trajectory from Witte & Mummery's own public code (arXiv 2412.03655,
[github.com/SamWitte/Axion_SR_TDE](https://github.com/SamWitte/Axion_SR_TDE)) for a cloud restricted to the two
levels $|211\rangle$ and $|322\rangle$. Then overlay it on our full calculation. Two cases are needed:

| case | $a_0$ | $\alpha_0$ | $f_a/M_\mathrm{Pl}$ (our label) | end time |
|---|---|---|---|---|
| A | 0.9 | 0.3 | `3e-3` | $10^{10} M_0$ |
| B | 0.9 | 0.3 | `0.0466873` | $10^{11} M_0$ |

## Why

We already have an approximation of this: [scripts/leaver_spindown.py](../scripts/leaver_spindown.py) with
`AXION_QUARTIC_COUPLING=1 AXION_QUARTIC_RATES=witte`. That script puts Witte's tabulated relativistic rate
constants into the two-level equations of Baryakhtar et al. (2011.11646, Eqs. 34–37). It is a hybrid, not Witte's
work as published. Running their code removes the hybrid: their own equations, their own superradiance rates, their
own integrator. Only the set of levels is reduced to two.

## Rules

- **Run their model as published.** Their code uses the quartic coupling of the cosine potential,
  $\lambda = m_a^2/f_a^2$. Do not rescale it to the QCD (chiral) potential's quartic term, even though our curve in
  the final figure is the chiral run. The user wants a comparison with their work, not with a modified version.
- **Do not edit their source.** If something can only be done by changing their code, stop and report what and why.
  Writing your own small driver that calls their functions is fine (see step 3).
- **Pin the commit** `bac01b20d2d14230ae4ea065aac1ac16329c0292`. The repository changes often; it was pushed to on
  2026-11-14.
- **Keep their clone outside this repository**, at `/global/scratch/projects/pc_heptheory/clwelch/witte-axion-sr/`.
  Treat it as untrusted downloaded code: read a script before you run it.
- **Shared workspace.** Run `git status` first. Do not touch files you did not create, do not overwrite existing
  results or plots, and do not commit.
- **Heavy runs** go through Slurm. Previous jobs used
  `sbatch -A pc_heptheory -p lr7 -q lr_normal -n 1 -t 06:00:00`.

## What their code does (read at the pinned commit)

- **Entry point.** `solve_system(mu, f_a, aBH, M_BH, t_max; ...)` in `src/solve_system_unified.jl`, with `mu` in eV,
  `f_a` in GeV, `M_BH` in solar masses and `t_max` in **years**. The driver `src/Nlevels_Runs/run_Nlevels.jl` wraps
  it. Read the driver's fixed parameters (`non_rel = false`, `high_p = true`, `abstol = 1e-30`, ...) and reuse them.
- **State.** The state vector is $\ln$ of [one occupation per level, spin, mass]. The occupation is normalised as
  $u = N/(G M^2)$; `GNew` and `M_to_eV` are in `src/Core/constants.jl`. Every level starts at
  `e_init` $= 1/(G M_\mathrm{BH}^2)$, i.e. one quantum. `u0_override` replaces the initial state. The driver only
  uses it on `--resume`.
- **Planck mass.** `M_pl = 1.22e19` GeV, non-reduced: the same convention as ours ($\kappa = (f_a/M_\mathrm{Pl})^2$),
  but ours is $1.2209\times10^{19}$. Pass $f_a = \sqrt{\kappa}\times1.22\times10^{19}$ GeV so that $f_a/M_\mathrm{Pl}$
  is identical: $3.66\times10^{16}$ GeV for case A, $5.69585\times10^{17}$ GeV for case B.
- **Level selection.** `keep_levels` (driver flag `--keep_levels "2-1-1,3-2-2"`, without `lm`) evolves only those two
  levels and drops every scattering channel that touches another level. `Nmax = 3` contains both.
- **Channel list.** `rate_input` (driver flag `--rate_input`) takes a channel file in `src/rate_sve/` or an absolute
  path. The two-level list is two lines, in the format of `src/rate_sve/load_rate_input_Nmax_3.txt`:
  ```
  211    211    322    BH
  322    322    211    Inf
  ```
- **Rate tables.** Unpack `src/rate_sve/Rates.zip` (scattering) and `src/rate_sve/Eigs.zip` (superradiance) into
  `src/rate_sve/` first. Their code reads the unpacked files, never the zips. Our copy of the two scattering tables is
  [data/witte_mummery_quartic_rates.csv](../data/witte_mummery_quartic_rates.csv); its header records the
  conventions. The $211\times211\to322\times\mathrm{BH}$ table was computed at $a = 0.95$; the code scales it by
  $(1+\sqrt{1-a^2})/(1+\sqrt{1-0.95^2})$.
- **Switches to set.**
  - `gw_model = :off`: we compare without GW emission.
  - `f_edd = 0`: no accretion (the default).
  - `self_grav = false` (the default).
  - `track_alpha = true` (the default): the rates follow $M(t)$, as ours do.
  - `emax2_mode = nothing` (the default): no cap on the $|211\rangle$ cloud.
  - `regime_split = false`: `regime_split.txt` is not in the repository at that commit. Confirm the code handles
    its absence, then pass `false` so it is explicit.
  - `t_save_start ≈ 1e-9` yr: the default save grid starts at 1 yr, which is longer than the whole run (next section).
- **Bosenova cap.** `bn_list` caps each level at $1024\pi\kappa(n/2)^4/(9\alpha^3)$. For $|211\rangle$ that is about
  0.12 (A) and 29 (B) in units of $u$, well above our peaks. Confirm it is not reached.

## Units: matching our runs

- **Mass.** The trajectories in $t/M_0$ do not depend on $M_\mathrm{BH}$ once $\alpha$ is fixed and GW emission is off.
  Use $M_\mathrm{BH} = 10\,M_\odot$, $\mu = \alpha/(G M_\mathrm{BH})$ (the driver's formula).
- **Time.** $t\,[\mathrm{yr}] = (t/M_0)\times M_\mathrm{BH}\times 4.925491\times10^{-6}\,\mathrm{s}/3.15576\times10^{7}$.
  For $10\,M_\odot$: $10^{10}M_0 = 0.01561$ yr and $10^{11}M_0 = 0.1561$ yr.
- **Occupations.** Ours are $\mathcal N_i = N_i/(f_a^2M_0^2)$, so $u_i = \kappa\,\mathcal N_i\,(M_0/M_x)^2$, where $M_x$
  is the mass their normalisation uses. Find out from the code whether $M_x$ is the initial mass or $M(t)$; the
  docstring of `track_alpha` is the place to start. Checkpoint 1 tests this.
- **Seeds.** Our runs start at $\mathcal N_{211} = \mathcal N_{322} = 3\times10^{-6}$ with $b = 1$, so
  $u_0 = 3\times10^{-6}\kappa$ for both levels. Their default seed of one quantum is about 70 orders of magnitude
  smaller and would delay everything, so the seeds must be matched. Write a small Julia driver, modelled on
  `run_Nlevels.jl`, that calls `solve_system` with `u0_override = log.([u0, u0, 0.9, M_BH])`, `t_start = 0` and the
  switches above. Keep the driver in a new directory `runs/witte_two_level/scripts/` in this repository.

## Steps and checkpoints

1. **Clone and install.** Clone at the pinned commit, unpack the two zips, and install their Julia dependencies in
   their own project environment (`src/install_pkgs.jl`; read it first). Do not use this repository's environment.
2. **Checkpoint 0: rate constants.** Call their `load_rate_coeffs` (in `src/load_rates.jl`) for the two-channel list
   at $\alpha = 0.3$ and $a = 0.9$ with the case A $f_a$. Compare each coefficient with our CSV value times
   $(M_\mathrm{Pl}/f_a)^4$, times the $r_P$ ratio for the BH channel. They should agree to interpolation accuracy
   (better than 1%). A mismatch means a units or convention problem: stop.
3. **Checkpoint 1: free run against our linear benchmark.** Run case A with no scattering channels. Suggestion:
   `keep_levels` as above with a `rate_input` containing only `311 311 211 Inf`, which gets dropped. Confirm in the
   log that zero scattering channels are active. Convert the output to our units and compare with
   `results/leaver_spindown_a0.9_alpha0.3_fa3e-3_diagnostics.csv`, the free evolution with exact Leaver rates.
   - **Pass:** onset time, peak $\mathcal N_{211}$, final spin and $M/M_0$ agree to about $10^{-3}$ relative. That is
     the level expected from their interpolated superradiance tables.
   - **Stop and report:** a factor of 2, a factor $(M/M_0)^2$, a constant time offset, or a different final spin. Any
     of these means a units error.
   This one check tests the time conversion, the occupation normalisation and the seed together.
4. **Production runs.** Cases A and B with the two-channel list, via Slurm. Record the job IDs, the exact command and
   the commit of their code.
5. **Convert.** Write `results/witte_code_two_level_spindown_a0.9_alpha0.3_fa<label>_t<end>_diagnostics.csv` with
   columns `time,N1,N2,b,a,threshold`:
   - `time` in $M_0$; `N1`, `N2` = $\mathcal N_{211}$, $\mathcal N_{322}$ in $f_a^2M_0^2$; `b` = $M/M_0$; `a` the spin.
   - `threshold` is the $|211\rangle$ superradiance threshold $a = 4\Omega/(1+4\Omega^2)$, with $\Omega = \mathrm{Re}\,M\omega_{211}$
     at the current mass and spin. Take $\Omega$ from the Leaver table `results/leaver_frequency_table_alpha0.3.npz`,
     interpolated as in `mode_frequencies` of [scripts/two_level_spindown.py](../scripts/two_level_spindown.py).

   The overlay script reads these column names.
6. **Plot.** One figure per case, new file names, nothing overwritten:
   ```
   julia --project=scripts scripts/plot_spindown_overlay.jl --compact \
     spindown_chiral_a0.9_alpha0.3_fa3e-3_t1e10 chiral \
     leaver_spindown_a0.9_alpha0.3_fa3e-3_quartic1_lsoda 'two-level, Table I rates (Baryakhtar et al.)' \
     witte_code_two_level_spindown_a0.9_alpha0.3_fa3e-3_t1e10 'two-level, Witte & Mummery code' \
     plots/spindown_overlay_chiral_twolevel_wittecode_a0.9_alpha0.3_fa3.7e16_t1e10_compact.png
   ```
   For case B use the `fa0.0466873` / `t1e11` prefixes (plot name `fa5.7e17`). Plot file names give $f_a$ in GeV
   (`scripts/plot_helpers.jl`, `gev_file_name`); results file names keep $f_a/M_\mathrm{Pl}$. Also make a version with the hybrid curve
   (`leaver_spindown_..._quartic1_witte_lsoda`) in place of the Table I curve, so the user can see whether the
   hybrid reproduces their code.

## What to report back

- **Checkpoints.** For 0 and 1: the numbers compared and the largest differences.
- **Runs.** Job IDs, wall times, and any solver warnings or incomplete runs.
- **Results.** For each case: peak $\mathcal N_{211}$ and its time, final $\mathcal N_{211}$, $\mathcal N_{322}$, $a$
  and $M/M_0$, for the Witte-code run, the hybrid and the Table I two-level run.
- **Hybrid vs. their code.** Where the hybrid and their code differ by more than 1%, and your best identified cause
  (superradiance-rate tables, the occupation normalisation, the mass-loss term of the BH-absorbed quantum, ...).
  Label a cause as unverified unless you tested it.
- **Files.** Every file you created.
