# Temporary Perlmutter setup

Prepared November 3, 2026 for NERSC user `clwelch`, project `m3166`.
Use `$SCRATCH/axion-sr-sim` as the checkout unless another destination is chosen.
Push the source commits to GitHub before transferring local artifacts.

## What to preserve

Git excludes all CSV tables and JLS checkpoints, campaign results, Dolan NPZ
files, generated figures/logs, and campaign `work/` directories. The latter
contain frozen source, bundled inputs, provenance, and standalone scripts.
Preserve the entire `runs/` tree except Julia depots and bytecode caches.
Also preserve `plots/`, `docs/multiphase-code-map.md`, and
`docs/diagrams/multiphase/`; the last two are excluded by `.git/info/exclude`.
Existing `recovery/` material belongs in a separate private archive.

The ignored staging directory is `recovery/perlmutter-20261103/`. Its data
archive includes the campaign recipes as they existed when packaged, so extract
it into a checkout of the matching source commit before making Perlmutter edits.
Keep this directory private: its recovery archive contains private backups.
The archives exclude credentials, SSH keys, and Julia package/compiler caches.

The staging directory also holds a Git bundle of all local branches and the
exact Julia 1.13.0 Linux x86-64 runtime. These preserve the source and runtime
independently of remote download availability. Consult `transfer-status.json`
there for actual preparation, push, and remote-validation status; preparation
alone does not establish that the data reached NERSC.

## Source and data transfer

After the commits have been pushed, open a Perlmutter shell and prepare the
transfer destination. The verified Git bundle carries the pushed source, so
obtaining the checkout does not need a second GitHub sign-in on Perlmutter:

```bash
mkdir -p "$SCRATCH/axion-transfer-20261103"
```

From the Perlmutter shell, pull the prepared directory through Lawrencium's
documented data-transfer node. Authenticate to Lawrencium in that terminal:

```bash
rsync -av --partial \
  clwelch@lrc-xfer.lbl.gov:/global/scratch/projects/pc_heptheory/clwelch/axion-sr-sim/recovery/perlmutter-20261103/ \
  "$SCRATCH/axion-transfer-20261103/"
```

This avoids the direct SSH connection from the Lawrencium login node to
Perlmutter, which timed out during preparation. Globus between `lbnl#lrc`
and `NERSC DTN` is another supported route if SSH transfer is unavailable.

On Perlmutter, verify the archive checksums before extraction:

```bash
cd "$SCRATCH/axion-transfer-20261103"
sha256sum -c SHA256SUMS
git clone --branch main "$PWD/axion-source.bundle" "$SCRATCH/axion-sr-sim"
git -C "$SCRATCH/axion-sr-sim" remote set-url origin \
  https://github.com/clwelch03/axion-sr-sim.git
tar -xzf axion-local-data.tar.gz -C "$SCRATCH/axion-sr-sim"
```

The separate private-recovery archive can be retained without extracting it.
The Git bundle can be checked with `git bundle verify axion-source.bundle` from
an existing repository. Use a new checkout directory if the default already
contains other work. A normal GitHub clone at the recorded source commit is
also suitable if NERSC already has GitHub authentication configured.

Use new output names on Perlmutter. Preserve existing checkpoint prefixes and
their frozen source. Never let two jobs write the same prefix. To refresh the
transfer after additional Lawrencium work, first let writers finish, update the
source commits, and recreate the archives and checksums.

## Recreate Python and Julia

`environments/BHsuper-lawrencium.yml` exports the current BHsuper environment
without Lawrencium's installation prefix. It pins Python 3.11.15 and package
versions, including pip-installed packages. Python 3.11 or newer is required by
recent campaign scripts using `tomllib`.

```bash
cd "$SCRATCH/axion-sr-sim"
module load conda
conda env create --name BHsuper --file environments/BHsuper-lawrencium.yml
conda activate BHsuper
export AXION_PYTHON="$(command -v python)"
```

If BHsuper already exists on NERSC, use a new environment name for this export.
Select Julia **1.13.0**, matching both checked-in manifests and the runtime used
to create the current checkpoints. Load NERSC's Julia module for site
preferences, then select the preserved pinned runtime:

```bash
module load julia
mkdir -p "$HOME/.local/opt"
tar -xzf "$SCRATCH/axion-transfer-20261103/julia-1.13.0-linux-x86_64.tar.gz" \
  -C "$HOME/.local/opt"
export AXION_JULIA="$HOME/.local/opt/julia-1.13.0+0.x64.linux.gnu/bin/julia"
export AXION_JULIA_DEPOT="$SCRATCH/axion-julia-depot"
bash scripts/hpc/setup_perlmutter.sh
```

The setup script instantiates and precompiles both Julia projects, installs
the Python checkout, and checks Python imports. Use the same depot for setup
and batch runs. Keep the exported `AXION_*` variables in the submission shell.
Package source/compiler caches are rebuilt on NERSC. Serialized objects can
depend on package and source definitions, so retain corresponding frozen source.

## Rehearse on a compute node

Submit from the checkout root, after setup has created the logs directory:

```bash
sbatch --account=m3166 --time=00:20:00 \
  scripts/hpc/perlmutter_job.sh "$PWD" scripts/hpc/smoke_perlmutter.jl
```

Success prints `PERLMUTTER SMOKE PASSED`. The check loads a transferred QCD
primary checkpoint, compares analytic phase averages and derivatives with
direct tensor sums, and runs a small nonlinear matching/continuation case.
This verifies startup, checkpoint compatibility, and the numerical paths used
by those regressions; it is not spatial convergence of a production campaign.

`perlmutter_job.sh` requests one serial process in the CPU shared QOS. It uses
one Julia/BLAS thread and verifies frozen source manifests before execution.
It accepts a checkout or frozen snapshot, a driver path relative to that root,
and driver arguments. Pass absolute paths for driver inputs and outputs.
For array drivers, supply the original `--array` range to `sbatch`.

For example, rerun one previously passing QCD comparison under a new output name:

```bash
campaign="$PWD/runs/adiabatic/multiphase_two_cloud"
sbatch --account=m3166 --array=1 \
  scripts/hpc/perlmutter_job.sh "$campaign/work/qcd_211_comparisons_20261030_v2" \
  runs/adiabatic/multiphase_two_cloud/scripts/qcd_211_comparison.jl \
  "$campaign/work/qcd_211_comparisons_20261030_v2/inputs/cases.toml" \
  '$SLURM_ARRAY_TASK_ID' "$campaign/results/perlmutter_qcd_comparison_20261103"
```

Drivers need the numeric array task index. The generic wrapper substitutes a
literal `$SLURM_ARRAY_TASK_ID` argument with the current task's index. It runs
the Julia/Python driver directly, bypassing the Lawrencium-only launch scripts.
Keep frozen snapshots unchanged so their hash checks remain meaningful.

NERSC references: [shared CPU jobs](https://docs.nersc.gov/jobs/examples/#shared),
[Julia modules and pinned runtimes](https://docs.nersc.gov/development/languages/julia/),
and [SSH connections](https://docs.nersc.gov/connect/).
Lawrencium references: [data-transfer node](https://scienceit-docs.lbl.gov/hpc/data-transfer-node/)
and [Globus collection](https://scienceit-docs.lbl.gov/data/globus-instructions/).
