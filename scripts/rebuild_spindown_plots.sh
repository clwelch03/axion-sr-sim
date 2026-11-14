#!/bin/bash
# Rebuild every spin-down figure in plots/ from the saved results/ (diagnostics CSVs and Slurm logs).
# One line per figure below: the plotting command and the output file name.
#
# Run from the package root:
#   bash scripts/rebuild_spindown_plots.sh [output directory, default plots] [parallel jobs, default 4] [name filter]
# The name filter is a grep pattern; only figures whose command line matches are rebuilt (e.g. overlay).
#
# The chiral runs and the cosine runs paired with them (the last twelve lines) were made with the sped-up code
# (branch speedup/integration-next at 057a167); see the first line of each Slurm log in results/.
#
# The twolevel overlays use scripts/leaver_spindown.py with AXION_QUARTIC_COUPLING=1 (the cosine potential's own quartic
# term) and AXION_ODE_METHOD=LSODA. Two-level models with the cosine rates rescaled to the QCD potential's quartic term
# were removed on 2026-11-13.
#
# The overlay labels are the short run names (cosine, chiral, free, ...); scripts/plot_helpers.jl turns them into the text
# shown on the figure.
output_directory=${1:-plots}
parallel_jobs=${2:-4}
name_filter=${3:-.}
mkdir -p "$output_directory"

single="julia --project=scripts scripts/plot_spindown.jl"
overlay="julia --project=scripts scripts/plot_spindown_overlay.jl"
cosine_3e3="spindown_a0.9_alpha0.3_fa3e-3_t1e10 cosine spindown_free_a0.9_alpha0.3_fa3e-3 free"
hydrogenic_3e3="hydrogenic_spindown_a0.9_alpha0.3_fa3e-3 hydrogenic"
leaver_3e3="leaver_spindown_a0.9_alpha0.3_fa3e-3"

while read -r command_line; do
    while [ "$(jobs -rp | wc -l)" -ge "$parallel_jobs" ]; do sleep 1; done
    sh -c "$command_line" &
done < <(grep -e "$name_filter" <<EOF
$single spindown_a0.9_alpha0.3_fa3e-3_t1e10 $output_directory/spindown_a0.9_alpha0.3_fa3e-3_t1e10.png
$single spindown_a0.9_alpha0.3_fa0.0466873_t1e11_restart $output_directory/spindown_a0.9_alpha0.3_fa0.0466873_t1e11_restart.png
$single --compact spindown_a0.9_alpha0.3_fa0.0466873_t1e11_restart $output_directory/spindown_a0.9_alpha0.3_fa0.0466873_t1e11_restart_compact.png
$single spindown_a0.9_alpha0.3_fa0.466873 $output_directory/spindown_a0.9_alpha0.3_fa0.466873.png
$single results/spindown_a0.9_alpha0.3_fa0.466873_t1e13-61102851.out $output_directory/spindown_a0.9_alpha0.3_fa0.466873_t1e13_partial.png
$single spindown_a0.99_alpha0.3_fa3e-3 $output_directory/spindown_a0.99_alpha0.3_fa3e-3.png
$single --compact results/spindown_a0.99_alpha0.3_fa1e-2-61102849.out $output_directory/spindown_a0.99_alpha0.3_fa1e-2_partial_compact.png
$single spindown_free_a0.9_alpha0.3_fa3e-3 $output_directory/spindown_free_a0.9_alpha0.3_fa3e-3.png
$single --compact spindown_free_a0.9_alpha0.3_fa3e-3 $output_directory/spindown_free_a0.9_alpha0.3_fa3e-3_compact.png
$overlay $cosine_3e3 $output_directory/spindown_overlay_cosine_free_a0.9_alpha0.3_fa3e-3.png
$overlay --compact $cosine_3e3 $output_directory/spindown_overlay_cosine_free_a0.9_alpha0.3_fa3e-3_compact.png
$overlay $cosine_3e3 $hydrogenic_3e3 $output_directory/spindown_overlay_cosine_free_hydrogenic_a0.9_alpha0.3_fa3e-3.png
$overlay --markers $cosine_3e3 $hydrogenic_3e3 $leaver_3e3 Leaver $output_directory/spindown_overlay_cosine_free_hydrogenic_Leaver_a0.9_alpha0.3_fa3e-3.png
$overlay --markers --compact $cosine_3e3 $leaver_3e3 'free (Leaver)' $output_directory/spindown_overlay_cosine_free_leaver_a0.9_alpha0.3_fa3e-3_compact.png
$overlay --markers --compact spindown_a0.9_alpha0.3_fa0.0466873_t1e11_restart cosine spindown_free_a0.9_alpha0.3_fa0.0466873_t1e11 free leaver_spindown_a0.9_alpha0.3_fa0.0466873 'free (Leaver)' $output_directory/spindown_overlay_cosine_free_leaver_a0.9_alpha0.3_fa0.0466873_compact.png
$overlay --markers --compact spindown_a0.9_alpha0.3_fa0.466873 cosine spindown_free_a0.9_alpha0.3_fa0.466873_t1.5e9 free leaver_spindown_a0.9_alpha0.3_fa0.466873 'free (Leaver)' $output_directory/spindown_overlay_cosine_free_leaver_a0.9_alpha0.3_fa0.466873_compact.png
$single --compact spindown_chiral_a0.9_alpha0.3_fa3e-3_t1e10 $output_directory/spindown_chiral_a0.9_alpha0.3_fa3e-3_t1e10_compact.png
$single --compact spindown_chiral_a0.9_alpha0.3_fa0.0466873_t1e11 $output_directory/spindown_chiral_a0.9_alpha0.3_fa0.0466873_t1e11_compact.png
$single --compact spindown_chiral_a0.9_alpha0.3_fa0.466873 $output_directory/spindown_chiral_a0.9_alpha0.3_fa0.466873_compact.png
$single --compact spindown_chiral_a0.99_alpha0.3_fa3e-3 $output_directory/spindown_chiral_a0.99_alpha0.3_fa3e-3_compact.png
$single --compact spindown_chiral_a0.99_alpha0.3_fa1e-2 $output_directory/spindown_chiral_a0.99_alpha0.3_fa1e-2_compact.png
$single --compact spindown_a0.9_alpha0.3_fa0.0466873_t1e11 $output_directory/spindown_a0.9_alpha0.3_fa0.0466873_t1e11_compact.png
$single --compact spindown_a0.99_alpha0.3_fa1e-2 $output_directory/spindown_a0.99_alpha0.3_fa1e-2_compact.png
$overlay --compact spindown_a0.9_alpha0.3_fa3e-3_t1e10 cosine spindown_chiral_a0.9_alpha0.3_fa3e-3_t1e10 chiral $output_directory/spindown_overlay_cosine_chiral_a0.9_alpha0.3_fa3e-3_t1e10_compact.png
$overlay --compact spindown_a0.9_alpha0.3_fa0.0466873_t1e11 cosine spindown_chiral_a0.9_alpha0.3_fa0.0466873_t1e11 chiral $output_directory/spindown_overlay_cosine_chiral_a0.9_alpha0.3_fa0.0466873_t1e11_compact.png
$overlay --compact spindown_a0.9_alpha0.3_fa0.466873 cosine spindown_chiral_a0.9_alpha0.3_fa0.466873 chiral $output_directory/spindown_overlay_cosine_chiral_a0.9_alpha0.3_fa0.466873_compact.png
$overlay --compact spindown_a0.99_alpha0.3_fa3e-3 cosine spindown_chiral_a0.99_alpha0.3_fa3e-3 chiral $output_directory/spindown_overlay_cosine_chiral_a0.99_alpha0.3_fa3e-3_compact.png
$overlay --compact spindown_a0.99_alpha0.3_fa1e-2 cosine spindown_chiral_a0.99_alpha0.3_fa1e-2 chiral $output_directory/spindown_overlay_cosine_chiral_a0.99_alpha0.3_fa1e-2_compact.png
$overlay --compact spindown_a0.9_alpha0.3_fa3e-3_t1e10 cosine leaver_spindown_a0.9_alpha0.3_fa3e-3_quartic1_lsoda 'two-level model, cosine coupling' $output_directory/spindown_overlay_cosine_twolevel_a0.9_alpha0.3_fa3e-3_t1e10_compact.png
$overlay --compact spindown_a0.9_alpha0.3_fa0.0466873_t1e11 cosine leaver_spindown_a0.9_alpha0.3_fa0.0466873_quartic1_lsoda 'two-level model, cosine coupling' $output_directory/spindown_overlay_cosine_twolevel_a0.9_alpha0.3_fa0.0466873_t1e11_compact.png
$overlay --compact spindown_a0.9_alpha0.3_fa0.466873 cosine leaver_spindown_a0.9_alpha0.3_fa0.466873_quartic1_lsoda 'two-level model, cosine coupling' $output_directory/spindown_overlay_cosine_twolevel_a0.9_alpha0.3_fa0.466873_compact.png
$overlay --compact spindown_a0.99_alpha0.3_fa3e-3 cosine leaver_spindown_a0.99_alpha0.3_fa3e-3_quartic1_lsoda 'two-level model, cosine coupling' $output_directory/spindown_overlay_cosine_twolevel_a0.99_alpha0.3_fa3e-3_compact.png
$overlay --compact spindown_a0.99_alpha0.3_fa1e-2 cosine leaver_spindown_a0.99_alpha0.3_fa1e-2_quartic1_lsoda 'two-level model, cosine coupling' $output_directory/spindown_overlay_cosine_twolevel_a0.99_alpha0.3_fa1e-2_compact.png
EOF
)
wait
