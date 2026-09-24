"""Summarize the frozen-problem experiments without modifying production data."""
import csv
from pathlib import Path

root = Path(__file__).resolve().parents[3]
directory = root / 'runs/adiabatic/radialsolvertest'
rows = []
for path in sorted(directory.glob('*.csv')):
    with path.open() as stream:
        rows.extend(csv.DictReader(stream))
lines = ['# Radial solver comparison', '',
         'All cases use the same frozen cosine seed, cubic prediction, and A0 increment 0.1.',
         'Setup includes profile construction and initial residual/Jacobian compilation.',
         'Solve time excludes that warmup and the independent production replay.',
         'These are single-target timings, not broad performance claims.', '',
         '| Method | Status | Setup s | Solve s | Corrections | Native residual | Shooting replay |',
         '|---|---|---:|---:|---:|---:|---:|']
for row in rows:
    lines.append(f"| {row['method']} | {row['status']} | {float(row['setup_seconds']):.1f} | "
                 f"{float(row['solve_seconds']):.1f} | {row['iterations']} | "
                 f"{float(row['residual']):.3g} | {float(row['replay_residual']):.3g} |")
lines.extend(['', f'{len(rows)} case summaries available (nine original cases plus any follow-ups).', '',
              'Collocation residuals measure radial defects; shooting residuals measure matching.',
              'A small native residual does not establish discretization accuracy.', ''])
reference = next((row for row in rows if row['method'] == 'shooting_ad' and
                  row['status'] == 'converged'), None)
if reference:
    lines.extend(['## Observable agreement with AD shooting', '',
                  '| Method | Relative energy difference | Relative flux difference | Relative Im(omega) difference |',
                  '|---|---:|---:|---:|'])
    for row in rows:
        if row['status'] != 'converged':
            continue
        differences = [abs(float(row[key])-float(reference[key])) / max(abs(float(reference[key])), 1e-300)
                       for key in ('native_energy', 'total_flux', 'omega_imag')]
        lines.append(f"| {row['method']} | " + ' | '.join(f'{value:.3g}' for value in differences) + ' |')
fine_reference = next((row for row in rows if row['method'] == 'shooting_ad_rk0025' and
                       row['status'] == 'converged'), None)
if fine_reference:
    lines.extend(['', '## Collocation versus refined shooting (RK4 spacing 0.025)', '',
                  '| Method | Relative native energy difference | Relative flux difference | Relative Im(omega) difference |',
                  '|---|---:|---:|---:|'])
    for row in rows:
        if not row['method'].startswith('collocation') or row['status'] != 'converged':
            continue
        differences = [abs(float(row[key])-float(fine_reference[key])) / abs(float(fine_reference[key]))
                       for key in ('native_energy', 'total_flux', 'omega_imag')]
        lines.append(f"| {row['method']} | " + ' | '.join(f'{value:.3g}' for value in differences) + ' |')
lines.extend(['', 'Deuflhard denotes a simplified-correction contraction controller inspired by NLEQ-ERR;',
              'it is not a full reproduction of that algorithm. Default native residual tolerance is',
              '1e-9 with a 20-correction budget. Suffixed cases override settings as listed below.', '',
              '| Method | Native tolerance | RK4 spacing | Spatial resolution |',
              '|---|---:|---:|---:|'])
for row in rows:
    lines.append(f"| {row['method']} | {row.get('tolerance','1e-9')} | "
                 f"{row.get('rk_step','0.1')} | {row.get('resolution','see method name')} |")
collocation = sorted((row for row in rows if row['method'].startswith('collocation') and
                     row['status'] == 'converged'),
                    key=lambda row: float(row.get('resolution',
                        {'collocation_050': '0.5', 'collocation_025': '0.25', 'collocation_0125': '0.125'}.get(row['method'], 'nan'))), reverse=True)
if len(collocation) >= 3:
    lines.extend(['', '## Collocation refinement', '',
                  'Ratios compare successive differences on meshes whose spacing halves.',
                  'Fourth-order asymptotic convergence predicts a ratio near 16.', ''])
    for index in range(len(collocation)-2):
        coarse, medium, fine = collocation[index:index+3]
        for key in ('native_energy', 'total_flux', 'omega_real'):
            first_difference = float(coarse[key])-float(medium[key])
            second_difference = float(medium[key])-float(fine[key])
            ratio = first_difference/second_difference if second_difference else float('nan')
            lines.append(f"- {coarse['method']} → {medium['method']} → {fine['method']}: {key} ratio {ratio:.5g}.")
lines.append('')
report = directory / 'SUMMARY.md'
report.write_text('\n'.join(lines))
print('\n'.join(lines))
