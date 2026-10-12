"""Collect all 100 cells, preserving excluded, screened, halted, and failed states."""
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from campaign import CAMPAIGN, ROOT, load_manifest


def main():
    manifest = load_manifest()
    records = []
    values = np.full((10, 10), np.nan)
    markers = []
    for case in manifest['cases']:
        path = CAMPAIGN / 'results' / case['case_id'] / 'summary.json'
        result = json.loads(path.read_text()) if path.exists() else {}
        record = {**case, 'status': result.get('status', 'missing' if case['included'] else 'excluded_alpha')}
        for name in ('final_time_years', 'max_bound_cloud_envelope_over_fa',
                     'max_primary_peak_over_fa', 'max_secondary_peak_over_fa',
                     'max_cubic_over_binding_scale', 'final_spin', 'error'):
            record[name] = result.get(name, '')
        records.append(record)
        if result.get('finished'):
            values[case['axion_index'], case['mass_index']] = result['max_bound_cloud_envelope_over_fa']
        markers.append((case['mass_index'], case['axion_index'], record['status']))
    with (CAMPAIGN / 'results/summary.csv').open('w') as output:
        writer = csv.DictWriter(output, fieldnames=records[0])
        writer.writeheader()
        writer.writerows(records)
    figure, axis = plt.subplots(figsize=(8, 6.5), constrained_layout=True)
    valid = values[np.isfinite(values)]
    image = axis.imshow(values, origin='lower', cmap='viridis',
        norm=LogNorm(vmin=max(float(np.min(valid)), 1e-16) if len(valid) else 1e-12,
                     vmax=max(float(np.max(valid)), 1e-11) if len(valid) else 1.))
    figure.colorbar(image, ax=axis, label=r'Maximum bound-cloud envelope estimate $|\Phi|/f_a$')
    used = set()
    for x, y, status in markers:
        if status == 'excluded_alpha':
            marker, color, label = 'x', '.5', r'Excluded: $\alpha_0>0.45$'
        elif status == 'age_limit_linear_screen':
            marker, color, label = 's', 'white', 'Small-alpha linear age screen'
        elif status == 'age_limit':
            continue
        elif status in ('failed', 'missing', 'running'):
            marker, color, label = 'x', 'red', 'Failed / unfinished'
        else:
            marker, color, label = '^', 'white', 'Stopped at model / track limit'
        axis.scatter(x, y, marker=marker, s=35, facecolors='none' if marker != 'x' else color,
                     edgecolors=color, label=label if label not in used else None)
        used.add(label)
    axis.set_xticks(range(10), [f'{mass:.1f}' for mass in np.geomspace(2, 100, 10)], rotation=45)
    axis.set_yticks(range(10), [f'{mass:.2g}' for mass in np.geomspace(1e-14, 1e-11, 10)])
    axis.set_xlabel(r'Initial black-hole mass [$M_\odot$]')
    axis.set_ylabel('Axion rest energy [eV]')
    axis.set_title('QCD two-mode survey: up to 10 million years\nPreliminary; halted cells show only the simulated interval')
    axis.legend(loc='upper left', fontsize=8, framealpha=.9)
    for suffix in ('png', 'pdf'):
        figure.savefig(ROOT / f'plots/qcd_grid_20261012.{suffix}', dpi=180)
    counts = {status: sum(row['status'] == status for row in records)
              for status in sorted({row['status'] for row in records})}
    print(json.dumps(counts, indent=2))


if __name__ == '__main__':
    main()
