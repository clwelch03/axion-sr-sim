"""Plot all three checkpointed Figure 2 branches; run from any directory."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

script_directory = Path(__file__).resolve().parent
root = script_directory.parents[3]
run_directory = script_directory.parent / "results"
plt.rcParams.update({'font.family': 'serif', 'mathtext.fontset': 'stix',
                     'font.size': 12})
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
potentials = (
    ('quartic', 'black', r'$\frac{1}{2}\phi^2-\frac{1}{4!}\phi^4$'),
    ('sextic', 'red', r'$\frac{1}{2}\phi^2-\frac{1}{4!}\phi^4+\frac{1}{6!}\phi^6$'),
    ('cosine', 'blue', r'$1-\cos\phi$'),
)
for potential, color, label in potentials:
    data = np.atleast_1d(np.genfromtxt(
        run_directory / f'{potential}_branch.csv', delimiter=',', names=True))
    for suffix in ('refined', 'gap'):
        refinement_path = run_directory / f'{potential}_{suffix}.csv'
        if not refinement_path.exists():
            continue
        refinement = np.atleast_1d(np.genfromtxt(
            refinement_path, delimiter=',', names=True))
        # Keep production values for shared amplitudes; add independently solved points.
        additional = refinement[~np.isin(refinement['amplitude'], data['amplitude'])]
        data = np.concatenate((data, additional))
        data = data[np.argsort(data['amplitude'])]
    for column in ('amplitude', 'energy', 'total_flux', 'residual'):
        if not np.all(np.isfinite(data[column])):
            raise ValueError(f'{potential}: nonfinite {column}')
    if np.any(np.diff(data['amplitude']) <= 0):
        raise ValueError(f'{potential}: amplitudes must increase')
    # Connect saved solutions in amplitude order, including the energy turn.
    axes[0].plot(data['amplitude'], data['energy'], color=color, lw=1.5,
                 label=label)
    axes[1].plot(data['energy'], data['total_flux'], color=color, lw=1.5,
                 label=label)
    print(f'{potential}: {len(data)} points; A0_max={data["amplitude"][-1]:.6g}; '
          f'E_last={data["energy"][-1]:.9g}; '
          f'max_residual={np.max(data["residual"]):.3g}')
axes[0].set_xlabel(r'$A_0$')
axes[0].set_ylabel(r'$E/[M(F_a/M_{\rm pl})^2]$')
axes[1].axhline(0, color='0.5', lw=0.8)
axes[1].set_xlabel(r'$E/[M(F_a/M_{\rm pl})^2]$')
axes[1].set_ylabel(r'$F_{\rm tot}/F_a^2$')
axes[0].set_xlim(0, 13.3)
axes[0].set_ylim(0, 1850)
axes[1].set_xlim(0, 1820)
axes[1].set_ylim(-0.016, 0.002)
axes[0].legend(loc='lower right', frameon=False)
axes[1].legend(loc='lower left', frameon=False)
for axis in axes:
    axis.minorticks_on()
    axis.tick_params(which='both', direction='in', top=True, right=True)
    axis.set_title(r'$\mu M=0.42,\ a/M=0.99$')
fig.suptitle('Preliminary Figure 2 reproduction — saved branches', fontsize=11)
output = root / 'plots/preliminary_fig2_20260922'
output.parent.mkdir(parents=True, exist_ok=True)
for extension in ('png', 'pdf'):
    fig.savefig(output.with_suffix('.' + extension), dpi=180)
print(output.with_suffix('.png'))
