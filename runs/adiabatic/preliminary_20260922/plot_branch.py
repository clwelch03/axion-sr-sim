"""Plot the checkpointed preliminary cosine branch; run from any directory."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

run_directory = Path(__file__).resolve().parent
root = run_directory.parents[2]
data = np.atleast_1d(np.genfromtxt(run_directory / 'cosine_branch.csv', delimiter=',', names=True))
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
axes[0].plot(data['amplitude'], data['energy'], 'o-', color='blue', ms=4)
axes[0].set_xlabel(r'$A_0$')
axes[0].set_ylabel(r'$E/[M(F_a/M_{\rm pl})^2]$')
axes[1].plot(data['energy'], data['total_flux'], 'o-', color='blue', ms=4,
             label=r'$1-\cos\phi$')
axes[1].axhline(0, color='0.5', lw=0.8)
axes[1].set_xlabel(r'$E/[M(F_a/M_{\rm pl})^2]$')
axes[1].set_ylabel(r'$F_{\rm tot}/F_a^2$')
axes[1].legend()
for axis in axes:
    axis.grid(alpha=0.2)
fig.suptitle(r'Preliminary cosine branch: $\mu M=0.42$, $a/M=0.99$'
             '\nSix modes; finite boundaries; not a converged Fig. 2 reproduction', fontsize=11)
output = root / 'plots/preliminary_cosine_branch_20260922'
output.parent.mkdir(parents=True, exist_ok=True)
for extension in ('png', 'pdf'):
    fig.savefig(output.with_suffix('.' + extension), dpi=180)
print(f'{len(data)} points; A0_max={data["amplitude"][-1]:.6g}; '
      f'E_last={data["energy"][-1]:.9g}; peak_field={data["peak_field"][-1]:.6g}')
print(output.with_suffix('.png'))
