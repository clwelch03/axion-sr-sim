# Nonlinear tower progress — 2026-09-15

We have extended the validated linear Dolan solver with the machinery to evolve a real axion field on fixed Kerr, obeying $\Box_g\Phi-\mu^2\sin\Phi=0$, where $\Phi=\phi/f_a$. The new implementation combines positive odd-$m$ towers, reconstructs the negative-$m$ contributions through the reality condition, and evaluates the nonlinear remainder $\Phi-\sin\Phi$ on an angular grid at each radial node. After multiplication by the required geometric factors, the projected source is added to the existing tower equations before the Thomas solve. A synchronized RK4 step recomputes this source from all towers at every intermediate stage. The existing linear operators, PML prescription, and boundary freezing are retained; the new source is evaluated accurately even at very small field amplitudes.

Incremental checks verified angular reconstruction/projection to roundoff, the analytic cubic source and generation of higher-$m$ components, and the source normalization against the existing mass terms (cancellation residual $3.4\times10^{-15}$ relative). Disabling nonlinearity reproduces the original linear steps exactly. Short nonlinear integrations show fourth-order timestep convergence, with a refinement ratio of 16.17, and recover linear evolution as the initial amplitude decreases: halving the initial peak $|\Phi|$ reduces the amplitude-normalized field discrepancy by approximately four, as expected from cubic forcing. These results validate the implementation and its small-amplitude limit on the tested grids; they do not yet establish accurate nonlinear cloud growth or saturation. Next are finite-amplitude angular and radial convergence checks, followed by physical energy, angular-momentum, and boundary-flux diagnostics. Long cloud evolutions and nonlinear growth-rate measurements remain to be done.

## Finite-amplitude angular sampling checkpoint

`tests/test_finite_amplitude_quadrature.py` now checks polar and azimuthal
quadrature refinement independently for a mixed six-channel field through
$\ell=m=5$, at sampled peak $|\Phi|=1$ and 3. At amplitude 3, the relative
source errors for 8, 16, 32 polar points (128 azimuthal points fixed) are
$1.70\times10^{-2}$, $1.61\times10^{-5}$, and $6.17\times10^{-14}$.
For 16, 32, 64 azimuthal points (64 polar points fixed), they are
$6.11\times10^{-3}$, $2.76\times10^{-6}$, and $7.24\times10^{-16}$.
The reference uses 64×128 points and is checked against 80×160.

A separate short evolution with initial peak $|\Phi|=3$, $\ell_{\max}=m_{\max}=3$,
$\Delta x=M$, and $\Delta t=0.025M$ runs to $t=0.5M$. Refining 6×12 to
12×24 angular points reduces the generated $m=3$ sector's relative error
from $1.03\times10^{-4}$ to roundoff against 32×64. These are sampling checks
at **fixed mode truncation**, not convergence in the number of retained modes
or evidence for accurate long-time evolution. Finite-amplitude mode-truncation
and radial convergence remain open.

## Runtime is now a priority

See [nonlinear-performance.md](nonlinear-performance.md) for explanations and
the subsequent temporary-array checkpoint: removing one redundant result copy
gave modest isolated gains with identical outputs; all 23 tests still pass.

The user reports hundreds of hours projected for $10^5M$ at $\Delta x=M/6$.
Before long evolutions, benchmark and improve the angular source; consider a
Julia kernel alongside a NumPy implementation using precomputed arrays and
reusable buffers. A full solver port has not been started.

`src/dolan/benchmark_nonlinear.py` provides a short, reproducible timing/profile probe:

```bash
module load miniconda3/22.11.1-gcc-11.4.0
OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  PYTHONPATH=src conda run -n BHsuper python -m dolan.benchmark_nonlinear --profile
```

Defaults: $a=0.99$, $\alpha=0.42$, $x/M\in[-200,511]$, 4267 radial nodes,
all six positive odd-$m$, odd-$\ell$ channels through 5, 32×64 angular points,
and $\Delta t=0.8\Delta x$. The seed is a Gaussian with peak $|\Phi|=1$.
These are assumed benchmark settings, not a recovered production configuration.

Local single-thread timings on 2026-09-15 (median of three warmed calls):

| Implementation | Nonlinear RK4 step | Projected 750,000 steps |
|---|---:|---:|
| Original source | 4.664 s | 972 hours |
| Replace `small_field**3` with `small_field * squared_field` | 1.712 s | 357 hours |

The original profile spent 4.09 of 4.56 seconds in `sine_remainder`; after
the change it spent 1.19 of 1.70 seconds there. Linear RK4 for the same towers
takes roughly 0.016–0.018 seconds. Thus source evaluation, rather than the
tower dictionaries or Thomas solve, is the measured first target. The source
already uses NumPy matrix products in radial blocks, but still constructs
temporary arrays and repeats static projection setup each RHS call.

These are fixed-state microbenchmarks on the current machine, excluding I/O
and setup; timings vary with hardware, field, quadrature, and thread settings.
Repeat on a compute node before allocating production time. All 23 tests pass
after the arithmetic change, including small-amplitude normalization and the
new finite-amplitude checks. Next performance checkpoint: compare an optimized
source kernel against this implementation for both accuracy and wall time,
before choosing whether to migrate the evolution machinery to Julia.

The subsequent polynomial-buffer checkpoint is complete: the small-field
expansion reuses one work array and preserves its terms and arithmetic order.
Paired local RK4 timings improved from 1.290 to 1.160 seconds per step (about
10% less time), with bitwise-identical output; all 23 tests pass. See the
performance note for isolated timings and allocation measurements. Production
evolutions remain unsubmitted.

A standalone Julia implementation of `sine_remainder` now passes direct
NumPy and high-precision checks. On identical 524,288-element inputs, the
allocating Julia kernel was about 9x faster for small fields, 8x for a sampled
Gaussian cloud slice, 2.6x for tiny tails, and 1.4x for mixed amplitudes in a
local benchmark. The reusable-output variant allocates no memory per call.
These timings exclude startup, compilation, and interoperability; Julia is
not yet connected to the Python evolution. Details and reproduction commands
are in `nonlinear-performance.md`.

The optional persistent Python–Julia bridge has now been tested end to end.
All 23 tests pass with the Julia remainder, and paired full-RK4 medians fell
from 1.481 to 0.948 seconds (1.56x speedup), including communication and array
copies. The benchmark temporarily substitutes the bridge; the production
solver still defaults to NumPy. This remains roughly 197 compute hours for
the benchmark's 1e5 M target, so further profiling is needed before long runs.
