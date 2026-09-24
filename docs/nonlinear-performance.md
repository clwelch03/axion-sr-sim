# Understanding the nonlinear runtime

> Dated development record. Later sections may supersede earlier next steps.
> For the implemented adiabatic approach and current work, see
> [adiabatic status](adiabatic-repro.md); time-domain campaign commands are in
> [Dolan Fig. 13](../runs/dolan/fig13/README.md).

These notes explain the measurements one small step at a time. Timings are
from the current machine on 2026-09-15, using the `BHsuper` environment
(NumPy 2.3.5), with BLAS/OpenMP thread counts set to one. They are not a
compute-node runtime guarantee.

## 1. Where is the time going?

The benchmark in [benchmark_nonlinear.py](../src/dolan/benchmark_nonlinear.py) uses
4267 radial nodes, six retained angular modes, and 32×64 angular points.
One RK4 step evaluates the nonlinear source four times: about 35 million
field values per step.

The original step took 4.66 seconds. A separate profiling call spent 4.09 of
4.56 seconds inside `sine_remainder`, which calculates

\[
R(\Phi)=\Phi-\sin\Phi.
\]

For small fields, subtracting these nearly equal numbers loses accuracy.
The code instead evaluates the series beginning with
$\Phi^3/6-\Phi^5/120$, through order $\Phi^{11}$.
We need to keep that accuracy when optimizing.

## 2. Why did changing the cube help?

The code already had `squared_field = small_field**2`. We changed
`small_field**3` to `small_field * squared_field`, reusing that result.

NumPy executes array operations separately. It does not automatically inspect
the whole function and recognize that an earlier square could be reused.
The power and multiplication expressions also invoke different operations.
The measurements establish that multiplication is much cheaper here; we
have not traced the compiled power implementation to establish its exact
machine-level cost.

For 524,288 values uniformly between −0.1 and 0.1, an isolated check measured
50.7 ms for cubing versus 0.56 ms for multiplication using the existing square.
Results agreed to floating-point roundoff. The full RK4 benchmark improved
from 4.66 to 1.71 seconds per step, since other work remains.

## 3. What is a temporary array?

Consider this expression:

```python
remainder = real_field - np.sin(real_field)
```

NumPy first makes an array containing the sine values, then makes another
array containing the subtraction results. The sine array can then be released.
Each operation traverses the data; writing and reading these arrays takes time
in addition to the arithmetic.

For one full radial block, the array shape is `(2048, 256)`: 2048 angular
points at each of 256 radii. Each float64 array occupies **4 MiB**.

Here are the main allocations in `sine_remainder` for an array input:

| Expression | What it creates |
|---|---|
| `np.abs(real_field) < 0.1` | Temporary absolute values and a Boolean mask |
| `real_field - np.sin(real_field)` | Sine values and a new result array |
| `real_field[small]` | A packed copy of the selected small values |
| `small_field**2` | Their squares |
| Polynomial expression | Intermediate arrays for the successive arithmetic operations |
| Former `np.array(remainder, copy=True)` | Another full copy of the result |

These arrays are not all alive simultaneously. Adding their sizes does not
give the peak memory usage.

## 4. One removable copy: measured and changed

The subtraction already produces a fresh result, so this subsequent line
was redundant for array inputs:

```python
remainder = np.array(remainder, copy=True)
```

It is now:

```python
remainder = np.asarray(remainder)
```

For an existing array, `asarray` reuses it. For a scalar result, it still
provides the zero-dimensional array needed for the later indexed assignment.
The caller's input remains untouched.

We compared the current function with a candidate differing only in that
line. Both were warmed up, then timed 15 times in alternating order, on
524,288 values generated with random seed 42:

| Input distribution | With copy | Without copy |
|---|---:|---:|
| Uniform between −0.099 and 0.099 | 30.84 ms | 28.09 ms |
| Uniform between −3 and 3 | 14.39 ms | 14.04 ms |
| Random signs, log-uniform magnitudes $10^{-160}$ to $10^{-10}$ | 19.85 ms | 19.36 ms |

These are modest improvements, roughly 2–9% in this run, not another factor
of three. Timing variability means the precise percentages should not be
treated as constants. Copying one 4 MiB array alone took about 0.35 ms.

`tracemalloc` measured the same peak, about 20.50 MiB, for both versions on
the all-small inputs. The mixed case decreased from 8.77 to 8.50 MiB.
Removing an allocation need not lower the peak: the later polynomial
temporaries can determine the maximum instead.

Outputs were bitwise identical on these inputs. Additional comparisons
covered scalars, empty arrays, and a strided array; inputs were preserved.
The full 23-test suite also passes after this change.

## 5. Reuse a polynomial work array: measured and changed

The small-field polynomial now uses one array, `polynomial`, updated in place.
Writing $q=\Phi^2$, we start with $q/39916800-1/362880$, then repeatedly
multiply by $q$ and add the next coefficient: $1/5040$, $-1/120$, and $1/6$.
Finally we multiply by $\Phi^3$. The selected `small_field` array is already
a copy made by Boolean indexing, so it can safely hold that cube. The input
field is not modified. The threshold and expansion through $\Phi^{11}$ are
unchanged, as is the arithmetic order.

A paired comparison against the immediately preceding implementation used
the same three 524,288-value distributions described above (seed 42), with
11 timings per implementation in alternating order and one BLAS/OpenMP thread:

| Input distribution | Previous expression | Reused buffer |
|---|---:|---:|
| Uniform between −0.099 and 0.099 | 30.11 ms | 22.90 ms |
| Uniform between −3 and 3 | 11.01 ms | 10.82 ms |
| Random signs, log-uniform magnitudes $10^{-160}$ to $10^{-10}$ | 29.70 ms | 23.16 ms |

`tracemalloc` peak fell from 20.50 to 16.50 MiB for the all-small and tiny
inputs; the mixed case remained at 8.50 MiB. Outputs were bitwise identical,
including additional scalar, empty-array, and strided-array checks, with
inputs preserved. All 23 tests pass.

For the full RK4 benchmark configuration in section 1, both implementations
were warmed and timed five times in alternating order on the same fixed
state. Median step time decreased from **1.290 s to 1.160 s**, about 10% less
wall time. Ranges were 1.269–1.311 s and 1.093–1.199 s respectively. Complete
RK4 output arrays were bitwise identical. These paired measurements are the
relevant comparison; absolute times differ from earlier runs on this machine.
This is a short compute-only benchmark, not a long evolution or a cluster
runtime guarantee.

## 6. What remains

This improvement does not resolve the long-run runtime problem. A subsequent
experiment could avoid calculating sine for entries whose result is replaced
by the polynomial; it must measure whether extra indexing costs outweigh the
saved sine evaluations. Angular matrix layout and a compiled source kernel
also remain candidates. None of those changes is included in this checkpoint.

## 7. Julia installed for the compiled-kernel experiment

The user chose Julia so the compiled loop remains readable to them. Julia
1.13.0 (official Linux x86_64 binary) is installed at
`/global/scratch/projects/pc_heptheory/clwelch/software/julia-1.13.0`.
The archive SHA256 matched the official release feed:
`8975da61c128a5e5ded3e719e868da8c8781deb7ad7913d37fb99be02a81904b`.
Startup and a small compiled arithmetic check passed with one thread.

To use it in a shell:

```bash
export PATH="/global/scratch/projects/pc_heptheory/clwelch/software/julia-1.13.0/bin:$PATH"
julia
```

Shell startup files were not modified. Verification used a temporary depot at
`/tmp/axion-julia-install-check`; no Julia packages or Python bridge have been
installed. No C kernel was written. The next step is the small Julia
`sine_remainder` loop and an accuracy/timing comparison, not a solver port.

## 8. Standalone Julia loop: accuracy and timing checkpoint

`src/dolan/sine_remainder.jl` now implements the same threshold and polynomial in a
compiled loop. Each iteration chooses the polynomial or sine subtraction,
then writes one result. No fast-math or explicit SIMD annotations are used.
`sine_remainder!` writes into a supplied output array; the allocating wrapper
`sine_remainder` creates an output like the NumPy function. The production
Python solver is unchanged and does not yet call Julia.

Reproduce the comparison (no additional Julia packages required):

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  PYTHONPATH=src python -m dolan.benchmark_julia_remainder \
  --julia /global/home/users/clwelch/.juliaup/bin/julia
```

The Python script creates identical float64 inputs for both languages via
temporary binary files, including 524,288-element small, mixed, and tiny-tail
distributions (seed 42), and a sampled real m=1 Gaussian field slice with
peak 0.7. Julia runs on one thread, with subnormal flushing disabled. Both
entry points are warmed before 11 alternating timings; NumPy is timed 11
times before and after the Julia process. Startup, compilation, file I/O,
and any future Python–Julia interoperability costs are excluded.

| Input | NumPy | Julia, allocating | Julia, reusable output |
|---|---:|---:|---:|
| Small fields | 27.099 ms | 3.026 ms | 1.609 ms |
| Mixed amplitudes | 11.153 ms | 8.163 ms | 7.796 ms |
| Tiny tails | 24.917 ms | 9.633 ms | 9.590 ms |
| Gaussian field slice | 15.455 ms | 1.907 ms | 1.923 ms |

Small-field and tiny-tail outputs matched NumPy exactly; largest absolute
differences were 4.44e-16 for mixed amplitudes and 1.11e-16 for the cloud
slice. All comparisons passed at relative tolerance 5e-13, with an absolute
allowance of eight minimum float64 subnormal units. Threshold-edge and
very-small-value checks also passed against 550-digit mpmath calculations.
Inputs were preserved; scalar-array, empty-array, and strided-view checks
passed. Julia's reusable-output calls reported zero allocated bytes; the
allocating calls used approximately the 4 MiB output size.

The representative cloud-slice kernel is about eight times faster here, but
this is **not** a full-RK4 speedup. The next decision is how to call Julia
without copying large arrays or starting a process for each call. That
integration and its end-to-end benchmark have not been implemented.

## 9. Persistent Python–Julia bridge: full-step checkpoint

`julia_remainder_bridge.py` and `julia_remainder_worker.jl` provide an
experimental connection using two shared memory-mapped files in a temporary
directory. The Julia process starts once and compiles before reporting ready.
Each call copies the input into the shared buffer, sends an element count,
waits for completion, and copies the output into an independently owned NumPy
array. Julia uses its standard-library `Mmap`; no additional packages are
required. This is a single-caller prototype, not a production parallel backend.
See [Julia's mapping semantics](https://docs.julialang.org/en/v1/stdlib/Mmap/).

The solver's default is unchanged. `benchmark_julia_bridge.py` temporarily
substitutes the bridge for `sine_remainder` only inside the experiment:

```bash
OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  PYTHONPATH=src python -m dolan.benchmark_julia_bridge --check-tests \
  --julia /global/home/users/clwelch/.juliaup/bin/julia
```

All 23 existing tests passed with Julia handling the remainder. Additional
checks covered scalar, empty and strided inputs, repeated buffer reuse,
preservation of inputs and prior outputs, and rejection of complex input.
The complete NumPy and Julia-backed RK4 states agreed at relative tolerance
2e-12 and absolute tolerance 1e-14.

For the same 4,267-node, six-channel, 32x64 angular benchmark, five warmed
calls per implementation were timed in alternating order:

| Full RK4 step | Median | Range |
|---|---:|---:|
| NumPy remainder | 1.481 s | 1.382–2.521 s |
| Julia bridge | 0.948 s | 0.883–1.973 s |

The median speedup was **1.56x** (36% less wall time), including buffer copies
and per-call communication. Both processes use one thread; Python waits for
Julia to complete each call. Startup/compilation took 3.13 s and is excluded
from steady-state timing. The broad timing ranges show local runtime
variability. The experiment excludes checkpoint I/O and uses repeated steps
from a fixed Gaussian state, not a full cloud evolution.

The resulting compute-only estimate is still about 197 hours to reach 1e5 M
at dx=M/6 and dt=0.8 dx. No production runs have been submitted. A fresh
profile of this configuration is the next performance checkpoint: the
remaining angular calculations, communication, and tiny-field arithmetic
must be measured before selecting another optimization.

## 10. Profile of the Julia-backed step

`benchmark_julia_bridge.py --profile` now profiles three additional warmed
steps with Python's `cProfile`. The numerical implementation is unchanged.
Using the same configuration and single-thread settings as section 9 gave
2.619 s total, or 0.873 s per profiled step:

| Work (exclusive categories) | Seconds per step | Share |
|---|---:|---:|
| NumPy source calculations outside the bridge | 0.441 | 50.5% |
| Julia bridge, including computation, waiting and copies | 0.410 | 46.9% |
| Remaining work, including radial RHS and RK4 bookkeeping | 0.022 | 2.6% |

The radial tower RHS alone took 0.017 s per step (1.9%). The bridge was
called 68 times per step: 17 radial blocks at each of four RK4 stages.
Waiting for Julia accounted for 0.329 s per step, but this includes Julia
computation and scheduling as well as communication; Python profiling
cannot separate these. The source's own time includes reconstruction,
projection, weights and array operations; this profile does not distinguish
those individual operations either.

The preceding unprofiled measurements gave medians 2.438 s for NumPy and
1.157 s for the bridge (three calls each, ranges 1.995–2.936 s and
0.911–1.384 s). Startup/compilation took 51.963 s, excluded throughout.
Complete-step state agreement passed at rtol=2e-12, atol=1e-14.
Timing variability remains substantial; these fixed-state measurements
are not a production evolution forecast.

The next small experiment should target angular matrix products: compare
the existing strided real/imaginary views with contiguous real arrays,
checking identical source values and timing before changing the solver.
This is a candidate optimization, not a demonstrated bottleneck within
the source's own time. No production runs were submitted.

## 11. Contiguous matrix experiment

The proposed layout change showed no measurable full-step benefit. The
optional `contiguous_arrays=True` argument to `nonlinear_source_on_grid`
packs harmonic and projection real/imaginary arrays once per source call,
and coefficient real/imaginary arrays once per radial block. All packing
costs are included; the default remains the original strided layout.

Reproduce with the section 9 environment and Julia executable, adding
`--compare-layouts --repeats 7` to `benchmark_julia_bridge.py`. Both variants
use the Julia remainder, with warmed steps timed in alternating order:

| Layout | Median RK4 step | Range |
|---|---:|---:|
| Strided | 1.527607 s | 1.344438–2.009898 s |
| Contiguous | 1.527888 s | 1.336625–1.646281 s |

The speedup was 1.000x. Random complex coefficients populated all six channels
at three coefficient scales (0.01, 0.7, 3.0; these are not peak field
amplitudes); the source arrays matched exactly in these checks. Complete
RK4 states agreed at rtol=2e-12, atol=1e-14. All 23 regression tests also
passed with the contiguous source and NumPy remainder. Julia startup took
2.214 s and was excluded. As before, this is a fixed Gaussian-state benchmark,
not a long cloud evolution.

This result does not support adopting the packing change. Separately timing
reconstruction, remainder evaluation and projection would narrow the next
optimization target further than the function-level profile.

## 12. Timing the individual source operations

`benchmark_julia_bridge.py --profile-source --repeats 5` uses targeted Python
line tracing on the existing source function, with a timer around the Julia
bridge. The bridge duration is subtracted from the enclosing weighting
expression, so the categories below do not double-count it. No numerical
implementation changed. Use the executable and single-thread environment
from section 9.

Five warmed, instrumented RK4 steps averaged 1.312 s:

| Operation | Mean seconds per step | Share of total |
|---|---:|---:|
| Reconstruct real Phi | 0.340 | 25.9% |
| Remainder bridge (Julia work, copies, communication and waiting) | 0.685 | 52.2% |
| Apply geometric weight to remainder | 0.150 | 11.4% |
| Project onto harmonics | 0.111 | 8.4% |
| Source setup and final assembly | 0.004 | 0.3% |
| Remaining solver work and instrumentation overhead | 0.023 | 1.8% |

Each step made 68 bridge calls. Instrumented and uninstrumented steps were
interleaved, and every resulting state agreed exactly with the warmed
reference. Untraced median/range: 1.270 s / 0.853–1.633 s; traced:
1.396 s / 0.857–1.584 s. Both paths include the bridge timer wrapper.
Timing drift is substantial, so the median difference does not isolate
tracing overhead. Julia startup took 2.055 s and was excluded.

These measurements identify reconstruction as the largest NumPy portion,
but the bridge remains the largest overall category. A timer inside the
Julia worker would distinguish kernel execution from the rest of the
bridge before choosing another implementation change. This profile uses
the same fixed Gaussian seed as previous benchmarks; it does not establish
performance for evolved or bound-mode clouds.
