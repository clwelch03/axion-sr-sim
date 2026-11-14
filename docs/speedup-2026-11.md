# Performance project, November 2026

This is the record of the speedup of the AxionCloud time evolution (`evolve_with_spindown`, `evolve_fixed_background`),
done on 2026-11-12 and 2026-11-13. It is written for the code's owner and for whoever later cleans the code up. It
says what changed, why, how each change was checked, and what is left. It was first written on 2026-11-13 and brought
up to date on 2026-11-14, after the pieces that were then unfinished had finished.

**Status on 2026-11-14.** The work was merged into `rewrite/axioncloud` on 2026-11-14 as pull request #1. Before the
merge, the branch took in a change that the owner had meanwhile made on `rewrite/axioncloud`, the conjugate-frequency
rule; section 4.1 says what that rule is and how the checks were redone for it. The comparison of time integrators is
finished (section 5). The numerical-parameter audit and the study of scaling with the number of bound modes are
finished and have their own documents, `docs/parameter-audit.md` and `docs/mode-scaling-study.md` (sections 7 and
8). The owner has answered the open questions of the first version; section 8 records the decisions and what came of
them.

The work was done by several Claude agent sessions, each named after its task (assembly, linear solver, Newton,
balance, integrator, and "polish" for the second round of cleanups), coordinated by one more session. Measurements are
attributed to the agent that made them, because they were taken on different nodes and at different stages.

How to read section 2. Each subsection describes one change in five parts: what this part of the calculation does;
why it was slow (or, for an accuracy item, what was wrong); what changed; what it bought, in numbers; and how we know
the answer did not change. Each computing term is explained where it first appears and then used. A closing line in
italics, "For the cleanup pass", carries the branch names, commit hashes, job IDs and function names; you can skip
those lines on a first reading.

*For the cleanup pass: merged as https://github.com/clwelch03/axion-sr-sim/pull/1, merge commit `6c0a9e6` on
`rewrite/axioncloud`. The merged branch is `speedup/integration-next`: `9e17100` was the state described in the first
version of this document, `057a167` its last commit before the merge with the target, `582b955` the merge with
`rewrite/axioncloud` at `6bec64a`, and `86c4b62` the rewritten regression reference. The conjugate-frequency rule is
commit `41fc2e7`. This document, the parameter audit and the scaling study are on `speedup/followup` (not merged, not
pushed). The coordinator's running record is `STATUS.md` in `axion-sr-sim-speedup/` (outside the repository); the
rules given to every agent are in `AGENT_BRIEF.md` there. This document is built from those two files, the agents'
reports and the code on `speedup/followup`, which contains the merged code.*

## If you read only this

**What was achieved.** The time evolution spends almost all its time computing the rates of change of the state
(amplitudes, black-hole mass and spin) at the points the time integrator asks for. One such computation now costs
0.19 to 0.22 s on one processor core instead of 2.99 s, a factor of about 15. On eight cores it costs about 0.07 s, a
factor of about 42; one earlier run on another node gave 0.11 s. A full rerun of your saved window
$t = 1.5\times10^9 \to 10^{11}$ took 295 s on eight cores. The original run took 13353 s, but on a different machine
with a different core layout, so the factor of about 45 is not a like-for-like comparison. With the stiff integrator
FBDF (section 5) the same window took 119 s.

| | Before | Now, 1 core | Now, 8 cores |
|---|---|---|---|
| one rate computation (plateau state) | 2.99 s | 0.189 to 0.224 s | 0.068 to 0.106 s |
| short evolution, 31 rate computations | 94.3 s | 6.2 to 7.7 s | 2.81 to 3.15 s |
| window $1.5\times10^9 \to 10^{11}$, Tsit5 at tolerance $10^{-7}$, 2233 to 2300 rate computations | 13353 s (other machine) | not measured | 266 to 296 s |
| same window, FBDF at tolerance $10^{-10}$, 1036 radial solves | | not measured | 119 s |

The ranges are over runs on different nodes and commits (section 1); the 8-core figures vary more between nodes than
the 1-core ones.

**Two things about your existing results.**

1. The old code computed the amplitude rates with a slightly out-of-date matrix. Individual rate evaluations were
   wrong by up to 1.1%. In the stored linear-growth window this left $A_{322}$ at $t = 1.5\times10^9$ off by
   $7.6\times10^{-5}$ relative. Every later window that was restarted from that state carries the same offset as its
   initial condition. Whether that offset changes late-time results has not been tested (section 3.1). You are
   rerunning everything with the new code (section 8).
2. At $\alpha = 0.187$ the radial solve cannot reach its default convergence tolerance, because of roundoff. With the
   old code every solve there ran the maximum 30 Newton iterations, about 90 s. The new code stops early and reports
   the solve as not converged; to run at $\alpha = 0.187$ you pass one extra keyword (section 3.2).

**Which integrator.** Tsit5 stays the default. On long stiff plateaus FBDF is much cheaper: use it at tolerance
$10^{-9}$ for plateaus and $10^{-10}$ through a $|322\rangle$ takeover. Neither Rosenbrock method is recommended
(section 5).

**More bound modes.** The speedup does not change how the cost grows with the number of bound modes $K$. Measured:
one radial solve takes 0.10 s at $K = 2$, 1.5 s at $K = 3$ and, from measured pieces, about 47 s at $K = 4$ (eight
cores). Of order ten fully nonlinear modes is out of reach even after every change the scaling study considered. Its
recommendation is two to four modes treated as now plus tens of weak modes, each treated as a separate linear problem
(section 7).

**Decisions** (section 8): the metric-work diagnostic stays off by default; the parameter audit was done and
recommends a cheaper phase and polar quadrature and a check of the angular truncation, its largest error; FBDF on
stiff runs, Tsit5 otherwise; you are rerunning old results; the remaining serial part of the computation was not
worked on.

## Terms used throughout

- **Rate computation**, or **right-hand-side (rhs) evaluation**. The time integrator advances the state
  ($\ln A_i$, the black-hole mass $b$ and angular momentum $j$). At each of its internal stages it asks for the rates
  of change of that state. Producing those rates for one state is one rhs evaluation: one `solve_radial` (Newton's
  method for the radial boundary-value problem at the current amplitudes, $b$ and $j$) followed by one
  `balance_terms` (the rates of $\ln A_i$, $b$ and $j$ from the fluxes). Total run cost is the number of evaluations
  times the cost of one.
- **Harness**: the benchmark and regression scripts in `benchmarks/`. They time each piece of an rhs evaluation on
  fixed stored states and compare the results with the baseline's (section 4.1).
- **A/B runs**: two runs that differ in one setting only, to measure the effect of that setting.
- **lr7**: the Lawrencium Slurm partition (Xeon Gold 6330 nodes) on which all quoted timings were taken, unless
  stated otherwise.
- **Thread** and **core** (2.1), **allocation** and **garbage collection** (2.1), **in-place function** and
  **workspace** (2.1), **race condition** (2.1), **dual number** and **automatic differentiation** (2.1),
  **roundoff** and **machine precision** (2.1).
- **Sparse matrix**, **LU factorisation**, **back-substitution**, **partial pivoting**, **almost-block-diagonal
  matrix** (2.2).
- **Chord (modified Newton) step**, **line search**, **implicit-function theorem for sensitivities**, **predictor**
  (2.3).
- **Stiffness**, **explicit, implicit, Rosenbrock and multistep methods**, **order**, **W-method**, **L-stable**
  (2.7 and 5); **work-precision comparison** (5).
- **Regression test**, **noise floor**, **unit test** (4.1, 4.2).

## 1. Summary

Cost of one rhs evaluation, and of a short evolution segment, on the benchmark's plateau state. The plateau state has
amplitudes $(2.8077\times10^{-3},\,3.396\times10^{-9})$: it is the $|211\rangle$ plateau at $t = 1.5\times10^9$. All
runs on lr7 with the matrix library restricted to one core (section 2.8). The baseline is the code before the project,
on one core.

| Piece | Baseline, 1 core | New, 1 core | New, 8 cores |
|---|---|---|---|
| spindown rhs call, `9e17100` (coordinator's harness run) | 2.99 s | 0.189 s | 0.106 s |
| evolution $1.5\times10^9 \to 1.9\times10^9$ (4 steps, 31 calls), `9e17100` | 94.3 s | 6.4 s | 3.15 s |
| spindown rhs call, `2246b4f` (polish agent's run) | 2.99 s | 0.206 s | 0.069 s |
| evolution segment, `2246b4f` (polish agent's run) | 94.3 s | 7.00 s | 2.84 s |
| spindown rhs call, `057a167` (last commit before the merge) | 2.99 s | 0.195 s | 0.071 s |
| evolution segment, `057a167` | 94.3 s | 6.23 s | 2.87 s |
| spindown rhs call, merged with the target (`582b955`) | not timed | not timed | 0.068 s |
| evolution segment, merged with the target (`582b955`) | not timed | not timed | 2.81 s |
| spindown rhs call, large mesh ($\alpha = 0.187$), `2246b4f` | 6.10 s | 0.402 s | 0.134 s |
| spindown rhs call, large mesh, `057a167` | 6.10 s | 0.381 s | 0.140 s |
| memory allocated per rhs call (2.1, 2.6), `057a167` | 6.6 GB | 56 MB | 90 MB |

At `057a167`, the last commit before the merge, the speedup per rhs call is 15 times on one core and 42 times on
eight. The 8-core figures differ between nodes: 0.106 s and 3.15 s in the coordinator's run of `9e17100`, 0.068 s to
0.071 s and 2.81 s to 2.87 s in the later runs. The 1-core figures agree between runs much better, but not perfectly:
a run at `ebf0ae9`, one commit before `057a167`, whose single-evaluation path does the same work, gave 0.224 s and
7.74 s on another node, and the `057a167` run gave 0.195 s and 6.23 s. That spread of about 15% is node-to-node
variation. The baseline of the merged code with the conjugate-frequency rule (`6bec64a`, the target before the merge)
was not timed; the rule changes the imaginary parts of some component frequencies, not the amount of work per
evaluation (my reading of the commit).

The regression test (168 quantities compared with the baseline's values, section 4.1) passes at one and at eight
cores at `057a167`, and on the merged code against a reference rewritten from the target's own original solver
(section 4.1). All 305 unit tests pass at `057a167`, and all 307 on the merged code. The costs of the separate pieces
(residual, Jacobian, factorisation, `balance_terms`) are in section 2 and in `benchmarks/README.md`.

End-to-end rerun of the saved restart window $t = 1.5\times10^9 \to 10^{11}$ at $f_a/M_{\mathrm{Pl}} = 0.0466873$
($f_a = 5.7\times10^{17}$ GeV), $a_0 = 0.9$, $\alpha = 0.3$, cosine potential, Tsit5. Every version took the same 2233
rhs evaluations.

| Code | Wall time | Hardware |
|---|---|---|
| original (Perlmutter job 61110567) | 13353 s | Perlmutter, 1 Julia thread, 4 BLAS threads |
| `01856b9`, first full merge (job 27820263) | 572 s | lr7, 8 Julia threads |
| `9e17100` (A/B runs, jobs 27821575 and 27821576) | 295 s | lr7, 8 Julia threads |

The 295 s figure is two runs: 294.6 s with a since-removed correction to chord steps (2.3.5) and 296.4 s without it
(`9e17100`). Both had 2233 calls, 371 accepted and 1 rejected step, and no failed solve. The first-merge run's
trajectory (61 saved times) agrees with the stored CSV to $1.4\times10^{-8}$ relative in $A_{211}$,
$7.4\times10^{-9}$ in $A_{322}$, $3\times10^{-11}$ in $b$ and $7\times10^{-11}$ in $j$. The integrator tolerance was
$10^{-7}$. The two A/B runs agree with each other in final state to $10^{-7}$ or better. The `9e17100` state
at $t = 10^{11}$ agrees with the stored final row to $1.4\times10^{-8}$ in $A_{211}$ and $1.2\times10^{-9}$ in
$A_{322}$. That is a comparison of the final state only, not of the 61 saved times.

The end-to-end ratio, 13353 s / 295 s or about 45, is not like-for-like: the original ran on different hardware with a
different thread layout. The harness numbers in the first table are like-for-like.

The integrator comparison (section 5) timed the same window again on the code with FBDF support (`a84f712`, the
`9e17100` solver plus the integrator fixes): Tsit5 at tolerance $10^{-7}$ took 266 s for 2300 radial solves (a few
more than 2233, because that script forces steps to land on its output times), and FBDF at $10^{-10}$ took 119 s for
1036. The end-to-end run was not repeated on the merged code with the conjugate-frequency rule.

Most important caveats:

1. **The old code's rates were wrong by up to about 1% in some evaluations** (section 3.1). The stored linear-window
   result has $A_{322}$ off by $7.6\times10^{-5}$ relative at $t = 1.5\times10^9$. Later windows restarted from it
   inherit that offset. The new code differs from a converged run by about $4\times10^{-7}$.
2. The speedup does not grow past about 8 cores (section 2.8). The factorisation, which runs on one core, and one
   boundary component set a floor.
3. The checks establish that the new code reproduces the old code's numbers within measured noise floors on four
   stored states and one short evolution, and that one full window reproduces; after the merge, that it reproduces
   the target's own solver, with the conjugate-frequency rule, on the same four states and short evolution. They do
   not establish discretisation accuracy, other parameter points, other potentials, the energy balance
   (`balance=:energy`), or more than two bound modes (section 4). The parameter audit addresses discretisation
   accuracy separately (section 8, item 2).
4. The full window was not rerun on the merged code.

*For the cleanup pass: harness jobs at `9e17100`: regression 27821631 (1 thread) and 27821632 (8 threads), unit
tests 27821630. At `ebf0ae9`: 27823848, 27823849. At `057a167`: regression and timing 27824651 (1 thread) and
27824652 (8 threads), logs in `integration-next/benchmarks/logs/` under `axion-sr-sim-speedup/`; unit tests 27824655.
Merged code (`582b955`): regression against the new reference and timing 27846902 (8 threads, log in
`pr-work/logs/`), unit tests 27846903. The 1-core and 8-core timings quoted are the plateau rows of these logs;
STATUS.md quotes the linear-state rows (0.188 s and 0.227 s at 1 thread), which differ by a few percent. Final-state
comparison against the stored CSV: jobs 27821575/27821576, logs `review-work/rescale_ab/ab_window_*.out`. First-merge
outputs in `review-work/e2e_01856b9/`. Integrator re-timing: `methods-work/logs/retime_*`.*

## 2. How the time is spent and what changed

### 2.0 Where the time went

One rhs evaluation is one radial solve followed by `balance_terms`. The radial solve is Newton's method. Each Newton
step does four things:

1. evaluate the residual, the list of how badly the current guess fails each discretised equation (2.1);
2. build the Jacobian, the matrix of derivatives of the residual with respect to every unknown (2.1);
3. factorise the Jacobian (2.2);
4. use the factors to solve for the Newton correction (2.2).

Baseline cost of each piece on the plateau state, one core:

| Piece | Time | Note |
|---|---|---|
| residual | 0.061 s | |
| Jacobian | 0.523 s | 1435 MB of memory allocated |
| factorisation | 0.122 s | |
| `solve_radial` after a small perturbation | 1.81 s | 4 Newton iterations |
| `solve_radial` at an already-solved point | 1.81 s | 1 iteration |
| `balance_terms` | 0.984 s | |
| whole rhs call | 2.99 s | |

Three things stand out. A radial solve was dominated by building and factorising the Jacobian. A solve at a point
whose solution was already known cost as much as a real solve. And `balance_terms` was a third of the call.

The changes below are grouped by what they touch. They were merged in the order assembly, Newton iteration, linear
solver, balance terms, integrator, followed by a second round of cleanups found by measuring the merged code.

*For the cleanup pass: baseline harness job 27815524. Merge order onto `speedup/integration`: `speedup/assembly`,
`speedup/newton`, `speedup/linear-solver`, `speedup/balance`, `speedup/integrator`, then `speedup/polish`.*

### 2.1 Building the residual and the Jacobian

**What this part of the calculation does.** The radial solve finds the cloud's profile on a mesh of 690 radial
points. At each point there are 8 field components (channels), and for each channel the code stores the complex
profile $R$ and its tortoise-coordinate slope $Q$, so 32 real numbers per point. Together with 32 parameters (the
complex mode frequencies and the leading series coefficients at the two ends) that makes 22112 unknowns. On each mesh
interval of width $h$ the unknowns must satisfy a fourth-order discretisation of the radial ODE $y' = f(y)$, Simpson's
rule with a Hermite midpoint:

$$y_R - y_L - \tfrac{h}{6}\,(f_L + 4 f_{\mathrm{mid}} + f_R) = 0, \qquad
y_{\mathrm{mid}} = \tfrac12 (y_L + y_R) - \tfrac{h}{8}\,(f_R - f_L).$$

Here $L$ and $R$ label the two ends of the interval and $f$ includes the self-interaction source, averaged over the
phases and polar angle. The two ends of the domain add matching conditions to the known series there. Newton's method
needs, at every step, the **residual** (the left-hand sides above, 22112 numbers) and the **Jacobian** (the derivative
of each residual entry with respect to each unknown). Building these two objects is called **assembly**. It is a loop
over the mesh points and intervals.

**Why it was slow.** Three separate reasons.

1. *Allocation.* At every mesh point the old code created fresh small arrays for intermediate results, used them once
   and dropped them. Creating an array is called **allocation**: the Julia runtime must find free memory for it. Memory
   that is no longer used is reclaimed later by the **garbage collector**, which scans memory to find what can be
   freed; while it runs, it pauses every thread of the program. The old code allocated 1435 MB of short-lived arrays
   per Jacobian, and finding and reclaiming that memory took a large share of the time.
2. *One core.* A processor has several **cores**, each of which can run its own stream of instructions. Julia calls
   the streams it can run at once **threads** and places each on a core (`julia -t 8` gives it 8 threads). The work
   at one mesh point does not need the result at any other point, so different threads could handle different
   stretches of the mesh at the same time. The old code used one thread and visited the points in turn.
3. *Differentiating with respect to more inputs than needed.* The old code obtained the derivatives of the outer
   boundary condition by **forward-mode automatic differentiation**, with **dual numbers**. A dual number carries a
   value together with its partial derivatives with respect to chosen inputs, $x + \sum_k \partial_k x\,\varepsilon_k$,
   and its arithmetic is defined so that every operation applies the chain rule to those partials. Running ordinary
   code on dual numbers therefore returns derivatives as well as values. Each operation costs roughly $1 + n$ times
   a plain one, where $n$ is the number of partials carried. The old code carried all 64 inputs of the boundary
   condition for every Fourier component. This alone took 0.35 s of the 0.58 s the Jacobian took (assembly agent's
   measurement).

**What changed.**

1. *In-place functions and workspaces.* The kernels (the functions doing the per-point work) are now **in-place**:
   instead of returning a new array they overwrite an array they are given. Each thread gets one set of scratch
   arrays, a **workspace**, allocated once and overwritten at every mesh point. The workspace types are
   `FieldWorkspace`, `RhsWorkspace` and `NodeWorkspace`; the Jacobian goes into a `JacobianBlocks`, which stores only
   its nonzero blocks (2.2). The inner loops no longer allocate. The old functions remain as thin wrappers with
   unchanged signatures.
2. *Threads.* The mesh is cut into consecutive stretches (chunks), and each chunk is a Julia **task**, a unit of work
   that the runtime schedules on any free thread. Each task has its own workspace. The rows of the boundary conditions
   run as one more task at the same time. If two tasks wrote into the same memory, each could overwrite the other's
   intermediate numbers, and the result would depend on timing. That is a **race condition**. The first threaded
   version had one: a variable captured by the parallel loop body was shared between tasks instead of private to each
   (a "closure-capture race"). It was found and fixed, and there is now a test that compares threaded results with
   sequential ones. Chunks are contiguous and every sum is taken in a fixed order, so the result does not depend on the
   number of threads (see "How we know" below).
3. *Differentiating each component separately.* Each Fourier component's far-field form depends only on its own
   complex frequency and its own leading coefficients: two real numbers for the frequency and two per channel of that
   component. Its boundary condition is therefore differentiated with dual numbers carrying only those partials. All
   other partials of that condition are exactly zero, so nothing is lost; only the work on zeros is skipped. The
   components run as concurrent tasks.
4. *Hand-coded chain rule for the source.* One more derivative is needed elsewhere (by `balance_terms`, 2.4): the
   derivative of the residual with respect to the black hole's mass and spin, $\partial F/\partial b$ and
   $\partial F/\partial j$. The old code ran the whole residual on dual numbers with 2 partials. The new code computes
   the self-interaction source's value with the ordinary kernel and its two partials from the analytic derivative of
   the source formula, instead of pushing dual numbers through every operation.

**What it bought.** Assembly agent's own scripts, lr7, before the benchmark harness existed:

| Piece | Before, 1 core | After, 1 core | After, 8 cores |
|---|---|---|---|
| Jacobian | 0.58 s | 0.102 s | 17.7 ms |
| residual | 0.0885 s | 0.0237 s | 4.1 ms |
| residual derivative with respect to $b$, $j$ | 0.407 s | 0.048 s | 7.9 ms |
| outer-boundary part of the Jacobian | 0.35 s | 0.033 s | |

**How we know the answer did not change.** First, a word on how computers store numbers. A 64-bit floating-point
number holds about 16 significant decimal digits; the relative spacing between neighbouring representable numbers,
the **machine precision**, is about $2\times10^{-16}$. Every arithmetic operation rounds its result to that grid, an
error called **roundoff**. Floating-point addition is not associative, so adding the same numbers in a different
order can change the last digit.

The old and the new code were fed identical inputs. The plain (non-dual) kernels agree to $4.6\times10^{-16}$ relative,
which is rounding in the last digit; the residual is identical bit for bit. The hand-coded mass-and-spin partials agree
with the old dual-number ones to $4.3\times10^{-15}$. The results are also identical bit for bit whether 1, 2, 4, 8 or
16 threads are used, which is the check that the threads do not interfere with one another.

*For the cleanup pass: branch `speedup/assembly`, commits `80b496f` (in-place kernels, threading, per-component
boundary Jacobian) and `83e8bb3` (dual-number source with analytic partials). Types `FieldWorkspace`, `RhsWorkspace`,
`NodeWorkspace`, `JacobianWorkspace`, `JacobianBlocks`; chunking in `_foreach_in_chunks` (`src/bvp.jl`); per-component
boundary derivative in `boundary_jacobian!`; files `src/bvp.jl`, `src/field.jl`, `src/radial.jl`. No `threadid()` is
used.*

### 2.2 Solving the linear system: block LU factorisation

**What this part of the calculation does.** Each Newton step solves the linear system $J\,\delta u = -F$ for the
correction $\delta u$, with $J$ the 22112-by-22112 Jacobian and $F$ the residual. The standard way is **LU
factorisation**: Gaussian elimination written as $J = P L U$, with $L$ lower triangular, $U$ upper triangular and $P$
a permutation of the rows (explained below). Computing the factors is the expensive step. Once they exist, solving for
any right-hand side is cheap: one forward substitution with $L$ and one **back-substitution** with $U$. The code
therefore factorises once and reuses the factors, both for the Newton correction and for the amplitude sensitivities
that feed the rates (2.3.4).

Most entries of $J$ are zero. A matrix stored with only its nonzero entries is a **sparse matrix**; the positions of
the nonzeros are its **sparsity pattern**. Here the pattern is very regular. The equations of mesh interval $j$
involve only the unknowns at nodes $j$ and $j+1$ (32 each) and the 4 real frequency unknowns. So the rows come in
blocks of 32, each touching two neighbouring node blocks: a block bidiagonal pattern, plus a horizon block at one end,
an outer block at the other, and the few frequency columns that touch every block. Matrices of this shape arise from
every discretised boundary-value problem and are called **almost block diagonal**.

**Why it was slow.** The old code passed the Jacobian to UMFPACK, a general sparse solver. A general solver does not
know the pattern in advance; it spends time analysing it and moving entries around. This took 0.12 to 0.13 s on the
production mesh and 0.259 s on the large mesh. Dense LU of a matrix this size would cost about $\tfrac23 n^3$, roughly
$10^{13}$ operations (my arithmetic), so a sparse method is essential; the question is how much the known pattern
helps.

**What changed.** A new solver, `BlockFactorization`, does Gaussian elimination in stages that follow the pattern. The
cost is linear in the number of mesh nodes times the cube of the block size, and nearly all of it is dense
matrix-matrix work on blocks of about 50 by 70.

How the stages work. The columns are eliminated in the order: horizon coefficients, node 1, node 2, ..., node $N$,
frequencies, outer coefficients.

1. Eliminating the 16 horizon-coefficient columns from the horizon block leaves 16 rows that involve only node 1 and
   the frequencies.
2. Stage $j$ stacks those 16 carried rows on the 32 rows of interval $j$, giving a 48-row panel whose columns are node
   $j$, node $j+1$ and the frequencies. It eliminates node $j$'s 32 columns from this panel with a routine from
   LAPACK, the standard library of dense linear-algebra routines. That leaves 16 rows that involve only node $j+1$ and
   the frequencies, and these are carried to the next stage.
3. The last stage is a dense 48-by-48 solve for the last node, the frequencies and the outer coefficients.

The 4 frequency columns touch every interval, so they are carried along as extra trailing columns of every panel.

**Partial pivoting, and why this is still full partial pivoting.** Gaussian elimination divides by the diagonal entry
(the pivot) of the current column. If that entry is small, the division magnifies roundoff in everything that follows.
**Partial pivoting** avoids this: before eliminating a column, it swaps rows so that the pivot is the entry of largest
magnitude in that column among the rows not yet used. The row swaps are the permutation $P$ above. At stage $j$ the
only rows with nonzero entries in node $j$'s columns are the 48 rows of the panel. A row with a zero in the column can
never be the largest entry. So choosing the pivot among the panel's rows makes exactly the same choice as choosing it
among all remaining rows of the whole matrix. The staged elimination is ordinary partial-pivoted LU of $J$ with a
particular column order, not an approximation.

**Why this matters for the radial equation.** A tempting alternative is to propagate from one end with one transfer
matrix per interval, $x_{j+1} = M_j x_j$, as in shooting. The radial equation has growing and decaying solutions, like
$e^{\pm\kappa r}$ for a bound component. The bound state is the decaying one, and any roundoff that excites the
growing solution is amplified across the domain; a transfer-matrix recursion also has no pivoting. The staged
elimination imposes both boundary conditions at once and pivots within each block, so it has the same stability as LU
of the whole matrix. The file's header comment says so, and the independent review checked it (section 4.4).

Two later improvements. The factorisation is now built directly from the Jacobian's blocks, so no general sparse
matrix is ever formed; that removed an 8 ms copy. The copying of blocks into the stage panels (with the Newton row and
column scaling weights of section 3.2 applied) runs on several threads. The elimination stages themselves run one
after another, because each needs the rows left by the previous one.

**What it bought.** Linear-solver agent's own runs, lr7, matrix library on one core:

| Piece | UMFPACK | Block LU |
|---|---|---|
| factorisation, 22112 unknowns | 0.123 to 0.133 s | 0.023 s |
| factorisation, 1355-node mesh | 0.25 s | 0.049 s |
| back-substitution, one right-hand side | 2.2 ms | 1.8 ms |
| back-substitution, 4 right-hand sides | 8.4 ms | 4.0 ms |

**How we know the answer did not change.** The **backward error** of a computed solution is the smallest relative
change of the matrix that would make that solution exact; it is the standard measure of a linear solver's accuracy.
On four states the block solver's backward error is at or below UMFPACK's. Newton iteration counts are unchanged, and
the unknowns differ by at most $4\times10^{-16}$ relative, which is roundoff.

Not established: behaviour near a stalled or ill-conditioned Newton iterate, and states other than the four tested. A
standard LAPACK routine for banded matrices was measured on a matrix of the same shape (0.063 s) and not pursued. The
sparse solver KLU was not tested, because it is a new dependency.

*For the cleanup pass: branch `speedup/linear-solver`, commits `05feccb`, `cf853b6`; `src/block_factorization.jl`,
behind `ScaledFactorization`. Building from blocks: `50ddeb9`. Threaded panel copy: `_foreach_in_chunks` at
`src/block_factorization.jl:196`. Per-stage routine: LAPACK `getrf` on 48-by-32 panels. Banded routine measured:
`gbtrf`.*

### 2.3 Newton iteration

**What this part of the calculation does.** `solve_radial` finds the unknowns $u$ that make the residual $F(u)$ zero
by Newton's method: $u \leftarrow u - J(u)^{-1} F(u)$. Near the solution Newton's method converges quadratically: the
error roughly squares at each step ($10^{-3}$, $10^{-6}$, $10^{-12}$). Two standard modifications save work.

- A **chord step** (the modified Newton method) reuses an older factorised Jacobian instead of building and
  factorising a new one. It costs one residual and one substitution instead of a full Jacobian, so it is much cheaper.
  Because the old matrix is only approximately right, the error falls by a roughly constant factor per step (linear
  convergence) instead of squaring. The code takes chord steps while each one at least halves the residual, and
  refactorises when one fails to.
- A **line search** guards against a step that makes things worse. If the full step does not reduce the residual, the
  code tries half the step, then a quarter, and so on (backtracking).

The evolution calls `solve_radial` at states close to the previous one, so the previous solution is a good starting
guess (a warm start) and its factorisation is a good chord matrix.

This subsection has five changes. Each has its own "why", "what changed", numbers and check.

#### 2.3.1 A solve at an already-solved point

**Why it was slow.** A solve started at the guess's own amplitudes, where the answer is already known, took 1.85 s,
as long as a real solve. Two reasons. The code records how far the point where the stored factorisation was computed
lies from the current iterate. For this case that distance was initialised to infinity instead of zero, which forced a
new factorisation at the end. And at the roundoff floor (2.3.3) no step can reduce the residual, so the line search
backtracked 13 times before giving up.

**What changed.** Both are fixed. Such a solve is now 2 residuals and one substitution.

**What it bought.** 1.85 s to 0.13 to 0.16 s (Newton agent's measurement). In the harness at `2246b4f`: 1.81 s to
0.048 s on one core.

**How we know the answer did not change.** The change affects only when the solve stops, not what it converges to.
The regression test, which compares re-solved quantities with the baseline's, passes.

#### 2.3.2 Starting Newton from a prediction

**Why it was slow.** Each evolution step started Newton from the previous solution, rescaled to the new amplitudes.
That start is off by first order in the change of the amplitudes, $b$ and $j$, so Newton needed several iterations
and often a new Jacobian.

**What changed.** Write $p$ for the parameters the solution depends on (the amplitudes, $b$ and $j$; my notation). The
solve now starts from $u(p) + (du/dp)\,\delta p$, a first-order Taylor extrapolation to the new state. The derivatives
$du/dp$ are already computed by `balance_terms` (2.3.4 explains how), so this needs no extra residual evaluations. The
neglected term is $O(\delta p^2)$, so the starting error is second order in the step. If the predicted start has a
residual above a threshold and the plain warm start has a smaller one, the plain warm start is used instead. The
predictor is on by default; `tangent_predictor=false` gives the old behaviour.

**What it bought.** Average residuals R and Jacobians J per solve over 10 adaptive Tsit5 steps (Newton agent):

| Regime | Before (R / J) | After (R / J) |
|---|---|---|
| plateau | 7.8 / 1.31 | 2.4 / 1.00 |
| transition | 11.6 / 3.23 | 4.7 / 1.05 |
| nonlinear | 23.8 / 2.67 | about 6 / 1.3 |
| linear | 16.5 / 3.45 | about 5 / 1.55 |

**How we know the answer did not change.** The prediction error was measured to be second order, as expected. Newton
converges to the same solution from any good start, to within its tolerance. The old path remains available as a
reference, and the unit tests compare the two.

#### 2.3.3 What to do at the roundoff floor

**What was wrong.** Near convergence the residual cannot fall below a floor set by roundoff (2.1), so a step that does
not reduce it is not a failure. The old logic treated it as one and backtracked.

**What changed.** A new decision rule handles this case. Write $\mathrm{tol}$ for the Newton tolerance
(default $10^{-10}$). It applies to the scaled residual, in which each row is divided by the expected size of its
channel (section 3.2).

- If the residual is below $\mathrm{tol}$ and the step is below $\sqrt{\mathrm{tol}}$, the solve has converged.
- For a chord step this also requires the step to be below $\mathrm{tol}$. Otherwise the large step may come from the
  out-of-date matrix, so the Jacobian is refactorised and a Newton step decides.
- A Newton step below $\mathrm{tol}$ that does not reduce a residual above $\mathrm{tol}$ means the floor lies above
  the tolerance. Repeating the step from the same point would give the same result, so the solve stops as not
  converged instead of backtracking.
- The opt-in keyword `floor_tol` relaxes this for a known floor (section 3.2).

**What it bought.** Part of the same-point saving above, and at $\alpha = 0.187$ a solve at the floor stops after about
7 iterations instead of 30 (3.2).

**How we know the answer did not change.** The rule changes only when the iteration stops. The unit tests cover the
control flow of `solve_radial`, including `floor_tol`.

#### 2.3.4 The factorisation used for the rates (accuracy fix)

**What was wrong.** The amplitude rates need the derivatives of the solution with respect to the amplitudes. These
come from the implicit-function theorem: if $F(u; p) = 0$ then

$$\frac{du}{dp} = -J^{-1}\,\frac{\partial F}{\partial p},$$

with $J = \partial F/\partial u$ evaluated at the solution. The code computes this with the factorisation returned with
the solution. The formula is exact only if that factorisation was taken at the converged point. The old code accepted
a factorisation taken anywhere within $\sqrt{\mathrm{tol}} \approx 10^{-5}$ of the converged point, for example the
one from an earlier chord sequence.

**What changed.** The returned factorisation must now be taken within $\mathrm{tol}/10$ of the converged point.

**What it bought.** Correct rates. Measured: rate errors are 0.1 to 6 times the distance of the factorisation from the
converged point. Individual rate evaluations with the old rule were wrong by up to 1.1% (section 3.1).

**How we know.** See section 3.1 for what the check does and does not show.

Not implemented: a cost-aware refactorisation rule, which would decide when to refactorise by comparing the cost of
further chord steps with the cost of a new factorisation. It was estimated to cut long chord sequences (7 to 15 steps)
by about 25%.

#### 2.3.5 A correction to chord steps that was added and then removed

**What it was.** The Newton work also added a correction to the frequency part of a chord step, for the case where the
stored factorisation came from a solve at a different unknown scale. A later test showed it made a chord step worse in
a linear-regime case.

**What was done.** Whole evolutions were run with and without it:

| Run | Calls (both) | Wall time with | Wall time without |
|---|---|---|---|
| linear window | 241 | 63.5 s | 64.0 s |
| long window | 2233 | 294.6 s | 296.4 s |

Final states agreed within $10^{-7}$. Wall times agreed within 1%. The correction was removed.

**How we know the answer did not change.** Chord steps change only how fast Newton converges, not the converged
solution, and the A/B runs confirm it.

*For the cleanup pass: branch `speedup/newton`, commits `02fda5c`, `753d59a`, `c5e5125`, `e1fca6d`; polish `eabe89e`.
Same-point fix and floor logic in `_floor_action` (`src/bvp.jl`). Predictor: `solve_radial(...; guess_change)` with
`SolutionTangents`; threshold keyword `chord_threshold`; at merge the coordinator took the tangents from
`BalanceTerms.amplitude_directions` and `background_directions`. Final-factorisation rule: validation script
`validation/newton_jacobian_accuracy.jl`. The removed correction was `_rescale_frequency_step!`, removed in `9e17100`;
A/B jobs 27821573 to 27821576.*

### 2.4 Balance terms

**What this part of the calculation does.** `balance_terms` turns a radial solution into the rates of $\ln A_i$, $b$
and $j$. It computes the cloud's charges and fluxes (the observables) and their derivatives with respect to the
background ($b$, $j$) and to the $K$ amplitudes.

**Why it was slow.** It took 0.984 s. It ran the observables several times. It took the derivative of the charges
with respect to the background by central differences, $[Q(b+\delta) - Q(b-\delta)]/2\delta$ with relative step
$\delta = 10^{-6}$, which needs extra evaluations and has an error from both truncation and roundoff.

**What changed.**

1. One pass of the observables on dual numbers (2.1) with $2 + K$ partials ($b$, $j$ and the $K$ amplitudes) replaces
   the separate passes. The derivative with respect to the background is therefore exact.
2. The phase-averaged potential on dual numbers uses an analytic partial instead of differentiating through every
   operation.
3. The integrals over mesh nodes run on several threads, with a fixed summation order so that the result does not
   depend on the number of threads (2.1). Later they were made allocation-free.
4. The metric-work term $w$ is needed only by the energy balance and for diagnostics. Under the action balance it is
   now opt-in, off by default. When it is off, `w` is `NaN` in the diagnostics and in the spindown CSV. Requesting it
   also computes $E$ and $J$. After the one-pass change the opt-in saves little; the owner decided to keep it off by
   default, because the action balance has no metric-work term (section 8, item 1).

**What it bought.** On its branch alone, `balance_terms` went from 1.07 s to 0.51 s. The remaining cost was then the
derivative of the residual with respect to $b$ and $j$, which was handed to the assembly work (2.1, item 4). With
everything merged (harness, at `2246b4f`):

| | Baseline, 1 core | Final, 1 core | Final, 8 cores |
|---|---|---|---|
| `balance_terms` | 0.984 s | 0.067 s | 0.025 s |
| extra cost of requesting $w$ (measured before later speedups) | | 30 to 65 ms | |

**How we know the answer did not change.** The exact background derivative differs from the old central difference
by $10^{-11}$ to $2\times10^{-9}$ relative, which is the size of the old method's own error. The action-balance rates
agree with the old ones to $2.7\times10^{-13}$ relative. The action-balance rates and background rates are identical
bit for bit with the metric work on and off.

*For the cleanup pass: branch `speedup/balance`, commit `cea6095`; `src/evolution.jl`, `src/charges.jl`. One-pass
observables: `observable_partials`. Dual potential average: `averaged_potential!`. Threaded node integrals:
`integrated_densities`, allocation-free in `a9754a4`. Opt-in: keyword `record_metric_work`, script variable
`AXION_RECORD_METRIC_WORK=1`. The residual derivative was in `background_tangents`.*

### 2.5 Series at the outer boundary

**What this part of the calculation does.** At large $r$ each propagating component is matched to a series solution of
the radial equation, an expansion in $1/r$ whose coefficients are found order by order (up to order 40) from the
radial operator. An inward integration then carries it to the edge of the mesh. This supplies the outer boundary
condition and its derivatives.

**Why it was slow.** The series was built from the whole operator at every evaluation. The polish agent found that the
allocation at the outer boundary came from this order-40 series construction, about 2.7 MB and 1 ms per propagating
component, and not from the inward integration.

**What changed.** The series coefficients are built from single operator entries, and the inward integration no longer
allocates.

**What it bought.** The outer-boundary part of the Jacobian went from 62 ms to 8 ms (polish agent).

**How we know the answer did not change.** No separate check is recorded for this change. It is covered by the
regression test and the unit tests passing at the final commit (section 4).

*For the cleanup pass: polish commit `5f8d0f6`; `src/series.jl`, type `SeriesTerm`; outer `boundary_jacobian!`.*

### 2.6 Less allocation

**What this part of the calculation does.** This is not one part of the calculation but a sweep across all of them,
to remove the remaining allocation (2.1).

**Why it was slow.** Allocation costs time directly and through garbage collection. With several threads the garbage
collector pauses all of them while it runs. At 8 to 16 threads garbage collection took 16 to 20% of wall time.

**What changed.**

1. The Jacobian's storage is created only when a Jacobian is actually needed, and its sparsity pattern is kept
   separately from the number buffers.
2. The factorisation is built from the Jacobian's blocks without an intermediate sparse matrix (2.2).
3. The observable densities are allocation-free, and one back-substitution handles all $K$ amplitude columns at once.
4. The per-chunk workspaces are reused across all the residuals and Jacobians of one solve.

The same round also changed behaviour in ways that do not change any compared number. The line search now reports the
damping it actually used. A chord step below $\mathrm{tol}$ that does not halve the residual now triggers a
refactorisation instead of counting as stalled. After a rejected prediction (2.3.2) the Newton scaling (the row and
column weights of 3.2) is rebuilt from the plain warm start. And `solve_radial` was split into three smaller functions.

**What it bought.** Memory allocated per rhs call went from 6.6 GB at baseline to 56 MB on one core and 87 MB on 8
cores (`2246b4f`).

**How we know the answer did not change.** The regression test passes at one and eight cores, and the unit tests pass.

*For the cleanup pass: polish commits `9ec3265` (Jacobian storage on demand, `JacobianPattern`), `50ddeb9`
(factorisation from `JacobianBlocks`), `a9754a4` (observable densities, one back-substitution for $K$ columns),
`5f8d0f6` (series), `2246b4f` (per-chunk `NodeWorkspace`s reused within a solve), `eabe89e` (`solve_radial` split into
`_newton_start`, `_line_search`, `_convergence_status`).*

### 2.7 Support for stiff integrators

**What this part of the calculation does.** The time integrator advances $y = (\ln A_i, b, j)$ under
$dy/dt = g(y)$, where each evaluation of $g$ is one rhs evaluation. Until now the only integrator was Tsit5.

**Stiffness.** Linearise $g$ about a state: small deviations obey $d\,\delta y/dt = \Lambda\,\delta y$ with
$\Lambda = \partial g/\partial y$. A problem is **stiff** when $\Lambda$ has an eigenvalue of large negative real
part, a fast relaxation, next to slow dynamics. On the $|211\rangle$ plateau this is the case: the fast eigenvalue is
$-3.27\times10^{-8}$ (a relaxation time of about $3\times10^7$ in units of $M_0$), while the slow eigenvalues have
absolute value below $5\times10^{-13}$ (times longer than $2\times10^{12}$). The fast direction is presumably the
relaxation of $A_{211}$ back to its plateau value after a small displacement.

**Explicit versus implicit.** An **explicit** method such as Tsit5 builds $y_{n+1}$ from evaluations of $g$ at known
states. Applied to $\dot y = \lambda y$ it multiplies $y$ by a polynomial in $h\lambda$ each step, and that polynomial
exceeds 1 in magnitude once $h|\lambda|$ is past the method's stability limit, about 3.5 for Tsit5 on the negative real
axis. So on the plateau its step is capped at about $3.5/3.27\times10^{-8} \approx 10^8$, however slowly the cloud
itself is changing. (My arithmetic: Tsit5's 758 solves on the plateau segment of length $1.35\times10^{10}$, at about 6
new evaluations per step, give steps of about $1.1\times10^8$, consistent with that cap.) An **implicit** method such as
backward Euler, $y_{n+1} = y_n + h\,g(y_{n+1})$, stays stable on decaying modes for any $h$. The price is an equation
for $y_{n+1}$ to solve at every step, by Newton's method, which needs $\Lambda = \partial g/\partial y$. A
**Rosenbrock** method replaces that Newton iteration by one linear solve with the matrix $I - h\gamma\Lambda$ per stage
($\gamma$ is a coefficient of the method), so it is called linearly implicit. It also needs $\Lambda$.

Note that $\Lambda$ is a second, much smaller Jacobian: $(K+2)$ by $(K+2)$, so 4 by 4 here, not the 22112-by-22112 one
of the radial solve.

**What changed.**

1. A keyword `algorithm` on both evolution drivers selects the integrator. The default is Tsit5.
2. `FiniteDifferenceJacobian` supplies $\Lambda$ for implicit and Rosenbrock methods. Column $k$ is a difference
   quotient, $[g(y + \epsilon e_k) - g(y)]/\epsilon$ (forward) or the central version. Each column needs its own
   radial solve at a perturbed state, started from the solution at the base state. The columns do not depend on one
   another, so they run as concurrent tasks on different threads. A perturbed solve that fails is retried with a
   smaller step. A reuse policy exists, off by default, because reusing the Jacobian gave erratic errors away from the
   plateau.
3. A record of radial solves by purpose: stage solves, Jacobian solves and solves at accepted states. The last kind is
   needed because of a property called **first same as last** (FSAL). Tsit5's last stage of each step is evaluated at
   the accepted new state and reused as the first stage of the next step, so the radial solution at the accepted state
   is a by-product. A method without this property, such as Rodas5P, never evaluates at the accepted state, so the code
   needs a separate solve there.
4. New package dependencies for the Rosenbrock methods and, later, for FBDF (section 5).

**What it bought.** The ability to use a stiff integrator; the solve counts per method are in section 5. The parallel
forward-difference Jacobian took 9.8 s on one thread and about 5.2 s on 2 to 8 threads, on the old radial code. It was
limited by the slowest column, by garbage collection and by a lock inside UMFPACK solves. On the merged solver it
takes 0.16 to 0.23 s on 8 threads (section 5).

**How we know the answer did not change.** With the default Tsit5 the sequence of rhs calls is the same as before. The
independent review found no wrong-number defect on the Tsit5, Rosenbrock23 and Rodas5P paths (section 4.4). The unit
tests cover the finite-difference Jacobian (steps, columns, failures, reuse) and the stiff integrators. The script
options are in section 6.

*For the cleanup pass: branch `speedup/integrator`, commits `b7e315e` (dependency: OrdinaryDiffEqRosenbrock 2.7.5),
`fcd98c7`, `1a6f91e`. Solve record: `RadialSolveCounts`, as `result.solve_counts`. Jacobian settings: forward or
central differences, step sizes, retry on failed solves (`step_reductions`), reuse policy (`max_reuse_steps`).*

### 2.8 What runs on several cores, and what limits scaling

On several threads: residual and Jacobian assembly (over chunks of mesh nodes and intervals), the per-component
outer-boundary derivatives and the boundary-row task, the copying of Jacobian blocks into the factorisation's panels,
the node integrals in the observables, and the perturbed solves of the finite-difference Jacobian. On one thread: the
elimination stages of a factorisation, and the logic that decides between chord and Newton steps.

Work that runs on one thread takes the same time however many threads are available, so it sets a floor on the total.
Scaling stops at about 8 threads: the evolution segment took 2.95 s at 8 threads and 2.91 s at 16 (polish agent). The
floor is the factorisation (16 ms per factorisation, its elimination stages being sequential) and one boundary
component (8 to 9 ms), and garbage collection is 16 to 20% of wall time at 8 to 16 threads. The agent reports this as
the remaining serial floor; I did not find a measurement of how much each part would give if removed.

The matrix routines (LAPACK) run on top of a library called BLAS, which can itself use several threads. Two layers of
threading competing for the same cores usually slows both down. The script therefore sets BLAS to one thread when Julia
has more than one. BLAS with more than one thread inside the threaded assembly was not tested.

## 3. Accuracy findings about the original code

### 3.1 Out-of-date factorisation in the rates (errors up to about 1%)

The old code computed the amplitude rates with a factorisation that could be up to $10^{-5}$ away from the converged
solution, and individual rate evaluations were wrong by up to 1.1%. In your stored linear-growth window this left
$A_{322}$ at $t = 1.5\times10^9$ off by $7.6\times10^{-5}$ relative (and $A_{211}$ by $8.5\times10^{-8}$); every
later window restarted from that state starts with the same offset. Nothing needs to change to use the new code.
You are rerunning old results with the new code (section 8, item 4); the test that would show whether the offset
matters at late times is to compare a late window restarted from a state produced by the new code with the old one.

The mechanism is in 2.3.4: the rates use $du/dp = -J^{-1}\,\partial F/\partial p$ with the returned factorisation,
and the old code returned one accepted within $\sqrt{\mathrm{tol}}$ (about $10^{-5}$) of the converged point. Errors in
the sensitivities are 0.1 to 6 times that distance. The Newton work found this by measurement, and the integrator work
found the same defect independently: the 1.1% error was in $d\ln A_{322}/dt$ at the takeover state. The rule is now
"within $\mathrm{tol}/10$".

Linear-window check, from $t = 0$ to $1.5\times10^9$ at $f_a/M_{\mathrm{Pl}} = 0.0466873$, state at the end:

| Run | $A_{211}$ | $A_{322}$ | rhs calls |
|---|---|---|---|
| stored (original code, tol $10^{-7}$) | $2.807735530\times10^{-3}$ | $3.3961332\times10^{-9}$ | 717 |
| new code, default tolerances (job 27821573, code at `2246b4f`) | $2.807734656\times10^{-3}$ | $3.396389002\times10^{-9}$ | 241 |
| new code, reltol = abstol = $10^{-10}$, `newton_tol` $10^{-12}$ (job 27821633) | $2.807735769\times10^{-3}$ | $3.396390155\times10^{-9}$ | 697 |

The tight run is the reference. The new default run is within $4\times10^{-7}$ relative of it in $A_{211}$ and
$3.4\times10^{-7}$ in $A_{322}$. The stored $A_{322}$ is $7.6\times10^{-5}$ away (and $A_{211}$ $8.5\times10^{-8}$).
For the code at `9e17100`, the default run (job 27821574, final state in its log) gives $A_{211} = 2.807734301\times10^{-3}$
and $A_{322} = 3.396388873\times10^{-9}$. From those logged values I computed differences from the tight run of
$5.2\times10^{-7}$ and $3.8\times10^{-7}$. That calculation is mine, not STATUS.md's.

What this establishes: the stored $A_{322}$ is not the converged value of the current discrete equations, and the new
code is much closer. What it does not establish: that the out-of-date factorisation is the cause. The check is
consistent with it, but the old code was not separately run with the fixed rule. The stored run used tolerance
$10^{-7}$, and the rate errors are of the size that would produce this offset.

Implications:

- Existing runs carry errors of this order in $A_{322}$ during linear growth. Later windows (such as
  $1.5\times10^9 \to 10^{11}$) were restarted from the stored state, so they inherit $7.6\times10^{-5}$ in $A_{322}$
  as an initial-condition offset (coordinator's note).
- Whether that offset changes the late-time result is not tested. The end-to-end comparison (section 1) started from
  the stored state in both versions, so it does not test this.
- Rate errors up to 1.1% were seen in single evaluations. Their effect on a whole run depends on how often the
  out-of-date case occurred, which was not measured.

*For the cleanup pass: validation script `validation/newton_jacobian_accuracy.jl`; A/B and tight-run logs in
`review-work/rescale_ab/`.*

### 3.2 Newton residual floor at $\alpha = 0.187$

On the large mesh used at $\alpha = 0.187$ (1355 nodes), the scaled Newton residual cannot go below about
$5\times10^{-10}$, so the default tolerance of $10^{-10}$ is never met. With the old code every solve there ran the
maximum 30 iterations, about 90 s. To run at $\alpha = 0.187$ with the new code, pass `newton_tol=1e-8`, or keep
the default `newton_tol` and pass `newton_floor_tol=3e-9` or larger, to `evolve_with_spindown`. The first version of
this document suggested `newton_floor_tol=1e-9`; the parameter audit found that this is not enough, because after a
small change of state the residual there stops at about $1.7\times10^{-9}$, above the $5\times10^{-10}$ of the stored
state, and every perturbed solve then fails. With `newton_tol=1e-8` the solve takes 4 instead of 7 iterations and
gives the same rates (`docs/parameter-audit.md`, section 5.1).

Why the floor exists. Each bound mode's amplitude is imposed at the outer end of the mesh, where its profile has
decayed like $r^\nu e^{-\kappa r}$ and can be many orders of magnitude below its peak. To keep that row from being
ignored, the code divides each channel by its expected size along the mesh (its envelope) before solving, so the outer
rows get large weights. The Newton agent's diagnosis is that roundoff in the outer-boundary row of channel 2, of order
machine precision relative to the numbers in that row, is multiplied by that row's weight of $1.1\times10^4$, and this
sets the floor. It is a roundoff limit, not a problem with the equations. I did not find an experiment in the sources
that removes the floor by changing the weight, so the cause is the agent's diagnosis, not a tested fix.

What the code does now. A solve at such a floor stops after about 7 iterations instead of 30 and reports not converged
(2.3.3). The keyword `floor_tol` of `solve_radial` accepts a Newton step at a residual below `floor_tol` as converged.
The evolution drivers pass it on as `newton_floor_tol`, whose default equals `newton_tol`, so default behaviour is
unchanged. The benchmark harness uses a tolerance of $10^{-8}$ on that state instead.

A related finding: a warm solve on the linear state after a $10^{-5}$ shift in $b$ and $j$ diverged. The harness
perturbs by $10^{-6}$.

*For the cleanup pass: logic in `_floor_action` (`src/bvp.jl`); weights from `channel_envelopes`; `newton_floor_tol`
added in `eabe89e` as a keyword of `EvolutionCache`, `evolve_fixed_background` and `evolve_with_spindown`.*

## 4. How correctness was checked

### 4.1 Regression test

A **regression test** runs the code on fixed inputs and compares the outputs with values saved from a trusted version,
here the code before the project. It answers "did anything change?", not "is it right?". The question is how close
counts as unchanged. Each compared quantity has a **noise floor**: the size of the changes it shows under variations
that should not change it at all.

**Inputs.** Four stored states: `linear`, `plateau`, `nonlinear`, and `large_mesh` at $\alpha = 0.187$. Each is a
converged unknown vector saved to a file, not rebuilt, so every run starts from the same numbers.

**Quantities.** 168, in two tiers. The evaluation tier is pure evaluations at the stored vector or at a fixed smooth
displacement of it: residual, Jacobian applied to a vector, charges, losses, amplitude response, background
derivatives, `balance_terms`, rates. The re-solve tier runs `solve_radial` from a fixed perturbation (frequencies, 12
sampled unknowns, charges) and `evolve_with_spindown` from the plateau over $1.5\times10^9 \to 1.9\times10^9$ (saved
times, amplitudes, $b$, $j$).

**Tolerances.** Each quantity was compared with itself under three variations that must not change it: 4 BLAS threads
instead of one, a roundoff-sized perturbation of the input, and a tighter Newton tolerance. The tolerance is a factor of
about 3 to 10 above the largest difference seen, with a lower limit of $10^{-12}$. Different CPU models agreed to
$5\times10^{-16}$ or better. The tolerance classes are tabulated in `benchmarks/README.md`.

**What it establishes.** The optimised code gives the baseline's numbers, within noise, on these four states.

**What it does not establish.** Accuracy or convergence of the discretisation, because the references are the
baseline's own output. Other states, potentials, the energy balance, or step rejections (none occur in the short
segment). Only 12 entries of long vectors are compared, so a change confined to other entries is seen only through
norms and derived quantities. The `derivative` class tolerance ($3\times10^{-6}$) is loose because the baseline's
background derivatives contained a central difference; the exact derivative of 2.4 passes it. The trajectory
tolerances are of the order of the integrator's tolerance, so a flipped accept-or-reject decision shows up in the
printed step count.

**Result before the merge.** 168 of 168 pass at 1 and at 8 threads at `9e17100` and again at `057a167`, the last
commit before the merge. The saved references belong to the baseline commit; writing new references on an optimised
branch would defeat the test.

**The merge with the conjugate-frequency rule, and a new reference.** While the speedup was in progress, the owner
changed the equations on `rewrite/axioncloud`. The physics: the field is a sum of Fourier components, each labelled by
an integer vector $q$ and oscillating at the combination $\sum_i q_i \omega_i$ of the bound-mode frequencies. A bound
mode $i$ has complex frequency $\omega_i + i\gamma_i$, where $\gamma_i$ is its growth rate. A driven component such as
$q = (2,-1)$ is produced by the product $\phi_1^2\phi_2^*$ of the mode fields, whose size grows like $A_1^2 A_2$, that
is at the rate $2\gamma_1 + \gamma_2$. The old rule gave that component the complex frequency
$\sum_i q_i(\omega_i + i\gamma_i)$, with imaginary part $2\gamma_1 - \gamma_2$, so a component containing a
complex-conjugated mode ($q_i < 0$) was treated as partly decaying. The new rule gives it the imaginary part
$\sum_i |q_i|\gamma_i$, so every driven component grows with the modes that produce it. Only components with a
negative entry in $q$ are affected.

The merge brought this rule into the optimised code in two places. In `rhs_derivatives!`, the function that computes
the derivatives of the radial right-hand side, it merged as text without conflict. In `boundary_jacobian!`, which the
speedup had rewritten by hand (2.1, item 3), the chain rule for the frequency derivative was changed by hand to
match. Because the rule changes the equations, the old references no longer apply: against them the merged branch
fails 23 of 168 quantities, as it should. So the references were rewritten from the target's own original solver.
The target commit before the merge (`6bec64a`), which has the new rule and none of the speedup, was checked out
together with the harness, and the harness wrote new references there. The merged branch passes 168 of 168 against
those at 8 threads, and its 307 unit tests pass. This is the same test as before, against a baseline that solves the
new equations.

What this adds to the limits above: the stored states were converged under the old rule, so on the merged code they
are fixed points at which to evaluate, not solutions of the present equations. The re-solve tier and the evolution
segment converge to new solutions from them, and those are compared with what the original solver gives from the same
starting points. On the merged code the regression was run at 8 threads only.

*For the cleanup pass: `benchmarks/`, branch `speedup/base`, commit `b5b3d67`; states in `benchmarks/states/*.bin`;
references in `reference.toml`, first written at `b4c4100`, rewritten on 2026-11-14 from `6bec64a` (commit
`86c4b62`; procedure in `benchmarks/README.md`); do not use `--write-reference` on an optimised branch. Merge commit
`582b955`; reference job 27846901 (`pr-work/reference_6bec64a.toml` under `axion-sr-sim-speedup/`), comparison
27846902, unit tests 27846903. In `boundary_jacobian!` the port is `abs(weight)` on the imaginary-part column.*

### 4.2 Unit tests

A **unit test** is a small automated check of one function, usually against an independent calculation or against the
function it replaces. 307 tests pass on the merged code and 305 passed at `057a167` (both run with 4 threads); 278
passed at `9e17100`, and 279 at `2246b4f`, before the chord-step correction of 2.3.5 and its test were removed. There
were 119 before the project. The tests added after `9e17100` cover FBDF, the per-attempt failure record and the
solution at the accepted state (4.4), the recording of $E$ and $J$ (section 8, item 7), and the conjugate-frequency
rule (from the target).

New tests cover: field kernels against the complex formulas; the dual-number source against finite differences; the
structured Jacobian against automatic differentiation; in-place assembly and boundary derivatives; the block
factorisation against general sparse LU; warm starts and the predictor; one-pass balance terms; the finite-difference
Jacobian (steps, columns, failures, reuse); stiff integrators; the control flow of `solve_radial` and `floor_tol`; and
threaded paths. The threaded-path test starts the code in a separate process with 1 and with 4 threads and checks that
concurrent evaluations equal sequential ones and that results do not depend on the thread count. So the threaded paths
are exercised even if the suite itself runs on one thread.

Does not establish: physical correctness. Most tests compare an optimised function with its reference, so they check
that the two agree, not that either is right.

*For the cleanup pass: `test/runtests.jl`; jobs 27821630 (`9e17100`), 27824655 (`057a167`), 27846903 (merged);
threaded-path subprocess `test/threads_check.jl`.*

### 4.3 Independence of the thread count

Assembly results were identical bit for bit at 1, 2, 4, 8 and 16 threads (assembly agent). Concurrent evaluations
were identical bit for bit to sequential ones (independent review). The unit test above compares 1 and 4 threads bit
for bit. This holds because chunks are contiguous and every sum is taken in a fixed order (2.1). It does not cover BLAS
with more than one thread inside the threaded assembly, which was not tested.

### 4.4 Independent review

A separate model session reviewed the first full merge. Confirmed by experiment: no wrong-number defect on the Tsit5,
Rosenbrock23 and Rodas5P paths; the pivoting of the block LU; concurrent evaluations identical to sequential ones; the
dual-number source partials. Findings and what happened to them:

- A flag that marks a failed solve can stay set with integrators that run their own nonlinear solve (FBDF, SDIRK), and
  with fixed steps (`adaptive=false`) a failure freezes the evolution. **Fixed** in the merged code; see section 8,
  item 6.
- For integrators without the FSAL property (2.7), the radial solution stored as "accepted", the checkpoint's solution
  and the final radial solution are a stage, not the accepted state, contrary to the documentation. **Fixed** in the
  merged code, docstring included; see section 8, item 6.
- Test gaps (threaded paths, the dual-number source against an independent check, control flow, block LU shapes).
  Addressed in the second round.
- Small code items. Several were fixed in the second round; the rest are in section 9.

Establishes: the combined code on the three integrator paths has no wrong-number defect that the review's experiments
could expose. Does not establish: absence of defects on paths not run (other states or potentials). FBDF was not part
of the review; it was run in the integrator comparison (section 5) and has a unit test, but no independent review.

*For the cleanup pass: review of `01856b9`, scripts in `review-work/claude-review/`. The failure flag `solve_failed`
is now `failed_attempt`; the accepted-state solution comes from `accepted_radial_solution!`; both in `2a1089a`. Test
gaps addressed in `eabe89e`.*

### 4.5 End-to-end comparison

See section 1: same steps and the same trajectory within tolerance for one run ($f_a/M_{\mathrm{Pl}} = 0.0466873$,
$a_0 = 0.9$, $\alpha = 0.3$, cosine, window $1.5\times10^9 \to 10^{11}$), and the final state of the code at `9e17100`
agrees with the stored final row to $1.4\times10^{-8}$. Does not establish: other $f_a$, other $\alpha$, other potentials, or the merged code with the
conjugate-frequency rule, on which the window was not rerun. The original ran on different hardware.

## 5. Time integrators

**What was compared.** The integrator agent ran four integrators, Tsit5, Rosenbrock23, Rodas5P and FBDF, at
tolerances from $10^{-6}$ to $10^{-10}$, and measured for each run its cost and its error. This is a
**work-precision comparison**: for each method, how much work buys how much accuracy. Work is counted in radial
solves. Every right-hand-side evaluation, every column of the rates Jacobian $\Lambda$ (2.7) and every solve at an
accepted state is one radial solve of similar cost, so the count measures cost independently of how fast one solve
is. Error is measured against a reference run at a much tighter tolerance. All runs are at $\alpha = 0.3$,
$a_0 = 0.9$, cosine potential, action balance, 8 threads:

- the **main window**, $f_a = 5.7\times10^{17}$ GeV ($f_a/M_{\mathrm{Pl}} = 0.0467$), $t = 1.5\times10^9 \to 10^{11}$:
  the $|211\rangle$ plateau, the $|322\rangle$ takeover and the decay of $|211\rangle$;
- a **long plateau**, $f_a = 5.7\times10^{18}$ GeV ($f_a/M_{\mathrm{Pl}} = 0.467$), $t = 1.5\times10^9 \to 10^{12}$;
- a **small $f_a$**, $3.7\times10^{16}$ GeV ($f_a/M_{\mathrm{Pl}} = 3\times10^{-3}$), $t = 10^{10} \to 2\times10^{11}$.

The solve counts and errors were measured on the code at `2a1089a`, which has the radial solver from before the
second round of cleanups (2.6) plus FBDF support. They do not depend on how fast a solve is, so they carry over to the
merged code. Wall times were measured again on the merged solver for the main window (below).

**Tolerance requested versus error achieved.** The integrator applies the requested tolerance to its own estimate of
the error made in each step, and chooses the step size to keep that estimate small. The error at the end of a run
accumulates over many steps and depends on how good the estimate is, so it can be much larger than the request; the
stiff methods below show this. The only reliable measure is a comparison with a tighter reference run.

**How the error is measured.** At every saved output time of a run (checkpoint), the agent took the largest of
$|\Delta \ln A_{211}|$, $|\Delta \ln A_{322}|$, $|\Delta b|$ and $|\Delta j|$ against the reference, and reports the
largest over all checkpoints. A difference in $\ln A$ is a relative error of the amplitude. Taking the largest over
the run, not the error at the end, is the conservative choice: at the end of the main window $A_{211}$ is decaying
fast, so a small error in timing appears as a large relative amplitude error. For example, FBDF at $10^{-9}$ ends
with an error of $4\times10^{-6}$ in $\ln A_{211}$, which is a shift of about $4\times10^4$ in $t$ (units of $M_0$).
The reference for the main window is Tsit5 at $10^{-11}$. It agrees with FBDF at $10^{-10}$ to $8\times10^{-8}$ and
with Rodas5P at $10^{-10}$ to $1.5\times10^{-9}$, so it is good to about $10^{-9}$.

**Where the problem is stiff.** In the main window only on the $|211\rangle$ plateau, with fast eigenvalue
$-3.27\times10^{-8}$ against slow ones of absolute value below $5\times10^{-13}$ (see 2.7 for what this means and why
it caps Tsit5's step). Tsit5 takes about 800 solves on the plateau segment whatever the tolerance, because its step is
held at $\Delta t \approx 1.1\times10^8 \approx 3.5/|\lambda|$. The long-plateau run is stiff for the same reason; the
small-$f_a$ run is not. The agent's reading is that the fast rate is set by the ratio of cloud mass to black-hole mass
on the plateau. Direct eigenvalue measurements made after its report fit this: $-3.26\times10^{-8}$ on the plateau of
the long-plateau run, the same as in the main window (the plateau amplitudes of the two runs, in units of $f_a$,
differ by the factor of 10 that compensates the factor of 10 in $f_a$, so, by my reading, the cloud masses are
equal), and nothing
faster than about $1.4\times10^{-10}$ in absolute value at the two times sampled on the small-$f_a$ run. So "small
$f_a$ means long stiff plateaus" did not hold here.

**Main window.** Total radial solves / largest error over all checkpoints:

| tolerance | Tsit5 | Rosenbrock23 | Rodas5P | FBDF |
|---|---|---|---|---|
| $10^{-6}$ | 2168 / $1.5\times10^{-5}$ | 682 / $1.4\times10^{-3}$ | 775 / $3.5\times10^{-5}$ | 334 / $1.3\times10^{-3}$ |
| $10^{-7}$ | 2300 / $2.3\times10^{-7}$ | 1332 / $3.5\times10^{-4}$ | 1195 / $5.8\times10^{-7}$ | 466 / $1.2\times10^{-4}$ |
| $10^{-8}$ | 2498 / $5.9\times10^{-8}$ | 2692 / $8.6\times10^{-5}$ | 2969 / $1.2\times10^{-7}$ | 687 / $8.8\times10^{-6}$ |
| $10^{-9}$ | 3062 / $3.1\times10^{-9}$ | 5544 / $1.3\times10^{-5}$ | 8031 / $2.0\times10^{-8}$ | 881 / $4.1\times10^{-6}$ |
| $10^{-10}$ | 4022 / $1.5\times10^{-9}$ | 11922 / $1.5\times10^{-6}$ | 24809 / $1.5\times10^{-9}$ | 1144 / $7.7\times10^{-8}$ |

Compared at equal achieved accuracy:

- about $10^{-6}$: Rodas5P at $10^{-7}$ or FBDF at $10^{-9}$ cost about half of Tsit5;
- about $10^{-7}$: FBDF at $10^{-10}$ (1144 solves) costs about half of Tsit5 at $10^{-8}$ (2498) and 0.4 times
  Rodas5P at $10^{-8}$;
- about $10^{-9}$: Tsit5 is cheapest; Rodas5P needs about 25000 solves;
- Rosenbrock23 is never competitive.

By segment: on the plateau ($1.5\times10^9 \to 1.5\times10^{10}$) FBDF at $10^{-10}$ takes 85 solves with error
$2\times10^{-9}$, about 10 times fewer than Tsit5. Through the takeover and decay ($1.5\times10^{10} \to 10^{11}$)
FBDF at $10^{-10}$ takes 1058 solves (error $7.7\times10^{-8}$), Tsit5 1536 at $10^{-7}$ ($2.3\times10^{-7}$) and
1728 at $10^{-8}$ ($6\times10^{-8}$). FBDF's error does not follow its tolerance closely: at $10^{-8}$ and at
$10^{-9}$ it ends near $4\times10^{-6}$ both times.

**Long plateau** ($f_a = 5.7\times10^{18}$ GeV). The reference is Rodas5P at $10^{-11}$; it agrees with FBDF at
$10^{-11}$ to $2\times10^{-10}$.

| run | solves | largest error |
|---|---|---|
| FBDF $10^{-8}$ | 52 | $1.2\times10^{-7}$ |
| FBDF $10^{-9}$ | 77 | $4.5\times10^{-9}$ |
| FBDF $10^{-11}$ | 167 | $1.7\times10^{-10}$ |
| Rosenbrock23 $10^{-9}$ | 180 | $1.0\times10^{-7}$ |
| Rodas5P $10^{-9}$ | 777 | $1.2\times10^{-10}$ |

Tsit5 first failed on this run: the library's rule for the first step chose a step so large that a stage landed at
negative spin, and constructing the Kerr background there threw an error instead of the step being rejected. The
merged code now treats such a trial state as a failed step attempt, so the integrator retries with a smaller step
(`29d8a6f`). Reruns with an explicit first step reached $t = 10^{11}$, not $10^{12}$: 5516 solves at tolerance
$10^{-7}$ and 5546 at $10^{-9}$, against 57 for FBDF at $10^{-9}$ over the same interval, about 100 times fewer. Their
errors were not analysed. Over the whole run to $10^{12}$ the agent estimated, from the production log, about
$6\times10^4$ Tsit5 solves, which would make FBDF 500 to 800 times cheaper; that is an extrapolation.

**Small $f_a$** ($3.7\times10^{16}$ GeV), not stiff. The reference is Rodas5P at $10^{-11}$; it agrees with FBDF at
$10^{-11}$ to $5\times10^{-8}$.

| run | solves | largest error |
|---|---|---|
| Tsit5 $10^{-7}$ | 248 | $2.2\times10^{-8}$ |
| Tsit5 $10^{-9}$ | 500 | $6.7\times10^{-10}$ |
| FBDF $10^{-9}$ | 542 | $1.6\times10^{-6}$ |
| FBDF $10^{-11}$ | 756 | $5.2\times10^{-8}$ |
| Rodas5P $10^{-7}$ | 599 | $2.6\times10^{-6}$ |

Tsit5 is the best method here.

**Other candidates on the main window.** QNDF, a close relative of FBDF, behaves like it (781 solves, error
$7.8\times10^{-6}$ at $10^{-9}$). The library's automatic-switching composites start with Tsit5 and switch to a stiff
method when they detect stiffness: AutoTsit5(Rodas5P) behaves like Rodas5P and AutoTsit5(Rosenbrock23) like
Rosenbrock23, and AutoTsit5(FBDF) fails at $10^{-9}$ (3830 solves, error $1.6\times10^{-2}$). KenCarp4, an implicit
Runge-Kutta method, was only smoke-tested: 474 solves where FBDF needed 23.

**Wall time on the merged solver code.** The main window was timed again on `a84f712`, the merge of the speedup
(`9e17100`) into the integrator branch, which has the solver of section 2 and FBDF support but not yet the
conjugate-frequency rule. Eight threads, lr7, one run each. The last column is the relative difference of the final
$A_{211}$ from the Tsit5 run at $10^{-8}$; I recomputed it from the logged final states.

| run | wall time | radial solves | final $A_{211}$ against Tsit5 $10^{-8}$ |
|---|---|---|---|
| Tsit5 $10^{-7}$ | 266 s | 2300 | $2.2\times10^{-7}$ |
| Tsit5 $10^{-8}$ | 275 s | 2492 | reference |
| FBDF $10^{-10}$ | 119 s | 1036 | $4.5\times10^{-7}$ |
| FBDF $10^{-8}$ | 77 s | 595 | $1.0\times10^{-5}$ |
| Rodas5P $10^{-8}$ | 306 s | 3225 | $1.9\times10^{-8}$ |
| Rosenbrock23 $10^{-8}$ | 209 s | 2692 | $8.6\times10^{-5}$ |

One radial solve costs about 0.1 s for every method, so wall time follows the solve count. FBDF at $10^{-10}$ takes
less than half of Tsit5's time. Reading the last column: the final state lies in the fast decay of $|211\rangle$,
the least favourable point for this measure; and the Tsit5 run at $10^{-8}$ is itself off the converged answer by
about $6\times10^{-8}$ (main-window table), so differences below that are not resolved. The FBDF error at
$10^{-10}$, $4.5\times10^{-7}$ at the final time, is larger than the $7.7\times10^{-8}$ (largest over all
checkpoints) that the same setting gave on `2a1089a`, and the solve count changed from 1144 to 1036. The two codes
give the same rates within the regression tolerances, but differences at the level of the Newton tolerance change
FBDF's sequence of steps, and at fixed tolerance its error moved by a factor of a few. The other methods hardly moved:
Tsit5 at $10^{-7}$ and Rosenbrock23 took the same number of solves on both codes, Tsit5 at $10^{-8}$ 2492 instead of
2498, and Rodas5P 3225 instead of 2969.

**Vocabulary for the methods.**

- The **order** $p$ of a method: the error after a fixed time scales as $h^p$, so halving the step reduces it by
  $2^p$. High-order methods reach small errors in fewer steps.
- A **multistep** method builds $y_{n+1}$ from several previous states instead of from stages within one step. The
  backward differentiation formulas (**BDF**) are implicit multistep methods made for stiff problems.
- A **W-method** is a Rosenbrock method whose accuracy does not require the exact $\Lambda$; any reasonable
  approximation keeps the order. That matters here because $\Lambda$ comes from finite differences.
- **L-stable**: as $h\lambda \to -\infty$ the method's amplification factor goes to zero, so a fast relaxation is
  damped out in one step instead of persisting as a slowly decaying oscillation.

**The four methods, and why they performed as they did.**

- **Tsit5**: explicit, order 5. It needs no $\Lambda$, but on a stiff plateau its step is held at the stability
  limit, whatever the tolerance.
- **Rosenbrock23**: order 2, linearly implicit, a W-method, so a finite-difference $\Lambda$ does not spoil it. Order
  2 is too low: it needs many steps for a small error, and each step needs a new $\Lambda$ (4 perturbed solves).
- **Rodas5P**: order 5, linearly implicit, not a W-method. Its accuracy depends on $\Lambda$ being accurate. Ours is
  a finite-difference approximation, so at tight tolerances it loses order and takes many small steps. It also
  computes a new $\Lambda$ at every step (245 in the re-timed run).
- **FBDF**: variable-order BDF (orders 1 to 5). It solves its implicit equation for $y_{n+1}$ by Newton's method, in
  which $\Lambda$ only steers the iteration: an inaccurate or old $\Lambda$ costs a few extra iterations, not
  accuracy. FBDF asks for a new $\Lambda$ only when its Newton iteration converges slowly or the step size changes a
  lot, so it computed 3 to 6 per run (3 in the re-timed runs). Each costs 4 perturbed solves plus one solve at the
  accepted state, because FBDF never evaluates the right-hand side there itself.

Settings of $\Lambda$ for FBDF were varied on the main window at $10^{-8}$: central differences, other step sizes, and
reuse of $\Lambda$ over 3 steps. They took 579 to 724 solves (default 687) with errors $8.8\times10^{-6}$ to
$2.8\times10^{-5}$ (default $8.8\times10^{-6}$). Only 12 to 32 of the roughly 700 solves are for $\Lambda$. So the
defaults stay: forward differences, the default steps, and no reuse, since FBDF already reuses $\Lambda$ itself.

Cost of one $\Lambda$ on the merged solver (8 threads, plateau and later states of the three runs): 0.16 to 0.23 s
with the four perturbed solves running concurrently (central differences 0.29 to 0.39 s), against 0.49 s on the old
code. One perturbed solve takes 0.08 to 0.11 s, about the cost of an ordinary solve one step away.

**Conclusion.**

- **Neither Rosenbrock method.** Rosenbrock23 lost at every accuracy. Rodas5P is competitive only at about $10^{-6}$
  and is slower than Tsit5 at $10^{-8}$.
- **FBDF for stiff runs.** On runs dominated by long stiff plateaus (here $f_a \gtrsim 6\times10^{17}$ GeV at this
  $\alpha$ and spin), FBDF at tolerance $10^{-9}$. Through a $|322\rangle$ takeover, FBDF at $10^{-10}$; expect errors
  of a few times $10^{-7}$, at less than half of Tsit5's cost. FBDF's error at fixed tolerance moves by a factor of a
  few with small changes in the sequence of solves, so treat its tolerance as a rough setting and check a new kind of
  run once against a tighter one.
- **Tsit5 stays the default.** It is safe when the problem is not stiff and is then the best method (the
  $f_a = 3.7\times10^{16}$ GeV run). It is also the choice when errors well below $10^{-8}$ are wanted.
- Making FBDF the default of `scripts/evolution_spindown.jl` is the owner's choice; the script accepts
  `AXION_ALGORITHM=fbdf` (section 6). The parameter audit adds a point (section 8, item 2): the right-hand side itself
  carries discretisation errors of a few times $10^{-6}$ relative in the rates, so tolerances much tighter than $10^{-7}$
  buy agreement with the same discrete equations, not physical accuracy.

**Limits of the evidence.**

- One $\alpha$ (0.3) and one spin ($a_0 = 0.9$); three values of $f_a$; cosine potential only.
- The linear-growth window ($t = 0 \to 1.5\times10^9$) was never compared: its jobs were held in the queue and never
  ran.
- Solve counts and errors are from `2a1089a`; wall times from `a84f712`, one run each on shared nodes. Neither code
  has the conjugate-frequency rule. Nothing was re-run on the merged code.
- FBDF's error does not always fall as its tolerance is tightened, and it moved by a factor of a few between two
  nearly identical codes.
- The Tsit5 cost of the full long plateau is extrapolated; the Tsit5 reruns there were not compared with the
  reference.
- Switching method by regime (Tsit5 off the plateau, FBDF on it) was not run end to end. At accuracy of $10^{-7}$ or
  looser FBDF beats Tsit5 in every segment of the stiff runs, so a switch could pay only below $10^{-8}$. A possible
  criterion, $|\lambda_{\mathrm{fast}}|\,\Delta t$ from the $\Lambda$ that FBDF already computes, is untested.

*For the cleanup pass: branch `speedup/methods` (`07287c0`, `2a1089a`, `a84f712`), merged through
`integration-next`. Report: `methods-work/REPORT_2026-11-13_1115.md` under `axion-sr-sim-speedup/`; scripts
`run_case.jl` (it places `tstops` at the checkpoints, which is why its Tsit5 runs take 2300 solves rather than 2233),
`analyze.jl`, `plot_wp.jl`, `jacobian_cost.jl`; tables in `methods-work/data/analysis_*.txt`. Re-timing:
`methods-work/logs/retime_*` (jobs 27823322 to 27823328). Cost of $\Lambda$ and eigenvalues: jobs 27823319 to 27823321
(`jaccost-*`). Long-plateau Tsit5 reruns: 27823295, 27823296. Linear-window jobs `wp_linear_*` (27820490 to
27820508): held, never ran. Trial-state fix: `29d8a6f` (`is_kerr_black_hole` in `_spindown_rhs!`). FBDF dependency
OrdinaryDiffEqBDF 2.4.14.*

## 6. How to run it now

All commands from the repository root. Julia is `~/.juliaup/bin/julia` (1.13.0 on Lawrencium).

Threads. Use Julia threads for the kernels (`julia -t 8 --project=. ...`) and one BLAS thread
(`OPENBLAS_NUM_THREADS=1`). `scripts/evolution_spindown.jl` sets BLAS to one thread itself when Julia has more than one,
so that the concurrent solves of the finite-difference Jacobian do not compete for cores (2.8). Scaling stops near 8
threads. Slurm: account `pc_heptheory`, partition `lr7`, QOS `lr_normal`.

The script `scripts/evolution_spindown.jl <f_a/M_Pl> [initial spin] [final t/M0] [restart CSV]` reads these environment
variables:

| Variable | Effect |
|---|---|
| `AXION_ALGORITHM` | `tsit5` (default), `rosenbrock23`, `rodas5p` or `fbdf`. Stiff integrators get the finite-difference Jacobian; output files get the suffix `_<algorithm>`. For long stiff plateaus use `fbdf` with `AXION_TOLERANCE=1e-9`, through a $|322\rangle$ takeover `1e-10` (section 5). |
| `AXION_TOLERANCE` | Sets reltol = abstol. Default is the library's $10^{-7}$. |
| `AXION_CHECKPOINT_MINUTES` | Writes `<prefix>_checkpoint.jls` after the first accepted step at least that many minutes after the last. If the file exists, the same command resumes from it. A resumed run agrees with an uninterrupted run to the tolerance, not bit for bit. |
| `AXION_RECORD_METRIC_WORK` | `1` computes $E$, $J$ and the metric work for the diagnostics (`w` column); otherwise `w` is `NaN`. Rates are unchanged. |
| `AXION_CHARGES_EVERY_EVALUATION` | `1` writes $E$ and $J$ in every diagnostics row. Otherwise they are written for the first evaluation and the last evaluation of each accepted step and are `NaN` in the other rows (section 8, item 7). Rates are unchanged. |
| `AXION_POTENTIAL` | `cosine` (default), `chiral` or `free` (from the owner's branch, unchanged by this work). |

New keywords on `evolve_with_spindown` (and, where it applies, `evolve_fixed_background`; both also take
`EvolutionCache` options):

| Keyword | Meaning |
|---|---|
| `tangent_predictor=true` | Start solves from the first-order prediction (2.3.2). `false` is the reference path. |
| `record_metric_work=false` | Also compute $E$, $J$ and metric work under the action balance (diagnostic only). |
| `newton_floor_tol=newton_tol` | `floor_tol` of each radial solve; raise it only at a known residual floor (3.2). |
| `charges_every_evaluation=false` | `true` computes $E$ and $J$ at every evaluation; otherwise only at the first and the last evaluation of each accepted step (spindown driver only). |
| `algorithm=Tsit5()` | Integrator; `Rosenbrock23()`, `Rodas5P()`, `FBDF()` and `QNDF()` are re-exported from the library. |
| `jacobian_settings=FiniteDifferenceJacobian()` | The rates Jacobian $\Lambda$ for implicit and Rosenbrock integrators (2.7). |

For a callback of your own (for example one that writes checkpoints), `accepted_radial_solution!(integrator)` returns
the radial solution at the accepted state; it costs one extra radial solve for an integrator that is not
first-same-as-last (FBDF, Rodas5P). `result.final_radial_solution` is the solution at the final state for every
integrator.

Benchmark harness (Slurm, lr7):

```
benchmarks/submit.sh <cpus> <julia-threads> [blas-threads=1] [walltime]    # regression, then timing, one job
julia --project=. benchmarks/regression.jl                                 # compare to reference.toml, about 4 min
julia --project=. benchmarks/regression.jl --skip-evolution                # about 2 min
julia --project=. benchmarks/timing.jl --repeats 3                         # about 10 min
```

Output goes to `benchmarks/logs/`. Do not write new references on an optimised branch (4.1).

Tests:

```
OPENBLAS_NUM_THREADS=1 julia -t 4 --project=. -e 'using Pkg; Pkg.test()'    # about 6 min in a Slurm job
```

This is the command of the coordinator's test jobs. The threaded-path test starts its own 1-thread and 4-thread
processes whatever the thread count of the suite.

Manifest. `Manifest.toml` is now written by Julia 1.13.0 (`julia_version = "1.13.0"`); before, it was written by
1.11.7. The README still gives Perlmutter instructions with `julia/1.11.7`. The two dependency changes are in their own
commits (`b7e315e`, adding OrdinaryDiffEqRosenbrock 2.7.5 with 35 Manifest entries; `07287c0`, adding
OrdinaryDiffEqBDF 2.4.14 with 14 new Manifest entries and no change to existing ones), so they can be redone on the
other cluster by adding the packages there and regenerating the Manifest with that cluster's Julia. Using the 1.13
Manifest with Julia 1.11.7 was not tested.

## 7. Scaling to more bound modes

You want the method to handle of order ten bound modes eventually. $K$ is the number of bound modes (in the code, the
length of the `ModeSet`). The first version of this section estimated the costs from counts in the code. A separate
scaling study has since measured them and tested remedies; its document, `docs/mode-scaling-study.md`, has the detail.
This section summarises it. None of the growth described here was changed by the speedup.

**What grows with $K$.** The field is expanded in Fourier components $e^{-iq\cdot\xi}$, where $\xi$ is the vector of
the $K$ phases and $q$ a vector of $K$ integers; the code keeps every $q$ with odd $\sum_i |q_i| \le 3$ (one of each
pair $\pm q$). Three things grow:

1. *The phase grid*, exponentially. The self-interaction term is computed by evaluating the field on a grid of phase
   points and averaging. The present grid has 16 points per phase, $16^K$ in all: 256, 4096, 65536 for
   $K = 2, 3, 4$.
2. *The number of components* $C = 2K + 2K(K-1) + 4\binom{K}{3}$, about $\tfrac23 K^3$: 8, 22, 48, 680 for
   $K = 2, 3, 4, 10$. Each component brings its own radial profile, so the dense blocks of the Jacobian have side
   $4C$, and the block LU of 2.2 costs $(4C)^3$ per mesh node: $K^9$ in time and $K^6$ in memory.
3. *The mesh*, for some choices of mode. Adding $|433\rangle$ at $\alpha = 0.3$ brings a combination component whose
   frequency lies just below the mass threshold, so its profile decays very slowly and the mesh must be six times
   longer (4256 nodes instead of 690).

**Measured cost** (8 threads, 690-node mesh, lr7). One warm radial solve, the closest single number to one
right-hand-side evaluation:

| $K$ | modes | components | phase points | warm solve | memory of Jacobian and factorisation |
|---|---|---|---|---|---|
| 2 | $|211\rangle, |322\rangle$ | 8 | 256 | 0.10 s | 37 MB |
| 3 | adding $|311\rangle$ | 22 | 4096 | 1.5 s | 268 MB |
| 4 | adding $|311\rangle, |433\rangle$ | 48 | 65536 | about 47 s, from measured pieces | 1.26 GB |

On one thread the $K = 3$ solve takes 7.7 s, 25 times the $K = 2$ solve. At this cost the production window (2233
evaluations) would take about an hour at $K = 3$, and at $K = 4$ on the 4256-node mesh that $|433\rangle$ needs,
about a week (extrapolated). From $K = 3$ to $K = 4$ the measured pieces follow the cost model (work proportional to
components times phase points, and $(4C)^3$ for the factorisation) to within 10 to 25%. A $K = 3$ case with
$|433\rangle$ on its own long mesh did not converge: the first Newton iterations hit an exactly singular pivot, cause
not investigated.

**A phase rule that is not a tensor grid.** The phase average only has to deliver a polynomial number of Fourier
coefficients, so a much smaller point set will do. A **rank-1 lattice rule** with $M$ points uses the points
$\xi_j = 2\pi (j z \bmod M)/M$, $j = 0, \dots, M - 1$, for one integer vector $z$: a straight line on the torus of
phases that wraps around. It averages every Fourier harmonic $e^{ik\cdot\xi}$ exactly except those with
$k\cdot z \equiv 0 \pmod M$, so $z$ is chosen to keep those (the *aliases*) at high order, where the harmonics of the
nonlinearity are negligible. At the field amplitudes of the clouds, about 100 lattice points at $K = 2$ and 700 at
$K = 3$ reproduce the source and its Jacobian to $10^{-10}$, against 256 and 4096 tensor points. The number needed
grows like $K^6$, not exponentially. The rule is in the library behind the `phase_grid` keyword (tensor grid by
default), and the parameter audit tested it in full solves at $K = 2$ (section 8, item 2).

**Which components matter.** The components whose frequency is near $3\mu$ ($\sum_i q_i = 3$) carry less than
$10^{-12}$ of the rates. What matters are the fundamentals and the components that are the outgoing fourth leg of a
two-to-two process between populated modes, such as $(2,-1)$ and $(-1,2)$: at $K = 3$, 6 of the 22 components give
the rates to $10^{-4}$. Size alone is a poor guide: $(-1,2)$ has a profile $10^{-5}$ of the fundamentals but carries
5 to 7% of the $|322\rangle$ rate, because it radiates to infinity.

**Linear algebra.** A Newton step can be solved without forming the dense blocks, by an iterative (Krylov) method that
only needs products of the Jacobian with vectors, preconditioned by the part of the Jacobian that does not couple
components. It converged in 2 to 5 iterations at the clouds' amplitudes, at $K = 2$ and 3 alike. It removes the
$(4C)^3$ cost but needs more residual evaluations, so it pays only beyond about $K = 5$.

**Strong and weak modes, and the recommendation.** Cannizzaro, Spieksma and Witte (arXiv:2611.15865, PDF in
`docs/papers/`) find two dominant states at the production parameters and at most four or five dominant ones at late
times elsewhere, with many weakly populated states that barely feed back. The study's estimate is that $K = 10$
fully nonlinear costs hours per right-hand-side evaluation even with the lattice, the component selection and the
Krylov solver, because the residual alone grows like $C$ times the number of phase points, about $K^3 \times K^6$. Its
recommendation is therefore not ten fully nonlinear modes but **two to four strong modes treated as now, plus tens of
weak modes**. A weak mode is one whose amplitude is small enough for the field equations to be linearised in it. Each
weak mode is then a separate linear problem on the background of the strong cloud, with $1 + 2S^2$ components for $S$
strong modes; the weak modes are independent of one another, so the cost is linear in their number and they can run
in parallel. This is the primary/secondary construction of Omiya et al. that the old branch implemented. A weak mode
is promoted to strong when the terms of second order in its amplitude reach the wanted accuracy. Suggested order: the
lattice in full solves (now done at $K = 2$ by the audit); a small-amplitude check of the $(2,-1)$ flux against the
paper's rates (needs their data release); the weak-mode solver; dropping the $\sum_i q_i = 3$ components; a $K = 3$
strong evolution; a Krylov solver only beyond about four strong modes.

What this does not establish: no accuracy or convergence study was done at $K = 3$; the $K = 4$ cost is from pieces,
not solves; the lattice was tested at single mesh nodes for $K = 3$; the component ranking was measured at small
field amplitudes and not near resonances; and the strong/weak scheme is a design, not code.

*For the cleanup pass: branch `speedup/scaling` (`e90db99`), library change `49ad331` (`PhaseGrid(points, mode_set)`,
`rank1_lattice`, `phase_grid` keyword; cherry-picked onto the audit branch as `b3bcdb4`), scratch
`axion-sr-sim-speedup/scaling-work/`. `ModeSet` still enumerates all $7^K$ integer vectors with entries from $-3$ to
$3$ and filters them (`_integer_vectors` in `src/modes.jl`), which is $2.8\times10^8$ candidates at $K = 10$;
generating only the vectors that pass the filter removes this.*

## 8. Decisions and their outcomes

The first version of this document ended with ten open items. You answered them on 2026-11-13; this section gives
each item, the decision, and what came of it.

**1. Default of the metric-work diagnostic: stays off.** You asked whether the metric work must be included in the
spindown and therefore recorded by default. It need not. Under the default action balance the rates come from the
balance of the wave actions $N_i$ of the modes. A wave action is an adiabatic invariant: when the background changes
slowly compared with the oscillation, the action does not change, and the work that the changing metric does on the
cloud shows up as a change of the cloud's frequency and energy at fixed action, not as a separate term. So the action
balance has no metric-work term, and `w` is only a diagnostic (the coordinator's answer). `record_metric_work` stays
off by default, and `w` is `NaN` in the output unless you set `AXION_RECORD_METRIC_WORK=1`. The energy balance
(`balance=:energy`) still uses the metric work when `include_metric_work=true`.

**2. Numerical-parameter audit: done.** The audit is `docs/parameter-audit.md`. For each numerical parameter it
measured how much the four rates the integrator uses ($d\ln A_{211}/dt$, $d\ln A_{322}/dt$, $db/dt$, $dj/dt$)
change when the parameter is made finer, at seven states from linear to strongly nonlinear, and what one evaluation
then costs. The two strongly nonlinear states come from a run at $f_a = 3.7\times10^{16}$ GeV, where $A_{211}$
reaches about $2\times10^{-2}$, seven times the plateau of the production run. The error measure is the change of the
amplitude rates divided by the faster of the two rates. Each table is self-convergence: it shows how far production is
from a finer discretisation of the same equations, not that the equations are right. Summary of its recommendations
(cost change is per evaluation, at 1 / 8 threads):

| Parameter | Production | Error at production, worst state | Recommendation | Cost change |
|---|---|---|---|---|
| phase quadrature | tensor grid, 16 points per phase | $3\times10^{-15}$ | 85-point rank-1 lattice (section 7) if the new code path is accepted, otherwise 12 points per phase | lattice $-42\%$ / $-19\%$; 12 points $-26\%$ / $-15\%$ |
| `n_polar` (polar quadrature nodes) | 12 | $2\times10^{-15}$ | 8 | $-22\%$ / $0\%$ |
| **lattice and `n_polar` 8 together (measured)** | | $1.8\times10^{-13}$ | **the recommended quadrature** | $-49\%$ / $-20\%$ |
| radial `step` | 0.4 | $2.5\times10^{-6}$ | keep | |
| cap on the stretched outer intervals (`max_step`) | 1.86 | about $3\times10^{-6}$ | keep, unless `n_angular` is raised; 0.8 removes this error | $+45\%$ at 8 threads |
| outer endpoint (`containment`) | $10^{-4}$ | about $4\times10^{-6}$; $2\times10^{-5}$ at the strong takeover on the merged code | keep | |
| inner endpoint (`horizon_offset_ratio`) | $10^{-6}$ | $\le 6\times10^{-8}$ on the merged code | keep | |
| `stretch_factor` | 1.03 | $3\times10^{-7}$ | keep | |
| `n_angular` (polar indices per component) | 1 | $2\times10^{-2}$ at the strong states, $7\times10^{-4}$ on the plateau and at the takeover | the largest error found: use 2 for strongly nonlinear runs, after checking 3 | $+91\%$ / $+148\%$ |
| `max_order` (Fourier truncation) | 3 | $\lesssim 3\times10^{-6}$ | keep | |
| `newton_tol` | $10^{-10}$ | $\le 2\times10^{-9}$ even at $10^{-7}$ | keep at $\alpha = 0.3$; $10^{-8}$ on the $\alpha = 0.187$ mesh | up to $-28\%$ at $\alpha = 0.187$ |
| integrator tolerance | $10^{-7}$ | | keep $10^{-7}$ with Tsit5, or FBDF at $10^{-9}$ to $10^{-10}$ (section 5) | |

Main findings:

- *The angular truncation is the largest error.* Each Fourier component carries one polar index by default. Allowing
  a second one changes $d\ln A_{322}/dt$ by 2% at the strong $|211\rangle$ peak and by $7\times10^{-4}$ on the
  plateau of the production run: $10^3$ to $10^4$ times the radial-step error. Whether two polar indices are enough
  is not known; three were not run.
- *The quadratures can be made cheaper at no cost in accuracy*: the recommended lattice with 8 polar nodes agrees with a
  much finer quadrature to $1.8\times10^{-13}$ and halves the one-thread cost.
- *On the weaker states the drift with the outer endpoint is a resolution error of the coarse outer intervals*: a cap
  of 0.8 instead of 1.86 reduces it about 30-fold, as expected for a fourth-order scheme.
- *The integrator tolerance is already tighter than the discretisation*: the radial step alone gives relative rate
  errors of a few times $10^{-6}$, so tightening below $10^{-7}$ does not make results more accurate.
- *At $\alpha = 0.187$*, `newton_floor_tol=1e-9` is not enough; use `newton_tol=1e-8` or a floor tolerance of at least
  $3\times10^{-9}$ (section 3.2).

The audit's tables were computed on the code before the merge (`ebf0ae9`), without the conjugate-frequency rule
(4.1). Its Section 7 reran the checks that involve the affected components on the merged code, at the strong states,
and this changed one conclusion. On the old code the rates at the strong states moved by $6.3\times10^{-4}$ when the
outer endpoint was moved, and by $2.3\times10^{-5}$ for the inner endpoint, together with the amount of radiation
(the components $(-1,2)$ and $(2,-1)$ travelling out or into the hole) held inside the domain. That looked like a
question of whether radiation in transit should count in the cloud's actions, and the first report to you asked you
to decide it. On the merged code the domain holds the same radiation, but the rates move by only $2.2\times10^{-5}$
and $6.3\times10^{-8}$: with the new rule those components' growth is consistent with the action balance, so the
large endpoint dependence was mostly a symptom of the old frequency rule, and that question is withdrawn. The polar
truncation error is unchanged on the merged code (1.5 to 2%). The rule itself moves the strong-state rates by up to
$1.2\times10^{-3}$. The weaker-state numbers were not rerun on the merged code; the audit expects them to be
unchanged, since the driven components are tiny there, but did not measure it.

**3. Integrator: FBDF for stiff runs, Tsit5 otherwise.** Section 5 has the comparison and the conclusion: neither
Rosenbrock method; FBDF at $10^{-9}$ on long stiff plateaus and $10^{-10}$ through a $|322\rangle$ takeover; Tsit5
stays the default and is best where the problem is not stiff or very small errors are wanted.

**4. Existing results: you are rerunning everything with the new code.** This removes the $7.6\times10^{-5}$ offset of
3.1 from new results. The end-to-end window has not been rerun on the merged code by this project, so your reruns are
also the first full runs with the conjugate-frequency rule and the new code together.

**5. Remaining serial floor: not worked on.** Beyond about 8 threads the speedup stops, because the elimination
stages of the factorisation, one outer-boundary component and garbage collection do not use more threads (2.8). You
asked what removing them would gain. The coordinator's estimate, not a measurement, is a factor of 1.3 to 1.6 on 8
threads and nothing on 1 thread, since on one thread there is nothing to parallelise. It was not done, and neither
was the cost-aware refactorisation rule of 2.3.4 (an estimated 25% fewer steps in long chord sequences).

**6. The two review findings: fixed.** Both concern how the evolution keeps track of radial solves when the
integrator is not Tsit5 (4.4).

- *Failure record per step attempt.* When a radial solve fails inside a step, the step must be rejected and retried
  with a smaller step. The old code set a flag on failure, and the check that the integrator makes at the end of each
  step attempt rejected the step and cleared the flag. An integrator that runs its own Newton iteration, such as
  FBDF, can abandon a step attempt by itself, and then skips that check, so the flag stayed set and spoiled the next
  attempt. With fixed steps (`adaptive=false`) the check is not made at all, so after one failure every later
  evaluation returned zero rates and the evolution froze. Now the integrator's
  own counter of step attempts, which increases at every attempt whether it is accepted or not, is recorded at the
  failure. The rest of that attempt returns zero rates and is rejected, and the next attempt has a new number, so
  there is nothing to clear. With fixed steps a failed solve now stops the run with an error. This was tested with a
  stand-in integrator and by reading the integrator library's source, not with a real integrator whose own Newton
  iteration fails.
- *Radial solution at the accepted state.* Checkpoints and the final result need the radial solution at the state
  the integrator accepted. For an integrator without the first-same-as-last property (2.7) the stored solution was
  that of a nearby stage: for Rodas5P its last stage, for FBDF its last Newton iterate, since FBDF never evaluates the
  right-hand side at the accepted state. A new function, `accepted_radial_solution!`, returns the solution at the
  accepted state, doing one extra radial solve when needed; the checkpoint writer of the script uses it, and the final
  radial solution of a run is now taken at the final state for every integrator. The docstrings were corrected.

**7. $E$ and $J$ under the action balance: recorded at each step, not each evaluation.** You wanted $E$ and $J$
recorded at each time step, not necessarily at every evaluation. They are now computed for the first evaluation and
for the last evaluation of each accepted step, and are `NaN` in the other rows of the diagnostics. For Tsit5 the last
evaluation of a step is at the accepted state itself; for FBDF and Rodas5P it is a nearby state. Setting
`charges_every_evaluation=true` (script: `AXION_CHARGES_EVERY_EVALUATION=1`) restores a value in every row. Computing
$E$ and $J$ took 19 ms of a 227 ms evaluation on 1 thread and 3 ms of 82 ms on 8. The rates do not use them and are
unchanged. Plotting code that reads those columns has to skip `NaN`. The fixed-background driver is unchanged.

**8. The unexplained $2.1\times10^{-8}$ difference: explained.** The finding was: on the nonlinear state, with fixed
background and fixed Tsit5 step, the old and the new code differed by $2.1\times10^{-8}$ in $\ln A_{211}$ after 6
steps, while each agreed with its own run at tighter Newton tolerance to $3\times10^{-11}$. The suspected cause was the
out-of-date factorisation of 3.1. An agent tested this by patching the old code to refactorise at the converged point.
In this run the old code already did so in every solve, and the patched code gave identical results bit for bit. So
the out-of-date factorisation is not the cause.

The cause is a property of the radial problem at that state. At two of the 37 right-hand-side evaluations (stages 5
and 6 of step 5) the Jacobian is nearly singular in one direction $d$ of the space of unknowns: moving the unknowns
along $d$ changes the residual very little, $|Jd|/|d| \approx 2\times10^{-3}$. The direction involves the horizon
coefficient of one driven channel, two outer coefficients and the unknowns at nodes near the horizon. Along such a
direction Newton's method cannot locate the solution better than the residual's roundoff floor allows: two points
$4.9\times10^{-9}$ apart both have residuals at that floor ($6\times10^{-13}$ to $6\times10^{-12}$), so both count as
the solution. (My arithmetic: a floor of $6\times10^{-12}$ divided by $2\times10^{-3}$ is $3\times10^{-9}$, consistent
with that separation.) The two points give values of $d\ln A_{211}/dt$ that differ by $3.5\times10^{-8}$ relative.
Which of them a solve returns depends on the path its Newton iteration takes, in either code, and the old and new code
take different paths. A tighter tolerance does not help, since both points already satisfy it, which is why each code
agrees with its own tight run (my reading). Tsit5 combines the stages of a step with fixed weights, $-3.29$ and $+2.32$
for those two stages, which turns the $3.5\times10^{-8}$ into the $2.1\times10^{-8}$ seen after the step. So at that
state the rate is determined only to about $3\times10^{-8}$ relative, by the problem, not by either code. That is below
the integrator tolerance ($10^{-7}$) and far below the discretisation errors of item 2, so it does not limit results.

What this does not establish: how common such states are along a run (2 of 37 evaluations here); whether the nearly
singular direction is physical or a weakness of how the unknowns are parametrised; and anything about the 1.1% rate
errors of the out-of-date factorisation elsewhere (3.1), which this test does not touch. A separate finding of the
same test: on the linear state the old code differs from its own tight run by $5.7\times10^{-7}$ (fixed background)
and $3.1\times10^{-6}$ (spindown) in $\ln A_{211}$, the new code by $4\times10^{-8}$ and $3\times10^{-9}$.

**9. Merge: done.** Merged into `rewrite/axioncloud` on 2026-11-14 as pull request #1, after the branch took in the
conjugate-frequency rule and passed the regression against the target's own solver (4.1, section 10).

**10. Reproduction on Perlmutter with Julia 1.11.7: still open.** Nothing in the sources records a run of the new code
there (section 6).

*For the cleanup pass: item 1, `record_metric_work`, `_computes_charges`; item 2, `validation/parameter_audit*.jl`,
results in `axion-sr-sim-speedup/audit-work/run1` to `run8_merged`, spot check job 27847112; item 6, `failed_attempt`,
`_step_attempt`, `_fail_step_attempt!`, `_reject_failed_step`, `accepted_radial_solution!`, `_accepted_evaluation!`,
`_final_spindown_solution!` in `src/evolution.jl`, commit `2a1089a`; item 7, commit `057a167`, `_recorded_charges`,
`_fill_last_charges!`, default `false` in `evolve_with_spindown` and `true` on a bare `EvolutionCache`; item 8,
scripts and outputs in `axion-sr-sim-speedup/item8-work/`, STATUS.md paragraph "Item 8 resolved".*

## 9. Notes for the cleanup pass

Each name below was checked with grep, first against `9e17100` and again on 2026-11-14 against the merged code (as
contained in `speedup/followup`). Items marked "already removed" were listed by an agent but are gone.

Already removed (nothing to do): `DEFAULT_FACTORIZATION_METHOD` (global switch, removed in `50ddeb9`),
`_rescale_frequency_step!` (removed in `9e17100`), `JacobianWorkspace.layout` (field no longer exists),
`node_rhs_blocks`, allocating `interval_residual` / `interval_jacobian`, `IntervalJacobian` (removed by assembly). The
name `jacobian!` appears in `src/evolution.jl` only as a local variable holding the integrator's Jacobian callback; it
is not a library function. The assembly-level in-place functions are `jacobian_blocks!`, `boundary_jacobian!` and
similar.

Reference paths kept for validation (still present):

- [ ] `tangent_predictor=false` (reference path for the predictor; used in tests).
- [ ] `resolve_background=true` in `balance_terms`, with `_resolved_balance_pieces`; used by one test and three
      scripts in `validation/`.
- [ ] `central_difference_background_derivatives` (`src/evolution.jl`).
- [ ] `boyer_lindquist_metric_work` (`src/evolution.jl`; used by the energy balance and validation scripts).
- [ ] `method=:sparse_lu` in `ScaledFactorization` (UMFPACK reference; used by tests and a validation script) and
      `benchmarks/linear_solver.jl` (exists at `9e17100`; it compares the two methods).
- [ ] `observable_jacobian`, `charge_jacobian`, `action_jacobian` (still used by the fixed-background
      `amplitude_rates(solution)`).
- [ ] `phase_averaged_potential` (non-dual fallback, `src/field.jl`) next to `averaged_potential!` in `src/charges.jl`;
      the allocating `averaged_potential` is referenced only by a test. The review suggests merging `averaged_potential`
      into `phase_averaged_potential`; `potential_slope` (`src/charges.jl`) belongs in `potential.jl`.
- [ ] Old allocating wrappers kept by assembly with unchanged signatures (for example `rhs_block`, `rhs_derivatives`,
      `radial_rhs`).

Functions that are now test-only or unused (from the review's list; each verified):

- [ ] `rhs_block` (`src/bvp.jl`): called only by `rhs_derivatives_autodiff`.
- [ ] `rhs_derivatives_autodiff` and `RhsDerivatives` (`src/bvp.jl`): used only by tests and the `rhs_derivatives`
      wrapper.
- [ ] `rhs_derivatives` (allocating wrapper): used by one test and one scratch script in `validation/logs/`.
- [ ] allocating `radial_rhs` (`src/radial.jl`): no caller found in `src`, `test`, `scripts` or `benchmarks`.
- [ ] `total_solves` (`src/evolution.jl`): defined and exported; used only by one test (added with the per-attempt
      failure record).
- [ ] `solution_tangents` outside tests (the review listed this as an unused export; it is in fact used internally,
      four references in `src`).

Duplicated formulas and code:

- [ ] `delta_B` is written out in `radial_rhs!` (`src/radial.jl:304`) and in `rhs_derivatives!` (`src/bvp.jl:475`), and
      a third time in `validation/nonlinear_boundary_terms.jl`.
- [ ] `_report_step` and `_report_line` (`src/evolution.jl`): the review proposed merging them.

Structure (deferred by the review to the cleanup):

- [ ] Split `EvolutionCache` into settings and state.
- [ ] Replace `NaN` "not computed" markers (for `w`, charges) with an explicit representation.
- [ ] Move test-only functions out of `src/`.
- [ ] Two pre-existing weak spots in the chord and line-search logic that the review asked polish to assess; the
      sources do not record the outcome.
- [ ] Docstring corrections: the review listed some; its list is not in the sources. The one known case, the
      `evolve_with_spindown` docstring on the accepted-state solution, was corrected in `2a1089a` (section 8, item 6).
- [ ] `EvolutionCache` has a field `integrator::Any` (the running integrator, read by `_step_attempt`). Julia compiles
      fast code when it knows the types of the values it works on; a field declared `Any` hides the type, so every
      access is resolved while the program runs. The cost is negligible here (one access per evaluation), but the
      field is untyped state on a struct that already mixes settings and state; it belongs in the split above.
- [ ] `u_modified!(integrator, false)` in `_record_accepted_solution!` (`src/evolution.jl`) and in the checkpoint
      callback of `scripts/evolution_spindown.jl` prints a deprecation warning in the test log: the integrator library
      has renamed it `derivative_discontinuity!(integrator, false)`. The call should be renamed before the old name is
      removed.
- [ ] `JacobianStore`'s reuse policy (`max_reuse_steps`, `_reuse_stored_jacobian`) is unused in practice. It is off by
      default because reuse gave erratic errors with the Rosenbrock methods (2.7), and FBDF, the recommended stiff
      integrator, already decides itself when to ask for a new $\Lambda$, so the store only serves a repeated request
      at the same state (section 5). Keep the store, consider removing the policy.

Validation scripts and scratch:

- [ ] `validation/newton_counts.jl` is stale (still so on the merged code): its wrappers hook `jacobian` and the
      sparse `ScaledFactorization`, which `solve_radial` no longer calls, so it no longer counts Jacobians and
      factorisations. Fix or retire.
- [ ] `benchmarks/` reference provenance. `reference.toml` was written from `b4c4100` and rewritten on 2026-11-14
      from `6bec64a`, the original solver with the conjugate-frequency rule (4.1); the README states this and the
      procedure for the next change of equations. Three things to keep in mind: the stored states in `states/*.bin`
      were converged under the earlier rule, so they are fixed evaluation points, not solutions of the present
      equations; the README's baseline timing table is still from `b4c4100`, which is fine for timing because the
      rule does not change the work; and if `6bec64a` and the unoptimised solver are ever gone, a future reference
      has to be written from the last commit before the change plus a separate argument for the change, as the
      README says.
- [ ] Scratch directories and worktrees under `axion-sr-sim-speedup/` that can be deleted once nothing is needed from
      them: worktrees `base/`, `assembly/`, `linear-solver/`, `newton/`, `balance/`, `integrator/`, `integration/`,
      `polish/`, `methods/`, `integration-next/` (merged), `pr-baseline/` (detached at `6bec64a`, used to write the
      new reference), `docs/` (superseded), and `audit/`, `scaling/`, `followup/` once their documents are merged;
      and the scratch directories `balance-scratch/`, `integrator-work/`, `linear-solver-work/`, `methods-work/`,
      `polish-work/`, `review-work/`, `audit-work/`, `scaling-work/`, `item8-work/`, `pr-work/`. Keep `review-work/`
      until the end-to-end outputs (`e2e_01856b9/`) and the A/B logs (`rescale_ab/`) cited here are no longer needed,
      and `methods-work/`, `audit-work/`, `pr-work/` while sections 5 and 8 cite their logs. The worktrees
      `chiral-runs/` and `conjugate-rule/` belong to other sessions of the owner's work, not to this project.
      Untracked scratch inside worktrees (`.scratch/` in `assembly/`, `newton/work/`) goes with the worktrees.
      Remove worktrees with `git worktree remove`, not `rm`.
- [ ] `benchmarks/logs/` is git-ignored, so harness logs are not in the repository; copy any you want to keep.

## 10. Branches and commits

The speedup was merged into `rewrite/axioncloud` on 2026-11-14 as pull request #1
(https://github.com/clwelch03/axion-sr-sim/pull/1), merge commit `6c0a9e6`, whose second parent is
`speedup/integration-next` at `86c4b62`. `speedup/integration-next` was pushed for the pull request; the other
branches below were not. All speedup branches start from `b4c4100`.

| Branch | Head | Contents |
|---|---|---|
| `speedup/base` | `b5b3d67` | Benchmark and regression harness (`benchmarks/`), four stored states, tolerances |
| `speedup/assembly` | `83e8bb3` | In-place kernels, threading over mesh nodes, per-component boundary Jacobian, dual-number source with analytic partials |
| `speedup/linear-solver` | `cf853b6` | `src/block_factorization.jl`, block elimination behind `ScaledFactorization` |
| `speedup/newton` | `e1fca6d` | Predictor, same-point cost, roundoff-floor logic, tightened final factorisation, `floor_tol`, tests and validation |
| `speedup/balance` | `cea6095` | One-pass `balance_terms`, exact background derivative, opt-in metric work, threaded node integrals |
| `speedup/integrator` | `1a6f91e` | `algorithm` keyword, `FiniteDifferenceJacobian`, `RadialSolveCounts`, script options, Rosenbrock dependency |
| `speedup/integration` | `229b09f` | Coordinator's merge branch of the above; superseded by `integration-next` |
| `speedup/polish` | `2246b4f` | Second round: workspaces, block LU from blocks, series, allocation, `solve_radial` split, `newton_floor_tol` in drivers, tests |
| `speedup/methods` | `a84f712` | OrdinaryDiffEqBDF dependency (`07287c0`); FBDF and QNDF, per-attempt failure record, solutions at accepted and final states (`2a1089a`); merge of `9e17100` (`a84f712`). Fast-forwarded into `integration-next` |
| `speedup/integration-next` | `86c4b62` | Everything above, then: trial states that are not a Kerr black hole fail the step attempt (`29d8a6f`); harness follows the renamed failure record (`ebf0ae9`); $E$ and $J$ once per accepted step (`057a167`); merge of `rewrite/axioncloud` at `6bec64a` (`582b955`); regression reference rewritten from `6bec64a` (`86c4b62`). Merged as PR #1 |
| `speedup/audit` | `7896c56` | Numerical-parameter audit (`docs/parameter-audit.md`, `validation/parameter_audit*.jl`), the lattice phase rule cherry-picked from `speedup/scaling` (`b3bcdb4`), merge of `6c0a9e6` (`9ec31f6`) |
| `speedup/scaling` | `e90db99` | Scaling study with the number of bound modes (`docs/mode-scaling-study.md`) and the library change `49ad331` (`PhaseGrid` from a point list, `rank1_lattice`, `phase_grid` keyword) |
| `speedup/docs` | `af45ebf` | Earlier drafts of this document; superseded by `speedup/followup` |
| `speedup/followup` | this document's commit | `speedup/audit` plus the scaling study and this document; candidate for a second pull request |

The 26 commits of `git log --oneline b4c4100..9e17100` are: the harness (`b5b3d67`); the dependency commit
(`b7e315e`); the four first-round merges and their commits (`80b496f`, `83e8bb3`, `05feccb`, `cf853b6`, `cea6095`,
`7b6cafd`, `02fda5c`, `753d59a`, `c5e5125`, `e1fca6d`, `fcd98c7`, `1a6f91e`); merge commits (`9950e6d`, `35f6928`,
`efd6772`, `229b09f`, `01856b9`); and the polish commits (`9ec3265`, `50ddeb9`, `a9754a4`, `5f8d0f6`, `eabe89e`,
`2246b4f`) followed by `9e17100`. From `9e17100` to the merged head `86c4b62` the branch adds `07287c0`, `2a1089a`,
`a84f712`, `29d8a6f`, `ebf0ae9`, `057a167`, `582b955` and `86c4b62`, plus the target's own commits brought in by the
merge.

## Points that the sources leave open

- **Which code the "new default run" in the linear-window table is.** STATUS.md's table and its quoted distances
  ($4\times10^{-7}$, $3.4\times10^{-7}$) are the run with the chord-step correction (job 27821573, code at `2246b4f`).
  The code without it is job 27821574; I computed $5.2\times10^{-7}$ and $3.8\times10^{-7}$ from its logged final
  state. The conclusion does not change.
- **Which code the 295 s end-to-end figure is.** STATUS.md calls it "9e17100-equivalent". The two runs are 294.6 s
  (with the correction, `2246b4f`) and 296.4 s (without, `9e17100`).
- **"Serial block factorisation".** STATUS.md calls the factorisation serial (16 ms floor). In the code only the
  copying of blocks into the stage panels is threaded; I read the elimination stages as serial from the code but did
  not time them separately.
- **The review's docstring list and the outcome of the two chord and line-search weak spots** are not recorded in the
  sources.

Settled since the first version: the harness logs at `9e17100` state that the reference was written at `b4c4100`;
the one test removed between `2246b4f` (279 tests) and `9e17100` (278) is the chord-step correction's test; the
8-core spread is resolved by later runs, which all gave 0.068 to 0.080 s (section 1).
