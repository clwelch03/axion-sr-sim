# Performance project, November 2026

This is the record of the speedup of the AxionCloud time evolution (`evolve_with_spindown`, `evolve_fixed_background`),
done on 2026-11-12 and 2026-11-13. It is written for the code's owner and for whoever later cleans the code up. It
says what changed, why, how each change was checked, and what is left.

All of the work is on one branch. Nothing has been merged into `rewrite/axioncloud` or pushed. Two pieces are still in
progress and are marked as such: the comparison of time integrators on the final code, and the audit of numerical
parameters.

The work was done by several Claude agent sessions, each named after its task (assembly, linear solver, Newton,
balance, integrator, and "polish" for the second round of cleanups), coordinated by one more session. Measurements are
attributed to the agent that made them, because they were taken on different nodes and at different stages.

How to read section 2. Each subsection describes one change in five parts: what this part of the calculation does;
why it was slow (or, for an accuracy item, what was wrong); what changed; what it bought, in numbers; and how we know
the answer did not change. Each computing term is explained where it first appears and then used. A closing line in
italics, "For the cleanup pass", carries the branch names, commit hashes, job IDs and function names; you can skip
those lines on a first reading.

*For the cleanup pass: the branch is `speedup/integration-next` at commit `9e17100` (26 commits on top of
`rewrite/axioncloud` at `b4c4100`). In progress: branch `speedup/methods` (integrator comparison) and the
numerical-parameter audit. The coordinator's running record is `STATUS.md` in `axion-sr-sim-speedup/` (outside the
repository); the rules given to every agent are in `AGENT_BRIEF.md` there. This document is built from those two
files and from the code at `9e17100`.*

## If you read only this

**What was achieved.** The time evolution spends almost all its time computing the rates of change of the state
(amplitudes, black-hole mass and spin) at the points the time integrator asks for. One such computation now costs
0.19 s on one processor core instead of 2.99 s, a factor of about 15. On eight cores it costs 0.07 to 0.11 s, a
factor of 28 to 43; the spread is between two runs on different nodes. A full rerun of your saved window
$t = 1.5\times10^9 \to 10^{11}$ took 295 s on eight cores. The original run took 13353 s, but on a different machine
with a different core layout, so the factor of about 45 is not a like-for-like comparison.

| | Before | Now, 1 core | Now, 8 cores |
|---|---|---|---|
| one rate computation (plateau state) | 2.99 s | 0.189 s | 0.069 to 0.106 s |
| short evolution, 31 rate computations | 94.3 s | 6.4 s | 2.84 to 3.15 s |
| window $1.5\times10^9 \to 10^{11}$, 2233 rate computations | 13353 s (other machine) | not measured | 295 s |

**Two things about your existing results.**

1. The old code computed the amplitude rates with a slightly out-of-date matrix. Individual rate evaluations were
   wrong by up to 1.1%. In the stored linear-growth window this left $A_{322}$ at $t = 1.5\times10^9$ off by
   $7.6\times10^{-5}$ relative. Every later window that was restarted from that state carries the same offset as its
   initial condition. Whether that offset changes late-time results has not been tested (section 3.1).
2. At $\alpha = 0.187$ the radial solve cannot reach its default convergence tolerance, because of roundoff. With the
   old code every solve there ran the maximum 30 Newton iterations, about 90 s. The new code stops early and reports
   the solve as not converged; to run at $\alpha = 0.187$ you pass one extra keyword (section 3.2).

**More bound modes.** The speedup does not change how the cost grows with the number of bound modes $K$. Three modes
are affordable and four are possible for single runs. Of order ten needs changes to the method: the phase average
uses $16^K$ grid points, and the number of Fourier components grows like $K^3$, which makes the linear solve grow like
$K^9$ (section 7).

**What is still pending.** The comparison of time integrators on the final code. The integrator numbers in section 5
were measured on the old code and are preliminary.

**What you need to decide** (details in section 8): whether to merge and push the work; whether to rerun results
computed with the old code; whether the metric-work diagnostic should be on or off by default; which integrator to
use, once the comparison is in; whether to work further on the parts that stop the speedup from growing beyond
eight cores; and how to approach of order ten bound modes (section 7).

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
  (2.7 and 5).
- **Regression test**, **noise floor**, **unit test** (4.1, 4.2).

## 1. Summary

Cost of one rhs evaluation, and of a short evolution segment, on the benchmark's plateau state. The plateau state has
amplitudes $(2.8077\times10^{-3},\,3.396\times10^{-9})$: it is the $|211\rangle$ plateau at $t = 1.5\times10^9$. All
runs on lr7 with the matrix library restricted to one core (section 2.8). The baseline is the code before the project,
on one core.

| Piece | Baseline, 1 core | Final, 1 core | Final, 8 cores |
|---|---|---|---|
| spindown rhs call, `9e17100` (coordinator's harness run) | 2.99 s | 0.189 s | 0.106 s |
| evolution $1.5\times10^9 \to 1.9\times10^9$ (4 steps, 31 calls), `9e17100` | 94.3 s | 6.4 s | 3.15 s |
| spindown rhs call, `2246b4f` (polish agent's run) | 2.99 s | 0.206 s | 0.069 s |
| evolution segment, `2246b4f` (polish agent's run) | 94.3 s | 7.00 s | 2.84 s |
| spindown rhs call, large mesh ($\alpha = 0.187$), `2246b4f` | 6.10 s | 0.402 s | 0.134 s |
| memory allocated per rhs call (2.1, 2.6), `2246b4f` | 6.6 GB | 56 MB | 87 MB |

At `2246b4f` the speedup per rhs call is 14.5 times on one core and 43 times on eight. The 8-core figures differ
between nodes: 0.106 s and 3.15 s in the coordinator's run of the final code, 0.069 s and 2.84 s in the polish agent's
run of `2246b4f`. Treat the 8-core numbers as good to a factor of about 1.5. The baseline and the 1-core numbers agree
between runs much better.

The regression test (168 quantities compared with the baseline's values, section 4.1) passes at one and at eight
cores on the final code. All 278 unit tests pass. The costs of the separate pieces (residual, Jacobian,
factorisation, `balance_terms`) are in section 2 and in `benchmarks/README.md`.

End-to-end rerun of the saved restart window $t = 1.5\times10^9 \to 10^{11}$ at $f_a/M_{\mathrm{Pl}} = 0.0466873$
($f_a = 5.7\times10^{17}$ GeV), $a_0 = 0.9$, $\alpha = 0.3$, cosine potential, Tsit5. Every version took the same 2233
rhs evaluations.

| Code | Wall time | Hardware |
|---|---|---|
| original (Perlmutter job 61110567) | 13353 s | Perlmutter, 1 Julia thread, 4 BLAS threads |
| `01856b9`, first full merge (job 27820263) | 572 s | lr7, 8 Julia threads |
| final code (A/B runs, jobs 27821575 and 27821576) | 295 s | lr7, 8 Julia threads |

The 295 s figure is two runs: 294.6 s with a since-removed correction to chord steps (2.3.5) and 296.4 s without it
(the final code). Both had 2233 calls, 371 accepted and 1 rejected step, and no failed solve. The first-merge run's
trajectory (61 saved times) agrees with the stored CSV to $1.4\times10^{-8}$ relative in $A_{211}$,
$7.4\times10^{-9}$ in $A_{322}$, $3\times10^{-11}$ in $b$ and $7\times10^{-11}$ in $j$. The integrator tolerance was
$10^{-7}$. The two final-code runs agree with each other in final state to $10^{-7}$ or better. The final code's state
at $t = 10^{11}$ agrees with the stored final row to $1.4\times10^{-8}$ in $A_{211}$ and $1.2\times10^{-9}$ in
$A_{322}$. That is a comparison of the final state only, not of the 61 saved times.

The end-to-end ratio, 13353 s / 295 s or about 45, is not like-for-like: the original ran on different hardware with a
different thread layout. The harness numbers in the first table are like-for-like.

Most important caveats:

1. **The old code's rates were wrong by up to about 1% in some evaluations** (section 3.1). The stored linear-window
   result has $A_{322}$ off by $7.6\times10^{-5}$ relative at $t = 1.5\times10^9$. Later windows restarted from it
   inherit that offset. The new code differs from a converged run by about $4\times10^{-7}$.
2. The speedup does not grow past about 8 cores (section 2.8). The factorisation, which runs on one core, and one
   boundary component set a floor.
3. The checks establish that the new code reproduces the old code's numbers within measured noise floors on four
   stored states and one short evolution, and that one full window reproduces. They do not establish discretisation
   accuracy, other parameter points, other potentials, the energy balance (`balance=:energy`), or more than two bound
   modes (section 4).
4. The integrator comparison (section 5) is preliminary and was measured on the old radial code.
5. Nothing is merged or pushed; that needs your sign-off.

*For the cleanup pass: harness jobs at `9e17100`: regression 27821631 (1 thread) and 27821632 (8 threads), unit
tests 27821630. Final-state comparison against the stored CSV: jobs 27821575/27821576, logs
`review-work/rescale_ab/ab_window_*.out` under `axion-sr-sim-speedup/`. First-merge outputs in
`review-work/e2e_01856b9/`.*

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
   also computes $E$ and $J$. After the one-pass change the opt-in saves little, so the default is an open decision
   (section 8).

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
4. A new package dependency for the Rosenbrock methods.

**What it bought.** The ability to use a stiff integrator; the solve counts per method are in section 5. The parallel
forward-difference Jacobian took 9.8 s on one thread and about 5.2 s on 2 to 8 threads, on the old radial code. It was
limited by the slowest column, by garbage collection and by a lock inside UMFPACK solves. It has not been re-measured
on the final code.

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
later window restarted from that state starts with the same offset. Nothing needs to change to use the new code;
whether to rerun old results is your decision (section 8, item 4), and the test that would show whether the offset
matters at late times is to restart a late window from a state produced by the new code.

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
For the final code, the default run (job 27821574, final state in its log) gives $A_{211} = 2.807734301\times10^{-3}$
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
maximum 30 iterations, about 90 s. To run at $\alpha = 0.187$ with the new code, pass for example
`newton_floor_tol=1e-9` to `evolve_with_spindown`.

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

**Result on the final code.** 168 of 168 pass at 1 and at 8 threads. The saved references belong to the baseline
commit; writing new references on an optimised branch would defeat the test.

*For the cleanup pass: `benchmarks/`, branch `speedup/base`, commit `b5b3d67`; states in `benchmarks/states/*.bin`;
references in `reference.toml`, written at `b4c4100`; do not use `--write-reference` on an optimised branch.*

### 4.2 Unit tests

A **unit test** is a small automated check of one function, usually against an independent calculation or against the
function it replaces. 278 tests pass on the final code (run with 4 threads); 279 passed at `2246b4f`, before the
chord-step correction of 2.3.5 was removed. There were 119 before the project.

New tests cover: field kernels against the complex formulas; the dual-number source against finite differences; the
structured Jacobian against automatic differentiation; in-place assembly and boundary derivatives; the block
factorisation against general sparse LU; warm starts and the predictor; one-pass balance terms; the finite-difference
Jacobian (steps, columns, failures, reuse); stiff integrators; the control flow of `solve_radial` and `floor_tol`; and
threaded paths. The threaded-path test starts the code in a separate process with 1 and with 4 threads and checks that
concurrent evaluations equal sequential ones and that results do not depend on the thread count. So the threaded paths
are exercised even if the suite itself runs on one thread.

Does not establish: physical correctness. Most tests compare an optimised function with its reference, so they check
that the two agree, not that either is right.

*For the cleanup pass: `test/runtests.jl`, job 27821630; threaded-path subprocess `test/threads_check.jl`.*

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
  with fixed steps (`adaptive=false`) a failure freezes the evolution. Assigned to the integrator-comparison work; not
  in the final code.
- For integrators without the FSAL property (2.7), the radial solution stored as "accepted", the checkpoint's solution
  and the final radial solution are a stage, not the accepted state, contrary to the documentation. Assigned to the
  integrator-comparison work, whose branch has an unmerged commit on this. The docstring of `evolve_with_spindown` in
  the final code still says the stored accepted solution is the radial solution at the integrator's current state.
- Test gaps (threaded paths, the dual-number source against an independent check, control flow, block LU shapes).
  Addressed in the second round.
- Small code items. Several were fixed in the second round; the rest are in section 9.

Establishes: the combined code on the three integrator paths has no wrong-number defect that the review's experiments
could expose. Does not establish: absence of defects on paths not run (FBDF, other states or potentials).

*For the cleanup pass: review of `01856b9`, scripts in `review-work/claude-review/`. Failure flag: `solve_failed`.
Accepted-state fields: `accepted_solution` (docstring refers to `integrator.p.accepted_solution` at `integrator.u`),
`final_radial_solution`. Fix on `speedup/methods`: `2a1089a`. Test gaps addressed in `eabe89e`.*

### 4.5 End-to-end comparison

See section 1: same steps and the same trajectory within tolerance for one run ($f_a/M_{\mathrm{Pl}} = 0.0466873$,
$a_0 = 0.9$, $\alpha = 0.3$, cosine, window $1.5\times10^9 \to 10^{11}$), and the final code's final state agrees with
the stored final row. Does not establish: other $f_a$, other $\alpha$, other potentials. The original ran on different
hardware.

## 5. Time integrators

**Preliminary.** This section reports the integrator agent's measurements on the **old radial code**, at one parameter
point ($a = 0.9$, $\alpha = 0.3$, $f_a = 5.7\times10^{17}$ GeV), with tolerance $10^{-7}$. The quantity compared is
the number of radial solves, which does not depend on how fast one solve is.

**Finding.** The problem is stiff only on the $|211\rangle$ plateau, with fast eigenvalue $-3.27\times10^{-8}$ against
slow ones of absolute value below $5\times10^{-13}$ (see 2.7 for what this means and why it caps Tsit5's step).
Elsewhere Tsit5 is efficient.

Radial solves per segment ("-" means not run):

| Segment | Tsit5 | Rosenbrock23 | Rodas5P | FBDF |
|---|---|---|---|---|
| plateau, $1.5\times10^9$ to $1.5\times10^{10}$ | 758 | 82 | 171 | 60 |
| takeover, $1.5\times10^{10}$ to $4.3\times10^{10}$ | 1141 | 640 | 500 | - |
| $|211\rangle$ decay, $4.3\times10^{10}$ to $10^{11}$ | 334 | 630 | 2700 | - |
| whole window $1.5\times10^9$ to $10^{11}$ | 2233 | 1354 | 3372 | 638 (about 510 with the final callback) |

Final-state error in $A_{211}$ against a tight reference: Tsit5 $2\times10^{-7}$, Rosenbrock23 $1.5\times10^{-4}$,
Rodas5P $9\times10^{-6}$, FBDF $1.1\times10^{-4}$.

**Tolerance requested versus error achieved.** All four runs requested $10^{-7}$. The integrator applies that
tolerance to its own estimate of the error made in each step, and chooses the step size to keep that estimate small.
The error at the end of a run accumulates over many steps and depends on how good the estimate is, so it can be much
larger than the request, as it is here for the stiff methods. The only reliable measure is a comparison with a tighter
reference run.

**Vocabulary for the methods.**

- The **order** $p$ of a method: the error after a fixed time scales as $h^p$, so halving the step reduces it by
  $2^p$.
- A **multistep** method builds $y_{n+1}$ from several previous states instead of from stages within one step. The
  backward differentiation formulas (**BDF**) are implicit multistep methods made for stiff problems.
- A **W-method** is a Rosenbrock method whose accuracy does not require the exact $\Lambda$; any reasonable
  approximation keeps the order. That matters here because $\Lambda$ comes from finite differences.
- **L-stable**: as $h\lambda \to -\infty$ the method's amplification factor goes to zero, so a fast relaxation is
  damped out in one step instead of persisting as a slowly decaying oscillation.

**The four methods.**

- **Tsit5**: explicit, order 5. No $\Lambda$ needed, but its step is held at the stability limit on the plateau.
- **Rosenbrock23**: order 2, linearly implicit. It is a W-method, so a finite-difference $\Lambda$ is acceptable. About
  6 radial solves per step. Low accuracy per step. (Its documentation describes it as L-stable.)
- **Rodas5P**: order 5, about 13 radial solves per step. Its accuracy degrades if $\Lambda$ is out of date.
- **FBDF**: a variable-order BDF method. It reuses $\Lambda$ over many steps by design. It had the fewest solves in
  every window it ran. It needs a new dependency (tested in scratch only) and a check at tighter tolerance, since its
  final error here is $1.1\times10^{-4}$.

**Agent's recommendation (preliminary).** A stiff method on plateaus and Tsit5 elsewhere; Rodas5P rather than
Rosenbrock23 when accuracy matters; FBDF is the best overall candidate. Reusing $\Lambda$ gave erratic errors away
from the plateau and stays off. On the final code the cost of one solve has changed but the number of solves per method
has not. The cost of one $\Lambda$ ($K + 2$ extra solves) now scales with the faster solve and with threads.

The work-precision comparison on the final code is pending. A work-precision comparison measures, for each
integrator, the cost (here radial solves) needed to reach a given error against a reference solution.

*For the cleanup pass: comparison on branch `speedup/methods`; FBDF dependency OrdinaryDiffEqBDF.*

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
| `AXION_ALGORITHM` | `tsit5` (default), `rosenbrock23` or `rodas5p`. Stiff integrators get the finite-difference Jacobian; output files get the suffix `_<algorithm>`. |
| `AXION_TOLERANCE` | Sets reltol = abstol. Default is the library's $10^{-7}$. |
| `AXION_CHECKPOINT_MINUTES` | Writes `<prefix>_checkpoint.jls` after the first accepted step at least that many minutes after the last. If the file exists, the same command resumes from it. A resumed run agrees with an uninterrupted run to the tolerance, not bit for bit. |
| `AXION_RECORD_METRIC_WORK` | `1` computes $E$, $J$ and the metric work for the diagnostics (`w` column); otherwise `w` is `NaN`. Rates are unchanged. |
| `AXION_POTENTIAL` | `cosine` (default) or `free` (existing). |

New keywords on `evolve_with_spindown` (and, where it applies, `evolve_fixed_background`; both also take
`EvolutionCache` options):

| Keyword | Meaning |
|---|---|
| `tangent_predictor=true` | Start solves from the first-order prediction (2.3.2). `false` is the reference path. |
| `record_metric_work=false` | Also compute $E$, $J$ and metric work under the action balance (diagnostic only). |
| `newton_floor_tol=newton_tol` | `floor_tol` of each radial solve; raise it only at a known residual floor (3.2). |
| `algorithm=Tsit5()` | Integrator; `Rosenbrock23()` and `Rodas5P()` are re-exported from the library. |
| `jacobian_settings=FiniteDifferenceJacobian()` | The rates Jacobian $\Lambda$ for implicit and Rosenbrock integrators (2.7). |

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
julia --project=. -e 'using Pkg; Pkg.test()'                    # default threads, about 3 min in a Slurm job
JULIA_NUM_THREADS=4 julia --project=. -e 'using Pkg; Pkg.test()'    # with 4 threads
```

The 4-thread line reflects how the coordinator ran the suite ("tests at 4 threads"). I did not find the exact command
line in the sources, so check how the thread count reaches the test process on your Julia version. The threaded-path
test starts its own 1-thread and 4-thread processes either way.

Manifest. `Manifest.toml` is now written by Julia 1.13.0 (`julia_version = "1.13.0"`); before, it was written by
1.11.7. The README still gives Perlmutter instructions with `julia/1.11.7`. The dependency change is in its own commit
(`b7e315e`, adding OrdinaryDiffEqRosenbrock 2.7.5 with 35 Manifest entries), so it can be redone on the other cluster
by adding the package there and regenerating the Manifest with that cluster's Julia. Using the 1.13 Manifest with Julia
1.11.7 was not tested.

## 7. Scaling to more bound modes

You want the method to handle of order ten bound modes eventually. $K$ is the number of bound modes (in the code, the
length of the `ModeSet`). This section says what grows how fast with $K$ in the code as it stands. The short version:
$K = 3$ is affordable, $K = 4$ is possible for single runs, and $K \approx 10$ needs changes to the method, not only to
the quadrature. None of the counts below was changed by the speedup work.

**What the calculation does with the modes.** The field is expanded in Fourier components $e^{-i q\cdot\xi}$, where
$\xi$ is the vector of the $K$ phases and $q$ is a vector of $K$ integers. The code keeps every $q$ with odd
$|q|_1 = \sum_i |q_i| \le 3$ and $\sum_i q_i > 0$: the $K$ fundamentals and all third-order combinations. Each
component has its own radial profile, so the number of components sets the number of unknowns at every mesh node.
The self-interaction term is obtained by reconstructing the field on a grid in $\xi$ (and in the polar angle),
applying the potential, and projecting back onto the components.

**Three things grow with $K$.**

1. *The phase grid: exponential.* The grid in $\xi$ is a tensor product with $n_{\mathrm{phase}} = 16$ points per phase,
   so it has $16^K$ points. Every evaluation of the self-interaction term at every mesh node does work proportional to
   that number, and the tables of sines and cosines at the grid points take memory proportional to it.
2. *The number of Fourier components: polynomial, but steep.* Counting the integer vectors above gives
   $2K + 2K(K-1) + 4\binom{K}{3}$ components, which grows like $\tfrac{2}{3}K^3$. (Checked against the code for
   $K = 2, 3, 4$.) With one angular channel per component the node block of the Jacobian has four reals per component.
   The block LU of section 2.2 costs the cube of the block size per mesh node, so the linear solve grows roughly like
   $K^9$, and storing the Jacobian blocks grows like $K^6$.
3. *Building the list of components: exponential, but trivially fixable.* The `ModeSet` constructor lists all $7^K$
   integer vectors with entries from $-3$ to $3$ and filters them. At $K = 10$ that is $2.8\times10^8$ candidates; a count
   I attempted for $K$ up to 10 did not finish in four minutes. Generating only the vectors that pass the filter removes
   this.

Everything else is mild: the rates Jacobian for a stiff integrator has $K + 2$ columns, and the derivative pass of the
balance terms carries $K + 2$ dual partials.

| $K$ | phase points $16^K$ | components | reals per mesh node | distinct wavevectors $q \pm q'$ | phase tables | stored Jacobian | node cost relative to $K = 2$ |
|---|---|---|---|---|---|---|---|
| 2 | 256 | 8 | 32 | 25 | 0.2 MB | 12 MB | 1 |
| 3 | 4096 | 22 | 88 | 116 | 12 MB | 88 MB | about 70 |
| 4 | 65536 | 48 | 192 | 417 | 590 MB | 415 MB | about 4000 |
| 10 | $1.1\times10^{12}$ | 680 | 2720 | about $5\times10^4$ | not feasible | about 80 GB | not feasible |

Sources and limits of this table: the first four columns for $K \le 4$ were read from the code on 2026-11-13 with the
example modes $|211\rangle, |322\rangle, |433\rangle, |544\rangle$; the $K = 10$ row is from the counting formula and my own
count of integer vectors with $|k|_1 \le 6$; memory figures and the "node cost" column are computed from the sizes and
from the assembly agent's cost formula for the source Jacobian (potential evaluations at every phase point and polar
node, plus a projection onto every distinct wavevector), not measured. The per-thread scratch arrays add about 1.6 MB
per thread at $K = 3$ and 25 MB at $K = 4$ (estimate).

**What this means in run time (extrapolated, not measured).** About 0.15 s of today's one-core right-hand-side
evaluation is work proportional to the phase grid. At $K = 3$ that becomes of order 10 s on one core, which roughly
cancels the speedup of this project. At $K = 4$ it becomes 4 to 10 minutes on one core and about a minute on eight,
with the factorisation at about 3.5 s; a run of 2000 evaluations would take a day or two on eight cores, so a single
run is feasible and a parameter scan is not. At $K = 10$ the tensor grid alone has $10^{12}$ points.

**What would remove each obstacle.**

- *Phase grid.* The $16^K$ points are used to compute a number of Fourier coefficients that grows only polynomially:
  the distinct wavevectors $q \pm q'$ in the table. A quadrature rule matched to those wavevectors, for example a rank-1
  lattice rule (points $\xi_n = 2\pi n z / P$ for one integer vector $z$, chosen so that no two needed wavevectors
  alias), needs a number of points of that order. The new node kernels take the phase points as a plain list with
  precomputed tables, so they do not need rewriting for a different point set. This changes the quadrature, so it needs
  a convergence comparison against the tensor grid at $K = 2$ and $3$. It is the first item of the numerical-parameter
  audit (section 8) and has not been started.
- *Component list.* Generate the vectors directly.
- *Number of components.* The quadrature does not help here. Treating all $K$ modes on an equal nonlinear footing is
  what produces the $K^3$ growth. A possible reformulation (a suggestion, not designed or tested): keep the two or
  three modes that are large at a given time fully nonlinear, as now, and treat each other mode as a linear
  perturbation on that nonlinear background, with a rule for promoting a mode to the nonlinear set when it has grown.
  The cost would then be one small nonlinear solve plus one linear solve per extra mode, linear in $K$. This is close in
  spirit to how Omiya et al. treat the secondary cloud. Whether it captures the physics you want from a ten-mode
  calculation is a physics judgement; if it does, it is a separate piece of method design on top of the present code.

*For the cleanup pass: `PhaseGrid` and its point-list tables `synthesis`, `analysis`, `jacobian_table` in `src/field.jl`;
`ModeSet` and `_integer_vectors` in `src/modes.jl` (`max_order`, `n_angular`); the agent's cost formula uses $P$ phase
points, $X$ polar nodes, $S$ distinct wavevectors up to sign, $L$ channels, and $h'$ the derivative of the nonlinearity.*

## 8. Open items and decisions for you

1. **Default of the metric-work diagnostic** (`record_metric_work`). Off now under the action balance, so `w` is
   `NaN` unless requested. Option off: saves 30 to 65 ms per evaluation on one core (measured before later speedups;
   the script comment says "a few percent"), at the price of a `NaN` column. Option on: simpler output, costs that
   time. The coordinator was to decide from the final cost per evaluation; no decision or recommendation is recorded.
2. **Numerical-parameter audit** ($n_{\mathrm{phase}}$, $n_{\mathrm{polar}}$, mesh step, tolerances): convergence
   tables for you to decide on. First item: a phase quadrature that is not $n_{\mathrm{phase}}^K$. Not started.
3. **Integrator choice.** Waiting for the work-precision comparison on the final code. The preliminary recommendation
   is in section 5.
4. **Existing results.** Whether to rerun anything computed with the old code, given 3.1: $A_{322}$ offsets of
   $7.6\times10^{-5}$ in the linear window, inherited by restarted windows. For the one parameter point measured, the
   new code took 64 s for the linear window and 295 s for the long window on 8 cores. The test of whether the offset
   matters is to restart a late window from a state produced by the new code. No recommendation is recorded.
5. **Remaining serial floor.** The factorisation (16 ms) and one boundary component (8 to 9 ms) limit scaling beyond 8
   threads, and garbage collection is 16 to 20% of wall time at 8 to 16 threads. Whether to work on this is not
   decided; the possible gain was not measured. A cost-aware refactorisation rule (about 25% fewer steps in long chord
   sequences) was not implemented.
6. **Review findings not yet in the code**, assigned to the integrator-comparison work: the failed-solve flag with
   integrators that run their own nonlinear solve; accepted state versus stage for integrators without FSAL (4.4).
7. **The action balance still computes $E$ and $J$ once per evaluation for diagnostics.** The review asks to document
   this or make it optional. Related to item 1.
8. **Unexplained difference** (Newton agent, preliminary). Nonlinear state, fixed background, fixed step: old and new
   code differ by $2.1\times10^{-8}$ in $\ln A_{211}$ after 6 steps, while each agrees with its own tighter-tolerance
   run to $3\times10^{-11}$. The agent's note named the out-of-date factorisation (3.1) as the likely cause, "to be
   checked". I did not find that it was checked.
9. **Merge and push.** Nothing is merged into `rewrite/axioncloud`; this needs your sign-off.
10. **Reproduction on Perlmutter** with Julia 1.11.7 (section 6).

## 9. Notes for the cleanup pass

Each name below was checked against `9e17100` with grep. Items marked "already removed" were listed by an agent but
are gone.

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
- [ ] `total_solves` (`src/evolution.jl`): defined and exported, no use found.
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
- [ ] Docstring corrections: the review listed some; its list is not in the sources. One known case: the
      `evolve_with_spindown` docstring says `integrator.p.accepted_solution` is the solution at `integrator.u`, which is
      not true for integrators that are not first-same-as-last (4.4).

Validation scripts and scratch:

- [ ] `validation/newton_counts.jl` is stale: its wrappers hook `jacobian` and the sparse `ScaledFactorization`,
      which `solve_radial` no longer calls, so it no longer counts Jacobians and factorisations. Fix or retire.
- [ ] Scratch directories and worktrees under `axion-sr-sim-speedup/` that can be deleted once nothing is needed from
      them: worktrees `base/`, `assembly/`, `linear-solver/`, `newton/`, `balance/`, `integrator/`, `integration/`,
      `polish/`, `methods/` and `docs/` (after merge), and the scratch directories `balance-scratch/`,
      `integrator-work/`, `linear-solver-work/`, `methods-work/`, `polish-work/`, `review-work/`. Keep `review-work/`
      until the end-to-end outputs (`e2e_01856b9/`) and the A/B logs (`rescale_ab/`) cited here are no longer needed.
      Untracked scratch inside worktrees (`.scratch/` in `assembly/`, `newton/work/`) goes with the worktrees.
      Remove worktrees with `git worktree remove`, not `rm`.
- [ ] `benchmarks/logs/` is git-ignored, so harness logs are not in the repository; copy any you want to keep.

## 10. Branches and commits

Nothing here is merged into `rewrite/axioncloud` or pushed. All branches start from `b4c4100`.

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
| `speedup/integration-next` | `9e17100` | Everything above merged (26 commits on `b4c4100`), chord-step correction removed. The branch to review or merge |
| `speedup/methods` | `2a1089a` | In progress, 2 commits beyond `9e17100`: OrdinaryDiffEqBDF dependency (`07287c0`), FBDF and QNDF, per-attempt failure record, solutions at accepted and final states (`2a1089a`). Work-precision comparison pending |
| `speedup/docs` | this document's commit | This file only |

The 26 commits of `git log --oneline b4c4100..9e17100` are: the harness (`b5b3d67`); the dependency commit
(`b7e315e`); the four first-round merges and their commits (`80b496f`, `83e8bb3`, `05feccb`, `cf853b6`, `cea6095`,
`7b6cafd`, `02fda5c`, `753d59a`, `c5e5125`, `e1fca6d`, `fcd98c7`, `1a6f91e`); merge commits (`9950e6d`, `35f6928`,
`efd6772`, `229b09f`, `01856b9`); and the polish commits (`9ec3265`, `50ddeb9`, `a9754a4`, `5f8d0f6`, `eabe89e`,
`2246b4f`) followed by `9e17100`.

## Unresolved while writing

- **Which code the "new default run" in the linear-window table is.** STATUS.md's table and its quoted distances
  ($4\times10^{-7}$, $3.4\times10^{-7}$) are the run with the chord-step correction (job 27821573, code at `2246b4f`).
  The final code is the run without it (job 27821574); I computed $5.2\times10^{-7}$ and $3.8\times10^{-7}$ from its
  logged final state. The conclusion does not change.
- **Which code the 295 s end-to-end figure is.** STATUS.md calls it "9e17100-equivalent". The two runs are 294.6 s
  (with the correction, `2246b4f`) and 296.4 s (without, `9e17100`).
- **8-thread timing spread** between the coordinator's run (0.106 s, 3.15 s) and the polish agent's (0.069 s, 2.84 s).
  I give both and did not choose.
- **Test count.** 278 at `9e17100`, 279 at `2246b4f`. I take the difference to be the removed correction's test, which
  I could not confirm from the sources.
- **"Serial block factorisation".** STATUS.md calls the factorisation serial (16 ms floor). In the code only the
  copying of blocks into the stage panels is threaded; I read the elimination stages as serial from the code but did
  not time them separately.
- **Whether the harness regression was run against the original `b4c4100` references at `9e17100`.** STATUS.md says
  "regression PASS"; the references are the baseline's by design.
- **The review's docstring list and the outcome of the two chord and line-search weak spots** are not recorded in the
  sources.
