# Numerical-parameter audit of AxionCloud.jl (2026-11-13)

Written for the owner before the production runs are repeated with the sped-up code. For each numerical parameter
this document says what it controls, how much the quantities that drive an evolution change when it is refined, what
one evaluation costs at each setting, and a recommendation. Status: complete for the stored states; the rows marked
*pending* are filled in as the last Slurm jobs finish.

## If you read only this

| Parameter | Production | Error at production (rate error / fastest rate, worst state) | Recommendation | Cost change per evaluation (1 / 8 threads) |
|---|---|---|---|---|
| `n_phase` (tensor phase grid) | 16 | $3\times10^{-15}$ | loosen to 12 | $-26\%$ / $-15\%$ |
| rank-1 lattice instead of the tensor grid | (not used) | 85 points: $2\times10^{-13}$ | adopt the 85-point lattice if the new code path is accepted (replaces the line above) | $-42\%$ / $-19\%$ |
| `n_polar` | 12 | $2\times10^{-15}$ | loosen to 8 | $-22\%$ / $0\%$ |
| radial `step` | 0.4 | $2.5\times10^{-6}$ | keep | |
| outer endpoint (`containment`) | $10^{-4}$ | $6\times10^{-4}$ (strong takeover), $4\times10^{-6}$ elsewhere; **does not converge** | not a parameter fix: see Section 3.2 | |
| inner endpoint (`horizon_offset_ratio`) | $10^{-6}$ | $\le 10^{-7}$ on the $\alpha = 0.3$ states (strong states pending) | keep | |
| `stretch_factor` | 1.03 | $3\times10^{-7}$ | keep | |
| `n_angular` (polar indices per component) | 1 | $7\times10^{-4}$ (plateau, takeover) | the largest error found; decide whether to run with 2 | $+91\%$ / $+148\%$ |
| `max_order` (Fourier truncation) | 3 | $\lesssim 3\times10^{-6}$ (mixed with the outer endpoint) | keep | |
| `newton_tol` | $10^{-10}$ | $\le 2\times10^{-9}$ even at $10^{-7}$ | keep on the $\alpha = 0.3$ meshes; use $10^{-8}$ on the $\alpha = 0.187$ mesh | $0$ to $-14\%$; $-28\%$ at $\alpha = 0.187$ |
| integrator `reltol`/`abstol` | $10^{-7}$ | see the methods report | keep $10^{-7}$ (Tsit5) or FBDF at $10^{-9}$-$10^{-10}$ on stiff plateaus | |

The two quadrature changes together (lattice or `n_phase` = 12, plus `n_polar` = 8) should save roughly a third to
half of an evaluation at 1 thread and 15-20% at 8 threads; the combination has not been run as one setting, so it
needs one confirmation run before production.

Three findings matter more than any saving:

1. **The angular truncation is the largest error.** Allowing a second polar index per Fourier component changes the
   amplitude rates by $7\times10^{-4}$ on the plateau and $6\times10^{-4}$ at the takeover. That is 300 times the
   radial-step error and $10^{4}$ times the integrator tolerance. Whether `n_angular` = 2 is itself converged is not
   known (3 was not run).
2. **The outer-endpoint dependence does not converge.** Moving the outer endpoint out changes the rates linearly in
   the endpoint radius, without levelling off, between $r_* = 428$ and $860$. At the strong takeover the cause is
   identified: outgoing radiation in the $q = (-1, 2)$ component inside the domain is counted in the phase actions.
3. **The integrator tolerance is already tighter than the discretisation.** The radial step alone gives relative rate
   errors of about $2\times10^{-6}$. Tightening `reltol` below $10^{-7}$ does not make the result more accurate.

## How the audit was done

**Quantities.** One spindown right-hand-side evaluation (`_spindown_rhs!`, action balance, as in
`scripts/evolution_spindown.jl`) gives the four numbers the integrator uses: $d\ln A_{211}/dt$, $d\ln A_{322}/dt$,
$db/dt$ and $dj/dt$, where $b = M/M_0$ and $j = J/M_0^2$ are the black hole's mass and angular momentum. The tables
also give the complex mode frequencies $\omega_{211}$, $\omega_{322}$ of the converged radial solution
($\mathrm{Im}\,\omega$ is the growth rate of the mode).

**States.** Seven states that span weak to strong nonlinearity (amplitudes $A$ are the leading outer coefficients of
the fundamental channels; the first five are at $f_a/M_\mathrm{Pl} = 0.0466873$, $f_a = 5.7\times10^{17}$ GeV):

| State | Where it comes from | $A_{211}$, $A_{322}$ |
|---|---|---|
| linear | harness, seed of the evolution | $9.9\times10^{-7}$, $3.1\times10^{-9}$ |
| plateau | harness, $t = 1.5\times10^9$ | $2.8\times10^{-3}$, $3.4\times10^{-9}$ |
| takeover | $t = 2.94\times10^{10}$ of the restart run's CSV, $|322\rangle$ rising fast (new) | $2.6\times10^{-3}$, $2.0\times10^{-6}$ |
| nonlinear | harness, $t = 10^{11}$ | $1.5\times10^{-4}$, $4.1\times10^{-6}$ |
| large_mesh | harness, $\alpha = 0.187$, 1355 nodes | $10^{-6}$, $10^{-6}$ |
| strong_peak | $f_a/M_\mathrm{Pl} = 3\times10^{-3}$ ($3.7\times10^{16}$ GeV) run, $|211\rangle$ peak (new) | $2.1\times10^{-2}$, $1.9\times10^{-6}$ |
| strong_takeover | same run, both modes grown (new) | $1.6\times10^{-2}$, $3.2\times10^{-5}$ |

The harness's "nonlinear" state is in fact weaker than the plateau (its $A_{211}$ is 20 times smaller), so the two
strong states were added from the $f_a/M_\mathrm{Pl} = 3\times10^{-3}$ runs, whose amplitudes are about 7 times the
plateau's. Without them the phase quadrature would have looked insensitive to everything.

**Re-solving at a new discretisation.** When only the quadrature changes (`n_phase`, `n_polar`, lattice), the unknown
vector has the same layout and the stored harness vector is the Newton starting point. When the mesh or the
truncation changes, the state is rebuilt from scratch the way `benchmarks/build_states.jl` does: a cold solve for the
two linear states, otherwise a cold solve at $10^{-3}$ of the amplitudes and 30 geometric amplitude steps. Every state
is then solved to a scaled residual of $10^{-12}$ (the $\alpha = 0.187$ mesh stops at its roundoff floor, between
$2\times10^{-10}$ and $2\times10^{-9}$), and the right-hand side is evaluated with all its radial solves at the same
tolerance, so that Newton error stays far below the changes being measured.

**Measures.** Each table gives every quantity's relative change from the finest setting of that table,
$|x - x_\mathrm{ref}|/|x_\mathrm{ref}|$. On the plateau $\mathrm{Im}\,\omega_{211}$, $d\ln A_{211}/dt$, $db/dt$ and
$dj/dt$ sit near zero, because the black hole has spun down to the superradiance threshold of $|211\rangle$; their
relative changes there look large (up to $10^{-1}$) while the absolute changes are of order $10^{-19}$ per unit time.
The column **rate error / fastest rate**, $\max_i |\Delta\, d\ln A_i/dt| / \max_j |d\ln A_j/dt|$, measures the error of
either amplitude rate against the faster of the two, and is the one to read. A relative rate error $\epsilon$ held
over an evolution in which an amplitude grows by $n$ e-folds gives an error of about $n\epsilon$ in $\ln A$, that is a
relative amplitude error of $n\epsilon$, or a shift of $n\epsilon$ e-folding times in when things happen. With $n$
of order 10, $\epsilon = 10^{-6}$ corresponds to $10^{-5}$ in the amplitudes, against about $2\times10^{-7}$ for the
integrator at `reltol` = $10^{-7}$.

**Cost.** The minimum wall time of one right-hand-side evaluation at a slightly perturbed point (the harness's fixed
perturbation: $A_{322}\times1.01$, $b$ and $j$ shifted by $10^{-6}$, so that the radial solve does 3-7 Newton
iterations as in an evolution), warm-started from the converged state, at production Newton tolerance, at 1 and at 8
Julia threads, on lr7 (Xeon Gold 6330), with BLAS on one thread. Times vary by about 10% between nodes and runs;
differences below that are not significant, and the 8-thread times are noisier than the 1-thread ones.

**Trajectory.** For every setting the harness's evolution segment $t = 1.5\times10^9 \to 1.9\times10^9$ (31
right-hand-side calls, Tsit5 at $10^{-7}$) was also run from the plateau state. It never changed the number of calls.
The first segment of each job includes compilation, so its wall time is not comparable.

**What these tables establish and what they do not.** Each table varies one parameter with all others at production,
so it shows that parameter's own error, not joint convergence; the "all finer" table gives one joint point. All of it
is self-convergence: it shows how far production is from a finer discretisation of the same equations, not that the
equations are right. One $\alpha$ (0.3, plus one linear state at 0.187), one spin, two values of $f_a$.

**Scripts and data.** `validation/parameter_audit.jl` (states, studies, solves, rates, timing, segment; it is in
`validation/` so that the owner can rerun it on the final code), `validation/parameter_audit_submit.sh` (one lr7 job
per group of studies), `validation/parameter_audit_tables.jl` (these tables). Results: CSV files and logs in
`../audit-work/run1` (all studies, five states), `run2_large` (the $\alpha = 0.187$ state with the corrected floor
tolerance), `run4_strong` (strong states), `run5_lattice`, `run6_tail`. Slurm jobs 27824736-27824738, 27824847,
27824972, 27825063, 27825217. The library defaults are unchanged. The branch also carries commit b3bcdb4 (the
`phase_grid` keyword and `rank1_lattice`, cherry-picked from `speedup/scaling`), which the lattice rows need.

## 1. Phase quadrature

### 1.1 `n_phase`: points per phase in the tensor grid

The nonlinear source and its Jacobian are averages over the two independent phases $\xi_1, \xi_2$ of the two bound
modes. The code evaluates the field on an $n \times n$ grid of equally spaced phases (a *tensor grid*, $n^2 = 256$
points at `n_phase` = 16) and averages there. This average is exact for every Fourier harmonic $e^{i k\cdot\xi}$
with $|k_i| < n$; a higher harmonic of the cosine nonlinearity is *aliased*, that is, the grid cannot tell it from a
lower one and adds it to the wrong place. The aliasing error therefore falls as fast as the high harmonics of the
field do, which depends on the field amplitude.

Summary over all seven states (full table in the appendix):

| `n_phase` | points | rate error / fastest rate | $db/dt$, $dj/dt$ | cost 1 thr | cost 8 thr |
|---|---|---|---|---|---|
| 8 | 64 | $2.8\times10^{-7}$ | $1.3\times10^{-7}$ | 0.53 | 0.85 |
| 12 | 144 | $1.5\times10^{-14}$ | $1.3\times10^{-15}$ | 0.74 | 0.85 |
| 16 (production) | 256 | $2.5\times10^{-15}$ | $6.7\times10^{-16}$ | 1.00 | 1.00 |
| 24 | 576 | $2.7\times10^{-14}$ | $1.2\times10^{-15}$ | 1.78 | 1.40 |
| 32 (reference) | 1024 | 0 | 0 | 2.86 | 1.97 |

The changes from 12 upwards are at roundoff (the differences between 12, 16 and 24 do not decrease, because they are
already noise). Only `n_phase` = 8 shows an aliasing error, $3\times10^{-7}$ at the two strong states and
$2\times10^{-10}$ on the plateau. The weak states do not see the phase grid at all.

**Recommendation: loosen to `n_phase` = 12**, saving 26% of an evaluation at 1 thread and about 15% at 8 threads,
with a margin of seven orders of magnitude between its error and the next smaller grid's at the strongest state
audited. If runs go to larger amplitudes than $A_{211} \approx 2\times10^{-2}$ (smaller $f_a$), recheck with
`n_phase` = 8 versus 12 at the new peak state; the 8 row is the early warning.

### 1.2 Rank-1 lattice rules instead of the tensor grid

A *rank-1 lattice rule* with $M$ points and generating vector $z = (1, g)$ (an idea from the scaling study on
`speedup/scaling`) uses the points $\xi_j = 2\pi\,(j z \bmod M)/M$, $j = 0, \dots, M-1$. It is exact for every
harmonic $k$ except those with $k \cdot z \equiv 0 \pmod M$. Unlike the tensor grid, it protects a diamond
$|k_1| + |k_2| \le L$ of harmonics, which is the shape that the nonlinearity actually fills, and its size does not
grow like $n^K$ with the number of modes $K$. Generating vectors are those of the scaling study's Korobov search
(aliases up to order 24): $(1,12)$, $(1,11)$, $(1,13)$, $(1,18)$, $(1,20)$ for $M$ = 31, 61, 85, 101, 151. The
85-point rule has no alias of order $\le 12$.

This is the first test of the lattice in full solves and right-hand-side evaluations (the scaling study tested it at
single mesh nodes). Summary over all seven states, against the 32-point tensor grid:

| phase rule | points | rate error / fastest rate | $db/dt$, $dj/dt$ | cost 1 thr | cost 8 thr |
|---|---|---|---|---|---|
| lattice 31 | 31 | $1.8\times10^{-2}$ | $2.1\times10^{-2}$ | 0.45 | 0.83 |
| lattice 61 | 61 | $1.2\times10^{-9}$ | $3.1\times10^{-11}$ | 0.52 | 0.85 |
| lattice 85 | 85 | $2.2\times10^{-13}$ | $5.3\times10^{-15}$ | 0.58 | 0.81 |
| lattice 101 | 101 | $8.7\times10^{-15}$ | $7.7\times10^{-16}$ | 0.63 | 0.84 |
| lattice 151 | 151 | $1.8\times10^{-15}$ | $1.3\times10^{-15}$ | 0.76 | 0.86 |
| tensor 16 (production) | 256 | $2.5\times10^{-15}$ | $6.7\times10^{-16}$ | 1.00 | 1.00 |
| tensor 32 (reference) | 1024 | 0 | 0 | 2.87 | 1.98 |

The plateau segment with the 61-151 point lattices ends within $1.4\times10^{-13}$ (in $\ln A_{322}$) of the tensor
grid's, with the same 31 calls. The 31-point lattice fails visibly ($10^{-2}$ at the strong takeover), which shows
that the test can fail. Every solve converged.

**Recommendation:** the 85-point lattice reproduces the tensor grid to $2\times10^{-13}$ at every state and costs
0.58 of production at 1 thread (0.81 at 8). It does as well as `n_phase` = 12 at 60% of its points. If the owner
accepts the `phase_grid` code path for production, use `phase_grid = PhaseGrid(rank1_lattice([1, 13], 85),
mode_set)`; otherwise use `n_phase` = 12. At $K = 2$ the gain over `n_phase` = 12 is modest (about 20% at 1 thread,
little at 8, because at 8 threads the phase average is no longer the largest cost); its real value is that it keeps
working when $K > 2$.

## 2. `n_polar`: polar quadrature nodes

The field is reconstructed at `n_polar` Gauss-Legendre nodes in $x = \cos\theta$ before the nonlinearity is applied and
projected back onto the polar harmonics. A Gauss-Legendre rule with $N$ nodes is exact for polynomials in $x$ up to
degree $2N - 1$; the error comes from the part of the nonlinear source that is not such a polynomial, which again
falls with the field amplitude.

| `n_polar` | rate error / fastest rate | $db/dt$, $dj/dt$ | cost 1 thr | cost 8 thr |
|---|---|---|---|---|
| 8 | $3.6\times10^{-14}$ | $1.8\times10^{-15}$ | 0.78 | 1.00 |
| 10 | $4.4\times10^{-15}$ | $7.7\times10^{-16}$ | 0.86 | 0.98 |
| 12 (production) | $1.7\times10^{-15}$ | $7.7\times10^{-16}$ | 1.00 | 1.00 |
| 16 | $1.2\times10^{-15}$ | $7.7\times10^{-16}$ | 1.21 | 1.18 |
| 24 (reference) | 0 | 0 | 1.69 | 1.43 |

Even 8 nodes are converged to $4\times10^{-14}$ at the strong states. Below 8 was not tried: the retained channels
reach $l = 6$, and products of three such harmonics need more nodes, so 8 is the smallest value with a clear
argument.

**Recommendation: loosen to `n_polar` = 8**, saving 22% at 1 thread; at 8 threads the saving is within the timing
noise.

## 3. Radial mesh

`RadialMesh(background, α, mode_set; step, stretch_factor, max_step, containment, horizon_offset_ratio)` exposes the
knobs used here. The mesh is uniform with spacing `step` in the tortoise coordinate $r_*$ out to twice the largest
bound-profile peak radius, then *stretched*: each interval is `stretch_factor` times the previous one, up to a cap
`max_step` (default: the smaller of 8 `step` and a quarter wavelength of the fastest propagating component, which is
1.86 here). The endpoints come from `domain_estimate`: the inner one where the horizon offset is
`horizon_offset_ratio` times the gap $r_+ - r_-$, the outer one where every bound component has fallen below
`containment` times its peak. Production: 0.4, 1.03, 1.86, $10^{-4}$, $10^{-6}$ (690 nodes, $r_* \in [-44.3, 517.9]$
at $\alpha = 0.3$).

### 3.1 `step`: radial step

The radial equations are discretised by a fourth-order MIRK (mono-implicit Runge-Kutta) scheme: the error of the
discrete solution falls like $h^4$ with the step $h$. This is the *discretisation error* of the radial problem.

| `step` | nodes | rate error / fastest rate | cost 1 thr | cost 8 thr |
|---|---|---|---|---|
| 0.8 | 444 | $2.1\times10^{-5}$ | 0.64 | 0.64 |
| 0.6 | 526 | $7.3\times10^{-6}$ | 0.74 | 0.79 |
| 0.4 (production) | 690 | $2.5\times10^{-6}$ | 1.00 | 1.00 |
| 0.3 | 851 | $1.7\times10^{-6}$ | 1.26 | 1.22 |
| 0.2 (reference) | 1197 | 0 | 1.81 | 1.76 |

On the linear states the convergence is a clean $h^4$: for $\mathrm{Im}\,\omega_{322}$ the changes from $h = 0.2$ are
$2.9\times10^{-4}$, $9.1\times10^{-5}$, $1.7\times10^{-5}$, $4.6\times10^{-6}$ at $h$ = 0.8, 0.6, 0.4, 0.3, which is
$C(h^4 - 0.2^4)$ with one constant $C$ to two digits. On the nonlinear states the errors at 0.4 and 0.3 are nearly the
same (plateau $2.2\times10^{-6}$ and $1.5\times10^{-6}$). The reason is that `step` = 0.2 also lowers the cap on the
stretched intervals from 1.86 to 1.6 (8 `step` becomes smaller than the quarter wavelength), so the reference differs
from the others also in the tail of the mesh; the tail resolution, not the uniform step, sets the remaining
difference. So the production error from the uniform step alone is about $10^{-6}$ (from the $h^4$ fit), and the rest
of the $2.5\times10^{-6}$ is the tail (see 3.2).

**Recommendation: keep 0.4.** Going to 0.3 costs 22-26% and removes less than half of the error; going to 0.6 saves
21-26% and triples it.

### 3.2 Outer endpoint (`containment`) and the tail of the mesh

| `containment` | $r_{*,\max}$ | nodes | rate error / fastest rate | cost 1 thr | cost 8 thr |
|---|---|---|---|---|---|
| $10^{-3}$ | 428 | 642 | $7.7\times10^{-4}$ | 0.90 | 0.92 |
| $10^{-4}$ (production) | 518 | 690 | $6.3\times10^{-4}$ | 1.00 | 1.00 |
| $10^{-5}$ | 605 | 737 | $4.7\times10^{-4}$ | 1.07 | 1.06 |
| $10^{-6}$ | 691 | 783 | $3.1\times10^{-4}$ | 1.09 | 1.09 |
| $10^{-8}$ (reference) | 860 | 873 | 0 | 1.25 | 1.28 |

The worst state is strong_takeover; on the other states the production error is $3.6\times10^{-6}$ (plateau),
$3.7\times10^{-6}$ (takeover), $4.2\times10^{-6}$ (strong_peak), $8.6\times10^{-6}$ (large_mesh, measured against
its $d\ln A_{211}/dt$) and below $10^{-8}$ on the linear $\alpha = 0.3$ state. **These numbers do not converge**: in
every state the change is proportional to the distance between the two endpoints (for the plateau about
$1.0\times10^{-8}$ per unit of $r_{*,\max}$ in $d\ln A_{322}/dt$; for strong_takeover about $2.5\times10^{-6}$ per unit
in $d\ln A_{322}/dt$) and shows no sign of levelling off at $r_{*,\max} = 860$. Refining this parameter would make the
numbers move further, not settle.

*Cause at the strong takeover (checked).* The component $q = (-1, 2)$, frequency $2\omega_{322} - \omega_{211} > \mu$,
propagates: it is radiation produced by the nonlinear coupling and travelling outwards. Its Klein-Gordon charge
$n_q$ inside the domain grows linearly with the endpoint, 0.042, 0.054, 0.066, 0.078, 0.102 for the five endpoints
(`../audit-work/radiation_in_box.jl`), because a steady outgoing wave fills the box. It enters the phase actions
$N_i = \sum_q q_i n_q$ with weights $-1$ (for $|211\rangle$) and $+2$ (for $|322\rangle$): radiation is 0.04% of
$N_{322}$ at production and 0.07% at $r_{*,\max} = 860$. The action balance treats $N_i$ as the cloud's action, so the
rates depend on how much radiation the box happens to contain. This is a definition question rather than a
numerical one: the radiation should probably be left out of the cloud actions (its flux through the outer boundary
is already in the loss terms). The owner should decide that before rerunning strongly nonlinear cases.

*At the weaker states the cause is not identified.* On the plateau the propagating components hold only $10^{-12}$ of
the actions, so radiation in the box is ruled out there; instead $n_{(1,0)}$ itself drifts down linearly, by
$4.6\times10^{-9}$ (relative) per unit of $r_{*,\max}$. A candidate is the discretisation error of the exponential
decay on the coarse stretched intervals (1.86 wide), which accumulates linearly with the length of the tail. The
decisive test, the same endpoint study with the stretched intervals capped at 0.8 (`outer_fine_tail`, job 27825217),
is *pending*: if the slope drops by about $(0.8/1.86)^4 \approx 0.03$, the tail discretisation is the cause and
`max_step` is the parameter to tighten.

**Recommendation:** keep `containment` = $10^{-4}$ for now, and report the rates with this endpoint error in mind:
about $4\times10^{-6}$ relative on the $f_a/M_\mathrm{Pl} = 0.047$ run and up to $10^{-3}$ at strong takeovers. Fix
the radiation accounting before relying on strongly nonlinear runs.

### 3.3 Inner endpoint (`horizon_offset_ratio`)

The radial functions are expanded in a Frobenius series at the inner endpoint; the series converges within the gap
$r_+ - r_-$, and `horizon_offset_ratio` says how far inside that radius the mesh starts.

| `horizon_offset_ratio` | nodes | rate error / fastest rate, five $f_a/M_\mathrm{Pl} = 0.047$ states | cost 1 thr | cost 8 thr |
|---|---|---|---|---|
| $10^{-4}$ | 652 | $1.1\times10^{-7}$ | 1.03 | 1.01 |
| $10^{-5}$ | 671 | $8.8\times10^{-8}$ | 1.01 | 1.00 |
| $10^{-6}$ (production) | 690 | $9.7\times10^{-8}$ | 1.00 | 1.00 |
| $10^{-8}$ (reference) | 727 | 0 | 1.08 | 0.99 |

The changes do not depend on the endpoint between $10^{-4}$ and $10^{-6}$; they are the change of the uniform mesh
itself (moving the inner endpoint shifts every node of the uniform region by a fraction of a step), not of the series.
Strong states: *pending*.

**Recommendation: keep.**

### 3.4 `stretch_factor`

| `stretch_factor` | nodes | rate error / fastest rate | cost 1 thr | cost 8 thr |
|---|---|---|---|---|
| 1.06 | 677 | $3.8\times10^{-7}$ | 1.04 | 1.08 |
| 1.03 (production) | 690 | $2.6\times10^{-7}$ | 1.00 | 1.00 |
| 1.015 (reference) | 715 | 0 | 1.09 | 1.06 |

(Five $f_a/M_\mathrm{Pl} = 0.047$ states; strong states *pending*.) **Recommendation: keep**; the factor mostly
changes how quickly the intervals reach the cap, and the cap is what matters (3.2).

## 4. Truncation

The field is expanded in Fourier components $q$ with $|q_1| + |q_2| \le$ `max_order` (odd orders only, for the even
cosine potential) and, per component, the `n_angular` lowest polar indices $l$ of the right parity. Leaving out
higher components or polar indices is a *truncation error*: it does not shrink when the mesh or quadrature is refined.

| setting | channels | rate error / fastest rate, against `n_angular` = 2 | cost 1 thr | cost 8 thr |
|---|---|---|---|---|
| production (`max_order` 3, `n_angular` 1) | 8 | $7.3\times10^{-4}$ (plateau), $6.1\times10^{-4}$ (takeover), $1.3\times10^{-5}$ (nonlinear), $7.5\times10^{-8}$ (linear) | 1.00 | 1.00 |
| `max_order` 5 | 18 | change from production $\le 2.7\times10^{-6}$ | 2.89 | 3.31 |
| `n_angular` 2 | 16 | reference | 1.91 | 2.48 |

`max_order` = 5 also enlarges the mesh (826 nodes: `domain_estimate` includes the new components), so its
$2.7\times10^{-6}$ includes an outer-endpoint change of the size found in 3.2; the Fourier truncation itself is at
most that. The polar truncation is the largest error in this audit. The large_mesh state shows $9\times10^{-5}$ with
`n_angular` = 2 although it is linear; at $\alpha = 0.187$ that would have to be the Kerr coupling between $l$ and
$l + 2$, and I have not checked why it is so much larger than at $\alpha = 0.3$ ($7.5\times10^{-8}$). Strong states:
*pending*. Whether `n_angular` = 2 is converged would need `n_angular` = 3, which was not run.

**Recommendation:** the owner should decide whether $7\times10^{-4}$ in the rates matters for the physics questions
of the reruns. If it does, run production with `n_angular` = 2 (about twice the cost per evaluation, 2.5 times at 8
threads) and first check `n_angular` = 3 at the plateau and takeover states. Keep `max_order` = 3.

## 5. Tolerances

### 5.1 Newton tolerance

`newton_tol` is the scaled residual (and step) at which each radial solve stops; it controls how far each solve is
from the exact discrete solution. Here the perturbed right-hand-side evaluation was repeated with `newton_tol` (and
the tolerance of all radial solves inside it) from $10^{-7}$ to $10^{-12}$, at production discretisation:

| state | largest relative change of the four rates, $10^{-7}$ vs $10^{-12}$ | Newton iterations ($10^{-7}$ / $10^{-10}$ / $10^{-12}$) | time 1 thr ($10^{-7}$ / $10^{-10}$) |
|---|---|---|---|
| linear | $8\times10^{-12}$ | 4 / 5 / 6 | 0.352 / 0.347 s |
| plateau | $1.8\times10^{-9}$ | 3 / 4 / 4 | 0.307 / 0.357 s |
| takeover | $2.1\times10^{-9}$ | 3 / 5 / 6 | 0.317 / 0.359 s |
| nonlinear | $1\times10^{-11}$ | 4 / 5 / 9 | 0.337 / 0.358 s |
| large_mesh | $0$ | 4 / 7 / 7 | 0.837 / 1.195 s |

Newton's method converges quadratically, so the last iteration usually lands far below the tolerance: at $10^{-7}$ the
final residuals were $4\times10^{-12}$ to $2\times10^{-10}$. Rates change by at most $2\times10^{-9}$, three orders below
the radial-step error. Loosening to $10^{-8}$ saves at most one iteration (0-14%, within the noise on most states).

On the $\alpha = 0.187$ mesh the residual cannot go below $1.7\times10^{-9}$ after the harness's small perturbation
(the documented floor of $5\times10^{-10}$ is that of the stored state itself; it ranges from $2\times10^{-10}$ to
$1\times10^{-9}$ with the step, see the step table's large_mesh rows). With the default `newton_tol` = $10^{-10}$ and
`newton_floor_tol` = $10^{-9}$ every perturbed solve there fails (the first audit run failed all its large_mesh timing
calls this way); the audit used `newton_floor_tol` = $3\times10^{-9}$. With `newton_tol` = $10^{-8}$ the solve stops
after 4 instead of 7 iterations with the same rates, 28% faster.

**Recommendation:** keep `newton_tol` = $10^{-10}$ on the $\alpha = 0.3$ meshes (the saving from loosening is within
noise). For $\alpha = 0.187$ and other large meshes use `newton_tol` = $10^{-8}$ (or keep $10^{-10}$ with
`newton_floor_tol` of at least $3\times10^{-9}$; $10^{-9}$ is not enough).

### 5.2 Integrator tolerance

The work-precision study of the methods agent (`methods-work/REPORT_2026-11-13_1115.md`) measured the trajectory
error of Tsit5, FBDF, Rodas5P and Rosenbrock23 at `reltol` = `abstol` from $10^{-6}$ to $10^{-10}$ on the
$f_a/M_\mathrm{Pl} = 0.047$ window $1.5\times10^9 \to 10^{11}$: Tsit5 at $10^{-7}$ gives a largest error over the
checkpoints of $2.3\times10^{-7}$ (in $\ln A$, $b$, $j$) for 2300 radial solves; FBDF at $10^{-10}$ gives
$7.7\times10^{-8}$ for 1144 solves, and FBDF at $10^{-9}$ is about 10 times cheaper than Tsit5 on stiff plateaus. Set
against this audit: the radial discretisation alone puts relative errors of about $2.5\times10^{-6}$ into the rates
(more from the outer endpoint and the polar truncation), which over 10 e-folds is a few times $10^{-5}$ in the
amplitudes. The integrator at $10^{-7}$ is already 100 times more accurate than the right-hand side it integrates.

**Recommendation:** do not tighten below $10^{-7}$ for accuracy; it only buys agreement with the same discrete
equations. The methods report's method choice stands (Tsit5 at $10^{-7}$ where the problem is not stiff, FBDF at
$10^{-9}$ to $10^{-10}$ on stiff plateaus).

## 6. All finer at once

| state | rate error of production against (n_phase 24, n_polar 16, step 0.3, stretch 1.015, containment $10^{-6}$, horizon $10^{-8}$) | cost of the finer setting 1 thr / 8 thr |
|---|---|---|
| linear | $5.1\times10^{-8}$ | 3.27 / 1.73 |
| plateau | $8.9\times10^{-7}$ | 3.42 / 2.26 |
| takeover | $9.0\times10^{-7}$ | 3.57 / 2.12 |
| nonlinear | $1.1\times10^{-7}$ | 3.37 / 2.20 |
| large_mesh | $4.3\times10^{-6}$ | (8 thr) 2.16 |
| strong_peak | $2.2\times10^{-5}$ | (8 thr) 2.32 |
| strong_takeover | $3.2\times10^{-4}$ | (8 thr) 2.60 |

The joint change is of the size of the larger single-parameter changes (step and outer endpoint), as expected if the
errors add; the strong takeover is dominated by the outer endpoint (3.2). This is one joint point, not a joint
convergence study, and it does not include the truncation.

## Appendix: full tables

Generated with `julia validation/parameter_audit_tables.jl <study> ../audit-work/run1 ../audit-work/run2_large
../audit-work/run4_strong ../audit-work/run5_lattice` (later directories replace a state's rows from earlier ones).
Settings are labelled by what differs from production: `np` `n_phase`, `npol` `n_polar`, `h` step, `sf`
stretch_factor, `c` containment, `hr` horizon_offset_ratio, `mo` max_order, `na` n_angular, `lat` lattice points,
`ms` max_step. "NaN" in a cost column: that timing was not taken (the 1-thread pass of the strong-state job had not
run when the table was made).

