# Numerical-parameter audit of AxionCloud.jl (2026-11-13, completed 2026-11-14)

Written for the owner before the production runs are repeated with the sped-up code. For each numerical parameter
this document says what it controls, how much the quantities that drive an evolution change when it is refined, what
one evaluation costs at each setting, and a recommendation.

**Code version.** All numbers below were computed on `speedup/integration-next` at ebf0ae9 (plus the lattice commit
b3bcdb4), that is, *before* the conjugate-frequency rule of commit 41fc2e7 on `rewrite/axioncloud`. That commit gives
a driven component $q$ the complex frequency $\sum_i q_i \omega_i + i \sum_i |q_i| \gamma_i$ (with
$	ilde\omega_i = \omega_i + i\gamma_i$ the complex frequency of bound mode $i$) instead of $q\cdot	ilde\omega$, so
components that contain a conjugated mode ($q_i < 0$, such as $(2,-1)$ and $(-1,2)$) now grow with the modes instead
of partly decaying. It changes only the imaginary parts of those components' frequencies, which are of the order of
the growth rates. I expect it to leave the quadrature conclusions (Sections 1, 2), the radial-step and stretch
conclusions (3.1, 3.4) and the tolerance conclusions (5) unchanged, since those concern how accurately the same
field is averaged and integrated. The conclusions that involve the $(2,-1)$ and $(-1,2)$ components directly should
be rechecked: the endpoint dependence at the strong states (3.2, 3.3), which comes from exactly those components'
waves in transit, and the size of the polar-truncation error at the strong states (4). Section 7 gives a spot check
on the merged code.

## If you read only this

| Parameter | Production | Error at production (rate error / fastest rate, worst state) | Recommendation | Cost change per evaluation (1 / 8 threads) |
|---|---|---|---|---|
| `n_phase` (tensor phase grid) | 16 | $3\times10^{-15}$ | loosen to 12 | $-26\%$ / $-15\%$ |
| rank-1 lattice instead of the tensor grid | (not used) | 85 points: $2\times10^{-13}$ | adopt the 85-point lattice if the new code path is accepted (replaces the line above) | $-42\%$ / $-19\%$ |
| `n_polar` | 12 | $2\times10^{-15}$ | loosen to 8 | $-22\%$ / $0\%$ |
| **85-point lattice and `n_polar` 8 together (measured)** | | $1.8\times10^{-13}$ | **the recommended quadrature** | $-49\%$ / $-20\%$ |
| radial `step` | 0.4 | $2.5\times10^{-6}$ (about $10^{-6}$ from the uniform step, the rest from the tail cap) | keep | |
| cap on stretched intervals (`max_step`) | 1.86 | about $3\times10^{-6}$ (weaker states) | keep, unless the polar truncation is raised; 0.8 removes it | $+45\%$ at 8 threads (1 thread not measured) |
| outer endpoint (`containment`) | $10^{-4}$ | $6\times10^{-4}$ (strong takeover), $4\times10^{-6}$ elsewhere; **does not converge** | not a parameter fix: see Section 3.2 | |
| inner endpoint (`horizon_offset_ratio`) | $10^{-6}$ | $\le 3\times10^{-8}$ on the weaker states, $2\times10^{-5}$ on the strong ones; **does not converge** (same cause as the outer endpoint) | keep; see Section 3.3 | |
| `stretch_factor` | 1.03 | $3\times10^{-7}$ | keep | |
| `n_angular` (polar indices per component) | 1 | $2\times10^{-2}$ (strong states), $7\times10^{-4}$ (plateau, takeover) | the largest error found; run with 2 at least for strongly nonlinear cases, after checking 3 | $+91\%$ / $+148\%$ |
| `max_order` (Fourier truncation) | 3 | $\lesssim 3\times10^{-6}$ (mixed with the outer endpoint) | keep | |
| `newton_tol` | $10^{-10}$ | $\le 2\times10^{-9}$ even at $10^{-7}$ | keep on the $\alpha = 0.3$ meshes; use $10^{-8}$ on the $\alpha = 0.187$ mesh | $0$ to $-14\%$; $-28\%$ at $\alpha = 0.187$ |
| integrator `reltol`/`abstol` | $10^{-7}$ | see the methods report | keep $10^{-7}$ (Tsit5) or FBDF at $10^{-9}$-$10^{-10}$ on stiff plateaus | |

The two quadrature changes together were run as one setting at all seven states (Section 1.3): the 85-point lattice
with `n_polar` = 8 agrees with a 32 x 32 tensor grid with 24 polar nodes to $1.8\times10^{-13}$ and costs 0.51 of
production at 1 thread and 0.80 at 8 threads.

Three findings matter more than any saving:

1. **The angular truncation is the largest error.** Allowing a second polar index per Fourier component changes
   $d\ln A_{322}/dt$ by 2% at the strong $|211\rangle$ peak, 1.3% at the strong takeover, $7\times10^{-4}$ on the
   plateau and $6\times10^{-4}$ at the takeover. That is $10^{3}$-$10^{4}$ times the radial-step error. Whether
   `n_angular` = 2 is itself converged is not known (3 was not run).
2. **The endpoint dependence does not converge, at either end.** Moving either endpoint changes the rates in
   proportion to the change of the domain length in $r_*$, without levelling off. The cause is identified for the
   strong states: the phase actions $N_i$ are integrals over the domain, and they include wave action that is only
   passing through, outgoing radiation near the outer endpoint (component $q = (-1,2)$) and waves falling into the
   horizon near the inner endpoint (component $q = (2,-1)$ and the $|211\rangle$ mode itself). That part grows
   linearly with the domain. It is a question of definition in the action balance, not of resolution.
3. **On the weaker states the outer-endpoint drift is a resolution error of the mesh tail.** Capping the stretched
   intervals at 0.8 instead of 1.86 reduces it about 30-fold, as expected for a fourth-order scheme (Section 3.2).
4. **The integrator tolerance is already tighter than the discretisation.** The radial step alone gives relative rate
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
tolerance), `run4_strong` (strong states), `run5_lattice`, `run6_tail` (finer tail cap), `run7_recommended`
(recommended quadrature); `radiation_in_box.jl` and `horizon_waves_in_box.jl` there compute the per-component charges
used in 3.2 and 3.3. Slurm jobs 27824736-27824738, 27824847, 27824972, 27825063, 27825217, 27825310. The library defaults are unchanged. The branch also carries commit b3bcdb4 (the
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

### 1.3 The recommended quadrature as one setting

The 85-point lattice and `n_polar` = 8 were run together at all seven states (study `recommended`, job 27825310),
against a reference that is finer in both (32 x 32 tensor grid, 24 polar nodes):

| setting | rate error / fastest rate | $db/dt$, $dj/dt$ | cost 1 thr | cost 8 thr |
|---|---|---|---|---|
| lattice 85, `n_polar` 8 | $1.8\times10^{-13}$ | $4.2\times10^{-15}$ | 0.51 | 0.80 |
| production | $5.5\times10^{-15}$ | $1.3\times10^{-15}$ | 1.00 | 1.00 |
| tensor 32, `n_polar` 24 (reference) | 0 | 0 | 5.63 | 3.16 |

The largest change is at strong_peak; every other state is below $4\times10^{-14}$. The plateau segment ends within
$4\times10^{-14}$ of the reference in $\ln A_{322}$, with the same 31 calls. The savings of the two changes nearly add
at 1 thread ($-49\%$ against $-42\%$ and $-22\%$ separately); at 8 threads the combination saves 20%, about what the
lattice alone saves, because the polar average is not a significant cost there.

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

*At the weaker states the cause is the resolution of the mesh tail (checked).* On the plateau the propagating
components hold only $10^{-12}$ of the actions, so radiation in the box is ruled out there; instead $n_{(1,0)}$
itself drifts down linearly, by $5\times10^{-9}$ (relative) per unit of $r_{*,\max}$. The candidate was the
discretisation error of the exponential decay on the coarse stretched intervals (1.86 wide), which accumulates
linearly with the length of the tail. The decisive test repeats the endpoint study with the stretched intervals
capped at `max_step` = 0.8 (study `outer_fine_tail`, job 27825217); for a fourth-order scheme the drift should fall by
about $(0.8/1.86)^4 \approx 0.03$. Rate error at `containment` = $10^{-4}$ against $10^{-8}$:

| state | default cap 1.86 | cap 0.8 | ratio |
|---|---|---|---|
| plateau | $3.6\times10^{-6}$ | $1.2\times10^{-7}$ | 0.033 |
| takeover | $3.7\times10^{-6}$ | $8.7\times10^{-8}$ | 0.024 |
| strong_peak | $4.2\times10^{-6}$ | $2.9\times10^{-7}$ | 0.07 |
| strong_takeover | $6.3\times10^{-4}$ | $8.6\times10^{-4}$ | (radiation, unchanged) |
| large_mesh | $8.6\times10^{-6}$ | $8.6\times10^{-6}$ | (unchanged, not explained) |

On the plateau and the takeover the drift falls by the predicted factor, so it is tail resolution. At the same
endpoint the cap itself changes $d\ln A_{322}/dt$ by $2.8\times10^{-6}$ (plateau), $2.9\times10^{-6}$ (takeover) and
$3.1\times10^{-6}$ (strong_peak): this is the part of the step study's error that did not follow $h^4$ (3.1). With the
finer cap the mesh has 941 nodes instead of 690, and one evaluation took 0.20-0.22 s at 8 threads against about 0.15 s
for production in the other jobs (about $+45\%$; different jobs and nodes, and no 1-thread timing). The strong
takeover is unchanged because its drift is the radiation counted in the box (above). The large_mesh drift
($\alpha = 0.187$, a linear state, measured on $d\ln A_{211}/dt$) is also unchanged, so it is neither tail resolution
nor nonlinear radiation; I have not found its cause.

**Recommendation:** keep `containment` = $10^{-4}$. Keep the default cap for now: with it the weaker states carry
about $3\times10^{-6}$ (cap) plus $4\times10^{-6}$ (endpoint) in the rates, and `max_step` = 0.8 would remove both
for about $+45\%$ per evaluation. That is worth paying only together with `n_angular` = 2, since the polar
truncation error ($7\times10^{-4}$) is 100 times larger. Fix the radiation accounting before relying on strongly
nonlinear runs.

### 3.3 Inner endpoint (`horizon_offset_ratio`)

The radial functions are expanded in a Frobenius series at the inner endpoint; the series converges within the gap
$r_+ - r_-$, and `horizon_offset_ratio` says how far inside that radius the mesh starts.

| `horizon_offset_ratio` | $r_{*,\min}$ | nodes | rate error / fastest rate, weaker five states | strong_peak | strong_takeover | cost 1 thr | cost 8 thr |
|---|---|---|---|---|---|---|---|
| $10^{-4}$ | $-29.2$ | 652 | $3.3\times10^{-8}$ | $4.6\times10^{-5}$ | $1.6\times10^{-5}$ | 1.03 | 1.01 |
| $10^{-5}$ | $-36.8$ | 671 | $2.4\times10^{-8}$ | $3.4\times10^{-5}$ | $1.2\times10^{-5}$ | 1.01 | 1.00 |
| $10^{-6}$ (production) | $-44.3$ | 690 | $1.6\times10^{-8}$ | $2.3\times10^{-5}$ | $7.9\times10^{-6}$ | 1.00 | 1.00 |
| $10^{-8}$ (reference) | $-59.5$ | 727 | 0 | 0 | 0 | 1.08 | 0.99 |

(Rate errors of the strong states are $d\ln A_{322}/dt$, the faster rate there. The weaker states' column is the
largest relative change of either amplitude rate.) For the weaker states the inner endpoint does not matter. For the
strong states the change is proportional to the shift of $r_{*,\min}$ (each decade of the offset moves $r_{*,\min}$
by about 7.6 and the rate by about $1.2\times10^{-5}$) and does not level off. The cause is checked
(`../audit-work/horizon_waves_in_box.jl` on the saved strong_peak states): the Klein-Gordon charge of the component
$q = (2,-1)$ in the domain grows linearly with $|r_{*,\min}|$, $2.585$, $2.638$, $2.692$, $2.800\times10^{-4}$, about
$7\times10^{-7}$ per unit of $r_*$. Near the horizon each component is a wave falling in, whose action density per
unit $r_*$ is constant, so the integral over the domain grows with its length. That component enters $N_{322}$ with
weight $-1$; the $|211\rangle$ mode's own charge also drifts ($-2.9\times10^{-5}$ per unit of $r_*$; it is
superradiant, so its infalling wave carries negative action). This is the same accounting question as in 3.2.

**Recommendation: keep $10^{-6}$**, and treat the inner and outer endpoint dependence together when the definition of
the cloud actions is revisited.

### 3.4 `stretch_factor`

| `stretch_factor` | nodes | rate error / fastest rate | cost 1 thr | cost 8 thr |
|---|---|---|---|---|
| 1.06 | 677 | $3.8\times10^{-7}$ | 1.04 | 1.08 |
| 1.03 (production) | 690 | $2.6\times10^{-7}$ | 1.00 | 1.00 |
| 1.015 (reference) | 715 | 0 | 1.09 | 1.06 |

(All seven states; the largest is strong_takeover, $3.6\times10^{-7}$ at production.) **Recommendation: keep**; the factor mostly
changes how quickly the intervals reach the cap, and the cap is what matters (3.2).

## 4. Truncation

The field is expanded in Fourier components $q$ with $|q_1| + |q_2| \le$ `max_order` (odd orders only, for the even
cosine potential) and, per component, the `n_angular` lowest polar indices $l$ of the right parity. Leaving out
higher components or polar indices is a *truncation error*: it does not shrink when the mesh or quadrature is refined.

| setting | channels | rate error / fastest rate, against `n_angular` = 2 | cost 1 thr | cost 8 thr |
|---|---|---|---|---|
| production (`max_order` 3, `n_angular` 1) | 8 | $2.0\times10^{-2}$ (strong_peak), $1.3\times10^{-2}$ (strong_takeover), $7.3\times10^{-4}$ (plateau), $6.1\times10^{-4}$ (takeover), $1.3\times10^{-5}$ (nonlinear), $7.5\times10^{-8}$ (linear) | 1.00 | 1.00 |
| `max_order` 5 | 18 | change from production $\le 2.7\times10^{-6}$ | 2.89 | 3.31 |
| `n_angular` 2 | 16 | reference | 1.91 | 2.48 |

`max_order` = 5 also enlarges the mesh (826 nodes: `domain_estimate` includes the new components), so its
$2.7\times10^{-6}$ includes an outer-endpoint change of the size found in 3.2; the Fourier truncation itself is at
most that. The polar truncation is the largest error in this audit. The large_mesh state shows $9\times10^{-5}$ with
`n_angular` = 2 although it is linear; at $\alpha = 0.187$ that would have to be the Kerr coupling between $l$ and
$l + 2$, and I have not checked why it is so much larger than at $\alpha = 0.3$ ($7.5\times10^{-8}$). At the strong
states `max_order` = 5 changes nothing beyond $10^{-6}$ while `n_angular` = 2 changes $d\ln A_{322}/dt$ by 1-2% and
$db/dt$ by up to 1.6%: the error grows roughly like the square of the field amplitude. Whether `n_angular` = 2 is converged would need `n_angular` = 3, which was not run.

**Recommendation:** at strong nonlinearity (the $f_a/M_\mathrm{Pl} = 3\times10^{-3}$ runs) a 1-2% error in the rates
is not acceptable for quantitative results; run those with `n_angular` = 2 after checking `n_angular` = 3 at the
strong states. For the $f_a/M_\mathrm{Pl} = 0.047$ runs the owner should decide whether $7\times10^{-4}$ matters. If it
does, run production with `n_angular` = 2 (about twice the cost per evaluation, 2.5 times at 8
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
| nonlinear | $4\times10^{-11}$ | 4 / 5 / 9 | 0.337 / 0.358 s |
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
../audit-work/run4_strong ../audit-work/run5_lattice ../audit-work/run6_tail ../audit-work/run7_recommended` (later
directories replace a state's rows from earlier ones; each study lives in one directory or is split by state).
Settings are labelled by what differs from production: `np` `n_phase`, `npol` `n_polar`, `h` step, `sf`
stretch_factor, `c` containment, `hr` horizon_offset_ratio, `mo` max_order, `na` n_angular, `lat` lattice points,
`ms` max_step. "NaN" in a cost column: that timing was not taken (the `outer_fine_tail` job ran no 1-thread pass, and
its 8-thread cost ratio is NaN because that job did not time the production setting).


### Study `n_phase`

| state | setting | nodes | Re ω (max) | Im ω₂₁₁ | Im ω₃₂₂ | d ln A₂₁₁/dt | d ln A₃₂₂/dt | db/dt | dj/dt | amplitude-rate error / fastest rate | t(1 thr) s | t(8 thr) s | cost vs prod (1/8 thr) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| linear | np8 | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.212 | 0.146 | 0.55 / 0.83 |
| linear | np12 | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.277 | 0.145 | 0.72 / 0.83 |
| linear | production | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.386 | 0.175 | 1.00 / 1.00 |
| linear | np24 | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.675 | 0.248 | 1.75 / 1.41 |
| linear | np32 (ref) | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1.129 | 0.363 | 2.92 / 2.07 |
| plateau | np8 | 690 | 3.8e-16 | 1.4e-07 | 1.3e-10 | 5.5e-08 | 1.8e-10 | 5.4e-08 | 5.3e-08 | 1.8e-10 | 0.197 | 0.152 | 0.55 / 0.93 |
| plateau | np12 | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.265 | 0.135 | 0.74 / 0.82 |
| plateau | production | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.358 | 0.164 | 1.00 / 1.00 |
| plateau | np24 | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.620 | 0.228 | 1.73 / 1.39 |
| plateau | np32 (ref) | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1.028 | 0.306 | 2.87 / 1.86 |
| takeover | np8 | 690 | 0 | 6.3e-11 | 1.3e-10 | 6.9e-11 | 1.3e-10 | 4.4e-11 | 3.0e-12 | 1.3e-10 | 0.196 | 0.159 | 0.54 / 0.93 |
| takeover | np12 | 690 | 0 | 1.1e-11 | 1.2e-13 | 3.3e-15 | 0 | 3.8e-16 | 6.0e-16 | 7.2e-16 | 0.281 | 0.156 | 0.77 / 0.91 |
| takeover | production | 690 | 0 | 1.2e-12 | 1.1e-13 | 9.8e-16 | 0 | 5.1e-16 | 6.0e-16 | 2.1e-16 | 0.367 | 0.171 | 1.00 / 1.00 |
| takeover | np24 | 690 | 0 | 4.5e-12 | 1.7e-13 | 3.9e-16 | 0 | 2.5e-16 | 0 | 8.4e-17 | 0.678 | 0.244 | 1.85 / 1.43 |
| takeover | np32 (ref) | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1.040 | 0.350 | 2.84 / 2.04 |
| nonlinear | np8 | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.189 | 0.138 | 0.49 / 0.81 |
| nonlinear | np12 | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.277 | 0.156 | 0.72 / 0.92 |
| nonlinear | production | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.385 | 0.170 | 1.00 / 1.00 |
| nonlinear | np24 | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.687 | 0.246 | 1.79 / 1.45 |
| nonlinear | np32 (ref) | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1.070 | 0.343 | 2.78 / 2.02 |
| large_mesh | np8 | 1355 | 0 | 3.3e-13 | 1.6e-09 | 2.3e-13 | 5.6e-13 | 6.2e-16 | 1.0e-15 | 2.3e-13 | 0.656 | 0.387 | 0.53 / 0.78 |
| large_mesh | np12 | 1355 | 0 | 0 | 0 | 0 | 3.5e-16 | 0 | 0 | 1.5e-21 | 0.950 | 0.381 | 0.77 / 0.77 |
| large_mesh | production | 1355 | 0 | 0 | 0 | 0 | 1.7e-16 | 0 | 0 | 7.5e-22 | 1.230 | 0.496 | 1.00 / 1.00 |
| large_mesh | np24 | 1355 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2.222 | 0.619 | 1.81 / 1.25 |
| large_mesh | np32 (ref) | 1355 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 3.563 | 0.854 | 2.90 / 1.72 |
| strong_peak | np8 | 690 | 1.5e-11 | 4.2e-08 | 2.8e-07 | 1.9e-08 | 2.8e-07 | 1.2e-08 | 4.3e-09 | 2.8e-07 | 0.209 | 0.125 | 0.60 / 0.82 |
| strong_peak | np12 | 690 | 0 | 9.4e-13 | 2.6e-13 | 7.1e-15 | 1.5e-14 | 1.8e-16 | 2.0e-16 | 1.5e-14 | 0.251 | 0.127 | 0.72 / 0.83 |
| strong_peak | production | 690 | 0 | 2.8e-13 | 4.8e-15 | 5.9e-16 | 9.7e-16 | 0 | 2.0e-16 | 9.7e-16 | 0.347 | 0.152 | 1.00 / 1.00 |
| strong_peak | np24 | 690 | 0 | 4.6e-13 | 1.3e-13 | 1.2e-15 | 2.7e-14 | 5.3e-16 | 2.0e-16 | 2.7e-14 | 0.615 | 0.230 | 1.77 / 1.51 |
| strong_peak | np32 (ref) | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1.008 | 0.322 | 2.90 / 2.11 |
| strong_takeover | np8 | 690 | 3.4e-12 | 1.2e-07 | 1.5e-07 | 1.3e-07 | 1.6e-07 | 1.3e-07 | 1.0e-09 | 1.3e-07 | 0.196 | 0.144 | 0.53 / 0.85 |
| strong_takeover | np12 | 690 | 0 | 1.2e-13 | 2.0e-14 | 1.8e-15 | 4.1e-15 | 1.3e-15 | 7.7e-16 | 3.0e-15 | 0.268 | 0.144 | 0.72 / 0.85 |
| strong_takeover | production | 690 | 0 | 1.4e-14 | 2.6e-15 | 9.7e-16 | 3.5e-15 | 6.7e-16 | 0 | 2.5e-15 | 0.370 | 0.170 | 1.00 / 1.00 |
| strong_takeover | np24 | 690 | 0 | 5.2e-14 | 2.4e-14 | 8.5e-16 | 1.6e-16 | 1.2e-15 | 0 | 8.5e-16 | 0.717 | 0.234 | 1.93 / 1.38 |
| strong_takeover | np32 (ref) | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1.072 | 0.335 | 2.89 / 1.97 |

Summary over states (largest change over the states listed; cost ratio to production, mean over those states):

| setting | states | amplitude-rate error / fastest rate | db/dt, dj/dt (relative) | cost 1 thr | cost 8 thr |
|---|---|---|---|---|---|
| np8 | 7 | 2.8e-07 | 1.3e-07 | 0.54 | 0.85 |
| np12 | 7 | 1.5e-14 | 1.3e-15 | 0.74 | 0.85 |
| production | 7 | 2.5e-15 | 6.7e-16 | 1.00 | 1.00 |
| np24 | 7 | 2.7e-14 | 1.2e-15 | 1.80 | 1.40 |
| np32 | 7 | 0 | 0 | 2.87 | 1.97 |

Segment 1.5e9 → 1.9e9 from the plateau state (relative change of the final state from the reference setting):

| setting | rhs calls | ln A₂₁₁ (abs) | ln A₃₂₂ (abs) | b (rel) | j (rel) | wall s (8 thr) |
|---|---|---|---|---|---|---|
| np8 | 31 | 3.6e-15 | 1.7e-11 | 1.1e-16 | 6.6e-16 | 9.05 |
| np12 | 31 | 0.0e+00 | 2.1e-14 | 0.0e+00 | 0.0e+00 | 3.07 |
| production | 31 | 0.0e+00 | 1.1e-14 | 0.0e+00 | 0.0e+00 | 3.47 |
| np24 | 31 | 0.0e+00 | 1.1e-14 | 0.0e+00 | 0.0e+00 | 4.81 |
| np32 | 31 | 0.0e+00 | 0.0e+00 | 0.0e+00 | 0.0e+00 | 6.74 |

### Study `lattice`

| state | setting | nodes | Re ω (max) | Im ω₂₁₁ | Im ω₃₂₂ | d ln A₂₁₁/dt | d ln A₃₂₂/dt | db/dt | dj/dt | amplitude-rate error / fastest rate | t(1 thr) s | t(8 thr) s | cost vs prod (1/8 thr) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| linear | lat31 | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.179 | 0.142 | 0.47 / 0.91 |
| linear | lat61 | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.194 | 0.132 | 0.51 / 0.85 |
| linear | lat85 | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.216 | 0.142 | 0.57 / 0.91 |
| linear | lat101 | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.232 | 0.139 | 0.61 / 0.89 |
| linear | lat151 | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.286 | 0.139 | 0.75 / 0.89 |
| linear | production | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.382 | 0.156 | 1.00 / 1.00 |
| linear | np32 (ref) | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1.066 | 0.321 | 2.79 / 2.06 |
| plateau | lat31 | 690 | 1.1e-15 | 1.9e-07 | 7.5e-07 | 1.2e-08 | 3.3e-07 | 3.3e-09 | 2.2e-09 | 3.3e-07 | 0.165 | 0.125 | 0.48 / 0.84 |
| plateau | lat61 | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.187 | 0.144 | 0.54 / 0.97 |
| plateau | lat85 | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.208 | 0.132 | 0.60 / 0.89 |
| plateau | lat101 | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.218 | 0.144 | 0.63 / 0.97 |
| plateau | lat151 | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.268 | 0.141 | 0.77 / 0.95 |
| plateau | production | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.346 | 0.148 | 1.00 / 1.00 |
| plateau | np32 (ref) | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.997 | 0.314 | 2.88 / 2.12 |
| takeover | lat31 | 690 | 7.2e-14 | 4.9e-06 | 1.6e-04 | 1.1e-04 | 1.1e-04 | 8.7e-05 | 3.0e-10 | 1.1e-04 | 0.166 | 0.142 | 0.43 / 0.86 |
| takeover | lat61 | 690 | 0 | 2.2e-12 | 3.3e-13 | 6.5e-15 | 3.9e-15 | 2.7e-15 | 1.8e-15 | 3.9e-15 | 0.194 | 0.145 | 0.50 / 0.88 |
| takeover | lat85 | 690 | 0 | 1.1e-12 | 4.6e-13 | 2.2e-15 | 2.5e-15 | 1.3e-15 | 0 | 2.5e-15 | 0.216 | 0.134 | 0.56 / 0.81 |
| takeover | lat101 | 690 | 0 | 2.2e-12 | 8.3e-14 | 5.9e-16 | 2.0e-15 | 5.1e-16 | 0 | 2.0e-15 | 0.232 | 0.138 | 0.60 / 0.83 |
| takeover | lat151 | 690 | 0 | 5.1e-13 | 6.5e-13 | 2.6e-15 | 6.8e-16 | 1.3e-16 | 0 | 6.8e-16 | 0.286 | 0.150 | 0.74 / 0.90 |
| takeover | production | 690 | 0 | 1.2e-12 | 1.1e-13 | 9.8e-16 | 0 | 5.1e-16 | 6.0e-16 | 2.1e-16 | 0.388 | 0.166 | 1.00 / 1.00 |
| takeover | np32 (ref) | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1.083 | 0.335 | 2.79 / 2.03 |
| nonlinear | lat31 | 690 | 1.3e-15 | 1.4e-09 | 2.0e-08 | 3.3e-08 | 1.4e-08 | 5.6e-08 | 3.9e-12 | 3.3e-08 | 0.166 | 0.143 | 0.45 / 0.80 |
| nonlinear | lat61 | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.194 | 0.133 | 0.52 / 0.74 |
| nonlinear | lat85 | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.218 | 0.131 | 0.59 / 0.73 |
| nonlinear | lat101 | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.232 | 0.141 | 0.62 / 0.78 |
| nonlinear | lat151 | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.288 | 0.148 | 0.77 / 0.82 |
| nonlinear | production | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.372 | 0.180 | 1.00 / 1.00 |
| nonlinear | np32 (ref) | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1.094 | 0.329 | 2.94 / 1.83 |
| large_mesh | lat31 | 1355 | 0 | 0 | 0 | 0 | 1.7e-16 | 0 | 0 | 7.5e-22 | 0.537 | 0.395 | 0.46 / 0.77 |
| large_mesh | lat61 | 1355 | 0 | 0 | 0 | 0 | 3.5e-16 | 0 | 0 | 1.5e-21 | 0.635 | 0.411 | 0.55 / 0.80 |
| large_mesh | lat85 | 1355 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.703 | 0.390 | 0.60 / 0.76 |
| large_mesh | lat101 | 1355 | 0 | 0 | 0 | 0 | 1.7e-16 | 0 | 0 | 7.5e-22 | 0.821 | 0.435 | 0.71 / 0.85 |
| large_mesh | lat151 | 1355 | 0 | 0 | 0 | 0 | 1.7e-16 | 0 | 0 | 7.5e-22 | 0.910 | 0.455 | 0.78 / 0.89 |
| large_mesh | production | 1355 | 0 | 0 | 0 | 0 | 1.7e-16 | 0 | 0 | 7.5e-22 | 1.163 | 0.512 | 1.00 / 1.00 |
| large_mesh | np32 (ref) | 1355 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 3.427 | 0.942 | 2.95 / 1.84 |
| strong_peak | lat31 | 690 | 1.9e-09 | 4.6e-05 | 1.3e-03 | 4.9e-05 | 8.6e-04 | 2.6e-05 | 1.9e-09 | 8.6e-04 | 0.157 | 0.126 | 0.45 / 0.75 |
| strong_peak | lat61 | 690 | 1.9e-16 | 1.1e-10 | 6.1e-10 | 6.2e-11 | 1.2e-09 | 3.1e-11 | 6.1e-16 | 1.2e-09 | 0.185 | 0.133 | 0.53 / 0.79 |
| strong_peak | lat85 | 690 | 0 | 2.9e-13 | 2.0e-13 | 1.0e-14 | 2.2e-13 | 5.3e-15 | 2.0e-16 | 2.2e-13 | 0.209 | 0.123 | 0.59 / 0.73 |
| strong_peak | lat101 | 690 | 0 | 9.4e-13 | 1.3e-13 | 1.2e-15 | 8.7e-15 | 5.3e-16 | 2.0e-16 | 8.7e-15 | 0.220 | 0.131 | 0.62 / 0.78 |
| strong_peak | lat151 | 690 | 0 | 1.9e-13 | 5.8e-14 | 5.9e-16 | 1.8e-15 | 1.8e-16 | 2.0e-16 | 1.8e-15 | 0.267 | 0.127 | 0.76 / 0.75 |
| strong_peak | production | 690 | 0 | 2.8e-13 | 4.8e-15 | 5.9e-16 | 9.7e-16 | 0 | 2.0e-16 | 9.7e-16 | 0.352 | 0.168 | 1.00 / 1.00 |
| strong_peak | np32 (ref) | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1.001 | 0.335 | 2.85 / 1.99 |
| strong_takeover | lat31 | 690 | 6.0e-09 | 1.6e-03 | 1.6e-02 | 1.8e-02 | 9.9e-03 | 2.1e-02 | 4.9e-06 | 1.8e-02 | 0.165 | 0.143 | 0.45 / 0.87 |
| strong_takeover | lat61 | 690 | 0 | 2.6e-11 | 8.9e-12 | 1.5e-11 | 1.8e-11 | 1.5e-11 | 1.5e-15 | 1.5e-11 | 0.194 | 0.148 | 0.53 / 0.90 |
| strong_takeover | lat85 | 690 | 0 | 1.4e-14 | 9.0e-15 | 9.7e-16 | 1.5e-15 | 6.7e-16 | 7.7e-16 | 1.1e-15 | 0.216 | 0.135 | 0.58 / 0.82 |
| strong_takeover | lat101 | 690 | 0 | 2.5e-14 | 6.1e-15 | 8.5e-16 | 6.6e-16 | 6.7e-16 | 7.7e-16 | 8.5e-16 | 0.231 | 0.128 | 0.63 / 0.78 |
| strong_takeover | lat151 | 690 | 0 | 1.3e-14 | 1.5e-14 | 9.7e-16 | 2.1e-15 | 1.3e-15 | 0 | 1.6e-15 | 0.284 | 0.138 | 0.77 / 0.84 |
| strong_takeover | production | 690 | 0 | 1.4e-14 | 2.6e-15 | 9.7e-16 | 3.5e-15 | 6.7e-16 | 0 | 2.5e-15 | 0.369 | 0.165 | 1.00 / 1.00 |
| strong_takeover | np32 (ref) | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1.078 | 0.328 | 2.92 / 1.99 |

Summary over states (largest change over the states listed; cost ratio to production, mean over those states):

| setting | states | amplitude-rate error / fastest rate | db/dt, dj/dt (relative) | cost 1 thr | cost 8 thr |
|---|---|---|---|---|---|
| lat31 | 7 | 1.8e-02 | 2.1e-02 | 0.45 | 0.83 |
| lat61 | 7 | 1.2e-09 | 3.1e-11 | 0.52 | 0.85 |
| lat85 | 7 | 2.2e-13 | 5.3e-15 | 0.58 | 0.81 |
| lat101 | 7 | 8.7e-15 | 7.7e-16 | 0.63 | 0.84 |
| lat151 | 7 | 1.8e-15 | 1.3e-15 | 0.76 | 0.86 |
| production | 7 | 2.5e-15 | 6.7e-16 | 1.00 | 1.00 |
| np32 | 7 | 0 | 0 | 2.87 | 1.98 |

Segment 1.5e9 → 1.9e9 from the plateau state (relative change of the final state from the reference setting):

| setting | rhs calls | ln A₂₁₁ (abs) | ln A₃₂₂ (abs) | b (rel) | j (rel) | wall s (8 thr) |
|---|---|---|---|---|---|---|
| lat31 | 31 | 1.7e-14 | 3.2e-08 | 9.0e-16 | 4.6e-15 | 8.19 |
| lat61 | 31 | 0.0e+00 | 1.4e-13 | 0.0e+00 | 0.0e+00 | 2.78 |
| lat85 | 31 | 0.0e+00 | 5.0e-14 | 0.0e+00 | 0.0e+00 | 2.81 |
| lat101 | 31 | 0.0e+00 | 6.0e-14 | 0.0e+00 | 0.0e+00 | 2.74 |
| lat151 | 31 | 0.0e+00 | 2.8e-14 | 0.0e+00 | 0.0e+00 | 2.93 |
| production | 31 | 0.0e+00 | 1.1e-14 | 0.0e+00 | 0.0e+00 | 3.31 |
| np32 | 31 | 0.0e+00 | 0.0e+00 | 0.0e+00 | 0.0e+00 | 6.33 |

### Study `recommended`

| state | setting | nodes | Re ω (max) | Im ω₂₁₁ | Im ω₃₂₂ | d ln A₂₁₁/dt | d ln A₃₂₂/dt | db/dt | dj/dt | amplitude-rate error / fastest rate | t(1 thr) s | t(8 thr) s | cost vs prod (1/8 thr) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| linear | npol8 lat85 | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.193 | 0.155 | 0.51 / 0.88 |
| linear | production | 690 | 0 | 0 | 0 | 6.7e-16 | 0 | 0 | 0 | 6.7e-16 | 0.379 | 0.177 | 1.00 / 1.00 |
| linear | np32 npol24 (ref) | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2.163 | 0.556 | 5.70 / 3.14 |
| plateau | npol8 lat85 | 690 | 0 | 0 | 0 | 8.5e-16 | 4.4e-16 | 0 | 0 | 4.4e-16 | 0.181 | 0.128 | 0.50 / 0.84 |
| plateau | production | 690 | 0 | 0 | 0 | 8.5e-16 | 4.4e-16 | 0 | 0 | 4.4e-16 | 0.358 | 0.153 | 1.00 / 1.00 |
| plateau | np32 npol24 (ref) | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1.933 | 0.513 | 5.40 / 3.35 |
| takeover | npol8 lat85 | 690 | 0 | 7.2e-13 | 4.4e-13 | 2.6e-15 | 1.2e-15 | 5.1e-16 | 6.0e-16 | 1.2e-15 | 0.191 | 0.127 | 0.50 / 0.73 |
| takeover | production | 690 | 0 | 1.2e-12 | 2.8e-14 | 2.8e-15 | 6.8e-16 | 1.3e-16 | 0 | 6.8e-16 | 0.379 | 0.174 | 1.00 / 1.00 |
| takeover | np32 npol24 (ref) | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2.235 | 0.540 | 5.90 / 3.11 |
| nonlinear | npol8 lat85 | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.195 | 0.120 | 0.51 / 0.70 |
| nonlinear | production | 690 | 0 | 0 | 0 | 0 | 8.2e-16 | 0 | 0 | 6.5e-18 | 0.385 | 0.171 | 1.00 / 1.00 |
| nonlinear | np32 npol24 (ref) | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2.109 | 0.556 | 5.48 / 3.26 |
| large_mesh | npol8 lat85 | 1355 | 0 | 0 | 0 | 0 | 1.7e-16 | 0 | 0 | 7.5e-22 | 0.594 | 0.390 | 0.51 / 0.72 |
| large_mesh | production | 1355 | 0 | 0 | 0 | 0 | 1.7e-16 | 0 | 0 | 7.5e-22 | 1.168 | 0.541 | 1.00 / 1.00 |
| large_mesh | np32 npol24 (ref) | 1355 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 6.526 | 1.556 | 5.59 / 2.87 |
| strong_peak | npol8 lat85 | 690 | 0 | 1.0e-13 | 6.4e-14 | 8.3e-15 | 1.8e-13 | 4.2e-15 | 2.0e-16 | 1.8e-13 | 0.180 | 0.132 | 0.52 / 0.86 |
| strong_peak | production | 690 | 0 | 2.9e-13 | 4.6e-14 | 5.9e-16 | 5.5e-15 | 1.8e-16 | 2.0e-16 | 5.5e-15 | 0.348 | 0.154 | 1.00 / 1.00 |
| strong_peak | np32 npol24 (ref) | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1.957 | 0.480 | 5.63 / 3.11 |
| strong_takeover | npol8 lat85 | 690 | 0 | 6.8e-15 | 6.4e-14 | 1.0e-14 | 5.1e-14 | 1.3e-15 | 7.7e-16 | 3.7e-14 | 0.186 | 0.141 | 0.50 / 0.85 |
| strong_takeover | production | 690 | 0 | 1.6e-14 | 1.5e-16 | 9.7e-16 | 2.1e-15 | 1.3e-15 | 7.7e-16 | 1.6e-15 | 0.369 | 0.166 | 1.00 / 1.00 |
| strong_takeover | np32 npol24 (ref) | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2.114 | 0.539 | 5.72 / 3.25 |

Summary over states (largest change over the states listed; cost ratio to production, mean over those states):

| setting | states | amplitude-rate error / fastest rate | db/dt, dj/dt (relative) | cost 1 thr | cost 8 thr |
|---|---|---|---|---|---|
| npol8 lat85 | 7 | 1.8e-13 | 4.2e-15 | 0.51 | 0.80 |
| production | 7 | 5.5e-15 | 1.3e-15 | 1.00 | 1.00 |
| np32 npol24 | 7 | 0 | 0 | 5.63 | 3.16 |

Segment 1.5e9 → 1.9e9 from the plateau state (relative change of the final state from the reference setting):

| setting | rhs calls | ln A₂₁₁ (abs) | ln A₃₂₂ (abs) | b (rel) | j (rel) | wall s (8 thr) |
|---|---|---|---|---|---|---|
| npol8 lat85 | 31 | 0.0e+00 | 3.9e-14 | 0.0e+00 | 0.0e+00 | 9.03 |
| production | 31 | 0.0e+00 | 1.1e-14 | 0.0e+00 | 0.0e+00 | 3.36 |
| np32 npol24 | 31 | 0.0e+00 | 0.0e+00 | 0.0e+00 | 0.0e+00 | 10.36 |

### Study `n_polar`

| state | setting | nodes | Re ω (max) | Im ω₂₁₁ | Im ω₃₂₂ | d ln A₂₁₁/dt | d ln A₃₂₂/dt | db/dt | dj/dt | amplitude-rate error / fastest rate | t(1 thr) s | t(8 thr) s | cost vs prod (1/8 thr) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| linear | npol8 | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.315 | 0.168 | 0.82 / 1.34 |
| linear | npol10 | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.339 | 0.127 | 0.88 / 1.01 |
| linear | production | 690 | 0 | 0 | 0 | 6.7e-16 | 0 | 0 | 0 | 6.7e-16 | 0.384 | 0.126 | 1.00 / 1.00 |
| linear | npol16 | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.461 | 0.143 | 1.20 / 1.14 |
| linear | npol24 (ref) | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.619 | 0.178 | 1.61 / 1.41 |
| plateau | npol8 | 690 | 0 | 0 | 0 | 8.5e-16 | 4.4e-16 | 0 | 0 | 4.4e-16 | 0.285 | 0.141 | 0.80 / 1.09 |
| plateau | npol10 | 690 | 0 | 0 | 0 | 8.5e-16 | 8.8e-16 | 0 | 0 | 8.8e-16 | 0.297 | 0.128 | 0.83 / 0.99 |
| plateau | production | 690 | 0 | 0 | 0 | 8.5e-16 | 4.4e-16 | 0 | 0 | 4.4e-16 | 0.357 | 0.129 | 1.00 / 1.00 |
| plateau | npol16 | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.419 | 0.153 | 1.17 / 1.19 |
| plateau | npol24 (ref) | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.601 | 0.162 | 1.68 / 1.26 |
| takeover | npol8 | 690 | 0 | 1.7e-13 | 6.8e-14 | 3.5e-15 | 6.8e-16 | 1.4e-15 | 1.8e-15 | 7.6e-16 | 0.284 | 0.134 | 0.77 / 0.99 |
| takeover | npol10 | 690 | 0 | 3.5e-14 | 6.7e-15 | 5.9e-16 | 6.8e-16 | 6.3e-16 | 0 | 6.8e-16 | 0.318 | 0.141 | 0.86 / 1.04 |
| takeover | production | 690 | 0 | 2.7e-13 | 1.4e-14 | 5.9e-16 | 6.8e-16 | 6.3e-16 | 6.0e-16 | 6.8e-16 | 0.369 | 0.136 | 1.00 / 1.00 |
| takeover | npol16 | 690 | 0 | 8.4e-14 | 7.3e-14 | 0 | 0 | 0 | 0 | 0 | 0.463 | 0.174 | 1.25 / 1.28 |
| takeover | npol24 (ref) | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.644 | 0.211 | 1.74 / 1.55 |
| nonlinear | npol8 | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.273 | 0.127 | 0.75 / 0.91 |
| nonlinear | npol10 | 690 | 0 | 0 | 0 | 5.0e-16 | 8.2e-16 | 0 | 0 | 5.0e-16 | 0.321 | 0.129 | 0.88 / 0.92 |
| nonlinear | production | 690 | 0 | 0 | 0 | 0 | 8.2e-16 | 0 | 0 | 6.5e-18 | 0.366 | 0.139 | 1.00 / 1.00 |
| nonlinear | npol16 | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.461 | 0.156 | 1.26 / 1.12 |
| nonlinear | npol24 (ref) | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.649 | 0.231 | 1.78 / 1.65 |
| large_mesh | npol8 | 1355 | 0 | 0 | 0 | 0 | 3.5e-16 | 0 | 0 | 1.5e-21 | 0.969 | 0.454 | 0.76 / 0.92 |
| large_mesh | npol10 | 1355 | 0 | 0 | 0 | 0 | 3.5e-16 | 0 | 0 | 1.5e-21 | 1.101 | 0.485 | 0.87 / 0.98 |
| large_mesh | production | 1355 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1.272 | 0.494 | 1.00 / 1.00 |
| large_mesh | npol16 | 1355 | 0 | 0 | 0 | 0 | 3.5e-16 | 0 | 0 | 1.5e-21 | 1.458 | 0.562 | 1.15 / 1.14 |
| large_mesh | npol24 (ref) | 1355 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2.059 | 0.638 | 1.62 / 1.29 |
| strong_peak | npol8 | 690 | 0 | 4.2e-13 | 1.9e-13 | 1.2e-15 | 1.7e-14 | 3.5e-16 | 2.0e-16 | 1.7e-14 | 0.290 | 0.120 | 0.84 / 0.92 |
| strong_peak | npol10 | 690 | 0 | 4.0e-13 | 1.9e-13 | 5.9e-16 | 4.4e-15 | 3.5e-16 | 0 | 4.4e-15 | 0.312 | 0.130 | 0.90 / 0.99 |
| strong_peak | production | 690 | 0 | 2.3e-13 | 6.3e-14 | 5.9e-16 | 1.7e-15 | 1.8e-16 | 0 | 1.7e-15 | 0.347 | 0.131 | 1.00 / 1.00 |
| strong_peak | npol16 | 690 | 0 | 9.0e-14 | 1.3e-13 | 5.9e-16 | 1.2e-15 | 1.8e-16 | 0 | 1.2e-15 | 0.432 | 0.166 | 1.25 / 1.27 |
| strong_peak | npol24 (ref) | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.621 | 0.189 | 1.79 / 1.44 |
| strong_takeover | npol8 | 690 | 0 | 8.3e-15 | 2.9e-14 | 8.2e-15 | 4.8e-14 | 0 | 7.7e-16 | 3.6e-14 | 0.282 | 0.123 | 0.77 / 0.83 |
| strong_takeover | npol10 | 690 | 0 | 1.2e-14 | 5.0e-15 | 2.8e-15 | 1.6e-16 | 6.7e-16 | 7.7e-16 | 2.8e-15 | 0.330 | 0.137 | 0.90 / 0.93 |
| strong_takeover | production | 690 | 0 | 8.0e-15 | 2.0e-15 | 0 | 8.2e-16 | 6.7e-16 | 7.7e-16 | 6.0e-16 | 0.365 | 0.148 | 1.00 / 1.00 |
| strong_takeover | npol16 | 690 | 0 | 5.8e-15 | 8.7e-15 | 9.7e-16 | 9.9e-16 | 6.7e-16 | 7.7e-16 | 9.7e-16 | 0.460 | 0.173 | 1.26 / 1.17 |
| strong_takeover | npol24 (ref) | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.635 | 0.211 | 1.74 / 1.43 |

Summary over states (largest change over the states listed; cost ratio to production, mean over those states):

| setting | states | amplitude-rate error / fastest rate | db/dt, dj/dt (relative) | cost 1 thr | cost 8 thr |
|---|---|---|---|---|---|
| npol8 | 7 | 3.6e-14 | 1.8e-15 | 0.79 | 1.00 |
| npol10 | 7 | 4.4e-15 | 7.7e-16 | 0.87 | 0.98 |
| production | 7 | 1.7e-15 | 7.7e-16 | 1.00 | 1.00 |
| npol16 | 7 | 1.2e-15 | 7.7e-16 | 1.22 | 1.18 |
| npol24 | 7 | 0 | 0 | 1.71 | 1.43 |

Segment 1.5e9 → 1.9e9 from the plateau state (relative change of the final state from the reference setting):

| setting | rhs calls | ln A₂₁₁ (abs) | ln A₃₂₂ (abs) | b (rel) | j (rel) | wall s (8 thr) |
|---|---|---|---|---|---|---|
| npol8 | 31 | 0.0e+00 | 1.1e-14 | 0.0e+00 | 0.0e+00 | 9.25 |
| npol10 | 31 | 0.0e+00 | 0.0e+00 | 0.0e+00 | 0.0e+00 | 2.73 |
| production | 31 | 0.0e+00 | 1.1e-14 | 0.0e+00 | 0.0e+00 | 2.87 |
| npol16 | 31 | 0.0e+00 | 0.0e+00 | 0.0e+00 | 0.0e+00 | 3.32 |
| npol24 | 31 | 0.0e+00 | 0.0e+00 | 0.0e+00 | 0.0e+00 | 3.96 |

### Study `step`

| state | setting | nodes | Re ω (max) | Im ω₂₁₁ | Im ω₃₂₂ | d ln A₂₁₁/dt | d ln A₃₂₂/dt | db/dt | dj/dt | amplitude-rate error / fastest rate | t(1 thr) s | t(8 thr) s | cost vs prod (1/8 thr) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| linear | h0.8 | 444 | 2.0e-10 | 2.2e-05 | 2.9e-04 | 1.2e-06 | 1.3e-05 | 2.2e-06 | 2.2e-06 | 1.2e-06 | 0.256 | 0.112 | 0.67 / 0.82 |
| linear | h0.6 | 526 | 6.4e-11 | 7.1e-06 | 9.1e-05 | 3.8e-07 | 4.0e-06 | 1.3e-06 | 1.3e-06 | 3.8e-07 | 0.288 | 0.143 | 0.75 / 1.05 |
| linear | production | 690 | 1.2e-11 | 1.3e-06 | 1.7e-05 | 7.0e-08 | 7.5e-07 | 9.0e-07 | 8.9e-07 | 7.0e-08 | 0.382 | 0.137 | 1.00 / 1.00 |
| linear | h0.3 | 851 | 3.2e-12 | 3.6e-07 | 4.6e-06 | 1.9e-08 | 2.0e-07 | 8.3e-07 | 8.3e-07 | 1.9e-08 | 0.463 | 0.196 | 1.21 / 1.43 |
| linear | h0.2 (ref) | 1197 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.652 | 0.282 | 1.71 / 2.06 |
| plateau | h0.8 | 444 | 1.7e-10 | 3.5e-02 | 4.1e-07 | 3.5e-02 | 1.6e-05 | 3.4e-02 | 3.3e-02 | 1.6e-05 | 0.231 | 0.089 | 0.66 / 0.60 |
| plateau | h0.6 | 526 | 4.9e-11 | 1.0e-02 | 7.8e-07 | 1.0e-02 | 6.0e-06 | 9.9e-03 | 9.7e-03 | 6.0e-06 | 0.272 | 0.135 | 0.78 / 0.92 |
| plateau | production | 690 | 3.7e-12 | 7.9e-04 | 1.2e-06 | 7.9e-04 | 2.2e-06 | 7.6e-04 | 7.4e-04 | 2.2e-06 | 0.349 | 0.147 | 1.00 / 1.00 |
| plateau | h0.3 | 851 | 3.7e-12 | 7.7e-04 | 1.3e-06 | 7.8e-04 | 1.5e-06 | 7.5e-04 | 7.3e-04 | 1.5e-06 | 0.439 | 0.197 | 1.26 / 1.34 |
| plateau | h0.2 (ref) | 1197 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.618 | 0.245 | 1.77 / 1.66 |
| takeover | h0.8 | 444 | 1.7e-10 | 8.5e-06 | 1.2e-06 | 1.1e-05 | 1.7e-05 | 9.4e-06 | 5.4e-06 | 1.7e-05 | 0.252 | 0.110 | 0.66 / 0.66 |
| takeover | h0.6 | 526 | 5.1e-11 | 2.5e-06 | 5.8e-07 | 3.7e-06 | 6.3e-06 | 3.6e-06 | 2.1e-06 | 6.3e-06 | 0.291 | 0.119 | 0.77 / 0.71 |
| takeover | production | 690 | 5.0e-12 | 1.8e-07 | 1.2e-06 | 9.8e-07 | 2.3e-06 | 1.4e-06 | 8.7e-07 | 2.3e-06 | 0.380 | 0.168 | 1.00 / 1.00 |
| takeover | h0.3 | 851 | 2.6e-12 | 1.9e-07 | 1.3e-06 | 5.3e-07 | 1.6e-06 | 1.0e-06 | 6.6e-07 | 1.6e-06 | 0.462 | 0.203 | 1.22 / 1.21 |
| takeover | h0.2 (ref) | 1197 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.643 | 0.280 | 1.69 / 1.67 |
| nonlinear | h0.8 | 444 | 1.9e-10 | 2.1e-05 | 3.1e-04 | 1.7e-06 | 2.0e-05 | 4.9e-05 | 9.2e-05 | 1.7e-06 | 0.255 | 0.109 | 0.68 / 0.66 |
| nonlinear | h0.6 | 526 | 6.2e-11 | 6.8e-06 | 9.7e-05 | 5.5e-07 | 6.4e-06 | 1.7e-05 | 3.1e-05 | 5.5e-07 | 0.291 | 0.122 | 0.77 / 0.74 |
| nonlinear | production | 690 | 1.1e-11 | 1.3e-06 | 1.8e-05 | 1.0e-07 | 1.2e-06 | 5.0e-06 | 7.5e-06 | 1.0e-07 | 0.377 | 0.166 | 1.00 / 1.00 |
| nonlinear | h0.3 | 851 | 3.0e-12 | 3.4e-07 | 4.9e-06 | 2.8e-08 | 3.2e-07 | 3.0e-06 | 3.7e-06 | 2.8e-08 | 0.465 | 0.180 | 1.23 / 1.09 |
| nonlinear | h0.2 (ref) | 1197 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.642 | 0.305 | 1.70 / 1.84 |
| large_mesh | h0.8 | 841 | 5.3e-12 | 3.4e-05 | 3.1e-04 | 1.8e-06 | 1.7e-05 | 3.1e-05 | 3.4e-05 | 1.8e-06 | 0.608 | 0.256 | 0.51 / 0.49 |
| large_mesh | h0.6 | 1013 | 1.7e-12 | 1.0e-05 | 9.6e-05 | 5.7e-07 | 5.4e-06 | 9.8e-06 | 1.1e-05 | 5.7e-07 | 0.740 | 0.335 | 0.62 / 0.64 |
| large_mesh | production | 1355 | 3.0e-13 | 2.0e-06 | 1.8e-05 | 1.1e-07 | 1.1e-06 | 1.9e-06 | 2.0e-06 | 1.1e-07 | 1.190 | 0.520 | 1.00 / 1.00 |
| large_mesh | h0.3 | 1765 | 7.9e-14 | 5.3e-07 | 4.9e-06 | 2.9e-08 | 3.1e-07 | 5.1e-07 | 5.6e-07 | 2.9e-08 | 1.655 | 0.630 | 1.39 / 1.21 |
| large_mesh | h0.2 (ref) | 2626 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2.600 | 0.972 | 2.18 / 1.87 |
| strong_peak | h0.8 | 444 | 1.1e-10 | 2.4e-05 | 1.1e-06 | 5.8e-07 | 1.8e-05 | 1.7e-06 | 2.1e-06 | 1.8e-05 | 0.238 | 0.097 | 0.68 / 0.62 |
| strong_peak | h0.6 | 526 | 2.3e-10 | 7.4e-06 | 1.4e-06 | 2.7e-07 | 6.8e-06 | 1.0e-06 | 1.2e-06 | 6.8e-06 | 0.270 | 0.106 | 0.77 / 0.68 |
| strong_peak | production | 690 | 2.7e-10 | 1.4e-06 | 1.5e-06 | 1.6e-07 | 2.5e-06 | 7.7e-07 | 8.1e-07 | 2.5e-06 | 0.349 | 0.155 | 1.00 / 1.00 |
| strong_peak | h0.3 | 851 | 2.7e-10 | 3.5e-07 | 1.5e-06 | 1.3e-07 | 1.7e-06 | 7.2e-07 | 7.4e-07 | 1.7e-06 | 0.427 | 0.168 | 1.22 / 1.08 |
| strong_peak | h0.2 (ref) | 1197 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.590 | 0.228 | 1.69 / 1.47 |
| strong_takeover | h0.8 | 444 | 4.9e-11 | 3.1e-06 | 1.3e-06 | 2.1e-05 | 2.7e-05 | 2.2e-05 | 2.1e-06 | 2.1e-05 | 0.245 | 0.096 | 0.60 / 0.62 |
| strong_takeover | h0.6 | 526 | 1.2e-10 | 1.6e-06 | 1.7e-06 | 7.3e-06 | 1.0e-05 | 8.2e-06 | 1.2e-06 | 7.3e-06 | 0.296 | 0.122 | 0.72 / 0.79 |
| strong_takeover | production | 690 | 1.7e-10 | 1.0e-06 | 1.9e-06 | 2.1e-06 | 3.5e-06 | 2.9e-06 | 8.4e-07 | 2.5e-06 | 0.411 | 0.154 | 1.00 / 1.00 |
| strong_takeover | h0.3 | 851 | 1.7e-10 | 9.2e-07 | 1.9e-06 | 1.3e-06 | 2.4e-06 | 1.9e-06 | 7.7e-07 | 1.7e-06 | 0.457 | 0.179 | 1.11 / 1.16 |
| strong_takeover | h0.2 (ref) | 1197 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.628 | 0.265 | 1.53 / 1.71 |

Summary over states (largest change over the states listed; cost ratio to production, mean over those states):

| setting | states | amplitude-rate error / fastest rate | db/dt, dj/dt (relative) | cost 1 thr | cost 8 thr |
|---|---|---|---|---|---|
| h0.8 | 7 | 2.1e-05 | 3.4e-02 | 0.64 | 0.64 |
| h0.6 | 7 | 7.3e-06 | 9.9e-03 | 0.74 | 0.79 |
| production | 7 | 2.5e-06 | 7.6e-04 | 1.00 | 1.00 |
| h0.3 | 7 | 1.7e-06 | 7.5e-04 | 1.23 | 1.22 |
| h0.2 | 7 | 0 | 0 | 1.75 | 1.76 |

Segment 1.5e9 → 1.9e9 from the plateau state (relative change of the final state from the reference setting):

| setting | rhs calls | ln A₂₁₁ (abs) | ln A₃₂₂ (abs) | b (rel) | j (rel) | wall s (8 thr) |
|---|---|---|---|---|---|---|
| h0.8 | 31 | 4.1e-09 | 1.5e-06 | 1.8e-10 | 7.3e-10 | 8.31 |
| h0.6 | 31 | 1.2e-09 | 5.7e-07 | 5.5e-11 | 2.2e-10 | 2.65 |
| production | 31 | 9.3e-11 | 2.1e-07 | 4.3e-12 | 1.7e-11 | 3.31 |
| h0.3 | 31 | 9.4e-11 | 1.4e-07 | 4.2e-12 | 1.7e-11 | 4.01 |
| h0.2 | 31 | 0.0e+00 | 0.0e+00 | 0.0e+00 | 0.0e+00 | 5.45 |

### Study `outer`

| state | setting | nodes | Re ω (max) | Im ω₂₁₁ | Im ω₃₂₂ | d ln A₂₁₁/dt | d ln A₃₂₂/dt | db/dt | dj/dt | amplitude-rate error / fastest rate | t(1 thr) s | t(8 thr) s | cost vs prod (1/8 thr) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| linear | c0.001 | 642 | 0 | 7.8e-16 | 3.1e-14 | 9.5e-15 | 2.2e-09 | 2.8e-06 | 2.8e-06 | 2.1e-13 | 0.339 | 0.166 | 0.92 / 1.01 |
| linear | production | 690 | 0 | 7.8e-16 | 2.8e-15 | 6.0e-15 | 4.6e-10 | 2.2e-06 | 2.2e-06 | 4.5e-14 | 0.367 | 0.164 | 1.00 / 1.00 |
| linear | c1e-05 | 737 | 0 | 0 | 1.4e-15 | 6.0e-15 | 3.8e-10 | 1.7e-06 | 1.7e-06 | 3.7e-14 | 0.397 | 0.171 | 1.08 / 1.05 |
| linear | c1e-06 | 783 | 0 | 7.8e-16 | 1.4e-15 | 9.5e-15 | 2.6e-10 | 1.1e-06 | 1.1e-06 | 2.5e-14 | 0.408 | 0.196 | 1.11 / 1.20 |
| linear | c1e-08 (ref) | 873 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.456 | 0.221 | 1.24 / 1.35 |
| plateau | c0.001 | 642 | 2.3e-11 | 4.8e-03 | 4.5e-06 | 4.8e-03 | 4.5e-06 | 4.6e-03 | 4.5e-03 | 4.5e-06 | 0.321 | 0.120 | 0.94 / 0.72 |
| plateau | production | 690 | 1.8e-11 | 3.8e-03 | 3.6e-06 | 3.8e-03 | 3.6e-06 | 3.7e-03 | 3.6e-03 | 3.6e-06 | 0.340 | 0.166 | 1.00 / 1.00 |
| plateau | c1e-05 | 737 | 1.4e-11 | 2.9e-03 | 2.7e-06 | 2.9e-03 | 2.7e-06 | 2.8e-03 | 2.7e-03 | 2.7e-06 | 0.369 | 0.185 | 1.08 / 1.11 |
| plateau | c1e-06 | 783 | 9.2e-12 | 1.9e-03 | 1.8e-06 | 1.9e-03 | 1.8e-06 | 1.9e-03 | 1.8e-03 | 1.8e-06 | 0.407 | 0.195 | 1.20 / 1.17 |
| plateau | c1e-08 (ref) | 873 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.433 | 0.207 | 1.27 / 1.24 |
| takeover | c0.001 | 642 | 1.9e-11 | 1.2e-06 | 4.7e-06 | 1.3e-06 | 4.0e-06 | 3.1e-06 | 2.0e-06 | 4.0e-06 | 0.340 | 0.162 | 0.92 / 0.98 |
| takeover | production | 690 | 1.5e-11 | 9.7e-07 | 3.8e-06 | 1.1e-06 | 3.7e-06 | 2.5e-06 | 1.6e-06 | 3.7e-06 | 0.369 | 0.164 | 1.00 / 1.00 |
| takeover | c1e-05 | 737 | 1.2e-11 | 7.3e-07 | 2.9e-06 | 7.9e-07 | 2.8e-06 | 1.9e-06 | 1.2e-06 | 2.8e-06 | 0.377 | 0.180 | 1.02 / 1.10 |
| takeover | c1e-06 | 783 | 7.7e-12 | 4.9e-07 | 1.9e-06 | 5.3e-07 | 1.9e-06 | 1.2e-06 | 8.2e-07 | 1.9e-06 | 0.440 | 0.185 | 1.19 / 1.12 |
| takeover | c1e-08 (ref) | 873 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.459 | 0.215 | 1.24 / 1.31 |
| nonlinear | c0.001 | 642 | 1.2e-13 | 4.9e-09 | 2.0e-07 | 7.7e-08 | 5.7e-07 | 7.6e-06 | 7.4e-06 | 7.7e-08 | 0.373 | 0.169 | 1.01 / 1.03 |
| nonlinear | production | 690 | 1.3e-13 | 4.1e-09 | 1.6e-07 | 6.3e-08 | 2.0e-07 | 6.1e-06 | 5.9e-06 | 6.3e-08 | 0.370 | 0.164 | 1.00 / 1.00 |
| nonlinear | c1e-05 | 737 | 1.2e-13 | 3.2e-09 | 1.2e-07 | 4.7e-08 | 1.5e-07 | 4.6e-06 | 4.5e-06 | 4.7e-08 | 0.402 | 0.170 | 1.09 / 1.04 |
| nonlinear | c1e-06 | 783 | 7.0e-14 | 2.3e-09 | 8.1e-08 | 3.2e-08 | 9.8e-08 | 3.1e-06 | 3.0e-06 | 3.2e-08 | 0.419 | 0.178 | 1.13 / 1.08 |
| nonlinear | c1e-08 (ref) | 873 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.464 | 0.209 | 1.25 / 1.28 |
| large_mesh | c0.001 | 1277 | 6.4e-14 | 6.3e-09 | 3.4e-08 | 1.1e-05 | 4.7e-05 | 1.2e-07 | 8.0e-08 | 1.1e-05 | 0.950 | 0.395 | 0.70 / 0.78 |
| large_mesh | production | 1355 | 5.3e-14 | 5.2e-09 | 2.1e-08 | 8.6e-06 | 8.0e-06 | 7.9e-08 | 5.4e-08 | 8.6e-06 | 1.351 | 0.507 | 1.00 / 1.00 |
| large_mesh | c1e-05 | 1430 | 4.1e-14 | 4.0e-09 | 1.4e-08 | 6.5e-06 | 5.6e-06 | 5.4e-08 | 3.8e-08 | 6.5e-06 | 1.437 | 0.503 | 1.06 / 0.99 |
| large_mesh | c1e-06 | 1504 | 2.8e-14 | 2.8e-09 | 1.1e-08 | 4.3e-06 | 3.7e-06 | 4.1e-08 | 2.8e-08 | 4.3e-06 | 1.125 | 0.419 | 0.83 / 0.83 |
| large_mesh | c1e-08 (ref) | 1650 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1.640 | 0.596 | 1.21 / 1.17 |
| strong_peak | c0.001 | 642 | 9.5e-10 | 7.1e-08 | 5.0e-06 | 4.5e-07 | 5.1e-06 | 2.4e-06 | 2.5e-06 | 5.1e-06 | 0.355 | 0.149 | 0.96 / 1.00 |
| strong_peak | production | 690 | 7.6e-10 | 5.7e-08 | 4.0e-06 | 3.7e-07 | 4.2e-06 | 1.9e-06 | 2.0e-06 | 4.2e-06 | 0.371 | 0.149 | 1.00 / 1.00 |
| strong_peak | c1e-05 | 737 | 5.7e-10 | 4.2e-08 | 3.0e-06 | 2.8e-07 | 3.1e-06 | 1.5e-06 | 1.5e-06 | 3.1e-06 | 0.381 | 0.157 | 1.03 / 1.06 |
| strong_peak | c1e-06 | 783 | 3.9e-10 | 2.8e-08 | 2.0e-06 | 1.8e-07 | 2.1e-06 | 9.8e-07 | 1.0e-06 | 2.1e-06 | 0.397 | 0.172 | 1.07 / 1.15 |
| strong_peak | c1e-08 (ref) | 873 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.430 | 0.192 | 1.16 / 1.29 |
| strong_takeover | c0.001 | 642 | 1.5e-09 | 7.5e-06 | 1.6e-05 | 1.6e-04 | 1.1e-03 | 1.4e-05 | 6.3e-06 | 7.7e-04 | 0.359 | 0.142 | 0.94 / 0.89 |
| strong_takeover | production | 690 | 2.1e-10 | 1.3e-06 | 2.1e-06 | 1.4e-04 | 8.6e-04 | 2.2e-06 | 8.8e-07 | 6.3e-04 | 0.383 | 0.159 | 1.00 / 1.00 |
| strong_takeover | c1e-05 | 737 | 3.6e-10 | 1.9e-06 | 3.8e-06 | 1.1e-04 | 6.4e-04 | 3.5e-06 | 1.6e-06 | 4.7e-04 | 0.418 | 0.167 | 1.09 / 1.05 |
| strong_takeover | c1e-06 | 783 | 3.0e-10 | 1.6e-06 | 3.2e-06 | 7.0e-05 | 4.3e-04 | 2.9e-06 | 1.3e-06 | 3.1e-04 | 0.466 | 0.172 | 1.22 / 1.08 |
| strong_takeover | c1e-08 (ref) | 873 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.497 | 0.210 | 1.30 / 1.32 |

Summary over states (largest change over the states listed; cost ratio to production, mean over those states):

| setting | states | amplitude-rate error / fastest rate | db/dt, dj/dt (relative) | cost 1 thr | cost 8 thr |
|---|---|---|---|---|---|
| c0.001 | 7 | 7.7e-04 | 4.6e-03 | 0.91 | 0.92 |
| production | 7 | 6.3e-04 | 3.7e-03 | 1.00 | 1.00 |
| c1e-05 | 7 | 4.7e-04 | 2.8e-03 | 1.06 | 1.06 |
| c1e-06 | 7 | 3.1e-04 | 1.9e-03 | 1.11 | 1.09 |
| c1e-08 | 7 | 0 | 0 | 1.24 | 1.28 |

Segment 1.5e9 → 1.9e9 from the plateau state (relative change of the final state from the reference setting):

| setting | rhs calls | ln A₂₁₁ (abs) | ln A₃₂₂ (abs) | b (rel) | j (rel) | wall s (8 thr) |
|---|---|---|---|---|---|---|
| c0.001 | 31 | 5.7e-10 | 4.3e-07 | 2.6e-11 | 1.0e-10 | 8.78 |
| production | 31 | 4.6e-10 | 3.4e-07 | 2.1e-11 | 8.3e-11 | 3.46 |
| c1e-05 | 31 | 3.5e-10 | 2.6e-07 | 1.6e-11 | 6.2e-11 | 3.69 |
| c1e-06 | 31 | 2.3e-10 | 1.7e-07 | 1.1e-11 | 4.2e-11 | 3.87 |
| c1e-08 | 31 | 0.0e+00 | 0.0e+00 | 0.0e+00 | 0.0e+00 | 4.21 |

### Study `outer_fine_tail`

| state | setting | nodes | Re ω (max) | Im ω₂₁₁ | Im ω₃₂₂ | d ln A₂₁₁/dt | d ln A₃₂₂/dt | db/dt | dj/dt | amplitude-rate error / fastest rate | t(1 thr) s | t(8 thr) s | cost vs prod (1/8 thr) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| plateau | ms0.8 | 941 | 6.2e-13 | 1.3e-04 | 1.2e-07 | 1.3e-04 | 1.2e-07 | 1.3e-04 | 1.2e-04 | 1.2e-07 | NaN | 0.216 | NaN / NaN |
| plateau | c1e-06 ms0.8 | 1158 | 3.1e-13 | 6.5e-05 | 6.2e-08 | 6.5e-05 | 6.2e-08 | 6.3e-05 | 6.2e-05 | 6.2e-08 | NaN | 0.240 | NaN / NaN |
| plateau | c1e-08 ms0.8 (ref) | 1368 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | NaN | 0.290 | NaN / NaN |
| takeover | ms0.8 | 941 | 5.3e-13 | 3.3e-08 | 1.3e-07 | 3.2e-08 | 8.7e-08 | 8.4e-08 | 5.5e-08 | 8.7e-08 | NaN | 0.224 | NaN / NaN |
| takeover | c1e-06 ms0.8 | 1158 | 2.7e-13 | 1.7e-08 | 6.5e-08 | 1.6e-08 | 5.1e-08 | 4.2e-08 | 2.8e-08 | 5.1e-08 | NaN | 0.274 | NaN / NaN |
| takeover | c1e-08 ms0.8 (ref) | 1368 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | NaN | 0.300 | NaN / NaN |
| large_mesh | ms0.8 | 2209 | 5.9e-16 | 4.2e-12 | 7.2e-10 | 8.6e-06 | 8.2e-06 | 2.3e-09 | 1.3e-09 | 8.6e-06 | NaN | 0.616 | NaN / NaN |
| large_mesh | c1e-06 ms0.8 | 2766 | 4.5e-16 | 1.5e-11 | 5.5e-10 | 4.3e-06 | 3.8e-06 | 1.8e-09 | 1.0e-09 | 4.3e-06 | NaN | 1.060 | NaN / NaN |
| large_mesh | c1e-08 ms0.8 (ref) | 3308 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | NaN | 1.223 | NaN / NaN |
| strong_peak | ms0.8 | 941 | 2.6e-11 | 1.9e-09 | 1.4e-07 | 1.2e-08 | 2.9e-07 | 6.6e-08 | 6.7e-08 | 2.9e-07 | NaN | 0.223 | NaN / NaN |
| strong_peak | c1e-06 ms0.8 | 1158 | 1.3e-11 | 9.6e-10 | 6.9e-08 | 5.3e-09 | 1.4e-07 | 3.3e-08 | 3.4e-08 | 1.4e-07 | NaN | 0.254 | NaN / NaN |
| strong_peak | c1e-08 ms0.8 (ref) | 1368 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | NaN | 0.293 | NaN / NaN |
| strong_takeover | ms0.8 | 941 | 2.6e-10 | 1.2e-06 | 2.9e-06 | 1.5e-04 | 8.6e-04 | 2.2e-06 | 1.1e-06 | 6.3e-04 | NaN | 0.204 | NaN / NaN |
| strong_takeover | c1e-06 ms0.8 | 1158 | 6.8e-11 | 3.6e-07 | 7.1e-07 | 7.1e-05 | 4.3e-04 | 6.4e-07 | 2.9e-07 | 3.2e-04 | NaN | 0.264 | NaN / NaN |
| strong_takeover | c1e-08 ms0.8 (ref) | 1368 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | NaN | 0.310 | NaN / NaN |

Summary over states (largest change over the states listed; cost ratio to production, mean over those states):

| setting | states | amplitude-rate error / fastest rate | db/dt, dj/dt (relative) | cost 1 thr | cost 8 thr |
|---|---|---|---|---|---|
| ms0.8 | 5 | 6.3e-04 | 1.3e-04 | NaN | NaN |
| c1e-06 ms0.8 | 5 | 3.2e-04 | 6.3e-05 | NaN | NaN |
| c1e-08 ms0.8 | 5 | 0 | 0 | NaN | NaN |

Segment 1.5e9 → 1.9e9 from the plateau state (relative change of the final state from the reference setting):

| setting | rhs calls | ln A₂₁₁ (abs) | ln A₃₂₂ (abs) | b (rel) | j (rel) | wall s (8 thr) |
|---|---|---|---|---|---|---|
| ms0.8 | 31 | 1.6e-11 | 1.2e-08 | 7.0e-13 | 2.8e-12 | 9.78 |
| c1e-06 ms0.8 | 31 | 7.8e-12 | 5.9e-09 | 3.5e-13 | 1.4e-12 | 5.26 |
| c1e-08 ms0.8 | 31 | 0.0e+00 | 0.0e+00 | 0.0e+00 | 0.0e+00 | 5.88 |

### Study `inner`

| state | setting | nodes | Re ω (max) | Im ω₂₁₁ | Im ω₃₂₂ | d ln A₂₁₁/dt | d ln A₃₂₂/dt | db/dt | dj/dt | amplitude-rate error / fastest rate | t(1 thr) s | t(8 thr) s | cost vs prod (1/8 thr) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| linear | hr0.0001 | 652 | 9.4e-14 | 1.0e-08 | 1.3e-07 | 7.5e-10 | 3.9e-09 | 9.0e-07 | 9.0e-07 | 7.5e-10 | 0.387 | 0.179 | 1.04 / 1.18 |
| linear | hr1e-05 | 671 | 9.6e-14 | 1.1e-08 | 1.4e-07 | 6.0e-10 | 6.1e-09 | 6.7e-07 | 6.7e-07 | 6.0e-10 | 0.367 | 0.150 | 0.98 / 0.99 |
| linear | production | 690 | 9.9e-14 | 1.1e-08 | 1.4e-07 | 6.6e-10 | 6.1e-09 | 4.5e-07 | 4.5e-07 | 6.6e-10 | 0.374 | 0.151 | 1.00 / 1.00 |
| linear | hr1e-08 (ref) | 727 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.390 | 0.151 | 1.04 / 1.00 |
| plateau | hr0.0001 | 652 | 8.1e-14 | 1.7e-05 | 7.9e-10 | 1.6e-05 | 2.1e-08 | 1.7e-05 | 1.6e-05 | 2.1e-08 | 0.346 | 0.132 | 0.95 / 0.87 |
| plateau | hr1e-05 | 671 | 8.4e-14 | 1.8e-05 | 8.3e-10 | 1.7e-05 | 1.4e-08 | 1.7e-05 | 1.7e-05 | 1.4e-08 | 0.345 | 0.155 | 0.94 / 1.02 |
| plateau | production | 690 | 8.6e-14 | 1.8e-05 | 8.4e-10 | 1.7e-05 | 6.7e-09 | 1.7e-05 | 1.7e-05 | 6.7e-09 | 0.366 | 0.151 | 1.00 / 1.00 |
| plateau | hr1e-08 (ref) | 727 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.370 | 0.156 | 1.01 / 1.03 |
| takeover | hr0.0001 | 652 | 9.0e-14 | 4.0e-09 | 1.2e-09 | 8.6e-09 | 1.1e-08 | 2.4e-09 | 2.0e-10 | 1.1e-08 | 0.380 | 0.150 | 1.02 / 0.97 |
| takeover | hr1e-05 | 671 | 8.6e-14 | 4.3e-09 | 1.3e-09 | 5.2e-09 | 6.2e-09 | 6.3e-10 | 5.7e-10 | 6.2e-09 | 0.367 | 0.146 | 0.98 / 0.95 |
| takeover | production | 690 | 8.7e-14 | 4.3e-09 | 1.3e-09 | 1.6e-09 | 1.5e-09 | 1.0e-09 | 1.2e-09 | 1.5e-09 | 0.374 | 0.154 | 1.00 / 1.00 |
| takeover | hr1e-08 (ref) | 727 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.402 | 0.158 | 1.07 / 1.03 |
| nonlinear | hr0.0001 | 652 | 9.2e-14 | 1.0e-08 | 1.4e-07 | 2.9e-09 | 9.5e-09 | 7.9e-09 | 1.0e-08 | 2.9e-09 | 0.378 | 0.138 | 1.00 / 0.86 |
| nonlinear | hr1e-05 | 671 | 9.4e-14 | 1.0e-08 | 1.5e-07 | 2.1e-09 | 1.0e-08 | 8.6e-10 | 2.1e-08 | 2.1e-09 | 0.370 | 0.154 | 0.98 / 0.96 |
| nonlinear | production | 690 | 9.7e-14 | 1.0e-08 | 1.5e-07 | 1.0e-09 | 1.0e-08 | 8.5e-09 | 2.9e-08 | 1.0e-09 | 0.378 | 0.161 | 1.00 / 1.00 |
| nonlinear | hr1e-08 (ref) | 727 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.396 | 0.148 | 1.05 / 0.92 |
| large_mesh | hr0.0001 | 1317 | 0 | 1.1e-09 | 1.2e-08 | 1.1e-10 | 3.3e-08 | 1.5e-08 | 9.7e-09 | 1.1e-10 | 1.417 | 0.662 | 1.13 / 1.26 |
| large_mesh | hr1e-05 | 1336 | 0 | 9.0e-10 | 8.2e-09 | 9.5e-12 | 2.4e-08 | 9.2e-09 | 4.9e-09 | 9.5e-12 | 1.434 | 0.638 | 1.15 / 1.22 |
| large_mesh | production | 1355 | 0 | 5.9e-10 | 5.3e-09 | 3.1e-11 | 1.6e-08 | 5.9e-09 | 3.0e-09 | 3.1e-11 | 1.250 | 0.524 | 1.00 / 1.00 |
| large_mesh | hr1e-08 (ref) | 1393 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1.558 | 0.603 | 1.25 / 1.15 |
| strong_peak | hr0.0001 | 652 | 8.9e-14 | 1.5e-08 | 1.9e-10 | 2.4e-06 | 4.6e-05 | 1.4e-06 | 7.9e-07 | 4.6e-05 | 0.376 | 0.146 | 1.09 / 0.97 |
| strong_peak | hr1e-05 | 671 | 9.1e-14 | 1.2e-08 | 2.3e-10 | 1.8e-06 | 3.4e-05 | 1.0e-06 | 6.0e-07 | 3.4e-05 | 0.339 | 0.137 | 0.98 / 0.91 |
| strong_peak | production | 690 | 9.3e-14 | 1.2e-08 | 2.1e-10 | 1.2e-06 | 2.3e-05 | 6.8e-07 | 4.0e-07 | 2.3e-05 | 0.345 | 0.151 | 1.00 / 1.00 |
| strong_peak | hr1e-08 (ref) | 727 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.364 | 0.135 | 1.05 / 0.89 |
| strong_takeover | hr0.0001 | 652 | 9.0e-14 | 1.4e-09 | 3.1e-10 | 1.3e-05 | 1.6e-05 | 2.0e-05 | 5.8e-06 | 1.3e-05 | 0.353 | 0.141 | 0.95 / 0.95 |
| strong_takeover | hr1e-05 | 671 | 9.2e-14 | 1.1e-09 | 3.5e-10 | 9.5e-06 | 1.2e-05 | 1.5e-05 | 4.4e-06 | 9.5e-06 | 0.361 | 0.143 | 0.97 / 0.96 |
| strong_takeover | production | 690 | 9.4e-14 | 1.0e-09 | 3.4e-10 | 6.3e-06 | 7.9e-06 | 1.0e-05 | 2.9e-06 | 6.3e-06 | 0.371 | 0.149 | 1.00 / 1.00 |
| strong_takeover | hr1e-08 (ref) | 727 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.383 | 0.140 | 1.03 / 0.94 |

Summary over states (largest change over the states listed; cost ratio to production, mean over those states):

| setting | states | amplitude-rate error / fastest rate | db/dt, dj/dt (relative) | cost 1 thr | cost 8 thr |
|---|---|---|---|---|---|
| hr0.0001 | 7 | 4.6e-05 | 2.0e-05 | 1.03 | 1.01 |
| hr1e-05 | 7 | 3.4e-05 | 1.7e-05 | 1.00 | 1.00 |
| production | 7 | 2.3e-05 | 1.7e-05 | 1.00 | 1.00 |
| hr1e-08 | 7 | 0 | 0 | 1.07 | 0.99 |

Segment 1.5e9 → 1.9e9 from the plateau state (relative change of the final state from the reference setting):

| setting | rhs calls | ln A₂₁₁ (abs) | ln A₃₂₂ (abs) | b (rel) | j (rel) | wall s (8 thr) |
|---|---|---|---|---|---|---|
| hr0.0001 | 31 | 2.0e-12 | 2.0e-09 | 8.8e-14 | 3.5e-13 | 8.39 |
| hr1e-05 | 31 | 2.0e-12 | 1.3e-09 | 8.9e-14 | 3.5e-13 | 3.17 |
| production | 31 | 2.1e-12 | 6.3e-10 | 9.4e-14 | 3.7e-13 | 3.24 |
| hr1e-08 | 31 | 0.0e+00 | 0.0e+00 | 0.0e+00 | 0.0e+00 | 3.39 |

### Study `stretch`

| state | setting | nodes | Re ω (max) | Im ω₂₁₁ | Im ω₃₂₂ | d ln A₂₁₁/dt | d ln A₃₂₂/dt | db/dt | dj/dt | amplitude-rate error / fastest rate | t(1 thr) s | t(8 thr) s | cost vs prod (1/8 thr) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| linear | sf1.06 | 677 | 6.9e-13 | 5.0e-11 | 3.7e-10 | 2.3e-09 | 6.9e-09 | 2.0e-07 | 2.0e-07 | 2.3e-09 | 0.407 | 0.164 | 1.15 / 1.08 |
| linear | production | 690 | 3.6e-13 | 1.4e-11 | 2.4e-10 | 7.7e-10 | 2.4e-09 | 1.4e-07 | 1.4e-07 | 7.7e-10 | 0.355 | 0.152 | 1.00 / 1.00 |
| linear | sf1.015 (ref) | 715 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.374 | 0.174 | 1.05 / 1.15 |
| plateau | sf1.06 | 677 | 1.4e-12 | 3.0e-04 | 3.1e-07 | 3.1e-04 | 3.2e-07 | 2.9e-04 | 2.9e-04 | 3.2e-07 | 0.331 | 0.126 | 1.02 / 0.85 |
| plateau | production | 690 | 1.1e-12 | 2.3e-04 | 2.2e-07 | 2.3e-04 | 2.2e-07 | 2.2e-04 | 2.1e-04 | 2.2e-07 | 0.326 | 0.149 | 1.00 / 1.00 |
| plateau | sf1.015 (ref) | 715 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.361 | 0.161 | 1.11 / 1.08 |
| takeover | sf1.06 | 677 | 1.2e-12 | 7.6e-08 | 3.3e-07 | 8.4e-08 | 3.3e-07 | 2.1e-07 | 1.4e-07 | 3.3e-07 | 0.358 | 0.176 | 1.02 / 1.08 |
| takeover | production | 690 | 9.0e-13 | 5.5e-08 | 2.3e-07 | 6.1e-08 | 2.3e-07 | 1.5e-07 | 1.0e-07 | 2.3e-07 | 0.352 | 0.163 | 1.00 / 1.00 |
| takeover | sf1.015 (ref) | 715 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.366 | 0.188 | 1.04 / 1.15 |
| nonlinear | sf1.06 | 677 | 6.5e-13 | 8.8e-10 | 1.1e-08 | 1.0e-09 | 1.1e-08 | 5.6e-07 | 5.8e-07 | 1.0e-09 | 0.354 | 0.168 | 1.00 / 1.09 |
| nonlinear | production | 690 | 3.5e-13 | 3.5e-10 | 8.6e-09 | 1.7e-12 | 4.9e-09 | 4.0e-07 | 4.1e-07 | 3.8e-11 | 0.356 | 0.154 | 1.00 / 1.00 |
| nonlinear | sf1.015 (ref) | 715 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.408 | 0.180 | 1.15 / 1.17 |
| large_mesh | sf1.06 | 1336 | 3.6e-14 | 8.6e-11 | 8.6e-10 | 1.0e-09 | 6.4e-09 | 2.9e-09 | 1.5e-09 | 1.0e-09 | 1.580 | 0.648 | 1.23 / 1.34 |
| large_mesh | production | 1355 | 2.1e-14 | 4.7e-11 | 6.0e-10 | 5.3e-10 | 4.6e-09 | 2.1e-09 | 1.1e-09 | 5.3e-10 | 1.288 | 0.482 | 1.00 / 1.00 |
| large_mesh | sf1.015 (ref) | 1393 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1.022 | 0.396 | 0.79 / 0.82 |
| strong_peak | sf1.06 | 677 | 7.0e-11 | 6.1e-09 | 3.7e-07 | 2.7e-08 | 3.7e-07 | 1.8e-07 | 1.8e-07 | 3.7e-07 | 0.387 | 0.149 | 1.12 / 1.11 |
| strong_peak | production | 690 | 4.9e-11 | 4.3e-09 | 2.6e-07 | 2.0e-08 | 2.6e-07 | 1.3e-07 | 1.3e-07 | 2.6e-07 | 0.347 | 0.135 | 1.00 / 1.00 |
| strong_peak | sf1.015 (ref) | 715 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.364 | 0.143 | 1.05 / 1.06 |
| strong_takeover | sf1.06 | 677 | 4.3e-11 | 2.1e-07 | 4.7e-07 | 2.2e-07 | 5.2e-07 | 3.9e-07 | 1.9e-07 | 3.8e-07 | 0.370 | 0.166 | 0.99 / 1.04 |
| strong_takeover | production | 690 | 3.1e-11 | 1.5e-07 | 3.3e-07 | 1.6e-07 | 3.6e-07 | 2.8e-07 | 1.3e-07 | 2.6e-07 | 0.375 | 0.160 | 1.00 / 1.00 |
| strong_takeover | sf1.015 (ref) | 715 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.401 | 0.160 | 1.07 / 1.00 |

Summary over states (largest change over the states listed; cost ratio to production, mean over those states):

| setting | states | amplitude-rate error / fastest rate | db/dt, dj/dt (relative) | cost 1 thr | cost 8 thr |
|---|---|---|---|---|---|
| sf1.06 | 7 | 3.8e-07 | 2.9e-04 | 1.07 | 1.08 |
| production | 7 | 2.6e-07 | 2.2e-04 | 1.00 | 1.00 |
| sf1.015 | 7 | 0 | 0 | 1.04 | 1.06 |

Segment 1.5e9 → 1.9e9 from the plateau state (relative change of the final state from the reference setting):

| setting | rhs calls | ln A₂₁₁ (abs) | ln A₃₂₂ (abs) | b (rel) | j (rel) | wall s (8 thr) |
|---|---|---|---|---|---|---|
| sf1.06 | 31 | 3.7e-11 | 3.0e-08 | 1.7e-12 | 6.5e-12 | 9.36 |
| production | 31 | 2.8e-11 | 2.1e-08 | 1.2e-12 | 4.9e-12 | 3.39 |
| sf1.015 | 31 | 0.0e+00 | 0.0e+00 | 0.0e+00 | 0.0e+00 | 3.45 |

### Study `truncation`

| state | setting | nodes | Re ω (max) | Im ω₂₁₁ | Im ω₃₂₂ | d ln A₂₁₁/dt | d ln A₃₂₂/dt | db/dt | dj/dt | amplitude-rate error / fastest rate | t(1 thr) s | t(8 thr) s | cost vs prod (1/8 thr) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| linear | production | 690 | 6.8e-11 | 7.5e-08 | 3.3e-09 | 7.5e-08 | 3.6e-09 | 1.2e-08 | 1.2e-08 | 7.5e-08 | 0.376 | 0.169 | 1.00 / 1.00 |
| linear | mo5 | 826 | 6.8e-11 | 7.5e-08 | 2.9e-09 | 7.6e-08 | 3.3e-09 | 1.6e-06 | 1.6e-06 | 7.6e-08 | 1.102 | 0.526 | 2.93 / 3.11 |
| linear | na2 (ref) | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.674 | 0.371 | 1.79 / 2.19 |
| plateau | production | 690 | 9.7e-10 | 1.7e-01 | 7.3e-04 | 1.7e-01 | 7.3e-04 | 1.7e-01 | 1.6e-01 | 7.3e-04 | 0.353 | 0.161 | 1.00 / 1.00 |
| plateau | mo5 | 826 | 9.6e-10 | 1.7e-01 | 7.2e-04 | 1.7e-01 | 7.2e-04 | 1.6e-01 | 1.6e-01 | 7.2e-04 | 1.006 | 0.510 | 2.85 / 3.18 |
| plateau | na2 (ref) | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.692 | 0.352 | 1.96 / 2.19 |
| takeover | production | 690 | 7.5e-10 | 3.3e-04 | 6.2e-04 | 3.5e-04 | 6.1e-04 | 2.3e-04 | 1.4e-05 | 6.1e-04 | 0.380 | 0.135 | 1.00 / 1.00 |
| takeover | mo5 | 826 | 7.4e-10 | 3.3e-04 | 6.2e-04 | 3.4e-04 | 6.0e-04 | 2.2e-04 | 1.5e-05 | 6.0e-04 | 1.058 | 0.521 | 2.78 / 3.85 |
| takeover | na2 (ref) | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.728 | 0.420 | 1.92 / 3.11 |
| nonlinear | production | 690 | 8.0e-11 | 1.3e-05 | 4.1e-05 | 1.3e-05 | 4.0e-05 | 1.9e-08 | 5.2e-07 | 1.3e-05 | 0.370 | 0.145 | 1.00 / 1.00 |
| nonlinear | mo5 | 826 | 8.0e-11 | 1.3e-05 | 4.1e-05 | 1.3e-05 | 4.0e-05 | 4.4e-06 | 5.0e-06 | 1.3e-05 | 1.105 | 0.517 | 2.98 / 3.57 |
| nonlinear | na2 (ref) | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.729 | 0.410 | 1.97 / 2.84 |
| large_mesh | production | 1355 | 2.1e-10 | 9.0e-05 | 5.4e-05 | 9.0e-05 | 7.4e-05 | 7.4e-07 | 7.1e-07 | 9.0e-05 | 1.401 | 0.546 | 1.00 / 1.00 |
| large_mesh | mo5 | 1578 | 2.1e-10 | 9.0e-05 | 5.4e-05 | 9.0e-05 | 7.4e-05 | 7.9e-07 | 7.4e-07 | 9.0e-05 | 3.095 | 1.547 | 2.21 / 2.83 |
| large_mesh | na2 (ref) | 1355 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2.017 | 1.134 | 1.44 / 2.08 |
| strong_peak | production | 690 | 7.1e-07 | 9.8e-04 | 2.0e-02 | 2.4e-05 | 2.0e-02 | 9.3e-04 | 4.1e-04 | 2.0e-02 | 0.394 | 0.160 | 1.00 / 1.00 |
| strong_peak | mo5 | 826 | 7.1e-07 | 9.8e-04 | 2.0e-02 | 2.3e-05 | 2.0e-02 | 9.3e-04 | 4.1e-04 | 2.0e-02 | 1.179 | 0.602 | 2.99 / 3.76 |
| strong_peak | na2 (ref) | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.703 | 0.395 | 1.78 / 2.46 |
| strong_takeover | production | 690 | 4.1e-07 | 1.5e-02 | 1.5e-02 | 1.5e-02 | 1.3e-02 | 1.6e-02 | 2.9e-04 | 1.5e-02 | 0.373 | 0.138 | 1.00 / 1.00 |
| strong_takeover | mo5 | 826 | 4.1e-07 | 1.5e-02 | 1.5e-02 | 1.5e-02 | 1.3e-02 | 1.6e-02 | 2.9e-04 | 1.5e-02 | 1.061 | 0.605 | 2.85 / 4.37 |
| strong_takeover | na2 (ref) | 690 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.683 | 0.401 | 1.83 / 2.90 |

Summary over states (largest change over the states listed; cost ratio to production, mean over those states):

| setting | states | amplitude-rate error / fastest rate | db/dt, dj/dt (relative) | cost 1 thr | cost 8 thr |
|---|---|---|---|---|---|
| production | 7 | 2.0e-02 | 1.7e-01 | 1.00 | 1.00 |
| mo5 | 7 | 2.0e-02 | 1.6e-01 | 2.80 | 3.53 |
| na2 | 7 | 0 | 0 | 1.81 | 2.54 |

Segment 1.5e9 → 1.9e9 from the plateau state (relative change of the final state from the reference setting):

| setting | rhs calls | ln A₂₁₁ (abs) | ln A₃₂₂ (abs) | b (rel) | j (rel) | wall s (8 thr) |
|---|---|---|---|---|---|---|
| production | 31 | 2.0e-08 | 6.9e-05 | 9.1e-10 | 3.6e-09 | 9.00 |
| mo5 | 31 | 2.0e-08 | 6.8e-05 | 8.9e-10 | 3.5e-09 | 10.21 |
| na2 | 31 | 0.0e+00 | 0.0e+00 | 0.0e+00 | 0.0e+00 | 7.46 |

### Study `combined`

| state | setting | nodes | Re ω (max) | Im ω₂₁₁ | Im ω₃₂₂ | d ln A₂₁₁/dt | d ln A₃₂₂/dt | db/dt | dj/dt | amplitude-rate error / fastest rate | t(1 thr) s | t(8 thr) s | cost vs prod (1/8 thr) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| linear | production | 690 | 8.6e-12 | 9.6e-07 | 1.2e-05 | 5.1e-08 | 5.4e-07 | 1.3e-06 | 1.3e-06 | 5.1e-08 | 0.396 | 0.174 | 1.00 / 1.00 |
| linear | np24 npol16 h0.3 sf1.015 c1e-06 hr1e-08 (ref) | 1027 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1.292 | 0.301 | 3.27 / 1.73 |
| plateau | production | 690 | 1.5e-11 | 3.2e-03 | 1.6e-06 | 3.2e-03 | 8.9e-07 | 3.1e-03 | 3.0e-03 | 8.9e-07 | 0.360 | 0.133 | 1.00 / 1.00 |
| plateau | np24 npol16 h0.3 sf1.015 c1e-06 hr1e-08 (ref) | 1027 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1.232 | 0.300 | 3.42 / 2.26 |
| takeover | production | 690 | 1.4e-11 | 7.9e-07 | 1.7e-06 | 5.3e-09 | 9.0e-07 | 6.7e-07 | 4.8e-07 | 9.0e-07 | 0.384 | 0.138 | 1.00 / 1.00 |
| takeover | np24 npol16 h0.3 sf1.015 c1e-06 hr1e-08 (ref) | 1027 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1.371 | 0.293 | 3.57 / 2.12 |
| nonlinear | production | 690 | 8.5e-12 | 9.1e-07 | 1.3e-05 | 1.1e-07 | 7.7e-07 | 5.5e-07 | 1.4e-06 | 1.1e-07 | 0.394 | 0.138 | 1.00 / 1.00 |
| nonlinear | np24 npol16 h0.3 sf1.015 c1e-06 hr1e-08 (ref) | 1027 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1.328 | 0.304 | 3.37 / 2.20 |
| large_mesh | production | 1355 | 2.3e-13 | 1.4e-06 | 1.3e-05 | 4.3e-06 | 5.1e-06 | 1.3e-06 | 1.5e-06 | 4.3e-06 | 1.332 | 0.521 | 1.00 / 1.00 |
| large_mesh | np24 npol16 h0.3 sf1.015 c1e-06 hr1e-08 (ref) | 2041 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 4.250 | 1.125 | 3.19 / 2.16 |
| strong_peak | production | 690 | 3.3e-10 | 1.0e-06 | 1.7e-06 | 1.3e-06 | 2.2e-05 | 1.5e-06 | 1.2e-06 | 2.2e-05 | 0.356 | 0.153 | 1.00 / 1.00 |
| strong_peak | np24 npol16 h0.3 sf1.015 c1e-06 hr1e-08 (ref) | 1027 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1.211 | 0.354 | 3.40 / 2.32 |
| strong_takeover | production | 690 | 1.2e-10 | 6.1e-07 | 1.5e-06 | 8.2e-05 | 4.4e-04 | 1.2e-05 | 3.5e-06 | 3.2e-04 | 0.368 | 0.151 | 1.00 / 1.00 |
| strong_takeover | np24 npol16 h0.3 sf1.015 c1e-06 hr1e-08 (ref) | 1027 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1.287 | 0.393 | 3.50 / 2.60 |

Summary over states (largest change over the states listed; cost ratio to production, mean over those states):

| setting | states | amplitude-rate error / fastest rate | db/dt, dj/dt (relative) | cost 1 thr | cost 8 thr |
|---|---|---|---|---|---|
| production | 7 | 3.2e-04 | 3.1e-03 | 1.00 | 1.00 |
| np24 npol16 h0.3 sf1.015 c1e-06 hr1e-08 | 7 | 0 | 0 | 3.39 | 2.20 |

Segment 1.5e9 → 1.9e9 from the plateau state (relative change of the final state from the reference setting):

| setting | rhs calls | ln A₂₁₁ (abs) | ln A₃₂₂ (abs) | b (rel) | j (rel) | wall s (8 thr) |
|---|---|---|---|---|---|---|
| production | 31 | 3.9e-10 | 8.4e-08 | 1.7e-11 | 6.9e-11 | 8.12 |
| np24 npol16 h0.3 sf1.015 c1e-06 hr1e-08 | 31 | 0.0e+00 | 0.0e+00 | 0.0e+00 | 0.0e+00 | 6.26 |
