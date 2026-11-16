# How the radial solve scales with the number of bound modes

Scoping study, November 13, 2026, branch `speedup/scaling`. Question: the method handles two bound modes
today (the $|211\rangle$ and $|322\rangle$ clouds); what does one radial solve cost with $K$ bound modes,
which parts make $K \approx 10$ impossible, and what would have to change? This is a study with small
measured prototypes, not an implementation of large $K$.

All timings are from Slurm jobs on lr7 (Intel Xeon Gold 6330), Julia 1.13.0, one BLAS thread; job numbers
are given with each table. Scripts and logs are in `axion-sr-sim-speedup/scaling-work/` (outside the
repository; listed at the end).

## If you read only this

* **Measured cost now.** One warm solve (about one right-hand-side evaluation) takes 0.10 s at $K = 2$, 1.5 s at
  $K = 3$ and, from measured pieces, about 47 s at $K = 4$ (8 threads, 690-node mesh). Adding $|433\rangle$ also
  makes the mesh six times longer, because a combination frequency ($2\omega_{433} - \omega_{322}$) lies just below
  the mass threshold.
* **Three walls.** The tensor phase grid ($16^K$ points), the dense block LU ($(4C)^3$ per node with $C \sim K^3$
  components, so $K^9$ time and $K^6$ memory), and the size of the component set itself.
* **Lattice quadrature works.** A rank-1 lattice rule reproduces the source and its Jacobian to $10^{-10}$ with
  about 100 points at $K = 2$ and 700 at $K = 3$ (against 256 and 4096), for the cosine and chiral potentials at
  the field amplitudes of the clouds ($\Phi \lesssim 0.1$). The number of points needed grows polynomially in $K$
  (about $K^6$), not exponentially. It is in the library behind the `phase_grid` keyword; the tensor grid stays
  the default.
* **Most components do not matter.** The components with frequency near $3\mu$ contribute below $10^{-12}$ of the
  rates; dropping them changes the rates by $10^{-8}$-$10^{-5}$. What matters are the fundamentals and the
  $2\to2$ fourth legs between populated modes: 6 of 22 components at $K = 3$ for $10^{-4}$ accuracy.
* **Newton–Krylov works.** GMRES preconditioned by the component-diagonal linear operator solves a Newton step in
  2-5 iterations at the clouds' field amplitudes (12 at a field six times larger), at $K = 2$ and 3 alike, with
  matrix-free products from one dual-number residual each. It removes the dense-LU wall but only pays beyond
  $K \approx 5$.
* **The paper (arXiv:2611.15865)** finds two dominant states at the production parameters ($\alpha = 0.2$-0.4,
  $f_a \approx 6\times10^{17}$ GeV) and at most a few weak ones; four or five dominant states at late times in
  richer regions, with tens of weakly populated states that barely feed back.
* **Recommendation.** Do not aim at order 10 fully nonlinear modes (hours per evaluation even with every change
  above). Aim at 2-4 strong modes treated as now, plus tens of weak modes, each a separate linear problem on the
  strong background (linear cost, parallel; the old branch's primary/secondary construction), chosen per $\alpha$
  from the paper's Fig. 2 and promoted when they grow. First steps: the lattice in full solves; a small-amplitude
  check of our $(2,-1)$ flux against the paper's $|211\rangle\times|211\rangle\to|322\rangle\times\mathrm{BH}$ rate;
  then the weak-mode solver.

## 1. What one radial solve does, and where $K$ enters

One radial solve is a Newton iteration on the discretised boundary-value problem. Each Newton iteration
needs three things:

* the **residual**: at every mesh node and interval midpoint, the right-hand side of the radial equations,
  including the nonlinear source $N_{l,q}$, which is the phase average
  $\langle e^{+i q\cdot\xi} h(\Phi)\rangle_\xi$ of the nonlinearity over the torus of independent phases
  $\xi = (\xi_1,\dots,\xi_K)$;
* the **Jacobian** (matrix of derivatives of the residual), whose nonlinear part needs the Fourier
  coefficients $\hat h'(k) = \langle e^{+ik\cdot\xi} h'(\Phi)\rangle_\xi$ at every wavevector $k = q \pm q'$
  formed from two retained components;
* the **linear solve** with that Jacobian (block LU factorisation, then back-substitution).

Three numbers grow with $K$:

* the number of retained Fourier components $C$ (with one polar harmonic per component, also the number
  of channels). With the default truncation (odd $|q|_1 \le 3$) $C = 2K + 2K(K-1) + 4\binom{K}{3}$:
  8, 22, 48, 152, 680 for $K$ = 2, 3, 4, 6, 10, cubic in $K$;
* the number of phase points $P$. The present rule is a tensor grid with 16 points per phase,
  $P = 16^K$, exponential in $K$;
* the number of wavevector slots $S$ (distinct $q \pm q'$ up to sign): 25, 116, 417, 3161 for
  $K$ = 2, 3, 4, 6, about $K^6$.

Per evaluation point (a mesh node or a midpoint) the residual costs about $4CPX$ multiply-adds
($X = 12$ polar nodes: synthesis of $\Phi$ from $2C$ real profiles, then projection back), and the source
Jacobian about $2SPX$ more. The interval blocks of the Jacobian are dense matrices of side $n = 4C$
(real and imaginary parts of $R$ and of its slope $Q$ per channel); forming them takes two
$n\times n$ by $n\times n$ products per interval, and the block LU eliminates one such block per node.

**Why dense block LU is cubic.** Gaussian elimination on a dense $n \times n$ block subtracts a multiple
of each pivot row from every row below it, for each of the $n$ pivots: about $n$ rows times $n$ columns
times $n$ pivots, so $n^3$ operations, and the block must be stored ($n^2$ numbers). Our Jacobian is
block-bidiagonal (each interval couples two neighbouring nodes), so the work is $N_\mathrm{nodes}$ such
eliminations: cost $\propto N_\mathrm{nodes} (4C)^3$ and memory $\propto N_\mathrm{nodes} (4C)^2$. With
$C \propto K^3$ that is $K^9$ in time and $K^6$ in memory.

A second, independent effect: adding a bound mode can lengthen the mesh. The mesh must contain the
slowest-decaying bound component, and combination frequencies such as $2\omega_{433} - \omega_{322}$ can
sit just below the mass threshold $\alpha$, with a very long decay length (Section 2).

## 2. Item 1: the present code at $K = 3$ and $K = 4$, measured

### Which third mode

At $\alpha = 0.3$, $a = 0.9$ the linear growth rates are $\mathrm{Im}\,\omega$ = 1.48e-8 ($|211\rangle$),
5.0e-9 ($|311\rangle$), 1.44e-12 ($|322\rangle$) and 1.5e-17 ($|433\rangle$) (from the linear solves below).
Two third modes were measured:

* $|311\rangle$: the fastest-growing state after $|211\rangle$ at this point (same $m = 1$, so it stops
  growing once $|211\rangle$ has spun the hole down). It keeps the $K = 2$ mesh (690 nodes), so it isolates
  the effect of $K$ on the cost at fixed mesh.
* $|433\rangle$: the next $\ell = m$ level, the one the literature (Section 4) finds populated by
  self-interactions at $\alpha \ge 0.2$ for $f_a \lesssim 10^{16}$ GeV. Its combination component
  $q = (0, -1, 2)$, frequency $2\omega_{433} - \omega_{322} \approx \alpha(1 - 0.007\alpha^2)$, is bound but
  close to threshold (decay constant $0.0106$ against $0.03$-$0.04$ for the others), so
  `RadialMesh(background, α, mode_set)` asks for a domain out to $r^* \approx 2466$ and 4256 nodes, six
  times the $K = 2$ mesh. On the $K = 2$ mesh the solve stalls at a residual of 1.8e-10, above the default
  tolerance 1e-10.

For $K = 4$ the modes were $|211\rangle, |322\rangle, |311\rangle, |433\rangle$ on the 690-node mesh
(pieces only, at the hydrogenic initial guess: the cost of a residual, a Jacobian or a factorisation does not
depend on the state) and $|211\rangle, |322\rangle, |433\rangle, |544\rangle$ on its own 4256-node mesh.

### Measured cost

Seconds per call (minimum of 3 calls after a warm-up; solves: minimum of 2), lr7. Jobs 27824080 (1 thread),
27824081 (8 threads), 27823936 ($K = 4$, 8 threads). "Warm solve" is `solve_radial` from the converged state
after multiplying the second amplitude by 1.01 (the harness's perturbed solve, without the background change);
it is the closest single number to the cost of one right-hand-side evaluation of the evolution. The nonlinear
state has the amplitudes of the harness's `nonlinear` state ($A_{211} = 1.47\times10^{-4}$,
$A_{322} = 4.1\times10^{-6}$; third mode $10^{-6}$), reached by 30 continuation steps on the $a = 0.9$ background.

| | $K = 2$ | $K = 3$ ($|311\rangle$) | $K = 3$ ($|433\rangle$, 690 nodes) | $K = 4$ (690 nodes) | $K = 4$ ($|433\rangle,|544\rangle$, 4256 nodes) |
|---|---|---|---|---|---|
| components $C$ / node block $4C$ | 8 / 32 | 22 / 88 | 22 / 88 | 48 / 192 | 48 / 192 |
| phase points $P = 16^K$ / slots $S$ | 256 / 25 | 4096 / 116 | 4096 / 116 | 65536 / 417 | 65536 / 417 |
| unknowns | 22112 | 60808 | 60808 | 132672 | 817344 |
| residual, 1 thread | 0.021 | 0.63 | | | |
| Jacobian blocks, 1 thread | 0.072 | 1.99 | | | |
| factorisation, 1 thread | 0.018 | 0.18 | | | |
| warm solve (nonlinear), 1 thread | 0.307 (4 it.) | 7.65 (4 it.) | | | |
| residual, 8 threads | 0.0039 | 0.102 | 0.099 | 3.25 | 19.6 |
| Jacobian blocks, 8 threads | 0.012 | 0.301 | 0.301 | 14.4 | 82.2 |
| factorisation, 8 threads | 0.016 | 0.149 | 0.147 | 1.10 | 6.71 |
| back-substitution, 8 threads | 0.0025 | 0.014 | 0.013 | 0.055 | 0.33 |
| warm solve (linear), 8 threads | 0.088 (3 it.) | 1.35 (3 it.) | 1.39 (3 it.)\* | | |
| warm solve (nonlinear), 8 threads | 0.101 (4 it.) | 1.47 (4 it.) | | | |
| cold solve (linear, hydrogenic start), 8 threads | 0.18 | 1.99 | 4.2, not converged\* | | |
| peak memory (whole process) | 1.7 GB | 2.7 GB | 2.8 GB | 4.5 GB | 14.2 GB |
| Jacobian workspace + factorisation | 37 MB | 268 MB | 268 MB | 1.26 GB | 7.7 GB |

\* On the 690-node mesh the $|433\rangle$ problem stalls at a scaled residual of 1.8e-10, above the default
tolerance $10^{-10}$, because its near-threshold component $(0,-1,2)$ is not contained (Section 2); the
continuation to the nonlinear state therefore stops at its first step. On its own 4256-node mesh (job 27824081, `scaling-work/logs/k8-27824081.out`) the cold solve from the hydrogenic guess failed: the block factorisation hit an exactly singular pivot (zero pivot in a triangular solve) in the first Newton iterations, after about 40 minutes. So there is no measured $K = 3$ case with $|433\rangle$ on an adequate mesh; the cause (the near-threshold component's scaling, or the hydrogenic start for $|433\rangle$, whose single-mode frequency solves also failed to converge) was not investigated.

What this establishes: the cost of the present code at $K = 3$ and $K = 4$ on this hardware, per piece, for the
default quadrature (tensor, $n_\mathrm{phase} = 16$) and truncation. What it does not establish: anything about
accuracy at $K = 3$ (no convergence study in $n_\mathrm{phase}$, mesh or truncation was done for the
$K = 3$ states), and the $K = 4$ rows are single pieces at the initial guess, not solves.

The one-thread numbers show the growth of the work: from $K = 2$ to $K = 3$ the residual costs 30 times more,
the Jacobian 28 times, the factorisation 10 times, a warm solve 25 times. Eight threads reduce the $K = 3$
warm solve from 7.7 s to 1.5 s (5.2 times). An evolution like the production one (2233 right-hand-side
calls for the window $1.5\times10^9 \to 10^{11}$, 295 s on 8 threads at $K = 2$) would take about an hour
at $K = 3$ on the same mesh (2233 times 1.5 s), and about a week at $K = 4$ on its 4256-node mesh (a warm
solve there is about five residuals and two Jacobians with factorisations, $5 \times 19.6 + 2 \times 89 \approx
275$ s; extrapolated, not measured).

### Against the cost model

The cost model of Section 1 predicts the ratios between measured cases (all 8 threads, 690 nodes):

| Piece | model | $K = 2 \to 3$ predicted | measured | $K = 3 \to 4$ predicted | measured |
|---|---|---|---|---|---|
| residual | $\propto CP$ | 44 | 26 | 35 | 32 |
| Jacobian | $\propto SP$ (+ dense $2(4C)^3$) | 74 | 25 | 57 | 48 |
| factorisation | $\propto (4C)^3$ | 21 | 9.5 | 10.4 | 7.4 |

From $K = 3$ to $K = 4$ the model holds to within 10-25%: the phase-point work ($CP$ for the residual, $SP$
for the source Jacobian) dominates and grows by 35-57 times per added mode, and the factorisation follows
$(4C)^3$ once the blocks are large enough for dense linear algebra to run at full speed. At $K = 2$ the
pieces are small enough that fixed overheads (boundary blocks, threading, allocation) hide part of the
growth, so the $2 \to 3$ ratios fall below the model. Taking the $K = 4$ measurement as the anchor, the cost
per mesh node is about $1.5\times10^{-9}$ s per unit of $CP$ for the residual, $7.6\times10^{-10}$ s per unit
of $SP$ for the Jacobian (8 threads), and $2.3\times10^{-10}$ s per unit of $(4C)^3$ for the factorisation
(serial). These constants are used for the extrapolations in Section 7.

The memory of the Jacobian blocks and the factorisation is $6 (4C)^2 \times 8$ bytes per node (three
$n \times n$ arrays in the Jacobian workspace and the factorisation panels of $1.5n \times 2n$): at $K = 4$
on 4256 nodes this gives 7.5 GB against 7.7 GB measured.

## 3. Item 2: a phase quadrature that is not a tensor grid

### The idea: what the phase average must deliver

The phase average is not needed as a function on the torus. The residual needs only the Fourier
coefficients $\hat h(q) = \langle e^{+iq\cdot\xi} h(\Phi)\rangle_\xi$ at the $C$ retained components, and the
Jacobian the coefficients $\hat h'(k)$ at the $S$ slot wavevectors. Any equal-weight point set that computes
those averages accurately will do; the kernels in `src/field.jl` already work on plain point lists
(`synthesis`, `analysis`, `jacobian_table`), so a different point set needs no change to them.

**The rank-1 lattice rule.** Choose an integer $M$ (the number of points) and an integer generating vector
$z = (z_1,\dots,z_K)$. The points are $\xi_j = 2\pi (j z \bmod M)/M$ for $j = 0,\dots,M-1$: one walks along a
straight line on the torus in direction $z$ and wraps around. The average of a single exponential over these
points is
$$\frac1M\sum_{j=0}^{M-1} e^{i k\cdot\xi_j} = \begin{cases}1 & k\cdot z \equiv 0 \pmod M\\ 0 & \text{otherwise,}\end{cases}$$
the discrete orthogonality of the $M$-th roots of unity, because $e^{ik\cdot\xi_j} = e^{2\pi i\, j (k\cdot z)/M}$.
So the rule integrates every exponential exactly except the wavevectors $k \neq 0$ with $k\cdot z \equiv 0
\pmod M$. The tensor grid with $n$ points per phase is the special case that is exact for every $k$ with all
$|k_i| < n$; the lattice instead protects a ball $|k|_1 \le L$ of wavevectors, which is the shape that matters
here, and it needs far fewer points to do so.

**Aliasing.** The source coefficient at target $q$ is computed as $\frac1M\sum_j e^{iq\cdot\xi_j}h(\Phi_j)$.
Expanding $h(\Phi) = \sum_{k'}\hat h(k') e^{-ik'\cdot\xi}$, the rule returns $\hat h(q)$ plus every
$\hat h(k')$ with $(k' - q)\cdot z \equiv 0 \pmod M$, $k' \ne q$: those harmonics are aliases of $q$, they are
indistinguishable from it on the point set. For a polynomial nonlinearity the harmonics stop at a finite order
(for $h = \Phi^3/6$ built from components with $|q|_1 \le 3$, at $|k'|_1 \le 9$; for $h' = \Phi^2/2$ at 6), so
a lattice whose smallest alias $d \neq 0$ (smallest $|d|_1$ with $d\cdot z\equiv 0$) exceeds 12 is exact.
The cosine and chiral nonlinearities have harmonics of every order, and the error of the rule is the size of
the harmonics that reach an alias. Call the smallest alias order of a lattice $L_a$ (my notation): the error
of the source is roughly the largest $|\hat h(k')|$ with $|k'|_1 \ge L_a - 3$, and that of the Jacobian table
the largest $|\hat h'(k')|$ with $|k'|_1 \ge L_a - 6$. How many points are needed therefore depends on how
fast the harmonics of $h(\Phi)$ decay at the field amplitudes that occur.

### How fast the harmonics decay

The prototype (`scaling-work/lattice_proto.jl`) takes the profiles at the mesh node where the field is
largest, rescales them so that the peak field is $\Phi_\mathrm{max}$, and computes the largest $|\hat h(k)|$
in each shell $|k|_1$ on a fine tensor grid. Relative to the largest coefficient (shell 1), $K = 2$ nonlinear
state, cosine potential (job 27824280):

| $\Phi_\mathrm{max}$ | shell 3 | 5 | 7 | 9 | 11 | 13 | 15 |
|---|---|---|---|---|---|---|---|
| 0.0034 (the state itself) | 0.3 | 9e-6 | 8e-11 | below roundoff | | | |
| 0.1 | 0.3 | 3e-5 | 2e-9 | 5e-14 | | | |
| 0.5 | 0.3 | 7e-4 | 1e-6 | 8e-10 | 4e-13 | | |
| 1 | 0.3 | 3e-3 | 2e-5 | 5e-8 | 1e-10 | 1e-13 | |
| 3 | 0.3 | 2e-2 | 1e-3 | 4e-5 | 7e-7 | 1e-8 | 1e-10 |

(Even shells vanish because $h$ is odd.) Shell 3 is the cubic term itself; beyond it each two shells cost a
factor that falls with the field, roughly $(\Phi_\mathrm{max}/15)^2$ (within a factor of about 3 over this range). The $K = 3$
state gives the same numbers to within a factor of two.

**Field amplitudes that occur.** The peak field is $\Phi_\mathrm{max} = 0.022$ on the $|211\rangle$ plateau
of the production run ($\alpha = 0.3$, $a = 0.9$, $f_a = 5.7\times10^{17}$ GeV) and $0.0034$ in its
two-cloud `nonlinear` state. The saturation amplitude is set by the balance between superradiant growth and
self-interaction losses, and in that balance $f_a$ drops out of $\Phi$ (the occupation number scales as
$f_a^2$ and $\Phi^2 \propto N/f_a^2$), so these values are typical for $\alpha \approx 0.3$ at any $f_a$; they
grow with the superradiant rate, that is with $\alpha$. At $\Phi \lesssim 0.1$ the nonlinearity is a cubic
polynomial plus corrections of relative size $10^{-4}$, and the harmonics beyond shell 7 are below $10^{-9}$.

### Accuracy against the number of points

Error of the source and of its Jacobian at one node, max-norm relative to the reference (tensor grid with 64
points per phase at $K = 2$, 32 at $K = 3$; tensor 48 against 64 at $K = 2$ and 24 against 32 at $K = 3$ agree
to roundoff except where noted). Lattice sizes from the list $M$ = 31, 61, 101, 151, 211, 307, 401, 503, 701,
1009 ($K = 2$) and 211, 401, 701, 1009, 1511, 2003, ... ($K = 3$); each with the Korobov generating vector
$z = (1, g, g^2, \dots) \bmod M$ that minimises a weighted count of aliases up to order 24 ($K = 2$) or 16
($K = 3$). Entries: the smallest $M$ in the list whose error is below the stated level; for comparison the
smallest tensor grid ($n^K$ points, $n$ = 8, 12, 16, 24, 32, 48) that does the same. Jobs 27824280 ($K = 2$),
27824279 ($K = 3$).

| $K$ | potential | $\Phi_\mathrm{max}$ | lattice, source $10^{-10}$ | lattice, Jacobian $10^{-10}$ | tensor, source $10^{-10}$ | tensor, Jacobian $10^{-10}$ |
|---|---|---|---|---|---|---|
| 2 | cosine | 0.0034 | 61 | 61 | 144 | 144 |
| 2 | cosine | 0.1 | 61 | 101 | 144 | 256 |
| 2 | cosine | 1 | 101 | 101 | 256 | 576 |
| 2 | cosine | 3 | 101 | 211 | 576 | 576 |
| 2 | chiral | 0.1 | 61 | 101 | 144 | 256 |
| 2 | chiral | 1 | 101 | 151 | 256 | 576 |
| 2 | chiral | 3 | 1009 | 1009 | 2304 | 2304 |
| 3 | cosine | 0.0034 | 211 | 401 | 1728 | 1728 |
| 3 | cosine | 0.1 | 211 | 701 | 1728 | 4096 |
| 3 | cosine | 1 | 701 | 1511 | 4096 | 4096 |
| 3 | cosine | 3 | 1511 | 2003 | 13824 | 13824 |
| 3 | chiral | 0.1 | 401 | 701 | 1728 | 4096 |
| 3 | chiral | 1 | 1009 | 1511 | 4096 | 13824 |
| 3 | chiral | 3 | not reached by 16001 | not reached | reference not converged | |

The $10^{-8}$ level needs the same or the next smaller lattice in every case. The default $n_\mathrm{phase} =
16$ ($256$ and $4096$ points) is accurate to $10^{-10}$ for $\Phi_\mathrm{max} \le 1$ (cosine); the
Jacobian at $\Phi_\mathrm{max} = 2$-3 is only good to $10^{-7}$-$10^{-6}$ there, which does not matter at
physical amplitudes.

At the amplitudes that occur ($\Phi_\mathrm{max} \le 0.1$) a lattice of about 100 points ($K = 2$) or 700
points ($K = 3$) does what the default tensor grid does with 256 and 4096 points, for both potentials: 2.5
and 6 times fewer points, so the residual and the source Jacobian, whose cost is proportional to the number
of points, become 2.5 and 6 times cheaper. The gain grows with $K$ (next subsection).

What this establishes: the lattice rule reproduces the tensor-grid source and Jacobian at one node of two
real states, for two potentials and a range of rescaled amplitudes, with the number of points given. What it
does not establish: a full solve or evolution with the lattice (not run); other nodes (the peak-field node is
the hardest, since the harmonics decay fastest where the field is small); the rescaled profiles at large
$\Phi_\mathrm{max}$ are not solutions of the field equations, only test inputs with the right harmonic content.

### Does the number of points grow polynomially with $K$?

A lattice with smallest alias order $L_a$ must map every wavevector with $|k|_1 \le L_a/2$ to a different
residue (two of them with equal residues would differ by an alias of order at most $L_a$), so it needs at least
as many points as there are such wavevectors:
$$M \ge B(K, r) = \sum_{j=0}^{\min(K,r)} 2^j\binom{K}{j}\binom{r}{j}, \qquad r = \lfloor L_a/2\rfloor,$$
the number of integer vectors in $K$ dimensions with $|k|_1 \le r$. For fixed $r$ this is a polynomial in $K$
of degree $r$ (about $(2K)^r/r!$ for $K \gg r$), whereas the tensor grid exact to the same order needs
$(L_a+1)^K$ points. The required $L_a$ is set by the decay of the harmonics, which depends on the field
amplitude and not on $K$: $L_a \approx 12$ for $10^{-10}$ in the Jacobian at $\Phi \lesssim 0.1$ (the
measured $K = 3$ case needed $L_a = 14$ for $10^{-10}$ in the Jacobian and $L_a = 9$ for the source).

The Korobov search finds lattices close to that bound (smallest $M$ in a scan with no alias of order $\le 12$;
`scaling-work/min_lattice.jl`, job 27824280 for $K = 4, 5$, login node for $K = 2, 3$):

| $K$ | bound $B(K,6)$ | smallest exact Korobov lattice found | tensor grid exact to order 12 ($13^K$) |
|---|---|---|---|
| 2 | 85 | 85 | 169 |
| 3 | 377 | 455 | 2197 |
| 4 | 1289 | 2398 | 28561 |
| 5 | 3653 | 8609 | 371293 |
| 6 | 8989 | (not searched) | $4.8\times10^6$ |
| 10 | 134245 | (not searched) | $1.4\times10^{11}$ |

(For $K = 4, 5$ only 4000 random multipliers $g$ were tried per $M$; a component-by-component construction
would do better.) So the answer is yes: under this approach the number of points grows polynomially with $K$,
like $K^6$ for the accuracy needed here, against $16^K$ now. The evidence is the measured $K = 2, 3$ accuracy
together with the lower bound and the constructed lattices for $K \le 5$. The polynomial is steep: at $K = 10$
a $10^{-10}$ Jacobian table needs of the order of $10^5$ points, which with $C = 680$ components and
$S = 52441$ slots is still far too expensive (Section 7). The lattice removes the exponential wall but not the
cost of a large component set.

**Alternatives not prototyped.** Sparse grids and hyperbolic-cross rules protect a different set of
wavevectors (products of one-dimensional degrees rather than the total degree $|k|_1$) and are a worse match
to harmonics that are organised by $|k|_1$. For the cosine potential, the Jacobi–Anger expansion
$e^{iA\cos\theta} = \sum_n i^n J_n(A)e^{in\theta}$ turns $\sin\Phi$ into a product over components of Bessel
series, so its Fourier coefficients can be built by successive convolutions on the wavevector lattice, one
component at a time, keeping only wavevectors with non-negligible coefficients; the cost is about
(components) × (kept wavevectors) × (Bessel orders), with the kept set again a ball of order $L_a$, so it scales
like the lattice and works only for the cosine. (The old branch's "analytic phase averaging" of November 2,
`HANDOFF.md`, was for the quadratic observables, not for $h(\Phi)$; it kept a quadrature for the potential.)

### Code added

`PhaseGrid(points, mode_set)` builds the tables for any equal-weight point list (the tensor constructor now
goes through it, with bitwise-identical tables), `rank1_lattice(z, M)` gives the lattice points, and
`FieldQuadrature` and `RadialProblem` accept `phase_grid = PhaseGrid(rank1_lattice(z, M), mode_set)`. The
tensor grid stays the default. Unit test (in `test/runtests.jl`, "Rank-1 lattice phase rule"): for the quartic
potential an 85-point lattice gives the same source and Jacobian as the 256-point tensor grid to $10^{-13}$.
The generating-vector search is in the scratch scripts, not in the library.

## 4. Which bound modes to evolve: Cannizzaro, Spieksma and Witte (arXiv:2611.15865)

All equation, figure and table numbers below are from the PDF in `docs/papers/2611.15865.pdf`.

**What the paper does.** It treats the quartic self-interaction to first order: each $2\to2$ process
$|a\rangle\times|b\rangle\to|c\rangle\times X$ has a rate computed from relativistic Kerr wavefunctions and the
sourced radial equation (their Eqs. 12-17), where the fourth leg $X$ is a forced oscillation absorbed by the
hole (BH) if $\omega_a + \omega_b - \omega_c < \mu$ and an axion radiated to infinity otherwise. The occupation
numbers then obey a nonlinear rate system (Eq. 21; its two-level version is Eqs. 22-25, with the processes
$|211\rangle\times|211\rangle\to|322\rangle\times\mathrm{BH}$ and $|322\rangle\times|322\rangle\to
|211\rangle\times\infty$). Frequency shifts and $3\to1$ processes are neglected (Sec. III A, items (ii)-(iii)),
and only superradiant ($m > 0$) states are kept. A state counts as populated if $\epsilon_{n\ell m} =
N_{n\ell m}/(GM^2)$ exceeds $\epsilon_\mathrm{th} = 10^{-20}$ at some time. The grid is $\alpha \in \{0.05,
0.2, 0.4, 0.6, 0.8, 1.0, 1.3, 1.6\}$, $f_a \in \{10^{12}, 10^{14}, 10^{16}, 10^{18}\}$ GeV, initial spin
$\tilde a = 0.95$, $M = 10\,M_\odot$ over $5\times10^7$ yr (and $10^4$, $10^8\,M_\odot$ over 1 Gyr). Black-hole
convergence (mass and spin within 1% up to a 30% time shift) is reached at $n_\mathrm{max} = 8$ everywhere
(Sec. I, result 2).

### (a) The populated states at $\alpha = 0.2$ and $0.4$, and what matters

The production runs here are at $\alpha = 0.3$ (and $0.187$), $a = 0.9$, $f_a/M_\mathrm{Pl} = 0.0467$, that
is $f_a = 5.7\times10^{17}$ GeV (non-reduced $M_\mathrm{Pl} = 1.22\times10^{19}$ GeV). The self-interaction
strength scales as $(M_\mathrm{Pl}/f_a)^4$ (their Eq. 19), so this is 10 times stronger than their
$f_a = 10^{18}$ GeV column and $10^7$ times weaker than $10^{16}$ GeV.

From Fig. 2 (states with $\epsilon > 10^{-20}$, $M = 10\,M_\odot$), read at the resolution of the figure:

| $\alpha$ | $f_a = 10^{18}$ GeV | $f_a = 10^{16}$ GeV | $f_a = 10^{14}$, $10^{12}$ GeV |
|---|---|---|---|
| 0.2 | $|211\rangle$, $|322\rangle$ | adds $|422\rangle$, $|544\rangle$, $|644\rangle$ (Fig. 4, left panel: these five and nothing beyond $n = 6$) | dozens of states, $|n66\rangle$ tower to $n = 32$ and beyond (Fig. 4, middle and right) |
| 0.4 | $|322\rangle$ only ($|211\rangle$ is not superradiant at $\tilde a = 0.95$, $\alpha = 0.4$) | about 20 states: $|422\rangle$, $|433\rangle$, $|533\rangle$, $|544\rangle$, $|633\rangle$, $|644\rangle$, $n = 7$-9 and $|n88\rangle$ to $n = 18$ | dozens, towers to $n = 45$ |

Fig. 25 counts the simultaneously populated states, $N_\mathrm{active}(t)$: at $f_a = 10^{18}$ GeV it is 1
at $\alpha = 0.4$, 0.6 and 0.8 apart from brief spikes to 2-3; at $10^{16}$ GeV it is 2-3 for most of the
evolution at $\alpha = 0.4$ and rises to 10-20 only after about $10^7$ yr; at $10^{14}$ and $10^{12}$ GeV it
reaches 30-60. $\alpha = 0.2$ is not in Fig. 25.

How many of these matter for the black hole rather than merely exceeding $10^{-20}$: the paper says so in
words, not in a table. Sec. III C 3: "the late-time equilibrium is still typically shared among only four or
five dominant modes" for $\alpha \simeq 0.2$-0.8; "eight or nine simultaneously populated modes" around
$\alpha \gtrsim 0.8$; for $\alpha \lesssim 0.15$ the two-level $|211\rangle$, $|322\rangle$ system holds.
Fig. 4 (left, $\alpha = 0.2$, $f_a = 10^{16}$ GeV) shows the hierarchy: $|211\rangle$ and $|322\rangle$ near
$10^{-5}$ of the hole's mass, $|544\rangle$ and $|422\rangle$ at $10^{-6}$-$10^{-7}$, $|644\rangle$ far
below. Witte and Mummery (their ref. 16, quoted in Sec. I) found $n \le 5$ sufficient for the $m = 1, 2$
spin-down, and the higher-$n$ states "barely feed back" on the dominant modes. The relative importance per
state is in the paper's data release (refs. 65, 66), which I have not used.

**So, at this $\alpha$ and $f_a$:** the paper's evidence points to two dominant modes ($|211\rangle$ and
$|322\rangle$) and at most a few weak ones ($|422\rangle$, $|433\rangle$, $|544\rangle$ appear at $f_a \le
10^{16}$ GeV; Fig. 12 shows $|433\rangle$ being kept from growing by GW transitions at $f_a = 10^{18}$ GeV,
$\alpha = 0.4$). For $f_a$ two or more decades lower, or $\alpha \gtrsim 0.6$, the picture is four to nine
dominant states plus a large weak population. Two caveats: their spin is 0.95 against our 0.9, and $\alpha =
0.3$ lies between their columns, so this reading is an interpolation.

### (b) Strong and weak modes: what it does to the cost

The component count is cubic in $K$ only if every bound mode is treated as fully nonlinear. Split the modes
into $S$ strong ones and $W$ weak ones (my notation), where a weak mode's amplitude is small enough that the
field equations may be expanded to first order in it. Then the retained components are all odd $|q|_1 \le 3$
combinations of the strong phases (as now), plus, for each weak mode $w$, the components with exactly one unit
of its phase and zero or two units of strong phases, $q = e_w + p$ with $|p|_1 \in \{0, 2\}$ over the strong
phases. That is $1 + 2S^2$ components per weak mode:

| | strong components | per weak mode | $W = 5$ | $W = 10$ | $W = 20$ | $W = 50$ |
|---|---|---|---|---|---|---|
| $S = 2$ | 8 | 9 | 53 | 98 | 188 | 458 |
| $S = 3$ | 22 | 19 | 117 | 212 | 402 | 972 |
| $S = 4$ | 48 | 33 | 213 | 378 | 708 | 1698 |

(`scaling-work/counts.jl`.) Linear in $W$ instead of cubic in $S + W$: $K = 12$ fully nonlinear would be 1168
components.

The structure matters more than the count. At first order in the weak amplitudes:

* the weak sector does not change the strong solution (the change is of second order), so the strong problem
  is solved first, exactly as now with $K = S$;
* each weak mode is then a separate linear problem: a quasi-bound state with its $2S^2$ satellite components
  on the background of the strong cloud, an eigenvalue problem for its complex frequency (growth rate) with
  $1 + 2S^2$ channels. The weak modes do not couple to one another at this order, so they are independent and
  can run in parallel;
* the phase average for a weak problem is over the strong phases only: the weak phase enters linearly, so its
  average is exact by picking the coefficient. A weak mode adds no dimension to the quadrature.

The cost is therefore (one strong solve at $K = S$) plus $W$ small linear solves, linear in $W$. This is the
asymmetric primary/secondary construction of Omiya et al. (`docs/papers/2211.01949v1.pdf`, Sec. III), which
the old branch implemented under `src/adiabatic/primary_secondary/` and used as "the small-secondary limit of
the joint cosine equations" (`HANDOFF.md`, October 30): a nonlinear $|211\rangle$ plus the coupled
$|322\rangle/(2,-1)$ response, with the $(-1,2)$ component omitted because it is second order in the secondary.

**What is lost.** Processes with two weak legs, which are second order in the weak amplitude: for example
$|322\rangle\times|322\rangle\to|211\rangle\times\infty$ (component $(-1,2)$) when $|322\rangle$ is weak, and the
back-reaction of a weak mode on the strong profiles (frequency shifts of order $|A_w|^2$). The $K = 2$
component table below quantifies this: on the plateau, where $A_{322}/A_{211} \approx 10^{-6}$, the $(-1,2)$
component's share of the rates is $10^{-7}$; in the `nonlinear` state ($A_{322}$ the larger cloud) it carries
2% of the $|211\rangle$ rate and is not negligible.

**When a weak mode must be promoted.** When terms of second order in its amplitude reach the accuracy wanted in
the rates. A usable criterion is the ratio of the flux in its two-weak-leg components to its one-leg flux,
which can be estimated from the linear solution (the second-order source is $h''$ times the square of the weak
field, computable without solving the second-order problem). Promotion changes the strong set from $S$ to
$S + 1$, which costs a factor of 15-25 per added strong mode in the present code (Section 2).

**Does the action balance close?** Yes, at the order kept. Each phase has an action $N_i = \sum_q q_i n_q$
(`phase_actions`). For a weak mode, $N_w$ and the fluxes of its components are all of order $|A_w|^2$, so
$d\ln A_w/dt$ is the growth rate of the linear problem: superradiant gain minus the losses of its satellite
components, functions of the strong amplitudes only. For a strong mode $i$, the weak components carry $q_i \ne
0$ (for example $(2,-1)$ in $|211\rangle\times|211\rangle\to|322\rangle\times\mathrm{BH}$, with $|322\rangle$
weak) and drain $N_i$ at a rate of order $|A_w|^2$; these fluxes come from the linear weak solution, so the
strong balance $\partial_A N\,\dot A = -L_N$ gets a correction $-\sum_w\sum_{q\ni w} q_i\,\Delta F_q$ without
any second-order solve. The strong profiles' own change at order $|A_w|^2$ (a frequency shift) is the part that
needs a second-order solve: one back-substitution with the strong Jacobian per weak mode, and it is the
quantity the paper neglects (Sec. III A, item iii).

**How the paper's rate system relates to this limit.** Their rates are first order in the quartic coupling
with bound-state wavefunctions as the zeroth-order fields: in our language, the fluxes of the $|q|_1 = 3$
components computed with $h(\Phi) \approx \Phi^3/6$ and unperturbed fundamentals, divided by the product of
occupations. That is the small-amplitude limit of the components in our strong/weak scheme. A decisive
consistency check: compute, at small amplitudes, the horizon flux of our $(2,-1)$ component and compare with
$\gamma^{322\times\mathrm{BH}}_{211\times211}$ (Eq. 22) at the same $\alpha$ and spin, after converting their
normalisation (Eq. 19, non-relativistic wavefunction normalisation, energy converted at $\mu$ rather than
$\omega$). The paper prints no numerical rate tables in the text (Fig. 9 shows coefficients of other channels
as a function of $n$; Table II lists super-resonant channels without values); the numbers are in their public
data (refs. 65, 66). This check needs that data and was not done. Two differences to expect: they include an
$\ell = 2$ part of the $m = 0$ fourth leg (Sec. III A), while our default keeps one polar harmonic per
component ($\ell = 0$ for $(2,-1)$; `n_angular = 2` would add $\ell = 2$); and they neglect frequency shifts,
which we include.

## 5. Item 3: which retained components matter

For each converged state, `scaling-work/comp_table.jl` lists every component's peak $|R|$ (relative to the
largest fundamental), its net action flux $\Delta F_q = F_q(\text{outer}) - F_q(\text{inner})$, and its
contribution to each mode's rate $d\ln A_i/dt$. The contributions are exact shares: the action balance
$[\partial_A N]\,\dot A = -L_N$ with $L_N = \sum_q q\,\Delta F_q$ is linear in the fluxes, so the shares add up
to the rate (checked against `amplitude_rates` to all printed digits). Then the state is re-solved with
reduced component sets and the rates compared; that second test also includes each component's feedback on
the others through the nonlinear source, which the shares do not. Fixed background, action balance, default
truncation and quadrature. $K = 2$ on the login node (single thread; not a timing), $K = 3$ in job 27824279.

### $K = 2$

Shares of $d\ln A_i/dt$, each relative to that mode's own rate (signed; they add up to 1), for the two-cloud `nonlinear` state
($A_{211} = 1.47\times10^{-4}$, $A_{322} = 4.1\times10^{-6}$, peak field 0.0034):

| $q$ | $m$ | $\omega_q/\alpha$ | peak $|R|$ | share of $|211\rangle$ rate | share of $|322\rangle$ rate |
|---|---|---|---|---|---|
| (1,0) | 1 | 0.988 bound | 0.42 | 1.02 | $-2\times10^{-3}$ |
| (0,1) | 2 | 0.995 bound | 1 | $-5\times10^{-6}$ | 1.07 |
| (2,−1) | 0 | 0.981 bound | 2.6e-5 | $1.5\times10^{-3}$ | $1.2\times10^{-3}$ |
| (−1,2) | 3 | 1.002 propagating | 9.2e-6 | $-0.022$ | $-0.072$ |
| (3,0), (2,1), (1,2), (0,3) | 3-6 | 2.96-2.98 propagating | 2e-9 - 4e-8 | below $3\times10^{-15}$ | below $10^{-16}$ |

The $(2,-1)$ component is the fourth leg of $|211\rangle\times|211\rangle\to|322\rangle\times\mathrm{BH}$
(bound, $m = 0$, absorbed by the hole) and $(-1,2)$ that of $|322\rangle\times|322\rangle\to|211\rangle\times
\infty$ (propagating, $m = 3$). The four components with $\sum q = 3$ oscillate near $3\mu$; they are the
"$3 \to 1$" processes that the paper also finds negligible.

Re-solving without them (keeping the four components with $\sum q = 1$, which I call the two-to-two set):

| state | peak field | change of $d\ln A_{211}/dt$ | change of $d\ln A_{322}/dt$ |
|---|---|---|---|
| plateau ($A_{322}/A_{211} = 10^{-6}$) | 0.022 | $3.0\times10^{-3}$ | $4.6\times10^{-6}$ |
| nonlinear | 0.0034 | $2.8\times10^{-9}$ | $9.5\times10^{-9}$ |
| nonlinear, amplitudes × 10 | 0.033 | $1.2\times10^{-5}$ | $1.2\times10^{-5}$ |

The $3\times10^{-3}$ on the plateau is not a flux: it is the frequency shift of $|211\rangle$ caused by the
$(3,0)$ harmonic ($\Delta\omega = 4\times10^{-12}$), amplified because the plateau sits at the superradiance
threshold, where the net rate is proportional to $m\Omega_H - \omega$. In an evolution it moves the plateau spin
by an amount of order $\Delta\omega$, about $10^{-11}$, so it is invisible in the trajectory; at fixed background
the relative rate change is the wrong measure there.

### $K = 3$ ($|211\rangle, |322\rangle, |311\rangle$, nonlinear state)

Of the 22 components, the ten with $\sum q = 3$ again contribute below $10^{-12}$ of every rate. The
largest shares, each relative to the mode's own rate:

| $q$ | $|211\rangle$ | $|322\rangle$ | $|311\rangle$ |
|---|---|---|---|
| (−1,2,0) | $2.2\times10^{-4}$ | $-0.049$ | $7\times10^{-9}$ |
| (2,−1,0) | $-1.5\times10^{-5}$ | $8.5\times10^{-4}$ | $-1\times10^{-9}$ |
| (−1,1,1) | $1\times10^{-6}$ | $-1.2\times10^{-4}$ | $-7\times10^{-6}$ |
| (1,−1,1) | $-3.5\times10^{-7}$ | $3.9\times10^{-5}$ | $-2.4\times10^{-6}$ |
| (0,−1,2) | $-1\times10^{-12}$ | $1.4\times10^{-5}$ | $-1.7\times10^{-6}$ |
| (−1,0,2) | $1\times10^{-8}$ | $3\times10^{-10}$ | $-1.6\times10^{-7}$ |
| (2,0,−1), (1,1,−1), (0,2,−1) | below $10^{-11}$ | below $2\times10^{-12}$ | below $2\times10^{-11}$ |

Re-solved rates: the two-to-two set (12 of 22 components) changes the rates by at most $7.5\times10^{-9}$. A
set chosen by rate share needs care in how the share is normalised: my first rule normalised every share by the
largest rate of the state ($|211\rangle$'s), kept only the three fundamentals at the $10^{-3}$ level, and changed
the $|322\rangle$ rate by 4.9% (it dropped $(-1,2,0)$, which carries 4.9% of that mode's own rate). With the same
global normalisation at $10^{-6}$ it kept six components and changed the rates by $5\times10^{-5}$ at most.
Normalising each share by the mode's own rate, the table above gives, for $10^{-4}$ accuracy, the three
fundamentals plus $(-1,2,0)$, $(2,-1,0)$, $(-1,1,1)$ (6 of 22), and for $10^{-6}$ three more (9 of 22).

### Selection rule

* **Drop the components with $\sum q = 3$** (frequency near $3\mu$). Evidence: their shares are below
  $10^{-12}$ at $K = 2$ and 3, and the re-solved rates change by $10^{-8}$ (two-cloud states) to $10^{-5}$ (field
  amplitude 0.03); the one larger number, on the plateau, is a frequency shift that does not affect the
  trajectory. This keeps $K^2 + 3\binom{K}{3}$ components (4, 12, 28, 96, 460 for $K$ = 2, 3, 4, 6, 10): a factor
  of about 1.5 at large $K$, still cubic. It also halves the slot count $S$ (11, 55, 214, 1817, 33851).
* **Keep a component only if its rate share, relative to each mode's own rate, exceeds the target accuracy.**
  At $K = 3$ this keeps 6 of 22 components for $10^{-4}$. The shares need a solve with the components present,
  so in practice the rule is applied once per parameter region (or re-checked occasionally along an evolution),
  not per step. The kept components are the fundamentals plus the $2\to2$ fourth legs between pairs of
  populated modes, which is the set of processes in the paper's rate system. Their number is set by how many
  pairs of modes interact appreciably, not by $K^3$.
* An amplitude-based rule (keep $q$ if $\prod_i A_i^{|q_i|}$ is large) is not enough on its own: the
  propagating $(-1,2)$ component has a peak $|R|$ of $10^{-5}$ of the fundamentals yet carries 5-7% of the
  $|322\rangle$ rate, because what matters is whether a component carries flux to infinity or into the horizon,
  not its size.

What this establishes: the importance ranking and the effect of dropping components on the rates at three
$K = 2$ states and one $K = 3$ state with small field amplitudes. It does not establish the ranking near
resonances (a fourth leg whose frequency lies close to a quasi-bound frequency, which the paper finds can
enhance a rate by orders of magnitude), at larger $\alpha$, or for the $|433\rangle$ variant (its nonlinear state
was not reached; see Section 2).

## 6. Item 4: the linear algebra at large channel count

### (i) The dense block LU against $K$

Per mesh node the factorisation eliminates a dense block of side $n = 4C$ and stores about $6n^2$ numbers
(Jacobian workspace and factorisation panels, Section 2). Measured on 690 nodes, 8 threads (the factorisation
itself runs on one thread), and extrapolated with $2.3\times10^{-10}$ s per $n^3$ per node and
$6n^2\times8$ bytes per node:

| $K$ | $C$ | $n = 4C$ | factorisation, 690 nodes | memory, 690 nodes | |
|---|---|---|---|---|---|
| 2 | 8 | 32 | 0.016 s | 34 MB | measured |
| 3 | 22 | 88 | 0.15 s | 0.26 GB | measured |
| 4 | 48 | 192 | 1.1 s | 1.2 GB | measured |
| 6 | 152 | 608 | 36 s | 12 GB | extrapolated |
| 10 | 680 | 2720 | 3200 s | 245 GB | extrapolated |

The interval blocks of the Jacobian cost two dense $n\times n$ products each, another $4n^3$ operations per
node, so assembly has the same $K^9$ growth on top of its phase-point work. With the meshes that larger mode
sets need (4256 nodes when $|433\rangle$ is included) multiply by six. Dense blocks are not viable beyond
$K \approx 5$ with the full truncation.

### (ii) Newton–Krylov

**What a Krylov method is.** To solve $Jx = b$ without factorising $J$, GMRES builds the solution from the
vectors $b, Jb, J^2b, \dots$ (the Krylov subspace): after $m$ steps it picks the combination of the first $m$
vectors that minimises the residual $|b - Jx|$. It only needs products $Jv$, never the matrix itself, and it
converges in few steps when the eigenvalues of $J$ are clustered.

**What a preconditioner is.** Our $J$ has eigenvalues spread over many orders of magnitude (the radial
operator), and plain GMRES makes almost no progress: no reduction below $10^{-3}$ in 150 steps (measured,
below). A preconditioner is an approximation $P \approx J$ that is cheap to invert; GMRES is applied to
$JP^{-1}$, whose eigenvalues cluster near 1 if $P$ captures the stiff part of $J$. Here the natural choice is
the linear part of the radial operator: without the nonlinear source it couples channels only within one
Fourier component (through the Kerr angular mixing), so it splits into $C$ independent small block-bidiagonal
systems with blocks of side $4 \times$ (channels per component) $= 4$, bordered only by the $2K$ frequency
columns. Factorising those costs of order $C N_\mathrm{nodes}$ small operations, nothing like $(4C)^3$.

**Matrix-free products.** $Jv$ is the directional derivative of the residual along $v$. The residual code
already accepts ForwardDiff dual numbers, so one residual evaluation at $x + \epsilon v$ with $\epsilon$ a dual
unit returns $Jv$ exactly, without assembling anything. Measured: it agrees with the assembled Jacobian to
$10^{-16}$ and costs about three residual evaluations (8 threads: 0.013-0.019 s against 0.005 s at $K = 2$,
0.27-0.37 s against 0.11-0.13 s at $K = 3$).

**The experiment** (`scaling-work/nk_proto.jl`; jobs 27824280 for $K = 2$, 27824279 for $K = 3$). From a
converged state, both amplitudes are multiplied by 1.01 and one Newton step is solved: exactly with the block LU,
and by GMRES (right-preconditioned, in the Newton scaling) with two preconditioners: *free*, the Jacobian of
the problem without the nonlinear source, and *diag*, the Jacobian with the source's Jacobian kept only within
each component (the $\hat h'(0)$ and $\hat h'(2q)$ terms) and dropped between components. Both are
block-diagonal across components apart from the frequency border; in the prototype they are factorised with
the existing block LU, so the prototype counts iterations, not time. Larger fields come from continuing the
`nonlinear` state's amplitudes upward.

| state | peak field $\Phi_\mathrm{max}$ | GMRES steps to $10^{-6}$ / $10^{-10}$, free | same, diag | no preconditioner |
|---|---|---|---|---|
| $K = 2$ linear | $8\times10^{-6}$ | 1 / 1 | 1 / 1 | not below $10^{-3}$ in 150 |
| $K = 2$ plateau | 0.022 | 2 / 3 | 2 / 2 | not below $10^{-3}$ in 150 |
| $K = 2$ nonlinear | 0.0034 | 2 / 3 | 1 / 3 | not below $10^{-3}$ in 150 |
| $K = 2$ nonlinear × 3 | 0.010 | 2 / 3 | 2 / 3 | |
| $K = 2$ nonlinear × 10 | 0.033 | 3 / 5 | 3 / 4 | |
| $K = 2$ nonlinear × 30 | 0.085 | 4 / 8 | 4 / 6 | |
| $K = 2$ nonlinear × 100 | 0.18 | 9 / 12 | 6 / 9 | |
| $K = 3$ linear | $7\times10^{-5}$ | 1 / 2 | 1 / 2 | not below $10^{-3}$ in 150 |
| $K = 3$ nonlinear | 0.0034 | 2 / 3 | 1 / 3 | not below $10^{-3}$ in 150 |

The GMRES step agrees with the exact Newton step to $10^{-10}$ or better whenever the relative residual
reaches $10^{-10}$. The preconditioner stays effective at the saturated state (the plateau, 3 steps) and the
step count grows slowly with the field (3 steps at $\Phi_\mathrm{max} = 0.01$, 12 at 0.18): the nonlinear coupling between
components, which the preconditioner leaves out, is of relative size $\Phi^2$ and appears only through the
nonlinear source. Keeping the within-component source terms (*diag*) saves a step or three at the larger
fields. At the field amplitudes that occur in the clouds ($\Phi \lesssim 0.03$) a Newton step costs about
3-5 matrix-free products, that is 10-15 residual evaluations, with no dense matrix at all; memory is a few
dozen vectors of the length of the unknown vector.

What this establishes: preconditioned GMRES solves the Newton system of real $K = 2$ and $K = 3$ states in a
few iterations, and the iteration count does not grow from $K = 2$ to $K = 3$. What it does not establish: the
time of a dedicated component-diagonal factorisation (not written; the operation count is small, but the
frequency border and the boundary rows have to be handled, for example by a Schur complement of size $2K$),
behaviour with many nearly degenerate frequencies or near resonances, and iteration counts at $K \ge 4$.

### (iii) Sensitivities

The implicit-function-theorem quantities need $J^{-1}$ applied to a few right-hand sides at the converged point:
`amplitude_response` needs $K$ (one per amplitude) and the spindown balance two more (mass and spin of the
hole). With a factorisation these cost one back-substitution each. With GMRES each is a separate solve of
3-12 iterations at the measured fields (or one block-GMRES solve for all), so the sensitivities cost about
$(K + 2) \times 3$-12 matrix-free products, comparable to or larger than the Newton solve itself; they need a
tight tolerance because the rates come from differences of nearly equal fluxes, which the measured convergence
to $10^{-10}$ in a dozen steps supports. Cost remains proportional to the residual cost, so it scales with the
quadrature and the component count, not with $(4C)^3$.

## 7. Item 5: synthesis

### Cost per right-hand-side evaluation

Time of one warm solve (the closest single number to one right-hand-side evaluation), 8 threads, on a
690-node mesh, and memory (Jacobian blocks, factorisation and phase tables). **M** = measured (Section 2);
all other entries are extrapolated with `scaling-work/cost_table.jl`:

* warm solve $\approx$ 5 residuals + 2 (Jacobian + factorisation) (checked: 1.41 s predicted from the measured
  $K = 3$ pieces against 1.47 s measured);
* residual $= 1.5\times10^{-9}\,\mathrm{s}\times C P N$, source Jacobian $= 7.6\times10^{-10}\,\mathrm{s}\times
  S P N$ plus the dense interval products, factorisation $= 2.3\times10^{-10}\,\mathrm{s}\times(4C)^3N$
  (constants from the measured $K = 4$ pieces); this omits a fixed cost of about 0.07 s per solve that dominates
  at $K = 2$ (boundary blocks, far-field series), so it underestimates $K = 2$ (0.03 s against 0.10 s measured);
* lattice: $P$ = 101, 701 (measured, Section 3), 2398 (constructed exact lattice), then $2.5\,B(6,6) = 22500$
  and $4\,B(10,6) = 5.4\times10^5$ (extrapolated from the trend of found-to-bound ratios 1.0, 1.2, 1.9, 2.4);
* selection: the two-to-two set ($\sum q = 1$; $C$ = 4, 12, 28, 96, 460, $S$ = 11, 55, 214, 1817, 33851); a
  rate-share selection would keep fewer (Section 5) but its count at large $K$ is not known;
* Krylov: per solve 2 Newton steps of 4 matrix-free products plus one residual, and $K + 2$ sensitivity solves of
  4 products, each product 3 residuals, so $(26 + 12(K+2))$ residuals; no Jacobian table, no dense blocks; the
  component-diagonal preconditioner taken as free (not measured).

| $K$ | present code | lattice quadrature | lattice + two-to-two selection | + Newton–Krylov |
|---|---|---|---|---|
| 2 | 0.10 s, 34 MB **M** | 0.02 s | 0.005 s, 9 MB | 0.03 s, 5 MB |
| 3 | 1.5 s, 0.27 GB **M** | 0.4 s | 0.12 s, 77 MB | 0.75 s, 16 MB |
| 4 | 47 s, 1.8 GB (pieces **M**) | 4 s, 1.2 GB | 1.4 s, 0.4 GB | 7 s, 40 MB |
| 6 | 19 h, 940 GB | 170 s, 13 GB | 74 s, 6 GB | 270 s, 0.2 GB |
| 10 | impossible ($16^{10}$ points) | 11 h, 710 GB | 6 h, 410 GB | 12 h, 9 GB |

Strong/weak scheme (strong modes on the lattice with the two-to-two set and the direct solver; each weak mode a
linear problem with $1 + 2S^2$ components on the same lattice, one Jacobian, factorisation and 3 residuals;
weak modes counted serially although they are independent):

| | strong solve | per weak mode | $W = 5$ | $W = 10$ | $W = 50$ |
|---|---|---|---|---|---|
| $S = 2$ | 0.005 s (0.10 s **M** with the present code) | 0.012 s | 0.07 s | 0.13 s | 0.6 s |
| $S = 3$ | 0.12 s (1.5 s **M** with the present code) | 0.14 s | 0.8 s | 1.5 s | 7 s |

All extrapolated except where marked. Multiply every time by about 6 when a near-threshold combination
frequency forces a 4256-node mesh (as $|433\rangle$ does at $\alpha = 0.3$); multiply the per-solve time by the
number of right-hand-side evaluations of an evolution (2233 for the production window) for the wall time.

What the table says:

* The tensor grid is the first wall: at $K = 6$ its tables alone are hundreds of GB. The lattice removes it.
* Dense block LU is the second wall: with the full truncation it reaches tens to hundreds of GB at $K = 6$-10.
  Krylov removes it, at the price of more residual evaluations, so it does not pay at $K \le 4$.
* Even with all three changes, $K = 10$ fully nonlinear costs hours per right-hand-side evaluation, because the
  residual itself scales like $C \times P \sim K^3 \times K^6$.
* The strong/weak scheme keeps every piece small: the cost of $S = 2$-3 strong modes plus tens of weak ones is
  comparable to today's $K = 2$-3 solve.

### Is "order 10 fully nonlinear modes" the right target?

No. The paper's evidence (Section 4) at the parameters of the production runs ($\alpha = 0.2$-0.4,
$f_a = 5.7\times10^{17}$ GeV) is two dominant states and at most a few weak ones; even in the richest regions
they report four or five dominant modes at late times for $\alpha = 0.2$-0.8 (eight or nine around $\alpha \gtrsim
0.8$), plus a large population of weak states that "barely feed back". The numerics point the same way: the
fully nonlinear cost grows like $K^9$ at best, while a weak mode costs a small linear solve. The design to aim
at is **2-4 strong modes plus tens of weak ones, chosen per $\alpha$ (and $f_a$) from the paper's Fig. 2, with
promotion of a weak mode to strong when its second-order terms reach the target accuracy**. Fully nonlinear
$K = 3$-4 is affordable with the lattice and the component selection; beyond that, weak modes.

### Recommended sequence

1. **Lattice quadrature in the solver** (behind `phase_grid`, already in the library). Run the harness and one
   evolution segment at $K = 2$ with a 101-point lattice against the tensor grid. Establishes that the lattice
   works in full solves and evolutions, not only at one node; at $K = 3$ it cuts the cost about 3 times.
2. **Small-amplitude consistency check against the paper's rates**: the horizon flux of $(2,-1)$ at small
   amplitudes against $\gamma^{322\times\mathrm{BH}}_{211\times211}$ from their data release, with
   `n_angular = 2` to include their $\ell = 2$ leg. Establishes that the two calculations agree in the limit where
   both apply, which is what makes the paper's state lists usable for our mode choice.
3. **Weak modes as linear problems on the strong background** (the old branch's primary/secondary construction,
   generalised to $S$ strong modes): the satellite components $e_w + p$, the linear eigenvalue problem per weak
   mode, the weak-sector fluxes in the strong action balance. Check: with $|322\rangle$ weak on the $|211\rangle$
   plateau, the weak-mode growth rate equals the full $K = 2$ solution's $d\ln A_{322}/dt$ to the size of the
   dropped second-order terms ($10^{-7}$ by Section 5).
4. **Component selection** for the strong sector: drop $\sum q = 3$ (cheap, established), then a rate-share rule
   re-checked along the evolution.
5. **A $K = 3$ strong evolution** ($|211\rangle, |322\rangle$ plus the next mode the paper lists at the chosen
   $\alpha$), with the promotion criterion monitored.
6. **Newton–Krylov** only if the strong sector has to grow beyond about 4 modes, or memory becomes the limit.

### Risks

* **Aliasing for non-polynomial potentials.** The number of lattice points is set by the decay of the harmonics
  of $h(\Phi)$, which is fast at the field amplitudes of the clouds ($\Phi \lesssim 0.03$) but slow for the chiral
  potential at $\Phi \gtrsim 2$. A lattice should be checked against a finer one at the largest field of a run.
* **Near-threshold components.** Combination frequencies close to the mass threshold (such as
  $2\omega_{433} - \omega_{322}$) need very long domains or a different outer treatment, and they multiply the
  mesh. This grows with the number of modes and is independent of the quadrature; it also governs which states
  can be added cheaply.
* **Resonances.** The paper finds fourth legs that land almost on a quasi-bound frequency and enhance a rate by
  orders of magnitude (super-resonances, Table II). A component-selection rule based on present shares can miss a
  component that becomes resonant as the hole spins down; such a component also makes the Newton system
  ill-conditioned (a nearly singular linear operator for that component), which would raise the Krylov iteration
  count and the condition number of the direct solve alike.
* **Nearly degenerate frequencies.** Modes such as $|311\rangle$ and $|322\rangle$ ($\omega/\alpha$ = 0.9946 and
  0.9948 here) give components with nearly equal frequencies and small frequency differences; the action balance
  matrix $\partial_A N$ and the frequency border of the Newton system can become ill-conditioned. Not seen in the
  $K = 3$ solves here, which converged in 3-4 iterations.
* **Energy versus action balance.** The energy/angular-momentum balance (`balance=:energy`) gives two equations
  and closes only for two modes; the action balance (the default) has one equation per phase and handles any $K$,
  including weak modes.
* **Polar truncation.** One polar harmonic per component (the default) leaves out the $\ell = 2$ part of $m = 0$
  fourth legs that the paper finds can dominate through resonances at larger $\alpha$.


## Files

In the repository (branch `speedup/scaling`): `src/field.jl` (`PhaseGrid(points, mode_set)`, `rank1_lattice`),
`src/radial.jl` (`phase_grid` keyword), `test/runtests.jl` ("Rank-1 lattice phase rule").

Scratch, `axion-sr-sim-speedup/scaling-work/`: `k_timing.jl` and `run_k.sh` (item 1), `lattice_proto.jl`,
`min_lattice.jl` (item 2), `comp_table.jl` (item 3), `nk_proto.jl` (item 4), `counts.jl`, `slot_counts.jl`,
`cost_table.jl` (counts and Section 7), `logs/` (Slurm output), `states/` (converged $K = 2, 3$ unknown vectors,
Julia serialization).
