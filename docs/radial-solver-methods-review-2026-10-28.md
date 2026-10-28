**Radial solver methods review — October 28, 2026**

**Follow-up result:** the isolated adaptive collocation experiment converged
directly at the target where the recorded shooting attempts stalled. Both
initial guesses required three Newton corrections, without amplitude
continuation. A tighter 881-point solve also passed independent residual and
short-integration checks. See the
[experiment and numerical results](../runs/adiabatic/radial_bvp/README.md).
This establishes a convergence improvement for this case, not a measured
warmed speedup. The implementation includes two tested experiment-local
corrections to the chosen solver library.

The review below motivated that experiment using the saved failures, code
inspection, and numerical literature. A secant or tangent continuation
predictor and multiple shooting remain useful alternatives, but neither has
been benchmarked here as a competing improvement.

This review used the recent local solver conversations, the
[target investigation](multiphase-target-solve-investigation-2026-10-28.md), the
[derivative verification](multiphase-ad-verification-2026-10-27.md), and the live
implementation. The conversations also record a request for twelve more
driven-seed corrections; this review does not assume their outcome. Live code
has concurrent edits, so any comparison should use a frozen source snapshot.
This research added no solver implementation, dependencies, or cluster jobs.

**What the evidence establishes**

The zero-driven target guess at join 60 gives a scaled residual of 56.5614.
The first full Newton step gives 6118.8183; a 1/32 step is the first to reduce
it. Subsequent corrections require still smaller steps. Six corrections reduce
the residual by only 4.75%. Equilibrated and 256-bit solutions of the *saved*
linear system agree with the original correction to about 9.5e-15 relative.
Consequently, inaccurate solution of that matrix is not the main obstacle.
These checks do not validate its entries at the target or prove spatial accuracy.

Moving the join from 5 to 60 reduces the original mismatch from 1.906e9 to
56.56. Adding driven boundary coefficients changes the Newton trajectory but
does not converge in the first six corrections. Both results make the
formulation and initialization worth investigating; neither proves that the
physical boundary-value problem has no solution.

The code explains the exposure to this difficulty:

- [Matching](../src/adiabatic/multiphase/profiles.jl) integrates from both
  boundaries to one join. Interior profiles are consequences of boundary
  guesses, not independently corrected unknowns.
- [Integration](../src/adiabatic/core/integration.jl) uses fixed-step RK4.
- [Newton](../src/adiabatic/core/newton.jl) recomputes a full matching Jacobian,
  takes its Newton direction, and halves the step until the norm decreases.
- The continuation predictor rescales each coefficient by a power of the
  amplitude ratio and leaves both frequencies unchanged. It does not extrapolate
  the measured branch or use the available implicit response.

My interpretation is that long shooting intervals are a strong candidate for
the small useful Newton steps. This remains a hypothesis until compared with
short intervals or a full-profile formulation. A well-conditioned BVP can have
badly conditioned initial-value integrations. Kramer and Mattheij explicitly
analyze growing modes, small Newton convergence domains in shooting, and the
benefits of shorter intervals and local collocation.
[Their paper](https://pure.tue.nl/ws/portalfiles/portal/2231026/353772.pdf)
is particularly close to the numerical issue we need to distinguish.

**Recommended order**

| Candidate | Expected benefit for this code | Priority and limitation |
|---|---|---|
| Secant/tangent continuation prediction | Predict frequencies and deformation from the actual branch; reduce corrections per step | Smallest useful experiment; still uses sensitive shooting |
| Adaptive full-profile collocation | Correct local profiles directly, exploit radial sparsity, concentrate resolution where needed | Leading replacement candidate; requires a BVP adapter and accuracy checks |
| Multiple shooting | Limit amplification within each integration interval | Useful intermediate control; may cost more per iteration |
| Better exterior boundary map | Improve outgoing/decaying boundary accuracy without evaluating the nonlinear source out to a very distant endpoint | Separate correctness task; exterior approximation must be justified |
| Multiple-domain spectral discretization | Potentially reduce expensive source evaluations further | Second-stage performance option after a robust BVP baseline |

**The smaller gap: continuation prediction**

For the real packed matching equations F(x,A)=0, a branch tangent at a
converged point satisfies

\[
 F_x x_A=-F_A,\qquad
 x_{\rm pred}=x+x_A\,\Delta A.
\]

This predicts the frequencies and all free boundary coefficients together.
The necessary mathematics is already implemented in
[`implicit_branch_response`](../src/adiabatic/core/sensitivity.jl). The current
helper rebuilds the derivatives; reusing a Jacobian/factorization would require
an explicit implementation change. In particular, Newton's final iteration
may have used a Jacobian at the previous point, so it should not be described
as an exact tangent at the accepted solution without checking or rebuilding it.

An even cheaper first experiment is a secant through two accepted solutions.
Keep a fixed packing and coefficient scale, and compare prediction in A2 and
log(A2). For the latter, the exact tangent is x_logA2 = A2*x_A2. Small changes
in weak-dependent coefficients can make a secant noisy, so inspect prediction
errors rather than assuming improvement. The existing power-law predictor can
remain a fallback. There is no need to compute E/J derivatives merely to
predict a profile.

This is also a detail of the original physics calculation: Omiya et al.
describe extrapolating previous solutions before iterative linearized
corrections in Appendix A. Their radial solve uses matching; their mention
of orthogonal collocation there concerns angular integration. Thus that paper
does not benchmark a radial collocation replacement.
[Omiya et al., 2201.04382, Appendix A](https://arxiv.org/html/2201.04382)

If the branch develops a turning point in A2, use pseudo-arclength continuation
to solve for the amplitude along with the profile. This addresses a bad branch
parameterization, which is different from sensitivity of radial shooting.
We have not established a turning point. BifurcationKit documents both secant
and bordered tangent methods, with an augmented Newton corrector for
pseudo-arclength continuation.
[BifurcationKit PALC](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/PALC/)

**The main candidate: solve the whole radial BVP**

For our 14 retained complex radial modes, set

\[
 y=(\Re R,\Im R,\Re Q,\Im Q)\in\mathbb R^{56},\qquad
 \nu=(\Re\omega_1,\Re\omega_2,\Im\omega_1,\Im\omega_2).
\]

Solve y'=f(r_*,y,nu) on the finite radial domain, with the current horizon
and outer boundary functions. The unknowns include y at mesh points and the
four real frequency parameters. The nonlinear source couples modes at the
same radius; it does not couple distant radial points. A local collocation
scheme therefore produces a sparse matrix with nearby radial blocks and a
small set of global parameter columns.

The simplest adapter can preserve the existing boundary coefficients as extra
unknown parameters: 28 real horizon coefficients, 24 free real outer
coefficients, and four real frequencies, totaling 56 parameters. Impose the
full 56-real-component state at each endpoint through the existing boundary
builders, inserting the two prescribed positive real outer amplitudes. There
are then 112 real boundary equations for 56 first-order equations and 56
parameters. With N mesh nodes, the count is
56(N-1)+112 = 56N+56 equations and unknowns. This avoids deriving a new boundary
condition for the first experiment.

A later reduction can eliminate the boundary coefficients using the allowed
boundary subspaces. It leaves 28 constraints at each end plus four real
amplitude/phase constraints, matching 56 differential equations and four
frequency parameters. Do not impose R and Q arbitrarily while also fixing
all boundary coefficients: that would overconstrain the eigenproblem.

This also changes how an initial guess can be used. We can provide the
*entire* primary, secondary, and forced radial profiles, including interior
values and slopes. The previous driven-seed test supplied their boundary
coefficients and then re-shot through the coupled nonlinear equations. Its
failure is not a test of a full-profile seed in collocation. Boundary
conventions must still be made consistent.

The existing
[`local_source_jacobian`](../src/adiabatic/multiphase/projection.jl) is directly
useful. Combine its R derivatives with the linear slope/operator terms to
assemble f_y locally. Differentiate frequency-dependent operators and endpoint
data as well. Real packing is necessary because reconstruction uses both
coefficients and their conjugates; treating the system as holomorphic would
be incorrect. SciPy's BVP documentation explicitly makes this distinction.
[SciPy solve_bvp](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_bvp.html)

For a Julia prototype, BoundaryValueDiffEq's MIRK4/MIRK6 are concrete candidates:
the documented methods have defect-control mesh adaptation and configurable
nonlinear solvers. Its two-point solver supports sparse differentiation.
[MIRK documentation](https://docs.sciml.ai/BoundaryValueDiffEq/v5.23/solvers/mirk/)
The documented `tune_parameters=true` interface jointly solves unknown
parameters with the BVP, provided enough boundary constraints are supplied.
[Unknown parameters](https://docs.sciml.ai/BoundaryValueDiffEq/v5.23/tutorials/unknown_parameters/)
These packages are absent from our Project.toml. This is a verified API option,
not an integration already tested against our custom AD specialization.

Use a sparse Jacobian or explicit local blocks; a dense ForwardDiff Jacobian
over every mesh unknown would discard the main structural advantage. Test the
chosen differentiation path explicitly, including its sparsity discovery and
our specialized Dual dispatch.
[BVP differentiation](https://docs.sciml.ai/BoundaryValueDiffEq/v5.23/basics/autodiff/)

**Why the larger system could be faster**

At the recorded step 0.1 over [-100,300], the two branches together take
4,000 RK4 steps and 16,000 source evaluations per residual. Each source
evaluation reconstructs 16*31*31 = 15,376 angular/phase samples: approximately
246 million potential evaluations per residual. These are counts derived
from the settings and code, not a new performance profile.

The saved target logs report roughly 8 seconds per residual and 50–77 seconds
per Jacobian. With nine trial residual evaluations, trial integrations alone
cost about 72 seconds per correction. Reducing rejected trials is therefore
as relevant as accelerating derivatives.

Collocation can evaluate sources on an adaptive mesh, assemble local
derivatives without propagating global boundary sensitivities, and reuse a
converged full profile as the next predictor. Its linear systems are larger,
and factorization/fill-in can become significant. No particular mesh size or
speedup is justified yet. Compare total time to a validated solution, not
only Jacobian time or iteration count.

After convergence, differentiating the discretized BVP gives the same form
G_U U_A = -G_A. A factorization at the accepted solution can serve multiple
amplitude-response right-hand sides. This is useful for evolution, but the
charges, fluxes, their derivatives, and profile corrections must all be included
in an eventual end-to-end timing.

**Multiple shooting: useful, with an implementation trap**

Introduce interior states y_j and enforce
Phi_j(y_j,nu)-y_(j+1)=0 on short intervals, alongside the boundary conditions.
Keep the interface unknowns in the sparse solve. Eliminating every interface
by multiplying transfer matrices across the full domain can recover the
sensitivity we intended to remove. The BVP literature discusses this problem
with conventional condensation.
[Kramer–Mattheij](https://pure.tue.nl/ws/portalfiles/portal/2231026/353772.pdf)

Start with a modest segmentation, then refine where segment sensitivities
grow. Merely dividing RK4 into chunks while immediately passing each endpoint
to the next chunk is still the same shooting calculation. Independent
interface variables and continuity equations are essential. More segments
also increase unknowns and derivative work; speedup is not automatic.

BoundaryValueDiffEq provides `MultipleShooting`, but its versioned
documentation says the grid is coarsened and a final single-shooting problem
is solved even for `grid_coarsening=false`. Before adopting it for this
purpose, inspect and test the pinned implementation's finalization behavior.
We need a retained multiple-shooting solution or a collocation solution, not
an obligatory return to a long unstable IVP.
[MultipleShooting options](https://docs.sciml.ai/BoundaryValueDiffEq/v5.23/solvers/shooting/)

**Boundary accuracy is a separate issue**

The investigation found a 0.283333 relative difference between outgoing
states constructed directly at 300 with order 7 and propagated from 5000
with order 16, holding the far-leading coefficients fixed. That is evidence
for a boundary check, not a measured 28% error in a solved profile. It also
does not establish the difference between the allowed subspaces after
refitting their free coefficients.

A useful comparison constructs the entire outer boundary basis for each
prescription, normalizes field and slope units consistently, and measures
subspace angles or projection residuals. For a linear exterior, propagate
the allowed solution subspace from a sufficiently distant radius and impose
its relation between values and slopes at the nonlinear domain edge. If a
matrix of field values U is invertible there, the relation is Q=Z R with
Z=U' U^{-1}; implement this with linear solves. QR-based subspace propagation
or changing Riccati charts is needed when that representation becomes singular.
Log-derivative propagation is an established multichannel scattering method.
[Johnson, 1973](https://doi.org/10.1016/0021-9991(73)90049-1)

For our problem this is a proposed exterior adapter. It requires checking that
omitted nonlinear sources in the exterior are small, including their integrated
effect, and accounting for the map's dependence on the trial frequencies.
The relation between finite-radius data and the prescribed asymptotic amplitudes
must also be preserved. A boundary map does not by itself solve the nonlinear
interior. For the first solver comparison, hold the existing boundaries fixed;
then vary them separately before interpreting physical differences.

**Other methods and their place**

A scaled trust-region or Levenberg–Marquardt step can change the direction
when a Newton step predicts poorly, whereas the current line search only
shortens one direction. This is worth a bounded comparison, especially on
the existing 56-variable problem. Choose meaningful variable and residual
scales, avoid explicitly forming poorly conditioned normal equations, and
require the original root residual to pass: least-squares stationarity alone
is insufficient. Julia already has these algorithms.
[NonlinearSolve methods](https://docs.sciml.ai/NonlinearSolve/stable/native/solvers/)
It remains an improvement to globalization, not a cure for radial sensitivity.

Multiple-domain Chebyshev or ultraspherical discretizations are plausible
next steps if source evaluations dominate after collocation works. Factor
known boundary oscillations/decay where helpful, and treat compactified
outgoing tails carefully. Olver–Townsend's almost-banded, preconditioned
spectral method is relevant to linearized radial operators; its stability
result is not a convergence theorem for our nonlinear coupled eigenproblem.
[Olver–Townsend](https://arxiv.org/abs/1202.1347)

There is adjacent physics precedent for spectral cloud calculations:
Chia et al. use KADATH for self-interacting complex-scalar clouds. However,
they neglect spin and fix the frequency to the free-field value in their
self-interacting calculations, explicitly identifying that approximation as
a limitation. Their result does not demonstrate a ready-made solver for our
two complex frequencies and dissipative real-field components.
[Chia et al., sections 2.2.3–2.2.4](https://arxiv.org/html/2212.11948)

A Green-function or block-preconditioned iteration could reuse our forced
solvers, but would still need a nonlinear convergence strategy and a bordered
treatment near fundamental eigenmodes, where an unprojected linear inverse
is singular or nearly singular. Likewise, multiplying the nonlinear source
by a homotopy parameter could give a path to the target at fixed amplitudes,
but requires a genuinely solved starting problem and is not guaranteed to
avoid folds. These are secondary options, not reasons to postpone the direct
BVP comparison.

An adaptive or implicit IVP integrator can reduce integration error/cost, but
does not remove amplification inherent in a long shooting map. Faster source
kernels, cached geometry, or parallel interval work may help after profiling.
Higher precision for the saved linear solve has already been tested without
changing its direction. None of these alone addresses the observed poor
finite-step prediction.

**A bounded experiment that would settle the next decision**

1. Check the target Jacobian along the actual saved Newton direction using a
   displacement sweep. Compare F(x+t*d) with F(x)+t*J*d and centered directional
   differences; inspect the expected small-t regime before roundoff. The
   previous full derivative verification was at weak-secondary states. Repeat
   a useful part of the check at smaller radial spacing if the derivative
   comparison exposes a discretization concern.
2. From identical accepted points, compare the existing continuation predictor
   with a secant and, if worthwhile, a tangent. Hold correction tolerances and
   residual weights fixed between methods at each amplitude. Count all
   prediction, Jacobian, and rejected-trial costs. Include larger proposed
   amplitude increments rather than judging only the current tiny increment.
3. Build an isolated full-profile BVP adapter using the existing equations,
   14-mode support, quadrature, finite boundaries, and amplitudes. First recover
   a known weak solution; then test the target using full reference profiles,
   with and without driven fields. Use adaptive collocation as the principal
   candidate; a short-interval multiple-shooting control can isolate whether
   shooting length is the dominant problem. Keep compilation separate from
   warmed timings, and bound both corrections and wall time.
4. Require independent accuracy evidence: off-mesh differential defects,
   endpoint constraints, finite fields, mesh refinement of frequencies and
   profiles, and convergence of the E/J/flux quantities we intend to use.
   Include componentwise scales and an absolute floor so weak modes are not
   hidden by the primary. A collocation tolerance of 1e-5 has different units
   and normalization from the current matching norm; it is not automatically
   an equivalent acceptance test. Report the old shooting residual where
   meaningful, but also use short local reintegrations: requiring a fresh
   unstable long shoot as the sole validator would defeat the comparison.
5. Once solver behavior is understood, vary outgoing boundary construction,
   radial endpoints, and angular/phase resolution independently. Track branch
   identity by overlap with accepted profiles and continuity of frequencies;
   a fast root on another branch is not the intended result.

The decision criterion is a reproducible reduction in total cost at comparable
verified accuracy, together with useful convergence at the target or much
larger accepted continuation steps. Small algebraic residuals alone cannot
establish spatial accuracy; residual control and solution error are related
through BVP conditioning.
[Kierzenka–Shampine on residual and error control](https://jnaiam.org/wp-content/uploads/2024/02/A-BVP-Solver-that-Controls-Residual-and-Error.pdf)

The subsequent target results support moving toward an optional collocation
backend, preserving the independently checked equations and boundary data.
The published methods provided a concrete test of whether a global profile
solve could improve convergence without changing the physical equations.
Library fixes, implementation overhead, and physical boundary-systematics
checks remain explicit parts of the integration work.
