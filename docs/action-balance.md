# The wave-action balance

This note explains the wave-action balance in `AxionCloud.jl`. It is an
alternative to the energy and angular-momentum balance for evolving the cloud
amplitudes. It was added on 2026-11-07 (commit `c1fc133`) and has been the
default since 2026-11-09. Pass `balance=:energy` to reproduce results computed
before then; the validation scripts that compare against the energy balance or
reproduce older logs pass it explicitly.

Code: [`src/charges.jl`](../src/charges.jl) (actions and their fluxes),
[`src/evolution.jl`](../src/evolution.jl) (balances, metric work).
Checks: [`validation/action_balance.jl`](../validation/action_balance.jl), and
the testset "Wave actions: identities and balances" in
[`test/runtests.jl`](../test/runtests.jl).
Background equations: [`docs/multiphase.tex`](multiphase.tex), cited below
by equation label.

## 1. Setup and notation

Notation follows `docs/multiphase.tex` and the code.

- The cloud has $K$ bound modes, labelled by their linear states $|nlm\rangle$.
  Mode $i$ has amplitude $A_i$, real frequency $\omega_i$, azimuthal number
  $m_i$, and phase $\xi_i = \vartheta_i - m_i\varphi$. Here $\varphi$ is the
  Boyer–Lindquist azimuth and $\vartheta_i = \int_0^t \omega_i\,dt'$.
- The field is a sum of Fourier components labelled by integer vectors
  $\mathbf q = (q_1,\dots,q_K)$. Each varies as
  $e^{-i\mathbf q\cdot\boldsymbol\xi}$, so it has frequency
  $\omega_{\mathbf q} = \mathbf q\cdot\boldsymbol\omega$ and azimuthal number
  $m_{\mathbf q} = \mathbf q\cdot\mathbf m$
  (`eq:spatial_profile_Fq`).
- Component $\mathbf q$ has spatial profile $F_{\mathbf q}(r,\theta)$, expanded
  in polar harmonics with radial functions $R_{l,\mathbf q}(r_*)$. In the code,
  $Q_{l,\mathbf q} = dR_{l,\mathbf q}/dr_*$ is called `slopes`.

For the production two-mode set ($|211\rangle$, $|322\rangle$; `max_order=3`,
even potential), the retained components are the fundamentals $(1,0)$ and
$(0,1)$, the second-order components $(2,-1)$ and $(-1,2)$, and the third-order
components $(3,0)$, $(2,1)$, $(1,2)$ and $(0,3)$. Only one of each $\pm\mathbf q$
pair is kept, the one with $\omega_{\mathbf q}>0$
([`src/modes.jl`](../src/modes.jl)).

## 2. What a wave action is

In Whitham's averaged-Lagrangian method, the Lagrangian is averaged over all
phases $\boldsymbol\xi$. The result, $\langle L\rangle(\boldsymbol\omega,
\mathbf A,\text{background})$, does not depend on the phases themselves. So each
shift $\xi_i \to \xi_i + c$ is a symmetry of the averaged problem. Its Noether
charge is the **wave action** of phase $i$:

$$
N_i = \frac{\partial\langle L\rangle}{\partial\omega_i}.
$$

### Computing it from the components

In $\langle L\rangle$, the frequencies enter only through the kinetic term of
each component. Differentiating component $\mathbf q$'s kinetic term with
respect to $\omega_{\mathbf q}$ gives its Klein–Gordon charge $n_{\mathbf q}$.
In the code (`component_action_densities`), its density per unit tortoise
coordinate is

$$
\frac{dn_{\mathbf q}}{dr_*} = 2\pi \sum_{x} w_x\,\rho^2\,\frac{\Delta}{r^2+a^2}\;
2|F_{\mathbf q}(x)|^2\,\bigl(-g^{tt}\omega_{\mathbf q} + g^{t\varphi} m_{\mathbf q}\bigr),
$$

where $x=\cos\theta$, the $w_x$ are Gauss–Legendre weights, and
$\rho^2 = r^2 + a^2x^2$. Since $\partial\omega_{\mathbf q}/\partial\omega_i = q_i$,
the chain rule gives

$$
N_i = \sum_{\mathbf q} q_i\, n_{\mathbf q} \qquad (\texttt{phase\_actions}).
$$

So $N_i$ counts quanta of phase $i$. A component with $q_i = 2$ counts twice. For
example, the $(2,-1)$ component (frequency $2\omega_1-\omega_2$, $m = 0$) adds
$+2\,n_{(2,-1)}$ to $N_1$ and $-n_{(2,-1)}$ to $N_2$.

### Why the sums are conserved but the individual charges are not

The nonlinear source moves charge between components. For example, two
$|211\rangle$ quanta can become one $|322\rangle$ quantum plus one $(2,-1)$
quantum. So no single $n_{\mathbf q}$ is conserved. The averaged potential is
nevertheless invariant under every phase shift. This makes the $q_i$-weighted
exchange rates sum to zero, and the transfer cancels out of each $N_i$. This
cancellation is what makes $N_i$ a balance variable.

## 3. Physical picture

Notation added in this section (mine, not the outline's): $\mu$ is the axion
mass, so $\mu M = \alpha$ and $\mu = \alpha$ in the code's $M=1$ units.
$\Omega_H = a/(r_+^2+a^2)$ is the horizon angular velocity.
$k_{\mathbf q} = \omega_{\mathbf q} - m_{\mathbf q}\Omega_H$ is a component's
wavenumber at the horizon.

### Actions count quanta

For a single linear mode, energy and angular momentum are fixed multiples of
the Klein–Gordon charge: $E = \omega N$ and $J = mN$. Each quantum carries
energy $\omega$ and angular momentum $m$. With $\hbar$ restored, $N/\hbar$ is the
number of axions in the mode. This is why superradiance always extracts energy
and angular momentum from the hole in the ratio $\omega/m$.

A self-interacting real scalar does not conserve particle number in general. A
quartic interaction allows $2\to2$, $1\to3$, $3\to1$ and $4\to0$ processes.
After the phase average, only processes whose frequencies are built from the
$K$ phases appear in the solution, and each of them conserves every $N_i$.
Two examples from the production set:

- **The harmonic $(3,0)$** (frequency $3\omega_1$) is three $|211\rangle$ quanta
  combined into one. It counts $3$ toward $N_1$, so generating it does not change
  $N_1$. In the action count, harmonics are part of the nonlinear mode, not
  separate excitations.
- **The component $(2,-1)$** is produced by $211 + 211 \to 322 + (2,-1)$.
  The energies match: $2\omega_1 = \omega_2 + (2\omega_1-\omega_2)$. The process
  removes two quanta from mode 1 and adds one to mode 2. The new $(2,-1)$
  quantum counts $+2$ toward $N_1$ and $-1$ toward $N_2$, so neither total changes.

Nonlinearly, then, $N_i$ is the number of quanta that belong to phase $i$,
wherever they are stored: in the fundamental, in harmonics, or in mixed
components.

The local density shows where the quanta sit. Its frequency factor is

$$
-g^{tt}\omega_{\mathbf q} + g^{t\varphi}m_{\mathbf q}
= -g^{tt}\bigl(\omega_{\mathbf q} - m_{\mathbf q}\,\Omega_\text{fd}\bigr),
\qquad \Omega_\text{fd} = \frac{g^{t\varphi}}{g^{tt}} = -\frac{g_{t\varphi}}{g_{\varphi\varphi}},
$$

where $\Omega_\text{fd}(r,\theta)$ (my notation) is the frame-dragging angular
velocity of zero-angular-momentum observers. Since $-g^{tt}>0$ outside the
horizon, the density is positive wherever the pattern speed
$\omega_{\mathbf q}/m_{\mathbf q}$ exceeds the local frame dragging. The frame
dragging tends to $\Omega_H$ at the horizon. So a superradiant component
($\omega_{\mathbf q} < m_{\mathbf q}\Omega_H$) has negative action density in a
thin layer next to the horizon. This is the local version of the horizon flux
sign below.

### Actions change only at the boundaries

The bulk exchange conserves every $N_i$, so $N_i$ can change only through
flux across the inner (horizon-side) and outer boundaries. Superradiant growth,
transfer between modes, and radiation all appear as boundary fluxes of
particular components.

**Horizon.** Each component is ingoing at the horizon,
$R_{l,\mathbf q} \propto e^{-ik_{\mathbf q}r_*}$. Then
$\operatorname{Im}(R^*Q) = -k_{\mathbf q}|R|^2$, and the component contributes

$$
+\,4\pi(r_+^2+a^2)\,k_{\mathbf q}\sum_l |R_{l,\mathbf q}|^2
$$

to its loss. The sign is that of $k_{\mathbf q}$:

- $k_{\mathbf q} > 0$: quanta fall into the hole (absorption).
- $k_{\mathbf q} < 0$, i.e. $\omega_{\mathbf q} < m_{\mathbf q}\Omega_H$: the loss
  is negative and quanta are created. This is superradiance, here stated for
  each component separately.

The hole receives energy $\omega_{\mathbf q}$ and angular momentum
$m_{\mathbf q}$ per absorbed quantum. So
$dM - \Omega_H\,dJ_\text{BH} = \sum_{\mathbf q} k_{\mathbf q}\,dn^\text{hole}_{\mathbf q}
\propto \sum_{\mathbf q} k_{\mathbf q}^2 |R_{\mathbf q}|^2 \ge 0$:
the area theorem holds component by component.

**Outer boundary.** Components with $\omega_{\mathbf q} < \mu$ are bound. They
decay as $e^{-\kappa r}$ and carry no flux out of a domain that contains them.
Components with $\omega_{\mathbf q} > \mu$ propagate. Their outgoing flux is
positive, so they only ever remove quanta.

**The transfer channels at $a=0.9$, $\alpha=0.3$** (component diagnostics in
`NOTES.md`):

| Component | $m_{\mathbf q}$ | Behaviour | Process | $(\Delta N_1, \Delta N_2)$ per quantum lost |
|---|---|---|---|---|
| $(2,-1)$ | 0 | bound; absorbed at the horizon ($k_{\mathbf q}=\omega_{\mathbf q}>0$ always) | $211+211 \to 322 + \text{BH}$ | $(-2, +1)$ |
| $(-1,2)$ | 3 | propagating; $\omega_{\mathbf q}$ just above $\mu$ (hydrogenic margin $\alpha^3/72$) | $322+322 \to 211 + \text{outgoing axion}$ | $(+1, -2)$ |
| $(3,0)$, $(2,1)$, $(1,2)$, $(0,3)$ | 3–6 | propagating, $\omega_{\mathbf q}\approx3\mu$ | harmonics radiated to infinity | fluxes $\le 10^{-12}$, negligible |

At the joint test point in `NOTES.md`, the horizon absorption of $(2,-1)$ is the
dominant loss channel. In the action balance this one flux both drains $|211\rangle$
and pumps $|322\rangle$. No separate interaction rate appears anywhere: the
transfer between modes is the $q_i$-weighted boundary flux of the driven
components.

**Relation to the perturbative rate equations.** At leading order in the
self-interaction, the $(2,-1)$ profile is proportional to $A_1^2A_2$, so its flux
scales as $N_1^2N_2$. Likewise the $(-1,2)$ flux scales as $N_1N_2^2$. The action
balance then reduces to

$$
\begin{aligned}
\dot N_1 &= \Gamma_1 N_1 - 2c_H N_1^2N_2 + c_\infty N_1N_2^2,\\
\dot N_2 &= \Gamma_2 N_2 + c_H N_1^2N_2 - 2c_\infty N_1N_2^2,
\end{aligned}
$$

with linear growth rates $\Gamma_i = 2\operatorname{Im}\omega_i$ (since
$N_i\propto A_i^2$) and non-negative coefficients $c_H$, $c_\infty$ from the horizon
and outer fluxes. These symbols are mine, for illustration. This is the structure
of the occupation-number equations in the perturbative literature, for example
Baryakhtar, Galanis, Lasenby & Simon (2021), arXiv:2011.11646. The code does
not make this expansion. It evaluates the fluxes on the full nonlinear solution,
so frequency shifts, saturation of the driven components and the higher harmonics
are all included.

### Why the actions are adiabatic invariants

The classical analogue is a pendulum whose length is changed slowly. Its energy
changes, but the action $E/\omega$ does not. Here the black-hole mass and spin
change on the spindown time scale, which is very long compared with $1/\omega$.
The mode frequencies $\omega_i(b,j)$ follow the background, but the number of
quanta changes only through the boundary fluxes. So the hole changes the energy
per quantum, not the number of quanta, and the cloud's energy changes by
$\sum_i N_i\,d\omega_i$. This is the metric work $w$ of §6 in the linear regime.

For the sign: at fixed $\mu$, a heavier hole binds more deeply (hydrogenically
$\omega \approx \mu(1-\alpha^2/2n^2)$ with $\alpha = \mu M$). So $w$ has the
opposite sign to $b'$, which is what the validation finds ($w/(E\,b')<0$).

The argument assumes the phases stay independent. If some retained
combination $\mathbf q\cdot\boldsymbol\omega$ passed through zero, the average
over $\boldsymbol\xi$ would not be valid near that point.

## 4. Relation to energy and angular momentum

Component $\mathbf q$'s energy and angular momentum are $\omega_{\mathbf q}$ and
$m_{\mathbf q}$ times its charge, both in the bulk and in the flux. This gives

$$
J = \sum_i m_i N_i \quad\text{(exact)},\qquad
L_E = \sum_i \omega_i L_{N_i},\qquad
L_J = \sum_i m_i L_{N_i}.
$$

Here $L_E$, $L_J$ and $L_{N_i}$ are the signed net losses through both
boundaries (§5).

For the energy itself, $E \ne \sum_i\omega_i N_i$ exactly. The difference is a
surface term at the outer endpoint plus potential-energy contributions. At the
linear test point it is $-3.6\times10^{-4}E$ at $r_{*,\max}=180$,
$-4\times10^{-5}E$ at $250$, and $-3\times10^{-6}E$ at $320$. On a fixed
background the differential form holds once the cloud is contained in the
domain:

$$
dE = \sum_i \omega_i\, dN_i .
$$

This is the first law of the averaged problem, the multi-phase version of
$E = \omega N$ for a single mode. Adding one quantum to mode $i$ costs energy
$\omega_i$, including that mode's nonlinear frequency shift.

**Where these relations come from.** In Whitham's theory the averaged energy is
the Legendre transform of the averaged Lagrangian,

$$
E = \sum_i \omega_i N_i - \mathcal L, \qquad \mathcal L = \int \langle L\rangle,
$$

the same way a Hamiltonian is $p\dot q - L$. Two consequences follow:

- **The first law.** Vary along the solution family:
  $dE = \sum_i\omega_i\,dN_i + \sum_i N_i\,d\omega_i - d\mathcal L$. Since
  $\partial\mathcal L/\partial\omega_i = N_i$ and $\mathcal L$ is stationary with
  respect to the profiles on a solution, $d\mathcal L = \sum_i N_i\,d\omega_i$
  plus boundary terms. So $dE = \sum_i\omega_i\,dN_i$ up to boundary terms,
  which vanish once the cloud is contained.
- **The offset $E - \sum_i\omega_iN_i = -\mathcal L$.** On a solution, integrating
  the kinetic term by parts with the field equation leaves a surface term plus a
  term of the form $\int(V - \tfrac12\Phi V'(\Phi))$. That term is zero for a free
  field and starts at quartic order in the amplitude otherwise. The free-field
  numbers above shrink as the outer endpoint moves out, as a surface term
  should.

The code does not compute $\mathcal L$. These are the standard relations that the
checks in §7 are consistent with.

## 5. Fluxes and the balance

The outward action flux of phase $i$ through the sphere at a mesh node is

$$
F_{N_i} = \sum_{\mathbf q} q_i\; 4\pi(r^2+a^2)\sum_{l}
\operatorname{Im}\!\bigl[R^*_{l,\mathbf q}\,Q_{l,\mathbf q}\bigr]
\qquad(\texttt{action\_flux}).
$$

Each term is a radial Wronskian. The net loss is
$L_{N_i} = F_{N_i}(\text{outer}) - F_{N_i}(\text{inner})$ (`action_losses`), and
the balance is

$$
\frac{dN_i}{dt} = -L_{N_i}.
$$

At each instant the cloud is the stationary multiphase solution at the current
amplitudes. So $N_i$ depends on $t$ only through $\mathbf A$, plus the background
$b = M/M_0$, $j = J_\text{BH}/M_0^2$ when the black hole spins down. The
derivatives $\partial N_j/\partial A_i$ include the response of the whole coupled
solution: profiles, frequencies and driven components. They are computed by
the implicit-function theorem (`action_jacobian`, the same machinery as
`charge_jacobian`).

**Fixed background:**

$$
\bigl[\partial_{\mathbf A} N\bigr]\,\mathbf A' = -L_N .
$$

**With spindown** (primes are $d/d\tau$, $\tau = t/M_0$; reference units
$N_\text{ref} = b^2 N$ like $J$, and losses $b$ times the current-unit values):

$$
\bigl[\partial_{\mathbf A} N\bigr]\,\mathbf A' = -L_N - \bigl(\partial_b N\, b' + \partial_j N\, j'\bigr).
$$

The background derivatives are taken at fixed amplitudes. For comparison, the
energy balance (`eq:evolving_background_balance`) is

$$
\begin{pmatrix}\partial_{\mathbf A}E\\ \partial_{\mathbf A}J\end{pmatrix}\mathbf A'
= -\begin{pmatrix}L_E\\ L_J\end{pmatrix} + \begin{pmatrix}w\\0\end{pmatrix}
- \begin{pmatrix}\partial_b E\,b' + \partial_j E\,j'\\ \partial_b J\,b' + \partial_j J\,j'\end{pmatrix}.
$$

## 6. Why use it

| | Energy balance (`:energy`) | Action balance (`:action`) |
|---|---|---|
| Equations | 2 ($E$, $J$) | $K$, one per phase |
| Bound modes allowed | exactly 2 | any number |
| Metric work $w$ | needed | not needed: actions are adiabatic invariants |
| Flux | needs the growth-rate term (below) | contains no time derivative |

### How the two balances are related

When $dE = \sum_i\omega_i\,dN_i$ holds, which requires a contained cloud, and
using $J = \sum_i m_iN_i$, the energy-balance quantities are linear combinations
of the action-balance quantities:

$$
\begin{pmatrix}\partial_{\mathbf A}E\\ \partial_{\mathbf A}J\end{pmatrix} = W\,\partial_{\mathbf A}N,
\qquad
\begin{pmatrix}L_E\\ L_J\end{pmatrix} = W\,L_N,
\qquad
W = \begin{pmatrix}\omega_1 & \cdots & \omega_K\\ m_1 & \cdots & m_K\end{pmatrix}.
$$

For two modes the energy balance is the action balance multiplied by the
$2\times2$ matrix $W$. The two have the same solution as long as
$\det W = \omega_1m_2 - \omega_2m_1 \neq 0$. For $|211\rangle$ and $|322\rangle$,
$\det W = 2\omega_1-\omega_2$, which is the frequency of the $(2,-1)$ component,
so it is far from zero. For $K>2$, $W$ is $2\times K$: energy and angular
momentum fix only two combinations of the $K$ action rates, and the remaining
$K-2$ directions are undetermined. The two balances differ only where
$dE = \sum_i\omega_i\,dN_i$ fails, which is on domains that do not contain the
cloud.

### Metric work

When the black hole changes slowly, it does work $w$ on the cloud. The energy
balance needs $w$ explicitly. The actions are adiabatic invariants, so they
change only through boundary fluxes, and their balance has no work term. The
same fact gives $w$ for the energy balance (`metric_work`):

$$
w = \Bigl(\partial_b E - \sum_i\omega_i\,\partial_b N_i\Bigr) b'
  + \Bigl(\partial_j E - \sum_i\omega_i\,\partial_j N_i\Bigr) j',
$$

with all derivatives at fixed amplitudes. The angular-momentum counterpart
vanishes identically because $J = \sum_i m_i N_i$, which matches the $(w, 0)$ in
`eq:evolving_background_balance`. In the linear regime ($E_i = \omega_i N_i$)
this reduces to $w = \sum_i N_i\, d\omega_i/d\tau$. This version needs no metric
derivative near the moving horizon, because the near-horizon densities of $E$
and $\sum_i\omega_iN_i$ cancel. The earlier Boyer–Lindquist expression
(`boyer_lindquist_metric_work`, kept for comparison) needed that derivative
and depended strongly on the inner endpoint.

### The energy-flux error

During the solve, each component has a complex trial frequency
$\tilde\omega_{\mathbf q}$ (my notation here: $\operatorname{Re}\tilde\omega =
\omega_{\mathbf q}$, and $\operatorname{Im}\tilde\omega$ is its growth rate). With
a complex frequency, the complete energy flux is

$$
F_E = 4\pi(r^2+a^2)\sum_{\mathbf q,l}\Bigl[\omega_{\mathbf q}\operatorname{Im}(R^*Q)
 - \operatorname{Im}\tilde\omega_{\mathbf q}\,\operatorname{Re}(R^*Q)\Bigr].
$$

`outward_flux` keeps only the first term. Since
$\operatorname{Re}(R^*Q) = \tfrac12\,d|R|^2/dr_*$, the dropped term is small only
where the profile has decayed. At an outer endpoint where the cloud is not
contained, it is as large as the losses. In the nonlinear adiabatic setting it
is also unclear which growth rate belongs in that term. The action flux is a
pure Wronskian and has no such term.

The difference comes from what each flux is built from. The action flux is the
radial component of the Klein–Gordon current,
$j^r \propto \operatorname{Im}(\Phi^*\nabla^r\Phi)$, which contains only a radial
derivative. The energy flux is $T^r{}_t \propto \operatorname{Re}(\nabla^r\Phi^*\,\partial_t\Phi)$.
Its time derivative brings in $-i\tilde\omega$, and with a complex frequency the
imaginary part produces the second term. Physically, the extra term is the
energy that a growing (or decaying) cloud stores in, or takes from, the region
inside the sphere. A stationary-flux formula leaves it out.

## 7. Validation

Run the full check with `julia --project=. validation/action_balance.jl`.
Parameters: $a=0.9$, $\alpha=0.3$, cosine potential, modes $|211\rangle$ and
$|322\rangle$.

| Check | Result | What it establishes |
|---|---|---|
| $J=\sum m_iN_i$, $L_E=\sum\omega_iL_{N_i}$, $L_J=\sum m_iL_{N_i}$ | round-off, at a linear and a nonlinear point | the action bookkeeping matches the existing charge and flux code |
| Linear point, default mesh $r_*\in[-80,180]$: $A'/(\gamma A)-1$ | energy balance as coded: $(7.4\times10^{-5},\,-4.2\times10^{-4})$; energy balance, complete flux: $(-1.5\times10^{-6},\,-1.2\times10^{-5})$; action balance: $(-1.5\times10^{-6},\,-1.2\times10^{-5})$ | the action balance recovers $A'=\gamma A$ and agrees with the corrected energy balance; the coded energy balance carries the flux error of §6 |
| Nonlinear point $\mathbf A=(10^{-2},10^{-4})$, default mesh | the three balances agree to $10^{-3}$ | consistency only; there is no exact answer at this point |
| Metric work vs inner endpoint, $r_{*,\min} = -45\ldots-120$ | $w/(E\,b') = -4.011009\times10^{-3}$, constant to 7 digits (the Boyer–Lindquist version runs from $-10^{-2}$ to $+5\times10^5$) | the action-based $w$ is insensitive to the inner endpoint |
| Metric work vs outer endpoint, $r_{*,\max}=180,250,320,400$ | $w/(E\,b') = -4.01, -9.52, -10.561, -10.660$ (all $\times10^{-3}$) vs linear prediction $-10.664\times10^{-3}$ | $w$ matches $\sum_iN_i\,d\omega_i/d\tau$ once the cloud is contained |
| Energy balance $+\,w$ vs action balance, contained domain | agree to $10^{-7}$ (without $w$: $10^{-4}$) | the two formulations agree when both are complete |

In $\gamma = \operatorname{Im}\omega$, $\gamma$ is the linear growth rate of each
mode, so $A'=\gamma A$ is the exact answer at the linear point.

`Pkg.test()` repeats the identities to $10^{-12}$, $A'=\gamma A$ to $10^{-4}$,
$dE=\sum\omega\,dN$ to $10^{-3}$ on a mesh with $r_{*,\max}=250$, and a short
log-amplitude evolution with `balance=:action` against $A\,e^{\gamma t}$.

**What these checks do not establish.** They test consistency and the linear
limit. They do not test convergence of nonlinear rates in radial resolution,
domain size, or harmonic truncation. The action balance has only been run with
two bound modes. The "any number of modes" claim is structural (one equation per
phase, no two-mode restriction in the code) and has not been tested with
$K\ge3$.

## 8. Using it

```julia
using AxionCloud

# Rates at one solution
rates = amplitude_rates(solution; balance=:action)

# Fixed-background evolution
result = evolve_fixed_background(problem, amplitudes, (0.0, final_time); balance=:action)

# With spindown (decay_constant_ratio = f_a / M_Pl)
result = evolve_with_spindown(problem, amplitudes, (0.0, final_time), decay_constant_ratio; balance=:action)

trajectory = amplitude_trajectory(result)
```

`balance=:energy` raises an `ArgumentError` for anything other than two bound
modes. $E$ and $J$ are still recorded in the evolution diagnostics, so they
remain available as checks.

## 9. Open questions

From `NOTES.md` (Day 3):

1. **Make `:action` the default?** Done on 2026-11-09. It is complete on a
   bounded domain, needs no metric work, removes the two-mode restriction, and
   avoids the energy-flux error.
2. **Long $a=0.9$ evolution.** It should be rerun on a contained domain with
   integrated outer data, preferably with `balance=:action`. The Day 2
   saturation numbers ($r_{*,\max}=180$, series boundary) should not be quoted
   until then.
3. **Drift terms on non-contained domains.** The background derivatives
   $\partial_b Q|_{\mathbf A}$ are large when the cloud is not contained,
   because the fixed outer coefficient rescales the whole tail. Using
   contained domains removes this problem.

## Code map

| Function | File | Role |
|---|---|---|
| `component_action_densities` | `src/charges.jl` | $dn_{\mathbf q}/dr_*$ per component |
| `phase_weights` | `src/charges.jl` | matrix $[q_i]$ mapping components to phases |
| `phase_actions` | `src/charges.jl` | $N_i$ (Simpson in $r_*$) |
| `action_flux`, `action_losses` | `src/charges.jl` | $F_{N_i}$ at a node, $L_{N_i}$ |
| `action_jacobian` | `src/charges.jl` | $\partial N_j/\partial A_i$ by the implicit-function theorem |
| `amplitude_rates(...; balance)` | `src/evolution.jl` | solves either balance |
| `metric_work` | `src/evolution.jl` | action-based $w$ |
| `reference_actions` | `src/evolution.jl` | $N_\text{ref} = b^2 N$ |
| `boyer_lindquist_metric_work` | `src/evolution.jl` | superseded $w$, kept for comparison |
