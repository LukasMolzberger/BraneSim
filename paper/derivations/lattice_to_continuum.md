# Lattice → continuum: long-wavelength elastic constants of the BraneSim cubic lattice

Sprint 1, subtask 1 of `paper-v4/validation_roadmap.md`.
Author: `physics-derivation` agent.

---

## 1. Layer interface

- **Microscopic layer (UV).** State: node positions $\{R_p\}_{p\in\mathbb{Z}^3}\subset\mathbb{R}^3$
  (only the in-brane / lateral components $\xi^i$; the $X^4$ amplitude direction is
  ignored by assumption A1 below). Energy:
  $$
    U_{\mathrm{lat}}=\tfrac12\sum_{n}\sum_{\delta\in\mathcal N} k_\delta\Big(|R_{n+\delta}-R_n|-\alpha a|\delta|\Big)^2,\qquad
    k_\delta=\frac{k_0}{|\delta|^2},
  $$
  with the 26-neighbor stencil
  $\mathcal N=\{\pm\hat e_i\}\cup\{\pm\hat e_i\pm\hat e_j:i<j\}\cup\{\pm\hat e_1\pm\hat e_2\pm\hat e_3\}$.
- **Continuum layer (IR).** State: vector field $u(x,t)\in\mathbb{R}^3$ on
  $\Omega\subset\mathbb{R}^3$ with Lagrangian density
  $\mathcal L=\tfrac12\rho_m|\partial_t u|^2-W(\nabla u)$ where $W$ is a quadratic form
  in $\partial u$.
- **Coarse-graining map.** Identify lattice site $n$ with continuum point $x=an$ and
  set $R_n(t)=an+u(an,t)$ (Eulerian/material near-identity around the **held**
  configuration; see assumption A2 — the prompt's $\alpha a n$ choice would expand
  around the stress-free state and is incompatible with the prestress story).
  Continuum fields are obtained by Taylor expansion in the small parameter $a|k|\ll1$.

---

## 2. Assumptions

1. **A1 — In-brane channel only.** We restrict to lateral displacements
   $u\in\mathbb{R}^3$ and ignore the $X^4$ (amplitude) channel. The geometric
   $\partial_i u\,\partial_j u$ coupling between lateral and amplitude appears at
   the same order in $\partial u$ but is orthogonal to the elastic-constant
   identification done here. (Stated in the task.)
2. **A2 — Reference configuration is the held (prestressed) lattice.** We expand
   $R_n=an+u(an)$. The prompt's $R_p=\alpha a n+u$ would expand around the
   stress-free configuration and would erase the prestress couplings the
   derivation is meant to expose; we therefore use the held configuration. This
   matches `BraneGrid3D.get_spatial_coordinates` (held spacing = `spacing`) and
   the convention of paper §3.
3. **A3 — Long-wavelength regime.** $a|k|\ll1$, equivalently $u$ varies slowly
   over $a$. Higher-order gradients ($a^2\partial^2 u$ etc.) are dropped.
4. **A4 — Quadratic in $\partial u$.** We retain the elastic order $(\partial u)^2$
   and discard $(\partial u)^3$ and higher. This is the linear-elasticity
   regime; geometric quartic terms (Skyrme-class) live at the next order and are
   the subject of a separate sprint.
5. **A5 — Bulk lattice / no boundary.** Surface terms drop after summing over
   $\pm\delta$ pairs.
6. **A6 — Wave-equation moduli.** When $\alpha<1$ the reference configuration is
   prestressed; the second derivative of $U$ at fixed $\partial u$ then mixes
   "pure elastic" and "geometric stiffness" contributions. Because the
   physically observable quantities (sound speeds) come directly from the
   **acoustic / wave-operator** tensor $A_{abcd}$, we report Lamé parameters
   from the wave dispersion. The split between "Lamé from static strain energy"
   and "geometric stiffness" is recorded in §3 for completeness but not used in
   §4.
7. **A7 — Inversion-symmetric stencil.** Each $\delta\in\mathcal N$ comes with
   $-\delta\in\mathcal N$ and equal weights; this kills the linear-in-$\partial u$
   terms in the bulk.

---

## 3. Derivation

### 3.1 Link strain to quadratic order in $\partial u$

Let $\Delta u\equiv u(an+a\delta)-u(an)$. Then
$|R_{n+\delta}-R_n|=|a\delta+\Delta u|$ and
$$
|a\delta+\Delta u|=a|\delta|\sqrt{1+\tfrac{2\delta\cdot\Delta u}{a|\delta|^2}+\tfrac{|\Delta u|^2}{a^2|\delta|^2}}
=a|\delta|+\tfrac{\delta\cdot\Delta u}{|\delta|}+\tfrac{|\Delta u|^2-(\hat\delta\cdot\Delta u)^2}{2a|\delta|}+\mathcal O\!\big(\tfrac{(\Delta u)^3}{a^2}\big),
$$
with $\hat\delta\equiv\delta/|\delta|$. The link strain is
$$
s_\delta=|R_{n+\delta}-R_n|-\alpha a|\delta|
=(1-\alpha)\,a|\delta|+\tfrac{\delta\cdot\Delta u}{|\delta|}+\tfrac{|\Delta u_\perp|^2}{2a|\delta|}+\cdots,
$$
where $|\Delta u_\perp|^2=|\Delta u|^2-(\hat\delta\cdot\Delta u)^2$ is the
component transverse to the bond.

Squaring and grouping by powers of $\Delta u$ (the $(1-\alpha)^2 a^2|\delta|^2$
zeroth-order piece is the constant background energy; the linear-in-$\Delta u$
piece sums to a total derivative by A7):
$$
\big(s_\delta^2\big)^{(2)}
=\big[\hat\delta\cdot\Delta u\big]^2+(1-\alpha)\,|\Delta u|^2
\;-\;(1-\alpha)\,(\hat\delta\cdot\Delta u)^2
=\alpha\,(\hat\delta\cdot\Delta u)^2+(1-\alpha)\,|\Delta u|^2.
$$
The second equality is the clean rewrite that exposes the role of $\alpha$:

- The longitudinal (along-bond) channel carries the full Hookean response with
  strength $\alpha$ (because the linear-in-$\Delta u$ stretch survives squared,
  and the $(1-\alpha)$-weighted "geometric" piece subtracts off its longitudinal
  part).
- The transverse (off-bond) channel is purely **geometric stiffness from
  prestress**, weighted by $(1-\alpha)$: a prestretched spring resists transverse
  motion of its endpoints to second order.

### 3.2 Continuum gradient expansion

Expand $\Delta u_a\simeq a\,\delta_b\partial_b u_a+\mathcal O(a^2)$. Then
$$
(\hat\delta\cdot\Delta u)^2=a^2|\delta|^2\,\hat\delta_a\hat\delta_b\hat\delta_c\hat\delta_d\,(\partial_b u_a)(\partial_d u_c),\qquad
|\Delta u|^2=a^2|\delta|^2\,\hat\delta_b\hat\delta_d\,\delta_{ac}(\partial_b u_a)(\partial_d u_c).
$$
With $k_\delta=k_0/|\delta|^2$ the per-link factor $k_\delta a^2|\delta|^2=k_0a^2$
becomes shell-independent, and the total quadratic energy density (per
$a^3$ volume, with the explicit $\tfrac12$ for bond double-counting) is
$$
\boxed{\;
W(\nabla u)=\frac{k_0}{4a}\sum_{\delta\in\mathcal N}\Big\{\alpha\,\hat\delta_a\hat\delta_b\hat\delta_c\hat\delta_d
+(1-\alpha)\,\hat\delta_b\hat\delta_d\,\delta_{ac}\Big\}\,(\partial_b u_a)(\partial_d u_c).
\;}
$$
Define the lattice tensor sums
$$
T^{(2)}_{ab}\equiv\sum_{\delta\in\mathcal N}\hat\delta_a\hat\delta_b,\qquad
T^{(4)}_{abcd}\equiv\sum_{\delta\in\mathcal N}\hat\delta_a\hat\delta_b\hat\delta_c\hat\delta_d.
$$

### 3.3 Lattice tensor sums for the 26-neighbor stencil

Shell-by-shell ($N_{\rm I}=6$ axial, $N_{\rm II}=12$ face-diagonal,
$N_{\rm III}=8$ body-diagonal):

**$T^{(2)}_{ab}$.** Each shell gives $T^{(2)}|_{\rm shell}\propto\delta_{ab}$:
$T^{(2)}|_{\rm I}=2\delta_{ab}$, $T^{(2)}|_{\rm II}=4\delta_{ab}$,
$T^{(2)}|_{\rm III}=\tfrac{8}{3}\delta_{ab}$. Total:
$$
T^{(2)}_{ab}=\Big(2+4+\tfrac{8}{3}\Big)\delta_{ab}=\tfrac{26}{3}\,\delta_{ab}.
$$

**$T^{(4)}_{abcd}$.** Define the cubic-anisotropy projector
$Q_{abcd}\equiv\sum_k\delta_{ka}\delta_{kb}\delta_{kc}\delta_{kd}$ ($=1$ iff $a{=}b{=}c{=}d$).
- Shell I: $T^{(4)}|_{\rm I}=2Q_{abcd}$ (only axial bonds, all four indices equal).
- Shell II: $T^{(4)}|_{\rm II}=(\delta_{ab}\delta_{cd}+\delta_{ac}\delta_{bd}+\delta_{ad}\delta_{bc})-Q_{abcd}$
  (the sign-sum kills odd powers, leaving the three pair-contractions and a
  $-Q$ correction so that the all-equal entry is $2$, not $3$).
- Shell III: $T^{(4)}|_{\rm III}=\tfrac{8}{9}\big(\delta_{ab}\delta_{cd}+\delta_{ac}\delta_{bd}+\delta_{ad}\delta_{bc}-2Q_{abcd}\big)$.

Sum:
$$
\boxed{\;
T^{(4)}_{abcd}=\frac{17}{9}\big(\delta_{ab}\delta_{cd}+\delta_{ac}\delta_{bd}+\delta_{ad}\delta_{bc}\big)\;-\;\frac{7}{9}\,Q_{abcd}.
\;}
$$
The $Q_{abcd}$ term is the **only** non-isotropic piece: a fully isotropic
4-tensor would consist solely of the symmetric pair-contractions. The negative
sign means the lattice is "stiffer along axes than off-axis" relative to the
isotropic baseline.

Numerical entries needed below:
$T^{(4)}_{1111}=3\cdot\tfrac{17}{9}-\tfrac{7}{9}=\tfrac{44}{9}$,
$T^{(4)}_{1122}=T^{(4)}_{1212}=\tfrac{17}{9}$.

### 3.4 Acoustic tensor and elastic constants

Write $W=\tfrac12 A_{abcd}(\partial_b u_a)(\partial_d u_c)$ with
$$
A_{abcd}=\frac{k_0}{2a}\Big[\alpha\,T^{(4)}_{abcd}+(1-\alpha)\,T^{(2)}_{bd}\,\delta_{ac}\Big].
$$
The wave equation $\rho_m\,\partial_t^2 u_a=\partial_b A_{abcd}\partial_d u_c$
gives the Christoffel matrix $M_{ac}(k)=A_{abcd}k_b k_d$ for plane waves
$u\propto e^{i(k\cdot x-\omega t)}$. For $k=k\hat e_1$:

- **Longitudinal** ($u\parallel\hat e_1$):
  $\displaystyle M_{11}=A_{1111}\,k^2=\frac{k_0}{2a}\Big[\alpha\,\tfrac{44}{9}+(1-\alpha)\,\tfrac{26}{3}\Big]k^2
  =\frac{k_0(78-34\alpha)}{18a}\,k^2.$

- **Transverse** ($u\perp\hat e_1$, e.g. $u\parallel\hat e_2$):
  $\displaystyle M_{22}=A_{2121}\,k^2=\frac{k_0}{2a}\Big[\alpha\,\tfrac{17}{9}+(1-\alpha)\,\tfrac{26}{3}\Big]k^2
  =\frac{k_0(78-61\alpha)}{18a}\,k^2.$

In the standard isotropic identification $\rho_m c_L^2=\lambda+2\mu$,
$\rho_m c_T^2=\mu$:

$$
\boxed{\;\;
\mu(\alpha)=\frac{k_0\,(78-61\alpha)}{18\,a},\qquad
\lambda(\alpha)+2\mu(\alpha)=\frac{k_0\,(78-34\alpha)}{18\,a},\qquad
\lambda(\alpha)=\frac{k_0\,(88\alpha-78)}{18\,a}.
\;\;}
$$

### 3.5 Mass density

The continuum kinetic term $\tfrac12\rho_m|\partial_t u|^2\,d^3x$ comes from the
discrete $\tfrac12 m|\dot u_n|^2$ on a lattice with one site per cell of volume
$a^3$:
$$
\boxed{\;\rho_m=\frac{m}{a^3}.\;}
$$
In `orchestration/configs/local_extensive.json` the field `mass_density` directly
sets $\rho_m$; with `spacing=1.0` the per-node mass equals $\rho_m$ numerically.

### 3.6 Static elastic tensor and Cauchy / isotropy

Under the symmetric-strain identification
$C_{ijkl}=\tfrac14(A_{ijkl}+A_{jikl}+A_{ijlk}+A_{jilk})$ the symmetrized form is
$$
C_{abcd}=\frac{k_0}{2a}\Big[\alpha\,T^{(4)}_{abcd}+(1-\alpha)\,\tfrac{13}{3}\big(\delta_{ac}\delta_{bd}+\delta_{ad}\delta_{bc}\big)\Big],
$$
whose cubic Voigt entries are
$$
\begin{aligned}
C_{1111}&=\frac{k_0(78-34\alpha)}{18a},\\
C_{1122}&=\frac{17\alpha\,k_0}{18a},\\
C_{1212}&=\frac{k_0(39-22\alpha)}{18a}.
\end{aligned}
$$
**Cauchy relation $C_{1122}\stackrel{?}{=}C_{1212}$:**
$17\alpha\stackrel{?}{=}39-22\alpha\Rightarrow\alpha\stackrel{?}{=}1$. Cauchy holds
**only at $\alpha=1$** (stress-free reference). For $\alpha<1$ the geometric
stiffness $(1-\alpha)|\Delta u|^2$ adds to $C_{1212}$ but not to $C_{1122}$,
which is the standard (and benign) signature of expanding a Hookean lattice
around a prestressed reference.

**Isotropy condition $C_{1111}-C_{1122}\stackrel{?}{=}2C_{1212}$:**
$$
C_{1111}-C_{1122}-2C_{1212}=\frac{k_0}{18a}\Big[(78-34\alpha)-17\alpha-2(39-22\alpha)\Big]
=\frac{k_0}{18a}\,\big(-7\alpha+(-7\alpha\cdot 0)\big)=-\,\frac{7\alpha\,k_0}{18\,a}.
$$
(Algebra: $78-34\alpha-17\alpha-78+44\alpha=-7\alpha$.) **Isotropy is violated**
for any $\alpha>0$. The anisotropy is **not** $O((ka)^2)$ or higher — it is
**leading-order** (zeroth in $ka$), with cubic dimensionless anisotropy index
$$
\eta_{\rm cub}\equiv\frac{C_{1111}-C_{1122}-2C_{1212}}{2C_{1212}}=-\frac{7\alpha}{2(39-22\alpha)}.
$$
At $\alpha=1$: $\eta_{\rm cub}=-7/34\approx-0.206$ (a $\sim20\%$ static
anisotropy). At $\alpha=0.2$: $\eta_{\rm cub}\approx-0.0207$.

### 3.7 Comparison with the principles claim

`principles.md` §1.1a and `paper-v4/backbone.md` point 15 assert that the
26-neighbor central-force lattice with $1/|\delta|^2$ shell weights satisfies
the cubic isotropy condition automatically, with anisotropy first appearing at
$O((ka)^4)$. **The derivation above does not support this claim.** What is
true is the textbook Cauchy relation $C_{1122}=C_{1212}$ at $\alpha=1$, which is
a consequence of central-force pair-wise interactions in a stress-free Bravais
lattice; this is **necessary but not sufficient** for cubic isotropy. The
extra condition $C_{1111}=3C_{1122}$ requires a specific tuning of shell
weights that the $1/|\delta|^2$ choice does not realise. We surface this as a
discrepancy that the project should resolve (either by retuning shell weights,
by accepting a small-but-finite static anisotropy, or by checking whether the
isotropy claim was intended only "after coarse-graining over many cells",
which would still leave the leading-order cubic correction in the elastic
tensor).

---

## 4. Result

For the BraneSim cubic lattice with 26 neighbors, shell weights
$k_\delta=k_0/|\delta|^2$, prestretch $\alpha$, lattice spacing $a$, node mass
$m$:

$$
\rho_m=\frac{m}{a^3},\qquad
\mu(\alpha)=\frac{k_0(78-61\alpha)}{18a},\qquad
\lambda(\alpha)=\frac{k_0(88\alpha-78)}{18a},\qquad
\lambda+2\mu=\frac{k_0(78-34\alpha)}{18a}.
$$

Cubic elastic tensor (in Voigt notation):
$C_{1111}=\tfrac{k_0(78-34\alpha)}{18a}$,
$C_{1122}=\tfrac{17\alpha k_0}{18a}$,
$C_{1212}=\tfrac{k_0(39-22\alpha)}{18a}$.

Wave speeds along a cubic axis (from the acoustic tensor $A_{abcd}$):
$$
c_L^2(\alpha)=\frac{k_0(78-34\alpha)}{18\,a\,\rho_m},\qquad
c_T^2(\alpha)=\frac{k_0(78-61\alpha)}{18\,a\,\rho_m}.
$$

In the unprestretched limit $\alpha\to1$:
$$
c_L^2=\frac{44\,k_0}{18\,a\rho_m}=\frac{22\,k_0}{9\,a\rho_m},\qquad
c_T^2=\frac{17\,k_0}{18\,a\rho_m},\qquad
\frac{c_L}{c_T}\Big|_{\alpha=1}=\sqrt{\tfrac{44}{17}}\approx1.609.
$$
For comparison, an isotropic medium with $\lambda=2\mu$ (Poisson ratio
$\nu=1/4$) gives $c_L/c_T=\sqrt{3}\approx1.732$; an incompressible medium
$\lambda\to\infty$ gives $c_L/c_T\to\infty$.

**Effect of $\alpha<1$.** Both speeds increase as $\alpha$ decreases, because
the prestress adds geometric stiffness ($+(1-\alpha)\,k_0\cdot\tfrac{26}{3}/(18a)$
to both $\lambda+2\mu$ and $\mu$). The transverse speed gains more than the
longitudinal speed (the $(1-\alpha)$ coefficient is the same in $A_{1111}$ and
$A_{2121}$, but $A_{2121}$ has a smaller $\alpha$-piece), so $c_L/c_T$
**decreases** monotonically from $1.609$ at $\alpha=1$ towards $1$ as
$\alpha\to0$ (where both speeds approach $\sqrt{26k_0/(54\,a\rho_m)}$).

---

## 5. Regime of validity

- **Small parameter:** $a|k|\ll1$ (long wavelength relative to lattice
  spacing). Leading neglected term is $O((ak)^4)$ in dispersion, $O(a^2\nabla^2u)$
  in the energy density.
- **Linear elasticity:** $|\nabla u|\ll1$. The geometric quartic
  $(\nabla u\cdot\nabla u)^2$ that gives Skyrme stabilization (paper-v4
  backbone point 17) is a separate sprint.
- **Static anisotropy.** The result $\eta_{\rm cub}\ne0$ found in §3.6 is
  exact within the long-wavelength expansion; it is **not** suppressed by any
  power of $ak$. Direction-dependence of $c_L,c_T$ along $\hat e_1$ vs along
  the body-diagonal $(\hat e_1+\hat e_2+\hat e_3)/\sqrt3$ is therefore
  $O(\eta_{\rm cub})$, **not** $O((ka)^2)$.
- **Prestress finiteness:** $0<\alpha\le1$. The result $\mu>0$ requires
  $\alpha<78/61\approx1.279$, which is automatic. $\lambda>0$ requires
  $\alpha>78/88\approx0.886$; for $\alpha<0.886$ the effective $\lambda$ is
  **negative**, i.e. $c_L^2<2c_T^2$, equivalently Poisson ratio $\nu<0$. The
  default config $\alpha=0.2$ is deep in the $\lambda<0$ regime; this is a
  legitimate regime for prestressed media (auxetic-like behavior coming from
  the geometric stiffness dominating over the bare Hookean stiffness) and is
  the operational regime of all current BraneSim runs.

---

## 6. Falsifiable numerical prediction

**Parameters from `orchestration/configs/local_extensive.json`:**
`spring_constant = k_0 = 1.0`, `spacing = a = 1.0`, `mass_density = ρ_m = 1.0`,
`rest_length = 0.2` so $\alpha=$ `rest_length/spacing` $=0.2$.

Plugging in:
$$
c_L^2=\frac{1\cdot(78-34\cdot0.2)}{18\cdot1\cdot1}=\frac{71.2}{18}=3.9556,\qquad
c_T^2=\frac{1\cdot(78-61\cdot0.2)}{18\cdot1\cdot1}=\frac{65.8}{18}=3.6556,
$$
$$
\boxed{\;\;c_L\approx1.989,\quad c_T\approx1.912,\quad \frac{c_L}{c_T}\approx1.040\;\;}
\qquad(\alpha=0.2,\ k_0=a=\rho_m=1).
$$
At the no-prestress reference $\alpha=1$ (same $k_0,a,\rho_m$):
$c_L=\sqrt{44/18}\approx1.563$, $c_T=\sqrt{17/18}\approx0.972$,
$c_L/c_T\approx1.609$.

For the body diagonal $k\parallel(1,1,1)$, the same acoustic tensor gives
$$
c_{L,[111]}^2
=\frac{k_0}{2a\rho_m}\left[(1-\alpha)\frac{26}{3}+\alpha\frac{146}{27}\right],
\qquad
c_{T,[111]}^2
=\frac{k_0}{2a\rho_m}\left[(1-\alpha)\frac{26}{3}+\alpha\frac{44}{27}\right].
$$
At $\alpha=0.2$ this predicts
$c_{L,[111]}\approx2.002$ and $c_{T,[111]}\approx1.905$, i.e. a small but
leading-order shift relative to the axis values ($+0.65\%$ longitudinal,
$-0.35\%$ transverse in speed).

**Measurement procedure** (executable in current `components/`):

1. Initialize a periodic 32³ lattice flat ($u\equiv0$) and superpose a small
   plane-wave seed $u_a(x,0)=\epsilon\,\hat p_a\cos(k\cdot x)$ with
   $\epsilon=10^{-3}\,a$, $k=k\hat e_1$, $|k|a\in\{0.1,0.2,0.4\}$. Run two
   polarizations: $\hat p=\hat e_1$ (longitudinal) and $\hat p=\hat e_2$
   (transverse).
2. Time-integrate for $t\in[0,T]$ with $T$ several wave periods. Record
   $u(x,t)$ at each saved step.
3. Fit $u(x,t)=\epsilon\,\hat p\cos(k\cdot x-\omega t)$ to extract $\omega$;
   compute $c=\omega/|k|$.
4. Extrapolate $c(k)\to c(0)$ by fitting $c^2(k)=c_0^2+c_2(ak)^2+\cdots$.

**Failure thresholds.**

- **Hard.** $|c_L^{\rm meas}-1.989|/1.989>5\%$ or
  $|c_T^{\rm meas}-1.912|/1.912>5\%$ at $|k|a=0.1$ → the derivation or the
  force code is wrong.
- **Anisotropy probe.** Repeat (1)–(3) with $k\parallel(\hat e_1+\hat e_2+
  \hat e_3)/\sqrt3$ and compare to the explicit $[111]$ predictions above
  after extrapolating $k\to0$. At $\alpha=0.2$ the expected effect in speed is
  sub-percent, so this probe needs either careful extrapolation or an
  additional $\alpha=1$ run, where the same bare anisotropy produces a much
  larger direction-dependent shift.
- **Cauchy probe.** Apply uniform shear $\partial_2 u_1=\gamma$ vs uniaxial
  $\partial_1 u_1=\gamma$ vs equibiaxial $\partial_1 u_1=\partial_2 u_2=\gamma$
  by setting up the corresponding initial $u$ field and reading off the
  potential energy. Compare $C_{1111},C_{1122},C_{1212}$ to §3.6. Failure of
  the predicted ratios at $\alpha=1$ by more than $1\%$ would indicate the
  symbolic algebra is wrong.

---

## 7. What remains open

1. **Shell-weight decision.** The old automatic-isotropy claim has been
   revised in the project principles, but the model still needs a decision:
   either retune the shell weights or carry the finite leading-order anisotropy
   explicitly. For effective shell weights $(w_{\rm I},w_{\rm II},w_{\rm III})$
   multiplying the axial, face-diagonal, and body-diagonal tensor sums after
   the $k_\delta|\delta|^2$ factor, cubic isotropy requires
   $2w_{\rm I}=w_{\rm II}+\tfrac{16}{9}w_{\rm III}$. The current
   $1/|\delta|^2$ implementation has $(1,1,1)$ and does not lie on this curve.
2. **Static vs acoustic moduli at $\alpha<1$.** The reported $\mu,\lambda$ are
   from the wave-equation tensor $A_{abcd}$ (asymmetric in the prestressed
   case). The static, symmetric $C_{ijkl}$ in §3.6 differs in $C_{1212}$. The
   wave-equation choice is the right one for sound speeds and dispersion
   tests, but for stress/strain analysis on slowly varying configurations
   the static $C_{ijkl}$ is the physically meaningful quantity. Both should
   coincide at $\alpha=1$.
3. **Negative $\lambda$ regime ($\alpha<0.886$).** For the default config
   $\alpha=0.2$ we get $\lambda\approx-3.36\,k_0/a<0$. This is mechanically
   stable (because $\lambda+2\mu>0$ and $\mu>0$ both hold), but it inverts
   the standard sign of cross-diagonal coupling and may have qualitative
   consequences for the gravity channel (which sources in-brane strain via
   the constitutive law). Worth checking whether other BraneSim experiments
   have implicitly assumed $\lambda>0$.
4. **Convention for the reference configuration.** I expanded around the
   **held** lattice $R_n=an+u(an)$, not the prompt's stress-free reference
   $R_n=\alpha a n+u$. Expanding around the stress-free reference would give
   a strain energy with no $(1-\alpha)$ term (because the linear-in-$\Delta u$
   stretch vanishes at the stress-free configuration only if there is no
   prestress). The prompt's choice would erase the prestress story and is
   inconsistent with the simulation's actual reference state. Surfaced under
   A2 for the reviewer.
5. **Higher-order anisotropy ($O((ka)^4)$).** The derivation here is at
   $O((ka)^0)$ in the long-wavelength limit; the claimed $O((ka)^4)$
   anisotropy of `principles.md` would presumably be the next non-trivial
   correction beyond a hypothetical leading-order isotropic theory. Since
   the leading order is **already** anisotropic in the present derivation,
   computing the $O((ka)^2)$ dispersion correction is the next natural
   sprint task, not something to defer to $O((ka)^4)$.
6. **X⁴ amplitude channel.** Excluded by A1. The full lattice→continuum map
   that includes amplitude (and hence the geometric quartic
   $(\partial_i u)^2(\partial_j u)^2$ relevant for Skyrme stabilization) is
   a separate sprint.
