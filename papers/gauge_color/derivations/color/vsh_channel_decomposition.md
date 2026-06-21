# VSH channel decomposition of the StVK + geometric quartic action at the soliton scale

**Sprint:** Soliton-hunt prerequisite (BACKBONE.md #20 verification).
**Author:** mathematical-physics agent.
**Operating point:** `α = 0.2`, `w/a ∈ {5, 10, 20}`, 6-neighbor axial-only cubic lattice.
**Question:** Is the vector-spherical-harmonic (VSH) description protocol of backbone #20
quantitatively sound at the soliton scale, given that the substrate point group is `O_h`
not `SO(3)`?

This brief stays in the descriptive language: substrate is and remains cubic. We ask only
whether the *emergent* approximate `SO(3)` at scales `≫ a` is good enough that channel
labels `(J, L, M)` separate the dynamics to the order needed for the soliton hunt.

---

## 0. Setup, conventions, master action

**Lateral triplet.** Following paper §3.7 we set
$\vecX(x,t) = (h_\star x + \boldsymbol\xi(x,t),\, u(x,t))$ with $\boldsymbol\xi\in\R^3$ the
lateral displacement and $u\in\R$ the embedding-amplitude direction. The induced metric
is
$$
g_{ij} = h_\star^2 \delta_{ij} + h_\star (\partial_i \xi_j + \partial_j \xi_i)
       + \partial_i\boldsymbol\xi\cdot\partial_j\boldsymbol\xi + \partial_i u\,\partial_j u.
$$
Define $\hat g = (g^0)^{-1} g$ with $g^0 = \ell_0^2 \delta_{ij}$. Using $\alpha = \ell_0/h_\star$,
$$
\hat g_{ij} = \alpha^{-2}\Big[\delta_{ij} + h_\star^{-1}(\partial_i\xi_j+\partial_j\xi_i)
            + h_\star^{-2}\,\partial_i\boldsymbol\xi\!\cdot\!\partial_j\boldsymbol\xi
            + h_\star^{-2}\,\partial_i u\,\partial_j u\Big].
$$
Set units $h_\star = 1$ and absorb the constant $\alpha^{-2}$ into the Lamé constants, so
the Green–Lagrange strain reads
$$
E_{ij} = \tfrac12(\partial_i\xi_j+\partial_j\xi_i)
       + \tfrac12\,\partial_i\boldsymbol\xi\!\cdot\!\partial_j\boldsymbol\xi
       + \tfrac12\,\partial_i u\,\partial_j u
       + (\text{const}).
$$
The constant offset shifts the reference and is dropped.

**StVK + geometric quartic.** The stored energy density is
$$
W = \mu\,\mathrm{tr}(E^2) + \tfrac{\lambda}{2}(\mathrm{tr} E)^2 .
$$
Expanding to quartic order in $(\partial\xi,\partial u)$ and discarding cubic cross terms
that vanish for the spherically symmetric `J=0, L=1` hedgehog (`tr` of an antisymmetric
piece is zero), the energy splits into a quadratic harmonic part `W_2` and a quartic
"Skyrme-class" piece `W_4`:
$$
W_2 = \mu\,\varepsilon_{ij}\varepsilon_{ij} + \tfrac{\lambda}{2}(\varepsilon_{kk})^2,
\qquad \varepsilon_{ij} = \tfrac12(\partial_i\xi_j+\partial_j\xi_i),
$$
$$
W_4 = \tfrac{\mu}{2}\big[\partial_i\boldsymbol\xi\!\cdot\!\partial_j\boldsymbol\xi\big]\!\big[\partial_i\boldsymbol\xi\!\cdot\!\partial_j\boldsymbol\xi\big]
    + \tfrac{\lambda}{4}\big[|\nabla\boldsymbol\xi|^2\big]^2 + \cdots
$$
(plus lateral–amplitude cross terms via `∂u`; these are part of the Skyrme-twist sector
in backbone #19 and are kept symbolically as a single $W_4^{(u)}$ piece in §1.3 below).

**Cubic anisotropy term.** From `lattice_to_continuum.md` the 6-neighbor axial-only
stencil gives a quadratic energy that, in continuum form, contains an extra cubic-anisotropy
operator
$$
W_{\rm cub} = \tfrac{1}{2}\,\eta_{\rm cub}\,\sum_a (\partial_a \xi_a)^2 \cdot a^2 \cdot \square,
$$
where the prefactor and structure are derived in §3 below from $D(k)$ along $[100]$ vs
$[111]$. The cubic correction is a $k^2$-times-leading operator (not $k^0$, in `D(k)/k^2`),
which is the crucial point for the splitting estimate.

---

## 1. VSH channel decomposition of the StVK + geometric quartic action

### 1.1 Layer interface

- **Source layer.** StVK + geometric quartic continuum on $\boldsymbol\xi(x)\in\R^3$.
- **Target layer.** Radial functions $\{f_{JLM}(r)\}$ in the VSH basis
  $\boldsymbol\xi^i(r,\theta,\varphi) = \sum_{JLM} f_{JLM}(r)\,Y^i_{JLM}(\theta,\varphi)$.
- **Projection.** $f_{JLM}(r) = \int d\Omega\,(Y^i_{JLM})^*\,\boldsymbol\xi^i$.

### 1.2 Assumptions

1. **A1 — Emergent SO(3) at quadratic order.** Cubic-anisotropy corrections to the
   isotropic continuum Lamé form scale as $(a/w)^2$ for a soliton of width $w$ (justified
   below in §3). Justifies using VSH for the *leading* term.
2. **A2 — Width ≫ a (Skyrme regime).** $w/a \gg 1$, so lattice corrections are subleading.
   (Backbone #17.)
3. **A3 — Linear elasticity at $W_2$ level.** $|\nabla\xi|\ll 1$. The quartic term
   $W_4$ contributes via Derrick balance but is small at each spatial point. **To be
   checked numerically** that the hedgehog amplitude `f(0)` satisfies this.
4. **A4 — Static problem.** We compute the static energy functional; time-dependent
   terms enter as $\rho_m\dot\xi^2/2$ and are diagonal in the VSH basis.
5. **A5 — Hedgehog channel definition.** "Canonical hedgehog" means the unique VSH
   with $J=0, L=1, M=0$ on the lateral triplet:
   $\boldsymbol\xi^i = f(r)\,\hat x^i = f(r)\,Y^i_{010}$ up to normalization.
6. **A6 — Drop $X^4$ Skyrme twist for §1, restore in §2.** For pure-lateral §1 we set
   $u\equiv 0$. The Skyrme-twist into $X^4$ enters Derrick balance and is restored in §2.

### 1.3 Quadratic-order channel structure

The quadratic part $W_2$ is built from $\varepsilon_{ij} = \tfrac12(\partial_i\xi_j+\partial_j\xi_i)$
plus the cubic correction $W_{\rm cub}$ (§3). For the isotropic part (the SO(3)-symmetric
limit of $W_2$):
$$
W_2^{\rm iso} = \mu\,\varepsilon_{ij}\varepsilon_{ij} + \tfrac\lambda 2(\varepsilon_{kk})^2.
$$

A VSH $Y^i_{JLM}$ is an eigenvector of $J^2$, $J_z$, $L^2$ (after gauge fixing)
simultaneously. The kinetic operator
$T_{ij}^{kl}\partial_l\partial_j$ derived from $W_2^{\rm iso}$ (the Christoffel matrix
$M_{ij}(k) = (\lambda+\mu)k_i k_j + \mu\,k^2\delta_{ij}$ rotated back to position space)
is rotationally invariant, so it conserves $J$ and $M$ and is **diagonal in the
$(J, M)$ blocks**.

However, **within a given $J$ block** the operator mixes the two `L = J±1` parity-even
VSH (and decouples from $L = J$, which has opposite parity). The decomposition is the
standard one for vector fields on $\R^3$:

- **For $J=0$:** Only $L=1$ exists (since $|J-L|=1$). Single channel, no mixing.
- **For $J\ge 1$:** Three channels per $(J,M)$ — $L\in\{J-1, J, J+1\}$ — with $L=J$
  decoupled (magnetic) and $L=J\pm 1$ mixed by the longitudinal term $(\lambda+\mu)\partial_i\partial_j$.
  In the Helmholtz language: the $L=J$ channel is the divergence-free transverse part
  ("magnetic" VSH), and the $L = J\pm 1$ channels span the gradient-of-scalar + radial-vector
  ("electric" + "longitudinal") sectors.

**At quadratic order:**
- Different $(J,M)$ blocks are **independent**.
- Inside $(J,M)$, the magnetic channel ($L=J$) decouples; the two electric/longitudinal
  channels ($L=J\pm 1$) mix into two normal modes per radial node.

### 1.4 Quartic interaction structure

$W_4$ is a polynomial in $\partial_i\boldsymbol\xi\!\cdot\!\partial_j\boldsymbol\xi$. In
the VSH basis this is a Clebsch–Gordan-controlled coupling between channels. Specifically
$\partial_i\boldsymbol\xi^a$ is an `L=1` tensor in spatial index $i$ and an `L=0` vector
in internal index $a$ (or vice versa under combined rotation). Schematically the
quartic vertex mixes four VSH modes obeying the standard angular-momentum addition rule:
$$
(J_1,L_1)\otimes(J_2,L_2)\otimes(J_3,L_3)\otimes(J_4,L_4) \supset (0,0)
\quad\Longleftrightarrow\quad
\sum L_i \text{ even, } J_{\rm total}=0.
$$

**Key consequence for the soliton hunt.** Working in the canonical hedgehog channel
$(J=0, L=1)$, the quartic vertex with four hedgehog legs is **non-vanishing** (since
$1\otimes 1\otimes 1\otimes 1 \supset 0$). Hence the geometric quartic provides a
self-coupling in the hedgehog channel, which is exactly what Derrick balance §2 needs.

### 1.5 Radial ODE in the canonical hedgehog channel

For the hedgehog ansatz $\boldsymbol\xi^i = f(r)\,\hat x^i$:
$$
\partial_i \xi_j = f'(r)\,\hat x_i\hat x_j + \frac{f(r)}{r}\big(\delta_{ij} - \hat x_i \hat x_j\big),
$$
$$
\mathrm{tr}(\nabla\boldsymbol\xi) = \partial_k \xi_k = f' + 2 f/r,
\qquad
\varepsilon_{ij}\varepsilon_{ij} = (f')^2 + 2(f/r)^2.
$$
Substituting into $W_2^{\rm iso}$:
$$
W_2^{\rm iso}[f] = \mu\left[(f')^2 + 2(f/r)^2\right] + \tfrac{\lambda}{2}\left(f' + 2f/r\right)^2.
$$
The radial energy is $E_2 = \int_0^\infty W_2^{\rm iso}[f]\,4\pi r^2\,dr$. The
Euler–Lagrange ODE is
$$
\boxed{
(\lambda+2\mu)\Big[f'' + \tfrac{2}{r}f' - \tfrac{2 f}{r^2}\Big] = 0 \qquad (W_2\text{-only}),
}
$$
i.e. the standard $L=1$ radial Laplacian for a vector field on $\R^3$. Without the
quartic, this admits only $f(r) = A r + B/r^2$ — no localized soliton. The quartic
term is essential.

---

## 2. Derrick balance for the $J=0, L=1$ hedgehog

### 2.1 Layer interface

- **Source.** Static energy functional $E[\xi]$ from $W_2 + W_4$.
- **Target.** Scalar function $E(\lambda)$ of the Derrick rescaling parameter $\lambda$
  applied to $f(r) \to f(r/\lambda)$. Stationary $\lambda^*$ gives the hedgehog size.

### 2.2 Assumptions

1. **B1 — Skyrme-twist coupling to $X^4$ provides the topological stabilizer.** Without
   the $X^4$ twist a pure lateral hedgehog is not topologically stable in the soliton
   sense; we treat the $J=0, L=1$ ansatz here as the angular-channel part and incorporate
   the radial $X^4$ profile $u(r)$ as a Skyrme angle $F(r)$ via
   $u(r) = u_0\,\cos F(r)$ and $f(r) = (1)\,\sin F(r)$ with $F(0)=\pi, F(\infty)=0$
   (standard Skyrme hedgehog ansatz).
2. **B2 — Isotropy approximation.** Cubic anisotropy enters as a separate small
   $(a/w)^2$ correction (§3), not as a Derrick destabilizer.
3. **B3 — StVK quartic is positive definite.** $\mu, \lambda+2\mu/3 > 0$ in the regime
   needed. For α=0.2 from the 6-neighbor stencil, $c_L^2 = 1 > 0$ and $c_T^2 = 0.8 > 0$,
   so both $\mu, \lambda+2\mu > 0$. (Note: $\lambda$ itself may be negative in the
   axial-only stencil, but $\lambda + 2\mu > 0$ is the stability condition.)

### 2.3 Derrick scaling

Let $E_n = \int (W_n[\xi]) d^3 x$ where $n$ is the number of derivatives. Under
$\xi(x) \to \xi(x/\lambda)$, in 3D:
$$
E_2 \to \lambda^{+1} E_2, \qquad E_4 \to \lambda^{-1} E_4.
$$
(Standard: count derivative powers minus spatial dimension.) Hence
$$
E(\lambda) = \lambda\,E_2 + \lambda^{-1} E_4,
$$
$$
\frac{dE}{d\lambda} = E_2 - \lambda^{-2} E_4 = 0
\quad\Longrightarrow\quad \lambda^* = \sqrt{E_4/E_2}.
$$

### 2.4 Second variation

$$
\frac{d^2 E}{d\lambda^2}\bigg|_{\lambda^*} = 2 \lambda^{*-3} E_4 > 0,
$$
so the extremum is a **minimum** in $\lambda$ (verified). This is the standard Skyrme
balance: $W_2$ wants the soliton to spread ($\lambda \to \infty$), $W_4$ wants it to
shrink ($\lambda \to 0$), and the balance is at finite $\lambda^*$.

**Topological lock.** Derrick stability in $\lambda$ is not absolute stability against
unwinding to vacuum; the Skyrme twist $F(r): [0,\infty]\to [\pi, 0]$ has winding number
$B=1$, which is conserved if $F$ does not pass through $0$ or $\pi$ inside the soliton
core. The Derrick balance gives the *size*; topology gives *existence*. Both are needed.
(Surfaced explicitly in §7.)

### 2.5 Predicted hedgehog radius

Run the Derrick balance with the **lattice-exact** densities
$W_2 = \frac{k_s(1-\alpha)}{2a}|\partial u_\perp|^2$ and
$W_4 = \frac{k_s\alpha}{8a}|\partial u_\perp|^4$ (the geometric quartic
$\propto k_s\alpha/a$; see `geometric_nonlinearity_alpha_scaling.md`). For the
Skyrme hedgehog $u_\perp = A\,\hat n$, $|\partial u_\perp|^2 = A^2[F'^2 + 2\sin^2F/r^2]$,
the Derrick energy is $E(b) = a_2 b + a_4 b^{-1}$ with
$a_2 = \frac{2\pi k_s(1-\alpha)A^2}{a}I_2$, $a_4 = \frac{\pi k_s\alpha A^4}{2a}I_4$,
giving $b_* = \sqrt{a_4/a_2} = \frac{A}{2}\sqrt{\frac{\alpha}{1-\alpha}\frac{I_4}{I_2}}$ and

$$
\boxed{\;\frac{R_h}{a} = \kappa\,\frac{A}{a}\,\sqrt{\frac{\alpha}{1-\alpha}}\;,\qquad
\kappa = \tfrac{\sigma}{2}\sqrt{I_4/I_2} = O(1).\;}
$$

**Key properties:** (i) $R_h \propto +A$ (linear in amplitude) — so a fat soliton
wants **large** $A$, not small. (ii) $R_h$ grows monotonically with α as
$\sqrt{\alpha/(1-\alpha)}$ — push α **up** toward the $\alpha\lesssim0.8$ edge
(above which $c_T^2=(1-\alpha)\to0$ kills confinement). The radius is set by
$\alpha/(1-\alpha)$; the $\alpha(1-\alpha)$ peak at 0.5 governs
binding-vs-confinement *quality*, a separate quantity.

**Falsifiable prediction** (fixed $A/a$, sweep α): $R_h(\alpha)/R_h(0.5)=\sqrt{\alpha/(1-\alpha)}$
= 0.65 / 1.00 / 1.53 / 2.00 at α = 0.3 / 0.5 / 0.7 / 0.8. A flat or *decreasing*
trend falsifies the prediction. **Operational seed target:** α=0.7, $A/a\approx10$
$\Rightarrow R_h/a\approx8$.

**Open risk:** $A/b_* = 2\sqrt{(1-\alpha)/\alpha}/\sqrt{I_4/I_2}$ is $O(1)$ and
α-set (independent of $A$) — single-scale Skyrme balance lands the soliton at
order-unity strain, the edge of the continuum premise. Whether the relaxed solution
stays continuum or pins at width $\sim a$ is undecided and is what the sweep measures.

### 2.6 Summary

- Derrick extremum exists at finite $\lambda^*$ when both $W_2$ and $W_4$ are positive.
- The extremum is a Derrick-minimum (second-derivative-positive).
- Absolute stability requires the topological winding $B=1$ from the Skyrme twist;
  this is an axiom for §2 (surfaced in §7).
- Hedgehog radius (§2.5): $R_h/a = \kappa\,(A/a)\sqrt{\alpha/(1-\alpha)}$
  — *linear* in the $X^4$ amplitude $A$ and *increasing* in α. **The width $w \gg a$
  regime is reachable at large amplitude and high α (≈0.7).**

---

## 3. Cubic-anisotropy splitting of SO(3) multiplets at soliton scale

This is the falsifiable test of the VSH protocol.

### 3.1 Layer interface

- **Source.** Closed-form $D(k)$ from `lattice_to_continuum.md` and §3 of the paper:
  on 6-neighbor axial-only, the dynamical matrix is diagonal in the Cartesian basis with
  eigenvalues that depend on direction $\hat k$.
- **Target.** Splitting $\Delta\omega$ within an SO(3) multiplet (e.g. $L=2\to E_g \oplus T_{2g}$)
  as a function of the small parameter $a/w$, where $w$ is the soliton width.

### 3.2 Closed-form $D(k)$ for 6-neighbor axial-only

From `lattice_to_continuum.md`: in units $k_{\rm axial}=a=\rho=1$, the dynamical matrix is
**diagonal** in the Cartesian basis at every $k$:
$$
D_{ab}(\mathbf k) = \delta_{ab}\,\Lambda_a(\mathbf k),
\qquad
\Lambda_a(\mathbf k) = 4 \sin^2(k_a/2) + (1-\alpha)\sum_{b\ne a} 4\sin^2(k_b/2).
$$
Expanding for small $|k|a$:
$$
\Lambda_a(\mathbf k) = k_a^2 + (1-\alpha)\sum_{b\ne a}k_b^2 \;-\; \tfrac{1}{12}\Big[k_a^4 + (1-\alpha)\sum_{b\ne a}k_b^4\Big] + O(k^6).
$$
Write $\Lambda_a = (1-\alpha)k^2 + \alpha k_a^2 - (a^2/12)\big[(1-\alpha)k^4_{\rm shell} + \alpha k_a^4\big] + \cdots$
where I have inserted lattice spacing $a$ in the dimensional sub-leading term to make
the small parameter explicit.

### 3.3 SO(3)-invariant + cubic-anisotropy split

At leading order:
$$
\Lambda_a(\mathbf k) = \underbrace{(1-\alpha)|\mathbf k|^2}_{SO(3)\text{-invariant, isotropic}}
                    + \underbrace{\alpha k_a^2}_{\text{cubic anisotropy}}
                    + O(k^4 a^2).
$$
The SO(3)-invariant piece is $(1-\alpha)k^2 \delta_{ab}$ — pure isotropic transverse-stiffness.
The cubic anisotropy is the operator $\alpha k_a^2 \delta_{ab}$, which **is not** an
$SO(3)$ scalar: it depends on direction.

Crucially the cubic anisotropy is at the **same order in $k^2$** as the isotropic
piece (not at $k^4$). What suppresses it at long wavelength is **not** a $k^2$ vs
$k^0$ hierarchy but rather the dimensionless ratio:
$$
\frac{\text{anisotropy strength}}{\text{isotropy strength}} = \frac{\alpha}{(1-\alpha)} \cdot \frac{k_a^2}{|\mathbf k|^2}.
$$
At α=0.2: $\alpha/(1-\alpha) = 0.25$, and the dispersion gap between $[100]$ and $[111]$
in $\omega^2$ is exactly $\alpha/3 \cdot 1/(1-\alpha + 2\alpha/3) \approx 0.083$, i.e.
**~8% in $\omega^2$, or ~4% in $\omega$**. This matches the ~10.6% in $c$ (L vs T along
[100]) reported in Sprint 1 to within the difference of two distinct anisotropy probes
(L–T gap at fixed $k$ vs L speed across directions).

### 3.4 Operator structure responsible for SO(3) splitting

Decompose the dispersion tensor $D_{ab}(\mathbf k)$ into SO(3)-irreducible pieces in
the internal index $a$:
$$
D_{ab}(\mathbf k) = D^{(0)}(\mathbf k)\,\delta_{ab} + D^{(2)}_{ab}(\mathbf k),
$$
with
$$
D^{(0)}(\mathbf k) = (1-\alpha)|\mathbf k|^2 + \tfrac{\alpha}{3}|\mathbf k|^2,
\qquad
D^{(2)}_{ab}(\mathbf k) = \alpha\Big(k_a^2 - \tfrac13 |\mathbf k|^2\Big)\delta_{ab}.
$$
The traceless piece $D^{(2)}$ is the entire $SO(3)$-breaking content at leading order
in $|\mathbf k|a$.

Two distinct things must be tracked:
(i) **band-structure anisotropy.** $D^{(2)}$ depends on direction $\hat k$ at order
$|\mathbf k|^2$. Plane waves at fixed $\mathbf k$ have $\omega^2_L([100]) - \omega^2_T([100]) = \alpha\,k^2$, the
~10.6% L–T gap measured in Sprint 1. **This is a property of the dispersion relation,
not of any localized mode**, and is *not* what splits the soliton's bound-state multiplets.
(ii) **multiplet splitting on a bound mode.** A soliton of width $w$ is a wavepacket of
plane waves with momentum distribution peaked at $|\mathbf k|\sim 1/w$. The splitting of
its $SO(3)$ multiplets is governed by the matrix elements of $D^{(2)}$ between
$O_h$-irreps within that multiplet.

### 3.5 Selection rule and $(a/w)^p$ power-counting

The cubic point group $O_h$ has the trivial irrep $A_{1g}$ appearing in $SO(3)$
representations at $L = 0, 4, 6, 8, \ldots$; it does *not* appear at $L=2$ (the $L=2$
multiplet of $SO(3)$ decomposes into $E_g \oplus T_{2g}$ of $O_h$, neither of which
is $A_{1g}$). Therefore the splitting of an $SO(3)$ multiplet into $O_h$-irreps requires
an angular operator with $L\ge 4$ cubic-harmonic content in $\hat k$.

The $|\mathbf k|^2$-level piece $D^{(2)}_{ab}(\mathbf k) = \alpha(k_a^2 - \tfrac13 k^2)\delta_{ab}$
contains only $L = 0$ and $L = 2$ angular content in $\hat k$ — neither contributes a
genuine splitting of an $L$-multiplet (the $L=0$ part is a uniform shift; the $L=2$
part has no diagonal cubic-invariant matrix element on an $L$-multiplet). The first
cubic-anisotropy operator that splits multiplets appears at the next order of the
Taylor expansion of $\Lambda_a$:
$$
\Lambda_a(\mathbf k) \supset -\tfrac{a^2}{12}\big[\alpha\,k_a^4 + (1-\alpha)\,k^4\big].
$$
Decompose $\sum_a k_a^4 = \tfrac35\,k^4 + K_4(\hat k)\,k^4$, where
$K_4(\hat k) \propto (\hat k_x^4 + \hat k_y^4 + \hat k_z^4 - \tfrac35)$ is the
**lowest cubic invariant** in the $L=4$ representation of $SO(3)$ — the first $A_{1g}$
of $O_h$ that lifts $SO(3)$-degeneracy. The isotropic $\tfrac35 k^4$ part shifts all
states equally; the $K_4$ part is the multiplet-splitting operator.

**Power counting.** For a localized mode with $|\mathbf k|\sim 1/w$ and $\omega^2\sim c^2 k^2$:
$$
\frac{\delta\omega^2_{\rm split}}{\omega^2} \;\sim\; \frac{\alpha\,a^2\,\langle K_4\rangle\,k^4 / 12}{c^2 k^2}
   \;\sim\; \alpha\cdot(a/w)^2,
$$
since one power of $k^2$ cancels between numerator and denominator. Hence
$$
\boxed{\;\; \frac{\Delta\omega_{\rm split}}{\omega} \;\sim\; \tfrac12\,\alpha\cdot(a/w)^2. \;\;}
$$
The exponent is **$p = 2$**, exactly at the protocol's validity threshold.

### 3.6 Result

$$
\boxed{
\frac{\Delta\omega_{\rm split}}{\omega_{\rm gap}} \;\approx\; C\cdot\alpha\cdot\Big(\frac{a}{w}\Big)^2,
}
$$
with $C$ a Clebsch–Gordan and cubic-harmonic prefactor of order unity (estimated $C\sim 0.05-0.2$
based on the $K_4$ normalization $\langle L=2|K_4|L=2\rangle$ matrix element table).
The exponent $p = 2 \ge 2$ — **at the threshold of the protocol's validity criterion**.

### 3.7 Numerical estimates at α=0.2

Using $\omega_{\rm gap}/\omega \approx 0.106$ (the L–T gap along $[100]$ at α=0.2,
backbone #15) as the inter-band scale, and taking $C \approx 0.1$ as a conservative
order-of-magnitude estimate:

| `w/a` | $(a/w)^2$ | $C\alpha(a/w)^2$ | $\Delta\omega/\omega_{\rm gap}$ | Verdict |
|:-:|:-:|:-:|:-:|:-:|
| 5  | 0.040 | 8.0e-4 | 7.5e-3 | OK |
| 10 | 0.010 | 2.0e-4 | 1.9e-3 | comfortable |
| 20 | 0.0025 | 5.0e-5 | 4.7e-4 | very comfortable |

The criterion $\Delta\omega/\omega_{\rm gap} < 0.1$ is **easily satisfied at $w/a = 10$**:
the predicted ratio is $\sim 2 \times 10^{-3}$, two orders of magnitude below the criterion.

**Even at $w/a = 5$** (smallest reasonable soliton width considered) the ratio is below 1%.

### 3.8 Sanity check via the L–T gap itself

The L–T gap at $[100]$ in $\omega^2$ is exactly $\alpha = 0.2$ (it equals $c_L^2 - c_T^2$
in our units). This is at $k\to 0$, so it is *not* an $(a/w)^p$-suppressed effect — it
is a leading-order anisotropy. **The L–T splitting is therefore at the band-structure
level, not a multiplet-splitting effect within a single soliton mode.**

The relevant question for VSH is not "does the lattice break SO(3)?" (it does, at
leading order, ~10%) but "do the $O_h$ corrections within a single soliton's bound-state
spectrum split SO(3) multiplets significantly?". The answer derived in §3.5–3.6 is no:
**$\Delta\omega_{\rm split}/\omega \sim \alpha(a/w)^2$**, which is small for $w/a \gtrsim 10$.

The L–T gap *bounds* the *carrier* frequency separation between modes that are in
different polarization branches; it does not split a single VSH multiplet.

---

## 4. Numerical estimate at α=0.2, w/a ∈ {5, 10, 20}

### 4.1 Per the table in §3.7

$\Delta\omega/\omega_{\rm gap} < 0.1$ holds at $w/a = 10$ by a factor of ~50.

### 4.2 Direct prediction

Taking $\omega_{\rm gap}$ as the lab L–T gap (~10.6% of carrier frequency), the
predicted absolute splitting at $w/a = 10$ is
$$
\Delta\omega \approx 0.106 \cdot 1.9\!\times\!10^{-3}\cdot\omega \;\approx\; 2\!\times\!10^{-4}\,\omega.
$$
For a soliton of natural mass scale $M_h$ and frequency $\omega \sim 1/R_h$, this is
a splitting two parts in $10^4$ of the soliton mass. **Numerically tiny — below the
resolution of typical Sprint 1 dispersion runs.**

### 4.3 Falsifiable test

**Measurement procedure (deferred to soliton-hunt simulations):**

1. Build a hedgehog seed $\boldsymbol\xi^i = f(r)\hat x^i$ with $f(r)$ a fitted Skyrme
   profile at width $w/a = 10$.
2. Evolve under the substrate dynamics for $T = 100\,R_h/c$.
3. Diagnose the angular structure of the localized mode via VSH projection
   $f_{JLM}(r,t) = \int Y^{i*}_{JLM}\,\boldsymbol\xi^i\,d\Omega$.
4. Measure cross-channel leakage: compute
   $\langle f_{010}|f_{210}\rangle$ (hedgehog into $L=2$) and other inter-multiplet
   overlaps as a function of time.
5. **Pass condition:** Inter-multiplet leakage from $J=0$ to $J\ne 0$ stays
   $< 0.01$ of the hedgehog norm over $T = 100\,R_h/c$.

**Failure threshold:** If leakage $> 0.1$ (10%), the cubic-anisotropy cross-channel
coupling is stronger than predicted and the VSH protocol is unusable in this regime.

---

## 5. (Stretch) Rotational zero modes of the hedgehog

### 5.1 Layer interface

Static hedgehog $\boldsymbol\xi_0^i = f(r)\hat x^i$ has zero modes from rigid rotation
of the internal frame (which equals spatial rotation in the hedgehog gauge).

### 5.2 Moment of inertia

Promote the hedgehog to a slowly-rotating configuration $\boldsymbol\xi^i(\mathbf x, t)
= R^i_j(t)\,f(r)\hat x^j$ with $R\in SO(3)$ a slowly varying rotation. The kinetic
energy comes from $\rho_m\,\dot\xi^2/2$:
$$
T = \tfrac12 \rho_m \int (\dot R^i_j \hat x^j f(r))^2 d^3x = \tfrac12 I_{\rm ab}\,\omega_a\omega_b
$$
where $\omega_a$ is the angular velocity vector and
$$
I_{ab} = \rho_m \int f^2(r) (\delta_{ab}\,r^2 - x_a x_b)\cdot \tfrac{1}{r^2}\,d^3x
     = \tfrac{4\pi\rho_m}{3}\,\delta_{ab}\,\int_0^\infty f^2(r)\,r^2 dr.
$$
**Scaling:** $I \sim \rho_m\,R_h^3\,f_0^2$, where $f_0$ is the typical hedgehog amplitude.

### 5.3 Rotor quantization

$$
E_J = \frac{J(J+1)}{2I}.
$$
The static mass is $M_h \sim W_2 \cdot R_h \sim \mu f_0^2\,R_h$ (since $E_2 \sim
\int \mu (\partial \xi)^2 d^3x \sim \mu (f_0/R_h)^2 \cdot R_h^3$).

Therefore:
$$
\frac{E_J}{M_h c^2} \sim \frac{J(J+1)}{2 I M_h c^2}
                \sim \frac{J(J+1)}{2\,(\rho_m R_h^3 f_0^2)(\mu f_0^2 R_h / c^2)}.
$$
Using $c^2 = \mu/\rho_m$ and $f_0 \sim 1$ (dimensionless angle):
$$
\frac{E_J}{M_h c^2} \sim \frac{J(J+1)}{f_0^2 (R_h/a)^4} \cdot \mathcal O(1).
$$
For $R_h/a = 10$ and $f_0 \sim 1$:
$$
E_{J=1/2}/M_h c^2 \sim 3/40000 \sim 10^{-4}.
$$
**The rotational band sits ~$10^{-4}$ below the static mass — well separated from radial
excitations** (which are at $\sim c/R_h$, i.e. comparable to the carrier frequency, so a
fraction $\sim 1$ of $M_h c^2$).

This is the QCD-like rotational baryon spectrum regime: rotational level spacing
$\ll$ radial-excitation scale $\ll$ band-gap scale. Cleanly separated.

### 5.4 Caveat

The full quantization of the soliton rotor requires accounting for the $SU(2)$
spin–isospin combined rotation in the Skyrme model. For our purely lateral hedgehog
without spinorial $X^4$-twist content, this is a *bosonic* rotor with $J = 0, 1, 2, \ldots$.
Half-integer spin requires the spinorial holonomy mechanism alluded to in principles
§4.4, which is *not* established by Derrick balance alone.

---

## 6. Regime of validity and go/no-go verdict

### 6.1 Validity envelope

The VSH description protocol is valid when:
- $w \gg a$ (Skyrme regime, backbone #17). For α=0.2 with the $X^4$-twist mechanism,
  reachable at $u_0/\ell_0 \approx 0.03$.
- The multiplet splitting $\Delta\omega_{\rm split}/\omega_{\rm gap} < 0.1$. By §3.6 this
  holds at $w/a \ge \sqrt{\alpha C \cdot 10/\omega_{\rm gap}} \approx 3$, i.e. *all*
  widths in $\{5, 10, 20\}$ are safe.
- The amplitude $f_0$ stays in the linear-elasticity range. To be checked numerically
  (A3).

### 6.2 Verdict: **GO** for the partial-wave description protocol at α=0.2,
**$w/a \in [5, 20]$**.

The cubic anisotropy in the splitting kernel enters at $(a/w)^2$ — at the boundary of
the protocol's needed scaling — and at α=0.2 the prefactor $\alpha C \sim 0.02$ keeps
the actual splitting two orders of magnitude below the L–T gap at $w/a = 10$.

**What remains to be verified at the simulation layer:**
- The hedgehog's Derrick-stable width *can* be tuned to $w/a \approx 10$ by choosing
  $u_0/\ell_0 \approx 0.03$ (deferred to baryon roadmap Phase 2).
- The cross-channel leakage from $J=0$ hedgehog into higher-$J$ channels stays below
  $1\%$ over 100 light-crossing times (concrete numerical test §4.3).

**If GO fails empirically:** the most plausible failure mode is *not* the cubic anisotropy
(which is too small) but the failure of the Derrick balance at width $w/a \approx 10$
because the geometric-quartic prefactor turns out to be sub-leading to other
constitutive contributions. In that case the right reformulation is *not* to abandon
VSH — the angular labels remain meaningful — but to add a microphysical mechanism for
the soliton stabilization (e.g. coupling to a deeper amplitude channel or adding a
$|\nabla u|^4$ term from a richer hyperelastic law than StVK).

---

## 7. What remains open

1. **AXIOM (B1).** Skyrme-twist topological stabilization is *assumed*, not derived.
   The $X^4$ profile $u(r) = u_0\cos F(r)$ with $F: [0,\infty]\to[\pi,0]$ has winding
   $B=1$, but the existence of a stable static solution of the full StVK+geometric-quartic
   equations of motion with this boundary condition is not proven here. **Deferred to
   simulation (baryon roadmap Phase 2).**

2. **AXIOM (§1.5 cubic-anisotropy prefactor $C$).** The dimensionless coefficient $C$
   in the splitting formula $\Delta\omega \sim C\alpha(a/w)^2 \omega$ is estimated as
   $C \approx 0.05$–$0.2$ from cubic-harmonic matrix-element tables (Stevens
   coefficients in crystal-field theory). A precise value would require the matrix
   element $\langle L=2, m'| K_4(\hat k)|L=2,m\rangle$ in the specific basis on the
   soliton's wavefunction. *Surfaced as: if the true $C$ is ~1, then the $w/a = 5$
   row of the §3.7 table flips to $\Delta\omega/\omega_{\rm gap} \sim 0.08$ — still under
   the criterion but no longer comfortable.*

3. **DERIVED:** the existence of finite Derrick extremum $\lambda^*$, the channel structure
   of $W_2$ (block-diagonal in $J$, mixing $L=J\pm 1$), the $(a/w)^2$ scaling of the
   multiplet splitting from the closed-form $D(k)$ Taylor expansion, the rotor inertia
   scaling, and the rotational-band placement below $M_h c^2$.

4. **CONSTITUTIVE-LAW DEPENDENCE.** The Derrick balance result (existence of a finite
   $\lambda^*$, ratio $\sim \sqrt{E_4/E_2}$) survives unchanged for any isotropic
   hyperelastic energy with a $(\nabla\xi)^4$ leading term, because both StVK and
   neo-Hookean and Mooney–Rivlin all reduce to the same quadratic-plus-quartic form in
   the small-strain expansion. The numerical prefactors change but the scaling $p=2$
   in §3.6 does not, since it depends only on the structure of $D(k)$ from the lattice,
   not on the constitutive law. **Robust result.**

5. **OPEN:** Half-integer spin assignment of rotational states. The bosonic rotor §5.2
   gives $J = 0, 1, 2, \ldots$. Promoting to fermionic $J = \tfrac12, \tfrac32, \ldots$
   needs a spinorial holonomy mechanism (principles §4.4) that is independent of the
   VSH protocol established here.

6. **OPEN:** Inclusion of the cross terms $W_4^{(u)}$ (lateral × amplitude). The Skyrme
   twist in $X^4$ contributes a quartic term $\sim (\nabla u\cdot\nabla \xi)^2$ that
   was treated only schematically in §1.3. A full derivation should include these
   cross channels — they are the load-bearing piece of the Skyrme stabilization in
   backbone #19.

7. **CRITIQUE LINK.** This derivation closes the gap "is the partial-wave / VSH
   description protocol actually compatible with the cubic substrate at the soliton
   scale?" by showing $\Delta\omega/\omega_{\rm gap}$ is small at $w/a = 10$. It does
   *not* close the gap "does a stable Skyrme-twisted hedgehog *exist* in the substrate
   dynamics?" — that is the simulation question of Sprint 4.

---

## 8. Falsifiable numerical predictions (consolidated)

| Test | Prediction | Failure threshold |
|------|------------|-------------------|
| Hedgehog Derrick balance at α=0.2 (§2.5) | $\lambda^*$ finite, $R_h/a = \kappa(A/a)\sqrt{\alpha/(1-\alpha)}$; increases with $\alpha$, linear in amplitude $A$ | No finite minimum, or flat/decreasing $R_h(\alpha)$ trend |
| Cubic-anisotropy multiplet splitting at $w/a=10$, α=0.2 | $\Delta\omega/\omega_{\rm gap} \approx 2\!\times\!10^{-3}$ | $\Delta\omega/\omega_{\rm gap} > 0.1$ |
| Cross-channel VSH leakage over $100 R_h/c$ | $< 0.01$ of hedgehog norm | $> 0.1$ |
| Rotational band E_{J=1/2} | $\sim 10^{-4} M_h c^2$ at $R_h/a = 10$ | Rotational level overlaps radial excitations |

**Bottom line.** At α=0.2 with width $w/a \ge 5$, the VSH description protocol is
mathematically sound to within the precision needed for the soliton hunt. **GO.**
