Below is the **Lagrangian (action) formulation** of the same brane model you just pinned down — **discrete first**, then **continuous** — with **pre-tension handled exactly as before**: it is a **stress state in material directions**, maintained by **boundary/cell constraints**, and **never** by singling out (X^4).

---

# 3) Discrete pre-tensioned brane: Lagrangian / action

## 3.1 Configuration space

* Graph (G=(\mathcal V,\mathcal E)), typically a regular (d)-dimensional grid (often (d=3)).
* Node positions (\mathbf R_p(t)\in\mathbb R^4), (p\in\mathcal V).

For each undirected edge ((p,q)\in\mathcal E),
[
\mathbf d_{pq} := \mathbf R_q-\mathbf R_p,\qquad
\ell_{pq}:= |\mathbf d_{pq}|_{\mathbb R^4}.
]

## 3.2 Discrete Lagrangian

Kinetic energy:
[
T(\dot{\mathbf R})=\sum_{p\in\mathcal V}\frac{m}{2}|\dot{\mathbf R}_p|^2.
]

Stretch energy (Hooke springs, “no bending”):
[
U(\mathbf R)=\sum_{(p,q)\in\mathcal E}\frac{k}{2}\big(\ell_{pq}-\ell_0\big)^2.
]

Optionally add an external potential (V_{\text{ext}}(\mathbf R,t)=\sum_p V_p(\mathbf R_p,t)).

Then the **Lagrangian** is
[
\boxed{
L(\mathbf R,\dot{\mathbf R},t)=T(\dot{\mathbf R})-U(\mathbf R)-V_{\text{ext}}(\mathbf R,t).
}
]

The **action** is
[
\boxed{
S[\mathbf R]=\int_{t_0}^{t_1} L(\mathbf R,\dot{\mathbf R},t),dt.
}
]

## 3.3 Euler–Lagrange equations (unconstrained form)

For “free” interior nodes (i.e. nodes not fixed by boundary conditions),
[
\frac{d}{dt}\Big(\frac{\partial L}{\partial \dot{\mathbf R}_p}\Big)-\frac{\partial L}{\partial \mathbf R_p}=0
\quad\Rightarrow\quad
m,\ddot{\mathbf R}_p
====================

\sum_{q:(p,q)\in\mathcal E}
k(\ell_{pq}-\ell_0)\frac{\mathbf R_q-\mathbf R_p}{\ell_{pq}}
-\frac{\partial V_p}{\partial \mathbf R_p}.
]
This is exactly your force law, now derived variationally.

## 3.4 How pre-tension is imposed in the discrete action

**Pre-tension is not an extra term in (L).** It is a choice of **admissible configurations** (constraints / BCs) such that the system is held at a base configuration (\mathbf R^\star) with (\ell^\star_{pq}\neq \ell_0).

There are three clean, standard ways to encode this variationally:

### (A) Dirichlet (clamped) boundary nodes

Pick a boundary subset (\partial\mathcal V\subset\mathcal V) and prescribe
[
\mathbf R_p(t)=\mathbf R_p^\star \quad (p\in\partial\mathcal V).
]
Then the variational principle (\delta S=0) is taken over variations (\delta\mathbf R_p) that vanish on (\partial\mathcal V).
No explicit constraint forces appear in the equations; they are implicit.

### (B) Periodic cell with fixed cell shape (uniform pre-stretch)

On a regular lattice, impose periodic identifications and keep the **cell deformation** fixed so that the affine base embedding
[
\mathbf R_n^\star = A (h n)+\mathbf b
]
is compatible with the periodic cell. This prevents relaxation and yields a uniform pre-stress state.

### (C) Holonomic constraints + Lagrange multipliers

Write constraints (C_\alpha(\mathbf R,t)=0) (e.g. fixed boundary positions, fixed average deformation gradient, fixed cell vectors, etc.). Use the augmented Lagrangian
[
\boxed{
L_c = L + \sum_\alpha \lambda_\alpha(t), C_\alpha(\mathbf R,t)
}
]
and vary w.r.t. (\mathbf R) and (\lambda_\alpha). Then the “boundary/constraint forces” are (\sum_\alpha \lambda_\alpha \nabla_{\mathbf R_p} C_\alpha).

**None of these distinguishes (X^4):** constraints are about the brane’s **material geometry** (the admissible (\mathbf R)), while (\mathbf R_p\in\mathbb R^4) is always treated as a 4-vector.

## 3.5 Discrete dictionary (Lagrangian objects)

* (S[\mathbf R]): action functional.
* (L(\mathbf R,\dot{\mathbf R},t)): Lagrangian.
* (T): kinetic energy.
* (U): spring stretch energy.
* (V_{\text{ext}}): external potential (optional).
* (C_\alpha), (\lambda_\alpha): constraints and Lagrange multipliers (optional).
* Canonical momentum (if/when you later go Hamiltonian):
  [
  \mathbf P_p := \frac{\partial L}{\partial \dot{\mathbf R}_p}=m\dot{\mathbf R}_p\in\mathbb R^4.
  ]

---

# 4) Continuous pre-tensioned brane: Lagrangian / action

## 4.1 Configuration field and induced metric

Material domain (\Omega\subset\mathbb R^d), coordinates (x=(x^1,\dots,x^d)).
Embedding field
[
\mathbf X:\Omega\times\mathbb R\to\mathbb R^4,\qquad \mathbf X(x,t)=(X^1,\dots,X^4).
]
Tangent vectors (\partial_i\mathbf X\in\mathbb R^4), induced metric
[
g_{ij}=\partial_i\mathbf X\cdot\partial_j\mathbf X.
]
Reference (stress-free) metric (g^0_{ij}) (e.g. (g^0_{ij}=\ell_0^2\delta_{ij})).

## 4.2 Lagrangian density and action

Kinetic energy density:
[
\mathcal T = \frac{\rho_m}{2}|\partial_t\mathbf X|^2.
]

Strain energy density (your stretch-only choice):
[
W(g;g^0)=\frac{K}{2}\sum_{a=1}^d(\sigma_a-1)^2,
\qquad \sigma_a=\sqrt{\lambda_a},\quad
\lambda_a=\text{eigs}\big(\widehat g\big),\ \widehat g=(g^0)^{-1}g.
]

Optional external potential density (\mathcal V_{\text{ext}}(\mathbf X,t)).

Then the **Lagrangian density** and **action** are:
[
\boxed{
\mathcal L(\mathbf X,\partial\mathbf X,t)
=========================================

## \frac{\rho_m}{2}|\partial_t\mathbf X|^2

## W(g;g^0)

\mathcal V_{\text{ext}}(\mathbf X,t),
}
]
[
\boxed{
S[\mathbf X]=\int_{t_0}^{t_1}\int_{\Omega}\mathcal L,dx,dt.
}
]

This treats all embedding components (X^A) identically.

## 4.3 Euler–Lagrange PDE

Vary (S) w.r.t. (X^A). The field Euler–Lagrange equations are
[
\frac{\partial}{\partial t}\Big(\frac{\partial\mathcal L}{\partial(\partial_t X^A)}\Big)
+
\partial_i\Big(\frac{\partial\mathcal L}{\partial(\partial_i X^A)}\Big)
-----------------------------------------------------------------------

\frac{\partial\mathcal L}{\partial X^A}
=0,
]
i.e.
[
\boxed{
\rho_m,\partial_{tt}X^A
=======================

## \partial_i P_{Ai}

\frac{\partial\mathcal V_{\text{ext}}}{\partial X^A},
\qquad A=1,\dots,4,
}
]
with the (first Piola–type) stress
[
\boxed{
P_{Ai}:=\frac{\partial W}{\partial(\partial_i X^A)}.
}
]

Using the chain rule through (g_{ij}),
[
\boxed{
P_{Ai}=2,\frac{\partial W}{\partial g_{ij}},\partial_j X^A.
}
]

## 4.4 Explicit (\partial W/\partial g_{ij}) for your eigen-stretch form

Let (\widehat g=(g^0)^{-1}g). Because (g^0) is symmetric positive definite, (\widehat g) is self-adjoint w.r.t. the (g^0)-inner product. Choose eigenpairs
[
\widehat g,e_a=\lambda_a e_a,\qquad \langle e_a,e_b\rangle_{g^0}=\delta_{ab}.
]
Define (\sigma_a=\sqrt{\lambda_a}). Then
[
\boxed{
\frac{\partial W}{\partial g_{ij}}
==================================

\sum_{a=1}^d
\frac{K}{4},\frac{\sigma_a-1}{\sigma_a},e_a^{,i}e_a^{,j}.
}
]
So the Piola-type stress is
[
\boxed{
P_{Ai}
======

\sum_{a=1}^d
\frac{K}{2},\frac{\sigma_a-1}{\sigma_a},
e_a^{,i}\Big(e_a^{,j}\partial_j X^A\Big).
}
]
This is the cleanest “closed form” without changing your chosen (W). (And again: no (X^4) exception anywhere.)

## 4.5 Boundary conditions in the variational principle

The usual boundary term from integration by parts is
[
\int_{\partial\Omega} P_{Ai},n_i,\delta X^A,dS.
]
So you get standard choices:

* **Dirichlet/clamped:** (\delta \mathbf X=0) on (\partial\Omega) (i.e. (\mathbf X) prescribed).
* **Neumann/traction:** prescribe (P_{Ai}n_i=\tau_A) on (\partial\Omega).
* **Periodic:** boundary contributions cancel pairwise.

## 4.6 How pre-tension appears in the continuous Lagrangian picture

Exactly as in your clean model statement:

* Pick a time-independent **equilibrium** embedding (\mathbf X_\star) satisfying the Euler–Lagrange equations with the chosen BCs.
* Ensure it is **not stress-free**: (g_\star\neq g^0) (equivalently (\sigma_a\neq 1)).
* Maintain it via **Dirichlet**, **periodic cell**, or **tractions**.

A common uniform pre-tension base is affine:
[
\mathbf X_\star(x)=Ax+\mathbf b,\qquad g_{\star}=A^\mathsf T A=\text{const}.
]
Then (P_{Ai}(\mathbf X_\star)) is constant and (\partial_i P_{Ai}=0) holds automatically, while the BCs/cell constraints prevent relaxation back to (g=g^0).

---

# 5) Continuous dictionary (Lagrangian objects)

* (\mathcal L): Lagrangian density.
* (S[\mathbf X]): action.
* (\rho_m): mass density in material coordinates.
* (W(g;g^0)): hyperelastic stretch energy density.
* (g_{ij}=\partial_i\mathbf X\cdot\partial_j\mathbf X): induced metric.
* (g^0_{ij}): reference (stress-free) metric.
* (\widehat g=(g^0)^{-1}g): relative mixed metric.
* (\lambda_a,\sigma_a): eigenvalues and principal stretches of (\widehat g).
* (P_{Ai}=\partial W/\partial(\partial_iX^A)): Piola-type stress (material index (i), embedding index (A)).
* Canonical momentum density:
  [
  \Pi_A:=\frac{\partial \mathcal L}{\partial(\partial_t X^A)}=\rho_m,\partial_t X^A.
  ]

---

## Where we can go next (still before Hamiltonian)

If you want the next step to be fully “paper-ready”, the usual sequence is:

1. Write (S[\mathbf X]) and the boundary conditions you actually use in the simulation (periodic cell vs clamped).
2. Choose whether you keep your **eigen-stretch** (W) *or* switch to an equivalent **tensor-quadratic** (W(g-g^0)) (much simpler algebra; same qualitative role).
3. Do the **second variation** of (S) around (\mathbf X_\star) to get the linearized operator (this is where the pre-stress enters the wave dynamics in a very explicit way).

If you say “use eigen-stretch (W) as written” vs “use quadratic (W=\frac{\mu}{8}|g-g^0|_{g0^{-1}}^2)”, I’ll write the linearized PDE in the cleanest closed form for that choice.



Perfect — with your choices

1. **Periodic cell (torus)**
2. **Keep the eigen-stretch (W)**
3. **Yes: do the second variation / linearized operator (paper-ready)**

here is the clean, insertable formulation (continuous first, then a short discrete periodic-cell note).

---

## Continuous model on a periodic cell (paper-ready)

### Domain and periodicity

Let the material domain be a (d)-torus,
[
\Omega=\mathbb T^d \equiv \prod_{i=1}^d [0,L_i] \quad\text{with endpoints identified.}
]
All “test variations” and perturbations are periodic in (x).

### Kinematics

Configuration:
[
\mathbf X:\Omega\times\mathbb R\to\mathbb R^4,\qquad \mathbf X(x,t)=(X^1,\dots,X^4).
]
Induced metric:
[
g_{ij}(\mathbf X)=\partial_i\mathbf X\cdot \partial_j\mathbf X.
]
Reference (stress-free) metric (g^0_{ij}) (often constant; e.g. (g^0_{ij}=\ell_0^2\delta_{ij})).

Define the relative mixed metric
[
\widehat g := (g^0)^{-1}g,
]
with eigenvalues (\lambda_a>0) and stretches (\sigma_a=\sqrt{\lambda_a}), (a=1,\dots,d).

### Constitutive law (kept as requested)

[
\boxed{
W(g;g^0)=\frac{K}{2}\sum_{a=1}^d(\sigma_a-1)^2,
\qquad \sigma_a=\sqrt{\lambda_a(\widehat g)}.
}
]

### Action (periodic ⇒ no boundary term)

[
\boxed{
S[\mathbf X]=\int_{t_0}^{t_1}\int_{\Omega}
\left(
\frac{\rho_m}{2}|\partial_t\mathbf X|^2
---------------------------------------

W(g(\mathbf X);g^0)
\right),dx,dt.
}
]

Because (\Omega) is periodic, integrations by parts generate **no boundary contributions**.

---

## Pre-tensioned periodic base state

To represent uniform pre-tension while staying periodic, use the standard decomposition:
[
\boxed{
\mathbf X(x,t)=\mathbf X_\star(x)+\boldsymbol\xi(x,t),
\qquad
\boldsymbol\xi(\cdot,t)\ \text{is periodic on }\Omega.
}
]

Choose an affine base embedding
[
\boxed{
\mathbf X_\star(x)=A x+\mathbf b,\qquad A\in\mathbb R^{4\times d},\ \mathrm{rank}(A)=d,
}
]
so the base metric is constant:
[
g_{\star,ij}=\partial_i\mathbf X_\star\cdot\partial_j\mathbf X_\star = (A^\mathsf T A)_{ij}.
]

Pre-tension means (g_\star\neq g^0) (equivalently (\sigma_{\star,a}\neq 1) for at least one (a)).
It is “maintained” because we **do not allow** the affine part (A) (the cell shape / macroscopic deformation) to relax; only the periodic displacement (\boldsymbol\xi) evolves.

Notation that will be useful below:
[
A_i := \partial_i\mathbf X_\star \in \mathbb R^4
\quad\text{(the (i)-th tangent vector of the base, constant).}
]

---

## First variation (Euler–Lagrange PDE)

Define
[
S^{ij}:=\frac{\partial W}{\partial g_{ij}},\qquad
P_{Ai}:=\frac{\partial W}{\partial(\partial_i X^A)}=2,S^{ij},\partial_j X^A.
]
Then the field equations are
[
\boxed{
\rho_m,\partial_{tt}X^A = \partial_i P_{Ai},\qquad A=1,\dots,4.
}
]

---

# Second variation about (\mathbf X_\star): quadratic action and linearized operator

This is the “pre-stress enters explicitly” step.

### Metric variations induced by (\boldsymbol\xi)

Expand (g(\mathbf X_\star+\boldsymbol\xi)) in (\boldsymbol\xi).

* First variation:
  [
  \boxed{
  \delta g_{ij}
  =
  A_i\cdot\partial_j\boldsymbol\xi

-

A_j\cdot\partial_i\boldsymbol\xi.
}
]

* Second variation:
  [
  \boxed{
  \delta^2 g_{ij}
  =
  \partial_i\boldsymbol\xi\cdot \partial_j\boldsymbol\xi.
  }
  ]

### Define the tangent moduli at the base

Let
[
\boxed{
H_\star^{ij,kl}:=
\left.\frac{\partial^2 W}{\partial g_{ij},\partial g_{kl}}\right|*{g=g*\star}.
}
]
This is a 4th-order tensor acting on symmetric ((\delta g_{ij})).

### Quadratic (second-variation) action

Because (\mathbf X_\star) is an equilibrium on the periodic domain, the linear term drops out, and the action expanded to second order is
[
\boxed{
S^{(2)}[\boldsymbol\xi]
=======================

\int_{t_0}^{t_1}\int_{\Omega}
\left[
\frac{\rho_m}{2}|\partial_t\boldsymbol\xi|^2
--------------------------------------------

\Big(
S_\star^{ij},\partial_i\boldsymbol\xi\cdot\partial_j\boldsymbol\xi
+
\frac12,H_\star^{ij,kl},(\delta g_{ij})(\delta g_{kl})
\Big)
\right]dx,dt.
}
]

**Interpretation (important):**

* The term (S_\star^{ij},\partial_i\boldsymbol\xi\cdot\partial_j\boldsymbol\xi) is the classic **geometric stiffness** term: it is *present only because the base is pre-stressed* ((S_\star\neq 0)).
* The term with (H_\star) is the **material tangent stiffness** around the base.

### Linearized PDE from (S^{(2)})

Taking (\delta S^{(2)}=0) gives the linear operator:
[
\boxed{
\rho_m,\partial_{tt}\xi^A
=========================

\partial_i\Big(
2,S_\star^{ij},\partial_j\xi^A
+
2,H_\star^{ij,kl},(\delta g_{kl}),A_j^A
\Big),
\qquad A=1,\dots,4,
}
]
where
[
\delta g_{kl}=A_k\cdot\partial_l\boldsymbol\xi + A_l\cdot\partial_k\boldsymbol\xi.
]

Since (A_j^A), (S_\star^{ij}), (H_\star^{ij,kl}) are constant for an affine base, you can also write it as an explicit constant-coefficient second-order system:
[
\boxed{
\rho_m,\partial_{tt}\xi^A
=========================

2,S_\star^{ij},\partial_i\partial_j\xi^A
+
2,H_\star^{ij,kl},A_j^A,
\big(A_k^B,\partial_i\partial_l\xi^B + A_l^B,\partial_i\partial_k\xi^B\big).
}
]
This form makes the **embedding-component coupling** ((A\leftrightarrow B)) completely explicit — without privileging (X^4).

---

## Spectral formulas for (S_\star) and (H_\star) for the eigen-stretch (W)

Let (\widehat g_\star=(g^0)^{-1}g_\star) have eigenpairs ((\lambda_{\star,a},\Pi_a)) (projectors (\Pi_a) in material index space; use (g^0) to define orthonormality). Define (\sigma_{\star,a}=\sqrt{\lambda_{\star,a}}).

### First derivative (base stress in material indices)

[
\boxed{
S_\star^{ij}
============

\sum_{a=1}^d
\frac{K}{4},\frac{\sigma_{\star,a}-1}{\sigma_{\star,a}},
\big(\Pi_a\big)^{ij}_{(g^0)}.
}
]

### Second derivative (tangent moduli)

Let (f(\lambda)=\frac{K}{2}(\sqrt{\lambda}-1)^2). Then
[
f'(\lambda)=\frac{K}{2}\frac{\sqrt{\lambda}-1}{\sqrt{\lambda}},\qquad
f''(\lambda)=\frac{K}{4}\lambda^{-3/2}.
]
For a symmetric perturbation (\delta \widehat g=(g^0)^{-1}\delta g), the bilinear form of the second derivative is the standard spectral expression:
[
\boxed{
D^2W_\star[\delta \widehat g,\delta \widehat g]
===============================================

\sum_a f''(\lambda_{\star,a})\big(\mathrm{tr}(\Pi_a,\delta \widehat g)\big)^2
+
2\sum_{a<b}
\frac{f'(\lambda_{\star,a})-f'(\lambda_{\star,b})}{\lambda_{\star,a}-\lambda_{\star,b}},
\mathrm{tr}\big(\Pi_a,\delta \widehat g,\Pi_b,\delta \widehat g\big).
}
]
From this you read off (H_\star^{ij,kl}) (it is the unique 4th-order tensor such that
(\tfrac12 H_\star^{ij,kl}\delta g_{ij}\delta g_{kl}= \tfrac12 D^2W_\star[(g^0)^{-1}\delta g,(g^0)^{-1}\delta g])).

**Note on degeneracy:** if (\lambda_{\star,a}=\lambda_{\star,b}), the “divided difference” is interpreted in the usual continuous limit and the expression remains well-defined when written with projectors onto the degenerate eigenspace (no (X^4) issue; it’s purely material-space linear algebra).

---

# Discrete periodic cell (short, matching statement)

On a regular (d)-grid with indices (n\in\prod_i{0,\dots,N_i-1}), periodicity means indices wrap mod (N_i). For uniform pre-tension you again use
[
\mathbf R_n(t)=\mathbf R^\star_n+\mathbf u_n(t),\qquad
\mathbf R^\star_n=A(hn)+\mathbf b,
]
with (\mathbf u_n) periodic in the indices.

**Implementation-critical detail:** for edges crossing a periodic boundary, the edge vector must include the appropriate cell shift so that (\ell^\star) matches the affine base (same concept as “affine base + periodic displacement” in the continuum).

---

## What I can write next, immediately, if you want the cleanest paper insertion

1. A self-contained LaTeX block titled something like **“Periodic pre-tension and linearization”** containing exactly:

* (S[\mathbf X]) on (\mathbb T^d),
* (\mathbf X=\mathbf X_\star+\xi) with (\xi) periodic,
* the quadratic action (S^{(2)}[\xi]),
* the linearized PDE in the constant-coefficient form above,
* the spectral formulas for (S_\star) and (H_\star).

If you’d like that now, I’ll format everything as LaTeX with consistent symbol choices matching your paper (e.g. whether you use (\rho) vs (\rho_m), (K) vs (\kappa), (d=3) fixed vs general (d)).
