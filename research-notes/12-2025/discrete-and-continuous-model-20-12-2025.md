Absolutely — here is a **clean, complete formal description** of the brane model **with pre-tension** in both the **discrete** and **continuous** versions, each with a **dictionary** and an explicit **assumptions list**.

I’ll also address your note about **“pre-tension in (X^4)”** up front:

* **Pre-tension is not an embedding-component concept.**
  It is a **tangential/material stress** living in the brane’s **(d)-dimensional material directions** (indices (i,j)).
  There is **no separate “pre-tension along (X^4)”** because (X^4) is just one component of the embedding (\mathbf X\in\mathbb R^4), while tension is defined along the brane’s internal directions.
* You *can* choose a base state whose image lies in some 3D subspace of (\mathbb R^4) (a “flat brane”). That does **not** make (X^4) special; it is just a choice of orientation in (\mathbb R^4).

---

# 1) Discrete brane model with pre-tension (spring–mass lattice in (\mathbb R^4), no bending)

## 1.1 Model statement

Let the brane be a graph (G=(\mathcal V,\mathcal E)). In practice (\mathcal V) is a regular (d)-dimensional grid of nodes (usually (d=3)), and (\mathcal E) connects nearest neighbors in the material grid.

Each node (p\in\mathcal V) has an embedding position
[
\mathbf R_p(t)\in\mathbb R^4.
]

For each edge ((p,q)\in\mathcal E), define
[
\mathbf d_{pq}(t)=\mathbf R_q(t)-\mathbf R_p(t),\qquad
\ell_{pq}(t)=|\mathbf d_{pq}(t)|_{\mathbb R^4}.
]

Each edge is a Hooke spring with rest length (\ell_0>0) and spring constant (k>0):
[
U({\mathbf R_p})
=\sum_{(p,q)\in\mathcal E}\frac{k}{2},\big(\ell_{pq}-\ell_0\big)^2.
]

Each node has mass (m>0). The equations of motion are
[
m,\ddot{\mathbf R}_p(t)=\mathbf F_p(t)+\mathbf F^{\text{bc}}_p(t),
\qquad
\mathbf F_p(t)=-\frac{\partial U}{\partial \mathbf R_p}.
]

Here (\mathbf F^{\text{bc}}_p) denotes **constraint forces** (if you use clamped boundaries) required to maintain pre-tension; for periodic boundaries, (\mathbf F^{\text{bc}}_p\equiv 0).

The internal spring force can be written explicitly:
[
\boxed{
\mathbf F_p
===========

\sum_{q:,(p,q)\in\mathcal E}
k,(\ell_{pq}-\ell_0),\frac{\mathbf R_q-\mathbf R_p}{\ell_{pq}}.
}
]

### Pre-tensioned base state (discrete)

A discrete brane is **pre-tensioned** if there exists a time-independent configuration ({\mathbf R^\star_p}) such that:

1. Each edge has a constant length
   [
   \ell^\star_{pq}:=|\mathbf R^\star_q-\mathbf R^\star_p| \neq \ell_0
   ]
   (typically (\ell^\star_{pq}>\ell_0), i.e. uniform extension), and

2. The configuration is in static equilibrium under the chosen boundary conditions:
   [
   \mathbf F_p({\mathbf R^\star}) + \mathbf F^{\text{bc}}_p({\mathbf R^\star}) = \mathbf 0
   \quad\text{for all }p.
   ]

A standard way to construct a **uniformly pre-tensioned** state on a regular grid is an **affine embedding**:
[
\boxed{
\mathbf R^\star_{n} = A,x_n + \mathbf b,
}
]
where

* (n) is a multi-index (n\in\mathbb Z^d) labeling grid points,
* (x_n = h,n\in\mathbb R^d) are material grid coordinates with spacing (h),
* (A\in\mathbb R^{4\times d}) has rank (d) (a linear embedding of the brane into (\mathbb R^4)),
* (\mathbf b\in\mathbb R^4) is a translation.

For a nearest-neighbor edge in material direction (e_i),
[
\ell^\star_i=|A(h e_i)|.
]
Uniform pre-tension corresponds to (\ell^\star_i) being constant (over the lattice) and (\ell^\star_i\neq \ell_0).

**Important:** this construction does not privilege (X^4). Any global rotation in (\mathbb R^4) (multiplying all (\mathbf R^\star_n) by an orthogonal (Q\in O(4))) produces an equivalent pre-tensioned brane.

---

## 1.2 Discrete dictionary

### Sets / indices

* (d): intrinsic/material dimension of the brane (typically (3)).
* (\mathcal V): set of nodes/vertices.
* (\mathcal E): set of undirected edges (springs).
* (p,q\in\mathcal V): node indices.
* (n\in\mathbb Z^d): multi-index for regular grids.
* (A\in{1,2,3,4}): embedding component index.

### Variables

* (\mathbf R_p(t)\in\mathbb R^4): position of node (p).
* (\dot{\mathbf R}_p(t)), (\ddot{\mathbf R}_p(t)): velocity, acceleration.
* (\mathbf d_{pq}=\mathbf R_q-\mathbf R_p): edge vector in (\mathbb R^4).
* (\ell_{pq}=|\mathbf d_{pq}|): current spring length in (\mathbb R^4).
* (\mathbf F_p): internal force from springs.
* (\mathbf F^{\text{bc}}_p): boundary/constraint force needed to maintain imposed pre-tension (if applicable).
* (\mathbf R^\star_p): pre-tensioned equilibrium configuration.

### Parameters

* (m>0): mass per node.
* (k>0): spring constant.
* (\ell_0>0): spring rest length (“rest_length” parameter).
* (h>0): material grid spacing (geometry of node index set).
* (A\in\mathbb R^{4\times d}), (\mathbf b\in\mathbb R^4): parameters defining an affine base embedding (one convenient way to impose uniform pre-tension).

---

## 1.3 Discrete assumptions (explicit)

1. **Embedding space:** positions live in Euclidean (\mathbb R^4) with the standard norm.
2. **No bending:** only rest-length spring stretching energy; no curvature or angle penalties.
3. **Hookean springs:** each edge contributes (\tfrac{k}{2}(\ell-\ell_0)^2).
4. **Mass points:** kinetic/inertial model (m\ddot{\mathbf R}=\mathbf F).
5. **Pre-tension requires constraints:** if (\ell^\star\neq\ell_0), you must enforce boundary conditions (periodic cell geometry or clamped boundaries) so the lattice cannot relax back to (\ell=\ell_0).
6. **Homogeneity (typical):** constant (m,k,\ell_0) and uniform neighbor pattern (unless explicitly varied).
7. **No dissipation / damping** unless you add it externally.
8. **Time integration** is an approximation (numerical solver), but the model is continuous-time ODE.

---

# 2) Continuous brane model with pre-tension (hyperelastic (d)-brane embedded in (\mathbb R^4), no bending)

This is a **true continuum** model: no discrete springs, no preferred lattice directions.

## 2.1 Model statement (kinematics)

Let (\Omega\subset\mathbb R^d) be the material domain with coordinates (x=(x^1,\dots,x^d)). The brane configuration is an embedding field
[
\mathbf X:\Omega\times\mathbb R\to\mathbb R^4,\qquad \mathbf X(x,t)=(X^1,\dots,X^4).
]

Define tangent vectors and induced metric:
[
\partial_i\mathbf X(x,t)\in\mathbb R^4,\qquad
g_{ij}(x,t) := \partial_i\mathbf X\cdot\partial_j\mathbf X.
]

### Rest geometry (reference metric)

A “rest-length” in the continuum is encoded by a **reference metric** (g^0_{ij}(x)). For a homogeneous, isotropic rest geometry with a single preferred length scale (\ell_0),
[
\boxed{
g^0_{ij} = \ell_0^2,\delta_{ij}.
}
]
This choice means: the brane is stress-free when (g=g^0).

(Using (\delta_{ij}) here avoids introducing an identity matrix symbol (I).)

Define the relative mixed metric
[
\widehat g^{i}{}*{j} := (g^0)^{ik},g*{kj},
]
where ((g^0)^{ik}) is the inverse of (g^0_{ik}).

Let (\lambda_a) be eigenvalues of (\widehat g^{i}{}_{j}) (all positive if the embedding is regular). Define relative principal stretches
[
\sigma_a := \sqrt{\lambda_a},\qquad a=1,\dots,d.
]
The stress-free condition is (\sigma_a=1) for all (a).

---

## 2.2 Constitutive law (stretch-only, isotropic, no bending)

Choose an isotropic strain energy density (W) minimized at (g=g^0). A direct “rest-length” analogue is:
[
\boxed{
W(g;g^0)
========

\frac{K}{2}\sum_{a=1}^d (\sigma_a-1)^2,
\qquad \sigma_a=\sqrt{\lambda_a(\widehat g)}.
}
]

* (K>0) is the stiffness scale (for (d=3): units Pa).
* (W) depends only on **first derivatives** of (\mathbf X) via (g_{ij}). (No bending.)

---

## 2.3 Dynamics (PDE)

Let (\rho_m>0) be mass density per material volume (kg/m(^d)). The evolution is:
[
\boxed{
\rho_m,\partial_{tt}X^A = \partial_i P_{Ai},
\qquad A=1,\dots,4,
}
]
where the first Piola-type stress is defined by
[
\boxed{
P_{Ai} := \frac{\partial W}{\partial(\partial_i X^A)}.
}
]

Using the chain rule through the metric,
[
\boxed{
P_{Ai} = 2,\frac{\partial W}{\partial g_{ij}},\partial_j X^A.
}
]

### Pre-tensioned base state (continuous)

A continuum brane is **pre-tensioned** if you choose a time-independent equilibrium embedding (\mathbf X_\star(x)) such that:

1. It is an equilibrium:
   [
   \partial_i P_{Ai}(\mathbf X_\star)=0 \quad \text{in }\Omega
   ]
   (with boundary conditions below), and

2. It is **not stress-free**:
   [
   g_\star \neq g^0 \quad\text{equivalently}\quad \sigma_a(\widehat g_\star)\neq 1
   ]
   so the resulting stress is nonzero.

A standard uniform pre-tension is obtained by an **affine base embedding**
[
\boxed{
\mathbf X_\star(x) = A x + \mathbf b,\qquad A\in\mathbb R^{4\times d},\ \mathrm{rank}(A)=d,
}
]
giving constant metric (g_{\star,ij}=(A^\mathsf T A)*{ij}).
Uniform pre-tension corresponds to (g*\star) being a constant multiple (or otherwise constant mismatch) relative to (g^0), e.g.
[
g_\star = s^2 g^0\quad\text{with } s\neq 1.
]

**How it is maintained:** You must impose boundary conditions (or external loads) that keep (\mathbf X) at/near (\mathbf X_\star). Typical choices:

* **Periodic** domain with a fixed cell shape consistent with (\mathbf X_\star).
* **Dirichlet/clamped:** (\mathbf X|*{\partial\Omega}=\mathbf X*\star|_{\partial\Omega}).
* **Traction** boundary: apply boundary tractions matching the desired pre-stress.

Finally, you evolve perturbations (\boldsymbol\xi) around the pre-tensioned base:
[
\mathbf X(x,t)=\mathbf X_\star(x)+\boldsymbol\xi(x,t).
]

Again: none of this treats (X^4) specially. The stress lives in material indices (i,j); the embedding indices (A) enter through (\mathbf X) and its derivatives.

---

## 2.4 Continuous dictionary

### Sets / indices

* (d): brane material dimension (typically (3)).
* (\Omega\subset\mathbb R^d): material domain.
* (x^i): material coordinates, (i=1,\dots,d).
* (t): time parameter.
* (A\in{1,2,3,4}): embedding-space component index.

### Unknown field

* (\mathbf X(x,t)\in\mathbb R^4): embedding map (brane configuration).
* (X^A(x,t)): components of (\mathbf X).

### Derived geometric objects

* (\partial_i \mathbf X): tangent vectors in (\mathbb R^4).
* (g_{ij}=\partial_i\mathbf X\cdot\partial_j\mathbf X): induced metric.
* (g^0_{ij}): reference (rest) metric.
* ((g^0)^{ij}): inverse reference metric.
* (\widehat g^{i}{}*{j}=(g^0)^{ik}g*{kj}): relative mixed metric.
* (\lambda_a): eigenvalues of (\widehat g).
* (\sigma_a=\sqrt{\lambda_a}): relative principal stretches.

### Material parameters

* (\rho_m): mass density per material volume.
* (K): stiffness/modulus scale in (W).
* (\ell_0): rest-length scale encoded in (g^0_{ij}=\ell_0^2\delta_{ij}).

### Constitutive quantities

* (W(g;g^0)): strain energy density (stretch-only).
* (P_{Ai}=\partial W/\partial(\partial_i X^A)): first Piola-type stress (maps material directions to embedding forces).

### Base state / pre-tension

* (\mathbf X_\star(x)): time-independent pre-tensioned equilibrium embedding.
* (g_{\star,ij}): induced metric of the base state.
* (s): uniform stretch factor if (g_\star=s^2 g^0).

### Standard symbols

* (\delta_{ij}): Kronecker delta.
* (\partial_i=\partial/\partial x^i), (\partial_{tt}=\partial^2/\partial t^2).

---

## 2.5 Continuous assumptions (explicit)

1. **Embedding space:** (\mathbf X) maps into Euclidean (\mathbb R^4) with standard inner product.
2. **No bending:** energy depends only on first derivatives (metric), not on curvature ((\partial_i\partial_j\mathbf X)).
3. **Hyperelastic:** stresses derive from a scalar energy density (W(g;g^0)).
4. **Isotropy (material):** (W) depends only on eigenvalues/invariants of ((g^0)^{-1}g), not on preferred directions.
5. **Regular embedding:** (g_{ij}) is positive definite (rank (d)) in the regime of interest (no local collapse).
6. **Pre-tension is an imposed equilibrium:** you choose (\mathbf X_\star) with (g_\star\neq g^0) and maintain it via boundary conditions/loads; otherwise the system relaxes toward (g=g^0).
7. **Homogeneous parameters** ((\rho_m,K,\ell_0) constant) unless specified otherwise.
8. **No dissipation** unless explicitly added.
9. **Classical continuum validity:** fields are sufficiently smooth for derivatives/eigenvalues to make sense (or interpreted weakly if you go to FEM/weak forms).

---

## A quick consistency remark about your “(X^4) pretension = 0” intuition

The correct statement is:

* Pre-tension is a **tangential stress state** on the brane (material indices).
* There is no “extra material direction” corresponding to the normal/4th embedding axis, so there is no independent “pretension along (X^4)”.
* Nevertheless, because spring lengths / induced metric are computed in (\mathbb R^4), **all embedding components are coupled** in the dynamics.
