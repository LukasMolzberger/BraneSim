## 1. Story in words

1. **There is a 4-dimensional brane** (3 space + 1 time) sitting inside a higher-dimensional flat space (here: 5D).

2. The brane is described by an **embedding map**
   [
   X^A(x^\mu)
   ]
   which tells you: for each brane coordinate (x^\mu) (“where you are in spacetime”), where that point lives in the surrounding flat 5D space (labelled by (X^A)).

3. The **geometry experienced by physics on the brane** is not the ambient 5D geometry, but the **induced metric**
   [
   g_{\mu\nu}(x) = \eta_{AB},\partial_\mu X^A,\partial_\nu X^B.
   ]
   This is just: “pull back the ambient metric to the brane”.

4. You pick a “static gauge” where, if there were no deformations, the brane would just be
   [
   X^\mu = x^\mu, \quad X^4 = 0.
   ]
   Any **physical deformation** is encoded as a **displacement field**
   [
   X^\mu(x) = x^\mu + u^\mu(x), \quad X^4(x) = \xi_\perp(x).
   ]

5. The **core idea** is then:

   > What we call “gravity” is simply the **slow, long-wavelength part** of the displacement field (u^\mu(x)):
   > *time dilation* is the time component (u^0(x)), and
   > *spatial curvature / contraction* is the spatial part (u^i(x)).
   > All matter and light see this through the induced metric (g_{\mu\nu}[u,\xi_\perp]).

   In particular, a **massive object** corresponds to a configuration of (\xi_\perp) (amplitude bulge) that, via the brane elasticity, induces a **signless inward lateral contraction** of the brane, encoded in (u^i(x)). That contraction pattern generates a metric that, in the weak field limit, looks like Einstein’s.

That’s the conceptual piece. Now let’s write it more precisely.

---

## 2. Formal setup (slowly)

### 2.1 Ambient space and brane coordinates

* Ambient 5D coordinates: (X^A), with (A=0,1,2,3,4).
* Brane coordinates (the ones we use for physics): (x^\mu), with (\mu =0,1,2,3).
* Ambient metric:
  [
  \eta_{AB} = \mathrm{diag}(-1,1,1,1,\sigma),
  ]
  where (\sigma = +1) if the transverse direction is Euclidean.

The brane is a map
[
X^A: \mathbb{R}^{1,3} \to \mathbb{R}^{1,4}, \quad X^A = X^A(x^\mu).
]

### 2.2 Static gauge + displacement fields

We choose a gauge where, without any deformations:

* (X^\mu = x^\mu) is the identity (flat brane),
* (X^4 = 0) (no transverse displacement).

We then describe deformations by:

[
\boxed{
\begin{aligned}
X^\mu(x) &= x^\mu + u^\mu(x), \
X^4(x)   &= \xi_\perp(x).
\end{aligned}
}
]

* (u^\mu(x)) is the **intrinsic 4D displacement field** (“how the brane stretches and shears within its own 4D directions”).
* (\xi_\perp(x)) is the **transverse / amplitude field** (“how the brane bulges in the extra dimension”).

### 2.3 Induced metric from displacements

The metric seen by any field living on the brane is

[
\boxed{
g_{\mu\nu}(x) = \eta_{AB},\partial_\mu X^A,\partial_\nu X^B.
}
]

Compute the derivatives:

* (\partial_\mu X^\nu = \delta_\mu^{\ \nu} + \partial_\mu u^\nu),
* (\partial_\mu X^4 = \partial_\mu \xi_\perp).

Plugging in and expanding for **small deformations**:

[
\begin{aligned}
g_{\mu\nu}
&= \eta_{\alpha\beta}(\delta_\mu^{\ \alpha} + \partial_\mu u^\alpha)(\delta_\nu^{\ \beta} + \partial_\nu u^\beta)

* \sigma,\partial_\mu \xi_\perp ,\partial_\nu \xi_\perp \
  &\simeq \eta_{\mu\nu}
* \partial_\mu u_\nu + \partial_\nu u_\mu
* \sigma,\partial_\mu \xi_\perp ,\partial_\nu \xi_\perp
* \mathcal O(u^2).
  \end{aligned}
  ]

So, to **linear order in small strains**:

[
\boxed{
h_{\mu\nu} := g_{\mu\nu} - \eta_{\mu\nu}
;\approx; \partial_\mu u_\nu + \partial_\nu u_\mu.
}
]

The amplitude field (\xi_\perp) only enters **quadratically** in this linear regime (through (\partial_\mu\xi_\perp,\partial_\nu\xi_\perp)). Its main role for gravity is: it creates energy density and stresses that drive (u^\mu) via the brane elasticity.

### 2.4 A simple “gravity” pattern: isotropic lateral contraction

We now pick a particular form of (u^\mu) that encodes a **gravitational potential** (\Phi_G(\mathbf x)) in the weak-field regime:

[
\boxed{
u^0(x) = \frac{\Phi_G(\mathbf x)}{c^2}, ct, \qquad
u^i(x) = -,\frac{\Phi_G(\mathbf x)}{c^2}, x^i.
}
]

This means:

* Time dilation: the local time coordinate is stretched by (1 + \Phi_G/c^2).
* Lateral contraction: the spatial coordinates are shrunk by (1 - \Phi_G/c^2).

Compute the relevant derivatives (schematically):

* (\partial_0 u^0 = \Phi_G/c^2),
* (\partial_i u^0 \approx 0) (static potential),
* (\partial_j u^i \approx -(\Phi_G/c^2),\delta^i_{\ j}), ignoring gradients of (\Phi_G) itself for the moment.

Insert into (h_{\mu\nu} \approx \partial_\mu u_\nu + \partial_\nu u_\mu), and you get:

[
g_{00} \simeq -\left(1 + \frac{2\Phi_G}{c^2}\right),
\qquad
g_{ij} \simeq \left(1 - \frac{2\Phi_G}{c^2}\right)\delta_{ij},
]

which is **exactly** the standard weak-field GR metric.

So the **metric of Einstein’s gravitational field** is literally:

> the symmetric gradient of your 4D displacement field
> (which encodes lateral contraction and time dilation).

---

## 3. The core assumption, now very explicitly

We can now state the **core assumption** in your model in three parts:

1. **Geometric part (what gravity *is*):**

   * The brane is embedded via (X^A(x) = (x^\mu + u^\mu(x),,\xi_\perp(x))).
   * The **physical spacetime metric** is the induced metric
     [
     g_{\mu\nu}(x) = \eta_{AB},\partial_\mu X^A,\partial_\nu X^B.
     ]
   * In the weak-field regime, the gravitational field is captured by the slow part of the displacement:
     [
     h_{\mu\nu} \approx \partial_\mu u_\nu + \partial_\nu u_\mu.
     ]
   * A mass distribution (localized energy in (\xi_\perp) etc.) induces a displacement (u^\mu) that, far away, looks like an **inward lateral contraction** plus time dilation, and therefore produces the **Einstein weak-field metric**.

2. **Matter coupling part (how gravity couples):**

   * All excitations (photon-like, electron-like solitons, phonons) live on the brane and have kinetic terms built from (g_{\mu\nu}).
   * After coarse graining, their effective actions use **minimal coupling**:
     [
     S_{\text{matter}}[\varphi^I,g] = \int d^4x \sqrt{-g},(\text{covariant kinetic terms and interactions}).
     ]
   * This realizes the **equivalence principle**: “gravity” is not an extra force; it is the geometry generated by (u^\mu).

3. **Source part (what generates the displacement):**

   * The brane has an elastic Lagrangian (tension, bending, saturation, …) which defines a microscopic stress tensor (\tau_{\mu\nu}[\mathbf X]).
   * Coarse-grained energy density
     [
     \rho_\text{eff}(x) \sim \langle \tau_{00}[\mathbf X] \rangle / c^2
     ]
     sources the slow-varying displacement (u^\mu), and hence (\Phi_G).
   * In the long-wavelength limit, the static equilibrium equations for (u^\mu) reduce to something Poisson-like:
     [
     \nabla^2 \Phi_G(\mathbf x) \approx 4\pi G \rho_\text{eff}(\mathbf x),
     ]
     reproducing Newtonian gravity; dynamically, the emergent evolution of (g_{\mu\nu}) is expected to be Einstein-like (induced gravity picture).

That’s “lateral contraction as gravity” in a clean, explicit form.

---

## 4. Symbol dictionary

Here’s a small dictionary of the main symbols we used, grouped by type.

### Indices and coordinates

* (A,B = 0,1,2,3,4): indices in the **ambient 5D space**.
* (\mu,\nu = 0,1,2,3): indices in the **brane (spacetime) coordinates**.
* (i,j = 1,2,3): **spatial** indices on the brane.
* (X^A): ambient coordinates.
* (x^\mu): intrinsic brane coordinates, used as physical spacetime labels.
* (t): time coordinate, so (x^0 = ct).
* (\mathbf x = (x^1,x^2,x^3)): spatial position on the brane.

### Fields and embeddings

* (X^A(x)): embedding of the brane into the ambient space.
* (u^\mu(x)): **intrinsic displacement field** of the brane in its own 4D directions.

  * (u^0(x)): time component (time dilation / time warping).
  * (u^i(x)): spatial components (lateral contraction / shear).
* (\xi_\perp(x) \equiv X^4(x)): **transverse / amplitude field** (brane height in the extra dimension).
* (\Phi_G(\mathbf x)): **emergent gravitational potential** in the Newtonian / weak-field limit.

### Geometric quantities

* (\eta_{AB}): ambient flat metric, typically (\mathrm{diag}(-1,1,1,1,\sigma)).
* (g_{\mu\nu}(x)): **induced spacetime metric** on the brane,
  [
  g_{\mu\nu} = \eta_{AB},\partial_\mu X^A,\partial_\nu X^B.
  ]
* (\eta_{\mu\nu}): flat Minkowski metric on the brane, (\mathrm{diag}(-1,1,1,1)).
* (h_{\mu\nu} = g_{\mu\nu} - \eta_{\mu\nu}): metric perturbation (weak-field).
* (R[g]): Ricci scalar of the metric (g_{\mu\nu}).
* (G_{\mu\nu}[g]): Einstein tensor built from (g_{\mu\nu}).

### Matter and energy

* (\mathcal L[\mathbf X]): **microscopic brane Lagrangian** (tension, bending, etc.).
* (\tau_{\mu\nu}[\mathbf X]): **microscopic stress tensor** derived from (\mathcal L).
* (T_{\mu\nu}^\text{(eff)}(x)): **coarse-grained (effective) stress–energy tensor**, e.g.
  [
  T_{\mu\nu}^\text{(eff)}(x)
  = \langle \tau_{\mu\nu}[\mathbf X](x) \rangle_{\text{cell}}.
  ]
* (\rho_\text{eff}(\mathbf x) = T_{00}^\text{(eff)}(\mathbf x)/c^2): effective mass density.
* (\Phi_G(\mathbf x)): potential that appears in the weak-field metric components:
  [
  g_{00} \simeq -\left(1 + 2\Phi_G/c^2 \right),\quad
  g_{ij} \simeq \left(1 - 2\Phi_G/c^2 \right)\delta_{ij}.
  ]

### Scales and constants

* (c): speed of light (also small-amplitude wave speed of the brane).
* (G): Newton’s gravitational constant (to be matched from brane parameters).
* (\ell_\ast): coarse-graining scale (many lattice spacings, small macroscopically).
* (\sigma): sign of the transverse metric component (usually (+1) for Euclidean).

