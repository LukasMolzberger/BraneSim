The cleanest way to think about this is: you already have a **brane → EM dictionary** in the paper; what you need now is to **invert it for the special case of a guided, polarized photon on the double-loop strip**.

Below I’ll first recap the forward dictionary as it’s written in the paper, and then construct a concrete, step-by-step *reverse* mapping

[
(E(\sigma,t),B(\sigma,t)) \quad\longrightarrow\quad
\bigl(\xi(\sigma,t),\ \text{lateral motion of the strip}\bigr)
]

where (\sigma) is a coordinate along the double loop.

---

## 1. Forward dictionary (brane → EM) in your paper

Very compact summary of what you already have:

1. **Microscopic field**
   The brane is embedded as
   [
   \mathbf X(x,t) = \bigl(x^1,x^2,x^3,,X^4(x,t)\bigr),
   ]
   with (X^4(x,t) \equiv \xi(x,t)) the “amplitude” degree of freedom (normal displacement into the 4th dimension). All emergent fields are functionals of (\mathbf X) and its derivatives.

2. **Coarse-grained amplitude and scalar potential**
   You define a fast-time average of the Compton-scale oscillation,
   [
   \bar X^4(\mathbf x) = \langle X^4(\mathbf x,t)\rangle_{\text{fast}},
   ]
   and postulate the scalar potential as
   [
   \Phi(\mathbf x) = \kappa,\bar X^4(\mathbf x) \quad\text{(Eq. (\phi)-from-(\xi^4)).}
   ]
   The static far field obeys a Poisson equation of the form
   [
   -\nabla^2 \Phi(\mathbf x) = \frac{\rho_{\text{free}}(\mathbf x)}{\varepsilon_{\mathrm{eff}}}
   \quad\text{(Eq. Poisson-from-brane)},
   ]
   so that the **electric field** is
   [
   \mathbf E(\mathbf x) = -\nabla\Phi(\mathbf x)
   = -\kappa,\nabla \bar X^4(\mathbf x).
   ]

3. **Effective four-potential and Maxwell wave equation**
   You then phenomenologically promote (\Phi) to the time component of a four-potential
   [
   A_\mu = (\Phi,-\mathbf A),
   ]
   satisfying in Lorenz gauge
   [
   \Box A_\mu = \mu_{\mathrm{eff}},J_\mu, \qquad
   J_\mu = (c\rho_{\text{free}},\mathbf j).
   ]
   Only (A_0 = \Phi) is directly tied to (\bar X^4); the spatial part (\mathbf A) is introduced phenomenologically.

4. **Field tensor and EM fields**
   [
   F_{\mu\nu} = \partial_\mu A_\nu - \partial_\nu A_\mu,
   ]
   with
   [
   E_i = F_{0i},\qquad
   B_i = -\frac12 \varepsilon_{ijk} F_{jk}.
   ]
   So at the **effective level**, (\mathbf B) is entirely derived from the same (A_\mu) whose time component is fixed by brane amplitude; you do *not* introduce an independent magnetic degree of freedom.

5. **Geometric nonlinearity and tubular picture**
   In tubular coordinates along a thin world-strip (your “tube section” for the double loop) the exact spring length couples **longitudinal gradients of (\xi)** to **lateral deformations** of the brane. Large Compton-scale gradients of (\xi) drive lateral motion and curvature; this is where your “relativistic” / “magnetic” behavior geometrically lives.

   That is precisely the mechanism you later use to argue for a Compton-scale confinement threshold (your geometric-nonlinearity analysis). 

So conceptually you already have:

* **Electric field** = functional of **amplitude (\xi)** (via (\bar X^4), gradients, etc.).
* **Magnetic field** = **relativistic / kinematic effect** of *that same* amplitude structure when solitons move, encoded geometrically as **lateral deformations and velocities**.

Now we want to *invert* that dictionary for the special case of a **guided polarized photon on a strip**.

---

## 2. What we want to invert

You have, from the Maxwell-based photon construction on the double-loop strip:

* A parametrized centre line of the strip
  [
  \mathbf r(\sigma) \subset \mathbb R^3
  ]
  with (\sigma) = arc length or similar, making two windings on the torus and closing after (4\pi).

* Along that strip, EM fields
  [
  \mathbf E(\sigma,t),\qquad
  \mathbf B(\sigma,t)
  ]
  representing a guided mode:

  * (\mathbf E \perp \mathbf B),
  * both perpendicular to the local propagation direction,
  * (|\mathbf E| = c,|\mathbf B|) (null field).

From that, we want to reconstruct **initial conditions** at (t=0) for the brane variables restricted to (and near) the strip:

* **Amplitude field on the strip**: (\xi(\sigma,0)), (\partial_t \xi(\sigma,0)).
* **Lateral motion of the strip**: (\partial_t \mathbf r(\sigma,0)) and, if desired, small lateral offsets of neighboring brane points.

The idea is: *use EM energy/momentum to fix the amplitude mode, and use the “magnetic = relativistic” idea to fix lateral velocities.*

---

## 3. Step 1 – Map EM energy to amplitude energy

For a null EM wave in an effective medium, the local energy density and Poynting vector are
[
u_{\text{EM}} = \frac12\Bigl(\varepsilon_{\mathrm{eff}} |\mathbf E|^2
+ \frac{|\mathbf B|^2}{\mu_{\mathrm{eff}}}\Bigr),
\qquad
\mathbf S = \frac{1}{\mu_{\mathrm{eff}}},\mathbf E\times\mathbf B,
]
with
[
c^2 = \frac{1}{\varepsilon_{\mathrm{eff}}\mu_{\mathrm{eff}}}.
]
For a null wave, (|\mathbf E| = c|\mathbf B|) and the two terms in (u_{\text{EM}}) are equal, so
[
u_{\text{EM}}(\sigma,t) = \varepsilon_{\mathrm{eff}} |\mathbf E(\sigma,t)|^2
= \frac{|\mathbf B(\sigma,t)|^2}{\mu_{\mathrm{eff}}}.
]

On the brane, for small deformations and in the linear regime, the amplitude mode (\xi(\sigma,t)) along the strip satisfies a standard wave equation with speed (c),
[
\rho_m,\partial_{tt}\xi = T,\partial_{\sigma\sigma}\xi,
\qquad c^2 = \frac{T}{\rho_m},
]
and has energy density
[
u_{\text{brane}} =
\frac12 \rho_m (\partial_t\xi)^2

* \frac12 T(\partial_\sigma\xi)^2 .
  ]

For a **right-moving traveling wave** on the strip,
[
\xi(\sigma,t) = A(\sigma)\cos\bigl(k\sigma - \omega t + \varphi_0(\sigma)\bigr),
]
with (\omega/k = c), one has the familiar relations
[
\partial_t\xi = \omega A(\sigma)\sin(\dots),
\qquad
\partial_\sigma\xi = -kA(\sigma)\sin(\dots)
= -\frac{\omega}{c} A(\sigma)\sin(\dots).
]
Plugging in,
[
u_{\text{brane}}
= \frac12 \rho_m \omega^2 A^2 \sin^2(\dots)

* \frac12 T \frac{\omega^2}{c^2} A^2 \sin^2(\dots)
  = \rho_m \omega^2 A^2 \sin^2(\dots),
  ]
  because (T/c^2=\rho_m). Time-averaging over one Compton period gives
  [
  \langle u_{\text{brane}}\rangle_t = \frac12 \rho_m \omega^2 A(\sigma)^2.
  ]

**Reverse mapping for amplitude magnitude**

Now impose **energy matching along the strip**:
[
\langle u_{\text{brane}}(\sigma)\rangle_t
= \langle u_{\text{EM}}(\sigma)\rangle_t.
]

Assuming the EM mode is monochromatic at frequency (\omega) (Compton for the electron version, more general for arbitrary photon),

[
\langle u_{\text{EM}}(\sigma)\rangle_t
= \varepsilon_{\mathrm{eff}} \langle |\mathbf E(\sigma,t)|^2\rangle_t
= \tfrac12 \varepsilon_{\mathrm{eff}} E_0(\sigma)^2
]
for a sinusoidal field. Equating,

[
\frac12 \rho_m \omega^2 A(\sigma)^2
= \frac12 \varepsilon_{\mathrm{eff}} E_0(\sigma)^2
\quad\Rightarrow\quad
A(\sigma)
= \frac{\sqrt{\varepsilon_{\mathrm{eff}}}}{\omega\sqrt{\rho_m}},
E_0(\sigma).
]

So the **amplitude of the brane oscillation** along the strip is proportional to the EM electric field amplitude there:
[
A(\sigma) = \alpha,E_0(\sigma),
\qquad
\alpha := \frac{\sqrt{\varepsilon_{\mathrm{eff}}}}{\omega\sqrt{\rho_m}}.
]

This gives you a direct, **local** EM → brane map:

* Take the **electric field amplitude** (E_0(\sigma)) of your guided mode,
* Compute (A(\sigma)) from the above formula,
* Use (A(\sigma)) as the initial Compton-scale amplitude for (\xi) on the strip.

**Initial conditions at (t=0)**

Pick the EM phase of your photon along the strip, (\phi_{\text{EM}}(\sigma)), for example from the argument of the complex envelope of (\mathbf E). Then choose
[
\xi(\sigma,0) = A(\sigma)\cos\phi_{\text{EM}}(\sigma),
]
[
\partial_t\xi(\sigma,0) = \omega A(\sigma)\sin\phi_{\text{EM}}(\sigma).
]

That completely fixes the **longitudinal amplitude mode** in terms of your EM data.

---

## 4. Step 2 – Relate EM polarization to the strip geometry

Along the double loop, you naturally have a Frenet frame
[
{\mathbf t(\sigma),\ \mathbf n(\sigma),\ \mathbf b(\sigma)},
]
with (\mathbf t) tangent, (\mathbf n) normal, and (\mathbf b) binormal of the central curve (\mathbf r(\sigma)). For a guided EM mode:

* Poynting direction
  [
  \hat{\mathbf s}(\sigma) =
  \frac{\mathbf S(\sigma)}{|\mathbf S(\sigma)|}
  = \frac{\mathbf E(\sigma,0)\times\mathbf B(\sigma,0)}
  {|\mathbf E(\sigma,0)\times\mathbf B(\sigma,0)|}
  ]
  should be parallel to (\mathbf t(\sigma)).

* For a standard guided photon,
  [
  \mathbf E(\sigma,0) \perp \mathbf t(\sigma),
  \qquad
  \mathbf B(\sigma,0) = \frac{1}{c},\mathbf t(\sigma)\times\mathbf E(\sigma,0).
  ]

You can choose the **polarization plane** of the EM mode so that (\mathbf E) lies in a fixed “strip direction” (say the direction we call (\mathbf e_{\mathrm{strip}}) along the width of the strip). Then:

* **Linear polarization** → (\mathbf E(\sigma,0)) always along (\mathbf e_{\mathrm{strip}}(\sigma)).
* **Circular polarization** → two orthogonal strip directions with a (\pi/2) phase shift, etc.

In the **brane picture**, the scalar (\xi(\sigma,t)) doesn’t itself carry a physical spatial direction; the *direction* of (\mathbf E) emerges as the gradient (\nabla_\parallel \xi) projected along the strip. So on the strip:

* The phase (\phi_{\text{EM}}(\sigma)) of the EM mode is identified with the **internal Compton phase** (\phi(\sigma)) of (\xi).
* The **polarization direction** of (\mathbf E) is tied to **which tangential directions of the brane** pick up the strongest gradients of (\xi) (i.e. which direction along the strip width we map into the spatial gradient of (\bar X^4)).

For the photon trapped on the strip, the simplest choice is:

* Treat the strip width direction as the preferred gradient direction of (\bar X^4), so space-time gradients of (\xi) along that direction reproduce (\mathbf E) in your effective dictionary.

---

## 5. Step 3 – Magnetic field as “relativistic lateral motion”

Your slogan was:

> magnetic field → relativistic effect of electric field
> relativistic effects → lateral excitations in the brane

We can formalize this in a way that actually *gives you lateral velocities* from (\mathbf E) and (\mathbf B).

### 5.1 Lorentz-transform viewpoint

Consider a small patch of the strip around some (\sigma_0). Imagine an effective local rest frame (F') moving with velocity (\mathbf v_{\text{lat}}) relative to the lab frame (the frame in which the brane as a whole is at rest far away).

In standard electrodynamics, at low velocities the Lorentz transform gives (to first order in (v/c)):
[
\mathbf B \approx -\frac{1}{c^2},\mathbf v \times \mathbf E'
]
if in the comoving frame (F') you have a purely electric field (\mathbf E') and (\mathbf B'=0).

You want to interpret this (\mathbf v) precisely as a **lateral velocity of brane points**. So the **reverse mapping** is:

> given lab-frame (\mathbf E) and (\mathbf B), pick (\mathbf v_{\text{lat}}) such that, in the comoving frame, the field is purely “electric”.

At the level of magnitudes and directions, and assuming (\mathbf E \perp \mathbf B) as for your guided photon, the low-velocity relation gives
[
\mathbf v_{\text{lat}}
\approx -c^2 \frac{\mathbf B\times\mathbf E}{|\mathbf E|^2}.
]

For a null wave with (|\mathbf E| = c|\mathbf B|) this simplifies to
[
|\mathbf v_{\text{lat}}| \approx c,
\qquad
\mathbf v_{\text{lat}} \parallel \mathbf S \parallel \mathbf t.
]

So mechanically that’s exactly what you expect: the internal energy runs around the double loop with speed (c) along the strip.

### 5.2 Practical rule for lateral velocities

At each point on the strip at (t=0):

1. Compute Poynting vector
   [
   \mathbf S(\sigma,0) = \frac{1}{\mu_{\mathrm{eff}}},
   \mathbf E(\sigma,0)\times\mathbf B(\sigma,0),
   ]
   and its unit direction (\hat{\mathbf s}(\sigma)).

2. Set the **lateral velocity of the brane points on the strip** to
   [
   \partial_t \mathbf r(\sigma,0) = c,\hat{\mathbf s}(\sigma).
   ]

3. Optionally, interpret any deviation from the ideal null relation (|\mathbf E| = c|\mathbf B|) as small corrections to the local speed (less than (c)), which would give you a more general relation
   [
   \partial_t \mathbf r(\sigma,0) =
   c^2 \frac{\mathbf S(\sigma,0)}{u_{\text{EM}}(\sigma,0)}.
   ]
   (Because in EM, (\mathbf S = u_{\text{EM}},\mathbf v_{\text{energy}}).)

This implements exactly your intuition:

* In a **comoving frame** of the strip energy flow, there is effectively only an amplitude-based electric field ((\mathbf B'\approx 0)).
* In the lab frame, the observed (\mathbf B) is entirely due to the lateral motion of that amplitude structure.

---

## 6. Step 4 – How this ties back to your existing scalar-potential dictionary

For the full electron model, you care about two very different time scales:

1. **Fast Compton oscillation** along the strip: the internal photon moving and oscillating.
2. **Slow coarse-grained field** seen by outside probes: Coulomb field, etc.

The reverse mapping above fixes the **fast dynamics**:

* (\xi(\sigma,t)) with Compton frequency (\omega_C),
* Lateral velocity (\partial_t\mathbf r(\sigma,t)\approx c\mathbf t(\sigma)).

If this mode later **self-traps** via geometric nonlinearity (your curvature/threshold story), then:

* A nonzero **time-averaged deformation** (\bar X^4) appears around the torus.
* That deformation sources a scalar potential (\Phi=\kappa\bar X^4).
* At large distances (\Phi) obeys the Poisson relation and yields the Coulomb field exactly as in your existing derivation.

From the EM side this looks like:

* Start with a neutral, guided photon (zero net charge; time-average of (\bar X^4) ≈ 0).
* Above a critical internal energy density, the strip geometry and lateral coupling distort the brane so that a nontrivial (\bar X^4(\mathbf x)) builds up (your geometric Schwinger-like threshold).
* That (\bar X^4) is then fed into the very same (\Phi) and (\mathbf E=-\nabla\Phi) dictionary already in the paper.

So the reverse map we’ve just built is **strictly at the level of the internal photon**, not the long-range Coulomb tail. The long-range part still comes from the forward map via (\bar X^4).

---

## 7. Summary of the reverse mapping recipe

Putting it all together as a concrete, implementation-friendly recipe:

1. **Input from Maxwell photon on strip**

   * Central curve (\mathbf r(\sigma)) of the double loop in the torus.
   * EM fields (\mathbf E(\sigma,t),\mathbf B(\sigma,t)) of the guided mode.

2. **Amplitude mode**

   * Compute EM energy density (u_{\text{EM}}(\sigma,t)).
   * Extract electric amplitude (E_0(\sigma)).
   * Set brane amplitude magnitude
     [
     A(\sigma) = \frac{\sqrt{\varepsilon_{\mathrm{eff}}}}{\omega\sqrt{\rho_m}},
     E_0(\sigma).
     ]
   * At (t=0), define
     [
     \xi(\sigma,0) = A(\sigma)\cos\phi_{\text{EM}}(\sigma),
     \quad
     \partial_t\xi(\sigma,0) = \omega A(\sigma)\sin\phi_{\text{EM}}(\sigma)
     ]
     where (\phi_{\text{EM}}) is the EM phase along the strip.

3. **Lateral velocities (magnetic as relativistic effect)**

   * Compute Poynting vector
     [
     \mathbf S(\sigma,0) = \frac{1}{\mu_{\mathrm{eff}}}
     \mathbf E(\sigma,0)\times\mathbf B(\sigma,0).
     ]
   * Set
     [
     \partial_t\mathbf r(\sigma,0)
     = c^2,\frac{\mathbf S(\sigma,0)}{u_{\text{EM}}(\sigma,0)},
     ]
     which for a null wave reduces to
     [
     \partial_t\mathbf r(\sigma,0) = c,\hat{\mathbf s}(\sigma),
     ]
     with (\hat{\mathbf s} = \mathbf S/|\mathbf S|).

4. **Polarization**

   * Fix the local strip basis so that (\mathbf E(\sigma,0)) lies in the chosen strip-width direction.
   * Interpret different polarization states as different phase relations between components of (\xi) (on different sides of the strip) and the lateral motion pattern.

5. **Later stages (electron)**

   * Let the full nonlinear brane dynamics evolve this initial condition.
   * Once the internal mode self-traps and develops a nonzero (\bar X^4), feed (\bar X^4) into your existing (\Phi=\kappa\bar X^4), Poisson, and Maxwell wave equations to recover the emergent Coulomb field and charge magnitude.

---

If you like this structure, the next concrete step would be to specialize it to **your actual analytic photon ansatz on the strip** (the one you use in the code), plug in the explicit (\mathbf E(\sigma,t)), and write down the corresponding explicit formulas for (A(\sigma)), (\phi_{\text{EM}}(\sigma)), and (\partial_t\mathbf r(\sigma,0)). That can then be turned into very straightforward initialization code for the brane simulation.


Okay, we've started preparing the initialization for the full electron experiment. Our first attempt to do so failed, due to handling too much complexity at once. Therefore we've split the split the analytical part of the initialization into two subgoals. One subgoal was the specification and implementation of the double looped strip within a torus topology as described in the W&vdM paper. The other subgoal was the specification and implementation of the polarized photon trapped on this double looped strip. These two parts are now completely separate, but build on one another. Before we can continue with the initialization of the full electron model, however, we probably need to take step backwards and have a closer look at the polarized photon again. The reason for that is that we've specified the polarized photon in terms of classical EM theory i.e. Maxwells equations. In our brane model however, we treat the EM fields as emergent from the mechanics of the underlying brane. In our paper we've already formalized the relation between our brane model and classical EM theory quite extensively. Basically we reduce the magnetic field to a relativistic effect of the electric field and the electric field to a property of the amplitude of the brane itself. The relativistic effects are reduced to lateral excitations within the brane. Now, to initialize our full electron model properly, we need to implement a reverse mapping from the EM fields in our polarized photon implementation to actual excitations of our brane, both in the amplitude direction as well as the lateral directions. Could you please give a detailed analysis on how such a mapping could look like? 



Good, this is exactly the right move: test the inverted mapping in the **simplest possible 3D setting** before you put it into the toroidal electron.

I’ll treat this as a design + “instructions for Claude code” task and keep three things clearly separated:

1. **EM mode generator** – `tubular_photon_mode.py` (unchanged)
2. **EM → brane mapping** – new module
3. **Experiment script** – `polarized_photon_experiment.py`, based on `photon_3d_experiment.py`

I’ll first sketch the architecture, then give concrete function signatures and code skeletons.

---

## 1. Architecture overview

Goal: A 3D experiment where

* The brane grid and solver are **exactly the same style** as in `photon_3d_experiment.py`.
* The “photon” is **not** a Gaussian packet in x, but a **straight tubular mode** (constant along the propagation direction, localized only in the transverse cross-section).
* The EM field is generated purely in **EM space** by `tubular_photon_mode.py`.
* A separate module applies the **inverted mapping**:

  * (\mathbf E, \mathbf B ;\Rightarrow; \xi) (amplitude) and lateral velocities of the brane.

No torus, no electron, no double-loop here – just a straight waveguide mode.

### Files and responsibilities

1. **`tubular_photon_mode.py` (existing, unchanged)**

   * Generates a *polarized* photon mode in a straight tubular waveguide.
   * From the point of view of the new experiment, it is a **black box** that can give you (\mathbf E(\mathbf x)) and (\mathbf B(\mathbf x)) at (t=0) for arbitrary 3D points (\mathbf x) along a straight tube.

2. **`branesim/physics/em_to_brane_mapping.py` (new)**

   * Implements the **inverted dictionary**:

     * Convert (\mathbf E,\mathbf B) in SI units to:

       * amplitude in the 4th coordinate (component 3 of `state.positions`),
       * lateral velocities of brane points (components 0–2 of `state.velocities`).
   * Encodes the formulas we discussed (energy matching and Poynting-based velocity).

3. **`experiments/polarized_photon_experiment.py` (new)**

   * Copy of `photon_3d_experiment.py` with three key differences:

     1. It does **not** call `initialize_waveguide_wave_shape_3d` (no longitudinal Gaussian).
     2. It instead:

        * queries `tubular_photon_mode` for (\mathbf E,\mathbf B) on the brane grid,
        * calls the EM→brane mapping to set `state.positions` and `state.velocities`.
     3. It configures the mode to be **straight** (propagation along x-axis) and **tubular in cross-section only**.

---

## 2. EM → brane mapping module

### 2.1. Location and imports

Create a new file:

* `branesim/physics/em_to_brane_mapping.py`

Suggested imports:

```python
# branesim/physics/em_to_brane_mapping.py
import torch
from dataclasses import dataclass
from typing import Tuple
```

### 2.2. Small helper dataclass

This keeps ε, μ, ρ etc. together and avoids “magical” constants spreading.

```python
@dataclass
class EMMaterialParams:
    epsilon_eff: float  # effective permittivity (≈ ε0)
    mu_eff: float       # effective permeability (≈ μ0)
    rho_mass: float     # effective mass density of the brane continuum [kg/m³]
```

In the experiment we will set

[
\rho_{\text{mass}} = \rho_D = \frac{m_\text{point}}{h_\text{phys}^3}
]

for 3D.

### 2.3. EM energy and Poynting

```python
def compute_em_energy_and_poynting(
    E: torch.Tensor,
    B: torch.Tensor,
    params: EMMaterialParams,
) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    Args:
        E: (N, 3) electric field in SI units [V/m]
        B: (N, 3) magnetic field in SI units [T]
        params: epsilon_eff, mu_eff

    Returns:
        u: (N,) EM energy density [J/m³]
        S: (N, 3) Poynting vector [W/m²]
    """
    epsilon = params.epsilon_eff
    mu = params.mu_eff

    E2 = (E ** 2).sum(dim=-1)           # |E|²
    B2 = (B ** 2).sum(dim=-1)           # |B|²

    # u = 1/2 (εE² + B²/μ)
    u = 0.5 * (epsilon * E2 + B2 / mu)

    # S = 1/μ E × B
    S = torch.cross(E, B, dim=-1) / mu

    return u, S
```

### 2.4. Amplitude mapping (energy matching)

This is the inverted formula from the previous message:

[
\langle u_{\text{brane}}\rangle_t = \frac12 \rho_{\text{mass}}\omega^2 A^2
\quad\stackrel!=\quad
\langle u_{\text{EM}}\rangle_t
]

Ignoring time-averaging factors (they cancel if both are sinusoidal), we take:

[
A(\mathbf x) = \sqrt{\frac{2 u_{\text{EM}}(\mathbf x)}
{\rho_{\text{mass}} \omega^2}}
]

and optionally clamp to a chosen fraction of (h).

```python
def compute_brane_amplitude_from_em_energy(
    u_em: torch.Tensor,
    omega_phys: float,
    rho_mass: float,
    max_amplitude: float = None,
) -> torch.Tensor:
    """
    Map EM energy density to brane amplitude magnitude.

    Args:
        u_em: (N,) EM energy density in [J/m³]
        omega_phys: physical angular frequency [rad/s]
        rho_mass: mass density of brane [kg/m³]
        max_amplitude: optional maximum amplitude in meters; if not None,
                       A is clamped to this value.

    Returns:
        A: (N,) amplitude magnitude [m]
    """
    # A = sqrt(2 u / (rho * omega^2))
    denom = rho_mass * (omega_phys ** 2)
    # regularization to avoid dividing by zero in empty regions
    denom = max(denom, 1e-60)
    A = torch.sqrt(2.0 * u_em / denom)

    if max_amplitude is not None:
        A = torch.clamp(A, max=max_amplitude)

    return A
```

### 2.5. Full initialization function

This function does **everything** the experiment needs:

* takes EM fields in physical units at (t=0),
* computes energy density and Poynting,
* computes amplitude (A(\mathbf x)),
* writes amplitude into `state.positions[:, field_component]` (component 3),
* writes lateral velocities from Poynting direction into `state.velocities[:, :3]`.

We deliberately *don’t* use tube logic here – it’s purely local and EM-based.

```python
def initialize_brane_from_em_fields(
    *,
    state,
    grid,
    mapper,
    m_point_phys: float,
    h_phys: float,
    omega_phys: float,
    E_field_phys: torch.Tensor,
    B_field_phys: torch.Tensor,
    epsilon_eff: float,
    mu_eff: float,
    c_light: float,
    field_component: int = 3,
    max_amplitude_fraction_of_h: float = 0.1,
    velocity_clip_to_c: bool = True,
) -> None:
    """
    Initialize brane positions and velocities from EM fields.

    Args:
        state: BraneState (positions and velocities will be modified in-place)
        grid: BraneGrid (used only to get N and spacing)
        mapper: DimensionalMapper (for SI -> sim units)
        m_point_phys: point mass [kg]
        h_phys: lattice spacing [m]
        omega_phys: photon angular frequency [rad/s]
        E_field_phys: (N, 3) E-field in SI units at t=0
        B_field_phys: (N, 3) B-field in SI units at t=0
        epsilon_eff, mu_eff: effective EM constants (≈ ε0, μ0)
        c_light: speed of light [m/s]
        field_component: index of amplitude dimension in state.positions (3)
        max_amplitude_fraction_of_h: clamp amplitude to this fraction of h_phys
        velocity_clip_to_c: if True, clip |v| ≤ c
    """

    # Number of grid points for sanity
    N = grid.num_points
    assert E_field_phys.shape == (N, 3)
    assert B_field_phys.shape == (N, 3)

    # 1) Effective mass density of brane continuum
    #    rho_D = m_point / h^3 for 3D
    rho_mass = m_point_phys / (h_phys ** 3)

    params = EMMaterialParams(
        epsilon_eff=epsilon_eff,
        mu_eff=mu_eff,
        rho_mass=rho_mass,
    )

    # 2) Energy density and Poynting vector
    u_em, S_phys = compute_em_energy_and_poynting(E_field_phys, B_field_phys, params)

    # 3) Amplitude magnitude from energy matching
    max_amp_phys = max_amplitude_fraction_of_h * h_phys
    A_phys = compute_brane_amplitude_from_em_energy(
        u_em=u_em,
        omega_phys=omega_phys,
        rho_mass=rho_mass,
        max_amplitude=max_amp_phys,
    )
    # Convert amplitude [m] to sim units (length)
    A_sim = mapper.to_sim_length(A_phys)  # in your mapping this is usually identity

    # 4) Write amplitude into 4th coordinate (positions[:,3])
    #    We choose phase φ_EM = 0 everywhere for now → pure cosine snapshot.
    #    If tubular_photon_mode exposes a spatial phase field, that can be
    #    incorporated here later.
    with torch.no_grad():
        state.positions[:, field_component] = A_sim

    # 5) Lateral velocities from Poynting flow:
    #    v_energy = S / u  (energy transport velocity)
    #    In EM, for a null wave, |v_energy| ≈ c.
    #    We map this directly to velocity of brane points in 3D coordinates.
    #    Add small epsilon to avoid division by zero in empty regions.
    eps = 1e-30
    u_expanded = (u_em + eps).unsqueeze(-1)  # (N,1)
    v_energy_phys = S_phys / u_expanded      # (N,3), units m/s

    if velocity_clip_to_c:
        # Clip |v| ≤ c_light
        v_norm = torch.linalg.norm(v_energy_phys, dim=-1, keepdim=True)
        factor = torch.clamp(v_norm, max=c_light) / (v_norm + eps)
        v_energy_phys = v_energy_phys * factor

    # Convert to sim units (velocity)
    v_energy_sim = mapper.to_sim_velocity(v_energy_phys)

    with torch.no_grad():
        # Overwrite only the lateral components (x,y,z)
        state.velocities[:, 0:3] = v_energy_sim
```

Notes:

* For now I set the *phase* of the amplitude to “all cosine” (φ=0). That’s enough to test the inverted mapping qualitatively:

  * amplitude ∝ local EM energy density,
  * lateral velocities ∝ Poynting vector.
* Later, if `tubular_photon_mode` can provide a per-point phase φ(σ), you can replace the simple assignment with:

  ```python
  state.positions[:, field_component] = A_sim * torch.cos(phi)
  state.velocities[:, field_component] = A_sim * omega_sim * torch.sin(phi)
  ```

  but for the first test I’d keep it simple and let the lateral motion be the main test.

---

## 3. New experiment: `polarized_photon_experiment.py`

### 3.1. File placement and structure

Place next to `photon_3d_experiment.py`, e.g.:

* `experiments/polarized_photon_experiment.py`

Start by copying `photon_3d_experiment.py` and then apply *surgical* changes.

### 3.2. High-level step list

In `main()`:

1. **Reuse all physical setup** from `photon_3d_experiment.py`:

   * `PhysicalConstants()`
   * Compton-based `h_phys`
   * `m_point`, `k_spring`, `rest_length_phys`
   * `DimensionalMapper`, `h_sim`, `m_sim`, `k_sim`, `c_sim`
   * `dt_phys`, `dt_sim`, CFL logic
   * `nx, ny, nz` and `BraneGrid`
   * `BraneState`, `SpringForceComputer`, `VelocityVerletSolver`
   * `state.initialize_flat_configuration(h_sim)`
   * `state.set_fixed_boundaries()` etc.
2. **Replace the wave initialization**:

   * Do *not* call `initialize_waveguide_wave_shape_3d`.
   * Do *not* use Gaussian envelopes in x.
3. **Instead**:

   * Get all spatial points from the brane grid:

     ```python
     coords_sim = grid.get_spatial_coordinates()  # (N, 3) in sim units
     coords_phys = mapper.to_phys_length(coords_sim)  # (N,3) in meters
     ```
   * Use `tubular_photon_mode.py` to compute (\mathbf E(\mathbf x)) and (\mathbf B(\mathbf x)) at `t=0` for those physical coordinates, with:

     * propagation direction = +x,
     * tubular cross-section (no Gaussian envelope along x),
     * chosen polarization (e.g. circular or linear).
   * Call `initialize_brane_from_em_fields(...)` from `em_to_brane_mapping.py` to fill `state.positions` and `state.velocities`.
4. **Run the solver loop** exactly as in `photon_3d_experiment.py`:

   * same diagnostic outputs,
   * same visualization of slices and lateral distortion.

### 3.3. Concrete code skeleton for the new initialization

Here’s a concrete sketch of the modified part of `main()`.

I’ll show only the section that differs from `photon_3d_experiment.py`; everything else (constants, mapper, solver, visualization) should be kept as is.

```python
# experiments/polarized_photon_experiment.py

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
import numpy as np

from branesim.grid import BraneGrid, Dimensionality
from branesim.state import BraneState
from branesim.physics.spring_force import SpringForceComputer
from branesim.physics.solver import VelocityVerletSolver
from branesim.physics.dimensional_mapping import DimensionalMapper
from branesim.physics.constants import PhysicalConstants
from branesim.physics.em_to_brane_mapping import (
    initialize_brane_from_em_fields,
)

# Import EM mode generator (do not modify its file)
import tubular_photon_mode  # or whatever the actual module name is


def main():
    print("=" * 70)
    print("Polarized Photon Experiment - Tubular Mode, Straight Waveguide")
    print("=" * 70)

    constants = PhysicalConstants()
    print(f"  c = {constants.c:.6e} m/s")
    print(f"  lambda_C = {constants.lambda_C:.6e} m")

    # --- PHYSICAL AND SIMULATION SETUP (COPY FROM photon_3d_experiment) ---

    # Physical lattice spacing, Compton-based
    h_phys = constants.lambda_C / 10.0  # example: 10 points per lambda_C
    D = 3  # 3D brane

    # Point mass and derived density
    m_point = constants.m_e  # or whatever you used in photon_3d_experiment
    rho_D = m_point / (h_phys ** D)
    T_D = rho_D * constants.c**2
    rest_length_phys = 0.0 * h_phys

    c_wave = constants.c
    k_spring = T_D * (h_phys ** (D - 2))

    mapper = DimensionalMapper(
        h_phys=h_phys,
        c_light=constants.c,
        mass_reference=m_point,
    )

    h_sim = mapper.to_sim_length(h_phys)
    m_sim = mapper.to_sim_mass(m_point)
    k_sim = mapper.to_sim_spring_constant(k_spring)
    c_sim = mapper.to_sim_velocity(c_wave)

    # CFL time step
    cfl_factor = 0.2
    dt_phys = cfl_factor * h_phys / c_wave
    dt_sim = mapper.to_sim_time(dt_phys)

    # Domain and grid
    nx = 100
    ny = 100
    nz = 100
    grid = BraneGrid((nx, ny, nz), spacing=h_sim, dimensionality=Dimensionality.THREE_D)

    state = BraneState((nx, ny, nz), Dimensionality.THREE_D, device=torch.device('cuda' if torch.cuda.is_available() else 'cpu'), dtype=torch.float32)
    state.initialize_flat_configuration(h_sim)
    initial_positions = state.positions.clone()

    state.set_fixed_boundaries()

    physics = SpringForceComputer(k_sim, mapper.to_sim_length(rest_length_phys))
    solver = VelocityVerletSolver(dt_sim, m_sim, physics, grid)

    # --- NEW: BUILD POLARIZED TUBULAR PHOTON MODE ---

    # Photon parameters (physical)
    wavelength_phys = constants.lambda_C  # or some multiple
    omega_phys = 2.0 * np.pi * constants.c / wavelength_phys

    # Get coordinates of grid points in physical units
    coords_sim = grid.get_spatial_coordinates()           # (N, 3), sim units
    coords_phys = mapper.to_phys_length(coords_sim)       # (N, 3), meters

    # Use tubular_photon_mode to compute EM fields on a STRAIGHT tube
    #
    # IMPORTANT:
    #   - propagation direction: +x
    #   - no Gaussian envelope in x
    #   - tubular cross-section in y-z
    #
    # The actual API depends on how tubular_photon_mode.py was implemented.
    # Implement a small wrapper here that calls into that module and returns:
    #
    #   E_field_phys, B_field_phys: (N, 3) torch tensors in SI units
    #
    # Below is a placeholder; Claude should replace it with the actual call.
    E_field_phys, B_field_phys = tubular_photon_mode.compute_em_fields_on_points(
        coords_phys,
        wavelength=wavelength_phys,
        polarization="circular",    # or "linear"
        time=0.0,
    )

    # Ensure EM fields are torch tensors on correct device/dtype
    device = state.positions.device
    dtype = state.positions.dtype
    E_field_phys = E_field_phys.to(device=device, dtype=dtype)
    B_field_phys = B_field_phys.to(device=device, dtype=dtype)

    # --- NEW: MAP EM FIELDS → BRANE INITIAL STATE ---

    initialize_brane_from_em_fields(
        state=state,
        grid=grid,
        mapper=mapper,
        m_point_phys=m_point,
        h_phys=h_phys,
        omega_phys=omega_phys,
        E_field_phys=E_field_phys,
        B_field_phys=B_field_phys,
        epsilon_eff=constants.epsilon0,
        mu_eff=constants.mu0,
        c_light=constants.c,
        field_component=3,
        max_amplitude_fraction_of_h=0.1,
        velocity_clip_to_c=True,
    )

    # After this call:
    #   - state.positions[:,3] encodes amplitude (ξ) based on EM energy
    #   - state.velocities[:,0:3] encode lateral motion based on Poynting vector

    # Initialize accelerations and apply boundaries
    solver.initialize_accelerations(state)
    state.apply_fixed_boundaries()

    # --- EVERYTHING BELOW: reuse time loop and visualization from photon_3d_experiment ---
    # ... (copy measurement loop, plotting, etc. from photon_3d_experiment.py) ...


if __name__ == "__main__":
    main()
```

Key points:

* **No Gaussian in x**: the tubular mode lives along a straight waveguide (x-axis) and is localized only radially in y–z. That property is enforced inside `tubular_photon_mode.compute_em_fields_on_points`.
* The experiment script **does not** know anything about torus geometry or double loops.
* The mapping logic is fully contained in `em_to_brane_mapping.py`.

---

## 4. How this satisfies your constraints

* **Complexity reduction**: We only test the inverted mapping in a straight waveguide geometry. No torus, no double loop, no electron self-trapping.
* **Separation of concerns**:

  * `tubular_photon_mode.py` is a pure EM-mode generator (unchanged).
  * `em_to_brane_mapping.py` is *the* mathematical implementation of the inverted dictionary.
  * `polarized_photon_experiment.py` just wires these together and uses the existing brane solver, like `photon_3d_experiment.py`.
* **No “tube logic” leakage** into the core brane physics:

  * The mapping is local in (\mathbf x) and uses only EM quantities and the global parameters (m_\text{point}, h_\text{phys}, \epsilon, \mu, c).
  * The brane core (grid, solver, spring forces) remains unchanged.

If you want, the next step after you’ve wired this up and run it is:

* Check visually whether the **lateral motion** encodes the polarization correctly (e.g. circular polarization shows up as a rotating distortion pattern in the cross-section, linear polarization as an oscillation in one transverse direction).
* Once that looks sane, we can think about how to pull this into the torus/double-loop geometry for the electron initialization.
