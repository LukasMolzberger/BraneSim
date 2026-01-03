Here’s the **complete Laplace–Beltrami eigenvalue problem on the 2-sphere** (with the geometry made explicit, and the operator written both invariantly and in coordinates).

---

## 1) Geometry of the sphere (S^2)

Take the **unit sphere**
[
S^2={,x\in\mathbb R^3:|x|=1,},
]
parameterized by standard spherical angles
[
(\theta,\phi)\in(0,\pi)\times(0,2\pi),
]
with embedding
[
\mathbf r(\theta,\phi)=
(\sin\theta\cos\phi,\ \sin\theta\sin\phi,\ \cos\theta).
]

### Induced Riemannian metric

The induced line element is
[
ds^2 = d\theta^2+\sin^2\theta, d\phi^2.
]
So the metric tensor and its inverse are
[
g_{ab}=
\begin{pmatrix}
1 & 0\
0 & \sin^2\theta
\end{pmatrix},
\qquad
g^{ab}=
\begin{pmatrix}
1 & 0\
0 & \csc^2\theta
\end{pmatrix}.
]
Determinant and area element:
[
|g|=\sin^2\theta,\qquad \sqrt{|g|}=\sin\theta,\qquad
d\Omega=\sin\theta, d\theta, d\phi.
]

(For a sphere of radius (R): (ds^2=R^2(d\theta^2+\sin^2\theta,d\phi^2)) and (d\Omega_R=R^2\sin\theta,d\theta,d\phi).)

---

## 2) Laplace–Beltrami operator on (S^2)

### Coordinate-free definition

For a smooth scalar function (f) on a Riemannian manifold,
[
\Delta f = \operatorname{div}(\nabla f).
]

### Coordinate expression (general)

In local coordinates (x^a),
[
\Delta f
========

\frac{1}{\sqrt{|g|}},
\partial_a!\Big(\sqrt{|g|}, g^{ab}, \partial_b f\Big).
]

### Specializing to (S^2) in ((\theta,\phi))

Plugging in (g^{ab}) and (\sqrt{|g|}) gives
[
\boxed{
\Delta_{S^2} f(\theta,\phi)
===========================

\frac{1}{\sin\theta}\frac{\partial}{\partial\theta}
!\left(\sin\theta,\frac{\partial f}{\partial\theta}\right)
+
\frac{1}{\sin^2\theta}\frac{\partial^2 f}{\partial\phi^2}.
}
]

For radius (R), the operator rescales:
[
\Delta_{S^2_R} f = \frac{1}{R^2},\Delta_{S^2} f.
]

---

## 3) The eigenvalue equation (spherical harmonics)

The **Laplace–Beltrami eigenvalue problem** on the unit sphere is:

[
\boxed{
\Delta_{S^2} Y(\theta,\phi) = -\lambda, Y(\theta,\phi).
}
]

Its complete set of smooth eigenfunctions is given by the **spherical harmonics** (Y_\ell^m), with
[
\boxed{
\Delta_{S^2} Y_\ell^m(\theta,\phi)
==================================

-\ell(\ell+1),Y_\ell^m(\theta,\phi),
\qquad
\ell=0,1,2,\dots,\quad m=-\ell,\dots,\ell.
}
]

So the eigenvalues are
[
\lambda_\ell=\ell(\ell+1),
]
and each eigenspace has dimension (2\ell+1).

---

## 4) Explicit separated form and the associated ODE

One common explicit representation is
[
Y_\ell^m(\theta,\phi)=N_{\ell m},P_\ell^{m}(\cos\theta),e^{i m\phi},
]
where (P_\ell^m) are associated Legendre functions and (N_{\ell m}) a normalization constant.

If you insert (Y(\theta,\phi)=\Theta(\theta)\Phi(\phi)) into the eigenvalue equation, you get:

* (\Phi(\phi)=e^{i m\phi}), (m\in\mathbb Z),
* and (\Theta(\theta)) satisfies
  [
  \frac{1}{\sin\theta}\frac{d}{d\theta}!\left(\sin\theta,\frac{d\Theta}{d\theta}\right)
  -\frac{m^2}{\sin^2\theta}\Theta

-

\ell(\ell+1)\Theta
=0,
]
whose regular solutions are (\Theta(\theta)\propto P_\ell^m(\cos\theta)) with (\ell\ge |m|).

---

## 5) Orthonormality and completeness (the “spectral decomposition”)

With the area element (d\Omega=\sin\theta,d\theta,d\phi),
[
\int_{S^2} Y_\ell^m(\theta,\phi), \overline{Y_{\ell'}^{m'}(\theta,\phi)}, d\Omega
=================================================================================

\delta_{\ell\ell'}\delta_{mm'}.
]

Any square-integrable function (f\in L^2(S^2)) can be expanded as
[
f(\theta,\phi)=\sum_{\ell=0}^{\infty}\sum_{m=-\ell}^{\ell} a_{\ell m},Y_\ell^m(\theta,\phi),
]
with coefficients
[
a_{\ell m}=\int_{S^2} f(\theta,\phi),\overline{Y_\ell^m(\theta,\phi)}, d\Omega.
]

---

If you want, I can also write the same eigenvalue equation **in fully intrinsic index notation** (using covariant derivatives ( \nabla_a \nabla^a Y )) or adapt it to your brane-paper conventions (e.g., distinguishing “material” vs “embedded” coordinates).
