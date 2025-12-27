### Local (U(1)) phase (gauge) transformation

Let (\alpha(x)) be a real-valued function on spacetime and (q) the charge/coupling of the Dirac field. Define the local phase rotation
[
\psi(x);\mapsto;\psi'(x)=e^{-i q \alpha(x)},\psi(x),
\qquad
\bar\psi(x);\mapsto;\bar\psi'(x)=\bar\psi(x),e^{+i q \alpha(x)}.
]

To keep the theory invariant, the electromagnetic 4-potential transforms as
[
A_\mu(x);\mapsto;A'*\mu(x)=A*\mu(x)+\partial_\mu \alpha(x).
]

With the covariant derivative
[
D_\mu := \partial_\mu + i q A_\mu,
]
one has the key covariance property
[
D'*\mu \psi'(x)=e^{-i q\alpha(x)},D*\mu \psi(x),
]
which is what makes the Lagrangian locally phase symmetric.

---

### Complete locally phase-symmetric Lagrangian (QED)

**Compact form**
[
\boxed{;\mathcal L
= \bar\psi\left(i\gamma^\mu D_\mu - m\right)\psi
;-;\frac14,F_{\mu\nu}F^{\mu\nu};}
]
with
[
F_{\mu\nu}:=\partial_\mu A_\nu-\partial_\nu A_\mu.
]

**Expanded into “Dirac term + coupling term + EM term”**
[
\boxed{;\mathcal L
==================

\underbrace{\bar\psi\left(i\gamma^\mu\partial_\mu-m\right)\psi}*{\text{Dirac term}}
;+;
\underbrace{\left(-,q,\bar\psi\gamma^\mu A*\mu\psi\right)}*{\text{coupling term}}
;+;
\underbrace{\left(-\frac14,F*{\mu\nu}F^{\mu\nu}\right)}_{\text{electromagnetic term}};}
]

This Lagrangian is invariant under the local transformations above (the (F_{\mu\nu}F^{\mu\nu}) part is invariant because (F_{\mu\nu}) itself is unchanged by (A_\mu\mapsto A_\mu+\partial_\mu\alpha)).

---

### Dictionary of variables (with types)

* (x^\mu): spacetime point (4 real coordinates), (\mu=0,1,2,3).
* (\partial_\mu := \frac{\partial}{\partial x^\mu}): spacetime derivative operator.
* (\alpha(x)): **real scalar function** on spacetime (gauge parameter), (\alpha:\mathbb R^{1,3}\to\mathbb R).
* (\psi(x)): **Dirac spinor field**, complex 4-component column vector (\psi:\mathbb R^{1,3}\to\mathbb C^{4}).
* (\bar\psi(x)): **Dirac adjoint**, (\bar\psi:=\psi^\dagger\gamma^0) (row spinor).
* (\gamma^\mu): **Dirac gamma matrices**, (4\times 4) complex matrices satisfying
  [
  {\gamma^\mu,\gamma^\nu}=2\eta^{\mu\nu}\mathbf 1_{4}.
  ]
* (\eta^{\mu\nu}): **Minkowski metric** (choose a convention, e.g. (\mathrm{diag}(+1,-1,-1,-1))).
* (m): **mass parameter** of the Dirac field (real constant).
* (q): **(U(1)) charge / coupling constant** of the Dirac field (real constant; for the electron (q=-e)).
* (A_\mu(x)): **electromagnetic 4-potential**, real 4-vector field (A:\mathbb R^{1,3}\to\mathbb R^{4}).
* (D_\mu): **gauge-covariant derivative**, (D_\mu=\partial_\mu+i q A_\mu).
* (F_{\mu\nu}(x)): **field strength tensor**, antisymmetric rank-2 tensor field,
  [
  F_{\mu\nu}=\partial_\mu A_\nu-\partial_\nu A_\mu.
  ]
* (\mathcal L): **Lagrangian density**, a real scalar density on spacetime.

If you want, I can also write the same Lagrangian with explicit (\hbar) and (c) factors (SI-friendly form), but the above is the standard natural-units presentation used in relativistic QFT.
