# Shell-weight retuning for leading-order cubic isotropy

Sprint 1, subtask 1b of `paper-v4/validation_roadmap.md`.

---

## 1. Definition of shell weights

Use the same 26-neighbor cubic stencil as `components/simulation/grid.py`, but
replace the current effective shell weights by
$$
k_\delta = \frac{k_0\,w_s}{|\delta|^2},
\qquad
s\in\{\mathrm{I},\mathrm{II},\mathrm{III}\},
$$
where shell I is axial, shell II is face-diagonal, and shell III is
body-diagonal. With this convention, the current implementation is
$$
(w_{\rm I},w_{\rm II},w_{\rm III})=(1,1,1).
$$

The weighted second-moment tensor is
$$
T^{(2)}_{ab}(w)=S(w)\,\delta_{ab},
\qquad
S(w)=2w_{\rm I}+4w_{\rm II}+\frac{8}{3}w_{\rm III}.
$$

The weighted fourth-moment tensor can be written as
$$
T^{(4)}_{abcd}(w)
=J(w)\big(\delta_{ab}\delta_{cd}
        +\delta_{ac}\delta_{bd}
        +\delta_{ad}\delta_{bc}\big)
+H(w)Q_{abcd},
$$
where
$$
J(w)=w_{\rm II}+\frac{8}{9}w_{\rm III},
\qquad
H(w)=2w_{\rm I}-w_{\rm II}-\frac{16}{9}w_{\rm III}.
$$

The $H(w)Q_{abcd}$ term is the leading-order cubic anisotropy.

---

## 2. Isotropy manifold

Cubic isotropy at leading order is equivalent to
$$
\boxed{\;H(w)=0
\quad\Longleftrightarrow\quad
2w_{\rm I}=w_{\rm II}+\frac{16}{9}w_{\rm III}.\;}
$$

Equivalently,
$$
w_{\rm I}=\frac12 w_{\rm II}+\frac89 w_{\rm III},
$$
with $w_{\rm II}>0$ and $w_{\rm III}>0$ giving positive spring weights.

The current weights $(1,1,1)$ fail this condition because
$$
H(1,1,1)=2-1-\frac{16}{9}=-\frac79.
$$

---

## 3. Acoustic speeds on the isotropy manifold

For in-brane linear waves, the acoustic tensor is
$$
A_{abcd}
=\frac{k_0}{2a}
\left[\alpha\,T^{(4)}_{abcd}(w)
+(1-\alpha)S(w)\,\delta_{bd}\delta_{ac}\right].
$$

On the isotropy manifold $H=0$, the direction-independent branch speeds are
$$
\boxed{\;
c_T^2=\frac{k_0}{2a\rho_m}\left[\alpha J+(1-\alpha)S\right],
\qquad
c_L^2=\frac{k_0}{2a\rho_m}\left[3\alpha J+(1-\alpha)S\right].
\;}
$$

At $\alpha=1$, this gives $c_L/c_T=\sqrt{3}$ because the central-force Cauchy
relation gives $\lambda=\mu$ in the isotropic limit. At $\alpha\to0$, the
prestress geometric stiffness dominates both branches and $c_L/c_T\to1$.

---

## 4. Recommended representative

A conservative representative should:

1. satisfy exact leading-order isotropy,
2. preserve the current second-moment stiffness scale $S=26/3$, and
3. stay as close as possible to the current effective weights $(1,1,1)$.

Minimizing
$$
(w_{\rm I}-1)^2+(w_{\rm II}-1)^2+(w_{\rm III}-1)^2
$$
subject to $H=0$ and $S=26/3$ gives
$$
\boxed{\;
w_{\rm I}=\frac{431}{345}\approx1.2493,\qquad
w_{\rm II}=\frac{334}{345}\approx0.9681,\qquad
w_{\rm III}=\frac{99}{115}\approx0.8609.
\;}
$$

These weights satisfy
$$
S=\frac{26}{3},\qquad
J=\frac{598}{345},\qquad
H=0.
$$

In terms of the code's per-neighbor `neighbor_weights = k_\delta/k_0`, this
would be
$$
\text{axis: }\frac{431}{345}\approx1.2493,\qquad
\text{face diagonal: }\frac{167}{345}\approx0.4841,\qquad
\text{body diagonal: }\frac{33}{115}\approx0.2870.
$$

The current code weights are axis $1$, face diagonal $1/2$, body diagonal
$1/3$.

---

## 5. Impact at the default BraneSim parameters

For `spring_constant = 1`, `spacing = 1`, `mass_density = 1`, and
`alpha = 0.2`, the recommended isotropic weights predict
$$
c_T^2=\frac{2093}{575}\approx3.6400,
\qquad
c_L^2=\frac{6877}{1725}\approx3.9867,
$$
so
$$
\boxed{\;
c_T\approx1.908,\qquad
c_L\approx1.997,\qquad
c_L/c_T\approx1.047.
\;}
$$

For comparison, the current unretuned weights predict
$$
c_T\approx1.912,\qquad c_L\approx1.989
$$
along a cubic axis at the same $\alpha$. The retuning therefore changes the
default branch speeds only mildly while removing the leading-order cubic
anisotropy.

---

## 6. Recommendation

Use the representative above for an explicit retuned-lattice experiment, but
do not silently replace the default model. The next dispersion sprint should
run both:

- **current weights** `(1,1,1)`: verifies the anisotropy predicted by
  `lattice_to_continuum.md`;
- **retuned weights** `(431/345,334/345,99/115)`: verifies that the same force
  code can recover a leading-order isotropic acoustic tensor when the shell
  moments are tuned.

If the retuned run matches the isotropic prediction and the current run matches
the anisotropic prediction, the numerical pipeline is validating the substrate
mechanics rather than hiding a preferred direction by analysis choices.
