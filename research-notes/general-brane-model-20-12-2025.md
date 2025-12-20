Because I mixed two different things:

* the **exact identity** you get under the *flat-base / constant-normal / pure-normal-displacement* assumptions, and
* the **general variation formula** for an arbitrary (possibly curved) embedded brane, where extra terms appear.

Let me be precise.

---

## 1) When (\delta g_{ij} = \partial_i w,\partial_j w) is **exact**

Take a base embedding (\mathbf X_\star(x)) that is **flat** in (\mathbb R^4), so you can choose material coordinates such that

* (\partial_i \mathbf X_\star) are constant tangent vectors,
* the unit normal (\mathbf n) is **constant** (does not depend on (x)),
* and consider a **pure normal perturbation**:
  [
  \mathbf X(x,t)=\mathbf X_\star(x)+ w(x,t),\mathbf n.
  ]

Then
[
\partial_i \mathbf X = \partial_i\mathbf X_\star + (\partial_i w),\mathbf n,
]
and the induced metric is
[
g_{ij}=\partial_i\mathbf X\cdot\partial_j\mathbf X
= g^\star_{ij}
+(\partial_i w)(\partial_j w),(\mathbf n\cdot\mathbf n)
+(\partial_i\mathbf X_\star\cdot \mathbf n)\partial_j w
+(\partial_j\mathbf X_\star\cdot \mathbf n)\partial_i w.
]

But (\mathbf n\cdot\mathbf n=1) and (\partial_i\mathbf X_\star\cdot \mathbf n=0) (normal is orthogonal to tangents), so the cross terms vanish and you get:

[
\boxed{g_{ij}=g^\star_{ij}+\partial_i w,\partial_j w}
]

This is **not an approximation** under those assumptions. It is exact.

So my “(\approx)” was misleading in that flat-base setting.

---

## 2) Why it becomes an **approximation** in the general case

For a **general** embedded brane (curved base, non-constant normal, and/or perturbations that include tangential components), you must keep extra terms.

The key point is: if (\mathbf n=\mathbf n(x)), then
[
\partial_i(w,\mathbf n)= (\partial_i w),\mathbf n + w,\partial_i \mathbf n.
]
That extra (w,\partial_i \mathbf n) piece generates additional contributions to (g_{ij}) that are **linear in (w)** (and proportional to curvature).

In standard differential geometry language (codimension 1 for your (d=3) in (\mathbb R^4)), write a general small displacement as
[
\boldsymbol\xi = u^k,\partial_k \mathbf X_\star + w,\mathbf n,
]
(tangential part (u), normal part (w)). Then the metric variation has the well-known structure
[
\boxed{
\delta g_{ij}
= \nabla_i u_j + \nabla_j u_i ;-; 2w,b_{ij} ;+; \partial_i w,\partial_j w ;+; \text{higher order terms}
}
]
where

* (\nabla) is the covariant derivative compatible with (g^\star),
* (b_{ij}) is the second fundamental form (curvature) of the base embedding.

Now you can see exactly what happened:

* On a **flat** base: (b_{ij}=0), and if you take **pure normal** perturbations: (u=0), then the leading change is indeed (\delta g_{ij}=\partial_i w,\partial_j w) (and it’s exact if (\mathbf n) is constant).
* On a **curved** base: the term (-2w,b_{ij}) is present, so (\delta g) has a **linear-in-(w)** contribution.
* If you allow **tangential** perturbations: (\nabla_i u_j+\nabla_j u_i) appears (linear in (u)).

So in the general case I should have written:
[
\delta g_{ij} = -2w,b_{ij} + \partial_i w,\partial_j w + \dots
]
and explained which terms I was dropping.

---

## 3) What remains true (and addresses your core concern)

Nothing here makes (X^4) special.

The special direction is the **geometric normal** (\mathbf n(x)), which is defined from (\mathbf X_\star) itself. In a flat base configuration you can choose coordinates so (\mathbf n) coincides with the 4th axis, but that’s just a coordinate convenience.

And the “linear restoring” statement I was aiming for is **about the linearized equation of motion** around a pre-stressed base state, not about putting a linear term into the energy.

If you want, I can write the pre-stress linearization cleanly *without ever introducing (w) as “(X^4)”*, i.e. entirely in terms of (\boldsymbol\xi), tangential/normal decomposition, (g_{ij}), and the base stress tensor.
