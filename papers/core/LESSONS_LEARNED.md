# Lessons Learned

A durable record of mistakes not to repeat and results we trust, so the knowledge
survives without being expensively re-encoded in tests.

## Methodology — mistakes we kept making

1. **Clamp the pretension (periodic spatial BC).** The substrate has prestress
   `α<1`: bonds are stretched (rest length `αa < a`), so the vacuum is under
   tension. With **open/free boundaries** an unpinned tensioned lattice just
   relaxes — it contracts globally toward `αa`. The tell: an `A=0` control (no
   soliton) collapses identically. Use periodic BC for every nonlinear run; the
   tensioned square lattice is then a stable vacuum (`A=0 → E_excess ≡ 0`).
2. **The discrete grid forbids Derrick collapse** (the core can't shrink below `a`)
   — but it does **not** prevent *unwinding* or *spreading*. "It can't collapse"
   was true and irrelevant to why solitons failed.
3. **Never state the consequence of a simplifying assumption as a no-go theorem.**
   The load-bearing nonlinearity is the Pythagorean link length (BACKBONE #3). Any
   analysis that linearizes it away — most importantly the **bare-vacuum
   Bloch/Hessian** `M_μ^{AB} = (1−α)δ^{AB} + α ê_μ ê_μ` — is only the quadratic
   tangent operator about a straight, axis-aligned background, and cannot decide
   questions that live *in* the nonlinearity. Concretely: the bare vacuum has no
   SU(3), but that is **not** a no-go for T1 — it only says the trivial background
   is the wrong place to look. The physically relevant object is the fluctuation
   Hessian about a finite-amplitude nonlinear background `R̄`,
   `H_{nμ}^{AB} = σ_μ k_μ[(1−r_μ/L̄)δ^{AB} + (r_μ/L̄) Q̂^A Q̂^B]`, a
   background-dependent frame tensor whose transverse sector rotates over the
   background and over `k`. Rule of thumb: **if a negative result depends on
   dropping the Euclidean square root, it is an artifact of the approximation, not
   a property of the substrate.**

