# Bell Bridge — Problem Statement

Strengthen the Bell bridge from "retrocausality is one consistent option" to "under the specific assumptions of this model, retrocausality is uniquely forced."

## The formal argument

Bell's theorem: empirical Bell-inequality violations require abandoning at least one of {local causality (L), determinism (D), measurement independence (MI)}.

This model commits to:
- Determinism: D = true (substrate is deterministic; all dynamics is root-finding ‖R‖ = 0 on the 4D block, no stochastic element).
- Measurement independence: MI = true (no superdeterminism — the model does not rely on fine-tuning of initial conditions correlated with measurement settings; the substrate's initial state is not correlated with measurement settings by assumption).

Given D = true and MI = true, the only option remaining is to abandon local causality (L).

Abandoning L allows two sub-options:
1. Nonlocal hidden variables (Bohmian-type).
2. Retrocausality.

This model is built on a local substrate (forces propagate at finite speed through the lattice; the Verlet stencil is purely local). Therefore nonlocal hidden variables are excluded by the substrate's own ontology.

CONCLUSION: Under D + MI + local substrate, retrocausality is uniquely forced, not merely preferred.

## Task

Make this argument rigorous and state it as a theorem with explicit premises:

**Theorem.** Let M be a deterministic (D), local (L_sub), measurement-independent (MI) completion of quantum mechanics. If the completion is local in the substrate sense (interaction range finite = lattice spacing a), then M must be retrocausal.

**Proof strategy:** (1) Show MI + D → not-L (from Bell). (2) Show local substrate → not-Bohmian (Bohmian requires nonlocal guidance equation; local substrate has no guidance equation). (3) Conclude retrocausality from the exclusive disjunction. State the no-signalling constraint explicitly and show it is compatible with retrocausality (operational no-signalling is preserved because the retrocausal influence propagates through the 4D worldtube, which respects the substrate's causal light cone on average).