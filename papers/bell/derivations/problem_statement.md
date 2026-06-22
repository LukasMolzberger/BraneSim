# Bell Bridge — Problem Statement

Strengthen the Bell bridge from "retrocausality is one consistent option" to "under the specific assumptions of this model, retrocausality is uniquely forced."

## The formal argument

Bell's theorem: empirical Bell-inequality violations require abandoning at least one of {determinism (D), local causality (L), measurement independence (MI)}.

This model commits to two of these as non-negotiable substrate facts (§2.2 `02_bell_constraint.tex`):
- Determinism: D = true (substrate is deterministic; all dynamics is root-finding ‖R‖ = 0 on the 4D block, no stochastic element).
- Locality: L = true (forces are nearest-neighbour spring contacts; information propagates at finite speed c_T/c_L; adding a superluminal channel would be a different model).

Given D = true and L = true, the assumption that must be relaxed is measurement independence (MI).

MI can be relaxed in two ways:
1. Superdeterminism — the hidden state λ is correlated with the settings through fine-tuned, lawless past data (settings effectively not freely chosen).
2. Retrocausality — the freely-chosen future setting back-propagates along the time-symmetric worldtube to the emission event, correlating λ with the setting via a uniform local action.

Superdeterminism is rejected on the mechanism-vs-fine-tuning discriminator (§2.5): it sources the correlation from lawless tuned past data rather than from the substrate's uniform local law. The separate give-up-L route — nonlocal hidden variables (Bohmian-type) — is excluded by the substrate's own ontology (local lattice, finite signal speed, no preferred foliation; §2.6).

CONCLUSION: Under D + L + local substrate, relaxing MI via retrocausality is the unique consistent completion, not merely preferred.

## Task

Make this argument rigorous and state it as a theorem with explicit premises:

**Theorem.** Let M be a deterministic (D), local (L_sub), measurement-independent (MI) completion of quantum mechanics. If the completion is local in the substrate sense (interaction range finite = lattice spacing a), then M must be retrocausal.

**Proof strategy:** (1) Show MI + D → not-L (from Bell). (2) Show local substrate → not-Bohmian (Bohmian requires nonlocal guidance equation; local substrate has no guidance equation). (3) Conclude retrocausality from the exclusive disjunction. State the no-signalling constraint explicitly and show it is compatible with retrocausality (operational no-signalling is preserved because the retrocausal influence propagates through the 4D worldtube, which respects the substrate's causal light cone on average).