---
description: Pick the right specialist agent for a physical-layer task. Usage: /layer <layer-name> <task description>
---

Map the layer name in `$ARGUMENTS` to the right BraneSim agent and dispatch the task.

Layer routing:
- `dispersion`, `linear`, `branch`, `isotropy` → `dispersion-analyst`
- `berry`, `holonomy`, `gauge`, `wilczek-zee`, `wz`, `triplet-phase` → `berry-validator`
- `soliton`, `baryon`, `confinement`, `localized`, `particle` → `soliton-hunter`
- `contraction`, `gravity`, `newtonian`, `phi`, `prestretch` → `contraction-channel`
- `derivation`, `math`, `lagrangian`, `eom`, `bridge` → `physics-derivation`
- `code`, `solver`, `pipeline`, `refactor`, `integrator` → `simulation-engineer`
- `paper`, `latex`, `tex`, `manuscript` → `paper-writer`
- `audit`, `principles`, `compliance` → `principles-auditor`

If the layer keyword is ambiguous, ask the user which one. Then spawn the agent with the rest of `$ARGUMENTS` as its task.