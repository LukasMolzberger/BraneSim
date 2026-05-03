# Claude Code multi-agent setup for BraneSim

The repository root `README.md` is the canonical instruction file (TL;DR: read `principles.md` first).

## Specialist agents (`.claude/agents/`)

Each agent has a narrow scope and only the context it needs. Spawn them via the `Agent` tool or via the `/layer` slash command.

| Agent | Use for | Model |
|---|---|---|
| `principles-auditor` | check any diff against `principles.md` and `paper-v4/backbone.md` | sonnet |
| `physics-derivation` | derive equations bridging two adjacent physical layers | opus |
| `dispersion-analyst` | linear-regime ω(k), branch speeds, isotropy (level 1) | sonnet |
| `berry-validator` | Berry / Wilczek–Zee holonomy on isolated bands (level 3) | sonnet |
| `soliton-hunter` | search for stable baryon-like triplet modes (level 5) | sonnet |
| `contraction-channel` | extract emergent gravity-like contraction field (level 6) | sonnet |
| `simulation-engineer` | non-physics code work in `components/` and `orchestration/` | sonnet |
| `paper-writer` | LaTeX edits in `paper-v4/` enforcing backbone | sonnet |

## Slash commands (`.claude/commands/`)

- `/principles` — read the canonical principles file
- `/audit` — run the principles-auditor on the current diff or named files
- `/pipeline` — run the 4-component pipeline (`orchestration/run_pipeline.py`)
- `/layer <name> <task>` — route a task to the right specialist agent

## When in doubt

1. Read `principles.md` (root).
2. Read `paper-v4/backbone.md`.
3. Pick the agent whose scope matches the change.
4. Run `/audit` before declaring a physics-core change done.