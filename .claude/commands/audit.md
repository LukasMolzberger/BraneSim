---
description: Audit the current diff or named files against principles.md non-negotiables. Spawns the principles-auditor agent.
---

Run the `principles-auditor` agent against the current changes.

If `$ARGUMENTS` is empty, audit `git diff` against the working tree and `git status -s`.
If `$ARGUMENTS` names files or directories, audit those.

Tell the agent which scope to use, then return its punch list verbatim.