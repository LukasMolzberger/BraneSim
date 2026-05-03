---
description: Run the 4-component pipeline (init → sim → viz → diag) using a JSON config. Defaults to orchestration/configs/local_extensive.json.
---

Run the BraneSim 4-component pipeline.

If `$ARGUMENTS` is empty, use `orchestration/configs/local_extensive.json`.
Otherwise, treat `$ARGUMENTS` as the path to a JSON config.

Execute:

```
PYTHONPATH=/Users/lukasmolzberger/PycharmProjects/BraneSim MPLBACKEND=Agg python3 -m orchestration.run_pipeline --config <config>
```

Report the resulting `pipeline_summary.json` and the path to the diagnostics summary.