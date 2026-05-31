# BraneSim — Cloud Deployment (AWS)

The 4D block/retrocausal solve is **memory-bound**, not FLOP-bound. This guide
covers running `branesim` experiments on AWS EC2 with enough RAM, using the
cost-safe scaffolding in `orchestration/aws/`.

Backend decision (2026-05-31): **high-memory CPU, numpy/scipy** — no GPU rewrite.
JFNK needs RAM more than FLOPs; high-memory CPU instances (up to multi-TB) fit
the working set without sharding the validated numpy core. GPU is a documented
future option (§5).

---

## 1. Memory sizing — the number that picks the instance

One world-volume state vector is

```
bytes_per_vector = (n_slices + 1) · N^dim · m_ambient · 8        (float64)
```

JFNK / GMRES holds roughly **~30** such vectors (Krylov subspace + iterates +
warm start). So plan for `~30 × bytes_per_vector` of RAM, plus headroom.

`branesim.run_experiment` prints this estimate at startup
(`JFNK working-set estimate ~X GB`). Sizing table for the canonical `m=4`:

| spatial `N³` | `n_slices` | GB / vector | ~working set (×30) | suggested instance (RAM) |
|---|---|---|---|---|
| 32³  | 128 | 0.13 | ~4 GB   | `r7i.xlarge` (32 GB) |
| 48³  | 256 | 0.91 | ~27 GB  | `r7i.2xlarge` (64 GB) |
| 64³  | 256 | 2.16 | ~65 GB  | `r7i.4xlarge` (128 GB) ← launcher default |
| 96³  | 512 | 14.5 | ~436 GB | `r7i.16xlarge` (512 GB) |
| 128³ | 512 | 34.4 | ~1.0 TB | `x2iedn.16xlarge` (1 TB) |
| 128³ | 1024| 68.7 | ~2.1 TB | `x2iedn.32xlarge` (2 TB) |

Instance families (memory-optimized, x86, numpy-friendly):
- **`r7i.*`** — up to `r7i.48xlarge` = 1.5 TB. General memory-optimized default.
- **`x2iedn.*`** — up to 4 TB. For large blocks.
- **`u-*` (High Memory)** — 6–24 TB. Only for extreme worldtube extents.

Reduce memory by: smaller `N` (spatial), fewer `n_slices`, or (future) the
matrix-free streaming/slab storage in ARCHITECTURE.md D4.

> Practical note: today's deployable BVP path is **Dirichlet**, which is
> ill-conditioned (`κ≈1e14`) and makes JFNK slow (minutes even for tiny grids).
> The well-posed **chiral** BC (`κ=1`) is still broken (ARCHITECTURE.md D2); once
> fixed, block solves become fast. For now, `mode="ivp"` is the cheap smoke and
> `bvp_dirichlet` is a correctness check, not a performance target.

---

## 2. Package the project

`pip install -e .` works because of the repo `pyproject.toml` (package
`branesim`, deps `numpy`+`scipy`). Build the upload archive:

```bash
orchestration/aws/package_project.sh /tmp/branesim-project.tar.gz
```

(Excludes `.git`, `test-runs/`, caches.)

## 3. Launch a job (uploads to S3 → EC2 bootstrap → run → sync → self-terminate)

```bash
python orchestration/aws/launch_branesim_job.py \
  --region us-east-1 \
  --ami-id ami-xxxxxxxx \
  --instance-type r7i.4xlarge \
  --subnet-id subnet-xxxxxxxx \
  --security-group-id sg-xxxxxxxx \
  --iam-instance-profile BraneSimEc2Profile \
  --s3-bucket my-branesim-bucket --s3-prefix block-runs \
  --project-archive /tmp/branesim-project.tar.gz \
  --volume-size-gb 200
```

Default remote command (override with `--remote-command`):

```bash
python -m branesim.run_experiment \
  --config orchestration/configs/branesim_ivp_smoke.json \
  --output-dir "$BRANESIM_RESULTS_DIR"
```

For a block solve, point `--config` at `branesim_bvp_dirichlet.json` (or your own).

## 4. Retrieve results

```bash
orchestration/aws/fetch_results.sh s3://my-branesim-bucket/block-runs/<job-id>/results ./out
```

Outputs per run: `worldvolume.zip` (solved slices + `manifest.json` with
`solver_report`), `summary.json`, and `user-data.log`.

## 5. Cost safety (built into the scaffolding)

- `instance-initiated-shutdown-behavior=terminate`; user-data always syncs +
  shuts down on the EXIT trap (even on failure).
- Root EBS `DeleteOnTermination=true`.
- `--spot` for one-time spot instances (big savings for long memory-bound runs).
- Size `--volume-size-gb` for the archive + result `worldvolume.zip`s (these can
  be large; consider writing sparse slice strides for big blocks).

## 6. GPU / acceleration (future, not enabled)

Deferred by the 2026-05-31 backend decision. If runs become FLOP-bound (not
memory-bound), a GPU backend would mean: porting `branesim/core` array ops to
JAX or CuPy behind a thin indirection (ARCHITECTURE.md D4), and targeting
`p4d`/`p5`. The catch is GPU memory (40–80 GB/GPU) caps the block size far below
high-memory CPU, so multi-GPU sharding of the residual would be required. Revisit
only with real memory-pressure/timing data from CPU runs.

---

## Prerequisites recap
- Local: AWS CLI authenticated; `package_project.sh` run.
- AWS: subnet / security group / Ubuntu AMI / instance profile with
  `s3:GetObject,PutObject,ListBucket`.
- The bootstrap (`ec2_user_data.sh.tmpl`) installs python+ffmpeg+awscli, makes a
  venv, `pip install -r requirements.txt` + `pip install -e .`, runs the command,
  syncs to S3 (and optional `rclone` Google Drive), and shuts down.
