# AWS Runbook — Breather Soliton Sweep

This runbook covers running the time-periodic breather eigen-BVP bracket
(`solve_breather(mode="topological")`) on AWS EC2.  The job is **CPU-bound
and tiny in memory** (~0.27 GB working set at 64^3); it does NOT use the
block-solver path (`branesim.run_experiment` / JFNK BVP).

---

## (a) Account Prerequisites

Pulled from `launch_branesim_job.py` and `orchestration/aws/README.md`:

| Item | Required value / action |
|---|---|
| AWS CLI | Installed and authenticated (`aws sts get-caller-identity` succeeds) |
| Region | `us-east-1` (or your preferred region; use consistently) |
| AMI | Ubuntu 22.04 LTS x86_64, e.g. `ami-0c7217cdde317cfec` (us-east-1); pick the current canonical Ubuntu AMI for your region |
| S3 bucket | `my-branesim-bucket` — create if absent; same region as instance |
| IAM instance profile | `BraneSimEc2Profile` — must have the following inline policy: |

```json
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Action": ["s3:GetObject", "s3:PutObject", "s3:ListBucket"],
      "Resource": [
        "arn:aws:s3:::my-branesim-bucket",
        "arn:aws:s3:::my-branesim-bucket/*"
      ]
    }
  ]
}
```

| Item | Required value / action |
|---|---|
| Security group | Allow outbound HTTPS (443) for apt/pip/S3; no inbound required for headless job |
| Subnet | Any subnet with outbound internet (NAT or IGW); note `subnet-xxxxxxxx` |
| Key pair | Optional; omit `--key-name` for headless jobs |

---

## (b) Exact commands

### Step 1 — Build the project archive

```bash
bash /path/to/BraneSim/orchestration/aws/package_project.sh /tmp/branesim-project.tar.gz
```

Verify the breather entrypoint is included:

```bash
tar -tzf /tmp/branesim-project.tar.gz | grep breather_sweep
# Expected output:
# ./branesim/experiments/breather_sweep.py
# ./branesim/experiments/__init__.py
```

### Step 2 — Launch: sequential (all 6 bracket points on one instance)

The 6 AWS bracket points are independent; the simplest approach is to run
them sequentially in one remote command, writing each result as it completes
so a partial run still yields data.

```bash
python /path/to/BraneSim/orchestration/aws/launch_branesim_job.py \
  --region us-east-1 \
  --ami-id ami-0c7217cdde317cfec \
  --instance-type c7i.4xlarge \
  --subnet-id subnet-xxxxxxxx \
  --security-group-id sg-xxxxxxxx \
  --iam-instance-profile BraneSimEc2Profile \
  --s3-bucket my-branesim-bucket \
  --s3-prefix breather-runs \
  --job-name breather-aws-bracket \
  --project-archive /tmp/branesim-project.tar.gz \
  --volume-size-gb 30 \
  --remote-command 'branesim-breather --aws-bracket --output-dir "$BRANESIM_RESULTS_DIR"'
```

### Step 2 (alternative) — Launch: 6 parallel processes on a fat-core instance

Run each bracket point in a separate background **process**, then `wait`. This is
the *only* effective parallelism here (see **Threading** below) and cuts walltime to
the slowest single point. Each process is pinned to 2 threads via env vars: the force
kernel is a Python neighbor-loop + `np.add.at` scatter (NOT BLAS), so more threads per
process barely help — and *uncapped*, each of the 6 processes' numpy/OpenBLAS would
grab all cores → 6× oversubscription → slower than sequential. 6 processes × 2 threads
= 12 fits a `c7i.4xlarge` (16 vCPU) with headroom. Thread-pinning also fixes the FP
reduction order, so the Floquet radii are reproducible run-to-run.

```bash
python /path/to/BraneSim/orchestration/aws/launch_branesim_job.py \
  --region us-east-1 \
  --ami-id ami-0c7217cdde317cfec \
  --instance-type c7i.4xlarge \
  --subnet-id subnet-xxxxxxxx \
  --security-group-id sg-xxxxxxxx \
  --iam-instance-profile BraneSimEc2Profile \
  --s3-bucket my-branesim-bucket \
  --s3-prefix breather-runs \
  --job-name breather-aws-bracket-parallel \
  --project-archive /tmp/branesim-project.tar.gz \
  --volume-size-gb 30 \
  --remote-command 'export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2; for i in 0 1 2 3 4 5; do branesim-breather --aws-idx $i --output-dir "$BRANESIM_RESULTS_DIR" & done; wait'
```

The bracket is a **SIZE-scan at fixed alpha=0.85** (re-scoped 2026-06-04: the local
32^3 study found Floquet rho flat at ~3.5 across alpha=0.80/0.85/0.90, so alpha is
NOT the stabilizer; the remaining lever is soliton size). `--aws-idx` runs one entry:
- 0 = 32^3, w=4   (local reference, rho=3.47)
- 1 = 48^3, w=6
- 2 = 48^3, w=7
- 3 = 64^3, w=7   (grid-control vs idx 2 — physical-size vs grid-resolution)
- 4 = 64^3, w=9
- 5 = 64^3, w=10  (fattest)

Decisive read: does rho fall monotonically with w toward 1 (stable)? A flat rho≈3.5
means the breathing ansatz is intrinsically unstable (→ pivot to the static
soliton + iso-rotation rest state).

### Step 3 — Retrieve results

The launcher prints `results_s3: s3://my-branesim-bucket/breather-runs/<job-id>/results`.

```bash
bash /path/to/BraneSim/orchestration/aws/fetch_results.sh \
  s3://my-branesim-bucket/breather-runs/<job-id>/results \
  ./breather-out/<job-id>
```

The retrieved directory will contain:
- `breather_sweep.csv` — one row per bracket point (appended as each completes)
- `<label>_result.json` — per-run diagnostics (harmonics, Floquet detail)
- `<label>_worldtube.npz` — converged worldtube (float32 slices + metadata)
- `user-data.log` — full bootstrap + solver stdout/stderr

The S3 sync in `ec2_user_data.sh.tmpl` uses `aws s3 sync --delete`, which
covers all file types including CSV/JSON/NPZ.  No path adjustments needed:
`$BRANESIM_RESULTS_DIR` is the S3-synced results directory.

---

## (c) Instance sizing — breather vs block solver

The breather solver is CPU-bound and low-memory (no JFNK Krylov subspace of
30 vectors; just one Newton-LGMRES pass per outer iteration).

Working-set estimate at 64^3, P=16, m_ambient=4:
```
state = P * N^3 * m * 8 bytes  (float64)
      = 16 * 262144 * 4 * 8   = 134 MB  (one state copy)
LGMRES inner_maxiter=2000: holds ~2000 Krylov vectors of size ~134 MB/4 = 33 MB each
    → peak ~67 GB  (worst case; LGMRES restarts limit actual footprint)
Practical peak at inner_maxiter=2000: ~4-8 GB  (LGMRES with restart=20 default)
```

Recommended instance: `c7i.4xlarge` (16 vCPU, 32 GB RAM, $0.68/hr on-demand).
For the parallel 6-process launch: `c7i.8xlarge` (32 vCPU, 64 GB RAM, $1.36/hr).

Do NOT use `r7i.*` (memory-optimized, expensive) for the breather sweep.
Memory-optimized instances are needed for the 4D block BVP (30-vector Krylov
subspace, see DEPLOYMENT.md §1), not for the breather eigen-BVP.

### Walltime estimates

Local timing baseline (from test-runs/sprint4b_skyrme_corrected/):
- 32^3, max_iter=100, inner_maxiter=1000: ~935 s/point
- 32^3 → 48^3 cost scale: (48/32)^3 = 3.375x in nodes per matvec
- 32^3 → 64^3 cost scale: (64/32)^3 = 8.0x in nodes per matvec
- inner_maxiter 1000 → 2000: ~2x LGMRES inner work per outer Newton step

| Grid | inner_maxiter | Estimated walltime/point (1 vCPU) |
|---|---|---|
| 32^3 | 1000 | ~935 s (~16 min) [measured] |
| 48^3 | 2000 | ~935 * 3.375 * 2 = ~6,300 s (~1.75 hr) |
| 64^3 | 2000 | ~935 * 8.0 * 2 = ~15,000 s (~4.2 hr) |

**Threading.** The hot path (`spacelike_force`) is a Python neighbor-loop with an
`np.add.at` scatter — not dense linear algebra — so per-process BLAS/OMP threading
buys almost nothing (~1.0–1.3×). The real speedup is **process-level across the 6
independent bracket points**: one process per point, each pinned to 1–2 threads
(`OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `NUMEXPR_NUM_THREADS`),
with the instance sized to (#parallel points) × (threads/point). Do NOT leave threads
uncapped under the parallel launch — 6 processes each grabbing all cores oversubscribes
and runs slower than sequential.

Sequential run (6 points, all cores to one point at a time): ~7–12 hr on `c7i.4xlarge`.
Parallel run (6 processes × 2 threads): ~2–4 hr total on `c7i.4xlarge` (limited by the
slowest 64^3 point; threading gives little, so don't expect < 2 hr).

### Cost estimate

| Mode | Instance | Walltime | On-demand cost |
|---|---|---|---|
| Sequential, 6 pts | `c7i.4xlarge` ($0.68/hr) | 10 hr | ~$6.80 |
| Parallel, 6 pts | `c7i.4xlarge` ($0.68/hr) | ~3 hr | ~$2.00 |
| Spot (parallel) | `c7i.4xlarge` spot (~$0.25/hr) | ~3 hr | ~$0.75 |

Add `--spot` to the launch command for spot pricing (safe for this workload:
if preempted, the EXIT trap still syncs whatever completed runs are in
`$BRANESIM_RESULTS_DIR` before shutdown).

---

## (d) S3 sync and self-terminate confirmation

- `ec2_user_data.sh.tmpl` sets `trap 'cleanup_and_exit $?' EXIT`.
- `cleanup_and_exit` calls `aws s3 sync "${RESULTS_DIR}" "${S3_RESULTS_URI}" --delete`
  unconditionally (even on failure / preemption).
- `instance-initiated-shutdown-behavior=terminate` is set in `launch_branesim_job.py`.
- Root EBS `DeleteOnTermination=true`.

The sync covers ALL files in `$BRANESIM_RESULTS_DIR` (CSV, JSON, NPZ, log).
Partial runs are safe: each bracket point appends its CSV row and writes its
JSON/NPZ independently inside `run_one()`, so results from completed points
are preserved even if a later point fails or the instance is preempted.

---

## Bracket summary

```
AWS_BRACKET (6 points, alpha in {0.80, 0.85, 0.90}, u0=10):
  Index  Label                    Grid   w    alpha
  0      aws_a0p80_u10_w8_64      64^3   8.0  0.80   box/w = 8.0 >= 6
  1      aws_a0p85_u10_w8_64      64^3   8.0  0.85   box/w = 8.0 >= 6
  2      aws_a0p90_u10_w8_64      64^3   8.0  0.90   box/w = 8.0 >= 6
  3      aws_a0p80_u10_w7_48      48^3   7.0  0.80   box/w = 6.86 >= 6
  4      aws_a0p85_u10_w7_48      48^3   7.0  0.85   box/w = 6.86 >= 6
  5      aws_a0p90_u10_w7_48      48^3   7.0  0.90   box/w = 6.86 >= 6

BreatherOpts: tol=1e-6, max_iter=200, inner_maxiter=2000, method="lgmres"
P=16, k_s=rho=a=1, m_ambient=4, analytic band_top override (no dense eigensystem)
```
