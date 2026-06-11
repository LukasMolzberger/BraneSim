# AWS Runbook — U(1) Vortex Seed (render + full diagnostics) at scale

This runbook covers running the U(1) spherical-harmonic vortex experiment
(`branesim.experiments.vortex_seed_render`) on AWS EC2 at higher spatial
resolution than the local 48³ pre-test.

**What this job is.** Inject the Y₁¹ vortex seed → **relax it along the
solver-iteration axis with the rotating-frame-periodic BC** → render both the
seed and the relaxed state (3D volume + 2D slice movies) → run the full 8-device
diagnostic suite (energy, confinement, winding, Berry, EM `A_μ`/`F_μν`,
per-colour SU(3), spectra, binding_probe) → write per-output `.md` docs. It produces
`iter_0000/` (the seed) and `iter_0015/` (after the JFNK relaxation), each
self-contained and fully diagnosed, under `$BRANESIM_RESULTS_DIR`.

**The science.** The relaxation uses the **rotating-frame-periodic** BC
(`PeriodicBC`): a closed cyclic time loop, all slices free, well-conditioned
(cond ~1e3–1e4 vs the old Dirichlet two-time κ~1e14 that froze the solve). The
brane genuinely relaxes toward ‖R‖=0 (residual drops ~100×) while the carrier
winding is preserved. **The key question this run answers: does the U(1) vortex
worldtube BIND (stay localized, winding intact) or RADIATE toward vacuum?** —
read off from the energy / confinement / winding diagnostics of `iter_0015/`.

**Cost profile.** Two parts: (a) the JFNK relaxation is **memory-bound** — it
holds ~`inner_maxiter` Krylov vectors of the full worldvolume; (b) the rendering
is **CPU-bound** (matplotlib 3D volume). Size for the Krylov working set:

```
GB per worldvolume vector = (n_slices+1)·N³·m·8 / 1e9        (m=4, float64)
JFNK working set ≈ inner_maxiter · (GB per vector)            (default inner=40)
```

To run **seed-only** (render + diagnostics, no solve, a few GB RAM), set
`BRANESIM_VORTEX_RELAX=0`. To run a fast **diagnostics-only** pass (relax +
diagnostics, no CPU-heavy movie render — e.g. an eigensolve pre-test), set
`BRANESIM_VORTEX_RENDER=0`. Both knobs are env vars on the **same module** — there
is no separate driver script (keeps the experiment single-sourced).

---

## (a) Account prerequisites

| Item | Required value / action |
|---|---|
| AWS CLI | Installed and authenticated (`aws sts get-caller-identity` succeeds) |
| Region | `us-east-1` (or your preferred region; use consistently) |
| AMI | Ubuntu 22.04 LTS x86_64, e.g. `ami-0c7217cdde317cfec` (us-east-1); pick the current canonical Ubuntu AMI for your region |
| S3 bucket | `branesim-runs-493652700851` — create if absent; same region as instance |
| IAM instance profile | `BraneSimEc2Profile` — must have the inline policy below |
| Security group | Allow outbound HTTPS (443) for apt/pip/S3; no inbound required for a headless job |
| Subnet | Any subnet with outbound internet (NAT or IGW); note `subnet-xxxxxxxx` |
| Key pair | Optional; omit `--key-name` for headless jobs |

```json
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Action": ["s3:GetObject", "s3:PutObject", "s3:ListBucket"],
      "Resource": [
        "arn:aws:s3:::branesim-runs-493652700851",
        "arn:aws:s3:::branesim-runs-493652700851/*"
      ]
    }
  ]
}
```

`ffmpeg` (for mp4 renders) is installed by the user-data bootstrap; `scipy>=1.15`
(for `sph_harm_y`) is pinned in `requirements.txt`. See `orchestration/aws/README.md`
for the launcher mechanics and `DEPLOYMENT.md` for the memory-sizing table.

---

## (b) Exact commands

### Step 1 — Build the project archive

```bash
bash /path/to/BraneSim/orchestration/aws/package_project.sh /tmp/branesim-project.tar.gz
# verify the vortex entrypoint is included:
tar -tzf /tmp/branesim-project.tar.gz | grep vortex_seed_render
# Expected: ./branesim/experiments/vortex_seed_render.py
```

### Step 2 — Launch the larger-scale run (recommended: 64³ × 32, relaxation ON)

The experiment is parameterized by env vars (so the same module scales). It
writes to `$BRANESIM_RESULTS_DIR`, which the launcher syncs to S3.  `n_slices=32`
keeps the periodic-operator conditioning and the Krylov memory modest while
giving a smooth carrier; bump it only if you need finer time resolution.

```bash
python /path/to/BraneSim/orchestration/aws/launch_branesim_job.py \
  --region us-east-1 \
  --ami-id ami-0c7217cdde317cfec \
  --instance-type r7i.4xlarge \
  --subnet-id subnet-xxxxxxxx \
  --security-group-id sg-xxxxxxxx \
  --iam-instance-profile BraneSimEc2Profile \
  --s3-bucket branesim-runs-493652700851 \
  --s3-prefix vortex-runs \
  --job-name vortex-seed-64 \
  --project-archive /tmp/branesim-project.tar.gz \
  --volume-size-gb 40 \
  --remote-command 'export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8; export BRANESIM_VORTEX_GRID=64 BRANESIM_VORTEX_NSLICES=32; python -m branesim.experiments.vortex_seed_render && { RUN=$(ls -dt "${BRANESIM_RESULTS_DIR:-./runs}"/vortex_seed_* | head -1); mv "$RUN" "${RUN}_aws"; }'
```

### Run provenance — encode AWS vs local in the run-dir name

The experiment writes its run folder as `vortex_seed_<ts>` regardless of where it
runs, so origin is marked **by renaming the run dir**, not by changing the
experiment code. The trailing `mv "$RUN" "${RUN}_aws"` in the remote-command above
appends an **`_aws`** suffix to the produced run dir *before* it is synced to S3, so
every AWS result is self-identifying:

```
vortex_seed_<ts>_aws/      ← produced on AWS (this runbook)
vortex_seed_<ts>           ← bare name = local pre-test (or append _local to be explicit)
```

The rename is safe: each run dir is self-contained (`config.json` + `iter_*/`), and
`run_measurements` / `fetch_results.sh` take the dir path as-is. For a **local**
run, leave the bare name (understood as local) or append `_local` yourself.

`r7i.4xlarge` (16 vCPU / 128 GB) covers the JFNK Krylov working set at 64³×32
(~`40 × 0.28 GB` ≈ 11 GB) with ample headroom for rendering. For a **seed-only**
render (no solve) a `c7i.2xlarge` is cheaper — add `BRANESIM_VORTEX_RELAX=0`.

### Scale knobs (env vars)

| env var | meaning | default |
|---|---|---|
| `BRANESIM_VORTEX_GRID` | spatial grid, `N` or `NX,NY,NZ` | `48,48,48` |
| `BRANESIM_VORTEX_NSLICES` | time-loop slices (carrier resolution & period P) | `32` |
| `BRANESIM_VORTEX_RELAX` | `1`/`0` — run the rotating-frame-periodic relaxation | `1` (ON) |
| `BRANESIM_VORTEX_RELAX_ITERS` | JFNK outer iterations | `15` |
| `BRANESIM_VORTEX_RENDER` | `1`/`0` — render movies+snapshot (`0` = fast diagnostics-only) | `1` (ON) |
| `BRANESIM_RESULTS_DIR` | output root (set by user-data) | `./runs` |

### Scale options (memory = JFNK Krylov set; render = matplotlib volume)

| grid × slices | GB/vector | JFNK set (×40) | render cost | instance |
|---|---|---|---|---|
| 48³ × 32 | 0.14 | ~5 GB | ~1 min/movie | local pre-test (done) |
| **64³ × 32** | **0.28** | **~11 GB** | **~3–5 min/movie** | **r7i.4xlarge (rec.)** |
| 96³ × 32 | 0.93 | ~37 GB | ~20–40 min/movie | r7i.4xlarge; `--volume-size-gb 60` |

Memory scales with `inner_maxiter` (default 40, set in `vortex_seed_render.py`
`RELAX_INNER_MAXITER`). Larger `n_slices` raises both the per-vector size and the
periodic-operator condition number (~P²), so prefer raising the spatial grid over
the slice count. At 96³ the 3D voxel render dominates wall-clock.

### Step 3 — Retrieve results

The launcher prints `results_s3: s3://branesim-runs-493652700851/vortex-runs/<job-id>/results`.

```bash
bash /path/to/BraneSim/orchestration/aws/fetch_results.sh \
  s3://branesim-runs-493652700851/vortex-runs/<job-id>/results \
  ./runs/<job-id>
```

Retrieved layout (`vortex_seed_<ts>_aws/` — the `_aws` suffix marks AWS origin):
- `config.json`, `manifest.json`, `README.md`
- `iter_0000/` — `config.json`, `world.npz`, `winding_closure.json`,
  `snapshot.png`, `renders/` (volume + slice mp4s, each + `.md`),
  `diagnostics/` (8 devices: CSV + PNG + `.md`) + `report.md`
- `user-data.log` — full bootstrap + run stdout/stderr

---

## (c) Pre-flight check (optional, local)

Confirm the pipeline before paying for the instance (tiny grid, ~seconds):

```bash
BRANESIM_VORTEX_GRID=16 BRANESIM_VORTEX_NSLICES=8 BRANESIM_VORTEX_OUTDIR=/tmp/vtest \
  python -m branesim.experiments.vortex_seed_render
```
