# AWS Runbook for BraneSim

> **See `DEPLOYMENT.md` (repo root) for the memory-sizing table and instance
> recommendations.** This file documents the launcher mechanics only.

This folder provides scripts to run memory-heavy `branesim` block-solve jobs on
AWS EC2 and automatically clean up compute resources. The 4D block solve is
memory-bound → use memory-optimized instances (`r7i.*` / `x2iedn.*`), not the
old compute-optimized default.

## What the launcher does

`launch_branesim_job.py` will:

1. Upload a project archive to S3 (or use an existing S3 archive).
2. Start an EC2 instance with user-data bootstrap.
3. Run your pipeline command on the instance.
4. Sync results/logs to S3.
5. Optionally sync results to Google Drive via `rclone`.
6. Shutdown the instance (instance-initiated shutdown behavior is `terminate`).

This avoids lingering compute costs after completion.

## Prerequisites

- Local machine: AWS CLI authenticated
- AWS: subnet/security group/AMI/instance profile prepared
- Instance profile permissions: `s3:GetObject`, `s3:PutObject`, `s3:ListBucket`
- Optional Google Drive sync: provide an `rclone` remote name and destination path

## 1) Build project archive

```bash
orchestration/aws/package_project.sh /tmp/branesim-project.tar.gz
```

## 2) Launch a job

```bash
python orchestration/aws/launch_branesim_job.py \
  --region us-east-1 \
  --ami-id ami-xxxxxxxx \
  --instance-type r7i.4xlarge \
  --subnet-id subnet-xxxxxxxx \
  --security-group-id sg-xxxxxxxx \
  --iam-instance-profile BraneSimEc2Profile \
  --s3-bucket my-branesim-bucket \
  --s3-prefix baryon-runs \
  --project-archive /tmp/branesim-project.tar.gz
```

Default remote command now runs the branesim entrypoint:

```bash
python -m branesim.run_experiment --config orchestration/configs/branesim_ivp_smoke.json --output-dir "$BRANESIM_RESULTS_DIR"
```

Point `--config` at `orchestration/configs/branesim_bvp_dirichlet.json` (or your
own) for a block solve.

Use `--remote-command` to run a custom workflow if needed.

## 3) Retrieve results

Use the `results_s3` URI printed by the launcher:

```bash
orchestration/aws/fetch_results.sh s3://my-branesim-bucket/baryon-runs/<job-id>/results ./downloaded-results
```

## Optional Google Drive sync

Add to launch command:

```bash
--enable-gdrive-sync --gdrive-remote mydrive --gdrive-dest BraneSim
```

`mydrive` must be a valid `rclone` remote configured on the instance image/bootstrap process.

## Cost safety

- EC2 instances are configured with `instance-initiated-shutdown-behavior=terminate`.
- User-data script always runs a result sync + shutdown path (even on failure).
- Root EBS volume is `DeleteOnTermination=true`.
