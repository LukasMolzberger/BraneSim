#!/usr/bin/env python3
"""Poll S3 tripwire markers written by ec2_user_data.sh.tmpl and auto-terminate
a DOA instance before it idles for hours.

Standalone usage:
    python orchestration/aws/watch_job.py \
        --region eu-north-1 \
        --instance-id i-0abc123 \
        --s3-bucket branesim-breather-493652700851 \
        --s3-prefix branesim/jobs \
        --job-id baryon-run-20260605-120000

Called programmatically from launch_branesim_job.py when --watch is passed.

Exit codes:
    0  job-complete marker seen (success)
    1  DOA tripwire fired; instance was terminated
    2  timeout waiting for job-complete (instance still running; no auto-action)
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from datetime import datetime, timezone


# ---------------------------------------------------------------------------
# Tripwire schedule: (marker_key, deadline_seconds_after_instance_running)
# ---------------------------------------------------------------------------
TRIPWIRES = [
    ("bootstrap-started", 10 * 60),   # 10 min — cloud-init must have started
    ("bootstrap-pip-ok",  20 * 60),   # 20 min — pip install -e . must complete
    ("run-started",       30 * 60),   # 30 min — command must have begun
]
JOB_COMPLETE_TIMEOUT = 4 * 60 * 60   # 4 h overall; adjust for long runs
POLL_INTERVAL = 30                    # seconds between S3 head-object checks


def _run(args: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(args, capture_output=True, text=True, check=True)


def _marker_exists(s3_bucket: str, s3_prefix: str, job_id: str, marker: str, region: str) -> bool:
    key = f"{s3_prefix}/{job_id}/markers/{marker}"
    result = subprocess.run(
        ["aws", "s3api", "head-object", "--region", region, "--bucket", s3_bucket, "--key", key],
        capture_output=True,
        text=True,
    )
    return result.returncode == 0


def _instance_state(instance_id: str, region: str) -> str:
    result = _run([
        "aws", "ec2", "describe-instances",
        "--region", region,
        "--instance-ids", instance_id,
        "--query", "Reservations[0].Instances[0].State.Name",
        "--output", "text",
    ])
    return result.stdout.strip()


def _terminate(instance_id: str, region: str) -> None:
    subprocess.run(
        ["aws", "ec2", "terminate-instances", "--region", region, "--instance-ids", instance_id],
        capture_output=True,
        text=True,
    )


def _wait_for_running(instance_id: str, region: str, timeout: int = 300) -> bool:
    """Wait up to timeout seconds for the instance to reach 'running'."""
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        state = _instance_state(instance_id, region)
        if state == "running":
            return True
        if state in ("terminated", "shutting-down", "stopped"):
            return False
        time.sleep(10)
    return False


def _ts() -> str:
    return datetime.now(timezone.utc).strftime("%H:%M:%S UTC")


def watch(
    region: str,
    instance_id: str,
    s3_bucket: str,
    s3_prefix: str,
    job_id: str,
) -> int:
    """
    Main watch loop.  Returns 0 on success, 1 on DOA, 2 on timeout.
    """
    print(f"[watch {_ts()}] Waiting for instance {instance_id} to reach 'running' ...")
    if not _wait_for_running(instance_id, region, timeout=300):
        state = _instance_state(instance_id, region)
        print(f"[watch {_ts()}] Instance never reached 'running' (state={state}). Aborting watch.")
        return 2
    print(f"[watch {_ts()}] Instance is running. Starting tripwire polling.")

    epoch = time.monotonic()

    # -----------------------------------------------------------------------
    # Phase 1: sequential tripwire checks with individual deadlines
    # -----------------------------------------------------------------------
    for marker, deadline_secs in TRIPWIRES:
        deadline = epoch + deadline_secs
        print(f"[watch {_ts()}] Waiting for marker '{marker}' (deadline T+{deadline_secs//60} min) ...")
        found = False
        while time.monotonic() < deadline:
            if _marker_exists(s3_bucket, s3_prefix, job_id, marker, region):
                elapsed = int(time.monotonic() - epoch)
                print(f"[watch {_ts()}] Marker '{marker}' seen at T+{elapsed}s. OK.")
                found = True
                break
            # Also bail early if the instance already terminated (success path)
            state = _instance_state(instance_id, region)
            if state in ("terminated", "shutting-down"):
                print(f"[watch {_ts()}] Instance {state} before all markers — checking job-complete.")
                if _marker_exists(s3_bucket, s3_prefix, job_id, "job-complete", region):
                    print(f"[watch {_ts()}] job-complete marker present. Job succeeded.")
                    return 0
                print(f"[watch {_ts()}] Instance terminated WITHOUT job-complete. Possible DOA or error.")
                return 1
            time.sleep(POLL_INTERVAL)

        if not found:
            elapsed = int(time.monotonic() - epoch)
            print(
                f"\n[TRIPWIRE {_ts()}] DOA — marker '{marker}' missing after {elapsed}s "
                f"(deadline was T+{deadline_secs}s)."
            )
            print(f"[TRIPWIRE {_ts()}] Auto-terminating instance {instance_id} ...")
            _terminate(instance_id, region)
            print(
                f"[TRIPWIRE {_ts()}] Instance termination requested. "
                "Check user-data.log in S3 results for the failure cause."
            )
            return 1

    # -----------------------------------------------------------------------
    # Phase 2: wait for job-complete or instance self-termination
    # -----------------------------------------------------------------------
    print(f"[watch {_ts()}] All bootstrap tripwires cleared. Waiting for job to finish ...")
    job_deadline = epoch + JOB_COMPLETE_TIMEOUT
    while time.monotonic() < job_deadline:
        if _marker_exists(s3_bucket, s3_prefix, job_id, "job-complete", region):
            elapsed = int(time.monotonic() - epoch)
            print(f"[watch {_ts()}] job-complete marker seen at T+{elapsed}s. Job succeeded.")
            return 0
        state = _instance_state(instance_id, region)
        if state in ("terminated", "shutting-down"):
            # Instance shut down — check whether job-complete marker is present
            # (it may have been written a moment before shutdown).
            time.sleep(15)
            if _marker_exists(s3_bucket, s3_prefix, job_id, "job-complete", region):
                print(f"[watch {_ts()}] Instance {state}; job-complete confirmed. Success.")
                return 0
            print(
                f"[watch {_ts()}] Instance {state} but job-complete marker absent. "
                "Job may have failed — check user-data.log."
            )
            return 1
        time.sleep(POLL_INTERVAL)

    elapsed = int(time.monotonic() - epoch)
    print(
        f"[watch {_ts()}] Timeout ({elapsed}s): job-complete not seen. "
        "Instance may still be running a long job. No auto-action taken."
    )
    return 2


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Poll S3 tripwire markers and auto-terminate a DOA EC2 instance."
    )
    parser.add_argument("--region", required=True)
    parser.add_argument("--instance-id", required=True)
    parser.add_argument("--s3-bucket", required=True)
    parser.add_argument("--s3-prefix", default="branesim/jobs")
    parser.add_argument("--job-id", required=True)
    args = parser.parse_args()

    rc = watch(
        region=args.region,
        instance_id=args.instance_id,
        s3_bucket=args.s3_bucket,
        s3_prefix=args.s3_prefix,
        job_id=args.job_id,
    )
    sys.exit(rc)


if __name__ == "__main__":
    main()
