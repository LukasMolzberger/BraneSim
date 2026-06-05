#!/usr/bin/env python3
"""Launch one AWS EC2 job for the BraneSim baryon pipeline.

The instance downloads a project archive from S3, runs a command, uploads results,
and self-terminates to avoid idle cost.
"""

from __future__ import annotations

import argparse
import base64
import json
import os
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path


def run_cmd(args: list[str], *, capture_output: bool = False) -> subprocess.CompletedProcess[str]:
    try:
        return subprocess.run(
            args,
            check=True,
            text=True,
            capture_output=capture_output,
        )
    except subprocess.CalledProcessError as e:
        # Surface the captured AWS error instead of swallowing it (was a silent
        # exit-254 with no message when capture_output=True).
        import sys
        if e.stderr:
            print(f"[run_cmd] command failed (exit {e.returncode}):\n{e.stderr}", file=sys.stderr)
        if e.stdout:
            print(f"[run_cmd] stdout:\n{e.stdout}", file=sys.stderr)
        raise


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Launch a BraneSim job on AWS EC2.")
    parser.add_argument("--region", default="us-east-1")
    parser.add_argument("--ami-id", required=True, help="Ubuntu AMI id, e.g. ami-xxxxxxxx")
    # Memory-optimized default: the 4D block BVP is memory-bound (state ~
    # n_slices*N^dim*m*8 bytes; JFNK holds ~30 such vectors). See DEPLOYMENT.md
    # for the sizing table; bump to r7i.8xlarge/16xlarge or x2iedn for big runs.
    parser.add_argument("--instance-type", default="r7i.4xlarge")
    parser.add_argument("--key-name", default=None)
    parser.add_argument("--subnet-id", default=None)
    parser.add_argument("--security-group-id", action="append", default=[], help="Repeat for multiple SGs")
    parser.add_argument("--iam-instance-profile", default=None, help="Instance profile name with S3 access")

    parser.add_argument("--s3-bucket", required=True)
    parser.add_argument("--s3-prefix", default="branesim/jobs")
    parser.add_argument("--job-name", default="baryon-run")

    parser.add_argument("--project-archive", default=None, help="Local .tar.gz project archive to upload")
    parser.add_argument("--project-s3-uri", default=None, help="Existing project archive in S3")

    parser.add_argument(
        "--remote-command",
        default=(
            "python -m branesim.run_experiment "
            "--config orchestration/configs/branesim_ivp_smoke.json "
            "--output-dir \"$BRANESIM_RESULTS_DIR\""
        ),
        help="Shell command executed on the EC2 instance inside the project directory. "
             "Default is the branesim IVP smoke; pass a bvp_dirichlet config for a block solve.",
    )

    parser.add_argument("--enable-gdrive-sync", action="store_true")
    parser.add_argument("--gdrive-remote", default="")
    parser.add_argument("--gdrive-dest", default="BraneSim")

    parser.add_argument(
        "--watch",
        action="store_true",
        help="After launch, run watch_job.py in-process to poll S3 markers and "
             "auto-terminate if a DOA tripwire fires.",
    )

    parser.add_argument("--spot", action="store_true", help="Launch as one-time spot instance")
    parser.add_argument("--volume-size-gb", type=int, default=120)

    return parser.parse_args()


def ensure_project_uri(args: argparse.Namespace, job_id: str) -> str:
    if args.project_s3_uri:
        return args.project_s3_uri

    if not args.project_archive:
        raise ValueError("Provide either --project-archive or --project-s3-uri")

    archive = Path(args.project_archive)
    if not archive.exists():
        raise FileNotFoundError(f"Project archive not found: {archive}")

    target_uri = f"s3://{args.s3_bucket}/{args.s3_prefix}/{job_id}/input/project.tar.gz"
    run_cmd(["aws", "s3", "cp", str(archive), target_uri])
    return target_uri


def render_user_data(template_path: Path, replacements: dict[str, str]) -> str:
    text = template_path.read_text(encoding="utf-8")
    for key, value in replacements.items():
        text = text.replace(key, value)
    return text


def main() -> None:
    args = parse_args()

    now = datetime.now(timezone.utc)
    job_id = f"{args.job_name}-{now.strftime('%Y%m%d-%H%M%S')}"
    s3_results_uri = f"s3://{args.s3_bucket}/{args.s3_prefix}/{job_id}/results"

    project_uri = ensure_project_uri(args, job_id)

    template_path = Path(__file__).with_name("ec2_user_data.sh.tmpl")
    command_b64 = base64.b64encode(args.remote_command.encode("utf-8")).decode("ascii")

    user_data = render_user_data(
        template_path,
        {
            "__JOB_ID__": job_id,
            "__S3_BUCKET__": args.s3_bucket,
            "__S3_PREFIX__": args.s3_prefix,
            "__S3_PROJECT_URI__": project_uri,
            "__S3_RESULTS_URI__": s3_results_uri,
            "__REMOTE_COMMAND_B64__": command_b64,
            "__ENABLE_GDRIVE_SYNC__": "true" if args.enable_gdrive_sync else "false",
            "__GDRIVE_REMOTE__": args.gdrive_remote,
            "__GDRIVE_DEST__": args.gdrive_dest,
        },
    )

    with tempfile.NamedTemporaryFile("w", suffix=".sh", delete=False) as tmp:
        tmp.write(user_data)
        user_data_path = tmp.name

    run_args = [
        "aws",
        "ec2",
        "run-instances",
        "--region",
        args.region,
        "--image-id",
        args.ami_id,
        "--instance-type",
        args.instance_type,
        "--count",
        "1",
        "--user-data",
        f"file://{user_data_path}",
        "--instance-initiated-shutdown-behavior",
        "terminate",
        "--block-device-mappings",
        json.dumps(
            [
                {
                    "DeviceName": "/dev/sda1",
                    "Ebs": {"VolumeSize": args.volume_size_gb, "VolumeType": "gp3", "DeleteOnTermination": True},
                }
            ]
        ),
        "--tag-specifications",
        json.dumps(
            [
                {
                    "ResourceType": "instance",
                    "Tags": [
                        {"Key": "Name", "Value": f"BraneSim-{job_id}"},
                        {"Key": "Project", "Value": "BraneSim"},
                        {"Key": "JobId", "Value": job_id},
                    ],
                }
            ]
        ),
    ]

    if args.key_name:
        run_args += ["--key-name", args.key_name]
    if args.subnet_id:
        run_args += ["--subnet-id", args.subnet_id]
    if args.security_group_id:
        run_args += ["--security-group-ids", *args.security_group_id]
    if args.iam_instance_profile:
        run_args += ["--iam-instance-profile", f"Name={args.iam_instance_profile}"]
    if args.spot:
        run_args += [
            "--instance-market-options",
            json.dumps({"MarketType": "spot", "SpotOptions": {"SpotInstanceType": "one-time"}}),
        ]

    output = run_cmd(run_args, capture_output=True)
    result = json.loads(output.stdout)
    instance_id = result["Instances"][0]["InstanceId"]

    os.unlink(user_data_path)

    print("EC2 job launched")
    print(f"  job_id:          {job_id}")
    print(f"  instance_id:     {instance_id}")
    print(f"  project_archive: {project_uri}")
    print(f"  results_s3:      {s3_results_uri}")
    print(f"  markers_s3:      s3://{args.s3_bucket}/{args.s3_prefix}/{job_id}/markers/")
    print("Instance will self-terminate after syncing results.")

    if args.watch:
        import importlib.util, sys as _sys
        watcher_path = Path(__file__).with_name("watch_job.py")
        spec = importlib.util.spec_from_file_location("watch_job", watcher_path)
        watcher_mod = importlib.util.module_from_spec(spec)  # type: ignore[arg-type]
        spec.loader.exec_module(watcher_mod)  # type: ignore[union-attr]
        watcher_mod.watch(
            region=args.region,
            instance_id=instance_id,
            s3_bucket=args.s3_bucket,
            s3_prefix=args.s3_prefix,
            job_id=job_id,
        )


if __name__ == "__main__":
    main()
