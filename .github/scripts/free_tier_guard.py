"""
Free-tier guard for the Cloud Run service (DEPLOY.md, "Staying at $0").

Rule: the service is public exactly when this month's usage, plus one more
hour of nonstop use (the time until the next hourly check), fits in the
free-tier allowance. Otherwise its public access (allUsers ->
roles/run.invoker) is removed. Cloud Run then rejects every request at
Google's front end, before it reaches the container, so nothing more is
billed. Public access comes back on its own once a new month resets usage.

This is the backstop. The app's own gate (BECKMANN_USAGE_GATE) already
refuses scans that could exceed the allowance; this also covers non-scan
traffic, such as a bot hammering the API.

Run by .github/workflows/free-tier-guard.yml (hourly) and at the end of
every deploy (a deploy re-adds public access, so the guard re-checks right
after). Needs gcloud authenticated as the deploy service account
(roles/run.admin + roles/monitoring.viewer).

Env: PROJECT_ID, REGION, SERVICE, CLOUD_RUN_CPU, CLOUD_RUN_MEMORY_GIB,
FREE_TIER_FRACTION.
"""
import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

# beckmann_pyscf/usage.py is standard-library only; load it straight from
# its file so this script needs no package install (no torch/pyscf).
_USAGE = Path(__file__).resolve().parents[2] / "packages/beckmann-pyscf/backend/beckmann_pyscf/usage.py"
_spec = importlib.util.spec_from_file_location("usage", _USAGE)
usage = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(usage)

CHECK_INTERVAL_SECONDS = 3600  # worst case: nonstop use until the next run


def gcloud(*args: str, check: bool = True) -> subprocess.CompletedProcess:
    return subprocess.run(["gcloud", *args], capture_output=True, text=True, check=check)


def is_public(service: str, region: str) -> bool:
    policy = json.loads(gcloud("run", "services", "get-iam-policy", service,
                               "--region", region, "--format", "json").stdout or "{}")
    return any(b.get("role") == "roles/run.invoker" and "allUsers" in b.get("members", [])
               for b in policy.get("bindings", []))


def summary(text: str) -> None:
    print(text)
    path = os.environ.get("GITHUB_STEP_SUMMARY")
    if path:
        with open(path, "a") as f:
            f.write(text + "\n")


def main() -> int:
    project = os.environ["PROJECT_ID"]
    region = os.environ["REGION"]
    service = os.environ.get("SERVICE", "beckmann-pyscf")
    cpu = float(os.environ.get("CLOUD_RUN_CPU", "1"))
    mem = float(os.environ.get("CLOUD_RUN_MEMORY_GIB", "4"))
    fraction = float(os.environ.get("FREE_TIER_FRACTION", "0.75"))

    exists = gcloud("run", "services", "describe", service, "--region", region,
                    "--format", "value(metadata.name)", check=False)
    if exists.returncode != 0:
        summary(f"Service {service} doesn't exist yet in {region}; nothing to guard.")
        return 0

    token = gcloud("auth", "print-access-token").stdout.strip()
    used = usage.month_usage(token, project, service)
    limit = usage.monthly_limit_seconds(cpu, mem, fraction)
    reason = usage.check(used, cpu, mem, fraction, reserve_seconds=CHECK_INTERVAL_SECONDS)
    public = is_public(service, region)

    summary(f"### Free-tier guard\n"
            f"- billable time this month: **{used['billable_seconds']:.0f} s** "
            f"of a {limit:.0f} s allowance ({fraction:.0%} of the free tier at {cpu:g} vCPU / {mem:g} GiB)\n"
            f"- outbound data this month: **{used['egress_bytes'] / 1e6:.1f} MB** "
            f"of {fraction * usage.FREE_EGRESS_BYTES / 1e6:.0f} MB\n"
            f"- service currently public: {public}")

    if reason and public:
        gcloud("run", "services", "remove-iam-policy-binding", service, "--region", region,
               "--member", "allUsers", "--role", "roles/run.invoker")
        summary(f"**Closed** public access: {reason}. It reopens automatically next month.")
    elif not reason and not public:
        gcloud("run", "services", "add-iam-policy-binding", service, "--region", region,
               "--member", "allUsers", "--role", "roles/run.invoker")
        summary("**Reopened** public access: usage is back within the allowance.")
    else:
        summary("No change.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
