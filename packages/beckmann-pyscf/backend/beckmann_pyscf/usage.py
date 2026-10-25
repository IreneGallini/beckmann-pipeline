"""
Free-tier accounting for the Cloud Run deployment: how much of this
month's Cloud Run free tier the service has used, read from Cloud
Monitoring's own billing metric. Used by two things (see DEPLOY.md,
"Staying at $0"):

- app.py's usage gate: before every scan, refuse if the scan could push
  this month's usage past the allowance.
- .github/scripts/free_tier_guard.py: the hourly backstop that makes the
  service private when usage nears the allowance (and public again once a
  new month resets it).

Standard library only (urllib), on purpose: the guard runs in a bare
GitHub Actions Python without torch/pyscf installed, and imports this file
directly by path.

The metric is `run.googleapis.com/container/billable_instance_time`: the
seconds Cloud Run actually bills (request time plus container startup and
shutdown, under request-based billing). Reading Google Cloud metrics
through the Monitoring API is free at this volume.
"""
import datetime as _dt
import json
import os
import urllib.parse
import urllib.request

# Cloud Run's monthly free tier (request-based billing, per billing
# account): https://cloud.google.com/run/pricing
FREE_VCPU_SECONDS = 180_000
FREE_GIB_SECONDS = 360_000
FREE_EGRESS_BYTES = 1 * 1024**3  # outbound data transfer, North America

_METADATA = "http://metadata.google.internal/computeMetadata/v1"
_MONITORING = "https://monitoring.googleapis.com/v3"


def monthly_limit_seconds(cpu: float, memory_gib: float, fraction: float) -> float:
    """Billable instance-seconds allowed per month: `fraction` of whichever
    of the CPU or memory free tiers runs out first at this instance size.
    E.g. 1 vCPU / 4 GiB / 0.75 -> 0.75 * min(180000, 90000) = 67500 s."""
    return fraction * min(FREE_VCPU_SECONDS / cpu, FREE_GIB_SECONDS / memory_gib)


def month_start(now: _dt.datetime | None = None) -> _dt.datetime:
    """Start of the current month, UTC. (Google's free-tier month follows
    the billing account's calendar; the safety fraction covers the few
    hours' difference.)"""
    now = now or _dt.datetime.now(_dt.timezone.utc)
    return now.replace(day=1, hour=0, minute=0, second=0, microsecond=0)


def metadata_token() -> tuple[str, str]:
    """(access_token, project_id) for the service's own runtime service
    account, from the Cloud Run metadata server."""
    def _get(path):
        req = urllib.request.Request(_METADATA + path, headers={"Metadata-Flavor": "Google"})
        with urllib.request.urlopen(req, timeout=5) as resp:
            return resp.read().decode()
    token = json.loads(_get("/instance/service-accounts/default/token"))["access_token"]
    return token, _get("/project/project-id")


def _sum_metric(token: str, project: str, metric: str, service: str,
                start: _dt.datetime, end: _dt.datetime) -> float:
    """Sum of a DELTA metric for one Cloud Run service over [start, end]:
    hourly ALIGN_SUM per series, REDUCE_SUM across revisions/instances,
    then the hourly buckets added up here."""
    params = {
        "filter": (f'metric.type = "{metric}" AND resource.type = "cloud_run_revision" '
                   f'AND resource.labels.service_name = "{service}"'),
        "interval.startTime": start.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "interval.endTime": end.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "aggregation.alignmentPeriod": "3600s",
        "aggregation.perSeriesAligner": "ALIGN_SUM",
        "aggregation.crossSeriesReducer": "REDUCE_SUM",
        "pageSize": "1000",
    }
    total = 0.0
    page_token = None
    while True:
        if page_token:
            params["pageToken"] = page_token
        url = f"{_MONITORING}/projects/{project}/timeSeries?" + urllib.parse.urlencode(params)
        req = urllib.request.Request(url, headers={"Authorization": f"Bearer {token}"})
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read())
        for series in data.get("timeSeries", []):
            for point in series.get("points", []):
                value = point["value"]
                total += float(value.get("doubleValue", value.get("int64Value", 0)))
        page_token = data.get("nextPageToken")
        if not page_token:
            return total


def month_usage(token: str, project: str, service: str,
                now: _dt.datetime | None = None) -> dict:
    """{"billable_seconds", "egress_bytes"} for `service` so far this month."""
    end = now or _dt.datetime.now(_dt.timezone.utc)
    start = month_start(end)
    return {
        "billable_seconds": _sum_metric(
            token, project, "run.googleapis.com/container/billable_instance_time", service, start, end),
        "egress_bytes": _sum_metric(
            token, project, "run.googleapis.com/container/network/sent_bytes_count", service, start, end),
    }


def check(usage: dict, cpu: float, memory_gib: float, fraction: float,
          reserve_seconds: float = 0.0) -> str | None:
    """None if `usage` plus `reserve_seconds` more billable time stays
    within `fraction` of the free tier; otherwise a short reason."""
    limit = monthly_limit_seconds(cpu, memory_gib, fraction)
    if usage["billable_seconds"] + reserve_seconds > limit:
        return (f"billable time {usage['billable_seconds']:.0f}s + {reserve_seconds:.0f}s reserved "
                f"> {limit:.0f}s allowance")
    if usage["egress_bytes"] > fraction * FREE_EGRESS_BYTES:
        return f"egress {usage['egress_bytes'] / 1e6:.0f} MB > {fraction * FREE_EGRESS_BYTES / 1e6:.0f} MB allowance"
    return None
