# Deploying beckmann-pyscf (Cloud Run API + GitHub Pages page), at $0

The site is split into three parts, each placed where it's free:

| Part | Where | Deployed by |
|---|---|---|
| The web page (HTML + the 24 MB Ketcher editor) | **GitHub Pages** | `.github/workflows/deploy-pages.yml` |
| The prediction API (`POST /predict`) | **Google Cloud Run** | `.github/workflows/deploy-cloud-run.yml` |
| The container image | **Docker Hub** (public repo) | the Cloud Run workflow |

Two more pieces keep Cloud Run inside its free tier:
- a **per-scan usage gate** in the app;
- an **hourly free-tier guard**, `.github/workflows/free-tier-guard.yml`.

GitHub authenticates to Google with **Workload Identity Federation**: short-lived tokens, and no Google key file anywhere.

This guide covers the one-time setup (steps 1–10), the first deploy (step 11), and how to check it (step 12).

---

## Staying at $0: what's guaranteed and what isn't

**Be clear on this first: Google Cloud has no hard spending cap.** A budget only sends email, and billing data lags by hours. So the setup below doesn't rely on billing alerts. It keeps *usage* below the free tier by construction, and shuts the service off automatically before it can leave the free tier.

| Possible charge | How it's kept at $0 |
|---|---|
| **Cloud Run compute.** Free each month: 180,000 vCPU-s and 360,000 GiB-s. At 1 vCPU / 4 GiB, memory runs out first, after 90,000 billable seconds (25 h) | The app allows only **75%** of that, 67,500 s ≈ **35–45 scans/month**. **Before every scan** it reads this month's billable time from Cloud Monitoring, and refuses the scan if the scan *at its worst case* (55 min timeout + 2 min startup) could go over. If usage can't be read, it refuses too, failing closed. |
| **Cloud Run, non-scan traffic** (e.g. a bot hammering the API) | The **hourly guard** removes public access once usage plus one more hour of nonstop use would pass the 75% line. Google's front end then rejects every request before it reaches the container, so nothing more is billed. Access reopens on its own when the new month resets usage. Deploys re-run the guard immediately. |
| **Cloud Run outbound data.** Free: 1 GiB/month in North America | The page and Ketcher are served by GitHub Pages, not Cloud Run. Cloud Run's `/` redirects to Pages and it serves no static files. A prediction result is a few KB. The guard also closes the service at 75% of the outbound allowance. |
| **Cloud Run requests.** Free: 2 million/month | One request per prediction. The browser sends it as a CORS "simple request", so there's no extra preflight request. |
| **Image storage** | Docker Hub hosts the image (free for public images), so there's no Artifact Registry repository at all. |
| Cloud Monitoring reads, Cloud Logging, IAM, Workload Identity Federation, GitHub Actions and Pages (public repo) | Free at this volume. |

**What remains, honestly:**
- **The 25% margin** absorbs a gap of a few hours between the UTC month the app uses and Google's billing month, a delayed hourly check, and cold starts.
- **The free tier is per billing account.** If you run other Cloud Run services on the *same billing account*, they share it, and this setup only counts its own service. Use a billing account with nothing else on Cloud Run.
- **Prices and free-tier sizes are Google's to change.** Check https://cloud.google.com/run/pricing once in a while. The numbers live in `packages/beckmann-pyscf/backend/beckmann_pyscf/usage.py`.
- **A $1 budget alert (step 3) is the tripwire.** It emails you if a single cent is ever charged, so a surprise can't grow quietly.
- **GitHub pauses scheduled workflows after 60 days without repo activity.** It emails you first. The per-scan gate keeps working, but the hourly backstop stops until you re-enable it from the Actions tab.

**The trade-off for $0:** capacity is roughly 35–45 scans a month. After that, the page says "This month's free compute allowance … is used up" until the 1st (UTC).

---

## Cloud Run settings, and why

All of these are in the workflow's `gcloud run deploy` command, so the service config is version-controlled.

| Setting | Value | Why |
|---|---|---|
| `--cpu` | `1` | PySCF, most of a scan, is pinned to 1 thread anyway. A second vCPU would mostly double the bill. |
| `--memory` | `4Gi` | Measured peak for α-tetralone: **2.6 GB** for the whole process tree. The PySCF step is capped at 2000 MB by `PYSCF_MAX_MEMORY`. |
| `--concurrency` / `--max-instances` | `1` / `1` | One request per instance, and at most one instance. The usage gate counts all instances, so raising max instances (`CLOUD_RUN_MAX_INSTANCES`) doesn't raise the monthly cap. It only lets two scans run at once. |
| `--min-instances` | `0` | Nothing billed when idle. |
| `--timeout` | `3600` | Cloud Run's maximum. The app's own scan timeout (`BECKMANN_JOB_TIMEOUT=3300`) fires first, with a clean "timeout" message. |
| `--cpu-throttling` | (request-based billing) | Billed only while a request runs, and the scan *is* the request. |
| `--execution-environment` | `gen2` | Full Linux process-group handling, used to kill a cancelled scan. |
| `BECKMANN_USAGE_GATE`, `BECKMANN_CLOUD_RUN_CPU`/`_MEMORY_GIB`, `BECKMANN_FREE_TIER_FRACTION` | `1`, `1`/`4`, `0.75` | The per-scan free-tier check (above). The CPU/memory values must match `--cpu`/`--memory`. |
| `BECKMANN_FRONTEND_URL`, `BECKMANN_CORS_ORIGINS` | the Pages URL / its origin | `/` redirects to Pages, no static files are served, and the Pages page may call the API. |
| `BECKMANN_API_MODE` | `sync` | The scan runs inside `/predict` (below). |
| `BECKMANN_TRUSTED_PROXIES` | `1` | The rate limit keys on the real client IP (below). |

### Why the scan runs inside the request
With request-based billing, Cloud Run throttles an instance's CPU to near zero whenever it isn't handling a request, so background work would stall.

So `POST /predict` runs the whole scan within the request. While it runs, it streams a space every 20 s, then sends the JSON result. The heartbeats keep the long connection alive. They also catch a closed tab: within about two heartbeats (~40 s) a write fails, and the scan's processes are killed. That matters for cost too: an abandoned scan doesn't keep billing.

### Why `BECKMANN_TRUSTED_PROXIES=1`
The rate limiter needs the real client IP. Google's front end **appends** the address of the client that connected to it to `X-Forwarded-For`, after anything the client already sent. So the header is `<anything the client sent>, <real client IP>`. Nothing inside the container adds another entry.

Trusting exactly one proxy (`ProxyFix(x_for=1)`) makes the app use the **right-most** entry, the one Google wrote. The other values are wrong:
- `0` would key everyone on Cloud Run's internal address, giving one shared limit.
- `2` would trust an entry the client can forge.

Check it on your deployment in step 12b.

---

## One-time setup

### Before you start
- **Use a personal Google account** (e.g. @gmail.com), not a Macalester/CMU account. University Google Workspace organizations often block public (`allUsers`) services, and some block billing-linked projects.
- **Use a billing account with no other Cloud Run usage** (see "What remains, honestly" above). A new account also usually comes with trial credit.
- **The repo must be public.** GitHub Pages and free Actions minutes depend on it.

Do steps 1–6 in **Google Cloud Shell**: open https://console.cloud.google.com and click `>_` at the top right. Use **one session**, since later steps use earlier variables. If it resets, re-run step 1, plus the `BILLING_ACCOUNT=` and `PROJECT_NUMBER=` lines from step 2.

### 1. Set variables

```bash
# Project IDs are globally unique: 6-30 chars, lowercase letters, digits, hyphens.
PROJECT_ID="beckmann-predictor-$RANDOM"
REGION="us-central1"                       # Tier-1 pricing; close to MN/PA
GITHUB_REPO="IreneGallini/beckmann-pipeline"
echo "PROJECT_ID=$PROJECT_ID"
```

### 2. Create the project and link billing

```bash
gcloud projects create "$PROJECT_ID" --name="Beckmann predictor"
gcloud config set project "$PROJECT_ID"

gcloud billing accounts list               # copy the ACCOUNT_ID of your billing account
BILLING_ACCOUNT="XXXXXX-XXXXXX-XXXXXX"     # <- paste it here
gcloud billing projects link "$PROJECT_ID" --billing-account="$BILLING_ACCOUNT"

PROJECT_NUMBER=$(gcloud projects describe "$PROJECT_ID" --format='value(projectNumber)')
echo "PROJECT_NUMBER=$PROJECT_NUMBER"
```

If `gcloud billing accounts list` is empty, create a billing account first at https://console.cloud.google.com/billing.

### 3. Enable the APIs and set a $1 "first cent" alert

```bash
gcloud services enable \
  run.googleapis.com \
  monitoring.googleapis.com \
  iam.googleapis.com \
  iamcredentials.googleapis.com \
  sts.googleapis.com \
  cloudresourcemanager.googleapis.com \
  billingbudgets.googleapis.com

gcloud billing budgets create \
  --billing-account="$BILLING_ACCOUNT" \
  --display-name="beckmann-predictor: any charge" \
  --budget-amount=1USD \
  --filter-projects="projects/$PROJECT_NUMBER" \
  --threshold-rule=percent=0.01 \
  --threshold-rule=percent=0.5 \
  --threshold-rule=percent=1.0
```

You get an email at $0.01, $0.50 and $1. Free-tier usage costs $0, so if this ever fires, something is being charged. Check the Billing → Reports page, and if in doubt stop the service (see "Turning it off").

### 4. Create the two service accounts

- **`beckmann-runtime`**: the identity the running app uses. It gets only **`roles/monitoring.viewer`**, enough for the usage gate to read its own billable time. It is not the default compute account, which has broad Editor rights.
- **`beckmann-deployer`**: the identity GitHub Actions uses. It deploys the service, runs it as `beckmann-runtime`, and reads usage for the hourly guard.

```bash
gcloud iam service-accounts create beckmann-runtime \
  --display-name="beckmann-pyscf Cloud Run runtime"
gcloud iam service-accounts create beckmann-deployer \
  --display-name="GitHub Actions deployer for beckmann-pyscf"

RUNTIME_SA="beckmann-runtime@${PROJECT_ID}.iam.gserviceaccount.com"
DEPLOY_SA="beckmann-deployer@${PROJECT_ID}.iam.gserviceaccount.com"

# Runtime: read Cloud Monitoring metrics (the per-scan usage gate).
gcloud projects add-iam-policy-binding "$PROJECT_ID" \
  --member="serviceAccount:$RUNTIME_SA" --role="roles/monitoring.viewer" --condition=None

# Deployer: deploy/update the service and change its public access. run.admin
# (not run.developer) because --allow-unauthenticated and the guard change
# the service's IAM policy.
gcloud projects add-iam-policy-binding "$PROJECT_ID" \
  --member="serviceAccount:$DEPLOY_SA" --role="roles/run.admin" --condition=None

# Deployer: read usage metrics (the hourly guard).
gcloud projects add-iam-policy-binding "$PROJECT_ID" \
  --member="serviceAccount:$DEPLOY_SA" --role="roles/monitoring.viewer" --condition=None

# Deployer: allowed to deploy services that run as beckmann-runtime, and only that SA.
gcloud iam service-accounts add-iam-policy-binding "$RUNTIME_SA" \
  --member="serviceAccount:$DEPLOY_SA" --role="roles/iam.serviceAccountUser"
```

### 5. Set up Workload Identity Federation for GitHub Actions

This lets workflow runs from **this repo's `main` branch, and nothing else**, act as `beckmann-deployer` without a key file.

```bash
gcloud iam workload-identity-pools create github \
  --location=global \
  --display-name="GitHub Actions"

gcloud iam workload-identity-pools providers create-oidc github-actions \
  --location=global \
  --workload-identity-pool=github \
  --display-name="GitHub Actions OIDC" \
  --issuer-uri="https://token.actions.githubusercontent.com" \
  --attribute-mapping="google.subject=assertion.sub,attribute.repository=assertion.repository,attribute.ref=assertion.ref" \
  --attribute-condition="assertion.repository=='${GITHUB_REPO}' && assertion.ref=='refs/heads/main'"

gcloud iam service-accounts add-iam-policy-binding "$DEPLOY_SA" \
  --role="roles/iam.workloadIdentityUser" \
  --member="principalSet://iam.googleapis.com/projects/${PROJECT_NUMBER}/locations/global/workloadIdentityPools/github/attribute.repository/${GITHUB_REPO}"
```

The attribute condition is the security boundary. Tokens from forks, other repos, or other branches are rejected. Scheduled runs (the guard) and `workflow_run` runs (Pages) execute on `main`, so they pass. Running a workflow by hand from another branch fails at the auth step, on purpose. IAM changes can take a few minutes to propagate; if the first run fails with a permission error right after setup, re-run it.

### 6. Print the values for GitHub

```bash
cat <<EOF

GitHub repository VARIABLES (not secrets):

GCP_PROJECT_ID    = $PROJECT_ID
GCP_REGION        = $REGION
GCP_WIF_PROVIDER  = projects/$PROJECT_NUMBER/locations/global/workloadIdentityPools/github/providers/github-actions
GCP_DEPLOY_SA     = $DEPLOY_SA
GCP_RUNTIME_SA    = $RUNTIME_SA
EOF
```

### 7. Create the Docker Hub repository and token

1. Sign up or log in at https://hub.docker.com (the free Personal plan is enough).
2. **Create repository**: name `beckmann-pyscf`, visibility **Public**. It must be public: private images aren't free, and Cloud Run couldn't pull them without extra setup.
3. Avatar → **Account settings → Personal access tokens → Generate new token**. Name it `github-actions-beckmann` and give it **Read & Write** access. Copy it; it's shown once.

### 8. Add the variables and the secret on GitHub

On github.com/IreneGallini/beckmann-pipeline, go to **Settings → Secrets and variables → Actions**:
- **Variables** tab → **New repository variable**, once for each:
  - the five values from step 6
  - `DOCKERHUB_USERNAME` = your Docker Hub username
- **Secrets** tab → **New repository secret**: `DOCKERHUB_TOKEN` = the token from step 7. This is the only secret. A leaked Docker Hub token could overwrite the public image, so keep it in Secrets, never in a variable or a file.
- Optional variables: `PAGES_URL`, only if your Pages address isn't `https://irenegallini.github.io/beckmann-pipeline/`; and `CLOUD_RUN_MAX_INSTANCES` (default `1`).
- If `HF_TOKEN` or `HF_SPACE` are still there from the earlier Hugging Face plan, delete them. Revoke that token on huggingface.co, and delete the Space if you made one.

Until `GCP_PROJECT_ID` and `DOCKERHUB_USERNAME` are set, the workflows skip themselves instead of failing.

### 9. Turn on GitHub Pages

**Settings → Pages → Build and deployment → Source: "GitHub Actions".** No branch selection is needed.

### 10. Make sure scheduled workflows are allowed

**Settings → Actions → General → Actions permissions**: "Allow all actions and reusable workflows" (the default). The guard runs on a schedule, and GitHub only runs schedules from the default branch, `main`, so the guard starts once the merge in step 11 lands.

---

## 11. First deploy

1. Merge the `deploy-hf-redis` branch into `main` through a PR.
2. **Deploy to Cloud Run** runs first. Its first build takes ~10–20 min, mostly installing torch, PySCF and the rest; later builds reuse cached layers. It's also the Dockerfile's first real build, so if it fails, the log says where. At the end it runs the free-tier guard, whose summary shows this month's usage.
3. **Deploy page to GitHub Pages** then runs automatically: it looks up the API URL and publishes the page. Your site is `https://irenegallini.github.io/beckmann-pipeline/`.
4. **Free-tier guard** runs every hour from then on. Each run's summary shows usage against the allowance.

You can re-run any of them from the **Actions** tab (select the workflow → Run workflow, on `main`).

## 12. Verify

Get the API URL from Cloud Shell or the deploy run's summary:
```bash
URL=$(gcloud run services describe beckmann-pyscf --region us-central1 --format='value(status.url)')
```

**a) Health and redirect:**
```bash
curl -s "$URL/health"                                  # -> {"status":"ok"}
curl -s -o /dev/null -w "%{http_code} %{redirect_url}\n" "$URL/"   # -> 302 https://irenegallini.github.io/beckmann-pipeline/
curl -s -o /dev/null -w "%{http_code}\n" "$URL/static/ketcher/standalone/index.html"   # -> 404 (not served from Cloud Run)
```
The first request after idle cold-starts the instance, which can take ~30–60 s for this large image.

**b) Rate-limit IP (checks `BECKMANN_TRUSTED_PROXIES`).** Temporarily turn on the debug endpoint:
```bash
gcloud run services update beckmann-pyscf --region us-central1 \
  --update-env-vars BECKMANN_DEBUG_CLIENT_IP=1
```
Then **from your own computer's terminal**, so you know which IP to expect:
```bash
URL=https://beckmann-pyscf-xxxxxxxx-uc.a.run.app      # your API URL
curl -s https://api.ipify.org; echo                    # your public IP
curl -s "$URL/debug/client-ip"; echo
curl -s -H "X-Forwarded-For: 1.2.3.4" "$URL/debug/client-ip"; echo
```

| What you see | Meaning | Action |
|---|---|---|
| `rate_limit_key` equals your IP in **both** calls, and the second call's `raw_x_forwarded_for` is `1.2.3.4, <your IP>` | Correct | Nothing |
| `rate_limit_key` is a Google address (e.g. `35.191.x.x`), with your IP just before it | Google added two hops | Set `BECKMANN_TRUSTED_PROXIES=2` in the workflow, redeploy, re-check |
| `rate_limit_key` is `1.2.3.4` in the second call | Trusting too many hops (forgeable) | Lower the value |
| `rate_limit_key` is `169.254.x.x`/`127.0.0.1`, and `trusted_proxies` is `0` | The env var isn't set | Check the deploy command |

Turn it off again:
```bash
gcloud run services update beckmann-pyscf --region us-central1 \
  --remove-env-vars BECKMANN_DEBUG_CLIENT_IP
```
The next GitHub deploy would remove it anyway, since the workflow replaces the full env var list.

**c) A real prediction.** Open the Pages site, run α-tetralone (`O=C1CCCC2=CC=CC=C21`), and keep the tab open. Expect **R** after roughly 20–35 min on 1 vCPU. Then look at:
- **Usage:** Actions → Free-tier guard → run it by hand → the summary shows billable seconds used. One scan should be ~1,200–2,100 s of the 67,500 s allowance.
- **Memory:** Cloud Run console → beckmann-pyscf → **Metrics** → *Memory utilization*. Expect a peak around 65%. If the logs show `Memory limit of 4096 MiB exceeded`, see Troubleshooting.

---

## Known limitations of this configuration

- **About 35–45 scans a month**, then predictions pause until the 1st (UTC). That's the price of staying inside the free tier.
- **One scan at a time.** A second visitor's prediction gets "The server is busy with another prediction" until the first finishes. Browsing the page still works, since it's on GitHub Pages.
- **Closing the tab cancels the scan.** There's no resume, and a dropped network connection does the same.
- **Cold starts.** After ~15 min idle, Cloud Run scales to zero, and the next prediction waits for the container to start (billed as part of that request).
- **Docker Hub pull limits.** Cloud Run pulls the image from Docker Hub when an instance starts. Docker Hub rate-limits anonymous pulls per IP, and Google's pull IPs are shared, so a cold start could occasionally fail with an image-pull error and work on retry. If that happens often, the fix (an Artifact Registry remote repository) costs a little storage money. Note it and decide then.

## Troubleshooting

| Symptom | Likely cause / fix |
|---|---|
| Workflows skipped | `GCP_PROJECT_ID` / `DOCKERHUB_USERNAME` variables not set (step 8) |
| Docker Hub login/push fails | `DOCKERHUB_TOKEN` secret missing or not Read & Write; the repo name must be `beckmann-pyscf` |
| `auth` step: "unauthorized_client" or "Permission 'iam.serviceAccounts.getAccessToken' denied" | WIF condition mismatch: run from `main`, check `GITHUB_REPO` in step 5, wait a few minutes |
| Deploy fails: "iam.serviceaccounts.actAs" | Step 4's `serviceAccountUser` binding on the runtime SA |
| Deploy fails setting IAM policy for allUsers | Org policy (domain restricted sharing). Use a personal-account project |
| Guard fails reading metrics (403) | Step 4's `monitoring.viewer` for the deployer |
| Every prediction says "temporarily unavailable" | The gate can't read usage: check `monitoring.viewer` for the runtime SA, and the Cloud Run logs for `usage gate:` lines |
| Page loads but predictions fail with a CORS error in the browser console | `BECKMANN_CORS_ORIGINS` doesn't match the Pages origin. Set the `PAGES_URL` variable and redeploy |
| "This took longer than expected" | Scan exceeded 3300 s. Raising `BECKMANN_JOB_TIMEOUT` (max ~3500) also raises each scan's reserved worst case, so fewer scans fit |
| "Memory limit exceeded" in logs | Raise `--memory` to `8Gi` **and** `CLOUD_RUN_MEMORY_GIB` to `8` in both workflows. That halves the monthly scan allowance, because memory is the binding free-tier limit |
| Budget email arrived | Billing → Reports shows what was charged. To stop everything immediately: "Turning it off" |

## Turning it off

```bash
# Make it private (reversible; the next deploy or guard run may reopen it if usage allows):
gcloud run services remove-iam-policy-binding beckmann-pyscf --region us-central1 \
  --member=allUsers --role=roles/run.invoker
# Delete the service (stops all Cloud Run charges):
gcloud run services delete beckmann-pyscf --region us-central1
# Or remove the whole project (recoverable for 30 days):
gcloud projects delete "$PROJECT_ID"
```
To stop it coming back, also disable the three workflows in the Actions tab.
