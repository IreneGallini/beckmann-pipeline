# Not built/tested in the environment this was written in (no Docker
# available there). Its first real build is the GitHub Actions deploy
# (.github/workflows/deploy-cloud-run.yml, which pushes it to Docker Hub)
# -- check that run's log, or
# build it yourself, before trusting it. Deployment guide: DEPLOY.md.
#
# Build (from the monorepo root, so both packages/beckmann-core/ and
# packages/beckmann-pyscf/ are in the build context):
#   docker build -f packages/beckmann-pyscf/backend/Dockerfile -t beckmann-pyscf .
#
# Run locally, with the same limits Cloud Run uses:
#   docker run -p 8080:8080 --memory=4g --cpus=1 beckmann-pyscf
#   -> http://localhost:8080
#
# Defaults to sync mode (BECKMANN_API_MODE=sync): one gunicorn process
# runs each scan inside its own /predict request -- what Cloud Run's
# request-based billing needs (see DEPLOY.md). Setting BECKMANN_REDIS_URL
# switches entrypoint.sh to the async multi-worker layout instead
# (redis-server + dispatcher + gunicorn -w WEB_CONCURRENCY), for hosts
# that keep CPU on between requests; set BECKMANN_API_MODE=async with it.

FROM python:3.11-slim

# redis-server: the shared job store / rate-limit storage (in-container,
#   localhost-only -- see entrypoint.sh).
# libgomp1: OpenMP runtime used under the hood by numpy/pyscf/torch's CPU
#   threading -- not always present on a slim base image.
# libxrender1/libxext6: rdkit's drawing code links against these.
RUN apt-get update && apt-get install -y --no-install-recommends \
    redis-server \
    libgomp1 \
    libxrender1 \
    libxext6 \
    && rm -rf /var/lib/apt/lists/*

# Run as an unprivileged user, not root.
RUN useradd -m -u 1000 user

# Where AIMNet2 caches its downloaded model weights (aimnet's
# model_registry.get_cache_dir()). Filled at build time below, so no job
# ever waits on -- or fails on -- a network download.
ENV AIMNET_CACHE_DIR=/opt/aimnet-cache \
    MPLCONFIGDIR=/tmp/matplotlib \
    KMP_DUPLICATE_LIB_OK=TRUE \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

WORKDIR /app

# CPU-only torch first, from PyTorch's own CPU wheel index: the default
# PyPI torch wheel bundles CUDA (several GB), useless on CPU-only hosting.
# Same version beckmann-core pins, so the install below is a no-op for it.
RUN pip install torch==2.12.0 --index-url https://download.pytorch.org/whl/cpu

COPY packages/beckmann-core /app/packages/beckmann-core
COPY packages/beckmann-pyscf /app/packages/beckmann-pyscf

# Real (non-editable) installs -- pulls in each package's pinned dependency
# versions from its own pyproject.toml. beckmann-core installed first so
# beckmann-pyscf's "beckmann-core==<pinned>" requirement resolves against
# the local copy just installed, not a PyPI lookup.
RUN pip install /app/packages/beckmann-core /app/packages/beckmann-pyscf

# Bake the AIMNet2 weights into the image (Auto3D's own AIMNET weights ship
# inside its wheel already).
RUN mkdir -p "$AIMNET_CACHE_DIR" \
    && python -c "from aimnet.calculators import AIMNet2Calculator; AIMNet2Calculator('aimnet2_2025')" \
    && chown -R user:user "$AIMNET_CACHE_DIR"

RUN chmod +x /app/packages/beckmann-pyscf/backend/entrypoint.sh

USER user
ENV HOME=/home/user
WORKDIR /app/packages/beckmann-pyscf/backend

# Cloud Run injects PORT (8080 by default) at runtime, overriding this;
# 8080 here just makes a plain `docker run` match it.
# PYSCF_MAX_MEMORY: PySCF's own per-calculation memory budget in MB
# (default 4000). Caps in-core ERI storage so a small molecule's SCF can't
# try to take 4 GB inside a 4 GiB container next to torch/AIMNet2.
ENV PORT=8080 \
    BECKMANN_API_MODE=sync \
    PYSCF_MAX_MEMORY=2000 \
    WEB_CONCURRENCY=2
EXPOSE 8080

HEALTHCHECK --interval=30s --timeout=5s --start-period=60s \
    CMD python -c "import os, urllib.request; urllib.request.urlopen(f'http://127.0.0.1:{os.environ[\"PORT\"]}/health', timeout=4)"

CMD ["./entrypoint.sh"]
