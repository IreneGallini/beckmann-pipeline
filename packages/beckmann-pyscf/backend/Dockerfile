# Not built/tested in the environment this was written in (no Docker
# available there) -- written to best practice, build and run this
# yourself before trusting it in production.
#
# Build (from the monorepo root, so both packages/beckmann-core/ and
# packages/beckmann-pyscf/ are in the build context):
#   docker build -f packages/beckmann-pyscf/backend/Dockerfile -t beckmann-pyscf .
#
# Run -- size the container to comfortably fit settings.MAX_CONCURRENT_JOBS
# simultaneous AIMNet2/PySCF jobs, not just one; resource limits are a
# deployment-time flag, not baked into the image:
#   docker run -p 5001:5001 --memory=4g --cpus=2 \
#     -e BECKMANN_MAX_CONCURRENT=2 beckmann-pyscf

FROM python:3.11-slim

# libgomp1: OpenMP runtime used under the hood by numpy/pyscf/torch's CPU
# threading -- not always present on a slim base image.
RUN apt-get update && apt-get install -y --no-install-recommends \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY packages/beckmann-core /app/packages/beckmann-core
COPY packages/beckmann-pyscf /app/packages/beckmann-pyscf

# Real (non-editable) installs -- pulls in each package's pinned dependency
# versions from its own pyproject.toml. beckmann-core installed first so
# beckmann-pyscf's "beckmann-core==<pinned>" requirement resolves against
# the local copy just installed, not a PyPI lookup.
RUN pip install --no-cache-dir \
    /app/packages/beckmann-core \
    /app/packages/beckmann-pyscf

# Required on macOS builds of this image; harmless elsewhere. Matches what
# every package's own Python modules already set at import time -- set
# here too so it's correct even before Python starts.
ENV KMP_DUPLICATE_LIB_OK=TRUE

WORKDIR /app/packages/beckmann-pyscf/backend

EXPOSE 5001

# See wsgi.py's own docstring: -w 1 is required, not a default left
# unconfigured -- more than one gunicorn worker breaks the in-memory job
# store's /submit-then-/status polling.
CMD ["gunicorn", "-w", "1", "-b", "0.0.0.0:5001", "wsgi:app"]
