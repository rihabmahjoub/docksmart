# --- Stage 1: compile fpocket from source -----------------------------------
# fpocket has no pip/PyPI distribution, but it's a small, fast-compiling C
# codebase (verified: ~1.9MB binary, builds cleanly with just gcc/g++/make —
# no BLAS/LAPACK/netcdf actually required despite what some docs imply).
# Compiling it ourselves on debian-slim avoids pulling in the full
# condaforge/mambaforge base image (several hundred MB before installing a
# single package), which matters a lot on Render's 512MB-RAM free tier.
FROM debian:bookworm-slim AS fpocket-builder
RUN apt-get update && apt-get install -y --no-install-recommends \
        gcc g++ make git ca-certificates \
    && rm -rf /var/lib/apt/lists/*
RUN git clone --depth 1 https://github.com/Discngine/fpocket.git /opt/fpocket-src \
    && cd /opt/fpocket-src && make

# --- Stage 2: the actual application image -----------------------------------
FROM python:3.11-slim

WORKDIR /app

# Runtime-only system deps: nothing heavy. libstdc++ is needed to run the
# fpocket binary we're copying in from the builder stage.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libstdc++6 \
    && rm -rf /var/lib/apt/lists/*

COPY --from=fpocket-builder /opt/fpocket-src/bin/fpocket /usr/local/bin/fpocket

COPY requirements.txt .
# All of these are verified-working pure pip wheels — no conda needed:
#   - rdkit, openmm, pdbfixer, meeko, vina: standard manylinux wheels
#   - openbabel: has its own pip wheel (via the 'openbabel' PyPI package)
#   - pymol-open-source-whl: a self-contained ~15MB wheel (no Qt/GStreamer),
#     pins numpy==1.26.4, which the rest of the stack is compatible with.
# meeko has two UNDECLARED runtime dependencies (scipy, gemmi) that pip
# won't pull in automatically — they're listed explicitly in requirements.txt
# to avoid a ModuleNotFoundError at job-run time instead of at build time.
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

RUN mkdir -p /app/data/jobs

EXPOSE 8000
CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000"]
