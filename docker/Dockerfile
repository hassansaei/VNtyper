# syntax=docker/dockerfile:1.7
# =============================================================================
# VNtyper application image - rebuilt on every commit
# =============================================================================
# Layers onto ghcr.io/<owner>/vntyper-base, which carries the conda environments,
# adVNTR and the reference genomes. Everything expensive lives there; this file only
# installs the application, so it builds in ~3 minutes instead of 40-70.
#
# BASE_IMAGE is passed by .github/workflows/docker-build.yml as a content-addressed
# tag derived from the base image's inputs, so a commit that changes conda/** cannot
# accidentally build against a stale base - the tag simply will not exist and the
# build fails loudly.
#
# Build locally with:
#   docker build -f docker/Dockerfile.base -t vntyper-base:local .
#   docker build -f docker/Dockerfile --build-arg BASE_IMAGE=vntyper-base:local -t vntyper:local .
# =============================================================================

ARG BASE_IMAGE=ghcr.io/hassansaei/vntyper-base:latest

FROM ${BASE_IMAGE} AS runtime

LABEL org.opencontainers.image.title="VNtyper" \
      org.opencontainers.image.description="VNtyper 2 - MUC1 VNTR genotyping pipeline for ADTKD-MUC1" \
      org.opencontainers.image.vendor="VNtyper Project" \
      org.opencontainers.image.licenses="BSD-3-Clause" \
      org.opencontainers.image.source="https://github.com/hassansaei/VNtyper" \
      org.opencontainers.image.url="https://github.com/hassansaei/VNtyper" \
      org.opencontainers.image.documentation="https://hassansaei.github.io/VNtyper/" \
      com.vntyper.type="bioinformatics" \
      com.vntyper.pipeline="genotyping" \
      com.vntyper.build="app"

WORKDIR /opt/vntyper

# --- application package -----------------------------------------------------
# Source comes from the build context, NOT a `git clone` of GitHub main. The old
# Dockerfile cloned the repo, which meant (a) PR builds validated main rather than the
# PR, and (b) the RUN layer's cache key never changed, so a cache hit could serve an
# arbitrarily stale checkout.
#
# --no-deps is required: dependencies are conda-managed, and pyproject.toml pins
# numpy>=1.26.0,<2.0.0 while conda/environment_vntyper.yml installs numpy=2.0.2. Without
# --no-deps, pip pulls a second numpy from PyPI over the conda-managed one.
COPY pyproject.toml README.md ./
COPY vntyper/ ./vntyper/
RUN conda run -n vntyper pip install --no-cache-dir --no-deps . && \
    conda run -n vntyper vntyper --version

# NOTE: deliberately NO `install-references` step here.
#
# The reference paths the pipeline actually reads - bwa_reference_hg19/hg38,
# advntr_reference_vntr_*, muc1_*, code_adVNTR_RUs - are already correct RELATIVE paths
# in the shipped config.json, and WORKDIR is /opt/vntyper with the references at
# /opt/vntyper/reference. install_references.py's `--config-path` only adds `ucsc_*` /
# `ncbi_*` absolute keys, which nothing in vntyper/ reads.
#
# Re-running it here is also actively harmful: `download_file` skips existing files but
# decompression and MD5 run unconditionally, so it re-downloaded 148 MB and rewrote
# ~508 MB of .fa on every application build - a 680 MB duplicated layer on top of the
# base, for no change in behaviour.
#
# This also matches the previous single-stage image, where install-references patched the
# source-tree config at /opt/vntyper/vntyper/config.json, never the installed copy that
# cli.py actually loads via pkg_resources.

# --- most volatile last ------------------------------------------------------
COPY --chown=appuser:appuser docker/app/ /opt/vntyper/app/
COPY docker/entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

EXPOSE 8000

USER appuser

HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
    CMD ["/usr/local/bin/entrypoint.sh", "health"]

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
