#!/usr/bin/env bash
# tcren bootstrap — source of truth for a reproducible install.
#
# Steps:
#   1. Create/update the `tcren` conda environment (python + mmseqs2 + toolchain).
#   2. pip install -e . — pulls in arda (pinned git tag) and builds its C++ extension.
#
# Flags:
#   --no-conda   Skip conda env creation (use the already-active environment).
#   --tests      After install, run the fast test suites.
#
# Usage:
#   bash setup.sh [--no-conda] [--tests]
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="tcren"
USE_CONDA=1
DO_TESTS=0

for arg in "$@"; do
  case "$arg" in
    --no-conda) USE_CONDA=0 ;;
    --tests)    DO_TESTS=1 ;;
    *) echo "Unknown flag: $arg" >&2; exit 2 ;;
  esac
done

log() { printf '\033[1;34m[tcren]\033[0m %s\n' "$*"; }

# --- 1. conda environment --------------------------------------------------
if [[ "$USE_CONDA" -eq 1 ]]; then
  if ! command -v conda >/dev/null 2>&1; then
    echo "conda not found on PATH; install miniconda/anaconda or pass --no-conda." >&2
    exit 1
  fi
  if conda env list | grep -qE "^${ENV_NAME}\s"; then
    log "updating conda env '${ENV_NAME}'"
    conda env update -n "$ENV_NAME" -f "$ROOT/environment.yml" --prune
  else
    log "creating conda env '${ENV_NAME}'"
    conda env create -f "$ROOT/environment.yml"
  fi
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "$ENV_NAME"
fi

# --- 2. tcren (+ arda, pinned git tag, from pyproject.toml) -----------------
log "installing tcren (editable); arda is pulled in as a pinned git dependency"
pip install -e "$ROOT"

# --- 3. reference data (HF) ------------------------------------------------
# Populate data/ with Native2026 (orientation refs) + Canonical2026 (the default
# `tcren superimpose` database). Skips folders already present. Set TCREN_NO_FETCH=1 to skip.
if [[ "${TCREN_NO_FETCH:-0}" -ne 1 ]]; then
  log "fetching reference structure sets into data/"
  tcren fetch-data
fi

# --- 4. tests --------------------------------------------------------------
if [[ "$DO_TESTS" -eq 1 ]]; then
  log "running fast tests"
  pytest "$ROOT/tests" -q
fi

log "done. Activate with: conda activate ${ENV_NAME}"
