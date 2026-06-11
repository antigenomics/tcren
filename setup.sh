#!/usr/bin/env bash
# tcren bootstrap — source of truth for a reproducible install.
#
# Steps:
#   1. Create/update the `tcren` conda environment (python + mmseqs2 + toolchain).
#   2. Install the sibling `arda` package in editable mode (TCR annotation backend).
#   3. pip install -e . (this package).
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
ARDA_DIR="${ARDA_DIR:-$(cd "$ROOT/../arda" 2>/dev/null && pwd || true)}"
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

# --- 2. arda (TCR annotation backend) --------------------------------------
# By default arda is installed from its dev branch via environment.yml (git+https).
# For local co-development, set ARDA_DIR to a checkout to install it editable instead.
if [[ -n "${ARDA_DIR:-}" && -f "$ARDA_DIR/pyproject.toml" ]]; then
  log "installing arda (editable) from $ARDA_DIR (overrides the git@dev install)"
  pip install -e "$ARDA_DIR"
else
  log "arda is installed from git@dev via environment.yml"
fi

# --- 3. tcren --------------------------------------------------------------
log "installing tcren (editable)"
pip install -e "$ROOT"

# --- 4. tests --------------------------------------------------------------
if [[ "$DO_TESTS" -eq 1 ]]; then
  log "running fast tests"
  pytest "$ROOT/tests" -q
fi

log "done. Activate with: conda activate ${ENV_NAME}"
