#!/usr/bin/env bash
# tcren bootstrap — source of truth for a reproducible install (uv, no conda).
#
# Steps:
#   1. Create a repo-local .venv with uv.
#   2. uv pip install -e . — builds the C++ extensions (scikit-build-core fetches
#      cmake+ninja) and pulls arda (PyPI: arda-mapper) + the rest from PyPI.
#   3. Fetch the reference structure sets into data/.
#
# mmseqs2 is NOT installed here: arda auto-fetches a static mmseqs binary on first use
# (override with $ARDA_MMSEQS / disable with $ARDA_NO_AUTO_FETCH). The only host requirement
# is a C++ compiler for tcren's own extensions.
#
# Flags:
#   --tests   After install, run the fast (non-slow) test suite.
#
# Usage (from zsh or bash):
#   bash setup.sh [--tests]
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DO_TESTS=0

for arg in "$@"; do
  case "$arg" in
    --tests) DO_TESTS=1 ;;
    *) echo "Unknown flag: $arg" >&2; exit 2 ;;
  esac
done

log() { printf '\033[1;34m[tcren]\033[0m %s\n' "$*"; }

# --- 0. prerequisites ------------------------------------------------------
if ! command -v uv >/dev/null 2>&1; then
  echo "uv not found on PATH. Install it, e.g.:" >&2
  echo "  curl -LsSf https://astral.sh/uv/install.sh | sh" >&2
  exit 1
fi
if ! { command -v cc >/dev/null 2>&1 || command -v clang >/dev/null 2>&1 \
       || command -v gcc >/dev/null 2>&1; }; then
  echo "No C++ compiler found (needed to build tcren's extensions)." >&2
  echo "  macOS:  xcode-select --install" >&2
  echo "  Debian: sudo apt-get install build-essential" >&2
  exit 1
fi

# --- 1. venv ---------------------------------------------------------------
if [[ ! -d "$ROOT/.venv" ]]; then
  log "creating .venv with uv"
  uv venv "$ROOT/.venv"
fi
# shellcheck disable=SC1091
source "$ROOT/.venv/bin/activate"

# --- 2. tcren (+ arda-mapper, from pyproject.toml) -------------------------
log "installing tcren (editable); arda (PyPI: arda-mapper) is pulled in as a dependency"
log "arda auto-fetches a static mmseqs binary on first use — no conda/bioconda needed"
uv pip install -e "$ROOT"

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
  pytest "$ROOT/tests" -m "not slow" -q
fi

log "done. Activate with: source .venv/bin/activate"
