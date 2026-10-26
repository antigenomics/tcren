#!/usr/bin/env bash
# Delegate a dataset-specific reproduction to an explicitly selected external driver.
set -euo pipefail
if [[ ${1-} == --help || ${1-} == -h ]]; then
  cat <<'EOF'
Run a dataset-specific reproduction with an immutable installed tcren wheel:
  TCREN_REPRO_PYTHON=/absolute/environment/bin/python \
  TCREN_REPRO_DRIVER=/absolute/benchmark/reproduce.py \
  ./reproduce.sh [driver arguments]

The external driver owns raw inputs, cohorts and publication outputs. It should
invoke tcren CLI commands and preserve wheel/input hashes. No dataset is assumed,
no fitted historical workflow is selected, and no calculation runs for --help.
EOF
  exit 0
fi
: "${TCREN_REPRO_PYTHON:?Set the absolute path to the pinned reproduction interpreter}"
: "${TCREN_REPRO_DRIVER:?Set the absolute path to the external reproduction driver}"
[[ $TCREN_REPRO_PYTHON == /* && -x $TCREN_REPRO_PYTHON ]] || { echo 'Invalid reproduction interpreter' >&2; exit 2; }
[[ $TCREN_REPRO_DRIVER == /* && -f $TCREN_REPRO_DRIVER ]] || { echo 'Invalid reproduction driver' >&2; exit 2; }
exec "$TCREN_REPRO_PYTHON" "$TCREN_REPRO_DRIVER" "$@"
