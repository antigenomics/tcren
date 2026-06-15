#!/usr/bin/env bash
# Set up the `tcren-nb` conda env for the Nat Comput Sci 2022 reproduction notebooks.
#
#   1. Create/update the conda env from environment.yml (python + mmseqs2 + pip).
#   2. Editable-install the sibling `arda` package and this `tcren` repo (with the
#      `notebooks` extra: jupyter, matplotlib, rapidfuzz, scikit-learn, logomaker).
#   3. Register a Jupyter kernel so VS Code / Jupyter can select it.
#
# Usage:  bash setup.sh
# Env overrides:  ENV_NAME (default tcren-nb), ARDA_DIR (default ../../../arda).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"                       # tcren-ms repo root
ARDA_DIR="${ARDA_DIR:-$(cd "$ROOT/../arda" && pwd)}"
ENV_NAME="${ENV_NAME:-tcren-nb}"

echo ">> Creating/updating conda env '$ENV_NAME' from environment.yml"
conda env create -n "$ENV_NAME" -f "$HERE/environment.yml" 2>/dev/null \
  || conda env update -n "$ENV_NAME" -f "$HERE/environment.yml"

echo ">> Editable-installing arda ($ARDA_DIR)"
conda run -n "$ENV_NAME" pip install -e "$ARDA_DIR"

echo ">> Editable-installing tcren + notebook deps ($ROOT)"
conda run -n "$ENV_NAME" pip install -e "$ROOT[notebooks]"

echo ">> Registering Jupyter kernel"
conda run -n "$ENV_NAME" python -m ipykernel install --user \
  --name "$ENV_NAME" --display-name "Python (tcren-nb)"

echo ">> Done. Select the 'Python (tcren-nb)' kernel in VS Code / Jupyter."
