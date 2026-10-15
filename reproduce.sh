#!/usr/bin/env bash
# Run the five fit-free benchmark blocks from raw inputs through the benchmark CLI.
# Requires bootstrap inputs, a pinned interpreter and a new output directory:
#   ./reproduce.sh --output /absolute/new/run --baseline-spec /absolute/baseline.json
# TCREN2_BENCH selects the benchmark checkout; TCREN_REPRO_PYTHON selects its interpreter.
# Use --help for the current driver. --historical explicitly selects the archived fitted
# manuscript workflow below; its --list/--library-only/subset flags are diagnostic only.
# No computed intermediates are resumed or reused by the final-model route.
set -uo pipefail

# Final reproduction has its own explicit raw-input driver. The old path below
# reproduces fitted historical evidence and is never a final-model default.
if [[ ${1-} != --historical ]]; then
  BENCH="${TCREN2_BENCH:-$HOME/vcs/projects/2026-tcren2-code}"
  PY="${TCREN_REPRO_PYTHON:-$BENCH/.venv/bin/python}"
  [[ -x "$PY" ]] || { echo "missing reproduction interpreter: $PY" >&2; exit 1; }
  exec "$PY" "$BENCH/bench/scripts/reproduce_fit_free.py" "$@"
fi
shift

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"                # the pytest node ids below are repo-relative
OUT="$ROOT/reproduce"
TCREN2_BENCH="${TCREN2_BENCH:-$HOME/vcs/projects/2026-tcren2-code}"
TCREN2_MS="${TCREN2_MS:-$HOME/vcs/manuscripts/2026-tcren2-ms}"

ALL_STAGES=(env potential model scores bench paper)
LIBRARY_STAGES=(env potential model scores)

log()  { printf '\033[1;34m[reproduce]\033[0m %s\n' "$*"; }
die()  { printf '\033[1;31m[reproduce] FAILED\033[0m %s\n' "$*" >&2; exit 1; }

usage() {
  cat <<'EOF'
stages:
  env         record the version, the descriptor-catalogue digest and the data root        1 s
  potential   re-derive TCRen2 from Native2026 and assert it reproduces the shipped      24 s
              matrix; emit the derivation's own contact table and count it
  model       validate the frozen object and catalogue (no fitting); the benchmark
              rebuilds the full hold-out from raw inputs and checks exact frozen identity
  scores      featurise all 374 Native2026 crystals and run the whole score set on       80 s
              them, live, through the two released commands
  bench       $TCREN2_BENCH/recompute.sh — every fixture, table and .dat the paper reads
  paper       $TCREN2_MS/make — the figures, then main.pdf and supplementary.pdf
EOF
}

WANT=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --list|-l)      usage; exit 0 ;;
    -h|--help)      awk 'NR>1 { if (/^#/) { sub(/^# ?/, ""); print } else { exit } }' "${BASH_SOURCE[0]}"; echo; usage; exit 0 ;;
    --library-only) WANT=("${LIBRARY_STAGES[@]}"); shift ;;
    -*)             die "unknown flag: $1 (try --help)" ;;
    *)
      # shellcheck disable=SC2076
      [[ " ${ALL_STAGES[*]} " == *" $1 "* ]] || die "no such stage: $1 (try --list)"
      WANT+=("$1"); shift ;;
  esac
done
[[ ${#WANT[@]} -gt 0 ]] || WANT=("${ALL_STAGES[@]}")

wants() { local s; for s in "${WANT[@]}"; do [[ "$s" == "$1" ]] && return 0; done; return 1; }

# --- preconditions, checked before anything runs ------------------------------------------------
# Every one of these is a `./bootstrap.sh` output. Checking them here rather than letting a stage
# die halfway through means the message names the missing input instead of a downstream symptom.
[[ -x "$ROOT/.venv/bin/python" ]] || die "no venv at $ROOT/.venv — run ./bootstrap.sh first"
# shellcheck disable=SC1091
source "$ROOT/.venv/bin/activate"
command -v tcren >/dev/null 2>&1 || die "the tcren console script is not on PATH — run ./bootstrap.sh"

DATA_DIR="$(python -c 'from tcren.paths import data_dir; print(data_dir())')" \
  || die "tcren does not import from $ROOT/.venv — run ./bootstrap.sh"
MHC_FASTA="$(python -c 'from tcren.mhc.reference import DATABASE_DIR; print(DATABASE_DIR / "alleles.aa.fasta")')"

if wants potential || wants scores; then
  [[ -d "$DATA_DIR/Native2026" ]] \
    || die "no Native2026 at $DATA_DIR/Native2026 — run ./bootstrap.sh (or tcren fetch-data)"
  [[ -s "$MHC_FASTA" ]] \
    || die "no MHC allele reference at $MHC_FASTA — run ./bootstrap.sh (or tcren build-mhc-ref)"
fi

mkdir -p "$OUT"
T0=$SECONDS
FAILED=()

run_stage() {           # run_stage <name> <function>
  local name="$1" fn="$2" t=$SECONDS
  log "--- $name ---"
  if "$fn"; then
    printf '\033[1;32m[reproduce]\033[0m %s ok (%ds)\n' "$name" "$((SECONDS - t))"
  else
    FAILED+=("$name")
    printf '\033[1;31m[reproduce]\033[0m %s FAILED (%ds)\n' "$name" "$((SECONDS - t))" >&2
    return 1
  fi
}

# --- env ----------------------------------------------------------------------------------------
# The two provenance stamps that decide whether any later number is comparable: the released
# version and the SHA-256 digest of the descriptor catalogue. A feature table is only valid under
# the catalogue that wrote it, and `tcren assess` refuses one written under a different digest.
stage_env() {
  {
    tcren info
    python - <<'PY'
from tcren.paths import data_dir, tcren_home
from tcren.provenance import registry_digest
from tcren.score import holdout_manifest
m = holdout_manifest()
print(f"catalogue digest: {registry_digest()}")
print(f"hold-out manifest: {m.height} structures, {m['epitope'].n_unique()} epitopes, "
      f"{int(m['y'].sum())} binders")
print(f"tcren_home: {tcren_home()}")
print(f"data_dir: {data_dir()}")
PY
  } | tee "$OUT/env.txt"
}

# --- potential ----------------------------------------------------------------------------------
# tests/regression/test_shipped_potentials.py re-derives every shipped matrix from the recipe in
# data/potentials.json and asserts it reproduces the file. That is the derivation reproduction; it
# is a test rather than a bare command so that a mismatch is an assertion and not a diff a reader
# has to eyeball. --contacts-out then emits the contact table the matrix was actually derived from,
# AFTER the alpha-beta filter and the withdrawn-entry filter, so the counts a downstream analysis
# quotes describe the derivation's own set.
stage_potential() {
  RUN_BENCHMARK=1 python -m pytest tests/regression/test_shipped_potentials.py -q -p no:cacheprovider \
    || return 1
  tcren derive-potential \
      --structure-dir "$DATA_DIR/Native2026" \
      --balance both \
      -o "$OUT/TCRen2_check.csv" \
      --contacts-out "$OUT/derivation_contacts.csv" \
    || return 1
  python - "$OUT/derivation_contacts.csv" "$OUT/TCRen2_check.csv" <<'PY' | tee "$OUT/potential.txt"
import sys
import polars as pl
from tcren.potential import tcren2, withdrawn_entries

contacts = pl.read_csv(sys.argv[1])
check = pl.read_csv(sys.argv[2])
sid = next(c for c in contacts.columns if "pdb" in c.lower() or c.lower() in ("id", "structure"))
print(f"derivation contacts: {contacts.height} rows over {contacts[sid].n_unique()} structures")
print(f"withdrawn entries dropped: {withdrawn_entries().height} listed")
print(f"re-derived matrix: {check.height} cells; shipped: {tcren2().matrix.height} cells")
if check.height != tcren2().matrix.height:
    sys.exit("re-derived matrix has a different number of cells than the shipped one")
PY
}

# --- model --------------------------------------------------------------------------------------
# Validate frozen inputs here. The benchmark rebuilds every hold-out descriptor from raw
# structures, then scripts/check_fresh_holdout.py fits a temporary object and compares it to
# the frozen arrays before any publication producers. No committed derived slice is reused.
stage_model() {
  python - <<'PY' | tee "$OUT/model.txt"
from tcren.score import holdout_model, holdout_manifest
m = holdout_model()
print(f"frozen model validated: {m.n_pos} binders, {m.n_neg} non-binders")
print(f"manifest: {holdout_manifest().height} structures")
print("No refit here. Full raw-input hold-out comparison is required in the benchmark stage.")
PY
}

# --- scores -------------------------------------------------------------------------------------
# The whole read-out, live, through the two released commands rather than from any cached table.
# `features` annotates the set in ONE arda call and ONE mmseqs MHC search — never one per
# structure — so this is an 80-second stage rather than an hour-long one.
stage_scores() {
  tcren features -s "$DATA_DIR/Native2026" \
      --all \
      -t 0 -o "$OUT/native_features.tsv" \
    || return 1
  tcren assess --features "$OUT/native_features.tsv" -o "$OUT/native_assess.tsv" \
    || return 1
  python - "$OUT/native_features.tsv" "$OUT/native_assess.tsv" <<'PY' | tee "$OUT/scores.txt"
import sys
import numpy as np
import polars as pl

feats = pl.read_csv(sys.argv[1], separator="\t", infer_schema_length=None)
sc = pl.read_csv(sys.argv[2], separator="\t", infer_schema_length=None)
print(f"features: {feats.height} structures x {feats.width} columns")
print(f"assess:   {sc.height} structures x {sc.width} columns")
required = ["pose_score", "binder_score", "peptide_score", "S",
            "channel_placement", "channel_interface", "channel_shape", "channel_energetics"]
missing = [c for c in required if c not in sc.columns]
if missing:
    sys.exit(f"assess table is missing read-outs: {', '.join(missing)}")
if sc.height != feats.height:
    sys.exit(f"assess wrote {sc.height} rows for {feats.height} featurised structures")
for c in required:
    finite = int(np.isfinite(sc[c].to_numpy().astype(float)).sum())
    print(f"  {c:<20} {finite:>4} finite, {sc.height - finite:>4} missing of {sc.height}")
    if not finite:
        sys.exit(f"assess read-out {c} has no finite observations")
PY
}

# --- bench --------------------------------------------------------------------------------------
stage_bench() {
  [[ -x "$TCREN2_BENCH/recompute.sh" ]] \
    || { echo "no recompute.sh at $TCREN2_BENCH — set \$TCREN2_BENCH to the benchmark repository" >&2
         return 1; }
  ( cd "$TCREN2_BENCH" && ./recompute.sh )
}

# --- paper --------------------------------------------------------------------------------------
stage_paper() {
  [[ -f "$TCREN2_MS/Makefile" ]] \
    || { echo "no Makefile at $TCREN2_MS — set \$TCREN2_MS to the manuscript repository" >&2
         return 1; }
  ( cd "$TCREN2_MS" && make )
}

for s in "${ALL_STAGES[@]}"; do
  wants "$s" || continue
  run_stage "$s" "stage_$s" || break
done

log "total $((SECONDS - T0))s; artifacts in $OUT"
if [[ ${#FAILED[@]} -gt 0 ]]; then
  die "stage(s) failed: ${FAILED[*]}"
fi
if [[ ${#WANT[@]} -eq ${#ALL_STAGES[@]} ]]; then
  log "every stage green"
else
  log "requested diagnostic subset passed; this is not a complete reproduction"
fi
