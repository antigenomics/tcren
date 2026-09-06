#!/usr/bin/env bash
# tcren — reproduce every number the TCRen2 manuscript prints. 2026-09-05
#
# Run ./bootstrap.sh FIRST. That script owns the environment and the ~620 MB of structures, and it
# is deliberately outside any time budget. This one assumes both are already on disk and fails
# immediately, by name, if they are not.
#
#   ./reproduce.sh                  every stage, in order
#   ./reproduce.sh --list           the stages and what each one costs
#   ./reproduce.sh --library-only   stop after the four stages this repository owns
#   ./reproduce.sh potential model  run these stages only (repeatable)
#
# THE THREE-REPOSITORY BOUNDARY. Every descriptor, potential, score and fit lives HERE. The
# evaluation (ROC/PR/CI, the cohort tables, the ledger) lives in the benchmark repository and the
# LaTeX lives in the manuscript repository. This script therefore runs four stages itself and then
# hands off, locating the other two by environment variable:
#
#   TCREN2_BENCH   default ~/vcs/projects/2026-tcren2-code   its ./recompute.sh
#   TCREN2_MS      default ~/vcs/manuscripts/2026-tcren2-ms  its make
#
# Point them somewhere else and the handoff follows; if either is absent the `bench` / `paper`
# stage EXITS NON-ZERO naming the path it looked for. It does not skip.
#
# THERE IS NO --skip-slow FLAG, on purpose. The benchmark repository has one and its own header
# claims it drops only stages "no reported table reads" — which is false for four of the eight it
# names. A flag that quietly changes which numbers a run produces is worse than a run that takes
# longer, so this script has one path through it. To run less, name the stages you want.
#
# WHAT IT COSTS. `--library-only` is 106 s end to end, measured 2026-09-05 on macOS arm64 against
# the 374-structure Native2026 set: env 1 s, potential 24 s, model 1 s, scores 80 s. The 30-minute
# target is therefore set by the `bench` stage, which is the benchmark repository's own recompute
# and is timed there. This script prints elapsed seconds for every stage it runs and a total at the
# end, so a stall is visible rather than inferred from silence.
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"                # the pytest node ids below are repo-relative
OUT="$ROOT/reproduce"
TCREN2_BENCH="${TCREN2_BENCH:-$HOME/vcs/projects/2026-tcren2-code}"
TCREN2_MS="${TCREN2_MS:-$HOME/vcs/manuscripts/2026-tcren2-ms}"

ALL_STAGES=(env potential model scores contacts bench paper)
LIBRARY_STAGES=(env potential model scores)

log()  { printf '\033[1;34m[reproduce]\033[0m %s\n' "$*"; }
die()  { printf '\033[1;31m[reproduce] FAILED\033[0m %s\n' "$*" >&2; exit 1; }

usage() {
  cat <<'EOF'
stages:
  env         record the version, the descriptor-catalogue digest and the data root        1 s
  potential   re-derive TCRen2 from Native2026 and assert it reproduces the shipped      24 s
              matrix; emit the derivation's own contact table and count it
  model       refit the frozen score model from the committed hold-out slice and          1 s
              assert it reproduces the shipped arrays within the optimiser's tolerance
  scores      featurise all 374 Native2026 crystals and run the whole score set on       80 s
              them, live, through the two released commands
  contacts    place and score a contact map from SEQUENCES ONLY -- no structure of the     30 s
              complex read -- and score its agreement with the crystal it was held out of
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
  python -m pytest tests/regression/test_shipped_potentials.py -q -p no:cacheprovider \
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
# The frozen score model refitted from the hold-out slice committed to the repository. It matches
# the shipped arrays to a relative 1e-5, which is what the optimiser earns and not a round number:
# each of the 55 Yeo-Johnson lambdas comes from a Brent search that stops at about 1.5e-8, and
# every array fitted through the transform inherits that. The FULL hold-out refit needs all 8,292
# structures featurised, which is hours; that test is opt-in behind TCREN_HOLDOUT_FEATURES, and if
# that variable is set here the file must exist. NOTE: that opt-in test asserts np.array_equal on
# arrays subject to the lambda term above, so it is expected to hold on macOS arm64 and not off it.
# The committed-slice test is the one with a defensible bound.
stage_model() {
  if [[ -n "${TCREN_HOLDOUT_FEATURES:-}" && ! -s "${TCREN_HOLDOUT_FEATURES}" ]]; then
    echo "TCREN_HOLDOUT_FEATURES is set to '$TCREN_HOLDOUT_FEATURES', which is not a file" >&2
    return 1
  fi
  if [[ -z "${TCREN_HOLDOUT_FEATURES:-}" ]]; then
    log "  full hold-out refit NOT RUN: set TCREN_HOLDOUT_FEATURES to a table from"
    log "  \`tcren features\` over the 8,292 structures score.holdout_manifest() names."
    log "  The committed-slice refit below runs either way and is the one under budget."
  fi
  python -m pytest tests/unit/test_score.py -q -p no:cacheprovider -ra \
    -k "refit or reproducible"
}

# --- scores -------------------------------------------------------------------------------------
# The whole read-out, live, through the two released commands rather than from any cached table.
# `features` annotates the set in ONE arda call and ONE mmseqs MHC search — never one per
# structure — so this is an 80-second stage rather than an hour-long one.
stage_scores() {
  tcren features -s "$DATA_DIR/Native2026" \
      -i placement,interface,topology,energetics \
      -t 0 -o "$OUT/native_features.tsv" \
    || return 1
  tcren assess --features "$OUT/native_features.tsv" -o "$OUT/native_assess.tsv" \
    || return 1
  python - "$OUT/native_features.tsv" "$OUT/native_assess.tsv" <<'PY' | tee "$OUT/scores.txt"
import sys
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
    print(f"  {c:<20} {sc[c].null_count():>4} null of {sc.height}")
PY
}

# --- contacts -----------------------------------------------------------------------------------
# The 3.7.x claim end to end through the released commands: a contact map placed from loop and
# peptide sequences with NO structure of the complex read, scored under every reference geometry,
# then compared against the crystal it was held out of. Two invocations, because they say different
# things -- the first proves the claim literally (a TSV of sequences goes in and a scored map comes
# out), the second gives it a number against ground truth.
stage_contacts() {
  local pdb=3qib
  # A leave-one-out prior, so what follows is not in-sample.  --keep-surfaces is not a tuning flag:
  # scoring one map off the MEAN surface reads under half the contacts, because a contacting pair
  # sits at the near end of the distribution being averaged and because per-pair averages need not
  # describe a chain any loop could adopt.
  tcren potts geom -s "$DATA_DIR/Native2026" --exclude "$pdb" --grid 32 --keep-surfaces \
      -o "$OUT/repro_geom.npz" \
    || return 1
  python - "$DATA_DIR/Native2026/$pdb.pdb.gz" "$OUT/repro_loops.tsv" <<'SNIP' || return 1
import sys
from tcren.annotation.batch import iter_annotated_set
from tcren.mhc.regions import annotate_mhc_batch
from tcren.potts.geom import loops_from_structure
from tcren.structure.model import PEPTIDE_TYPE

s = next(iter(iter_annotated_set(sys.argv[1])))
annotate_mhc_batch([s])
loops = loops_from_structure(s)
pep = "".join(r.aa for c in s.chains
              if getattr(c, "chain_type", None) == PEPTIDE_TYPE for r in c.residues if r.aa)
cls = "MHCII" if any(getattr(c, "mhc_class", None) == "MHCII"
                     for c in getattr(s, "mhc_calls", []) or []) else "MHCI"
cols = sorted(loops)
with open(sys.argv[2], "w") as fh:
    fh.write("\t".join(["id", "peptide", "mhc_class", *cols]) + "\n")
    fh.write("\t".join([s.pdb_id, pep, cls, *[loops[c] for c in cols]]) + "\n")
print(f"wrote {sys.argv[2]}: {len(cols)} loops, peptide {pep} ({cls})")
SNIP
  # (1) sequences only.  Nothing below this line opens a structure.
  tcren potts predict --sequences "$OUT/repro_loops.tsv" -m potts_tcr_peptide_sc --ensemble \
      --geom-prior "$OUT/repro_geom.npz" -o "$OUT/repro_pred_seq.tsv" \
    || return 1
  # (2) the same map off the structure, which joins the crystal's own sigma back on so the table
  #     can be scored against what the crystal actually did.
  tcren potts predict -s "$DATA_DIR/Native2026/$pdb.pdb.gz" -m potts_tcr_peptide_sc --ensemble \
      --geom-prior "$OUT/repro_geom.npz" -o "$OUT/repro_pred.tsv" \
    || return 1
  tcren potts contacts -s "$DATA_DIR/Native2026/$pdb.pdb.gz" -m potts_tcr_peptide_sc \
      -o "$OUT/repro_contacts.tsv" \
    || return 1
  tcren potts agreement -i "$OUT/repro_pred.tsv" -o "$OUT/repro_agree_pred.tsv" || return 1
  tcren potts agreement -i "$OUT/repro_contacts.tsv" -o "$OUT/repro_agree_cry.tsv" || return 1
  python - "$OUT/repro_pred_seq.tsv" "$OUT/repro_agree_pred.tsv" "$OUT/repro_agree_cry.tsv" \
      <<'SNIP' | tee "$OUT/contacts.txt"
import sys
import polars as pl

seq = pl.read_csv(sys.argv[1], separator="\t", infer_schema_length=None)
pred = pl.read_csv(sys.argv[2], separator="\t").row(0, named=True)
cry = pl.read_csv(sys.argv[3], separator="\t").row(0, named=True)
print(f"from sequences alone: {seq.height} sites, "
      f"expected contacts {seq['p_model'].sum():.2f} -- no structure was read")
hdr = ("geometry", "sites", "contacts", "expected", "ratio", "js", "auc")
print(f"{hdr[0]:<20s}{hdr[1]:>7s}{hdr[2]:>10s}{hdr[3]:>10s}{hdr[4]:>7s}{hdr[5]:>7s}{hdr[6]:>7s}")
for tag, r in (("the crystal's own", cry), ("predicted", pred)):
    print(f"{tag:<20s}{r['sites']:7d}{r['contacts']:10d}{r['expected']:10.2f}"
          f"{r['ratio']:7.2f}{r['js']:7.3f}{r['auc']:7.3f}")
if seq.height == 0 or seq["p_model"].sum() <= 0:
    sys.exit("the sequence-only prediction scored nothing")
# js is scale-free in the level, so the crystal's own geometry is what the predicted one is read
# against -- not 1, which would need a deterministic map.
if not cry["js"] > pred["js"]:
    sys.exit(f"predicted geometry read js={pred['js']:.3f} against the crystal's "
             f"{cry['js']:.3f}; the two model rows of Figure 2c are not ordered as reported")
SNIP
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
log "every stage green"
