#!/usr/bin/env bash
# tcren bootstrap — the environment and the data. 2026-09-05
#
# THIS SCRIPT IS OUTSIDE ANY TIME BUDGET, deliberately. It downloads about 620 MB of structures
# from Hugging Face and builds the ~28k-allele IMGT/HLA MHC reference plus its mmseqs prefilter
# index. How long that takes is a property of the network and the machine, not of this project,
# so it is not something `reproduce.sh` can be asked to fit inside 30 minutes. Run this once;
# run `./reproduce.sh` as often as you like afterwards.
#
#   ./bootstrap.sh                    venv, editable install, every structure set, MHC reference
#   ./bootstrap.sh --sets Native2026,Canonical2026,tcrvdb
#                                     only these sets (see `tcren fetch-data --list`)
#   ./bootstrap.sh --tests            also run the fast test suite at the end
#   ./bootstrap.sh --tables           also fetch the derived-table deposit (SEE THE NOTE BELOW)
#
# WHAT IT LEAVES BEHIND, and where:
#   $TCREN_HOME       the on-disk root. $TCREN_HOME, else this checkout, else $XDG_CACHE_HOME/tcren
#   $TCREN_DATA_DIR   the structure sets. Defaults to <root>/data
#   <root>/database/mhc/alleles.aa.fasta   the built MHC allele reference
#
# --tables IS EXPECTED TO FAIL TODAY. It forwards to `tcren fetch-data --tables`, which reads the
# derived-table deposit `isalgo/tcren2_reproduce` — the outputs of the corpus-scale featurisations
# that do not fit in a reproduction budget. That repository does not exist on the Hub yet, so the
# flag is off by default and nothing in `reproduce.sh` depends on it. It is wired up so that the
# day the deposit is published, one flag reaches it.
#
# FAILS LOUDLY. Every step is checked afterwards against the filesystem: a set that did not land,
# an MHC reference that was not written, an import that does not resolve. Nothing is skipped
# silently, and a missing input is an exit code, not a warning.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SETS="all"
DO_TESTS=0
DO_TABLES=0
FORCE_MHC=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sets)    SETS="${2:?--sets needs a comma-separated list, or 'all'}"; shift 2 ;;
    --sets=*)  SETS="${1#*=}"; shift ;;
    --tests)   DO_TESTS=1; shift ;;
    --tables)  DO_TABLES=1; shift ;;
    --force-mhc) FORCE_MHC=1; shift ;;
    -h|--help) awk 'NR>1 { if (/^#/) { sub(/^# ?/, ""); print } else { exit } }' "${BASH_SOURCE[0]}"; exit 0 ;;
    *) echo "bootstrap: unknown flag: $1 (try --help)" >&2; exit 2 ;;
  esac
done

log()  { printf '\033[1;34m[bootstrap]\033[0m %s\n' "$*"; }
die()  { printf '\033[1;31m[bootstrap] FAILED\033[0m %s\n' "$*" >&2; exit 1; }

# --- 1. venv, editable install, and the two orientation reference sets ------------------------
# setup.sh is the source of truth for the install: it checks for uv and a C++ compiler, creates
# .venv, runs `uv pip install -e .` (which builds the C++ extensions and pulls arda-mapper), and
# calls `tcren fetch-data` for Native2026 + Canonical2026. Reimplementing any of that here would
# give the project two install recipes that drift apart.
log "environment and the orientation reference sets (setup.sh)"
bash "$ROOT/setup.sh"

# shellcheck disable=SC1091
source "$ROOT/.venv/bin/activate"

command -v tcren >/dev/null 2>&1 \
  || die "the tcren console script is not on PATH after setup.sh; the install did not complete"

TCREN_ROOT="$(python -c 'from tcren.paths import tcren_home; print(tcren_home())')"
DATA_DIR="$(python -c 'from tcren.paths import data_dir; print(data_dir())')"
MHC_FASTA="$(python -c 'from tcren.mhc.reference import DATABASE_DIR; print(DATABASE_DIR / "alleles.aa.fasta")')"
log "root      $TCREN_ROOT"
log "data dir  $DATA_DIR"
log "mhc ref   $MHC_FASTA"

# --- 2. the structure sets --------------------------------------------------------------------
# The listing is read from the Hub, so `--sets all` reaches whatever is deposited rather than a
# hard-coded tuple that names a third of it. Sets already present are skipped without a network
# round trip; pass --force to `tcren fetch-data` yourself to re-fetch one.
log "structure sets: --sets $SETS"
tcren fetch-data --sets "$SETS" || die "tcren fetch-data --sets $SETS exited non-zero"

if [[ "$DO_TABLES" -eq 1 ]]; then
  log "derived-table deposit (expected to fail until isalgo/tcren2_reproduce is published)"
  tcren fetch-data --tables --sets Native2026 \
    || die "tcren fetch-data --tables failed. The deposit is not published yet; drop --tables."
fi

# --- 3. the MHC allele reference ----------------------------------------------------------------
# Built from IMGT/HLA + the UniProt mouse H-2 records, not downloaded ready-made and not bundled
# in the wheel. Without it every structure fails to annotate, which is the single most common
# reason a fresh checkout produces an empty table.
if [[ "$FORCE_MHC" -eq 1 || ! -s "$MHC_FASTA" ]]; then
  log "building the MHC allele reference (IMGT/HLA + UniProt mouse) — minutes, not seconds"
  tcren build-mhc-ref || die "tcren build-mhc-ref exited non-zero"
else
  log "MHC allele reference already built, leaving it alone (--force-mhc to rebuild)"
fi

# --- 4. verify, loudly ---------------------------------------------------------------------------
log "verifying what landed"

[[ -s "$MHC_FASTA" ]] \
  || die "no MHC reference at $MHC_FASTA. Re-run with --force-mhc."
MHC_N="$(grep -c '^>' "$MHC_FASTA")"
[[ "$MHC_N" -gt 1000 ]] \
  || die "the MHC reference at $MHC_FASTA holds only $MHC_N alleles; expected tens of thousands"
log "  MHC reference: $MHC_N alleles"

# Every set that was asked for must exist on disk and hold at least one non-hidden file. "Holds a
# .pdb" is the wrong check and was a real defect: several sets ship a single .tar.gz and one ships
# a table, so requiring a loose structure file re-downloaded them on every run.
python - "$DATA_DIR" "$SETS" <<'PY' || die "one or more requested structure sets did not land"
import sys
from pathlib import Path

data_dir, sets = Path(sys.argv[1]), sys.argv[2]
from tcren.paper.bootstrap import hf_listing

listing = hf_listing()
want = sorted(listing) if sets.strip().lower() == "all" else [s.strip() for s in sets.split(",") if s.strip()]

bad = []
for name in want:
    d = data_dir / name
    files = [p for p in d.glob("*") if p.is_file() and not p.name.startswith(".")] if d.is_dir() else []
    print(f"  {name:<26} {len(files):>5} files")
    if not files:
        bad.append(name)
for name in bad:
    print(f"MISSING: {data_dir / name} is absent or empty "
          f"(re-run: tcren fetch-data --sets {name} --force)", file=sys.stderr)
sys.exit(1 if bad else 0)
PY

# The two provenance stamps a reader has to record before any number is comparable.
tcren info || die "tcren info exited non-zero"
python -c 'from tcren.provenance import registry_digest; print("catalogue digest:", registry_digest())' \
  || die "the descriptor catalogue would not hash; the install is broken"

# --- 5. tests -------------------------------------------------------------------------------------
if [[ "$DO_TESTS" -eq 1 ]]; then
  log "fast test suite"
  pytest "$ROOT/tests" -m "not slow" -q || die "the fast test suite is red"
fi

log "done. Next: ./reproduce.sh"
