#!/bin/bash
#SBATCH --job-name=tcren-fold-bench
#SBATCH --partition=medium
#SBATCH --constraint=compute
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=fold_bench_%j.out
#SBATCH --error=fold_bench_%j.err
#
# Full-scale peptide-recovery benchmark: the native tcren engines against their reference oracles
# (OpenMM, ProMod3, PyRosetta FlexPepDock) over the whole Native2026 set.
#
# Run it locally at --limit 6 first; this exists because FlexPepDock is minutes per structure, so
# n = 374 is a cluster job rather than a laptop one.
#
#   aldan3 slurm submit scripts/fold_benchmark.sbatch --env tcren-fold --partition medium
#   aldan3 slurm queue && aldan3 slurm log <jobid>
#   aldan3 pull /projects/<...>/fold_benchmark.csv scratch/
#
# --cpus-per-task=16 is for the per-structure parallelism, not for one engine: OpenMM's
# LocalEnergyMinimizer and FlexPepDock are both single-threaded here. Right-size it to the actual
# structure count before scaling up (CLAUDE.md §10).
set -euo pipefail

REPO="${TCREN_REPO:-$HOME/tcren}"
DATASET="${DATASET:-$REPO/data/Native2026}"
OUT="${OUT:-$PWD/fold_benchmark.csv}"

cd "$REPO"
echo "host=$(hostname) cpus=${SLURM_CPUS_PER_TASK:-?} dataset=$DATASET"
python -c "import tcren, sys; print('tcren', tcren.__version__, sys.version)"
python -c "from tcren.refine.engines import available_engines; print('engines:', available_engines())"

RUN_BENCHMARK=1 python scripts/fold_benchmark.py \
    --dataset "$DATASET" \
    --oracle \
    --out "$OUT" \
    "$@"

echo "wrote $OUT"
