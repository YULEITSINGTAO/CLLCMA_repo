#!/usr/bin/env bash
# Re-run microenvironment associations for each deconvolution method.
set -euo pipefail
cd "$(dirname "$0")/../.."
export CLL_TME_DIR="${CLL_TME_DIR:-$PWD}"

RSCRIPT="${RSCRIPT:-$HOME/miniconda3/envs/deconvolution/bin/Rscript}"
LOG_DIR="CLLCMA_alpha/results/deconvolution/cll_broad_2022/association_rerun_logs"
mkdir -p "$LOG_DIR"

METHODS=(cibersort bayesprism mcp_counter quantiseq xcell)

run_one() {
  local method="$1"
  local log="$LOG_DIR/${method}.log"
  {
    echo "==== $(date -Is) START $method ===="
    export DECONV_METHOD="$method"

    echo "[1/2] mutation_microenv_association -> alteration_microenv_association/"
    "$RSCRIPT" CLLCMA_alpha/scripts/mutation_microenv_association.R

    echo "[2/2] genomic_microenv_association -> alteration_microenv_association/"
    "$RSCRIPT" CLLCMA_alpha/scripts/genomic_microenv_association.R

    echo "==== $(date -Is) DONE $method ===="
  } >"$log" 2>&1
}

echo "Logs -> $LOG_DIR"
for m in "${METHODS[@]}"; do
  echo "Running $m (sequential to limit memory) ..."
  run_one "$m"
done
echo "All methods finished."
ls -la CLLCMA_alpha/results/deconvolution/cll_broad_2022/alteration_microenv_association/*/ 2>/dev/null | head -80
