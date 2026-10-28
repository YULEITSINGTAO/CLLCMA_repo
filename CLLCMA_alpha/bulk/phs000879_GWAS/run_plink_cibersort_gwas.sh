#!/usr/bin/env bash
# PLINK linear GWAS for CIBERSORT INT phenotypes (phs000879)
# Sensitivity / comparison run vs GEMMA LMM.
#
# Usage:
#   nohup bash run_plink_cibersort_gwas.sh > plink_cibersort_linear/nohup_gwas.log 2>&1 &
#
# Uses the same sample-kept genotypes as the relaxed GEMMA run:
#   gemma_cibersort_relaxed/phs000879_keep.{bed,bim,fam}

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

CONDA_BASE="${CONDA_BASE:-${HOME}/miniconda3}"
# shellcheck disable=SC1091
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate gwas

command -v plink >/dev/null
command -v python >/dev/null

BFILE_IN="gemma_cibersort_relaxed/phs000879_keep"
OUTDIR="plink_cibersort_linear"
# Same traits as GEMMA (skip ILC: zero variance). Include DC for completeness.
PHENOS=(CLL CD4_T CD8_T DC Monocyte NK Progenitor RS)

# Optional PLINK SNP filters (0 = none, to mirror relaxed GEMMA LMM)
MAF="0"
GENO="1"

mkdir -p "${OUTDIR}/output" "${OUTDIR}/logs"
export OUTDIR

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }

log "=== PLINK linear CIBERSORT GWAS start ==="
log "cwd=${SCRIPT_DIR}"
log "plink=$(plink --version 2>&1 | head -1)"
log "bfile=${BFILE_IN}"
log "params: MAF=${MAF} GENO=${GENO}"

if [[ ! -s "${BFILE_IN}.bed" ]]; then
  echo "ERROR: missing ${BFILE_IN}.bed — run GEMMA relaxed prep first" >&2
  exit 1
fi

# ---------------------------------------------------------------------------
# 1) Phenotype + covariate files for PLINK
# ---------------------------------------------------------------------------
log "Writing phenotype and covariate files from multi-pheno .fam"
python - <<'PY'
from pathlib import Path
import os
import pandas as pd

outdir = Path(os.environ["OUTDIR"])
features = ["CLL", "CD4_T", "CD8_T", "DC", "ILC", "Monocyte", "NK", "Progenitor", "RS"]
fam = pd.read_csv(
    "gemma_cibersort_relaxed/phs000879_keep.fam",
    sep=r"\s+",
    header=None,
    names=["FID", "IID", "father", "mother", "sex"] + features,
)

# PLINK phenotype file: FID IID + traits; missing as -9
pheno = fam[["FID", "IID"] + features].copy()
for f in features:
    pheno[f] = pd.to_numeric(pheno[f], errors="coerce")
pheno = pheno.fillna(-9)
pheno_path = outdir / "phenotypes.txt"
pheno.to_csv(pheno_path, sep="\t", index=False, header=False)
# also keep a headered copy for humans
pheno.to_csv(outdir / "phenotypes.tsv", sep="\t", index=False)

# Covariates: FID IID sex (1/2)
cov = fam[["FID", "IID", "sex"]].copy()
cov_path = outdir / "covariates_sex.txt"
cov.to_csv(cov_path, sep="\t", index=False, header=False)

# phenotype column index map (1-based within pheno file after FID IID)
legend = ["FID IID then phenotype columns (PLINK --mpheno):", ""]
for i, f in enumerate(features, start=1):
    run = "yes" if f != "ILC" else "no"
    legend.append(f"{i}\t{f}\t{run}")
legend.append(f"\nn_samples\t{len(fam)}")
(outdir / "phenotype_legend.txt").write_text("\n".join(legend) + "\n")

print(f"Wrote {pheno_path} and {cov_path} for {len(fam)} samples")
PY

# Map phenotype name -> --mpheno index (1-based among phenotype columns)
declare -A MPHENO=(
  [CLL]=1 [CD4_T]=2 [CD8_T]=3 [DC]=4 [ILC]=5
  [Monocyte]=6 [NK]=7 [Progenitor]=8 [RS]=9
)

# ---------------------------------------------------------------------------
# 2) PLINK linear association per phenotype
# ---------------------------------------------------------------------------
FILTER_ARGS=()
if [[ "${MAF}" != "0" ]]; then
  FILTER_ARGS+=(--maf "${MAF}")
fi
if [[ "${GENO}" != "1" ]]; then
  FILTER_ARGS+=(--geno "${GENO}")
fi

FAILED=()
for name in "${PHENOS[@]}"; do
  n="${MPHENO[$name]}"
  out="${OUTDIR}/output/assoc_${name}"
  log "PLINK --linear --mpheno ${n} (${name})"
  set +e
  plink \
    --bfile "${BFILE_IN}" \
    --pheno "${OUTDIR}/phenotypes.txt" \
    --mpheno "${n}" \
    --covar "${OUTDIR}/covariates_sex.txt" \
    --linear hide-covar \
    --allow-no-sex \
    --ci 0.95 \
    "${FILTER_ARGS[@]}" \
    --out "${out}" \
    > "${OUTDIR}/logs/plink_linear_${name}.log" 2>&1
  rc=$?
  set -e
  if [[ "${rc}" -ne 0 || ! -s "${out}.assoc.linear" ]]; then
    log "WARNING: PLINK failed for ${name} (rc=${rc}); continuing"
    FAILED+=("${name}")
    tail -30 "${OUTDIR}/logs/plink_linear_${name}.log" || true
    continue
  fi
  nsnp=$(awk 'END{print NR-1}' "${out}.assoc.linear")
  log "  -> ${out}.assoc.linear (${nsnp} SNPs)"
done

# ---------------------------------------------------------------------------
# 3) Top hits summary (by P)
# ---------------------------------------------------------------------------
python - <<'PY'
from pathlib import Path
import os
import pandas as pd

outdir = Path(os.environ["OUTDIR"])
phenos = ["CLL", "CD4_T", "CD8_T", "DC", "Monocyte", "NK", "Progenitor", "RS"]
rows = []
for name in phenos:
    f = outdir / "output" / f"assoc_{name}.assoc.linear"
    if not f.exists():
        continue
    df = pd.read_csv(f, sep=r"\s+")
    # PLINK additive test rows only if TEST column present
    if "TEST" in df.columns:
        df = df[df["TEST"] == "ADD"]
    if "P" not in df.columns:
        continue
    df = df.dropna(subset=["P"])
    top = df.nsmallest(10, "P").copy()
    top.insert(0, "feature", name)
    rows.append(top)

summary = outdir / "top_hits_summary.tsv"
if rows:
    out = pd.concat(rows, ignore_index=True)
    out.to_csv(summary, sep="\t", index=False)
    print(f"Wrote {summary} ({len(out)} rows)")
else:
    print("No association files to summarize")
PY

if ((${#FAILED[@]})); then
  log "Completed with failures: ${FAILED[*]}"
else
  log "=== PLINK linear CIBERSORT GWAS done ==="
fi
log "Results under: ${SCRIPT_DIR}/${OUTDIR}/"
log "Nohup log:     ${SCRIPT_DIR}/${OUTDIR}/nohup_gwas.log"
