#!/usr/bin/env bash
# Reproducible GEMMA GWAS for CIBERSORT INT phenotypes (phs000879)
#
# This run:
#   - PLINK: sample keep only (no SNP MAF/GENO filters)
#   - GEMMA kinship (-gk): moderate SNP QC for a stable K
#   - GEMMA LMM: relaxed -miss/-maf so chr3–22 SNPs are analyzed
#     (GEMMA defaults are -miss 0.05 -maf 0.01, which wiped later chromosomes)
#
# Usage:
#   nohup bash run_gemma_cibersort_gwas.sh > gemma_cibersort_relaxed/nohup_gwas.log 2>&1 &
#
# Inputs (same directory as this script):
#   phs000879.{bed,bim,fam}
#
# Outputs:
#   gemma_cibersort_relaxed/

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------
CONDA_BASE="${CONDA_BASE:-${HOME}/miniconda3}"
# shellcheck disable=SC1091
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate gwas

command -v gemma >/dev/null
command -v plink >/dev/null
command -v python >/dev/null

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
BFILE_IN="phs000879"
OUTDIR="gemma_cibersort_relaxed"
PREFIX="${OUTDIR}/phs000879_keep"

MAX_SAMPLE_MISS="0.95"   # drop near-empty genotype samples (JB-0008)
FILTER_SNPS=0            # PLINK SNP filters off

# GEMMA SNP QC
# Kinship: keep well-typed SNPs so K is numerically stable
GK_MISS="0.05"
GK_MAF="0.05"
# Association: analyze nearly all SNPs (overrides GEMMA defaults 0.05 / 0.01)
LMM_MISS="1"
LMM_MAF="0"

GK_TYPE=1
LMM_TYPE=1
# CLL=1 CD4_T=2 CD8_T=3 DC=4 ILC=5 Monocyte=6 NK=7 Progenitor=8 RS=9
PHENOS=(CLL CD4_T CD8_T DC Monocyte NK Progenitor RS)
PHENOS_N=(1 2 3 4 6 7 8 9)

export OUTDIR PREFIX

mkdir -p "${OUTDIR}/output" "${OUTDIR}/logs"

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }

log "=== GEMMA CIBERSORT GWAS start (relaxed LMM SNP QC) ==="
log "cwd=${SCRIPT_DIR}"
log "gemma=$(command -v gemma); $(gemma -h 2>&1 | head -1 || true)"
log "plink=$(plink --version 2>&1 | head -1)"
log "params: MAX_SAMPLE_MISS=${MAX_SAMPLE_MISS} FILTER_SNPS=${FILTER_SNPS}"
log "  GK:  -miss ${GK_MISS} -maf ${GK_MAF} -gk ${GK_TYPE}"
log "  LMM: -miss ${LMM_MISS} -maf ${LMM_MAF} -lmm ${LMM_TYPE}"

# ---------------------------------------------------------------------------
# 1) Sample keep list
# ---------------------------------------------------------------------------
KEEP="${OUTDIR}/keep_pheno.txt"
MISS_PREFIX="${OUTDIR}/missingness"
log "Computing sample missingness for phenotype samples..."
awk '{
  ok=0
  for (i=6; i<=NF; i++) if ($i != "NA" && $i != "nan" && $i != ".") ok=1
  if (ok) print $1, $2
}' "${BFILE_IN}.fam" > "${OUTDIR}/keep_pheno_all.txt"

plink \
  --bfile "${BFILE_IN}" \
  --keep "${OUTDIR}/keep_pheno_all.txt" \
  --missing \
  --out "${MISS_PREFIX}" \
  > "${OUTDIR}/logs/plink_missing.log" 2>&1

awk -v maxm="${MAX_SAMPLE_MISS}" 'NR>1 && ($6+0) <= (maxm+0) {print $1, $2}' \
  "${MISS_PREFIX}.imiss" > "${KEEP}"
NSAMP=$(wc -l < "${KEEP}")
log "Samples kept (phenotype + F_MISS<=${MAX_SAMPLE_MISS}): ${NSAMP}"
log "Excluded for missingness:"
awk -v maxm="${MAX_SAMPLE_MISS}" 'NR>1 && ($6+0) > (maxm+0) {printf "  %s %s F_MISS=%s\n",$1,$2,$6}' \
  "${MISS_PREFIX}.imiss" || true

if [[ "${NSAMP}" -lt 5 ]]; then
  echo "ERROR: too few samples after missingness filter (${NSAMP})" >&2
  exit 1
fi

# ---------------------------------------------------------------------------
# 2) PLINK subset: samples only
# ---------------------------------------------------------------------------
if [[ "${FILTER_SNPS}" -eq 1 ]]; then
  log "PLINK: --keep + SNP filters (unexpected in this run)"
  exit 1
else
  log "PLINK: --keep only (NO SNP filters)"
  plink \
    --bfile "${BFILE_IN}" \
    --keep "${KEEP}" \
    --make-bed \
    --out "${PREFIX}" \
    > "${OUTDIR}/logs/plink_qc.log" 2>&1
fi

NSNP=$(wc -l < "${PREFIX}.bim")
NIND=$(wc -l < "${PREFIX}.fam")
log "Genotypes for GEMMA: ${NIND} individuals, ${NSNP} SNPs -> ${PREFIX}.{bed,bim,fam}"

# ---------------------------------------------------------------------------
# 3) Restore all INT phenotype columns into .fam
# ---------------------------------------------------------------------------
log "Restoring multi-phenotype columns into .fam"
python - <<'PY'
from pathlib import Path
import os
import pandas as pd

outdir = Path(os.environ["OUTDIR"])
prefix = Path(os.environ["PREFIX"])
features = ["CLL", "CD4_T", "CD8_T", "DC", "ILC", "Monocyte", "NK", "Progenitor", "RS"]

orig = pd.read_csv(
    "phs000879.fam",
    sep=r"\s+",
    header=None,
    names=["FID", "IID", "father", "mother", "sex"] + features,
)
qc = pd.read_csv(
    f"{prefix}.fam",
    sep=r"\s+",
    header=None,
    names=["FID", "IID", "father", "mother", "sex", "pheno"],
)

merged = qc[["FID", "IID", "father", "mother", "sex"]].merge(
    orig[["IID"] + features], on="IID", how="left", validate="one_to_one"
)
assert len(merged) == len(qc), "fam phenotype restore length mismatch"
assert merged[features].notna().all(axis=1).all(), "Missing phenotypes after restore"

lines = []
for _, r in merged.iterrows():
    parts = [str(r["FID"]), str(r["IID"]), "0", "0", str(int(r["sex"]))]
    for f in features:
        parts.append(f"{float(r[f]):.10g}")
    lines.append(" ".join(parts))
Path(f"{prefix}.fam").write_text("\n".join(lines) + "\n")

legend = [
    f"{prefix.name}.fam phenotype column legend",
    "Phenotypes: CIBERSORT fractions, rank-based INT (Blom).",
    "ILC present but skipped in association (zero variance).",
    "PLINK SNP filters: none.",
    "GEMMA GK SNP QC: -miss 0.05 -maf 0.05",
    "GEMMA LMM SNP QC: -miss 1 -maf 0",
    "",
    "col\tfeature\tGEMMA_-n\trun",
]
run = {"CLL", "CD4_T", "CD8_T", "DC", "Monocyte", "NK", "Progenitor", "RS"}
for i, f in enumerate(features, start=1):
    legend.append(f"{i+5}\t{f}\t{i}\t{'yes' if f in run else 'no'}")
legend.append(f"\nn_samples\t{len(merged)}")
(outdir / "phenotype_legend.txt").write_text("\n".join(legend) + "\n")
print(f"Restored {len(features)} phenotypes for {len(merged)} samples")
PY

# ---------------------------------------------------------------------------
# 4) Covariates: intercept + sex
# ---------------------------------------------------------------------------
COV="${OUTDIR}/covariates_sex.txt"
awk '{print 1, $5}' "${PREFIX}.fam" > "${COV}"
log "Covariates (intercept + sex): ${COV}"

# ---------------------------------------------------------------------------
# 5) Relatedness matrix (moderate SNP QC)
# ---------------------------------------------------------------------------
log "GEMMA relatedness matrix (-gk ${GK_TYPE} -miss ${GK_MISS} -maf ${GK_MAF})"
gemma \
  -bfile "${PREFIX}" \
  -gk "${GK_TYPE}" \
  -miss "${GK_MISS}" \
  -maf "${GK_MAF}" \
  -outdir "${OUTDIR}/output" \
  -o "phs000879_relatedness" \
  > "${OUTDIR}/logs/gemma_gk.log" 2>&1

KFILE="${OUTDIR}/output/phs000879_relatedness.cXX.txt"
if [[ ! -s "${KFILE}" ]]; then
  echo "ERROR: relatedness matrix not found: ${KFILE}" >&2
  tail -50 "${OUTDIR}/logs/gemma_gk.log" >&2 || true
  exit 1
fi
log "Relatedness matrix: ${KFILE}"
grep -E 'analyzed SNPs|total SNPs' "${OUTDIR}/output/phs000879_relatedness.log.txt" \
  >> "${OUTDIR}/logs/gemma_gk.log" || true

# ---------------------------------------------------------------------------
# 6) Association: LMM Wald per phenotype (relaxed SNP QC)
# ---------------------------------------------------------------------------
FAILED=()
for i in "${!PHENOS[@]}"; do
  name="${PHENOS[$i]}"
  n="${PHENOS_N[$i]}"
  out="assoc_${name}"
  log "GEMMA LMM -n ${n} (${name}) -miss ${LMM_MISS} -maf ${LMM_MAF}"
  set +e
  gemma \
    -bfile "${PREFIX}" \
    -k "${KFILE}" \
    -c "${COV}" \
    -lmm "${LMM_TYPE}" \
    -n "${n}" \
    -miss "${LMM_MISS}" \
    -maf "${LMM_MAF}" \
    -outdir "${OUTDIR}/output" \
    -o "${out}" \
    > "${OUTDIR}/logs/gemma_lmm_${name}.log" 2>&1
  rc=$?
  set -e
  if [[ "${rc}" -ne 0 || ! -s "${OUTDIR}/output/${out}.assoc.txt" ]]; then
    log "WARNING: association failed for ${name} (rc=${rc}); continuing"
    FAILED+=("${name}")
    tail -20 "${OUTDIR}/logs/gemma_lmm_${name}.log" || true
    continue
  fi
  # log how many SNPs analyzed
  n_anal=$(grep -E 'number of analyzed SNPs' "${OUTDIR}/output/${out}.log.txt" | tail -1 || true)
  log "  -> ${OUTDIR}/output/${out}.assoc.txt (${n_anal})"
done

# ---------------------------------------------------------------------------
# 7) Summary of top hits
# ---------------------------------------------------------------------------
python - <<'PY'
from pathlib import Path
import os
import pandas as pd

outdir = Path(os.environ["OUTDIR"]) / "output"
phenos = ["CLL", "CD4_T", "CD8_T", "DC", "Monocyte", "NK", "Progenitor", "RS"]
rows = []
for name in phenos:
    f = outdir / f"assoc_{name}.assoc.txt"
    if not f.exists():
        continue
    df = pd.read_csv(f, sep="\t")
    if "p_wald" not in df.columns:
        continue
    top = df.nsmallest(10, "p_wald").copy()
    top.insert(0, "feature", name)
    rows.append(top)

summary = Path(os.environ["OUTDIR"]) / "top_hits_summary.tsv"
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
  log "=== GEMMA CIBERSORT GWAS done (relaxed LMM SNP QC) ==="
fi
log "Results under: ${SCRIPT_DIR}/${OUTDIR}/"
log "Logs under:    ${SCRIPT_DIR}/${OUTDIR}/logs/"
log "Nohup log:     ${SCRIPT_DIR}/${OUTDIR}/nohup_gwas.log"
