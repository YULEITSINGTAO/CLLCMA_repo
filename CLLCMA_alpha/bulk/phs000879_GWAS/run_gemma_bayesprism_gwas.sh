#!/usr/bin/env bash
# GEMMA LMM GWAS: BayesPrism cell fractions + covariates (sex, age/10, IGHV)
#
# Usage:
#   nohup bash run_gemma_bayesprism_gwas.sh > gemma_bayesprism_covar/nohup_gwas.log 2>&1 &
#
# Pipeline:
#   1) prepare_bayesprism_gwas_inputs.py  (BayesPrism INT phenos + covars)
#   2) PLINK --keep analysis-ready samples
#   3) GEMMA -gk then -lmm with covariates: intercept + sex + age10 + IGHV

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

CONDA_BASE="${CONDA_BASE:-${HOME}/miniconda3}"
# shellcheck disable=SC1091
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate gwas

command -v gemma >/dev/null
command -v plink >/dev/null
command -v python >/dev/null

BFILE_IN="phs000879"
INPUT_DIR="bayesprism_gwas_inputs"
OUTDIR="gemma_bayesprism_covar"
PREFIX="${OUTDIR}/phs000879_keep"

MAX_SAMPLE_MISS="0.95"
FILTER_SNPS=0

GK_MISS="0.05"
GK_MAF="0.05"
LMM_MISS="1"
LMM_MAF="0"

GK_TYPE=1
LMM_TYPE=1

export OUTDIR PREFIX INPUT_DIR

mkdir -p "${OUTDIR}/output" "${OUTDIR}/logs"

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }

log "=== GEMMA BayesPrism GWAS start (sex + age10 + IGHV covariates) ==="
log "cwd=${SCRIPT_DIR}"
log "gemma=$(command -v gemma); $(gemma -h 2>&1 | head -1 || true)"
log "plink=$(plink --version 2>&1 | head -1)"

# ---------------------------------------------------------------------------
# 0) Prepare phenotypes / covariates
# ---------------------------------------------------------------------------
log "Preparing BayesPrism phenotypes and covariates"
python prepare_bayesprism_gwas_inputs.py \
  > "${OUTDIR}/logs/prepare_inputs.log" 2>&1
cp -f "${INPUT_DIR}/phenotype_legend.txt" "${OUTDIR}/phenotype_legend.txt"
cp -f "${INPUT_DIR}/sample_qc_summary.tsv" "${OUTDIR}/sample_qc_summary.tsv"
cp -f "${INPUT_DIR}/phenotype_variance.tsv" "${OUTDIR}/phenotype_variance.tsv"
cp -f "${INPUT_DIR}/phenotypes_to_run.txt" "${OUTDIR}/phenotypes_to_run.txt"

mapfile -t PHENOS < "${INPUT_DIR}/phenotypes_to_run.txt"
if ((${#PHENOS[@]} == 0)); then
  echo "ERROR: no phenotypes to run" >&2
  exit 1
fi
log "Phenotypes (${#PHENOS[@]}): ${PHENOS[*]}"

# ---------------------------------------------------------------------------
# 1) Keep analysis-ready samples (BayesPrism + complete covariates)
# ---------------------------------------------------------------------------
KEEP="${OUTDIR}/keep_analysis_ready.txt"
awk -F'\t' 'NR>1 && $5=="True" {print $1, $2}' \
  "${INPUT_DIR}/sample_qc_summary.tsv" > "${KEEP}"
# pandas may write True/False; also accept TRUE/true/1
if [[ ! -s "${KEEP}" ]]; then
  awk -F'\t' 'NR>1 && ($5=="TRUE" || $5=="true" || $5=="1") {print $1, $2}' \
    "${INPUT_DIR}/sample_qc_summary.tsv" > "${KEEP}"
fi

# Prefer python for robust boolean parsing
python - <<'PY'
from pathlib import Path
import os
import pandas as pd

outdir = Path(os.environ["OUTDIR"])
inp = Path(os.environ["INPUT_DIR"])
qc = pd.read_csv(inp / "sample_qc_summary.tsv", sep="\t")
ready = qc[qc["analysis_ready"].astype(str).str.lower().isin(["true", "1"])]
keep = outdir / "keep_analysis_ready.txt"
keep.write_text(
    "\n".join(f"{r.FID} {r.IID}" for r in ready.itertuples()) + ("\n" if len(ready) else "")
)
print(f"Keep list: {len(ready)} samples -> {keep}")
PY

NSAMP=$(wc -l < "${KEEP}")
log "Analysis-ready samples: ${NSAMP}"
if [[ "${NSAMP}" -lt 5 ]]; then
  echo "ERROR: too few analysis-ready samples (${NSAMP})" >&2
  exit 1
fi

# Optional extra filter: drop near-empty genotype samples
MISS_PREFIX="${OUTDIR}/missingness"
plink \
  --bfile "${BFILE_IN}" \
  --keep "${KEEP}" \
  --missing \
  --out "${MISS_PREFIX}" \
  > "${OUTDIR}/logs/plink_missing.log" 2>&1

awk -v maxm="${MAX_SAMPLE_MISS}" 'NR>1 && ($6+0) <= (maxm+0) {print $1, $2}' \
  "${MISS_PREFIX}.imiss" > "${OUTDIR}/keep_final.txt"
KEEP="${OUTDIR}/keep_final.txt"
NSAMP=$(wc -l < "${KEEP}")
log "After genotype missingness filter (F_MISS<=${MAX_SAMPLE_MISS}): ${NSAMP}"

# ---------------------------------------------------------------------------
# 2) PLINK subset
# ---------------------------------------------------------------------------
log "PLINK: --keep only (NO SNP filters)"
plink \
  --bfile "${BFILE_IN}" \
  --keep "${KEEP}" \
  --make-bed \
  --out "${PREFIX}" \
  > "${OUTDIR}/logs/plink_qc.log" 2>&1

NSNP=$(wc -l < "${PREFIX}.bim")
NIND=$(wc -l < "${PREFIX}.fam")
log "Genotypes for GEMMA: ${NIND} individuals, ${NSNP} SNPs"

# ---------------------------------------------------------------------------
# 3) Restore BayesPrism INT phenotypes + write covariates in fam order
# ---------------------------------------------------------------------------
log "Writing multi-phenotype .fam and covariate file (INT recomputed on final keep)"
python - <<'PY'
from pathlib import Path
import os
import numpy as np
import pandas as pd
from scipy.stats import norm

outdir = Path(os.environ["OUTDIR"])
prefix = Path(os.environ["PREFIX"])
inp = Path(os.environ["INPUT_DIR"])


def blom_int(v: pd.Series) -> pd.Series:
    x = pd.to_numeric(v, errors="coerce")
    ok = x.notna()
    out = pd.Series(np.nan, index=v.index, dtype=float)
    if ok.sum() < 3:
        return out
    r = x[ok].rank(method="average")
    n = float(ok.sum())
    out.loc[ok] = norm.ppf((r - 3.0 / 8.0) / (n + 1.0 / 4.0))
    return out


pheno = pd.read_csv(inp / "phenotypes_raw_and_int.tsv", sep="\t")
features = [c[4:] for c in pheno.columns if c.startswith("raw_")]
qc_fam = pd.read_csv(
    f"{prefix}.fam",
    sep=r"\s+",
    header=None,
    names=["FID", "IID", "father", "mother", "sex", "pheno"],
)

merged = qc_fam[["FID", "IID", "father", "mother"]].merge(
    pheno, on=["FID", "IID"], how="left", validate="one_to_one"
)
assert len(merged) == len(qc_fam)

# Recompute Blom INT on the final kept sample set
for f in features:
    raw_col = f"raw_{f}"
    if raw_col not in merged.columns:
        raise SystemExit(f"Missing raw phenotype column {raw_col}")
    if merged[raw_col].isna().any():
        raise SystemExit(f"Missing raw phenotype {f} after keep")
    merged[f"int_{f}"] = blom_int(merged[raw_col])

if merged[["sex", "age10", "IGHV_mutated"]].isna().any().any():
    raise SystemExit("Missing covariates after keep")

lines = []
for _, r in merged.iterrows():
    parts = [str(r["FID"]), str(r["IID"]), "0", "0", str(int(r["sex"]))]
    for f in features:
        parts.append(f"{float(r[f'int_{f}']):.10g}")
    lines.append(" ".join(parts))
Path(f"{prefix}.fam").write_text("\n".join(lines) + "\n")

cov_lines = []
for _, r in merged.iterrows():
    cov_lines.append(
        f"1 {int(r['sex'])} {float(r['age10']):.10g} {int(r['IGHV_mutated'])}"
    )
cov_path = outdir / "covariates_sex_age_ighv.txt"
cov_path.write_text("\n".join(cov_lines) + "\n")

plink_cov = merged[["FID", "IID", "sex", "age10", "IGHV_mutated"]].copy()
plink_cov.to_csv(outdir / "covariates_sex_age_ighv.tsv", sep="\t", index=False)

run = [
    ln.strip()
    for ln in (outdir / "phenotypes_to_run.txt").read_text().splitlines()
    if ln.strip()
]
map_path = outdir / "pheno_n_map.tsv"
map_rows = ["name\tgemma_n\tfam_col\trun"]
for i, f in enumerate(features, start=1):
    map_rows.append(f"{f}\t{i}\t{i+5}\t{'yes' if f in run else 'no'}")
map_path.write_text("\n".join(map_rows) + "\n")

print(f"Restored {len(features)} phenotypes for {len(merged)} samples (INT on final keep)")
print(f"Covariates -> {cov_path}")
print("Feature order:", ", ".join(features))
PY

COV="${OUTDIR}/covariates_sex_age_ighv.txt"
log "Covariates (intercept + sex + age10 + IGHV): ${COV}"

# ---------------------------------------------------------------------------
# 4) Relatedness matrix
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

# ---------------------------------------------------------------------------
# 5) Association per phenotype
# ---------------------------------------------------------------------------
FAILED=()
while IFS=$'\t' read -r name gemma_n fam_col runflag; do
  [[ "${name}" == "name" || -z "${name}" ]] && continue
  [[ "${runflag}" != "yes" ]] && continue
  out="assoc_${name}"
  log "GEMMA LMM -n ${gemma_n} (${name})"
  set +e
  gemma \
    -bfile "${PREFIX}" \
    -k "${KFILE}" \
    -c "${COV}" \
    -lmm "${LMM_TYPE}" \
    -n "${gemma_n}" \
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
  n_anal=$(grep -E 'number of analyzed SNPs' "${OUTDIR}/output/${out}.log.txt" | tail -1 || true)
  log "  -> ${OUTDIR}/output/${out}.assoc.txt (${n_anal})"
done < "${OUTDIR}/pheno_n_map.tsv"

# ---------------------------------------------------------------------------
# 6) Top hits summary
# ---------------------------------------------------------------------------
python - <<'PY'
from pathlib import Path
import os
import pandas as pd

outdir = Path(os.environ["OUTDIR"])
run = [
    ln.strip()
    for ln in (outdir / "phenotypes_to_run.txt").read_text().splitlines()
    if ln.strip()
]
rows = []
for name in run:
    f = outdir / "output" / f"assoc_{name}.assoc.txt"
    if not f.exists():
        continue
    df = pd.read_csv(f, sep="\t")
    if "p_wald" not in df.columns:
        continue
    top = df.nsmallest(10, "p_wald").copy()
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
  log "=== GEMMA BayesPrism GWAS done ==="
fi
log "Results under: ${SCRIPT_DIR}/${OUTDIR}/"
log "Logs under:    ${SCRIPT_DIR}/${OUTDIR}/logs/"
log "Nohup log:     ${SCRIPT_DIR}/${OUTDIR}/nohup_gwas.log"
