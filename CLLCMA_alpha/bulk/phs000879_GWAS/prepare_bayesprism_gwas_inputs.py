#!/usr/bin/env python3
"""Prepare BayesPrism INT phenotypes + sex/age/IGHV covariates for phs000879 GWAS.

Writes under bulk/phs000879_GWAS/bayesprism_gwas_inputs/:
  - phenotypes_raw_and_int.tsv
  - phs000879_bayesprism.fam   (multi-pheno; NA if incomplete)
  - covariates_complete.tsv    (IID-level sex, age, age10, IGHV)
  - phenotype_legend.txt
  - sample_qc_summary.tsv

Does not modify genotype .bed/.bim. The GEMMA runner subsets to complete cases.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parents[1]  # CLLCMA_alpha
BP_TSV = ROOT / "results/deconvolution/cll_broad_2022/bayesprism_sample_by_feature.tsv"
CLIN_TSV = SCRIPT_DIR / "srr_clinical.txt"
FAM_IN = SCRIPT_DIR / "phs000879.fam"
OUTDIR = SCRIPT_DIR / "bayesprism_gwas_inputs"


def blom_int(x: pd.Series) -> pd.Series:
    """Rank-based inverse-normal transform (Blom)."""
    v = pd.to_numeric(x, errors="coerce")
    ok = v.notna()
    out = pd.Series(np.nan, index=x.index, dtype=float)
    if ok.sum() < 3:
        return out
    r = v[ok].rank(method="average")
    n = float(ok.sum())
    out.loc[ok] = norm.ppf((r - 3.0 / 8.0) / (n + 1.0 / 4.0))
    return out


def encode_sex_plink(sex: object, gender: object) -> float:
    for raw in (sex, gender):
        s = str(raw).strip().lower()
        if s.startswith("m") or s in {"1", "1.0"}:
            return 1.0
        if s.startswith("f") or s in {"2", "2.0"}:
            return 2.0
    return np.nan


def encode_ighv(status: object) -> float:
    s = str(status).strip().lower()
    if s == "mutated":
        return 1.0
    if s == "unmutated":
        return 0.0
    return np.nan


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    fam = pd.read_csv(
        FAM_IN,
        sep=r"\s+",
        header=None,
        usecols=[0, 1, 4],
        names=["FID", "IID", "sex_fam"],
    )
    fam["IID"] = fam["IID"].astype(str)
    fam["FID"] = fam["FID"].astype(str)

    bp = pd.read_csv(BP_TSV, sep="\t")
    bp["SAMPLE_ID"] = bp["SAMPLE_ID"].astype(str).str.strip('"')
    features = [c for c in bp.columns if c != "SAMPLE_ID"]
    bp = bp.set_index("SAMPLE_ID")

    clin = pd.read_csv(CLIN_TSV, sep="\t")
    clin = clin.drop_duplicates(subset=["submitted_subject_id"], keep="first")
    clin["submitted_subject_id"] = clin["submitted_subject_id"].astype(str)
    clin = clin.set_index("submitted_subject_id")

    rows = []
    for _, fr in fam.iterrows():
        iid = fr["IID"]
        row: dict[str, object] = {
            "FID": fr["FID"],
            "IID": iid,
            "sex_fam": fr["sex_fam"],
            "has_bayesprism": iid in bp.index,
        }
        if iid in clin.index:
            c = clin.loc[iid]
            row["sex"] = encode_sex_plink(c.get("sex"), c.get("gender"))
            age = pd.to_numeric(c.get("age_at_reference_sample"), errors="coerce")
            if pd.isna(age):
                age = pd.to_numeric(c.get("age_at_diagnosis"), errors="coerce")
            row["age"] = float(age) if pd.notna(age) else np.nan
            row["IGHV_mut_raw"] = c.get("IGHV_mut")
            row["IGHV_mutated"] = encode_ighv(c.get("IGHV_mut"))
        else:
            row["sex"] = np.nan
            row["age"] = np.nan
            row["IGHV_mut_raw"] = np.nan
            row["IGHV_mutated"] = np.nan

        # Prefer clinical sex; fall back to fam sex if coded 1/2
        if pd.isna(row["sex"]):
            sf = pd.to_numeric(fr["sex_fam"], errors="coerce")
            if sf in (1, 2):
                row["sex"] = float(sf)

        row["age10"] = (row["age"] / 10.0) if pd.notna(row["age"]) else np.nan

        if iid in bp.index:
            for f in features:
                row[f"raw_{f}"] = float(bp.loc[iid, f])
        else:
            for f in features:
                row[f"raw_{f}"] = np.nan
        rows.append(row)

    tab = pd.DataFrame(rows)
    tab["cov_complete"] = (
        tab["sex"].notna() & tab["age"].notna() & tab["IGHV_mutated"].notna()
    )
    tab["analysis_ready"] = tab["has_bayesprism"] & tab["cov_complete"]

    # INT within analysis-ready samples only (Blom)
    ready = tab["analysis_ready"].to_numpy()
    for f in features:
        raw = tab[f"raw_{f}"]
        int_all = pd.Series(np.nan, index=tab.index, dtype=float)
        int_all.loc[ready] = blom_int(raw.loc[ready]).to_numpy()
        tab[f"int_{f}"] = int_all

    # Phenotype variance among analysis-ready
    var_rows = []
    run_features: list[str] = []
    for f in features:
        v = tab.loc[ready, f"raw_{f}"].astype(float)
        sd = float(v.std(ddof=1)) if len(v) > 1 else 0.0
        n_nz = int((v > 1e-12).sum())
        keep = sd > 1e-12 and v.nunique(dropna=True) >= 3
        var_rows.append(
            {
                "CELL_TYPE": f,
                "N_READY": int(ready.sum()),
                "MEAN_RAW": float(v.mean()),
                "SD_RAW": sd,
                "N_NONZERO": n_nz,
                "RUN_GWAS": keep,
            }
        )
        if keep:
            run_features.append(f)
    var_df = pd.DataFrame(var_rows)

    # Multi-pheno FAM for all genotype samples (NA if not analysis-ready)
    fam_lines = []
    for _, r in tab.iterrows():
        sex_i = int(r["sex"]) if pd.notna(r["sex"]) else 0
        parts = [str(r["FID"]), str(r["IID"]), "0", "0", str(sex_i)]
        for f in features:
            val = r[f"int_{f}"]
            parts.append("NA" if pd.isna(val) else f"{float(val):.10g}")
        fam_lines.append(" ".join(parts))
    fam_out = OUTDIR / "phs000879_bayesprism.fam"
    fam_out.write_text("\n".join(fam_lines) + "\n")

    # Covariate table (IID order = genotype fam)
    cov = tab[["FID", "IID", "sex", "age", "age10", "IGHV_mutated", "analysis_ready"]].copy()
    cov.to_csv(OUTDIR / "covariates_complete.tsv", sep="\t", index=False)

    # Full phenotype table
    tab.to_csv(OUTDIR / "phenotypes_raw_and_int.tsv", sep="\t", index=False)
    var_df.to_csv(OUTDIR / "phenotype_variance.tsv", sep="\t", index=False)

    sample_qc = tab[
        [
            "FID",
            "IID",
            "has_bayesprism",
            "cov_complete",
            "analysis_ready",
            "sex",
            "age",
            "IGHV_mut_raw",
            "IGHV_mutated",
        ]
    ]
    sample_qc.to_csv(OUTDIR / "sample_qc_summary.tsv", sep="\t", index=False)

    legend = [
        "BayesPrism GWAS phenotype / covariate preparation",
        f"Source fractions: {BP_TSV}",
        f"Clinical: {CLIN_TSV}",
        "Phenotypes: BayesPrism cell fractions, Blom INT within analysis-ready samples",
        "Covariates for association: intercept + sex(1=male,2=female) + age10 + IGHV_mutated(1/0)",
        "Age: age_at_reference_sample (fallback age_at_diagnosis), divided by 10",
        "Analysis-ready: has BayesPrism fractions AND non-missing sex, age, IGHV",
        "",
        f"n_genotype_samples\t{len(tab)}",
        f"n_with_bayesprism\t{int(tab['has_bayesprism'].sum())}",
        f"n_cov_complete\t{int(tab['cov_complete'].sum())}",
        f"n_analysis_ready\t{int(tab['analysis_ready'].sum())}",
        f"n_phenotypes_total\t{len(features)}",
        f"n_phenotypes_run\t{len(run_features)}",
        "",
        "fam_col\tfeature\tGEMMA_-n\trun",
    ]
    for i, f in enumerate(features, start=1):
        legend.append(
            f"{i + 5}\t{f}\t{i}\t{'yes' if f in run_features else 'no'}"
        )
    legend.append("")
    legend.append("Excluded from analysis-ready (examples):")
    for _, r in tab.loc[~tab["analysis_ready"]].iterrows():
        reasons = []
        if not r["has_bayesprism"]:
            reasons.append("no_BayesPrism")
        if not r["cov_complete"]:
            reasons.append("missing_sex_age_or_IGHV")
        legend.append(f"  {r['IID']}\t{','.join(reasons)}")
    (OUTDIR / "phenotype_legend.txt").write_text("\n".join(legend) + "\n")

    # Machine-readable run list for the shell script
    (OUTDIR / "phenotypes_to_run.txt").write_text("\n".join(run_features) + "\n")

    print(f"Wrote {OUTDIR}")
    print(f"Analysis-ready samples: {int(tab['analysis_ready'].sum())} / {len(tab)}")
    print(f"Phenotypes to run: {len(run_features)} / {len(features)}")
    print("Ready IIDs:", ", ".join(tab.loc[tab["analysis_ready"], "IID"].tolist()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
