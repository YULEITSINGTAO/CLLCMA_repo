#!/usr/bin/env python3
"""Splicing-outlier gene burden and per-gene × TME fraction association.

For each alternative-splicing class (SplAdder simplified matrices):

  1. Flag outlier PSI values per event: two-sided extreme tenth percentile
     (≤P10 or ≥P90 across analysis samples) AND |PSI − median| > 0.10.
  2. Collapse to genes: a gene is an outlier in a sample if any of its events
     is an outlier in that sample.
  3. Sum the gene×sample 0/1 table → outlier-gene count per sample; associate
     counts with TME fractions.
  4. Per-gene association: for each gene passing prevalence filters, test
     outlier (0/1) vs each TME fraction with LM controlling for IGHV, sex,
     and age at sampling / 10 (plus Wilcoxon descriptive p / Cliff's δ).

ID map: sequences_and_participants.xlsx  Run → participant_id
        participant_id must match deconvolution SAMPLE_ID.

Env:
  CLL_TME_DIR, DECONV_METHOD (default bayesprism)
"""

from __future__ import annotations

import os
import re
import sys
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

AS_TYPES = [
    "alt_3prime",
    "alt_5prime",
    "exon_skip",
    "intron_retention",
    "mult_exon_skip",
    "mutex_exons",
]

MIN_N_FINITE = 20
PCT_LO = 10.0
PCT_HI = 90.0
MEDIAN_DELTA = 0.10
CHUNKSIZE = 5000
# Per-gene association: require enough outlier and non-outlier samples
MIN_GENE_OUTLIER = 15
MIN_GENE_WT = 15


def project_root() -> Path:
    env = os.environ.get("CLL_TME_DIR")
    if env:
        return Path(env).expanduser().resolve()
    return Path(__file__).resolve().parents[2]


def load_run_participant_map(xlsx: Path) -> pd.DataFrame:
    """Read Run / participant_id via R readxl (no openpyxl dependency)."""
    import subprocess

    rscript = os.environ.get(
        "RSCRIPT",
        str(Path.home() / "miniconda3/envs/deconvolution/bin/Rscript"),
    )
    cmd = [
        rscript,
        "-e",
        (
            "library(readxl); "
            f'd=read_excel("{xlsx.as_posix()}"); '
            'write.table(d[,c("Run","participant_id")], stdout(), '
            "sep='\\t', row.names=FALSE, quote=FALSE)"
        ),
    ]
    out = subprocess.check_output(cmd, text=True)
    df = pd.read_csv(StringIO(out), sep="\t")
    df["Run"] = df["Run"].astype(str)
    df["participant_id"] = df["participant_id"].astype(str)
    return df.drop_duplicates(subset=["Run"], keep="first")


def parse_srr_from_column(col: str) -> str:
    col = str(col)
    m = re.match(r"(SRR\d+)", col)
    if m:
        return m.group(1)
    return re.sub(r"Aligned.*$", "", col).replace(":psi", "")


def strip_ensembl_version(gene: str) -> str | float:
    if gene is None or (isinstance(gene, float) and np.isnan(gene)):
        return np.nan
    g = str(gene).strip()
    if g in {"", "nan", "None", "NA"}:
        return np.nan
    return re.sub(r"\.\d+$", "", g)


def read_cbioportal_clinical(path: Path) -> pd.DataFrame:
    # Skip "#..." meta lines used by cBioPortal exports.
    return pd.read_csv(path, sep="\t", comment="#")


def encode_ighv_mutated(status: pd.Series) -> pd.Series:
    s = status.astype(str)
    out = pd.Series(np.nan, index=status.index, dtype=float)
    out[s.str.contains("unmut", case=False, na=False)] = 0.0
    mut = s.str.contains(r"^mut|mutated", case=False, na=False) & ~s.str.contains(
        "unmut", case=False, na=False
    )
    out[mut] = 1.0
    return out


def encode_sex_male(sex: pd.Series) -> pd.Series:
    s = sex.astype(str).str.upper().str.strip()
    out = pd.Series(np.nan, index=sex.index, dtype=float)
    out[s == "MALE"] = 1.0
    out[s == "FEMALE"] = 0.0
    return out


def attach_covariates(micro: pd.DataFrame, sample_clin: Path, patient_clin: Path) -> pd.DataFrame:
    samp = read_cbioportal_clinical(sample_clin)[["SAMPLE_ID", "PATIENT_ID"]]
    pat = read_cbioportal_clinical(patient_clin)
    keep = [
        c
        for c in ["PATIENT_ID", "SEX", "AGE", "AGE_SAMPLING", "IGHV_MUTATION_STATUS"]
        if c in pat.columns
    ]
    pat = pat[keep]
    meta = samp.merge(pat, on="PATIENT_ID", how="left")
    meta["IGHV_MUTATED"] = encode_ighv_mutated(meta["IGHV_MUTATION_STATUS"])
    meta["SEX_MALE"] = encode_sex_male(meta["SEX"])
    age = pd.to_numeric(meta.get("AGE_SAMPLING"), errors="coerce")
    age2 = pd.to_numeric(meta.get("AGE"), errors="coerce")
    meta["AGE_YEARS"] = age.fillna(age2)
    return micro.merge(meta[["SAMPLE_ID", "IGHV_MUTATED", "SEX_MALE", "AGE_YEARS"]], on="SAMPLE_ID", how="left")


def flag_outliers_chunk(X: np.ndarray) -> np.ndarray:
    """Return boolean outlier matrix same shape as X (NaN positions False)."""
    ok = np.isfinite(X)
    n_ok = ok.sum(axis=1)
    out = np.zeros_like(X, dtype=bool)
    enough = n_ok >= MIN_N_FINITE
    if not np.any(enough):
        return out

    Xe = X[enough]
    oke = ok[enough]
    # nan-aware stats along samples
    med = np.nanmedian(Xe, axis=1)
    p10 = np.nanpercentile(Xe, PCT_LO, axis=1)
    p90 = np.nanpercentile(Xe, PCT_HI, axis=1)
    delta = np.abs(Xe - med[:, None]) > MEDIAN_DELTA
    extreme = (Xe <= p10[:, None]) | (Xe >= p90[:, None])
    out[enough] = oke & extreme & delta
    return out


def accumulate_gene_outliers(
    event_ids: np.ndarray,
    outlier: np.ndarray,
    event_to_gene: pd.Series,
    gene_flags: dict[str, np.ndarray],
) -> None:
    genes = event_to_gene.reindex(event_ids).to_numpy()
    # Process gene-wise within chunk for speed
    order = np.argsort(np.where(pd.isna(genes), "", genes.astype(str)))
    event_ids = event_ids  # noqa: F841 — keep parallel arrays clear
    genes_sorted = genes[order]
    out_sorted = outlier[order]
    i = 0
    n = len(genes_sorted)
    while i < n:
        g = genes_sorted[i]
        if not isinstance(g, str) or g == "" or g == "nan":
            i += 1
            continue
        j = i + 1
        while j < n and genes_sorted[j] == g:
            j += 1
        block = out_sorted[i:j]
        gene_or = np.any(block, axis=0)
        if gene_or.any():
            if g in gene_flags:
                gene_flags[g] |= gene_or
            else:
                gene_flags[g] = gene_or.copy()
        i = j


def counts_from_gene_flags(gene_flags: dict[str, np.ndarray], n_samp: int) -> np.ndarray:
    counts = np.zeros(n_samp, dtype=np.int32)
    for arr in gene_flags.values():
        counts += arr.astype(np.int32)
    return counts


def process_as_type(
    matrix_path: Path,
    info_path: Path,
    use_cols: list[str],
    sample_ids: list[str],
) -> tuple[np.ndarray, dict[str, np.ndarray], dict]:
    info = pd.read_csv(info_path, usecols=["event_id", "gene_name"])
    info["event_id"] = info["event_id"].astype(str)
    info["gene"] = info["gene_name"].map(strip_ensembl_version)
    event_to_gene = info.set_index("event_id")["gene"]

    gene_flags: dict[str, np.ndarray] = {}
    n_events = 0
    n_events_tested = 0
    n_event_sample_flags = 0

    usecols = ["event_id"] + use_cols
    for chunk in pd.read_csv(matrix_path, usecols=usecols, chunksize=CHUNKSIZE):
        ev = chunk["event_id"].astype(str).to_numpy()
        X = chunk[use_cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=np.float32)
        out = flag_outliers_chunk(X)
        n_events += X.shape[0]
        tested = np.isfinite(X).sum(axis=1) >= MIN_N_FINITE
        n_events_tested += int(tested.sum())
        n_event_sample_flags += int(out.sum())
        accumulate_gene_outliers(ev, out, event_to_gene, gene_flags)
        if n_events % 50000 == 0:
            print(
                f"    … {n_events:,} events | outlier genes so far {len(gene_flags):,}",
                flush=True,
            )

    counts = counts_from_gene_flags(gene_flags, len(sample_ids))
    meta = {
        "n_events": n_events,
        "n_events_tested": n_events_tested,
        "n_event_sample_outlier_flags": n_event_sample_flags,
        "n_outlier_genes_any_sample": len(gene_flags),
        "count_min": int(counts.min()) if len(counts) else 0,
        "count_median": float(np.median(counts)) if len(counts) else np.nan,
        "count_mean": float(counts.mean()) if len(counts) else np.nan,
        "count_max": int(counts.max()) if len(counts) else 0,
    }
    return counts, gene_flags, meta


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    """Cliff's δ (a vs b); positive => values in a tend larger."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) == 0 or len(b) == 0:
        return np.nan
    u, _ = stats.mannwhitneyu(a, b, alternative="two-sided")
    return float(2.0 * u / (len(a) * len(b)) - 1.0)


def write_gene_binary(
    path: Path,
    gene_flags: dict[str, np.ndarray],
    sample_ids: list[str],
) -> None:
    genes = sorted(gene_flags.keys())
    if not genes:
        return
    mat = np.vstack([gene_flags[g].astype(np.uint8) for g in genes])
    df = pd.DataFrame(mat, index=genes, columns=sample_ids)
    df.index.name = "GENE"
    df.to_csv(path, sep="\t", compression="gzip")


def per_gene_tme_association(
    gene_flags: dict[str, np.ndarray],
    sample_ids: list[str],
    micro: pd.DataFrame,
    feature_cols: list[str],
    as_type: str,
    min_outlier: int = MIN_GENE_OUTLIER,
    min_wt: int = MIN_GENE_WT,
) -> pd.DataFrame:
    """Associate each gene's outlier status (0/1) with each TME fraction.

    Primary model: fraction ~ outlier + IGHV + sex + age10
    Also report Wilcoxon p and Cliff's δ (outlier vs non-outlier).
    """
    m = micro.set_index("SAMPLE_ID").reindex(sample_ids)
    cov_ok = (
        m["IGHV_MUTATED"].notna()
        & m["SEX_MALE"].notna()
        & m["AGE_YEARS"].notna()
        & m[feature_cols].notna().all(axis=1)
    )
    keep = cov_ok.to_numpy()
    m = m.loc[keep]
    keep_idx = np.flatnonzero(keep)
    if len(m) < 40:
        return pd.DataFrame()

    Y = m[feature_cols].to_numpy(dtype=float)
    ighv = m["IGHV_MUTATED"].to_numpy(dtype=float)
    sex = m["SEX_MALE"].to_numpy(dtype=float)
    age10 = m["AGE_YEARS"].to_numpy(dtype=float) / 10.0
    n = len(m)
    k = 5
    df_resid = n - k

    rows: list[dict] = []
    genes_tested = 0
    for gene, flags_full in gene_flags.items():
        alt = flags_full[keep_idx].astype(float)
        n_out = int(alt.sum())
        n_wt = int(n - n_out)
        if n_out < min_outlier or n_wt < min_wt:
            continue
        genes_tested += 1
        design = np.column_stack([np.ones(n), alt, ighv, sex, age10])
        try:
            beta, *_ = np.linalg.lstsq(design, Y, rcond=None)
            xtx_inv = np.linalg.inv(design.T @ design)
        except np.linalg.LinAlgError:
            continue
        resid = Y - design @ beta
        sigma2 = np.sum(resid**2, axis=0) / df_resid
        se_alt = np.sqrt(np.maximum(sigma2 * xtx_inv[1, 1], 0.0))
        t_alt = np.divide(
            beta[1],
            se_alt,
            out=np.full_like(beta[1], np.nan),
            where=se_alt > 0,
        )
        p_adj = 2.0 * stats.t.sf(np.abs(t_alt), df_resid)

        alt_mask = alt == 1.0
        wt_mask = ~alt_mask
        n_a = int(alt_mask.sum())
        n_w = int(wt_mask.sum())
        for j, ct in enumerate(feature_cols):
            y = Y[:, j]
            a = y[alt_mask]
            w = y[wt_mask]
            try:
                u_res = stats.mannwhitneyu(a, w, alternative="two-sided")
                p_w = float(u_res.pvalue)
                delta = float(2.0 * u_res.statistic / (n_a * n_w) - 1.0)
            except ValueError:
                p_w = np.nan
                delta = np.nan
            rows.append(
                {
                    "AS_TYPE": as_type,
                    "GENE": gene,
                    "CELL_TYPE": ct,
                    "N_OUTLIER": n_out,
                    "N_WT": n_wt,
                    "N_ADJ": n,
                    "MEDIAN_OUTLIER": float(np.median(a)),
                    "MEDIAN_WT": float(np.median(w)),
                    "CLIFFS_DELTA": delta,
                    "P_WILCOX": p_w,
                    "BETA_ADJ": float(beta[1, j]),
                    "SE_ADJ": float(se_alt[j]),
                    "P_ADJ": float(p_adj[j]),
                }
            )
        if genes_tested % 2000 == 0:
            print(
                f"    per-gene {as_type}: tested {genes_tested:,} genes …",
                flush=True,
            )

    print(
        f"  per-gene {as_type}: {genes_tested:,} genes × {len(feature_cols)} cell types",
        flush=True,
    )
    return pd.DataFrame(rows)


def bh_fdr(p: np.ndarray) -> np.ndarray:
    """Benjamini–Hochberg FDR; NaNs preserved."""
    p = np.asarray(p, dtype=float)
    q = np.full_like(p, np.nan, dtype=float)
    ok = np.isfinite(p) & (p >= 0) & (p <= 1)
    if not np.any(ok):
        return q
    pv = p[ok]
    n = len(pv)
    order = np.argsort(pv)
    ranked = pv[order]
    ranks = np.arange(1, n + 1)
    qvals = ranked * n / ranks
    qvals = np.minimum.accumulate(qvals[::-1])[::-1]
    qvals = np.clip(qvals, 0, 1)
    out = np.empty(n, dtype=float)
    out[order] = qvals
    q[ok] = out
    return q


def ols_count_coef(y: np.ndarray, design: np.ndarray) -> tuple[float, float, float]:
    """Return (beta, se, p) for design column 1 (col 0 = intercept)."""
    n, k = design.shape
    if n <= k:
        return np.nan, np.nan, np.nan
    beta, *_ = np.linalg.lstsq(design, y, rcond=None)
    resid = y - design @ beta
    df = n - k
    if df <= 0:
        return np.nan, np.nan, np.nan
    sigma2 = float(np.sum(resid**2) / df)
    try:
        xtx_inv = np.linalg.inv(design.T @ design)
    except np.linalg.LinAlgError:
        return np.nan, np.nan, np.nan
    se = float(np.sqrt(max(sigma2 * xtx_inv[1, 1], 0.0)))
    if not np.isfinite(se) or se == 0:
        return float(beta[1]), se, np.nan
    t = float(beta[1] / se)
    p = float(2 * stats.t.sf(abs(t), df))
    return float(beta[1]), se, p


def spearman_assoc(counts: pd.Series, fractions: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for ct in fractions.columns:
        y = fractions[ct].astype(float)
        ok = counts.notna() & y.notna()
        if ok.sum() < 30:
            continue
        rho, p = stats.spearmanr(counts[ok], y[ok])
        rows.append(
            {
                "CELL_TYPE": ct,
                "N": int(ok.sum()),
                "SPEARMAN_RHO": float(rho),
                "P_SPEARMAN": float(p),
            }
        )
    return pd.DataFrame(rows)


def covariate_adjusted_assoc(
    counts: pd.Series,
    fractions: pd.DataFrame,
    cov: pd.DataFrame,
) -> pd.DataFrame:
    """fraction ~ scale(count) + IGHV + sex + age10."""
    rows = []
    x_count = counts.astype(float)
    for ct in fractions.columns:
        y = fractions[ct].astype(float)
        df = pd.DataFrame(
            {
                "y": y,
                "count": x_count,
                "ighv": cov["IGHV_MUTATED"],
                "sex": cov["SEX_MALE"],
                "age10": cov["AGE_YEARS"] / 10.0,
            }
        ).dropna()
        if len(df) < 40 or df["count"].nunique() < 2:
            continue
        sd = float(df["count"].std(ddof=0))
        if not np.isfinite(sd) or sd == 0:
            continue
        count_z = ((df["count"] - df["count"].mean()) / sd).to_numpy(dtype=float)
        design = np.column_stack(
            [
                np.ones(len(df)),
                count_z,
                df["ighv"].to_numpy(dtype=float),
                df["sex"].to_numpy(dtype=float),
                df["age10"].to_numpy(dtype=float),
            ]
        )
        beta, se, p = ols_count_coef(df["y"].to_numpy(dtype=float), design)
        rows.append(
            {
                "CELL_TYPE": ct,
                "N_ADJ": int(len(df)),
                "BETA_COUNT_Z": beta,
                "SE_COUNT_Z": se,
                "P_ADJ": p,
            }
        )
    return pd.DataFrame(rows)


def main() -> int:
    root = project_root()
    os.chdir(root)

    method = os.environ.get("DECONV_METHOD", "bayesprism").lower()
    splice_dir = root / "CLLCMA_alpha/bulk/cll_broad_2022_splicing/spladder_output"
    map_xlsx = (
        root
        / "CLLCMA_alpha/bulk/cll_broad_2022_splicing/clinical_outcome/sequences_and_participants.xlsx"
    )
    deconv_file = (
        root
        / f"CLLCMA_alpha/results/deconvolution/cll_broad_2022/{method}_sample_by_feature.tsv"
    )
    clin_sample = root / "CLLCMA_alpha/bulk/cll_broad_2022/data_clinical_sample.txt"
    clin_patient = root / "CLLCMA_alpha/bulk/cll_broad_2022/data_clinical_patient.txt"
    out_dir = (
        root
        / "CLLCMA_alpha/results/deconvolution/cll_broad_2022"
        / "splicing_outlier_tme_association"
        / method
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Method: {method} -> {out_dir}", flush=True)

    if not deconv_file.exists():
        print(f"Missing deconvolution matrix: {deconv_file}", file=sys.stderr)
        return 1

    micro = pd.read_csv(deconv_file, sep="\t")
    micro["SAMPLE_ID"] = micro["SAMPLE_ID"].astype(str).str.strip('"')
    feature_cols = [c for c in micro.columns if c != "SAMPLE_ID"]
    micro = attach_covariates(micro, clin_sample, clin_patient)

    run_map = load_run_participant_map(map_xlsx)
    # Build column map from first available matrix header
    example = splice_dir / f"simplified_{AS_TYPES[0]}.csv.gz"
    header = pd.read_csv(example, nrows=0).columns.tolist()
    psi_cols = [c for c in header if c != "event_id"]
    col_df = pd.DataFrame(
        {
            "col": psi_cols,
            "Run": [parse_srr_from_column(c) for c in psi_cols],
        }
    ).merge(run_map, on="Run", how="left")
    col_df = col_df[
        col_df["participant_id"].notna()
        & col_df["participant_id"].isin(micro["SAMPLE_ID"])
    ].copy()
    # one column per participant
    col_df = col_df.drop_duplicates(subset=["participant_id"], keep="first")
    use_cols = col_df["col"].tolist()
    sample_ids = col_df["participant_id"].tolist()
    print(
        f"Mapped splicing→deconv samples: {len(sample_ids)} "
        f"(map runs={len(run_map)}; deconv n={len(micro)})",
        flush=True,
    )
    col_df.to_csv(out_dir / "srr_to_sample_map.tsv", sep="\t", index=False)

    reuse = os.environ.get("SPLICING_OUTLIER_REUSE_BINARY", "").lower() in {
        "1", "true", "yes",
    }
    count_path = out_dir / "outlier_gene_counts_per_sample.tsv"
    summary_path = out_dir / "outlier_count_summary.tsv"

    count_table = pd.DataFrame({"SAMPLE_ID": sample_ids})
    type_meta_rows = []
    union_flags: dict[str, np.ndarray] = {}
    flags_by_type: dict[str, dict[str, np.ndarray]] = {}

    def load_binary(path: Path) -> dict[str, np.ndarray]:
        df = pd.read_csv(path, sep="\t", index_col=0)
        # align columns to sample_ids
        missing = [s for s in sample_ids if s not in df.columns]
        if missing:
            raise RuntimeError(f"{path.name} missing {len(missing)} sample columns")
        df = df.loc[:, sample_ids]
        return {str(g): df.loc[g].to_numpy(dtype=bool) for g in df.index}

    if reuse and all(
        (out_dir / f"outlier_gene_binary_{t}.tsv.gz").exists() for t in AS_TYPES
    ) and (out_dir / "outlier_gene_binary_all_types_union.tsv.gz").exists() and count_path.exists():
        print("Reusing saved gene-binary matrices + counts", flush=True)
        count_table = pd.read_csv(count_path, sep="\t")
        # keep only current sample order
        count_table = count_table.set_index("SAMPLE_ID").reindex(sample_ids).reset_index()
        for as_type in AS_TYPES:
            print(f"=== load {as_type} ===", flush=True)
            flags_by_type[as_type] = load_binary(
                out_dir / f"outlier_gene_binary_{as_type}.tsv.gz"
            )
        union_flags = load_binary(
            out_dir / "outlier_gene_binary_all_types_union.tsv.gz"
        )
        flags_by_type["all_types_union"] = union_flags
        if summary_path.exists():
            summary = pd.read_csv(summary_path, sep="\t")
        else:
            summary = pd.DataFrame()
        as_cols = [c for c in AS_TYPES if c in count_table.columns]
    else:
        for as_type in AS_TYPES:
            mat = splice_dir / f"simplified_{as_type}.csv.gz"
            info = splice_dir / f"simplified_{as_type}_event_info.csv.gz"
            if not mat.exists() or not info.exists():
                print(f"Skipping {as_type}: missing files", flush=True)
                continue
            print(f"=== {as_type} ===", flush=True)
            counts, gene_flags, meta = process_as_type(mat, info, use_cols, sample_ids)
            count_table[as_type] = counts
            meta["AS_TYPE"] = as_type
            type_meta_rows.append(meta)
            flags_by_type[as_type] = gene_flags
            write_gene_binary(
                out_dir / f"outlier_gene_binary_{as_type}.tsv.gz",
                gene_flags,
                sample_ids,
            )
            print(
                f"  tested events={meta['n_events_tested']:,}/{meta['n_events']:,}; "
                f"outlier genes={meta['n_outlier_genes_any_sample']:,}; "
                f"count median={meta['count_median']:.1f} "
                f"[{meta['count_min']}-{meta['count_max']}]",
                flush=True,
            )
            for g, arr in gene_flags.items():
                if g in union_flags:
                    union_flags[g] |= arr
                else:
                    union_flags[g] = arr.copy()

        count_table["all_types_union"] = counts_from_gene_flags(
            union_flags, len(sample_ids)
        )
        as_cols = [c for c in AS_TYPES if c in count_table.columns]
        count_table["all_types_sum"] = count_table[as_cols].sum(axis=1)
        count_table.to_csv(count_path, sep="\t", index=False)
        flags_by_type["all_types_union"] = union_flags
        write_gene_binary(
            out_dir / "outlier_gene_binary_all_types_union.tsv.gz",
            union_flags,
            sample_ids,
        )

        summary = pd.DataFrame(type_meta_rows)
        summary.loc[len(summary)] = {
            "AS_TYPE": "all_types_union",
            "n_events": np.nan,
            "n_events_tested": np.nan,
            "n_event_sample_outlier_flags": np.nan,
            "n_outlier_genes_any_sample": len(union_flags),
            "count_min": int(count_table["all_types_union"].min()),
            "count_median": float(count_table["all_types_union"].median()),
            "count_mean": float(count_table["all_types_union"].mean()),
            "count_max": int(count_table["all_types_union"].max()),
        }
        summary.to_csv(summary_path, sep="\t", index=False)

    # Merge TME for association
    merged = count_table.merge(micro, on="SAMPLE_ID", how="inner")
    print(f"Association samples: {len(merged)}", flush=True)

    assoc_rows = []
    burden_cols = as_cols + ["all_types_union", "all_types_sum"]
    for burden in burden_cols:
        sp = spearman_assoc(merged[burden], merged[feature_cols])
        adj = covariate_adjusted_assoc(
            merged[burden],
            merged[feature_cols],
            merged[["IGHV_MUTATED", "SEX_MALE", "AGE_YEARS"]],
        )
        tab = sp.merge(adj, on="CELL_TYPE", how="outer")
        tab.insert(0, "BURDEN", burden)
        assoc_rows.append(tab)

    assoc = pd.concat(assoc_rows, ignore_index=True)
    # FDR within burden measure across cell types (primary: adjusted p)
    assoc["FDR_ADJ"] = np.nan
    assoc["FDR_SPEARMAN"] = np.nan
    for burden, idx in assoc.groupby("BURDEN").groups.items():
        ii = list(idx)
        p_adj = assoc.loc[ii, "P_ADJ"].to_numpy()
        p_sp = assoc.loc[ii, "P_SPEARMAN"].to_numpy()
        assoc.loc[ii, "FDR_ADJ"] = bh_fdr(p_adj)
        assoc.loc[ii, "FDR_SPEARMAN"] = bh_fdr(p_sp)

    assoc = assoc.sort_values(["BURDEN", "P_ADJ"], na_position="last")
    assoc.to_csv(out_dir / "association_all.tsv", sep="\t", index=False)
    sig = assoc[assoc["FDR_ADJ"] < 0.05].copy()
    sig.to_csv(out_dir / "association_significant_fdr05.tsv", sep="\t", index=False)

    # ---------------- Per-gene associations ----------------
    print("=== Per-gene outlier × TME associations ===", flush=True)
    per_gene_parts = []
    for as_type, gflags in flags_by_type.items():
        part = per_gene_tme_association(
            gflags,
            sample_ids,
            micro,
            feature_cols,
            as_type=as_type,
        )
        if len(part):
            per_gene_parts.append(part)
    if per_gene_parts:
        per_gene = pd.concat(per_gene_parts, ignore_index=True)
        per_gene["FDR_ADJ"] = np.nan
        per_gene["FDR_WILCOX"] = np.nan
        for as_type, idx in per_gene.groupby("AS_TYPE").groups.items():
            ii = list(idx)
            per_gene.loc[ii, "FDR_ADJ"] = bh_fdr(per_gene.loc[ii, "P_ADJ"].to_numpy())
            per_gene.loc[ii, "FDR_WILCOX"] = bh_fdr(
                per_gene.loc[ii, "P_WILCOX"].to_numpy()
            )
        per_gene = per_gene.sort_values(
            ["AS_TYPE", "P_ADJ"], na_position="last"
        )
        per_gene_path = out_dir / "association_per_gene_all.tsv.gz"
        per_gene.to_csv(per_gene_path, sep="\t", index=False, compression="gzip")
        per_sig = per_gene[per_gene["FDR_ADJ"] < 0.05].copy()
        per_sig.to_csv(
            out_dir / "association_per_gene_significant_fdr05.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )
        # Best cell-type hit per gene within each AS type
        best = (
            per_gene.sort_values("P_ADJ")
            .groupby(["AS_TYPE", "GENE"], as_index=False)
            .first()
            .sort_values(["AS_TYPE", "P_ADJ"])
        )
        best.to_csv(
            out_dir / "association_per_gene_best_per_gene.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )
        per_gene_summary = (
            per_gene.groupby("AS_TYPE", as_index=False)
            .agg(
                N_TESTS=("P_ADJ", "size"),
                N_GENES=("GENE", "nunique"),
                N_SIG_FDR05=("FDR_ADJ", lambda s: int((s < 0.05).sum())),
                N_SIG_FDR10=("FDR_ADJ", lambda s: int((s < 0.10).sum())),
            )
        )
        per_gene_summary.to_csv(
            out_dir / "association_per_gene_summary.tsv", sep="\t", index=False
        )
        n_per_sig = int(len(per_sig))
        n_per_tests = int(len(per_gene))
    else:
        n_per_sig = 0
        n_per_tests = 0
        per_gene = pd.DataFrame()
        per_sig = pd.DataFrame()

    # Key findings
    lines = [
        "SPLICING OUTLIER GENE BURDEN × TME ASSOCIATION",
        f"Deconvolution method: {method}",
        f"Samples (splicing ∩ deconv): {len(merged)}",
        (
            "Outlier rule: PSI ≤P10 or ≥P90 (across analysis samples) "
            f"AND |PSI−median|>{MEDIAN_DELTA}; min finite n/event={MIN_N_FINITE}"
        ),
        "Gene collapse: gene outlier if any event outlier; count = sum of genes",
        "Burden association: fraction ~ z(outlier_gene_count) + IGHV + sex + age/10",
        f"Burden covariate-adjusted FDR<0.05 hits: {len(sig)}",
        (
            f"Per-gene association: fraction ~ outlier(0/1) + IGHV + sex + age/10; "
            f"min outlier/WT samples = {MIN_GENE_OUTLIER}/{MIN_GENE_WT}"
        ),
        f"Per-gene tests: {n_per_tests:,}; FDR_ADJ<0.05 hits: {n_per_sig:,}",
        "",
        "=== Count summary (median [min-max]) ===",
    ]
    for _, r in summary.iterrows():
        lines.append(
            f"{r['AS_TYPE']}: {r['count_median']:.0f} [{r['count_min']}-{r['count_max']}] "
            f"(genes ever outlier={r['n_outlier_genes_any_sample']})"
        )
    lines.append("")
    lines.append("=== Top burden associations by adjusted p ===")
    top = assoc.nsmallest(15, "P_ADJ", keep="all")
    for _, r in top.iterrows():
        lines.append(
            f"{r['BURDEN']} ~ {r['CELL_TYPE']} | rho={r.get('SPEARMAN_RHO', np.nan):+.3f} "
            f"| beta_z={r.get('BETA_COUNT_Z', np.nan):+.3g} "
            f"| adj-p={r.get('P_ADJ', np.nan):.3g} FDR={r.get('FDR_ADJ', np.nan):.3g}"
        )
    if len(per_gene):
        lines.append("")
        lines.append("=== Per-gene summary (FDR within AS type) ===")
        for _, r in per_gene_summary.iterrows():
            lines.append(
                f"{r['AS_TYPE']}: genes={r['N_GENES']:,} tests={r['N_TESTS']:,} "
                f"FDR<0.05={r['N_SIG_FDR05']:,} FDR<0.10={r['N_SIG_FDR10']:,}"
            )
        lines.append("")
        lines.append("=== Top per-gene associations by adjusted p ===")
        topg = per_gene.nsmallest(25, "P_ADJ", keep="all")
        for _, r in topg.iterrows():
            lines.append(
                f"{r['AS_TYPE']} | {r['GENE']} ~ {r['CELL_TYPE']} | "
                f"n_out={r['N_OUTLIER']} | delta={r['CLIFFS_DELTA']:+.2f} | "
                f"beta={r['BETA_ADJ']:+.3g} | adj-p={r['P_ADJ']:.3g} "
                f"FDR={r['FDR_ADJ']:.3g}"
            )
    (out_dir / "KEY_FINDINGS.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[:30]), flush=True)
    print(f"Wrote {out_dir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
