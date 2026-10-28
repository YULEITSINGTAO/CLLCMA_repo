#!/usr/bin/env python3
"""Per-cell-type gene ranking and enrichment for splicing-outlier × TME.

For each cell type (default: all_types_union associations):
  score = Cliff's delta * (-log10(P_ADJ))
  rank genes by score (rank 1 = highest score)
  plot rank vs score
  preranked GSEA (local GMT) + ORA on top/bottom tails

Env:
  CLL_TME_DIR, DECONV_METHOD (default bayesprism)
  HEATMAP_AS_TYPE (default all_types_union)
  ORA_TOP_N (default 200)
"""

from __future__ import annotations

import os
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

AS_DEFAULT = "all_types_union"
ORA_TOP_N_DEFAULT = 200
MIN_SET_SIZE = 10
MAX_SET_SIZE = 500


def project_root() -> Path:
    env = os.environ.get("CLL_TME_DIR")
    if env:
        return Path(env).expanduser().resolve()
    return Path(__file__).resolve().parents[2]


def map_ensembl_to_symbol(ensembl_ids: list[str]) -> dict[str, str]:
    """Map ENSG → SYMBOL via R org.Hs.eg.db."""
    import subprocess
    import tempfile

    ids = sorted({i for i in ensembl_ids if isinstance(i, str) and i.startswith("ENSG")})
    if not ids:
        return {}
    rscript = os.environ.get(
        "RSCRIPT",
        str(Path.home() / "miniconda3/envs/deconvolution/bin/Rscript"),
    )
    with tempfile.TemporaryDirectory() as td:
        td_path = Path(td)
        id_path = td_path / "ids.txt"
        out_path = td_path / "map.tsv"
        r_path = td_path / "map_ids.R"
        id_path.write_text("\n".join(ids) + "\n")
        r_path.write_text(
            "\n".join(
                [
                    "suppressPackageStartupMessages(library(org.Hs.eg.db))",
                    f'ids <- scan("{id_path.as_posix()}", what = character(), quiet = TRUE)',
                    "m <- AnnotationDbi::select(",
                    "  org.Hs.eg.db, keys = ids, keytype = 'ENSEMBL', columns = 'SYMBOL'",
                    ")",
                    "m <- unique(m[!is.na(m$SYMBOL) & m$SYMBOL != '', c('ENSEMBL', 'SYMBOL')])",
                    f'write.table(m, "{out_path.as_posix()}", sep = "\\t", row.names = FALSE, quote = FALSE)',
                    "",
                ]
            )
        )
        subprocess.check_call([rscript, str(r_path)])
        df = pd.read_csv(out_path, sep="\t")
    return dict(zip(df["ENSEMBL"].astype(str), df["SYMBOL"].astype(str)))


def read_gmt(path: Path) -> dict[str, set[str]]:
    sets: dict[str, set[str]] = {}
    with path.open() as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            name = parts[0]
            genes = {g.strip().upper() for g in parts[2:] if g.strip()}
            if genes:
                sets[name] = genes
    return sets


def ora_one_list(
    hit_genes: set[str],
    background: set[str],
    gene_sets: dict[str, set[str]],
    direction: str,
    cell_type: str,
    library: str,
) -> pd.DataFrame:
    bg = {g.upper() for g in background}
    hits = {g.upper() for g in hit_genes} & bg
    n_bg = len(bg)
    n_hit = len(hits)
    rows = []
    if n_hit < 5 or n_bg < 50:
        return pd.DataFrame()
    for term, members in gene_sets.items():
        s = {g.upper() for g in members} & bg
        k = len(s)
        if k < MIN_SET_SIZE or k > MAX_SET_SIZE:
            continue
        a = len(hits & s)
        if a < 2:
            continue
        # Fisher: [[in_set & hit, hit & not set],[set & not hit, neither]]
        table = np.array(
            [
                [a, n_hit - a],
                [k - a, n_bg - n_hit - (k - a)],
            ]
        )
        if table.min() < 0:
            continue
        oddsr, p = fisher_exact(table, alternative="greater")
        rows.append(
            {
                "CELL_TYPE": cell_type,
                "DIRECTION": direction,
                "LIBRARY": library,
                "TERM": term,
                "OVERLAP": a,
                "SET_SIZE_IN_BG": k,
                "HIT_SIZE": n_hit,
                "BG_SIZE": n_bg,
                "ODDS_RATIO": float(oddsr),
                "P_VALUE": float(p),
                "GENES": ";".join(sorted(hits & s)),
            }
        )
    if not rows:
        return pd.DataFrame()
    out = pd.DataFrame(rows)
    # BH within this cell_type × direction × library
    p = out["P_VALUE"].to_numpy()
    order = np.argsort(p)
    n = len(p)
    q = np.empty(n)
    ranked = p[order]
    qvals = ranked * n / np.arange(1, n + 1)
    qvals = np.minimum.accumulate(qvals[::-1])[::-1]
    qvals = np.clip(qvals, 0, 1)
    q[order] = qvals
    out["FDR"] = q
    return out.sort_values("P_VALUE")


def plot_rank_curves(
    ranked: pd.DataFrame,
    out_pdf: Path,
    out_png: Path,
    method_label: str,
    as_type: str,
) -> None:
    """Rank vs score panels; positive scores in vermillion, negative in blue."""
    # Diverging red / blue (positive / negative), high contrast on white
    col_pos = "#B2182B"  # red
    col_neg = "#2166AC"  # blue
    col_zero = "#666666"
    col_spine = "#333333"

    ct_order = [
        "CLL",
        "CD4_T_naive", "CD4_T_memory",
        "CD8_T_naive", "CD8_T_memory", "CD8aa", "DN_T",
        "T_proliferating", "Treg", "gdT", "MAIT",
        "NK_CD56dim", "NK_CD56bright", "NK_proliferating",
        "Monocyte_CD14", "Monocyte_CD16", "Monocyte_intermediate",
        "cDC1", "cDC2", "ASDC", "ILC", "Progenitor",
    ]
    present = list(ranked["CELL_TYPE"].unique())
    cts = [c for c in ct_order if c in present] + sorted(
        set(present) - set(ct_order)
    )

    n = len(cts)
    ncols = 4
    nrows = int(np.ceil(n / ncols))
    fig_w = 180 / 25.4  # ~180 mm wide provisional general figure
    fig_h = max(140 / 25.4, 32 / 25.4 * nrows)
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(fig_w, fig_h),
        layout="constrained",
        sharex=False,
        sharey=False,
    )
    axes = np.atleast_1d(axes).ravel()

    def _split_signed_curve(x: np.ndarray, y: np.ndarray):
        """Return (x_pos, y_pos), (x_neg, y_neg) with zero-crossing interpolated."""
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        xp, yp, xn, yn = [], [], [], []
        for i in range(len(y) - 1):
            x0, x1 = x[i], x[i + 1]
            y0, y1 = y[i], y[i + 1]
            if y0 >= 0:
                xp.append(x0)
                yp.append(y0)
            else:
                xn.append(x0)
                yn.append(y0)
            # crossing
            if y0 == 0 or y1 == 0 or (y0 > 0) == (y1 > 0):
                continue
            t = -y0 / (y1 - y0)
            xc = x0 + t * (x1 - x0)
            xp.append(xc)
            yp.append(0.0)
            xn.append(xc)
            yn.append(0.0)
        # last point
        if y[-1] >= 0:
            xp.append(x[-1])
            yp.append(y[-1])
        else:
            xn.append(x[-1])
            yn.append(y[-1])
        return (
            np.asarray(xp),
            np.asarray(yp),
            np.asarray(xn),
            np.asarray(yn),
        )

    for i, ct in enumerate(cts):
        ax = axes[i]
        d = ranked[ranked["CELL_TYPE"] == ct].sort_values("RANK")
        x = d["RANK"].to_numpy()
        y = d["SCORE"].to_numpy()
        xp, yp, xn, yn = _split_signed_curve(x, y)

        ax.axhline(0, color=col_zero, lw=0.7, zorder=1)
        if len(xp):
            ax.fill_between(xp, 0, yp, color=col_pos, alpha=0.18, linewidth=0, zorder=2)
            ax.plot(xp, yp, color=col_pos, lw=1.1, solid_capstyle="round", zorder=3)
        if len(xn):
            ax.fill_between(xn, 0, yn, color=col_neg, alpha=0.18, linewidth=0, zorder=2)
            ax.plot(xn, yn, color=col_neg, lw=1.1, solid_capstyle="round", zorder=3)

        # Label extremes (redundant to color)
        top = d.nsmallest(2, "RANK")
        bot = d.nlargest(2, "RANK")
        for _, r in top.iterrows():
            ax.annotate(
                str(r["SYMBOL"])[:10],
                xy=(r["RANK"], r["SCORE"]),
                xytext=(4, 2),
                textcoords="offset points",
                fontsize=5.5,
                color=col_pos,
                ha="left",
                va="bottom",
            )
        for _, r in bot.iterrows():
            ax.annotate(
                str(r["SYMBOL"])[:10],
                xy=(r["RANK"], r["SCORE"]),
                xytext=(-4, -2),
                textcoords="offset points",
                fontsize=5.5,
                color=col_neg,
                ha="right",
                va="top",
            )

        ax.set_title(ct.replace("_", " "), fontsize=8.5, fontweight="bold", pad=3)
        ax.set_xlabel("Gene rank", fontsize=7)
        ax.set_ylabel(r"Score ($\delta \times -\log_{10}P$)", fontsize=7)
        ax.tick_params(labelsize=6, length=2.5, width=0.6)
        for spine in ax.spines.values():
            spine.set_color(col_spine)
            spine.set_linewidth(0.6)
        ax.set_facecolor("white")
        ax.margins(x=0.02)

    for j in range(len(cts), len(axes)):
        axes[j].axis("off")

    # Shared legend
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    legend_elems = [
        Patch(facecolor=col_pos, edgecolor=col_pos, alpha=0.35, label="Positive score"),
        Patch(facecolor=col_neg, edgecolor=col_neg, alpha=0.35, label="Negative score"),
        Line2D([0], [0], color=col_zero, lw=0.8, label="Score = 0"),
    ]
    fig.legend(
        handles=legend_elems,
        loc="lower center",
        ncol=3,
        fontsize=8,
        frameon=False,
        bbox_to_anchor=(0.5, -0.02),
    )

    fig.suptitle(
        f"Splicing-outlier gene ranking by cell type ({method_label}; {as_type})",
        fontsize=11,
        fontweight="bold",
    )
    fig.text(
        0.5,
        0.965,
        r"score = Cliff's $\delta \times (-\log_{10} P_{\mathrm{adj}})$; "
        "rank 1 = highest score (higher fraction in gene-outlier samples)",
        ha="center",
        va="top",
        fontsize=8,
        color="#333333",
    )

    fig.savefig(out_pdf, facecolor="white")
    fig.savefig(out_png, dpi=400, facecolor="white")
    plt.close(fig)

    # Provenance sidecar
    prov = out_pdf.with_name(out_pdf.stem + "_provenance.txt")
    prov.write_text(
        "\n".join(
            [
                "Provisional general figure (journal/phase not specified).",
                f"Source: gene_scores_ranked_by_celltype.tsv.gz ({as_type})",
                "Encoding: x=gene rank; y=Cliff's delta * (-log10 P_adj)",
                "Color: red (#B2182B) positive; blue (#2166AC) negative (ColorBrewer RdBu ends)",
                "Fill: signed area between curve and y=0 (alpha 0.18)",
                "Labels: top-2 and bottom-2 genes by rank per panel",
                f"Exported: {pd.Timestamp.now()}",
                "",
            ]
        )
    )


def plot_enrichment_dots(enr: pd.DataFrame, out_pdf: Path, title: str) -> None:
    if enr is None or enr.empty:
        return
    # top terms overall after filtering
    d = enr[enr["FDR"] < 0.05].copy() if "FDR" in enr.columns else enr.copy()
    if d.empty:
        d = enr.nsmallest(30, "P_VALUE" if "P_VALUE" in enr.columns else "FDR q-val")
    # handle ORA vs GSEA column names
    if "TERM" in d.columns:
        d = d.sort_values("FDR").groupby("CELL_TYPE", group_keys=False).head(5)
        d["NEGLOG10"] = -np.log10(d["FDR"].clip(lower=1e-300))
        d["TERM_SHORT"] = d["TERM"].str.slice(0, 60)
        fig_h = max(4, 0.28 * len(d) + 1.5)
        fig, ax = plt.subplots(figsize=(10, fig_h))
        # order
        d = d.sort_values(["CELL_TYPE", "FDR"], ascending=[True, True])
        d["Y"] = np.arange(len(d))
        colors = d["DIRECTION"].map({"positive": "#D55E00", "negative": "#0072B2"}).fillna("#333")
        ax.scatter(d["NEGLOG10"], d["Y"], s=30 + 8 * d["OVERLAP"].astype(float), c=colors, alpha=0.85)
        ax.set_yticks(d["Y"])
        ax.set_yticklabels(
            [f"{r.CELL_TYPE} | {r.DIRECTION[:3]} | {r.TERM_SHORT}" for r in d.itertuples()],
            fontsize=6,
        )
        ax.set_xlabel(r"$-\log_{10}$(FDR)")
        ax.set_title(title, fontsize=10)
        fig.tight_layout()
        fig.savefig(out_pdf, bbox_inches="tight")
        plt.close(fig)
    elif "Term" in d.columns:
        # GSEA res2d
        d = d.sort_values("FDR q-val").groupby("CELL_TYPE", group_keys=False).head(5)
        d["NEGLOG10"] = -np.log10(d["FDR q-val"].clip(lower=1e-300))
        d["TERM_SHORT"] = d["Term"].astype(str).str.slice(0, 60)
        fig_h = max(4, 0.28 * len(d) + 1.5)
        fig, ax = plt.subplots(figsize=(10, fig_h))
        d = d.sort_values(["CELL_TYPE", "FDR q-val"])
        d["Y"] = np.arange(len(d))
        ax.scatter(
            d["NES"],
            d["Y"],
            s=40,
            c=np.where(d["NES"] >= 0, "#D55E00", "#0072B2"),
            alpha=0.85,
        )
        ax.axvline(0, color="#888", lw=0.6)
        ax.set_yticks(d["Y"])
        ax.set_yticklabels(
            [f"{r.CELL_TYPE} | {r.TERM_SHORT}" for r in d.itertuples()],
            fontsize=6,
        )
        ax.set_xlabel("NES")
        ax.set_title(title, fontsize=10)
        fig.tight_layout()
        fig.savefig(out_pdf, bbox_inches="tight")
        plt.close(fig)


def main() -> int:
    root = project_root()
    os.chdir(root)
    method = os.environ.get("DECONV_METHOD", "bayesprism").lower()
    as_type = os.environ.get("HEATMAP_AS_TYPE", AS_DEFAULT)
    ora_top_n = int(os.environ.get("ORA_TOP_N", str(ORA_TOP_N_DEFAULT)))
    method_label = {
        "bayesprism": "BayesPrism",
        "cibersort": "CIBERSORT",
    }.get(method, method)

    in_dir = (
        root
        / "CLLCMA_alpha/results/deconvolution/cll_broad_2022"
        / "splicing_outlier_tme_association"
        / method
    )
    out_dir = in_dir / "rank_enrichment"
    out_dir.mkdir(parents=True, exist_ok=True)
    rnk_dir = out_dir / "rnk"
    rnk_dir.mkdir(exist_ok=True)

    assoc_file = in_dir / "association_per_gene_all.tsv.gz"
    print(f"Reading {assoc_file}", flush=True)
    pg = pd.read_csv(assoc_file, sep="\t")
    pg = pg[pg["AS_TYPE"] == as_type].copy()
    pg = pg[np.isfinite(pg["CLIFFS_DELTA"]) & np.isfinite(pg["P_ADJ"])]
    if pg.empty:
        raise SystemExit(f"No rows for AS_TYPE={as_type}")

    # Map symbols
    print("Mapping Ensembl → symbol …", flush=True)
    sym = map_ensembl_to_symbol(pg["GENE"].astype(str).unique().tolist())
    pg["SYMBOL"] = pg["GENE"].map(sym)
    # drop unmapped for enrichment; keep ENSG as last resort for ranking plots
    pg["SYMBOL_PLOT"] = pg["SYMBOL"].fillna(pg["GENE"])

    # Score and rank within cell type
    pg["NEGLOG10_P"] = -np.log10(np.clip(pg["P_ADJ"].to_numpy(), 1e-300, None))
    pg["SCORE"] = pg["CLIFFS_DELTA"] * pg["NEGLOG10_P"]
    # one row per gene per cell type already
    parts = []
    for ct, sub in pg.groupby("CELL_TYPE"):
        s = sub.sort_values("SCORE", ascending=False).copy()
        # if duplicate symbols, keep strongest |score|
        s = s.sort_values("SCORE", key=lambda x: x.abs(), ascending=False)
        s = s.drop_duplicates(subset=["GENE"], keep="first")
        s = s.sort_values("SCORE", ascending=False).reset_index(drop=True)
        s["RANK"] = np.arange(1, len(s) + 1)
        parts.append(s)
    ranked = pd.concat(parts, ignore_index=True)
    ranked.to_csv(out_dir / "gene_scores_ranked_by_celltype.tsv.gz", sep="\t", index=False, compression="gzip")
    print(f"Wrote ranked scores ({len(ranked):,} rows)", flush=True)

    plot_rank_curves(
        ranked.assign(SYMBOL=ranked["SYMBOL_PLOT"]),
        out_dir / "rank_vs_score_by_celltype.pdf",
        out_dir / "rank_vs_score_by_celltype.png",
        method_label,
        as_type,
    )
    print("Wrote rank-vs-score plots", flush=True)

    # GMT libraries (local, offline)
    gmt_paths = {
        "GO_BP": root
        / "CLLCMA_alpha/results/bcell_stage_deseq2/publication_figures/stage_progression/GO_Biological_Process_2026_retrieved_2026-09-08.gmt",
        "Reactome": root
        / "CLLCMA_alpha/results/bcell_stage_deseq2/publication_figures/stage_progression/Reactome_Pathways_2024_retrieved_2026-09-08.gmt",
    }
    # fallback to gseapy cache
    cache = Path.home() / ".cache/gseapy"
    if not gmt_paths["GO_BP"].exists():
        alt = cache / "Enrichr.GO_Biological_Process_2026.gmt"
        if alt.exists():
            gmt_paths["GO_BP"] = alt
    if not gmt_paths["Reactome"].exists():
        alt = cache / "Enrichr.Reactome_Pathways_2024.gmt"
        if alt.exists():
            gmt_paths["Reactome"] = alt

    gene_sets_all = {}
    for lib, path in gmt_paths.items():
        if path.exists():
            gene_sets_all[lib] = read_gmt(path)
            print(f"Loaded {lib}: {len(gene_sets_all[lib])} sets from {path.name}", flush=True)
        else:
            print(f"WARNING: missing GMT {path}", flush=True)

    # ---- ORA on positive / negative tails ----
    ora_rows = []
    for ct, sub in ranked.groupby("CELL_TYPE"):
        sub_sym = sub.dropna(subset=["SYMBOL"]).copy()
        sub_sym["SYMBOL"] = sub_sym["SYMBOL"].str.upper()
        # collapse duplicate symbols by max |score|
        sub_sym = (
            sub_sym.assign(ABS=sub_sym["SCORE"].abs())
            .sort_values("ABS", ascending=False)
            .drop_duplicates("SYMBOL")
        )
        bg = set(sub_sym["SYMBOL"])
        pos = set(sub_sym.nlargest(ora_top_n, "SCORE")["SYMBOL"])
        neg = set(sub_sym.nsmallest(ora_top_n, "SCORE")["SYMBOL"])
        for lib, gsets in gene_sets_all.items():
            ora_rows.append(ora_one_list(pos, bg, gsets, "positive", ct, lib))
            ora_rows.append(ora_one_list(neg, bg, gsets, "negative", ct, lib))
    ora = pd.concat([r for r in ora_rows if r is not None and len(r)], ignore_index=True)
    if len(ora):
        ora.to_csv(out_dir / "ora_by_celltype.tsv.gz", sep="\t", index=False, compression="gzip")
        ora_sig = ora[ora["FDR"] < 0.05].copy()
        ora_sig.to_csv(out_dir / "ora_significant_fdr05.tsv.gz", sep="\t", index=False, compression="gzip")
        plot_enrichment_dots(
            ora_sig if len(ora_sig) else ora.nsmallest(40, "P_VALUE"),
            out_dir / "ora_top_terms_by_celltype.pdf",
            f"ORA on score tails (top/bottom {ora_top_n}; {as_type})",
        )
        print(f"ORA: {len(ora):,} tests; FDR<0.05: {len(ora_sig):,}", flush=True)
    else:
        print("ORA produced no rows", flush=True)

    # ---- Preranked GSEA via gseapy (local GMT) ----
    try:
        import gseapy as gp
    except ImportError:
        print("gseapy not available; skipping preranked GSEA", flush=True)
        gp = None

    gsea_rows = []
    if gp is not None and gene_sets_all:
        # merge libraries into one dict for prerank, keep library tag in term via prefix
        for ct, sub in ranked.groupby("CELL_TYPE"):
            sub_sym = sub.dropna(subset=["SYMBOL"]).copy()
            sub_sym["SYMBOL"] = sub_sym["SYMBOL"].str.upper()
            sub_sym = (
                sub_sym.assign(ABS=sub_sym["SCORE"].abs())
                .sort_values("ABS", ascending=False)
                .drop_duplicates("SYMBOL")
            )
            rnk = sub_sym.set_index("SYMBOL")["SCORE"].sort_values(ascending=False)
            # write rnk
            safe_ct = re.sub(r"[^A-Za-z0-9._-]+", "_", ct)
            rnk.to_csv(rnk_dir / f"{safe_ct}.rnk", sep="\t", header=False)
            for lib, gsets in gene_sets_all.items():
                # gseapy wants dict of list
                gs = {k: sorted(v) for k, v in gsets.items() if MIN_SET_SIZE <= len(v) <= MAX_SET_SIZE}
                try:
                    pre = gp.prerank(
                        rnk=rnk,
                        gene_sets=gs,
                        min_size=MIN_SET_SIZE,
                        max_size=MAX_SET_SIZE,
                        permutation_num=200,
                        seed=1,
                        threads=4,
                        outdir=None,
                        verbose=False,
                    )
                    res = pre.res2d.copy()
                    res.insert(0, "CELL_TYPE", ct)
                    res.insert(1, "LIBRARY", lib)
                    gsea_rows.append(res)
                    print(f"  GSEA {ct} / {lib}: {len(res)} sets", flush=True)
                except Exception as e:
                    print(f"  GSEA failed {ct}/{lib}: {e}", flush=True)

    if gsea_rows:
        gsea = pd.concat(gsea_rows, ignore_index=True)
        # normalize column names
        gsea.to_csv(out_dir / "gsea_prerank_by_celltype.tsv.gz", sep="\t", index=False, compression="gzip")
        fdr_col = "FDR q-val" if "FDR q-val" in gsea.columns else "FDR"
        gsea_sig = gsea[pd.to_numeric(gsea[fdr_col], errors="coerce") < 0.25].copy()
        gsea_sig.to_csv(
            out_dir / "gsea_prerank_fdr025.tsv.gz", sep="\t", index=False, compression="gzip"
        )
        plot_enrichment_dots(
            gsea_sig if len(gsea_sig) else gsea.head(40),
            out_dir / "gsea_top_terms_by_celltype.pdf",
            f"Preranked GSEA (score=δ×-log10P; {as_type}; FDR displayed <0.25)",
        )
        print(f"GSEA rows: {len(gsea):,}; FDR q<0.25: {len(gsea_sig):,}", flush=True)

    # Key findings
    lines = [
        "SPLICING-OUTLIER GENE RANKING + ENRICHMENT",
        f"Method: {method_label}",
        f"AS matrix: {as_type}",
        "Score: Cliff's delta * (-log10 P_ADJ)",
        "Rank 1 = highest score (higher TME fraction in gene-outlier samples)",
        f"ORA tails: top/bottom {ora_top_n} genes by score",
        f"Samples/genes follow association_per_gene_all.tsv.gz",
    ]
    if len(ora):
        lines.append(f"ORA tests: {len(ora)}; FDR<0.05: {(ora['FDR']<0.05).sum()}")
        top_ora = ora.nsmallest(15, "P_VALUE")
        lines.append("=== Top ORA terms ===")
        for _, r in top_ora.iterrows():
            lines.append(
                f"{r['CELL_TYPE']} | {r['DIRECTION']} | {r['LIBRARY']} | {r['TERM'][:70]} | "
                f"overlap={r['OVERLAP']} FDR={r['FDR']:.3g}"
            )
    if gsea_rows:
        gsea = pd.concat(gsea_rows, ignore_index=True)
        fdr_col = "FDR q-val" if "FDR q-val" in gsea.columns else "FDR"
        gsea[fdr_col] = pd.to_numeric(gsea[fdr_col], errors="coerce")
        lines.append("=== Top GSEA terms (by FDR q) ===")
        for _, r in gsea.nsmallest(15, fdr_col).iterrows():
            term = str(r.get("Term", r.get("TERM", "")))[:70]
            lines.append(
                f"{r['CELL_TYPE']} | {r['LIBRARY']} | {term} | NES={r.get('NES', np.nan):+.3f} "
                f"FDR={r[fdr_col]:.3g}"
            )
    (out_dir / "KEY_FINDINGS.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[:25]), flush=True)
    print(f"Done -> {out_dir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
