#!/usr/bin/env python3
"""Combined rank-score + preranked GSEA bubble figures per cell type.

For each cell type with significant GSEA terms (default FDR q < 0.05):
  one row = [rank-vs-score curve | GSEA bubble plot].
Cell types without significant terms are omitted.
Multiple pages if needed.

Env:
  CLL_TME_DIR, DECONV_METHOD
  GSEA_FDR_CUTOFF (default 0.05)
  GSEA_TOP_N_TERMS (default 10)
  COMBINED_ROWS_PER_PAGE (default 5)
"""

from __future__ import annotations

import os
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


COL_POS = "#B2182B"
COL_NEG = "#2166AC"
COL_ZERO = "#666666"
COL_SPINE = "#333333"

CT_ORDER = [
    "CLL",
    "CD4_T_naive",
    "CD4_T_memory",
    "CD8_T_naive",
    "CD8_T_memory",
    "CD8aa",
    "DN_T",
    "T_proliferating",
    "Treg",
    "gdT",
    "MAIT",
    "NK_CD56dim",
    "NK_CD56bright",
    "NK_proliferating",
    "Monocyte_CD14",
    "Monocyte_CD16",
    "Monocyte_intermediate",
    "cDC1",
    "cDC2",
    "ASDC",
    "ILC",
    "Progenitor",
]


def project_root() -> Path:
    env = os.environ.get("CLL_TME_DIR")
    if env:
        return Path(env).expanduser().resolve()
    return Path(__file__).resolve().parents[2]


def split_signed_curve(x: np.ndarray, y: np.ndarray):
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
        if y0 == 0 or y1 == 0 or (y0 > 0) == (y1 > 0):
            continue
        t = -y0 / (y1 - y0)
        xc = x0 + t * (x1 - x0)
        xp.extend([xc])
        yp.extend([0.0])
        xn.extend([xc])
        yn.extend([0.0])
    if y[-1] >= 0:
        xp.append(x[-1])
        yp.append(y[-1])
    else:
        xn.append(x[-1])
        yn.append(y[-1])
    return map(np.asarray, (xp, yp, xn, yn))


def short_term(term: str, max_len: int = 55) -> str:
    t = re.sub(r"\s+", " ", str(term)).strip()
    # drop trailing GO ids for readability when long
    t = re.sub(r"\s*\(GO:\d+\)\s*$", "", t)
    if len(t) > max_len:
        t = t[: max_len - 1] + "…"
    return t


def draw_rank_panel(ax, d: pd.DataFrame, cell_type: str) -> None:
    d = d.sort_values("RANK")
    x = d["RANK"].to_numpy()
    y = d["SCORE"].to_numpy()
    xp, yp, xn, yn = split_signed_curve(x, y)

    ax.axhline(0, color=COL_ZERO, lw=0.7, zorder=1)
    if len(xp):
        ax.fill_between(xp, 0, yp, color=COL_POS, alpha=0.18, linewidth=0, zorder=2)
        ax.plot(xp, yp, color=COL_POS, lw=1.15, solid_capstyle="round", zorder=3)
    if len(xn):
        ax.fill_between(xn, 0, yn, color=COL_NEG, alpha=0.18, linewidth=0, zorder=2)
        ax.plot(xn, yn, color=COL_NEG, lw=1.15, solid_capstyle="round", zorder=3)

    top = d.nsmallest(2, "RANK")
    bot = d.nlargest(2, "RANK")
    for _, r in top.iterrows():
        ax.annotate(
            str(r["SYMBOL"])[:10],
            xy=(r["RANK"], r["SCORE"]),
            xytext=(3, 2),
            textcoords="offset points",
            fontsize=6,
            color=COL_POS,
            ha="left",
            va="bottom",
        )
    for _, r in bot.iterrows():
        ax.annotate(
            str(r["SYMBOL"])[:10],
            xy=(r["RANK"], r["SCORE"]),
            xytext=(-3, -2),
            textcoords="offset points",
            fontsize=6,
            color=COL_NEG,
            ha="right",
            va="top",
        )

    ax.set_xlabel("Gene rank", fontsize=7.5)
    ax.set_ylabel(r"Score ($\delta \times -\log_{10}P$)", fontsize=7.5)
    # Cell-type identity as left spine label (no panel title)
    ax.text(
        -0.18,
        0.5,
        cell_type.replace("_", " "),
        transform=ax.transAxes,
        fontsize=9,
        fontweight="bold",
        va="center",
        ha="right",
        rotation=90,
    )
    ax.tick_params(labelsize=6.5, length=2.5, width=0.6)
    for spine in ax.spines.values():
        spine.set_color(COL_SPINE)
        spine.set_linewidth(0.6)
    ax.margins(x=0.02)


def draw_gsea_bubble(ax, gsea_ct: pd.DataFrame, fdr_cut: float) -> None:
    d = gsea_ct.copy()
    d["FDR"] = pd.to_numeric(d["FDR q-val"], errors="coerce")
    d["NES"] = pd.to_numeric(d["NES"], errors="coerce")
    d = d[np.isfinite(d["FDR"]) & np.isfinite(d["NES"])]
    d = d[d["FDR"] < fdr_cut].sort_values("FDR")
    # unique terms (prefer better FDR if duplicated across libraries)
    d = d.sort_values("FDR").drop_duplicates(subset=["Term"], keep="first")
    top_n = int(os.environ.get("GSEA_TOP_N_TERMS", "10"))
    d = d.head(top_n).copy()
    if d.empty:
        ax.axis("off")
        ax.text(0.5, 0.5, "No significant terms", ha="center", va="center", fontsize=8)
        return

    d["TERM_SHORT"] = d["Term"].map(short_term)
    d["NEGLOG10_FDR"] = -np.log10(np.clip(d["FDR"].to_numpy(), 1e-300, None))
    # order by NES for visual grouping
    d = d.sort_values("NES")
    d["Y"] = np.arange(len(d))

    # size by -log10 FDR
    s = 35 + 55 * (d["NEGLOG10_FDR"] / max(d["NEGLOG10_FDR"].max(), 1e-6))
    colors = np.where(d["NES"] >= 0, COL_POS, COL_NEG)

    ax.axvline(0, color=COL_ZERO, lw=0.7, zorder=1)
    ax.scatter(
        d["NES"],
        d["Y"],
        s=s,
        c=colors,
        alpha=0.85,
        edgecolors="white",
        linewidths=0.4,
        zorder=3,
    )
    ax.set_yticks(d["Y"])
    ax.set_yticklabels(d["TERM_SHORT"], fontsize=6.5)
    ax.set_xlabel("NES", fontsize=7.5)
    ax.tick_params(axis="x", labelsize=6.5, length=2.5, width=0.6)
    for spine in ax.spines.values():
        spine.set_color(COL_SPINE)
        spine.set_linewidth(0.6)
    ax.set_ylim(-0.6, len(d) - 0.4)
    # library tag as small text on right of points? optional: marker style by library
    # Use different edge for Reactome vs GO
    for _, r in d.iterrows():
        if r.get("LIBRARY") == "Reactome":
            ax.scatter(
                [r["NES"]],
                [r["Y"]],
                s=[35 + 55 * (r["NEGLOG10_FDR"] / max(d["NEGLOG10_FDR"].max(), 1e-6))],
                facecolors="none",
                edgecolors="#111111",
                linewidths=0.9,
                zorder=4,
            )


def main() -> int:
    root = project_root()
    os.chdir(root)
    method = os.environ.get("DECONV_METHOD", "bayesprism").lower()
    fdr_cut = float(os.environ.get("GSEA_FDR_CUTOFF", "0.05"))
    rows_per_page = int(os.environ.get("COMBINED_ROWS_PER_PAGE", "5"))
    method_label = {"bayesprism": "BayesPrism", "cibersort": "CIBERSORT"}.get(
        method, method
    )

    base = (
        root
        / "CLLCMA_alpha/results/deconvolution/cll_broad_2022"
        / "splicing_outlier_tme_association"
        / method
        / "rank_enrichment"
    )
    ranked = pd.read_csv(base / "gene_scores_ranked_by_celltype.tsv.gz", sep="\t")
    if "SYMBOL" not in ranked.columns:
        ranked["SYMBOL"] = ranked["GENE"]
    else:
        ranked["SYMBOL"] = ranked["SYMBOL"].fillna(ranked["GENE"])

    gsea = pd.read_csv(base / "gsea_prerank_by_celltype.tsv.gz", sep="\t")
    gsea["FDR q-val"] = pd.to_numeric(gsea["FDR q-val"], errors="coerce")
    gsea["NES"] = pd.to_numeric(gsea["NES"], errors="coerce")

    sig = gsea[gsea["FDR q-val"] < fdr_cut]
    cts_sig = set(sig["CELL_TYPE"].unique())
    cts = [c for c in CT_ORDER if c in cts_sig]
    cts += sorted(cts_sig - set(cts))
    if not cts:
        raise SystemExit(f"No cell types with GSEA FDR q < {fdr_cut}")

    print(
        f"Cell types with significant GSEA (FDR q < {fdr_cut}): {len(cts)} -> {cts}",
        flush=True,
    )

    out_dir = base / "combined_rank_gsea"
    out_dir.mkdir(parents=True, exist_ok=True)

    # chunk into pages
    pages = [
        cts[i : i + rows_per_page] for i in range(0, len(cts), rows_per_page)
    ]

    for page_i, page_cts in enumerate(pages, start=1):
        n_rows = len(page_cts)
        fig_w = 190 / 25.4
        # Taller rows for readability of GSEA term labels
        fig_h = max(100 / 25.4, 65 / 25.4 * n_rows)
        fig, axes = plt.subplots(
            n_rows,
            2,
            figsize=(fig_w, fig_h),
            layout="constrained",
            gridspec_kw={"width_ratios": [1.05, 1.35], "hspace": 0.22},
        )
        if n_rows == 1:
            axes = np.array([axes])

        for row, ct in enumerate(page_cts):
            draw_rank_panel(axes[row, 0], ranked[ranked["CELL_TYPE"] == ct], ct)
            draw_gsea_bubble(axes[row, 1], gsea[gsea["CELL_TYPE"] == ct], fdr_cut)

        legend_elems = [
            Patch(facecolor=COL_POS, edgecolor=COL_POS, alpha=0.4, label="Positive NES / score"),
            Patch(facecolor=COL_NEG, edgecolor=COL_NEG, alpha=0.4, label="Negative NES / score"),
            Line2D(
                [0],
                [0],
                marker="o",
                color="none",
                markerfacecolor="white",
                markeredgecolor="#111111",
                markersize=7,
                label="Reactome term (black ring)",
            ),
            Line2D([0], [0], color=COL_ZERO, lw=0.8, label="Zero line"),
        ]
        fig.legend(
            handles=legend_elems,
            loc="lower center",
            ncol=4,
            fontsize=7.5,
            frameon=False,
            bbox_to_anchor=(0.5, -0.01),
        )

        stem = f"rank_gsea_combined_fdr{fdr_cut:g}_page{page_i:02d}"
        pdf = out_dir / f"{stem}.pdf"
        png = out_dir / f"{stem}.png"
        fig.savefig(pdf, facecolor="white")
        fig.savefig(png, dpi=400, facecolor="white")
        plt.close(fig)
        print(f"Wrote {pdf}", flush=True)

    # Also one multipage-friendly list of included cell types
    summary = (
        sig.groupby("CELL_TYPE", as_index=False)
        .agg(N_SIG_TERMS=("Term", "nunique"), MIN_FDR=("FDR q-val", "min"))
        .sort_values("MIN_FDR")
    )
    summary.to_csv(out_dir / f"celltypes_with_gsea_fdr{fdr_cut:g}.tsv", sep="\t", index=False)

    (out_dir / f"rank_gsea_combined_fdr{fdr_cut:g}_provenance.txt").write_text(
        "\n".join(
            [
                "Provisional general figure (journal/phase not specified).",
                f"Method: {method_label}",
                f"GSEA FDR q cutoff: {fdr_cut}",
                f"Cell types included ({len(cts)}): {', '.join(cts)}",
                "Left panel: rank vs score (red=positive, blue=negative)",
                "Right panel: preranked GSEA bubble (x=NES, size=-log10 FDR q)",
                "Black ring = Reactome; filled only = GO BP",
                f"Source ranked: {base / 'gene_scores_ranked_by_celltype.tsv.gz'}",
                f"Source GSEA: {base / 'gsea_prerank_by_celltype.tsv.gz'}",
                f"Exported: {pd.Timestamp.now()}",
                "",
            ]
        )
    )
    print(f"Done -> {out_dir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
