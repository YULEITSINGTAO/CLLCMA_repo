#!/usr/bin/env python3
"""Preranked GSEA + combined figures using CellMarker gene sets.

Builds blood/immune Normal CellMarker gene sets from
CLLCMA_alpha/references/human_cell_marker.zip, runs preranked GSEA on
existing splicing-outlier rank scores, and writes combined
rank-score | GSEA bubble figures (same layout as GO/Reactome combined).

Env:
  CLL_TME_DIR, DECONV_METHOD
  GSEA_FDR_CUTOFF (default 0.05)
  GSEA_TOP_N_TERMS (default 10)
  COMBINED_ROWS_PER_PAGE (default 5)
  CELLMARKER_MAX_SET (default 2000)
  CELLMARKER_MIN_SET (default 10)
"""

from __future__ import annotations

import os
import zipfile
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# Reuse plotting helpers from combined GO/Reactome script
import importlib.util


COL_POS = "#B2182B"
COL_NEG = "#2166AC"
COL_ZERO = "#666666"
COL_SPINE = "#333333"

KEEP_TISSUE = {
    "Blood",
    "Bone marrow",
    "Lymph node",
    "Spleen",
    "Tonsil",
    "Mesenteric lymph node",
    "Peripheral blood",
    "Umbilical cord blood",
}

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


def load_plot_helpers():
    root = project_root()
    path = root / "CLLCMA_alpha/scripts/plot_rank_gsea_combined.py"
    spec = importlib.util.spec_from_file_location("rank_gsea_combined", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def build_cellmarker_gmt(
    zip_path: Path,
    gmt_path: Path,
    summary_path: Path,
    min_set: int,
    max_set: int,
) -> dict[str, list[str]]:
    """Return gene_sets dict and write GMT for blood/immune Normal markers."""
    sets: dict[str, set[str]] = defaultdict(set)
    with zipfile.ZipFile(zip_path) as zf:
        names = zf.namelist()
        member = names[0]
        with zf.open(member) as fh:
            for ch in pd.read_csv(
                fh,
                sep="\t",
                usecols=["tissue_class", "disease", "cell_name", "symbol"],
                chunksize=250_000,
                low_memory=False,
            ):
                m = (
                    (ch["disease"].astype(str) == "Normal")
                    & (ch["tissue_class"].isin(KEEP_TISSUE))
                    & ch["symbol"].notna()
                )
                sub = ch.loc[m]
                for cell, sym in zip(sub["cell_name"], sub["symbol"]):
                    s = str(sym).strip().upper()
                    if s and s != "NAN":
                        sets[str(cell)].add(s)

    gene_sets: dict[str, list[str]] = {}
    rows = []
    for cell, genes in sorted(sets.items()):
        n = len(genes)
        use = min_set <= n <= max_set
        rows.append({"CELL_NAME": cell, "N_GENES": n, "USED": use})
        if use:
            gene_sets[cell] = sorted(genes)

    gmt_path.parent.mkdir(parents=True, exist_ok=True)
    with gmt_path.open("w") as out:
        for name, genes in gene_sets.items():
            out.write(name + "\tna\t" + "\t".join(genes) + "\n")
    summary = pd.DataFrame(rows).sort_values("N_GENES", ascending=False)
    summary.to_csv(summary_path, sep="\t", index=False)
    print(
        f"CellMarker GMT: {len(gene_sets)} sets "
        f"(min={min_set}, max={max_set}) -> {gmt_path}",
        flush=True,
    )
    return gene_sets


def main() -> int:
    root = project_root()
    os.chdir(root)
    method = os.environ.get("DECONV_METHOD", "bayesprism").lower()
    fdr_cut = float(os.environ.get("GSEA_FDR_CUTOFF", "0.05"))
    top_n = int(os.environ.get("GSEA_TOP_N_TERMS", "10"))
    rows_per_page = int(os.environ.get("COMBINED_ROWS_PER_PAGE", "5"))
    min_set = int(os.environ.get("CELLMARKER_MIN_SET", "10"))
    max_set = int(os.environ.get("CELLMARKER_MAX_SET", "2000"))
    method_label = {"bayesprism": "BayesPrism", "cibersort": "CIBERSORT"}.get(
        method, method
    )
    os.environ["GSEA_TOP_N_TERMS"] = str(top_n)

    zip_path = root / "CLLCMA_alpha/references/human_cell_marker.zip"
    if not zip_path.exists():
        raise SystemExit(f"Missing {zip_path}")

    base = (
        root
        / "CLLCMA_alpha/results/deconvolution/cll_broad_2022"
        / "splicing_outlier_tme_association"
        / method
        / "rank_enrichment"
    )
    out_dir = base / "cell_marker_gsea"
    out_dir.mkdir(parents=True, exist_ok=True)
    combined_dir = out_dir / "combined_rank_gsea"
    combined_dir.mkdir(parents=True, exist_ok=True)

    gmt_path = out_dir / "CellMarker_blood_immune_Normal.gmt"
    summary_path = out_dir / "CellMarker_gene_set_sizes.tsv"
    gene_sets = build_cellmarker_gmt(zip_path, gmt_path, summary_path, min_set, max_set)

    ranked = pd.read_csv(base / "gene_scores_ranked_by_celltype.tsv.gz", sep="\t")
    if "SYMBOL" not in ranked.columns:
        ranked["SYMBOL"] = ranked["GENE"]
    else:
        ranked["SYMBOL"] = ranked["SYMBOL"].fillna(ranked["GENE"])

    import gseapy as gp

    gsea_rows = []
    for ct, sub in ranked.groupby("CELL_TYPE"):
        sub_sym = sub.dropna(subset=["SYMBOL"]).copy()
        sub_sym["SYMBOL"] = sub_sym["SYMBOL"].astype(str).str.upper()
        sub_sym = (
            sub_sym.assign(ABS=sub_sym["SCORE"].abs())
            .sort_values("ABS", ascending=False)
            .drop_duplicates("SYMBOL")
        )
        rnk = sub_sym.set_index("SYMBOL")["SCORE"].sort_values(ascending=False)
        try:
            pre = gp.prerank(
                rnk=rnk,
                gene_sets=gene_sets,
                min_size=min_set,
                max_size=max_set,
                permutation_num=200,
                seed=1,
                threads=4,
                outdir=None,
                verbose=False,
            )
            res = pre.res2d.copy()
            res.insert(0, "CELL_TYPE", ct)
            res.insert(1, "LIBRARY", "CellMarker_blood_immune")
            gsea_rows.append(res)
            print(f"  CellMarker GSEA {ct}: {len(res)} sets", flush=True)
        except Exception as e:
            print(f"  CellMarker GSEA failed {ct}: {e}", flush=True)

    if not gsea_rows:
        raise SystemExit("No CellMarker GSEA results")

    gsea = pd.concat(gsea_rows, ignore_index=True)
    gsea.to_csv(out_dir / "gsea_prerank_by_celltype.tsv.gz", sep="\t", index=False, compression="gzip")
    gsea["FDR q-val"] = pd.to_numeric(gsea["FDR q-val"], errors="coerce")
    gsea["NES"] = pd.to_numeric(gsea["NES"], errors="coerce")
    gsea_sig = gsea[gsea["FDR q-val"] < fdr_cut].copy()
    gsea_sig.to_csv(
        out_dir / f"gsea_prerank_fdr{fdr_cut:g}.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    print(
        f"Significant terms FDR q < {fdr_cut}: {len(gsea_sig)} "
        f"across {gsea_sig['CELL_TYPE'].nunique()} cell types",
        flush=True,
    )

    # Combined figures
    helpers = load_plot_helpers()
    draw_rank_panel = helpers.draw_rank_panel
    draw_gsea_bubble = helpers.draw_gsea_bubble

    cts_sig = set(gsea_sig["CELL_TYPE"].unique())
    cts = [c for c in CT_ORDER if c in cts_sig] + sorted(cts_sig - set(CT_ORDER))
    if not cts:
        print("No cell types with significant CellMarker GSEA; skip combined plots")
        return 0

    summary_ct = (
        gsea_sig.groupby("CELL_TYPE", as_index=False)
        .agg(N_SIG_TERMS=("Term", "nunique"), MIN_FDR=("FDR q-val", "min"))
        .sort_values("MIN_FDR")
    )
    summary_ct.to_csv(
        combined_dir / f"celltypes_with_gsea_fdr{fdr_cut:g}.tsv", sep="\t", index=False
    )

    pages = [cts[i : i + rows_per_page] for i in range(0, len(cts), rows_per_page)]
    for page_i, page_cts in enumerate(pages, start=1):
        n_rows = len(page_cts)
        fig_w = float(os.environ.get("COMBINED_FIG_WIDTH_MM", "260")) / 25.4
        fig_h = max(100 / 25.4, 65 / 25.4 * n_rows)
        fig, axes = plt.subplots(
            n_rows,
            2,
            figsize=(fig_w, fig_h),
            layout="constrained",
            gridspec_kw={"width_ratios": [1.0, 1.55], "hspace": 0.22},
        )
        if n_rows == 1:
            axes = np.asarray([axes])

        for row, ct in enumerate(page_cts):
            draw_rank_panel(axes[row, 0], ranked[ranked["CELL_TYPE"] == ct], ct)
            draw_gsea_bubble(axes[row, 1], gsea[gsea["CELL_TYPE"] == ct], fdr_cut)

        legend_elems = [
            Patch(facecolor=COL_POS, edgecolor=COL_POS, alpha=0.4, label="Positive NES / score"),
            Patch(facecolor=COL_NEG, edgecolor=COL_NEG, alpha=0.4, label="Negative NES / score"),
            Line2D([0], [0], color=COL_ZERO, lw=0.8, label="Zero line"),
        ]
        fig.legend(
            handles=legend_elems,
            loc="lower center",
            ncol=3,
            fontsize=7.5,
            frameon=False,
            bbox_to_anchor=(0.5, -0.01),
        )

        stem = f"rank_cellmarker_gsea_fdr{fdr_cut:g}_page{page_i:02d}"
        fig.savefig(combined_dir / f"{stem}.pdf", facecolor="white")
        fig.savefig(combined_dir / f"{stem}.png", dpi=400, facecolor="white")
        plt.close(fig)
        print(f"Wrote {combined_dir / stem}.pdf", flush=True)

    (combined_dir / f"rank_cellmarker_gsea_fdr{fdr_cut:g}_provenance.txt").write_text(
        "\n".join(
            [
                "Provisional general figure (journal/phase not specified).",
                f"Method: {method_label}",
                "Gene sets: CellMarker (human_cell_marker.zip)",
                "Filter: Human, disease=Normal, tissue in Blood/Bone marrow/Lymph node/Spleen/Tonsil/…",
                f"Set size filter: {min_set}-{max_set} genes (larger sets truncated to max)",
                f"GSEA FDR q cutoff: {fdr_cut}",
                f"Cell types included ({len(cts)}): {', '.join(cts)}",
                f"GMT: {gmt_path}",
                f"Exported: {pd.Timestamp.now()}",
                "",
            ]
        )
    )

    # KEY findings
    lines = [
        "CELLMARKER GSEA ON SPLICING-OUTLIER RANKS",
        f"Method: {method_label}",
        f"Sets: {len(gene_sets)} blood/immune Normal CellMarker cell types",
        f"FDR q < {fdr_cut}: {len(gsea_sig)} terms / {len(cts)} TME cell types",
        "=== Top terms ===",
    ]
    for _, r in gsea.nsmallest(20, "FDR q-val").iterrows():
        lines.append(
            f"{r['CELL_TYPE']} | {r['Term']} | NES={r['NES']:+.3f} FDR={r['FDR q-val']:.3g}"
        )
    (out_dir / "KEY_FINDINGS.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[:25]), flush=True)
    print(f"Done -> {out_dir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
