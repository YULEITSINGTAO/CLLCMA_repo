#!/usr/bin/env python3
"""Preranked GO:BP or Reactome GSEA across B -> MBL -> CLL -> RT."""

from __future__ import annotations

import argparse
import csv
import math
import sys
import textwrap
from collections import defaultdict
from datetime import date
from pathlib import Path

import os

import gseapy as gp
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure_io import export_figure, style_context  # noqa: E402

PROJECT = Path(os.environ.get("CLL_TME_DIR", Path.home() / "research" / "CLL_TME")).resolve()
INPUT_DIR = PROJECT / "CLLCMA_alpha/results/bcell_stage_deseq2"
OUTPUT_DIR = INPUT_DIR / "publication_figures/stage_progression"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

DATABASES = {
    "go": {
        "library": "GO_Biological_Process_2026",
        "prefix": "GO_BP",
        "term_axis": "GO biological-process term",
        "figure_title": "Rank-based GO biological-process enrichment across B-cell disease progression",
        "shared_title": "Top shared GO biological processes across disease-stage transitions",
        "caption_name": "GO biological-process",
    },
    "reactome": {
        "library": "Reactome_Pathways_2024",
        "prefix": "Reactome",
        "term_axis": "Reactome pathway",
        "figure_title": "Rank-based Reactome pathway enrichment\nacross B-cell disease progression",
        "shared_title": "Top shared Reactome pathways across disease-stage transitions",
        "caption_name": "Reactome pathway",
    },
}
SEED = 42
PERMUTATIONS = 1000
MIN_SIZE = 15
MAX_SIZE = 500
FDR_CUTOFF = 0.05
MAX_TERMS = 10
MAX_KEY_GENES = 4

CONTRASTS = [
    ("MBL_vs_B", "B \u2192 MBL"),
    ("CLL_vs_MBL", "MBL \u2192 CLL"),
    ("RT_vs_CLL", "CLL \u2192 RT"),
]


def write_gmt(gene_sets: dict[str, list[str]], path: Path, library: str) -> None:
    """Freeze the retrieved gene-set definitions used by this analysis."""
    with path.open("w") as handle:
        for term in sorted(gene_sets):
            genes = sorted(set(map(str.upper, gene_sets[term])))
            handle.write("\t".join([term, library, *genes]) + "\n")


def ranked_genes(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t")
    table = table.loc[
        table["gene"].notna() & np.isfinite(table["stat"]),
        ["gene", "stat"],
    ].copy()
    table["gene"] = table["gene"].astype(str).str.upper()
    # Deterministic tie-breaking by symbol keeps the preranked order reproducible.
    table = table.sort_values(
        ["stat", "gene"],
        ascending=[False, True],
        kind="mergesort",
    )
    table = table.drop_duplicates("gene", keep="first")
    return table


def normalize_result(result: pd.DataFrame, stem: str, transition: str) -> pd.DataFrame:
    out = result.copy()
    out.columns = [str(c).strip() for c in out.columns]
    aliases = {
        "Name": "name",
        "Term": "term",
        "ES": "es",
        "NES": "nes",
        "NOM p-val": "nominal_p",
        "FDR q-val": "fdr",
        "FWER p-val": "fwer",
        "Tag %": "tag_percent",
        "Gene %": "gene_percent",
        "Lead_genes": "leading_edge",
    }
    out = out.rename(columns=aliases)
    for col in ("es", "nes", "nominal_p", "fdr", "fwer"):
        if col in out:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    out.insert(0, "contrast", stem)
    out.insert(1, "transition", transition)
    return out


def leading_set(value: object) -> set[str]:
    if pd.isna(value):
        return set()
    return {x.strip().upper() for x in str(value).split(";") if x.strip()}


def key_genes_for_term(significant: pd.DataFrame) -> str:
    """Rank leading-edge genes by transition recurrence then mean edge rank."""
    gene_ranks: dict[str, list[int]] = defaultdict(list)
    for value in significant["leading_edge"]:
        if pd.isna(value):
            continue
        genes = [x.strip().upper() for x in str(value).split(";") if x.strip()]
        for rank, gene in enumerate(genes, start=1):
            gene_ranks[gene].append(rank)
    ranked = sorted(
        gene_ranks,
        key=lambda gene: (
            -len(gene_ranks[gene]),
            float(np.mean(gene_ranks[gene])),
            gene,
        ),
    )
    return ";".join(ranked[:MAX_KEY_GENES])


def select_shared_terms(all_results: pd.DataFrame) -> pd.DataFrame:
    """Select shared terms and collapse highly overlapping leading edges."""
    records: list[dict[str, object]] = []
    for term, group in all_results.groupby("term", sort=False):
        significant = group.loc[group["fdr"] < FDR_CUTOFF]
        if len(significant) < 2:
            continue
        records.append(
            {
                "term": term,
                "n_significant_transitions": len(significant),
                "min_fdr": float(significant["fdr"].min()),
                "mean_abs_nes_significant": float(significant["nes"].abs().mean()),
                "max_abs_nes": float(group["nes"].abs().max()),
                "key_leading_edge_genes": key_genes_for_term(significant),
            }
        )
    metrics = pd.DataFrame(records)
    if metrics.empty:
        raise RuntimeError("No terms are significant in at least two transitions")
    metrics = (
        metrics.sort_values(
            ["n_significant_transitions", "min_fdr", "mean_abs_nes_significant"],
            ascending=[False, True, False],
            kind="mergesort",
        )
        .reset_index(drop=True)
    )

    edge_by_term: dict[str, set[str]] = defaultdict(set)
    significant_rows = all_results.loc[all_results["fdr"] < FDR_CUTOFF]
    for row in significant_rows.itertuples():
        edge_by_term[row.term].update(leading_set(row.leading_edge))

    selected_indices: list[int] = []
    selected_terms: list[str] = []
    for index, row in metrics.iterrows():
        genes = edge_by_term[row["term"]]
        redundant = False
        for kept in selected_terms:
            other = edge_by_term[kept]
            union = genes | other
            if union and len(genes & other) / len(union) >= 0.5:
                redundant = True
                break
        if not redundant:
            selected_indices.append(index)
            selected_terms.append(row["term"])
        if len(selected_indices) >= MAX_TERMS:
            break
    return metrics.loc[selected_indices].reset_index(drop=True)


def clean_term(term: str) -> str:
    label = term.replace("_", " ")
    if label.endswith(")") and "(GO:" in label:
        label = label[: label.rfind("(GO:")].strip()
    return "\n".join(textwrap.wrap(label.title(), width=46))


def make_trajectory_figure(
    all_results: pd.DataFrame,
    selected_metrics: pd.DataFrame,
    spec: dict[str, str],
) -> None:
    selected = selected_metrics["term"].tolist()
    gene_labels = dict(
        zip(
            selected_metrics["term"],
            selected_metrics["key_leading_edge_genes"],
            strict=True,
        )
    )
    transition_order = [label for _, label in CONTRASTS]
    block = all_results.loc[all_results["term"].isin(selected)].copy()
    lookup = {
        (row.term, row.transition): row
        for row in block.itertuples()
    }

    # Group terms by the transition with their strongest absolute NES.
    def term_sort_key(term: str) -> tuple[int, float]:
        vals = [
            abs(float(lookup[(term, tr)].nes))
            if (term, tr) in lookup and np.isfinite(lookup[(term, tr)].nes)
            else -1
            for tr in transition_order
        ]
        peak = int(np.argmax(vals))
        signed = float(lookup[(term, transition_order[peak])].nes)
        return peak, -signed

    ordered_terms = sorted(selected, key=term_sort_key)
    max_abs_nes = max(2.0, float(np.nanmax(np.abs(block["nes"]))))
    norm = mpl.colors.TwoSlopeNorm(vmin=-max_abs_nes, vcenter=0, vmax=max_abs_nes)
    cmap = mpl.colormaps["RdBu_r"]

    with style_context("default", palette_name="okabe_ito_on_white"):
        height = max(4.8, 0.36 * len(ordered_terms) + 1.6)
        fig, ax = plt.subplots(
            figsize=(180 / 25.4, height),
            layout="constrained",
        )

        for y, term in enumerate(ordered_terms):
            for x, transition in enumerate(transition_order):
                row = lookup.get((term, transition))
                if row is None or not np.isfinite(row.nes):
                    continue
                fdr = max(float(row.fdr), np.finfo(float).tiny)
                strength = min(-math.log10(fdr), 10.0)
                size = 14 + 13 * strength
                significant = fdr < FDR_CUTOFF
                ax.scatter(
                    x,
                    y,
                    s=size,
                    c=[cmap(norm(float(row.nes)))],
                    marker="o" if significant else "x",
                    linewidths=0.8 if significant else 1.0,
                    edgecolors="black" if significant else None,
                    zorder=3,
                )

        ax.set_xticks(range(len(transition_order)), transition_order)
        ax.set_yticks(
            range(len(ordered_terms)),
            [
                f"{clean_term(term)}\nLeading edge: {gene_labels[term].replace(';', ', ')}"
                for term in ordered_terms
            ],
        )
        ax.invert_yaxis()
        ax.set_xlabel("Ordered disease-stage transition (not elapsed time)")
        ax.set_ylabel(spec["term_axis"])
        ax.set_title(
            spec["shared_title"],
            loc="left",
            fontweight="bold",
        )
        ax.grid(axis="both", color="0.90", linewidth=0.5)
        ax.set_axisbelow(True)

        colorbar = fig.colorbar(
            mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
            ax=ax,
            fraction=0.035,
            pad=0.02,
        )
        colorbar.set_label("Normalized enrichment score (NES)\nred: later-stage enriched; blue: earlier-stage enriched")

        size_handles = []
        for fdr in (0.05, 0.01, 0.001):
            strength = min(-math.log10(fdr), 10.0)
            size_handles.append(
                ax.scatter(
                    [], [], s=14 + 13 * strength, facecolor="white",
                    edgecolor="black", linewidth=0.8, label=f"FDR={fdr:g}"
                )
            )
        nonsig = ax.scatter(
            [], [], s=35, marker="x", color="black", linewidth=1.0,
            label="FDR\u22650.05"
        )
        ax.legend(
            handles=[*size_handles, nonsig],
            title="Point size: -log10(FDR)",
            loc="lower left",
            bbox_to_anchor=(1.01, 0.02),
            frameon=False,
        )

        export_figure(
            fig,
            OUTPUT_DIR / f"{spec['prefix']}_transition_GSEA_trajectory",
            formats=["pdf", "png"],
            dpi=400,
            bbox_inches=None,
            overwrite=True,
            write_manifest=True,
            provenance={
                "inputs": [str(INPUT_DIR / f"deseq2_{stem}.tsv") for stem, _ in CONTRASTS],
                "ranking": "DESeq2 Wald statistic; all finite tested genes; duplicate uppercase symbols keep highest statistic",
                "method": "GSEApy prerank",
                "gene_sets": spec["library"],
                "gene_set_retrieval_date": str(date.today()),
                "parameters": {
                    "permutations": PERMUTATIONS,
                    "seed": SEED,
                    "min_size": MIN_SIZE,
                    "max_size": MAX_SIZE,
                    "fdr_cutoff": FDR_CUTOFF,
                },
                "term_selection": "Up to 10 terms with FDR<0.05 in at least two transitions; ranked by recurrence, minimum FDR, then mean absolute NES; leading-edge Jaccard >=0.5 collapsed",
                "key_gene_selection": "Up to four leading-edge genes per term; ranked by recurrence across significant transitions, mean leading-edge rank, then symbol",
                "interpretation": "Columns are separate stage contrasts, not longitudinal elapsed-time measurements",
                "publisher_status": "provisional general figure; target journal not specified",
            },
        )
        plt.close(fig)


def main(database: str, plot_only: bool = False) -> None:
    spec = DATABASES[database]
    library = spec["library"]
    prefix = spec["prefix"]
    combined_file = OUTPUT_DIR / f"{prefix}_GSEA_all_transitions.tsv"
    if plot_only:
        combined = pd.read_csv(combined_file, sep="\t")
    else:
        gene_sets = gp.get_library(name=library, organism="Human")
        write_gmt(
            gene_sets,
            OUTPUT_DIR / f"{library}_retrieved_{date.today()}.gmt",
            library,
        )

        results: list[pd.DataFrame] = []
        rank_summary: list[dict[str, object]] = []
        for stem, transition in CONTRASTS:
            rank = ranked_genes(INPUT_DIR / f"deseq2_{stem}.tsv")
            rank.to_csv(OUTPUT_DIR / f"rank_{stem}.tsv", sep="\t", index=False)
            pre = gp.prerank(
                rnk=rank,
                gene_sets=gene_sets,
                min_size=MIN_SIZE,
                max_size=MAX_SIZE,
                permutation_num=PERMUTATIONS,
                seed=SEED,
                threads=8,
                outdir=None,
                verbose=True,
            )
            table = normalize_result(pre.res2d, stem, transition)
            table.to_csv(
                OUTPUT_DIR / f"{prefix}_GSEA_{stem}.tsv",
                sep="\t",
                index=False,
            )
            results.append(table)
            rank_summary.append(
                {
                    "contrast": stem,
                    "transition": transition,
                    "n_ranked_genes": len(rank),
                    "n_tested_terms": len(table),
                    "n_fdr_lt_0.05": int((table["fdr"] < FDR_CUTOFF).sum()),
                    "n_positive_fdr_lt_0.05": int(
                        ((table["fdr"] < FDR_CUTOFF) & (table["nes"] > 0)).sum()
                    ),
                    "n_negative_fdr_lt_0.05": int(
                        ((table["fdr"] < FDR_CUTOFF) & (table["nes"] < 0)).sum()
                    ),
                }
            )

        combined = pd.concat(results, ignore_index=True)
        combined.to_csv(combined_file, sep="\t", index=False)
        pd.DataFrame(rank_summary).to_csv(
            OUTPUT_DIR / f"{prefix}_GSEA_summary.tsv",
            sep="\t",
            index=False,
        )

    selected_metrics = select_shared_terms(combined)
    selected_metrics.to_csv(
        OUTPUT_DIR / f"{prefix}_GSEA_representative_terms.tsv",
        sep="\t",
        index=False,
    )
    make_trajectory_figure(combined, selected_metrics, spec)

    (OUTPUT_DIR / f"{prefix}_transition_GSEA_caption.txt").write_text(
        "\n".join(
            [
                f"Figure Y. Rank-based {spec['caption_name']} enrichment across ordered B-cell disease transitions.",
                "Preranked GSEA used the complete DESeq2 Wald-statistic ranking for MBL versus B, CLL versus MBL, and RT versus CLL.",
                "Positive NES (red) denotes enrichment toward genes increased in the later stage; negative NES (blue) denotes enrichment toward the earlier stage.",
                "Point size encodes -log10(FDR), capped at 10; circles have FDR<0.05 and crosses have FDR\u22650.05.",
                "The figure shows up to 10 terms significant at FDR<0.05 in at least two transitions, ranked by recurrence, minimum FDR, and mean absolute NES; terms with leading-edge Jaccard overlap \u22650.5 were collapsed.",
                "Up to four key leading-edge genes are shown per term, ranked by recurrence across significant transitions and then mean leading-edge rank.",
                f"Gene sets: {library}, retrieved {date.today()}; sizes {MIN_SIZE}-{MAX_SIZE}; {PERMUTATIONS} permutations; seed {SEED}.",
                "The columns are separate disease-stage contrasts and must not be interpreted as repeated measurements at uniform elapsed times.",
                "MBL versus B is partly cohort-confounded. Results require biological validation and are not evidence of temporal causality.",
            ]
        )
        + "\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--database",
        choices=sorted(DATABASES),
        default="go",
        help="Gene-set database to analyze (default: go)",
    )
    parser.add_argument(
        "--plot-only",
        action="store_true",
        help="Regenerate selection and figure from existing combined GSEA table",
    )
    args = parser.parse_args()
    main(args.database, plot_only=args.plot_only)
