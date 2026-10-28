#!/usr/bin/env python3
"""Summarize inferCNV as donor-level chromosome means.

Publisher/journal, article type, and submission phase are not specified;
this is a provisional general figure. Visual choices are not journal-compliant
by default and must be re-checked against the target journal.
"""

from __future__ import annotations

import csv
import sys
from collections import defaultdict
from pathlib import Path

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure_io import OKABE_ITO, export_figure, style_context  # noqa: E402

PROJECT = Path(os.environ.get("CLL_TME_DIR", Path.home() / "research" / "CLL_TME")).resolve()
DATA_DIR = (
    PROJECT
    / "CLLCMA_alpha/results/bcell_stage_deseq2/publication_figures/"
    / "infercnv_bcell_stage/chromosome_summary"
)
OUT_BASE = DATA_DIR / "infercnv_chromosome_donor_summary"

STAGE_ORDER = ("Normal_B", "MBL", "CLL", "RT")
OBS_STAGES = ("MBL", "CLL", "RT")
CHR_LEVELS = tuple(f"chr{i}" for i in range(1, 23))
STAGE_COLORS = {
    "Normal_B": OKABE_ITO["blue"],
    "MBL": OKABE_ITO["bluish_green"],
    "CLL": OKABE_ITO["orange"],
    "RT": OKABE_ITO["vermillion"],
}
STAGE_MARKERS = {
    "Normal_B": "o",
    "MBL": "s",
    "CLL": "^",
    "RT": "D",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def donor_matrix(rows: list[dict[str, str]]) -> tuple[np.ndarray, list[str], list[str]]:
    values: dict[tuple[str, str, str], float] = {}
    cells: dict[tuple[str, str], int] = {}
    for row in rows:
        stage = row["stage"]
        donor = row["donor"]
        chr_name = row["chr"]
        values[(stage, donor, chr_name)] = float(row["mean_log2_dev"])
        cells[(stage, donor)] = int(row["n_cells"])

    ordered_donors: list[tuple[str, str]] = []
    for stage in OBS_STAGES:
        donors = sorted({donor for (st, donor) in cells if st == stage})
        chr12 = {
            donor: values.get((stage, donor, "chr12"), 0.0) for donor in donors
        }
        donors = sorted(donors, key=lambda d: chr12[d], reverse=True)
        ordered_donors.extend((stage, donor) for donor in donors)

    mat = np.full((len(ordered_donors), len(CHR_LEVELS)), np.nan)
    labels = []
    stages = []
    for i, (stage, donor) in enumerate(ordered_donors):
        labels.append(f"{donor}")
        stages.append(stage)
        for j, chr_name in enumerate(CHR_LEVELS):
            mat[i, j] = values[(stage, donor, chr_name)]
    return mat, labels, stages


def stage_summary_matrix(rows: list[dict[str, str]], field: str) -> np.ndarray:
    lookup = {(r["stage"], r["chr"]): float(r[field]) for r in rows}
    mat = np.full((len(STAGE_ORDER), len(CHR_LEVELS)), np.nan)
    for i, stage in enumerate(STAGE_ORDER):
        for j, chr_name in enumerate(CHR_LEVELS):
            mat[i, j] = lookup[(stage, chr_name)]
    return mat


def write_key_findings(
    path: Path,
    stage_rows: list[dict[str, str]],
    rank_rows: list[dict[str, str]],
    donor_n: list[dict[str, str]],
    abs_limit: float,
) -> None:
    n_map = {r["stage"]: r for r in donor_n}
    lines = [
        "inferCNV chromosome summary (RNA-based; HMM was not run).",
        "Unit of replication: donor. Estimator: chromosome-mean log2(inferCNV expr.data),",
        "then the mean of those cell values within each donor. Donors with <5 cells excluded.",
        "Normal_B is the inferCNV reference group. Values are relative expression deviations,",
        "not genomic copy-number calls.",
        "",
        "Donors retained:",
    ]
    for stage in STAGE_ORDER:
        rec = n_map[stage]
        lines.append(
            f"  {stage}: n_donors={rec['n_donors']} "
            f"(median {rec['median_cells']} cells, range {rec['min_cells']}-{rec['max_cells']})"
        )
    lines += [
        "",
        f"Heatmap color limits: log2 deviation { -abs_limit:.2f} to {abs_limit:.2f} "
        "(symmetric about 0; 0 = diploid-centered).",
        "",
        "Largest systematic shifts vs Normal_B (mean of observation-donor chromosome means):",
    ]
    by_shift = sorted(
        rank_rows, key=lambda r: abs(float(r["mean_minus_normal_b"])), reverse=True
    )
    for rec in by_shift[:8]:
        lines.append(
            f"  {rec['chr']}: mean-Normal_B={float(rec['mean_minus_normal_b']):+.4f} "
            f"(observation SD={float(rec['sd_observation_donors']):.4f}, "
            f"MAD={float(rec['mad_observation_donors']):.4f})"
        )
    lines += [
        "",
        "Largest between-donor heterogeneity among MBL+CLL+RT (SD of donor means):",
    ]
    by_sd = sorted(
        rank_rows, key=lambda r: float(r["sd_observation_donors"]), reverse=True
    )
    for rec in by_sd[:8]:
        lines.append(
            f"  {rec['chr']}: SD={float(rec['sd_observation_donors']):.4f} "
            f"mean-Normal_B={float(rec['mean_minus_normal_b']):+.4f}"
        )
    lines += [
        "",
        "Stage-specific chromosome means minus Normal_B (selected chromosomes):",
    ]
    interesting = ["chr12", "chr13", "chr8", "chr4", "chr22", "chr19", "chr11", "chr17"]
    lookup = {(r["stage"], r["chr"]): r for r in stage_rows}
    for chr_name in interesting:
        bits = []
        for stage in OBS_STAGES:
            rec = lookup[(stage, chr_name)]
            bits.append(
                f"{stage} {float(rec['mean_minus_normal_b']):+.4f} "
                f"(SD {float(rec['sd_of_donor_means']):.4f}, n={rec['n_donors']})"
            )
        lines.append(f"  {chr_name}: " + "; ".join(bits))
    lines += [
        "",
        "Notes:",
        "- chr12 gain is subset-driven (high SD relative to MAD): compatible with trisomy 12 in a minority of donors.",
        "- Whole-chromosome means dilute focal events (e.g. 13q, 17p, 11q). See donor_chromosome_focal_bin.tsv.",
        "- chr19/chr22 are gene-dense; RNA-based inferCNV can produce non-copy-number expression shifts.",
        "- Matched genomic CNA data should be used before interpreting these as cytogenetic lesions.",
    ]
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    donor_rows = read_tsv(DATA_DIR / "donor_chromosome_mean_log2_deviation.tsv")
    stage_rows = read_tsv(DATA_DIR / "chromosome_stage_summary.tsv")
    rank_rows = read_tsv(DATA_DIR / "chromosome_ranked_by_donor_mad.tsv")
    donor_n = read_tsv(DATA_DIR / "donor_counts_by_stage.tsv")
    gene_n = {r["chr"]: int(r["N"]) for r in read_tsv(DATA_DIR / "chromosome_gene_counts.tsv")}

    obs_mat, donor_labels, _donor_stages = donor_matrix(donor_rows)
    diff_mat = stage_summary_matrix(stage_rows, "mean_minus_normal_b")
    abs_limit = 0.10

    n_by_stage = {r["stage"]: int(r["n_donors"]) for r in donor_n}
    write_key_findings(
        DATA_DIR / "KEY_FINDINGS.txt",
        stage_rows,
        rank_rows,
        donor_n,
        abs_limit,
    )

    chr_labels = [c.replace("chr", "") for chr_name in CHR_LEVELS for c in [chr_name]]
    cmap = mpl.colormaps["RdBu_r"].with_extremes(bad="#777777")
    norm = mpl.colors.TwoSlopeNorm(vmin=-abs_limit, vcenter=0.0, vmax=abs_limit)

    with style_context("default", palette_name="okabe_ito_on_white"):
        fig = plt.figure(figsize=(180 / 25.4, 215 / 25.4), layout="constrained")
        fig.set_facecolor("white")
        mosaic = fig.subplot_mosaic(
            [
                ["stage", "heat", "cbar"],
                ["sdlab", "sd", "cbar_blank_sd"],
                ["difflab", "diff", "cbar_blank_diff"],
            ],
            width_ratios=[0.001, 1.0, 0.025],
            height_ratios=[2.6, 0.85, 1.05],
        )
        ax_stage = mosaic["stage"]
        ax_heat = mosaic["heat"]
        ax_cbar = mosaic["cbar"]
        ax_sd = mosaic["sd"]
        ax_diff = mosaic["diff"]
        mosaic["sdlab"].axis("off")
        mosaic["difflab"].axis("off")
        mosaic["cbar_blank_sd"].axis("off")
        mosaic["cbar_blank_diff"].axis("off")
        ax_stage.axis("off")

        im = ax_heat.imshow(
            obs_mat,
            aspect="auto",
            cmap=cmap,
            norm=norm,
            interpolation="nearest",
        )
        ax_heat.set_xticks(np.arange(len(CHR_LEVELS)))
        ax_heat.set_xticklabels(chr_labels)
        ax_heat.set_yticks([])
        ax_heat.set_xlabel("Chromosome")
        ax_heat.set_title(
            "Donor-mean inferCNV log2 deviation by chromosome "
            f"(MBL n={n_by_stage['MBL']}, CLL n={n_by_stage['CLL']}, RT n={n_by_stage['RT']})",
            loc="left",
        )
        boundaries = np.cumsum(
            [n_by_stage["MBL"], n_by_stage["CLL"], n_by_stage["RT"]]
        )[:-1]
        for y in boundaries:
            ax_heat.axhline(y - 0.5, color="black", linewidth=0.6)
        ax_heat.tick_params(axis="x", length=0)
        # Rows ordered within stage by chr12 mean (descending), disclosed in caption.

        sd_lookup = {r["chr"]: float(r["sd_observation_donors"]) for r in rank_rows}
        sd_vals = np.array([sd_lookup[c] for c in CHR_LEVELS])
        ax_sd.bar(
            np.arange(len(CHR_LEVELS)),
            sd_vals,
            color=OKABE_ITO["black"],
            width=0.8,
        )
        ax_sd.set_xlim(-0.5, len(CHR_LEVELS) - 0.5)
        ax_sd.set_xticks(np.arange(len(CHR_LEVELS)))
        ax_sd.set_xticklabels(chr_labels)
        ax_sd.set_ylabel("SD of donor\nmeans")
        ax_sd.set_title(
            "Between-donor variation among MBL+CLL+RT (SD of chromosome-mean log2 deviation)",
            loc="left",
        )
        ax_sd.set_ylim(0, max(sd_vals) * 1.15)

        offsets = {"MBL": -0.22, "CLL": 0.0, "RT": 0.22}
        for i, chr_name in enumerate(CHR_LEVELS):
            for stage in OBS_STAGES:
                rec = next(
                    r
                    for r in stage_rows
                    if r["stage"] == stage and r["chr"] == chr_name
                )
                x = i + offsets[stage]
                y = float(rec["mean_minus_normal_b"])
                sd = float(rec["sd_of_donor_means"])
                ax_diff.errorbar(
                    x,
                    y,
                    yerr=sd,
                    fmt=STAGE_MARKERS[stage],
                    color=STAGE_COLORS[stage],
                    ecolor=STAGE_COLORS[stage],
                    elinewidth=0.7,
                    capsize=1.2,
                    markersize=4.5,
                    markeredgecolor="white",
                    markeredgewidth=0.3,
                    zorder=3,
                )
        ax_diff.axhline(0.0, color="0.3", linewidth=0.6, linestyle="-")
        ax_diff.set_xlim(-0.5, len(CHR_LEVELS) - 0.5)
        ax_diff.set_xticks(np.arange(len(CHR_LEVELS)))
        ax_diff.set_xticklabels(chr_labels)
        ax_diff.set_ylabel("Mean minus\nNormal_B")
        ax_diff.set_xlabel("Chromosome (autosomes; gene counts in source table)")
        ax_diff.set_title(
            "Stage mean of donor chromosome means versus Normal_B\n"
            f"(n={n_by_stage['Normal_B']} Normal_B donors). Error bars: SD of donor means.",
            loc="left",
        )
        ax_diff.legend(
            handles=[
                Line2D(
                    [0],
                    [0],
                    marker=STAGE_MARKERS[st],
                    color=STAGE_COLORS[st],
                    markeredgecolor="white",
                    linestyle="None",
                    markersize=5,
                    label=f"{st} (n={n_by_stage[st]})",
                )
                for st in OBS_STAGES
            ],
            loc="lower left",
            ncol=3,
            frameon=False,
        )

        ax_cbar.set_box_aspect(15)
        cbar = fig.colorbar(im, cax=ax_cbar)
        cbar.set_label("inferCNV log2 deviation\n(chromosome mean; 0 = diploid-centered)")

        caption = DATA_DIR / "infercnv_chromosome_donor_summary_caption.txt"
        caption.write_text(
            "\n".join(
                [
                    "Figure X. Chromosome-level summary of inferCNV signal across donors.",
                    "The official inferCNV heatmap is cell-by-gene. Here each value is the chromosome-wide mean of log2(inferCNV expr.data) for one donor (cells averaged within donor; donors with fewer than 5 cells excluded).",
                    f"Panel (top): one row per MBL, CLL, or RT donor, ordered within stage by chromosome-12 mean (descending). Color is a diverging map centered at 0 with limits ±{abs_limit:.2f}. Missing values would be gray; none are present.",
                    f"Panel (middle): between-donor SD among the {n_by_stage['MBL']+n_by_stage['CLL']+n_by_stage['RT']} observation donors; this is sample-to-sample spread, not a standard error.",
                    f"Panel (bottom): stage mean of those donor means minus the Normal_B donor mean (n={n_by_stage['Normal_B']}). Error bars are SD of donor means. Markers encode stage in addition to color.",
                    "HMM copy-number calling was not run. Whole-chromosome means under-represent focal events; companion table donor_chromosome_focal_bin.tsv stores the most extreme 80-gene bin per donor and chromosome.",
                    "These are RNA-expression-based inferCNV estimates relative to Normal_B, not genomic CNA calls. Publisher formatting was not applied.",
                    "Gene counts per chromosome: "
                    + ", ".join(f"{c.replace('chr','')}={gene_n[c]}" for c in CHR_LEVELS)
                    + ".",
                ]
            )
            + "\n"
        )

        export_figure(
            fig,
            OUT_BASE,
            formats=["pdf", "png"],
            dpi=600,
            bbox_inches=None,
            overwrite=True,
            mkdir=True,
            write_manifest=True,
            provenance={
                "raw_data": str(
                    PROJECT
                    / "CLLCMA_alpha/results/bcell_stage_deseq2/publication_figures/"
                    / "infercnv_bcell_stage/infercnv_obj_final.rds"
                ),
                "transformations": [
                    "log2(max(inferCNV expr.data, 1e-6))",
                    "chromosome-wide mean across genes",
                    "donor mean of cells; exclude donors with <5 cells",
                    "no HMM CNV calling",
                    "heatmap rows: observation donors ordered by stage then chr12 mean",
                ],
                "uncertainty": "SD of donor-level chromosome means; unit of replication = donor",
                "missing_data": "donors with <5 cells excluded; no missing chromosomes in retained donors",
                "color_scale": f"RdBu_r TwoSlopeNorm vcenter=0, limits ±{abs_limit:.2f}",
                "publisher_status": "pending verification; provisional general figure",
            },
        )
        plt.close(fig)


if __name__ == "__main__":
    main()
