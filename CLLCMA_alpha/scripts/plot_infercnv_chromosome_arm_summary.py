#!/usr/bin/env python3
"""Polished donor-level inferCNV summary at chromosome-arm resolution."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure_io import export_figure, style_context  # noqa: E402

PROJECT = Path(os.environ.get("CLL_TME_DIR", Path.home() / "research" / "CLL_TME")).resolve()
DATA_DIR = (
    PROJECT
    / "CLLCMA_alpha/results/bcell_stage_deseq2/publication_figures/"
    / "infercnv_bcell_stage/chromosome_summary"
)
OUT_BASE = DATA_DIR / "infercnv_chromosome_arm_donor_summary"

STAGE_ORDER = ("Normal_B", "MBL", "CLL", "RT")
OBS_STAGES = ("MBL", "CLL", "RT")
ARM_LEVELS = tuple(
    arm for chromosome in range(1, 23) for arm in (f"{chromosome}p", f"{chromosome}q")
)
STAGE_COLORS = {
    "Normal_B": "#0072B2",
    "MBL": "#009E73",
    "CLL": "#E69F00",
    "RT": "#D55E00",
}
STAGE_MARKERS = {"MBL": "s", "CLL": "^", "RT": "D"}
LOSS_COLOR = "#0072B2"
NEUTRAL_COLOR = "#F7F7F7"
GAIN_COLOR = "#D55E00"
MISSING_COLOR = "#9E9E9E"
ABS_LIMIT = 0.10


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def donor_arm_matrix(
    rows: list[dict[str, str]],
) -> tuple[np.ndarray, list[str]]:
    values: dict[tuple[str, str, str], float] = {}
    donors_by_stage: dict[str, set[str]] = {stage: set() for stage in OBS_STAGES}
    for row in rows:
        stage, donor, arm = row["stage"], row["donor"], row["arm"]
        values[(stage, donor, arm)] = float(row["mean_log2_dev"])
        if stage in donors_by_stage:
            donors_by_stage[stage].add(donor)

    ordered: list[tuple[str, str]] = []
    for stage in OBS_STAGES:
        donors = sorted(donors_by_stage[stage])
        donors.sort(
            key=lambda donor: np.nanmean(
                [
                    values.get((stage, donor, "12p"), np.nan),
                    values.get((stage, donor, "12q"), np.nan),
                ]
            ),
            reverse=True,
        )
        ordered.extend((stage, donor) for donor in donors)

    matrix = np.full((len(ordered), len(ARM_LEVELS)), np.nan)
    stages: list[str] = []
    for row_index, (stage, donor) in enumerate(ordered):
        stages.append(stage)
        for column_index, arm in enumerate(ARM_LEVELS):
            matrix[row_index, column_index] = values.get(
                (stage, donor, arm),
                np.nan,
            )
    return matrix, stages


def add_chromosome_boundaries(ax: mpl.axes.Axes) -> None:
    for boundary in np.arange(1.5, len(ARM_LEVELS) - 0.5, 2):
        ax.axvline(boundary, color="black", linewidth=0.35, alpha=0.55)


def bold_axis_text(ax: mpl.axes.Axes) -> None:
    for label in [*ax.get_xticklabels(), *ax.get_yticklabels()]:
        label.set_fontweight("bold")
    ax.xaxis.label.set_fontweight("bold")
    ax.yaxis.label.set_fontweight("bold")
    ax.title.set_fontweight("bold")


def write_findings(
    stage_rows: list[dict[str, str]],
    rank_rows: list[dict[str, str]],
    gene_counts: dict[str, int],
) -> None:
    by_shift = sorted(
        rank_rows,
        key=lambda row: abs(float(row["mean_minus_normal_b"])),
        reverse=True,
    )
    by_sd = sorted(
        rank_rows,
        key=lambda row: float(row["sd_observation_donors"]),
        reverse=True,
    )
    unavailable = [arm for arm in ARM_LEVELS if gene_counts.get(arm, 0) < 10]
    lines = [
        "inferCNV chromosome-arm summary (RNA-based; HMM was not run).",
        "Unit: donor. Values are arm-wide means of log2(inferCNV expr.data).",
        "Arms require at least 10 mapped genes; unavailable arms are gray.",
        "",
        "Largest systematic arm shifts versus Normal_B:",
    ]
    for row in by_shift[:10]:
        lines.append(
            f"  {row['arm']}: mean-Normal_B={float(row['mean_minus_normal_b']):+.4f}; "
            f"SD={float(row['sd_observation_donors']):.4f}"
        )
    lines.append("")
    lines.append("Largest between-donor arm heterogeneity:")
    for row in by_sd[:10]:
        lines.append(
            f"  {row['arm']}: SD={float(row['sd_observation_donors']):.4f}; "
            f"mean-Normal_B={float(row['mean_minus_normal_b']):+.4f}"
        )
    lines += [
        "",
        "Unavailable (<10 mapped genes): " + ", ".join(unavailable),
        "",
        "These are expression-derived inferCNV estimates, not genomic CNA calls.",
        "Arm averaging can still dilute focal lesions; validate against matched genomic CNA data.",
    ]
    (DATA_DIR / "CHROMOSOME_ARM_KEY_FINDINGS.txt").write_text("\n".join(lines) + "\n")


def main() -> None:
    donor_rows = read_tsv(DATA_DIR / "donor_chromosome_arm_mean_log2_deviation.tsv")
    stage_rows = read_tsv(DATA_DIR / "chromosome_arm_stage_summary.tsv")
    rank_rows = read_tsv(DATA_DIR / "chromosome_arm_ranked_by_donor_mad.tsv")
    donor_counts = read_tsv(DATA_DIR / "donor_counts_by_stage.tsv")
    gene_count_rows = read_tsv(DATA_DIR / "chromosome_arm_gene_counts.tsv")
    gene_counts = {row["arm"]: int(row["N"]) for row in gene_count_rows}
    usable = {row["arm"]: row["usable"].upper() == "TRUE" for row in gene_count_rows}
    n_by_stage = {row["stage"]: int(row["n_donors"]) for row in donor_counts}

    obs_matrix, donor_stages = donor_arm_matrix(donor_rows)
    stage_lookup = {
        (row["stage"], row["arm"]): row
        for row in stage_rows
    }
    rank_lookup = {row["arm"]: row for row in rank_rows}
    write_findings(stage_rows, rank_rows, gene_counts)

    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        "infercnv_blue_white_orange",
        [LOSS_COLOR, NEUTRAL_COLOR, GAIN_COLOR],
    ).with_extremes(bad=MISSING_COLOR)
    norm = mpl.colors.TwoSlopeNorm(vmin=-ABS_LIMIT, vcenter=0, vmax=ABS_LIMIT)
    arm_positions = np.arange(len(ARM_LEVELS))

    with style_context("default", palette_name="okabe_ito_on_white"):
        mpl.rcParams.update(
            {
                "font.weight": "bold",
                "axes.labelweight": "bold",
                "axes.titleweight": "bold",
            }
        )
        fig = plt.figure(
            figsize=(220 / 25.4, 210 / 25.4),
            layout="constrained",
        )
        fig.set_facecolor("white")
        axes = fig.subplot_mosaic(
            [
                ["stage", "heat", "cbar"],
                ["blank_sd", "sd", "cbar_blank_sd"],
                ["blank_effect", "effect", "cbar_blank_effect"],
            ],
            width_ratios=[0.055, 1.0, 0.022],
            height_ratios=[2.55, 0.82, 0.82],
        )
        for key in ("blank_sd", "blank_effect", "cbar_blank_sd", "cbar_blank_effect"):
            axes[key].axis("off")

        stage_index = np.array(
            [{"MBL": 0, "CLL": 1, "RT": 2}[stage] for stage in donor_stages]
        )[:, None]
        stage_cmap = mpl.colors.ListedColormap(
            [STAGE_COLORS["MBL"], STAGE_COLORS["CLL"], STAGE_COLORS["RT"]]
        )
        ax_stage = axes["stage"]
        ax_stage.imshow(
            stage_index,
            aspect="auto",
            cmap=stage_cmap,
            vmin=0,
            vmax=2,
            interpolation="nearest",
        )
        ax_stage.set_xticks([])
        ax_stage.set_yticks([])
        ax_stage.set_ylabel("Donors")

        ax_heat = axes["heat"]
        image = ax_heat.imshow(
            obs_matrix,
            aspect="auto",
            cmap=cmap,
            norm=norm,
            interpolation="nearest",
        )
        ax_heat.set_xticks(arm_positions)
        ax_heat.set_xticklabels(ARM_LEVELS, rotation=90, fontsize=6)
        ax_heat.set_yticks([])
        ax_heat.set_xlabel("Chromosome arm (GRCh38)")
        ax_heat.set_title(
            "Donor-mean inferCNV signal by chromosome arm "
            f"(MBL n={n_by_stage['MBL']}, CLL n={n_by_stage['CLL']}, RT n={n_by_stage['RT']})",
            loc="left",
            fontsize=10,
        )
        add_chromosome_boundaries(ax_heat)

        stage_boundaries = np.cumsum(
            [n_by_stage["MBL"], n_by_stage["CLL"], n_by_stage["RT"]]
        )[:-1]
        for boundary in stage_boundaries:
            for ax in (ax_stage, ax_heat):
                ax.axhline(boundary - 0.5, color="black", linewidth=0.7)
        stage_centers = {
            "MBL": n_by_stage["MBL"] / 2 - 0.5,
            "CLL": n_by_stage["MBL"] + n_by_stage["CLL"] / 2 - 0.5,
            "RT": (
                n_by_stage["MBL"]
                + n_by_stage["CLL"]
                + n_by_stage["RT"] / 2
                - 0.5
            ),
        }
        for stage, center in stage_centers.items():
            ax_stage.text(
                -1.7,
                center,
                f"{stage}\nn={n_by_stage[stage]}",
                ha="right",
                va="center",
                fontsize=7,
                fontweight="bold",
                color=STAGE_COLORS[stage],
                clip_on=False,
            )

        ax_sd = axes["sd"]
        sd_values = np.array(
            [
                float(rank_lookup[arm]["sd_observation_donors"])
                if arm in rank_lookup and usable.get(arm, False)
                else np.nan
                for arm in ARM_LEVELS
            ]
        )
        bar_colors = [
            "#4C566A" if np.isfinite(value) else MISSING_COLOR
            for value in sd_values
        ]
        ax_sd.bar(
            arm_positions,
            np.nan_to_num(sd_values, nan=0),
            color=bar_colors,
            width=0.78,
        )
        ax_sd.set_xlim(-0.5, len(ARM_LEVELS) - 0.5)
        ax_sd.set_xticks(arm_positions)
        ax_sd.set_xticklabels(ARM_LEVELS, rotation=90, fontsize=6)
        ax_sd.set_ylabel("SD of donor\narm means")
        ax_sd.set_title(
            "Between-donor variation among MBL, CLL, and RT",
            loc="left",
            fontsize=9,
        )
        ax_sd.set_ylim(0, np.nanmax(sd_values) * 1.14)
        add_chromosome_boundaries(ax_sd)

        ax_effect = axes["effect"]
        offsets = {"MBL": -0.22, "CLL": 0, "RT": 0.22}
        for arm_index, arm in enumerate(ARM_LEVELS):
            if not usable.get(arm, False):
                continue
            for stage in OBS_STAGES:
                row = stage_lookup.get((stage, arm))
                if row is None:
                    continue
                ax_effect.errorbar(
                    arm_index + offsets[stage],
                    float(row["mean_minus_normal_b"]),
                    yerr=float(row["sd_of_donor_means"]),
                    fmt=STAGE_MARKERS[stage],
                    color=STAGE_COLORS[stage],
                    ecolor=STAGE_COLORS[stage],
                    elinewidth=0.7,
                    capsize=1,
                    markersize=4,
                    markeredgecolor="white",
                    markeredgewidth=0.35,
                    zorder=3,
                )
        ax_effect.axhline(0, color="black", linewidth=0.6)
        ax_effect.set_xlim(-0.5, len(ARM_LEVELS) - 0.5)
        ax_effect.set_xticks(arm_positions)
        ax_effect.set_xticklabels(ARM_LEVELS, rotation=90, fontsize=6)
        ax_effect.set_ylabel("Mean minus\nNormal_B")
        ax_effect.set_xlabel("Chromosome arm; gray denotes <10 mapped genes")
        ax_effect.set_title(
            f"Stage mean versus Normal_B (n={n_by_stage['Normal_B']}); error bars: donor SD",
            loc="left",
            fontsize=9,
        )
        add_chromosome_boundaries(ax_effect)
        ax_effect.legend(
            handles=[
                Line2D(
                    [0],
                    [0],
                    marker=STAGE_MARKERS[stage],
                    color=STAGE_COLORS[stage],
                    linestyle="None",
                    markeredgecolor="white",
                    markersize=5,
                    label=f"{stage} (n={n_by_stage[stage]})",
                )
                for stage in OBS_STAGES
            ],
            loc="lower left",
            ncol=3,
            frameon=False,
            prop={"weight": "bold", "size": 7},
        )

        colorbar = fig.colorbar(image, cax=axes["cbar"])
        colorbar.set_label(
            "Arm-mean inferCNV log2 deviation\n"
            "blue: relative loss; orange: relative gain",
            fontweight="bold",
        )
        for label in colorbar.ax.get_yticklabels():
            label.set_fontweight("bold")

        for ax in (ax_stage, ax_heat, ax_sd, ax_effect):
            bold_axis_text(ax)

        unavailable = [arm for arm in ARM_LEVELS if not usable.get(arm, False)]
        caption_lines = [
            "Figure X. Donor-level inferCNV signal summarized by chromosome arm.",
            "Genes were assigned to GRCh38 p or q arms by midpoint using UCSC hg38 cytoband boundaries.",
            "Each heatmap value is the arm-wide mean log2 deviation of final inferCNV expression, averaged across cells within a donor.",
            "Rows are grouped by stage and ordered within stage by mean chromosome-12 arm signal.",
            f"The diverging color scale is centered at zero and fixed to \u00b1{ABS_LIMIT:.2f}; values outside this range are color-saturated.",
            "The middle panel shows between-donor SD. The lower panel shows stage means minus Normal_B with donor SD.",
            "Arms with fewer than 10 mapped genes are gray and excluded from arm-level estimates: "
            + ", ".join(unavailable)
            + ".",
            "HMM calling was not run; these are RNA-expression-based inferCNV estimates, not genomic CNA calls.",
            "Arm averaging can dilute focal lesions and should be validated against matched genomic CNA data.",
        ]
        (DATA_DIR / "infercnv_chromosome_arm_donor_summary_caption.txt").write_text(
            "\n".join(caption_lines) + "\n"
        )

        export_figure(
            fig,
            OUT_BASE,
            formats=["pdf", "png"],
            dpi=600,
            bbox_inches=None,
            overwrite=True,
            write_manifest=True,
            provenance={
                "raw_data": str(
                    DATA_DIR.parent / "infercnv_obj_final.rds"
                ),
                "coordinate_system": "GRCh38/hg38",
                "arm_boundaries": str(
                    PROJECT / "CLLCMA_alpha/references/hg38_chromosome_arm_boundaries.tsv"
                ),
                "arm_boundary_source": "UCSC hg38 cytoBand; retrieved 2026-09-08",
                "transformations": [
                    "log2(max(inferCNV expr.data, 1e-6))",
                    "gene assigned by midpoint to p or q arm",
                    "arm-wide gene mean per cell",
                    "cell mean within donor; donors with <5 cells excluded",
                    "arms with <10 mapped genes shown as unavailable",
                    "no HMM CNV calling",
                ],
                "uncertainty": "SD of donor-level arm means",
                "missing_data": "unavailable arms colored gray",
                "color_scale": "Okabe-Ito blue-white-vermillion, centered 0, limits \u00b10.10",
                "publisher_status": "provisional general figure; journal not specified",
            },
        )
        plt.close(fig)


if __name__ == "__main__":
    main()
