#!/usr/bin/env python3
"""Network figures for splicing-outlier × LR prioritization (Idea A).

Reads existing tables under rank_enrichment/ligand_receptor_prioritization/
and writes:
  - lr_network_bipartite.{pdf,png}     ligand | receptor bipartite graph
  - lr_network_celltype.{pdf,png}      sender→receiver cell-type graph
  - lr_network_multipartite.{pdf,png}  sender → ligand → receptor → receiver

Env:
  CLL_TME_DIR, DECONV_METHOD
  LR_NET_TOP (default 35) — edges kept per figure
"""

from __future__ import annotations

import os
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

COL_LIGAND = "#E69F00"
COL_RECEPTOR = "#56B4E9"
COL_CELL = "#009E73"
COL_POS = "#B2182B"
COL_NEG = "#2166AC"
COL_EDGE = "#666666"


def project_root() -> Path:
    env = os.environ.get("CLL_TME_DIR")
    if env:
        return Path(env).expanduser().resolve()
    return Path(__file__).resolve().parents[2]


def _short(s: str, n: int = 14) -> str:
    s = str(s)
    return s if len(s) <= n else s[: n - 1] + "…"


def _cell_colors(names: list[str]) -> dict[str, str]:
    # Okabe-Ito-ish cycle
    palette = [
        "#009E73",
        "#E69F00",
        "#56B4E9",
        "#CC79A7",
        "#0072B2",
        "#D55E00",
        "#F0E442",
        "#000000",
        "#999999",
        "#882255",
        "#44AA99",
        "#117733",
    ]
    return {n: palette[i % len(palette)] for i, n in enumerate(sorted(set(names)))}


def plot_bipartite_lr(within: pd.DataFrame, out_stem: Path, top_n: int = 35) -> None:
    """Ligands (left) — receptors (right); edge width = pair strength."""
    d = within.nlargest(top_n, "PAIR_STRENGTH").copy()
    if d.empty:
        return

    ligands = list(dict.fromkeys(d["LIGAND"].tolist()))
    receptors = list(dict.fromkeys(d["RECEPTOR"].tolist()))
    # y positions
    def _ys(nodes: list[str]) -> dict[str, float]:
        if len(nodes) == 1:
            return {nodes[0]: 0.5}
        return {n: i / (len(nodes) - 1) for i, n in enumerate(nodes)}

    yL, yR = _ys(ligands), _ys(receptors)
    ct_cols = _cell_colors(d["CELL_TYPE"].tolist())

    fig_h = max(6.0, 0.28 * max(len(ligands), len(receptors)) + 1.5)
    fig, ax = plt.subplots(figsize=(10, fig_h), layout="constrained")

    smax = d["PAIR_STRENGTH"].max()
    for _, r in d.iterrows():
        x0, y0 = 0.0, yL[r["LIGAND"]]
        x1, y1 = 1.0, yR[r["RECEPTOR"]]
        lw = 0.6 + 3.5 * (r["PAIR_STRENGTH"] / smax)
        color = ct_cols.get(r["CELL_TYPE"], COL_EDGE)
        # sign of product for linestyle
        ls = "-" if r["PAIR_PRODUCT"] >= 0 else "--"
        ax.plot([x0, x1], [y0, y1], color=color, lw=lw, alpha=0.75, solid_capstyle="round", linestyle=ls, zorder=1)

    ax.scatter(
        [0] * len(ligands),
        [yL[n] for n in ligands],
        s=120,
        c=COL_LIGAND,
        edgecolors="white",
        linewidths=0.6,
        zorder=3,
    )
    ax.scatter(
        [1] * len(receptors),
        [yR[n] for n in receptors],
        s=120,
        c=COL_RECEPTOR,
        edgecolors="white",
        linewidths=0.6,
        zorder=3,
    )
    for n in ligands:
        ax.text(-0.03, yL[n], n, ha="right", va="center", fontsize=7.5)
    for n in receptors:
        ax.text(1.03, yR[n], n, ha="left", va="center", fontsize=7.5)

    ax.set_xlim(-0.45, 1.45)
    ax.set_ylim(-0.05, 1.05)
    ax.axis("off")
    ax.set_title(
        f"Bipartite LR network (top {len(d)} within–cell-type pairs)\n"
        "Edge color = TME cell type; width ∝ |score_L|+|score_R|; dashed = opposite-sign scores",
        fontsize=10,
        loc="left",
    )
    # legends
    ct_handles = [
        Line2D([0], [0], color=ct_cols[c], lw=2.5, label=c)
        for c in sorted(d["CELL_TYPE"].unique())
    ]
    type_handles = [
        mpatches.Patch(color=COL_LIGAND, label="Ligand"),
        mpatches.Patch(color=COL_RECEPTOR, label="Receptor"),
        Line2D([0], [0], color="#333", lw=1.5, linestyle="-", label="Same-sign scores"),
        Line2D([0], [0], color="#333", lw=1.5, linestyle="--", label="Opposite-sign scores"),
    ]
    leg1 = ax.legend(handles=type_handles, loc="upper left", fontsize=7, frameon=False, bbox_to_anchor=(-0.02, 1.02))
    ax.add_artist(leg1)
    ax.legend(handles=ct_handles, loc="lower left", fontsize=6.5, frameon=False, ncol=2, bbox_to_anchor=(-0.02, -0.02))

    fig.savefig(f"{out_stem}.pdf")
    fig.savefig(f"{out_stem}.png", dpi=350, facecolor="white")
    plt.close(fig)
    print(f"Wrote {out_stem}.pdf/.png", flush=True)


def plot_celltype_network(cross: pd.DataFrame, out_stem: Path, top_n: int = 40) -> None:
    """Cell-type nodes; edge = aggregated cross LR strength."""
    if cross.empty:
        return
    d = cross.nlargest(max(top_n * 3, 50), "PAIR_STRENGTH").copy()
    edges = (
        d.groupby(["SENDER_CELL", "RECEIVER_CELL"], as_index=False)
        .agg(
            N_PAIRS=("PAIR", "nunique"),
            MAX_STRENGTH=("PAIR_STRENGTH", "max"),
            MEAN_STRENGTH=("PAIR_STRENGTH", "mean"),
            TOP_PAIR=("PAIR", "first"),
        )
        .sort_values("MAX_STRENGTH", ascending=False)
        .head(top_n)
    )
    nodes = sorted(set(edges["SENDER_CELL"]) | set(edges["RECEIVER_CELL"]))
    # circular layout
    n = len(nodes)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    pos = {name: (np.cos(a), np.sin(a)) for name, a in zip(nodes, angles)}
    cols = _cell_colors(nodes)

    fig, ax = plt.subplots(figsize=(8.5, 8.5), layout="constrained")
    smax = edges["MAX_STRENGTH"].max()
    for _, r in edges.iterrows():
        x0, y0 = pos[r["SENDER_CELL"]]
        x1, y1 = pos[r["RECEIVER_CELL"]]
        lw = 0.5 + 4.0 * (r["MAX_STRENGTH"] / smax)
        # slight curve via midpoint offset
        xm, ym = (x0 + x1) / 2, (y0 + y1) / 2
        # perpendicular nudge
        dx, dy = x1 - x0, y1 - y0
        nx, ny = -dy, dx
        norm = np.hypot(nx, ny) + 1e-9
        xm += 0.08 * nx / norm
        ym += 0.08 * ny / norm
        ax.plot([x0, xm, x1], [y0, ym, y1], color="#555555", lw=lw, alpha=0.55, zorder=1)
        # arrow head near receiver
        ax.annotate(
            "",
            xy=(x1, y1),
            xytext=(xm, ym),
            arrowprops=dict(arrowstyle="->", color="#555555", lw=lw * 0.7, mutation_scale=10),
            zorder=2,
        )

    for name in nodes:
        x, y = pos[name]
        ax.scatter([x], [y], s=600, c=cols[name], edgecolors="white", linewidths=1.2, zorder=3)
        ax.text(x * 1.18, y * 1.18, name, ha="center", va="center", fontsize=8)

    ax.set_xlim(-1.55, 1.55)
    ax.set_ylim(-1.55, 1.55)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(
        "Cell–cell LR network (sender → receiver)\n"
        "Edge width ∝ max pair strength; aggregated from cross–cell-type LR table",
        fontsize=10,
    )
    # annotate a few top edges with pair name
    for _, r in edges.head(8).iterrows():
        x0, y0 = pos[r["SENDER_CELL"]]
        x1, y1 = pos[r["RECEIVER_CELL"]]
        ax.text(
            (x0 + x1) / 2,
            (y0 + y1) / 2,
            _short(r["TOP_PAIR"], 16),
            fontsize=5.5,
            color="#222",
            ha="center",
            va="center",
            bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.7),
            zorder=4,
        )

    fig.savefig(f"{out_stem}.pdf")
    fig.savefig(f"{out_stem}.png", dpi=350, facecolor="white")
    plt.close(fig)
    print(f"Wrote {out_stem}.pdf/.png", flush=True)


def plot_multipartite(cross: pd.DataFrame, out_stem: Path, top_n: int = 25) -> None:
    """Four columns: sender cell | ligand | receptor | receiver cell."""
    d = cross.nlargest(top_n, "PAIR_STRENGTH").copy()
    if d.empty:
        return

    senders = list(dict.fromkeys(d["SENDER_CELL"]))
    ligands = list(dict.fromkeys(d["LIGAND"]))
    receptors = list(dict.fromkeys(d["RECEPTOR"]))
    receivers = list(dict.fromkeys(d["RECEIVER_CELL"]))

    def _ys(nodes: list[str]) -> dict[str, float]:
        if len(nodes) == 1:
            return {nodes[0]: 0.5}
        return {n: i / (len(nodes) - 1) for i, n in enumerate(nodes)}

    yS, yL, yR, yT = _ys(senders), _ys(ligands), _ys(receptors), _ys(receivers)
    xs = [0.0, 0.33, 0.66, 1.0]
    smax = d["PAIR_STRENGTH"].max()

    fig_h = max(7.0, 0.25 * max(len(senders), len(ligands), len(receptors), len(receivers)) + 1.8)
    fig, ax = plt.subplots(figsize=(12, fig_h), layout="constrained")

    for _, r in d.iterrows():
        path = [
            (xs[0], yS[r["SENDER_CELL"]]),
            (xs[1], yL[r["LIGAND"]]),
            (xs[2], yR[r["RECEPTOR"]]),
            (xs[3], yT[r["RECEIVER_CELL"]]),
        ]
        lw = 0.5 + 2.8 * (r["PAIR_STRENGTH"] / smax)
        col = COL_POS if r["PAIR_PRODUCT"] >= 0 else COL_NEG
        xs_p, ys_p = zip(*path)
        ax.plot(xs_p, ys_p, color=col, lw=lw, alpha=0.65, zorder=1)

    layers = [
        (xs[0], senders, yS, COL_CELL, "Sender"),
        (xs[1], ligands, yL, COL_LIGAND, "Ligand"),
        (xs[2], receptors, yR, COL_RECEPTOR, "Receptor"),
        (xs[3], receivers, yT, COL_CELL, "Receiver"),
    ]
    for x, nodes, ymap, color, title in layers:
        ax.scatter(
            [x] * len(nodes),
            [ymap[n] for n in nodes],
            s=90,
            c=color,
            edgecolors="white",
            linewidths=0.5,
            zorder=3,
        )
        for n in nodes:
            ha = "right" if x == 0 else ("left" if x == 1.0 else "center")
            dx = -0.02 if x == 0 else (0.02 if x == 1.0 else 0.0)
            va = "bottom" if x in (0.33, 0.66) else "center"
            ax.text(
                x + dx,
                ymap[n] + (0.015 if x in (0.33, 0.66) else 0),
                _short(n, 16),
                ha=ha,
                va=va,
                fontsize=6.5,
                zorder=4,
            )
        ax.text(x, 1.06, title, ha="center", va="bottom", fontsize=9, fontweight="bold")

    ax.set_xlim(-0.2, 1.2)
    ax.set_ylim(-0.05, 1.12)
    ax.axis("off")
    ax.set_title(
        f"Multipartite LR network (top {len(d)} cross pairs)\n"
        "Red = same-sign score product; blue = opposite-sign; width ∝ strength",
        fontsize=10,
        loc="left",
    )
    fig.savefig(f"{out_stem}.pdf")
    fig.savefig(f"{out_stem}.png", dpi=350, facecolor="white")
    plt.close(fig)
    print(f"Wrote {out_stem}.pdf/.png", flush=True)


def main() -> int:
    root = project_root()
    method = os.environ.get("DECONV_METHOD", "bayesprism").lower()
    top_n = int(os.environ.get("LR_NET_TOP", "35"))
    out_dir = (
        root
        / "CLLCMA_alpha/results/deconvolution/cll_broad_2022"
        / "splicing_outlier_tme_association"
        / method
        / "rank_enrichment"
        / "ligand_receptor_prioritization"
    )
    within_path = out_dir / "lr_pairs_within_celltype.tsv.gz"
    cross_path = out_dir / "lr_pairs_cross_celltype.tsv.gz"
    if not within_path.exists() or not cross_path.exists():
        raise SystemExit(f"Missing pair tables under {out_dir}")

    within = pd.read_csv(within_path, sep="\t")
    cross = pd.read_csv(cross_path, sep="\t")

    plot_bipartite_lr(within, out_dir / "lr_network_bipartite", top_n=top_n)
    plot_celltype_network(cross, out_dir / "lr_network_celltype", top_n=min(top_n, 30))
    plot_multipartite(cross, out_dir / "lr_network_multipartite", top_n=min(top_n, 25))

    # also a focused CLL-centric bipartite from cross table
    cll = cross[cross["SENDER_CELL"].eq("CLL") | cross["RECEIVER_CELL"].eq("CLL")].copy()
    if len(cll):
        # reshape to within-like for bipartite using pair strength
        fake = cll.rename(columns={"PAIR_PRODUCT": "PAIR_PRODUCT"})
        fake["CELL_TYPE"] = fake["SENDER_CELL"] + "→" + fake["RECEIVER_CELL"]
        plot_bipartite_lr(
            fake,
            out_dir / "lr_network_bipartite_CLL_focus",
            top_n=min(30, len(fake)),
        )

    print(f"Done -> {out_dir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
