#!/usr/bin/env python3
"""Idea A: splicing-outlier rank × ligand–receptor prioritization.

Uses existing gene ranks (Cliff's δ × −log10 P_adj) from
splicing_outlier_tme_association and OmniPath ligrecextra (+ optional
CellPhoneDB simple pairs) to:

  1) Test enrichment of ligands / receptors in positive and negative score tails
  2) Prioritize LR pairs whose ligand and receptor are highly ranked
  3) Build cross–cell-type sender→receiver pair tables for key TME compartments
  4) Write heatmaps / network-style summary figures

No network download in this script — expects local files under
CLLCMA_alpha/references/ligand_receptor/ (approved separately).

Env:
  CLL_TME_DIR, DECONV_METHOD (default bayesprism)
  LR_TOP_N (default 500) — tail size for ORA
  LR_PAIR_TOP (default 50) — pairs kept per direction / sender-receiver
  LR_MIN_ABS_SCORE (default 0.5) — min |score| to keep a gene in pair tables
"""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

DEFAULT_TOP_N = 500
DEFAULT_PAIR_TOP = 50
DEFAULT_MIN_ABS = 0.5

KEY_CELL_TYPES = [
    "CLL",
    "CD8aa",
    "T_proliferating",
    "Monocyte_CD16",
    "Monocyte_CD14",
    "Monocyte_intermediate",
    "CD8_T_naive",
    "CD8_T_memory",
    "CD4_T_naive",
    "NK_proliferating",
    "Treg",
    "ASDC",
]


def project_root() -> Path:
    env = os.environ.get("CLL_TME_DIR")
    if env:
        return Path(env).expanduser().resolve()
    return Path(__file__).resolve().parents[2]


def bh_fdr(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    n = len(p)
    q = np.full(n, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return q
    idx = np.where(ok)[0]
    order = idx[np.argsort(p[idx])]
    ranked = p[order]
    m = len(ranked)
    raw = ranked * m / np.arange(1, m + 1)
    raw = np.minimum.accumulate(raw[::-1])[::-1]
    q[order] = np.clip(raw, 0, 1)
    return q


def load_omnipath_lr(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df = df.rename(
        columns={
            "source_genesymbol": "LIGAND",
            "target_genesymbol": "RECEPTOR",
        }
    )
    df["LIGAND"] = df["LIGAND"].astype(str).str.upper().str.strip()
    df["RECEPTOR"] = df["RECEPTOR"].astype(str).str.upper().str.strip()
    df = df[
        (df["LIGAND"] != "")
        & (df["RECEPTOR"] != "")
        & (df["LIGAND"] != "NAN")
        & (df["RECEPTOR"] != "NAN")
    ]
    df["PAIR"] = df["LIGAND"] + "→" + df["RECEPTOR"]
    # collapse duplicate directed edges
    keep_cols = ["LIGAND", "RECEPTOR", "PAIR"]
    agg = df[keep_cols].drop_duplicates()
    agg["DATABASE"] = "OmniPath_ligrecextra"
    if "is_stimulation" in df.columns:
        stim = df.groupby(["LIGAND", "RECEPTOR"])["is_stimulation"].any()
        inhib = df.groupby(["LIGAND", "RECEPTOR"])["is_inhibition"].any()
        agg = agg.merge(
            stim.rename("IS_STIM").reset_index(), on=["LIGAND", "RECEPTOR"], how="left"
        )
        agg = agg.merge(
            inhib.rename("IS_INHIB").reset_index(), on=["LIGAND", "RECEPTOR"], how="left"
        )
    return agg


def load_cellphone_simple_pairs(gene_csv: Path, int_csv: Path) -> pd.DataFrame:
    """Keep CellPhoneDB interactions where both partners are single UniProt proteins."""
    genes = pd.read_csv(gene_csv)
    genes["uniprot"] = genes["uniprot"].astype(str)
    genes["hgnc_symbol"] = genes["hgnc_symbol"].astype(str).str.upper().str.strip()
    u2sym = (
        genes.dropna(subset=["hgnc_symbol"])
        .drop_duplicates("uniprot")
        .set_index("uniprot")["hgnc_symbol"]
        .to_dict()
    )
    inter = pd.read_csv(int_csv)
    rows = []
    for _, r in inter.iterrows():
        a, b = str(r["partner_a"]), str(r["partner_b"])
        # skip named complexes
        if "complex" in a.lower() or "complex" in b.lower():
            continue
        if a not in u2sym or b not in u2sym:
            continue
        la, rb = u2sym[a], u2sym[b]
        if not la or not rb or la == "NAN" or rb == "NAN":
            continue
        rows.append(
            {
                "LIGAND": la,
                "RECEPTOR": rb,
                "PAIR": f"{la}→{rb}",
                "DATABASE": "CellPhoneDB_simple",
                "CLASSIFICATION": r.get("classification", ""),
            }
        )
    if not rows:
        return pd.DataFrame(columns=["LIGAND", "RECEPTOR", "PAIR", "DATABASE"])
    out = pd.DataFrame(rows).drop_duplicates(["LIGAND", "RECEPTOR"])
    return out


def fisher_enrich(
    hits: set[str],
    background: set[str],
    gene_set: set[str],
    label: str,
    cell_type: str,
    direction: str,
) -> dict:
    bg = {g.upper() for g in background}
    hit = {g.upper() for g in hits} & bg
    s = {g.upper() for g in gene_set} & bg
    n_bg, n_hit, k = len(bg), len(hit), len(s)
    if n_hit < 5 or k < 5:
        return {
            "CELL_TYPE": cell_type,
            "DIRECTION": direction,
            "SET": label,
            "OVERLAP": 0,
            "HIT_SIZE": n_hit,
            "SET_SIZE_IN_BG": k,
            "BG_SIZE": n_bg,
            "ODDS_RATIO": np.nan,
            "P_VALUE": np.nan,
            "GENES": "",
        }
    a = len(hit & s)
    table = np.array([[a, n_hit - a], [k - a, n_bg - n_hit - (k - a)]])
    if table.min() < 0:
        oddsr, p = np.nan, np.nan
    else:
        oddsr, p = fisher_exact(table, alternative="greater")
    return {
        "CELL_TYPE": cell_type,
        "DIRECTION": direction,
        "SET": label,
        "OVERLAP": int(a),
        "HIT_SIZE": n_hit,
        "SET_SIZE_IN_BG": int(k),
        "BG_SIZE": n_bg,
        "ODDS_RATIO": float(oddsr) if np.isfinite(oddsr) else np.nan,
        "P_VALUE": float(p) if np.isfinite(p) else np.nan,
        "GENES": ";".join(sorted(hit & s)),
    }


def prepare_ranks(ranked: pd.DataFrame) -> pd.DataFrame:
    d = ranked.copy()
    d["SYMBOL"] = d["SYMBOL"].fillna(d.get("SYMBOL_PLOT", d["GENE"])).astype(str).str.upper()
    d = d[d["SYMBOL"].notna() & (d["SYMBOL"] != "") & (d["SYMBOL"] != "NAN")]
    d = d[np.isfinite(d["SCORE"])]
    # collapse duplicate symbols within cell type by max |score|
    d = (
        d.assign(ABS=d["SCORE"].abs())
        .sort_values("ABS", ascending=False)
        .drop_duplicates(["CELL_TYPE", "SYMBOL"])
    )
    return d.drop(columns=["ABS"])


def pair_scores_within_celltype(
    ranks_ct: pd.DataFrame,
    lr: pd.DataFrame,
    min_abs: float,
) -> pd.DataFrame:
    score = ranks_ct.set_index("SYMBOL")["SCORE"]
    delta = ranks_ct.set_index("SYMBOL")["CLIFFS_DELTA"] if "CLIFFS_DELTA" in ranks_ct.columns else None
    rows = []
    for _, e in lr.iterrows():
        L, R = e["LIGAND"], e["RECEPTOR"]
        if L not in score.index or R not in score.index:
            continue
        sL, sR = float(score[L]), float(score[R])
        if abs(sL) < min_abs and abs(sR) < min_abs:
            continue
        # pair priority: both high |score|; prefer same-sign concordance or L+/R+ etc.
        pair_score = sL * sR  # positive if same direction
        pair_strength = abs(sL) + abs(sR)
        rows.append(
            {
                "LIGAND": L,
                "RECEPTOR": R,
                "PAIR": e["PAIR"],
                "DATABASE": e.get("DATABASE", ""),
                "SCORE_L": sL,
                "SCORE_R": sR,
                "DELTA_L": float(delta[L]) if delta is not None else np.nan,
                "DELTA_R": float(delta[R]) if delta is not None else np.nan,
                "PAIR_PRODUCT": pair_score,
                "PAIR_STRENGTH": pair_strength,
                "SAME_SIGN": np.sign(sL) == np.sign(sR) and sL != 0 and sR != 0,
            }
        )
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).sort_values("PAIR_STRENGTH", ascending=False)


def pair_scores_cross(
    ranks: pd.DataFrame,
    lr: pd.DataFrame,
    sender: str,
    receiver: str,
    min_abs: float,
) -> pd.DataFrame:
    """Ligand ranked in sender cell type; receptor ranked in receiver cell type."""
    s_send = ranks[ranks["CELL_TYPE"] == sender].set_index("SYMBOL")["SCORE"]
    s_recv = ranks[ranks["CELL_TYPE"] == receiver].set_index("SYMBOL")["SCORE"]
    rows = []
    for _, e in lr.iterrows():
        L, R = e["LIGAND"], e["RECEPTOR"]
        if L not in s_send.index or R not in s_recv.index:
            continue
        sL, sR = float(s_send[L]), float(s_recv[R])
        if abs(sL) < min_abs or abs(sR) < min_abs:
            continue
        rows.append(
            {
                "SENDER_CELL": sender,
                "RECEIVER_CELL": receiver,
                "LIGAND": L,
                "RECEPTOR": R,
                "PAIR": e["PAIR"],
                "DATABASE": e.get("DATABASE", ""),
                "SCORE_L_SENDER": sL,
                "SCORE_R_RECEIVER": sR,
                "PAIR_STRENGTH": abs(sL) + abs(sR),
                "PAIR_PRODUCT": sL * sR,
            }
        )
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).sort_values("PAIR_STRENGTH", ascending=False)


def plot_enrichment_heatmap(enr: pd.DataFrame, out_pdf: Path, out_png: Path) -> None:
    d = enr.copy()
    d = d[np.isfinite(d["P_VALUE"])]
    if d.empty:
        return
    d["NEGLOG10"] = -np.log10(np.clip(d["P_VALUE"], 1e-300, None))
    # facet-like: SET × DIRECTION as columns
    d["COL"] = d["SET"] + " | " + d["DIRECTION"]
    mat = d.pivot_table(index="CELL_TYPE", columns="COL", values="NEGLOG10", aggfunc="max")
    # order cell types
    order = [c for c in KEY_CELL_TYPES if c in mat.index] + sorted(
        set(mat.index) - set(KEY_CELL_TYPES)
    )
    mat = mat.reindex(order)
    fig_h = max(4.0, 0.35 * len(mat) + 1.5)
    fig_w = max(6.0, 0.55 * mat.shape[1] + 2)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), layout="constrained")
    im = ax.imshow(mat.to_numpy(), aspect="auto", cmap="YlOrRd", vmin=0)
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels(mat.columns, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(mat.shape[0]))
    ax.set_yticklabels(mat.index, fontsize=8)
    ax.set_title("Ligand / receptor enrichment in splicing-outlier score tails\n(−log10 Fisher P)", fontsize=10)
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label=r"$-\log_{10}P$")
    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=300, facecolor="white")
    plt.close(fig)


def plot_top_pairs_heatmap(pairs: pd.DataFrame, out_pdf: Path, out_png: Path, top_n: int = 40) -> None:
    if pairs is None or pairs.empty:
        return
    d = pairs.nlargest(top_n, "PAIR_STRENGTH").copy()
    d = d.sort_values("PAIR_STRENGTH")
    fig, ax = plt.subplots(figsize=(8.5, max(4.0, 0.28 * len(d) + 1.2)), layout="constrained")
    y = np.arange(len(d))
    colors = np.where(d["PAIR_PRODUCT"] >= 0, "#B2182B", "#2166AC")
    ax.barh(y, d["PAIR_STRENGTH"], color=colors, alpha=0.85, height=0.7)
    labels = [
        f"{r.CELL_TYPE} | {r.PAIR}"
        for r in d.itertuples()
    ]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("|score_L| + |score_R|")
    ax.set_title("Top within–cell-type LR pairs by splicing-outlier rank strength", fontsize=10)
    ax.axvline(0, color="#666", lw=0.5)
    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=300, facecolor="white")
    plt.close(fig)


def plot_cross_chord_style(cross: pd.DataFrame, out_pdf: Path, out_png: Path, top_n: int = 30) -> None:
    """Simple sender→receiver edge plot (not a full chord diagram)."""
    if cross is None or cross.empty:
        return
    d = cross.nlargest(top_n, "PAIR_STRENGTH").copy()
    # aggregate edge weight by sender-receiver
    edges = (
        d.groupby(["SENDER_CELL", "RECEIVER_CELL"], as_index=False)
        .agg(N_PAIRS=("PAIR", "nunique"), MEAN_STRENGTH=("PAIR_STRENGTH", "mean"), MAX_STRENGTH=("PAIR_STRENGTH", "max"))
        .sort_values("MAX_STRENGTH", ascending=False)
    )
    nodes = sorted(set(edges["SENDER_CELL"]) | set(edges["RECEIVER_CELL"]))
    pos = {n: i for i, n in enumerate(nodes)}
    fig, axes = plt.subplots(1, 2, figsize=(11, max(4.0, 0.35 * len(nodes) + 1.5)), layout="constrained")

    ax = axes[0]
    for _, r in edges.iterrows():
        y0, y1 = pos[r["SENDER_CELL"]], pos[r["RECEIVER_CELL"]]
        ax.plot([0, 1], [y0, y1], color="#555", lw=0.5 + 1.5 * r["MAX_STRENGTH"] / edges["MAX_STRENGTH"].max(), alpha=0.7)
    ax.set_yticks(range(len(nodes)))
    ax.set_yticklabels(nodes, fontsize=8)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Sender\n(ligand rank)", "Receiver\n(receptor rank)"], fontsize=8)
    ax.set_title("Cross–cell-type LR links", fontsize=10)
    ax.set_xlim(-0.1, 1.1)

    ax2 = axes[1]
    show = d.nlargest(min(top_n, 25), "PAIR_STRENGTH").sort_values("PAIR_STRENGTH")
    y = np.arange(len(show))
    ax2.barh(y, show["PAIR_STRENGTH"], color="#542788", alpha=0.85)
    ax2.set_yticks(y)
    ax2.set_yticklabels(
        [f"{r.SENDER_CELL}→{r.RECEIVER_CELL} | {r.PAIR}" for r in show.itertuples()],
        fontsize=6.5,
    )
    ax2.set_xlabel("|score_L(sender)| + |score_R(receiver)|")
    ax2.set_title("Top cross LR pairs", fontsize=10)
    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=300, facecolor="white")
    plt.close(fig)


def main() -> int:
    root = project_root()
    os.chdir(root)
    method = os.environ.get("DECONV_METHOD", "bayesprism").lower()
    top_n = int(os.environ.get("LR_TOP_N", str(DEFAULT_TOP_N)))
    pair_top = int(os.environ.get("LR_PAIR_TOP", str(DEFAULT_PAIR_TOP)))
    min_abs = float(os.environ.get("LR_MIN_ABS_SCORE", str(DEFAULT_MIN_ABS)))

    ref = root / "CLLCMA_alpha/references/ligand_receptor"
    op_path = ref / "omnipath_ligrecextra.tsv"
    if not op_path.exists():
        raise SystemExit(
            f"Missing {op_path}. Download OmniPath first (user-approved script)."
        )

    base = (
        root
        / "CLLCMA_alpha/results/deconvolution/cll_broad_2022"
        / "splicing_outlier_tme_association"
        / method
        / "rank_enrichment"
    )
    rank_path = base / "gene_scores_ranked_by_celltype.tsv.gz"
    if not rank_path.exists():
        raise SystemExit(f"Missing ranked scores: {rank_path}")

    out_dir = base / "ligand_receptor_prioritization"
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Loading OmniPath LR …", flush=True)
    lr_op = load_omnipath_lr(op_path)
    lr_parts = [lr_op]
    cp_gene = ref / "cellphone_gene_input.csv"
    cp_int = ref / "cellphone_interaction_input.csv"
    if cp_gene.exists() and cp_int.exists():
        print("Loading CellPhoneDB simple pairs …", flush=True)
        lr_cp = load_cellphone_simple_pairs(cp_gene, cp_int)
        print(f"  CellPhoneDB simple pairs: {len(lr_cp)}", flush=True)
        lr_parts.append(lr_cp[["LIGAND", "RECEPTOR", "PAIR", "DATABASE"]])
    lr = pd.concat(lr_parts, ignore_index=True).drop_duplicates(["LIGAND", "RECEPTOR"])
    lr.to_csv(out_dir / "lr_pairs_combined.tsv.gz", sep="\t", index=False, compression="gzip")

    ligands = set(lr["LIGAND"])
    receptors = set(lr["RECEPTOR"])
    lr_genes = ligands | receptors
    print(
        f"LR database: {len(lr)} unique pairs; "
        f"{len(ligands)} ligands; {len(receptors)} receptors",
        flush=True,
    )

    print("Loading gene ranks …", flush=True)
    ranked = prepare_ranks(pd.read_csv(rank_path, sep="\t"))
    # prefer all_types_union if present
    if "AS_TYPE" in ranked.columns and (ranked["AS_TYPE"] == "all_types_union").any():
        ranked = ranked[ranked["AS_TYPE"] == "all_types_union"].copy()

    # ---- ORA enrichment ----
    enr_rows = []
    for ct, sub in ranked.groupby("CELL_TYPE"):
        bg = set(sub["SYMBOL"])
        pos = set(sub.nlargest(top_n, "SCORE")["SYMBOL"])
        neg = set(sub.nsmallest(top_n, "SCORE")["SYMBOL"])
        for direction, hits in (("positive", pos), ("negative", neg)):
            enr_rows.append(fisher_enrich(hits, bg, ligands, "Ligand", ct, direction))
            enr_rows.append(fisher_enrich(hits, bg, receptors, "Receptor", ct, direction))
            enr_rows.append(fisher_enrich(hits, bg, lr_genes, "LR_gene", ct, direction))
    enr = pd.DataFrame(enr_rows)
    enr["FDR"] = bh_fdr(enr["P_VALUE"].to_numpy())
    enr = enr.sort_values(["P_VALUE", "CELL_TYPE"])
    enr.to_csv(out_dir / "lr_set_enrichment_by_celltype.tsv.gz", sep="\t", index=False, compression="gzip")
    enr_sig = enr[enr["FDR"] < 0.05].copy()
    enr_sig.to_csv(out_dir / "lr_set_enrichment_fdr05.tsv.gz", sep="\t", index=False, compression="gzip")
    print(f"Enrichment tests: {len(enr)}; FDR<0.05: {len(enr_sig)}", flush=True)

    plot_enrichment_heatmap(
        enr,
        out_dir / "lr_enrichment_heatmap.pdf",
        out_dir / "lr_enrichment_heatmap.png",
    )

    # ---- Within–cell-type pair prioritization ----
    within_rows = []
    for ct, sub in ranked.groupby("CELL_TYPE"):
        p = pair_scores_within_celltype(sub, lr, min_abs=min_abs)
        if p.empty:
            continue
        p.insert(0, "CELL_TYPE", ct)
        within_rows.append(p.head(pair_top * 2))
    within = pd.concat(within_rows, ignore_index=True) if within_rows else pd.DataFrame()
    if len(within):
        within.to_csv(
            out_dir / "lr_pairs_within_celltype.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )
        # top overall for figure
        plot_top_pairs_heatmap(
            within,
            out_dir / "lr_pairs_within_top_heatmap.pdf",
            out_dir / "lr_pairs_within_top_heatmap.png",
            top_n=40,
        )
        print(f"Within-cell-type pairs scored: {len(within)}", flush=True)

    # ---- Cross–cell-type (selected senders/receivers) ----
    senders = [c for c in ["CLL", "ASDC", "Monocyte_CD14", "Monocyte_CD16", "cDC1", "cDC2"] if c in ranked["CELL_TYPE"].unique()]
    receivers = [
        c
        for c in [
            "CD8aa",
            "CD8_T_naive",
            "CD8_T_memory",
            "CD4_T_naive",
            "T_proliferating",
            "Treg",
            "NK_proliferating",
            "Monocyte_CD16",
            "CLL",
        ]
        if c in ranked["CELL_TYPE"].unique()
    ]
    cross_rows = []
    for s in senders:
        for r in receivers:
            if s == r:
                continue
            p = pair_scores_cross(ranked, lr, s, r, min_abs=min_abs)
            if p.empty:
                continue
            cross_rows.append(p.head(pair_top))
    cross = pd.concat(cross_rows, ignore_index=True) if cross_rows else pd.DataFrame()
    if len(cross):
        cross.to_csv(
            out_dir / "lr_pairs_cross_celltype.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )
        plot_cross_chord_style(
            cross,
            out_dir / "lr_pairs_cross_summary.pdf",
            out_dir / "lr_pairs_cross_summary.png",
            top_n=40,
        )
        print(f"Cross-cell-type pairs scored: {len(cross)}", flush=True)

    # ---- KEY FINDINGS ----
    lines = [
        "SPLICING-OUTLIER RANK × LIGAND–RECEPTOR PRIORITIZATION (Idea A)",
        f"Method: {method}",
        f"Ranks: {rank_path.name} (AS_TYPE filtered to all_types_union if present)",
        f"LR DB: OmniPath ligrecextra"
        + (" + CellPhoneDB simple pairs" if cp_gene.exists() else ""),
        f"Unique LR pairs: {len(lr)}; ligands={len(ligands)}; receptors={len(receptors)}",
        f"ORA tails: top/bottom {top_n} genes by score",
        f"Pair min |score|: {min_abs}",
        f"Enrichment FDR<0.05: {len(enr_sig)} / {len(enr)}",
        "",
        "=== Top enrichment hits ===",
    ]
    for _, r in enr.nsmallest(20, "P_VALUE").iterrows():
        lines.append(
            f"{r['CELL_TYPE']} | {r['DIRECTION']} | {r['SET']} | "
            f"overlap={r['OVERLAP']}/{r['HIT_SIZE']} OR={r['ODDS_RATIO']:.2f} "
            f"P={r['P_VALUE']:.3g} FDR={r['FDR']:.3g}"
        )
    if len(within):
        lines.append("")
        lines.append("=== Top within–cell-type pairs ===")
        for _, r in within.nlargest(15, "PAIR_STRENGTH").iterrows():
            lines.append(
                f"{r['CELL_TYPE']} | {r['PAIR']} | "
                f"sL={r['SCORE_L']:+.2f} sR={r['SCORE_R']:+.2f} "
                f"strength={r['PAIR_STRENGTH']:.2f}"
            )
    if len(cross):
        lines.append("")
        lines.append("=== Top cross–cell-type pairs ===")
        for _, r in cross.nlargest(15, "PAIR_STRENGTH").iterrows():
            lines.append(
                f"{r['SENDER_CELL']}→{r['RECEIVER_CELL']} | {r['PAIR']} | "
                f"sL={r['SCORE_L_SENDER']:+.2f} sR={r['SCORE_R_RECEIVER']:+.2f} "
                f"strength={r['PAIR_STRENGTH']:.2f}"
            )
    lines.append("")
    lines.append(
        "Interpretation: correlative prioritization only; bulk ranks + LR databases "
        "do not prove cell–cell signaling."
    )
    (out_dir / "KEY_FINDINGS.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[:35]), flush=True)
    print(f"Done -> {out_dir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
