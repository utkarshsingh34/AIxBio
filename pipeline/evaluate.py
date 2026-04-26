"""
evaluate.py

Four-figure analysis with mutation-rate sweep support:

  FIGURE 1 (per-fragment, the fundamental measurement):
      For each fragment length T and mutation rate M, what fraction of
      HAZARDOUS fragments and HONEST fragments produce a BLAST hit at
      e-value <= cutoff?

  FIGURE 2 (per-order, the deployment scenario):
      ROC panels at each (T, mutation rate). Per-order score = min e-value
      across the order's fragments (any-flag aggregation).

  FIGURE 3 (NEW, the headline visual):
      Two-panel heatmap of TPR @ FPR=1% across (T, mutation rate). Left =
      pure evasion, right = dilute evasion. The "evasion gap" between panels
      is the visible argument that dilute attacks defeat any-flag aggregation
      where pure attacks don't.

  FIGURE 4 (NEW, the operating-point visual):
      Bar chart comparing TPR @ FPR=1% at three representative scenarios
      (naive / realistic / sophisticated adversary), at the OSTP-mandated
      50 bp threshold. Designed for slides.

USAGE (from your pipeline folder, with blast-env active):
    python evaluate.py
    python evaluate.py --T 20,30,50,75,100,150,200
    python evaluate.py --evalue-cutoffs 1e-3,1e-9,1e-30

Outputs (to ../results/):
    fig1_fragment_sensitivity.png  - per-fragment TPR/FPR vs T x mutation
    fig2_order_roc_panels.png      - per-order ROC, full grid
    fig3_evasion_heatmap.png       - HEADLINE: pure vs dilute heatmaps
    fig4_operating_points.png      - SLIDE: bar chart at T=50
    fragment_sensitivity.csv       - tidy data behind Fig 1
    roc_data.csv                   - tidy data behind Fig 2
    summary.csv                    - per-(T, mut, evasion type) AUC + TPR @ FPR
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


DIR_RE = re.compile(r"^T(?P<T>\d+)(?:_mut(?P<mut>\d+))?$")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--corpora-root", default="../data/corpora")
    p.add_argument("--T", default=None)
    p.add_argument("--mutation-rates", default=None)
    p.add_argument("--out-dir", default="../results")
    p.add_argument("--fpr-targets", default="0.01,0.05,0.10")
    p.add_argument("--evalue-cutoffs", default="1e-3,1e-9,1e-30")
    p.add_argument("--ostp-T", type=int, default=50,
                   help="The OSTP-mandated screening threshold to highlight in fig4 (default: 50)")
    return p.parse_args()


def discover_corpus_dirs(root):
    out = {}
    for p in sorted(root.iterdir()):
        if not p.is_dir():
            continue
        m = DIR_RE.match(p.name)
        if not m:
            continue
        T = int(m.group("T"))
        mut = int(m.group("mut")) if m.group("mut") is not None else 0
        out[(T, mut)] = p
    return out


# -----------------------------------------------------------------------------
# FIGURE 1: per-fragment sensitivity
# -----------------------------------------------------------------------------

def per_fragment_sensitivity(corpus_dirs, evalue_cutoffs):
    rows = []
    for (T, mut), T_dir in corpus_dirs.items():
        for label, corpus_name in [("hazardous", "evasion_pure_fragments.csv"),
                                   ("honest", "honest_fragments.csv")]:
            csv = T_dir / corpus_name
            if not csv.exists():
                continue
            frags = pd.read_csv(csv)
            if label == "hazardous":
                frags = frags[frags["truth_hazard"] == True]
            else:
                frags = frags[frags["truth_hazard"] == False]
            if frags.empty:
                continue
            mev = frags["min_evalue"].fillna(np.inf).values
            for cutoff in evalue_cutoffs:
                n_hits = int(np.sum(mev <= cutoff))
                rows.append({
                    "T": T, "mut_pct": mut, "label": label, "cutoff": cutoff,
                    "n": len(frags), "n_hits": n_hits,
                    "detection_rate": n_hits / len(frags),
                })
    return pd.DataFrame(rows)


def plot_figure_1(sens_df, out_path):
    cutoffs = sorted(sens_df["cutoff"].unique(), reverse=True)
    muts = sorted(sens_df["mut_pct"].unique())

    nrows = len(cutoffs)
    fig, axes = plt.subplots(nrows, 2, figsize=(12, 4 * nrows),
                             sharex=True, sharey=True, squeeze=False)

    cmap = plt.get_cmap("viridis")
    color_for = {m: cmap(i / max(1, len(muts) - 1)) for i, m in enumerate(muts)}

    for row_i, cutoff in enumerate(cutoffs):
        for col_i, (label, ylabel) in enumerate([
            ("hazardous", "TPR (hazardous fragments hit)"),
            ("honest", "FPR (honest fragments hit)"),
        ]):
            ax = axes[row_i][col_i]
            sub = sens_df[(sens_df["cutoff"] == cutoff) & (sens_df["label"] == label)]
            for mut in muts:
                row = sub[sub["mut_pct"] == mut].sort_values("T")
                if row.empty:
                    continue
                ax.plot(row["T"], row["detection_rate"],
                        marker="o", linewidth=2,
                        color=color_for[mut],
                        label=f"{mut}% mutation")
            ax.set_xscale("log")
            ax.set_ylim(-0.02, 1.05)
            ax.grid(alpha=0.3)
            if row_i == nrows - 1:
                ax.set_xlabel("Fragment length T (bp)")
            ax.set_ylabel(ylabel)
            ax.set_title(f"{'Hazardous' if label == 'hazardous' else 'Honest (E. coli)'}"
                         f" — cutoff e ≤ {cutoff:g}")
            if col_i == 1 and row_i == 0:
                ax.legend(title="Mutation rate", loc="best", fontsize=9)

    fig.suptitle("Figure 1 — BLAST per-fragment sensitivity vs fragment length and mutation rate",
                 fontsize=12, y=1.00)
    plt.tight_layout()
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)


# -----------------------------------------------------------------------------
# FIGURE 2: per-order ROC
# -----------------------------------------------------------------------------

def compute_roc(honest_scores, evasion_scores):
    def to_score(s):
        s = np.asarray(s, dtype=float)
        with np.errstate(divide="ignore"):
            sc = -np.log10(np.where(s > 0, s, 1e-300))
        sc[np.isinf(s)] = -np.inf
        return sc
    h = to_score(honest_scores)
    e = to_score(evasion_scores)
    all_scores = np.concatenate([h, e])
    finite = all_scores[np.isfinite(all_scores)]
    if finite.size == 0:
        return pd.DataFrame({"threshold": [np.inf], "fpr": [0.0], "tpr": [0.0]})
    thresholds = np.unique(np.concatenate([[-np.inf], np.sort(finite), [np.inf]]))
    rows = []
    for t in thresholds[::-1]:
        rows.append({"threshold": float(t),
                     "fpr": float(np.mean(h >= t)),
                     "tpr": float(np.mean(e >= t))})
    return pd.DataFrame(rows)


def auc_trapezoid(fpr, tpr):
    order = np.argsort(fpr)
    return float(np.trapezoid(np.asarray(tpr)[order], np.asarray(fpr)[order]))


def tpr_at_fpr(roc_df, target_fpr):
    eligible = roc_df[roc_df["fpr"] <= target_fpr]
    return float(eligible["tpr"].max()) if not eligible.empty else 0.0


def plot_figure_2(roc_data, T_values, mut_values, out_path):
    nrows = len(T_values)
    ncols = len(mut_values)
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.5 * ncols, 3.5 * nrows),
                             squeeze=False)
    color_for = {"evasion_pure": "tab:blue", "evasion_dilute": "tab:orange"}

    for i, T in enumerate(T_values):
        for j, mut in enumerate(mut_values):
            ax = axes[i][j]
            for ev_type in ["evasion_pure", "evasion_dilute"]:
                sub = roc_data[(roc_data["T"] == T) &
                               (roc_data["mut_pct"] == mut) &
                               (roc_data["evasion_type"] == ev_type)]
                if sub.empty:
                    continue
                sub_sorted = sub.sort_values("fpr")
                auc = auc_trapezoid(sub_sorted["fpr"].values, sub_sorted["tpr"].values)
                ax.plot(sub_sorted["fpr"], sub_sorted["tpr"],
                        color=color_for[ev_type], linewidth=2,
                        label=f"{ev_type.replace('evasion_', '')} ({auc:.2f})")
            ax.plot([0, 1], [0, 1], "k--", alpha=0.3, linewidth=1)
            ax.set_xlim(0, 1); ax.set_ylim(0, 1.02)
            if i == nrows - 1:
                ax.set_xlabel("FPR")
            if j == 0:
                ax.set_ylabel(f"T={T}\nTPR")
            if i == 0:
                ax.set_title(f"Mutation = {mut}%")
            ax.legend(loc="lower right", fontsize=7)
            ax.grid(alpha=0.3)

    fig.suptitle("Figure 2 — Per-order ROC under any-flag aggregation (rows: T, cols: mutation)",
                 fontsize=12, y=1.00)
    plt.tight_layout()
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)


# -----------------------------------------------------------------------------
# FIGURE 3: HEADLINE heatmap (NEW)
# -----------------------------------------------------------------------------

def plot_figure_3(summary_df, out_path, metric="tpr_at_fpr0.01"):
    """Two-panel heatmap: pure (left) vs dilute (right). T on Y, mutation on X.
    Color = TPR at FPR=1%. The visible 'evasion gap' between panels IS the headline."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), sharey=True)

    ev_types = [("evasion_pure", "Pure evasion"),
                ("evasion_dilute", "Dilute evasion (1 hazardous + benign filler)")]

    for ax, (ev_type, title) in zip(axes, ev_types):
        sub = summary_df[summary_df["evasion_type"] == ev_type]
        if sub.empty:
            continue
        pivot = sub.pivot_table(index="T", columns="mut_pct", values=metric)
        pivot = pivot.sort_index(ascending=False)  # large T at top

        im = ax.imshow(pivot.values, aspect="auto", cmap="RdYlGn",
                       vmin=0, vmax=1, interpolation="nearest")

        # Annotate every cell
        for i in range(pivot.shape[0]):
            for j in range(pivot.shape[1]):
                v = pivot.values[i, j]
                # Choose readable text color
                txt_color = "white" if (v < 0.35 or v > 0.85) else "black"
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        color=txt_color, fontsize=9)

        ax.set_xticks(range(len(pivot.columns)))
        ax.set_xticklabels([f"{m}%" for m in pivot.columns])
        ax.set_yticks(range(len(pivot.index)))
        ax.set_yticklabels([f"{T}" for T in pivot.index])
        ax.set_xlabel("Adversary mutation rate")
        ax.set_ylabel("Fragment length T (bp)")
        ax.set_title(title)

    cbar = fig.colorbar(im, ax=axes, fraction=0.04, pad=0.02)
    cbar.set_label(f"TPR at FPR = 1%  (catch rate at strict threshold)")

    fig.suptitle("Figure 3 — The evasion gap: pure vs dilute attacks under any-flag aggregation",
                 fontsize=13, y=1.02)
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)


# -----------------------------------------------------------------------------
# FIGURE 4: SLIDE bar chart (NEW)
# -----------------------------------------------------------------------------

def plot_figure_4(summary_df, out_path, ostp_T=50):
    """Bar chart at T=ostp_T showing TPR@FPR=1% for representative adversaries."""
    sub = summary_df[summary_df["T"] == ostp_T].copy()
    if sub.empty:
        print(f"  [skip] fig4: no data at T={ostp_T}")
        return

    # Three representative adversary scenarios
    scenarios = [
        ("Naive\n(0% mutation)", 0),
        ("Realistic\n(10% mutation)", 10),
        ("Sophisticated\n(20% mutation)", 20),
    ]
    metrics = ["tpr_at_fpr0.01", "tpr_at_fpr0.05", "tpr_at_fpr0.1"]
    metric_labels = ["FPR ≤ 1%", "FPR ≤ 5%", "FPR ≤ 10%"]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    width = 0.25
    x = np.arange(len(scenarios))

    for ax, ev_type, panel_title in zip(
        axes, ["evasion_pure", "evasion_dilute"],
        ["Pure evasion", "Dilute evasion"]
    ):
        for i, (metric, mlabel) in enumerate(zip(metrics, metric_labels)):
            vals = []
            for _, mut in scenarios:
                row = sub[(sub["evasion_type"] == ev_type) & (sub["mut_pct"] == mut)]
                vals.append(float(row[metric].iloc[0]) if not row.empty else 0.0)
            offset = (i - 1) * width
            bars = ax.bar(x + offset, vals, width, label=mlabel,
                          edgecolor="black", linewidth=0.5)
            # Annotate bars
            for bar, v in zip(bars, vals):
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
                        f"{v:.2f}", ha="center", va="bottom", fontsize=9)

        ax.set_xticks(x)
        ax.set_xticklabels([s[0] for s in scenarios], fontsize=10)
        ax.set_ylim(0, 1.15)
        ax.set_ylabel("TPR (catch rate)")
        ax.set_title(f"{panel_title}\n(at T = {ostp_T} bp, the OSTP threshold)")
        ax.grid(alpha=0.3, axis="y")
        ax.legend(title="Operating point", loc="upper right", fontsize=9)
        ax.axhline(0.5, color="red", linestyle=":", alpha=0.5, linewidth=1)

    fig.suptitle(f"Figure 4 — Catch rates at the OSTP {ostp_T} bp threshold across adversary types",
                 fontsize=13, y=1.02)
    plt.tight_layout()
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    args = parse_args()
    root = Path(args.corpora_root).resolve()
    if not root.exists():
        sys.exit(f"Corpora root not found: {root}")

    corpus_dirs = discover_corpus_dirs(root)
    if not corpus_dirs:
        sys.exit(f"No T*/ or T*_mut*/ subdirs found under {root}.")

    if args.T:
        wanted_T = set(int(x) for x in args.T.split(",") if x.strip())
        corpus_dirs = {k: v for k, v in corpus_dirs.items() if k[0] in wanted_T}
    if args.mutation_rates is not None:
        wanted_mut = set(int(x) for x in args.mutation_rates.split(",") if x.strip())
        corpus_dirs = {k: v for k, v in corpus_dirs.items() if k[1] in wanted_mut}
    if not corpus_dirs:
        sys.exit("No matching corpus directories after filtering.")

    Ts = sorted({k[0] for k in corpus_dirs})
    muts = sorted({k[1] for k in corpus_dirs})
    print(f"[evaluate] T values: {Ts}")
    print(f"[evaluate] mutation rates: {muts}")

    fpr_targets = [float(x) for x in args.fpr_targets.split(",") if x.strip()]
    evalue_cutoffs = [float(x) for x in args.evalue_cutoffs.split(",") if x.strip()]

    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # ----- Figure 1 -----
    print(f"\n[fig1] computing per-fragment sensitivity")
    sens_df = per_fragment_sensitivity(corpus_dirs, evalue_cutoffs)
    if sens_df.empty:
        print("  [warn] no per-fragment data")
    else:
        sens_df.to_csv(out_dir / "fragment_sensitivity.csv", index=False)
        plot_figure_1(sens_df, out_dir / "fig1_fragment_sensitivity.png")
        print(f"  [done] fig1_fragment_sensitivity.png")

    # ----- Figure 2 -----
    print(f"\n[fig2] computing per-order ROC")
    all_roc = []
    summary_rows = []
    for (T, mut), T_dir in corpus_dirs.items():
        honest_csv = T_dir / "honest_orders.csv"
        if not honest_csv.exists():
            continue
        honest = pd.read_csv(honest_csv)
        honest_scores = honest["order_min_evalue"].values
        for ev_type in ["evasion_pure", "evasion_dilute"]:
            ev_csv = T_dir / f"{ev_type}_orders.csv"
            if not ev_csv.exists():
                continue
            evasion = pd.read_csv(ev_csv)
            evasion_scores = evasion["order_min_evalue"].values

            roc = compute_roc(honest_scores, evasion_scores)
            roc.insert(0, "T", T)
            roc.insert(1, "mut_pct", mut)
            roc.insert(2, "evasion_type", ev_type)
            all_roc.append(roc)

            row = {
                "T": T, "mut_pct": mut, "evasion_type": ev_type,
                "n_honest": len(honest_scores),
                "n_evasion": len(evasion_scores),
                "auc": auc_trapezoid(roc["fpr"].values, roc["tpr"].values),
            }
            for fpr_t in fpr_targets:
                row[f"tpr_at_fpr{fpr_t:g}"] = tpr_at_fpr(roc, fpr_t)
            summary_rows.append(row)

    summary = None
    if not all_roc:
        print("  [warn] no ROC data computed.")
    else:
        roc_data = pd.concat(all_roc, ignore_index=True)
        roc_data.to_csv(out_dir / "roc_data.csv", index=False)
        summary = pd.DataFrame(summary_rows)
        summary.to_csv(out_dir / "summary.csv", index=False)
        plot_figure_2(roc_data, Ts, muts, out_dir / "fig2_order_roc_panels.png")
        print(f"  [done] fig2_order_roc_panels.png")

    # ----- Figure 3 (NEW): headline heatmap -----
    if summary is not None:
        print(f"\n[fig3] plotting headline evasion-gap heatmap")
        plot_figure_3(summary, out_dir / "fig3_evasion_heatmap.png")
        print(f"  [done] fig3_evasion_heatmap.png")

    # ----- Figure 4 (NEW): operating-point bar chart -----
    if summary is not None:
        print(f"\n[fig4] plotting operating-point bars at T={args.ostp_T}")
        plot_figure_4(summary, out_dir / "fig4_operating_points.png", ostp_T=args.ostp_T)
        print(f"  [done] fig4_operating_points.png")

    # ----- Console outputs -----
    if not sens_df.empty:
        print("\n=== Per-fragment sensitivity (Figure 1 data) ===")
        wide = (sens_df.pivot_table(index=["T", "mut_pct", "label"],
                                    columns="cutoff",
                                    values="detection_rate")
                       .round(3))
        print(wide.to_string())

    if summary is not None:
        print("\n=== Per-order summary (Figure 2 data) ===")
        print(summary.round(3).to_string(index=False))

    print(f"\n[done] all outputs in {out_dir}")


if __name__ == "__main__":
    main()