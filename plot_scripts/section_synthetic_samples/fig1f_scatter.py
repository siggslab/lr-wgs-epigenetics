"""Figure 1F: recovery of known cell-type proportions in synthetic mixtures.

Reads the leave-one-out deconvolution table and plots, per compartment,
the deconvolution estimate (y) against the expected proportion (x) of each
synthetic sample.

The x-axis is NOT taken from the table; it is derived. Each deficiency sample
is built by removing one cell type's reads and leaving the other populations at
their original read counts, so the remaining cells' expected proportions are
just the base weightings reweighted for the missing cell:

    expected(cell) = weight(cell) / (1 - weight(dropped)),   0 for the dropped cell.

Worked example (synthetic blood, "w/out CD4s"): neutrophils have a base
fraction of 0.5. Removing CD4 T cells (base 0.175) leaves 1 - 0.175 = 0.825 of
the reads, so neutrophils now make up 0.5 / 0.825 = 0.606 of that sample -
they were not added to, they are simply a larger slice of a smaller pool.

y (the deconvolution estimate) is taken directly from the table. Points on the
1:1 line mean the method recovered the known composition; the deficient-lineage
diamonds sit near the origin (expected ~= 0, recovered ~= 0).

Usage:  python3 fig1f_scatter_final.py [input.csv] [output.png]

------------------------------------------------------------------------------
FIGURE LEGEND (encodings, for the manuscript caption)
------------------------------------------------------------------------------
Two panels: synthetic blood (6 cell types) and synthetic PBMCs (5 cell types).
Each point is one cell type in one synthetic sample.
  x-axis  Expected fraction: the known mixing proportion of that cell type in
          that sample, i.e. the base weighting reweighted for the removed
          lineage, weight / (1 - weight_removed).
  y-axis  Deconvolution estimate: the methylation-based predicted proportion.
  dashed grey line  1:1 identity (perfect recovery).
  circles  cell types present in the sample.
  diamonds (dark outline)  the deficient lineage in each leave-one-out sample
           (expected fraction ~= 0); their clustering at the origin shows the
           method correctly assigns ~0 to an absent cell type.
  r, RMSE  Pearson correlation and root-mean-square error across all points in
           that panel.
Colours match Figure 1 panels B/D/E (COL dict below; its key order sets the
legend order):
  Monocytes blue, B-cells orange, CD4 T green, NK red, CD8 T purple,
  Neutrophils brown.

Suggested caption sentence:
  "F) Deconvolution estimates versus known (reweighted) cell-type proportions
  for each synthetic blood and PBMC mixture. Diamonds mark the depleted lineage
  in each leave-one-out sample (expected proportion ~= 0). Colours as in B-E;
  dashed line, 1:1 identity."
------------------------------------------------------------------------------
"""
import csv
import sys
import numpy as np
import matplotlib.pyplot as plt
CSV = sys.argv[1] if len(sys.argv) > 1 else "fig1f_leave_one_out.csv"
OUT = sys.argv[2] if len(sys.argv) > 2 else "fig1f_scatter.png"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 12, "axes.linewidth": 0.8,
})

# Matched to Figure 1 panels B/D/E (matplotlib default cycle, legend order:
# Monocytes, B-cells, CD4, NK, CD8, Neutrophils). Keys match CSV headers.
COL = {"Monocytes": "#1f77b4", "B Cells": "#ff7f0e", "CD4 T Cells": "#2ca02c",
       "NK Cells": "#d62728", "CD8 T Cells": "#9467bd", "Neutrophils": "#8c564b"}
SHORT = {"Neutrophils": "Neutrophils", "Monocytes": "Monocytes",
         "CD4 T Cells": "CD4 T", "CD8 T Cells": "CD8 T",
         "B Cells": "B", "NK Cells": "NK"}
# keyword in a "w/out ..." label  ->  column it refers to
DROP_KEY = {"neutrophil": "Neutrophils", "monocyte": "Monocytes",
            "cd4": "CD4 T Cells", "cd8": "CD8 T Cells",
            "b cell": "B Cells", "nk": "NK Cells"}


def dropped_cell(label):
    if "w/out" not in label.lower():
        return None
    low = label.lower()
    for key, col in DROP_KEY.items():
        if key in low:
            return col
    raise ValueError(f"cannot parse dropped cell from {label!r}")


def parse(path):
    """Return list of blocks; each block = (cells, weighting, [(label, dropped, preds)])."""
    with open(path) as fh:
        rows = [r for r in csv.reader(fh) if any(c.strip() for c in r)]
    header = rows[0][1:]
    blocks, cur = [], None
    for r in rows[1:]:
        label, vals = r[0].strip(), r[1:]
        present = [c for c, v in zip(header, vals) if v.strip() != ""]
        num = {c: float(v) for c, v in zip(header, vals) if v.strip() != ""}
        if label.upper() == "WEIGHTING":
            cur = {"cells": present, "W": num, "samples": []}
            blocks.append(cur)
        else:
            cur["samples"].append((label, dropped_cell(label), num))
    return blocks


def expected(W, dropped):
    if dropped is None:
        return dict(W)
    denom = 1.0 - W[dropped]
    return {c: (0.0 if c == dropped else w / denom) for c, w in W.items()}


def points(block):
    pts = []  # (cell, expected, predicted, is_deficient)
    for _, dropped, pred in block["samples"]:
        exp = expected(block["W"], dropped)
        for c in block["cells"]:
            pts.append((c, exp[c], pred[c], c == dropped))
    return pts


def draw(ax, pts, title):
    ax.plot([0, 0.62], [0, 0.62], color="#bbb", lw=1, ls=(0, (4, 4)), zorder=1)
    for c, e, p, defc in pts:
        ax.scatter(e, p, s=70, color=COL[c], zorder=3,
                   marker="D" if defc else "o",
                   edgecolor="#333" if defc else "white",
                   linewidth=1.4 if defc else 0.7, alpha=0.95)
    e = np.array([p[1] for p in pts]); p_ = np.array([p[2] for p in pts])
    rmse = np.sqrt(np.mean((e - p_) ** 2))
    r = np.corrcoef(e, p_)[0, 1]
    ax.set_xlim(-0.02, 0.62); ax.set_ylim(-0.02, 0.62); ax.set_aspect("equal")
    ax.set_xlabel("Expected fraction (known synthetic mixture)")
    ax.set_title(title, fontweight="bold")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.text(0.03, 0.58, f"$r$ = {r:.3f}\nRMSE = {rmse:.3f}", fontsize=12, va="top")
    return r, rmse

if __name__ == '__main__':
    blocks = parse(CSV)
    titles = {"Blood": "Synthetic blood", "PBMCs": "Synthetic PBMCs"}
    # name each block by its full-reference (dropped is None) sample label
    named = []
    for b in blocks:
        ref = next(lbl for lbl, d, _ in b["samples"] if d is None)
        named.append((titles.get(ref, ref), b))

    fig, axes = plt.subplots(1, len(named), figsize=(5.25 * len(named), 5.4),
                             gridspec_kw={"wspace": 0.12})
    axes = np.atleast_1d(axes)
    for ax, (title, b) in zip(axes, named):
        draw(ax, points(b), title)
    axes[0].set_ylabel("Deconvolution estimate")
    for ax in axes[1:]:
        ax.tick_params(labelleft=False)

    cells_order = list(COL.keys())
    handles = [plt.Line2D([], [], marker="o", ls="", mfc=COL[c], mec="white",
                          ms=8, label=SHORT[c]) for c in cells_order]
    handles.append(plt.Line2D([], [], marker="D", ls="", mfc="#999", mec="#333",
                              ms=8, label="Deficient lineage\n(expected ≈ 0)"))
    legend = axes[-1].legend(handles=handles, loc="lower right", fontsize=10,
                    handletextpad=0.4, labelspacing=0.6)
    legend.get_frame().set_facecolor("#e8e8e8")
    legend.get_frame().set_alpha(0.7)

    fig.savefig(OUT, dpi=300, bbox_inches="tight", facecolor="white")
    print("saved", OUT)
