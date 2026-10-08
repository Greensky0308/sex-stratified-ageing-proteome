#!/usr/bin/env python3
"""Figure S1 -- the circular 'sarcopenia index' (creatinine/cystatin C) at the CST3 locus.

a: regional association of the CST3 pQTL (SomaScan) and the sarcopenia index (SI)
   along chr20.  b: per-SNP coupling of the two z-scores (arithmetic, not biology).
Out: figures/FigS1.{pdf,tif}
"""
import csv, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(ROOT, "figures")
plt.rcParams.update({
    "font.family": "Arial", "font.size": 9, "axes.linewidth": 0.8,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "xtick.major.size": 3.5, "ytick.major.size": 3.5,
    "axes.labelsize": 9, "xtick.labelsize": 9, "ytick.labelsize": 9,
    "legend.fontsize": 9, "font.weight": "normal",
    "axes.spines.top": False, "axes.spines.right": False,
})
PC, FC, MC = "#2E8B8B", "#C05299", "#3E7CB1"
LEAD_PROT, LEAD_SI = 23.641629, 23.632100


def load(p):
    out = {}
    for row in csv.reader(open(p), delimiter="\t"):
        pos, b, se = int(row[1]), float(row[2]), float(row[3])
        if se > 0:
            out[row[0]] = (pos, b / se)
    return out


def main():
    prot = load(os.path.join(ROOT, "out/cst3_coloc_somascan.tsv"))
    siF = load(os.path.join(ROOT, "out/cst3_coloc_sarcF.tsv"))
    siM = load(os.path.join(ROOT, "out/cst3_coloc_sarcM.tsv"))

    fig = plt.figure(figsize=(7.6, 3.6))
    gs = fig.add_gridspec(1, 2, wspace=0.32)

    # ---- a: regional
    ax = fig.add_subplot(gs[0, 0])
    lo, hi = 23.55, 23.72
    ax.axvspan(23.6261, 23.6326, color="#EAF3EA", zorder=0)  # CST3 gene (approx)
    for d, col, lab, al in ((prot, PC, "CST3 pQTL", 0.5), (siF, FC, "Sarcopenia index (F)", 0.5),
                            (siM, MC, "Sarcopenia index (M)", 0.5)):
        xs = [p / 1e6 for (p, z) in d.values() if lo <= p / 1e6 <= hi]
        zs = [z for (p, z) in d.values() if lo <= p / 1e6 <= hi]
        ax.scatter(xs, zs, s=4, color=col, alpha=al, zorder=2, label=lab)
    ax.axvline(LEAD_PROT, color=PC, lw=0.8, ls=":")
    ax.axvline(LEAD_SI, color=FC, lw=0.8, ls=":")
    ax.axhline(0, color="0.7", lw=0.6, ls="--")
    ax.text(LEAD_PROT + 0.002, ax.get_ylim()[0] * 0.9, "CST3 lead", color=PC, fontsize=7)
    ax.text(LEAD_SI - 0.002, ax.get_ylim()[1] * 0.9, "SI lead", color=FC, fontsize=7, ha="right")
    ax.set_xlabel("chr20 (Mb, GRCh38)")
    ax.set_ylabel("z")
    ax.legend(loc="upper right", frameon=False, fontsize=8, handletextpad=0.3)
    ax.text(-0.22, 1.02, "a", transform=ax.transAxes, fontsize=12, va="bottom")

    # ---- b: per-SNP coupling
    ax = fig.add_subplot(gs[0, 1])
    sh = set(prot) & set(siF)
    zp = np.array([prot[r][1] for r in sh])
    zs = np.array([siF[r][1] for r in sh])
    ax.scatter(zp, zs, s=4, color="0.4", alpha=0.3, zorder=2)
    k, b0 = np.polyfit(zp, zs, 1)
    xr = np.linspace(zp.min(), zp.max(), 50)
    ax.plot(xr, k * xr + b0, color=FC, lw=1.2, ls="--", zorder=3)
    r = np.corrcoef(zp, zs)[0, 1]
    ax.text(0.05, 0.92, f"r = {r:+.3f}", transform=ax.transAxes, fontsize=9, color=FC)
    ax.axhline(0, color="0.7", lw=0.6, ls="--"); ax.axvline(0, color="0.7", lw=0.6, ls="--")
    ax.set_xlabel("CST3 pQTL z"); ax.set_ylabel("Sarcopenia index z (female)")
    ax.text(-0.28, 1.02, "b", transform=ax.transAxes, fontsize=12, va="bottom")

    for ext in ("pdf", "tif"):
        fig.savefig(os.path.join(FIG, f"FigS1.{ext}"), dpi=600, bbox_inches="tight")
    plt.close(fig)
    print("FigS1 done; r =", r, "n =", len(sh))


if __name__ == "__main__":
    main()
