#!/usr/bin/env python3
"""Main figures Fig1-Fig4 (each a-d) -> figures/Fig*.pdf + .tif.

Style per project standard: Arial 11 pt (not bold), neutral data points #303030,
emphasis #C00000, category colours #5B9BD5 / #D08A6A, legend below the axis and
horizontally centred, no in-image title (panel letter only), TIFF at 600 dpi.
"""
import csv, math, os, sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(ROOT, "figures")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

plt.rcParams.update({
    "font.family": "Arial", "font.size": 11, "axes.linewidth": 0.8,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "xtick.major.size": 3.5, "ytick.major.size": 3.5,
    "axes.labelsize": 11, "xtick.labelsize": 11, "ytick.labelsize": 11,
    "legend.fontsize": 11, "font.weight": "normal", "axes.titleweight": "normal",
    "axes.spines.top": False, "axes.spines.right": False,
})
BX = -0.377711            # instrument: A -> SULT2A1 (per allele)
NEU = "#303030"                       # neutral data points / lines
EMPH = "#C00000"                      # emphasis / highlight
F_COL, M_COL = "#5B9BD5", "#D08A6A"   # category colours: female / male


_PANELS = []


def panel(ax, letter):
    _PANELS.append((ax, letter))


def _place_panels(fig):
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    rows = []
    for ax, letter in _PANELS:
        if ax not in fig.axes:
            continue
        pos = ax.get_position()
        ylx = None
        yl = ax.yaxis.get_label()
        if yl.get_text():
            ylx = yl.get_window_extent(r).x0
        text_xs = ([ylx] if ylx is not None else []) +                   [t.get_window_extent(r).x0 for t in ax.get_yticklabels() if t.get_text()]
        rows.append([letter, ylx, (min(text_xs) if text_xs else None), pos.x0, pos.y1])
    for row in rows:
        col = [o[2] for o in rows if o[2] is not None and abs(o[3] - row[3]) < 0.2]
        if row[1] is not None:
            anchor = row[1]
        elif col:
            anchor = min(col)
        else:                       # axis-off panel (no ylabel/tick text): sit left of the column
            anchor = row[3] - 0.045
        row.append(anchor)
    for letter, _, _, _, y1, anchor in rows:
        fig.text(anchor / fig.bbox.width, y1 + 0.012, letter, fontsize=12, ha="left", va="bottom")
    _PANELS.clear()


def save(fig, name):
    os.makedirs(FIG, exist_ok=True)
    _place_panels(fig)
    for ext in ("pdf", "tif"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(name, "done")


def leg(ax, ncol=2):
    ax.legend(frameon=False, ncol=ncol, loc="upper center",
              bbox_to_anchor=(0.5, -0.30), handletextpad=0.4, columnspacing=1.2)


def load_row(f, prot="SULT2A1"):
    for i, l in enumerate(open(f)):
        if i == 0:
            continue
        c = l.rstrip("\n").split("\t")
        if c[0] == prot:
            return [float(x) for x in c[1:5]]
    return None


# ---------------------------------------------------------------- panel set 1
def fig1():
    fig = plt.figure(figsize=(8.0, 6.4))
    gs = fig.add_gridspec(2, 2, hspace=0.34, wspace=0.35)

    ax = fig.add_subplot(gs[0, 0]); ax.axis("off"); panel(ax, "a")
    ax.set_xlim(0, 12); ax.set_ylim(0, 6)

    def box(x, y, t, fc):
        ax.text(x, y, t, ha="center", va="center", fontsize=9.5,
                bbox=dict(boxstyle="round,pad=0.45", fc=fc, ec=NEU, lw=0.8))

    box(2.0, 3.9, "1,884\ncis-pQTLs", "#CFE2F3")
    box(6.0, 3.9, "MR\nby sex", "#FFF2CC")
    box(10.0, 3.9, "ageing\noutcomes", "#E4D7F5")
    ax.text(4.0, 3.9, "→", ha="center", va="center", fontsize=14, color=NEU)
    ax.text(8.0, 3.9, "→", ha="center", va="center", fontsize=14, color=NEU)
    ax.text(6.0, 2.6, "outcomes: brain IDPs · muscle mass", ha="center", fontsize=8.5)
    ax.text(6.0, 1.3, "z = (b$_F$ − b$_M$)/SE · FDR · colocalization",
            ha="center", va="center", fontsize=9.5, style="italic")

    ax = fig.add_subplot(gs[0, 1]); panel(ax, "b")
    d = list(csv.DictReader(open(ROOT + "/out/sexdiff_ukbppp2/cis.tsv"), delimiter="\t"))
    z = np.array([float(r["z"]) for r in d])
    pag = next(i for i, r in enumerate(d) if r["protein"] == "PAG1")
    ax.scatter(np.arange(len(z)), z, s=4, color=NEU)
    ax.scatter([pag], [z[pag]], s=45, color=EMPH, zorder=3)
    ax.axhline(0, color="0.7", lw=0.7, ls="--")
    ax.set_ylabel("sex-difference z"); ax.set_xlabel("cis proteins (brain IDPs)")
    ax.annotate("PAG1\ncoloc H4 ≤ 0.06", (pag, z[pag]), (0.42, 0.12),
                textcoords="axes fraction", fontsize=9,
                arrowprops=dict(arrowstyle="->", lw=0.7))

    ax = fig.add_subplot(gs[1, 0]); panel(ax, "c")
    ys = [1, 0]; vals = [4.52, -0.04]
    ax.axvline(0, color=NEU, lw=0.7)
    ax.axvline(1.96, color="0.5", lw=0.6, ls=":")
    ax.axvline(-1.96, color="0.5", lw=0.6, ls=":")
    ax.hlines(ys, 0, vals, color="0.82", lw=1.0, zorder=1)
    ax.scatter(vals, ys, s=45, color=[EMPH, "0.55"], zorder=3)
    ax.text(4.90, 1, "+4.52", va="center", fontsize=10)
    ax.text(0.18, 0, "-0.04", va="center", fontsize=10)
    ax.text(0.15, 1.48, "CST3 = the index denominator", fontsize=8.5, color=NEU, ha="left")
    ax.set_yticks(ys)
    ax.set_yticklabels(["sarcopenia index\n(creatinine / cystatin C)", "true ALM\n(DXA)"],
                       fontsize=8.5)
    ax.set_xlabel("CST3 causal effect (sex-difference z)")
    ax.set_xlim(-0.6, 6.2); ax.set_ylim(-0.6, 1.9)
    ax.tick_params(axis="y", length=0)

    ax = fig.add_subplot(gs[1, 1]); panel(ax, "d")
    rows = [l.rstrip("\n").split("\t") for l in
            open(ROOT + "/out/alm_ukbppp/alm_sexdiff.tsv")][1:]
    zz = np.array([float(r[5]) for r in rows]); prots = [r[0] for r in rows]
    s_i = prots.index("SULT2A1")
    ax.scatter(np.arange(len(zz)), zz, s=4, color=NEU)
    ax.scatter([s_i], [zz[s_i]], s=45, color=EMPH, zorder=3)
    ax.axhline(0, color="0.7", lw=0.7, ls="--")
    ax.set_ylabel("sex-difference z"); ax.set_xlabel("cis proteins (ALM)")
    ax.annotate("SULT2A1\nq = 8.6×10$^{-4}$", (s_i, zz[s_i]), (0.30, 0.12),
                textcoords="axes fraction", fontsize=9,
                arrowprops=dict(arrowstyle="->", lw=0.7))
    save(fig, "figure1")



# ---------------------------------------------------------------- panel set 2
MASS = [
    ("Appendicular lean mass", ROOT + "/out/alm_ukbppp/alm_sexdiff.tsv"),
    ("Whole-body fat-free mass", ROOT + "/out/neale_23101/ffm_sexdiff.tsv"),
    ("Leg fat-free mass", ROOT + "/out/neale_23117/ffm_sexdiff.tsv"),
    ("Arm fat-free mass", ROOT + "/out/neale_23121/ffm_sexdiff.tsv"),
]
CAT = {"Whole-body fat-free mass": "lean", "Whole-body water mass": "lean",
       "Arm fat-free mass (L)": "lean", "Leg fat-free mass (R)": "lean",
       "Basal metabolic rate": "lean", "Whole-body fat mass": "size",
       "Arm fat mass (L)": "size", "Leg fat mass (R)": "size", "BMI": "size",
       "Weight": "size", "Height": "size", "Heel BMD": "size",
       "Grip strength (L)": "strength", "Grip strength (R)": "strength",
       "Testosterone": "hormone", "SHBG": "hormone", "Oestradiol": "hormone",
       "IGF-1": "hormone", "C-reactive protein": "metabolic", "HbA1c": "metabolic",
       "AST": "metabolic", "HDL cholesterol": "metabolic", "LDL cholesterol": "metabolic",
       "Urate": "metabolic", "Glucose": "metabolic", "Total cholesterol": "metabolic",
       "Albumin": "metabolic", "Total bilirubin": "metabolic", "Cystatin C": "renal",
       "Creatinine": "renal", "Diastolic BP": "cv", "Systolic BP": "cv",
       "Reaction time": "cognitive"}
CATCOL = {"lean": EMPH, "size": "0.75", "strength": "#3E7CB1", "hormone": "#5B9BD5",
          "metabolic": "#8E7CC3", "renal": "#70AD47", "cv": "#D08A6A", "cognitive": "0.5"}


# ---------------------------------------------------------------- panel set 2
def fig2():
    fs = 9
    fig = plt.figure(figsize=(8.6, 7.8))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.32], hspace=0.36, wspace=0.42)

    # a: forest F vs M (95% CI); legend top-left, y-axis extended for it
    ax = fig.add_subplot(gs[0, 0]); panel(ax, "a")
    rows = [(l, load_row(f)) for l, f in MASS]
    rows = [d for d in rows if d[1]]
    y = np.arange(len(rows))[::-1]
    for i, (lab, v) in enumerate(rows):
        bF, sF, bM, sM = v[0] * BX, v[1] * abs(BX), v[2] * BX, v[3] * abs(BX)
        ax.errorbar(bF, y[i] + 0.17, xerr=1.96 * sF, fmt="o", color="#C05299", ms=5, capsize=3,
                    lw=1.0, label="Female" if i == 0 else None)
        ax.errorbar(bM, y[i] - 0.17, xerr=1.96 * sM, fmt="s", color="#3E7CB1", ms=5, capsize=3,
                    lw=1.0, label="Male" if i == 0 else None)
        if abs(bF / sF) > 1.96:
            ax.text(bF + 1.96 * sF + 0.004, y[i] + 0.17, "*", ha="left", va="center", fontsize=12, color=EMPH)
        if abs(bM / sM) > 1.96:
            ax.text(bM + 1.96 * sM + 0.004, y[i] - 0.17, "*", ha="left", va="center", fontsize=12, color=EMPH)
    ax.axvline(0, color="0.6", lw=0.7, ls="--")
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=fs)
    ax.set_ylim(-0.7, len(rows) + 0.9)
    ax.set_xlabel("causal effect on mass (per allele A)", fontsize=fs)
    ax.tick_params(axis="x", labelsize=fs)
    ax.legend(frameon=False, fontsize=8, loc="upper right", handletextpad=0.4)

    # b: female vs male effect across lean-mass traits (opposite quadrants)
    ax = fig.add_subplot(gs[0, 1]); panel(ax, "b")
    ph = list(csv.DictReader(open(ROOT + "/out/sult2a1_neale_phewas_all.tsv"), delimiter="\t"))
    lean = [r for r in ph if CAT.get(r["trait"]) == "lean"]
    TAG = {"Whole-body fat-free mass": "Whole-body FFM", "Whole-body water mass": "Water mass",
           "Arm fat-free mass (L)": "Arm FFM", "Leg fat-free mass (R)": "Leg FFM",
           "Basal metabolic rate": "BMR"}
    xF = [float(r["bF"]) for r in lean]; yM = [float(r["bM"]) for r in lean]
    ax.axhline(0, color="0.7", lw=0.7, ls="--"); ax.axvline(0, color="0.7", lw=0.7, ls="--")
    ax.scatter(xF, yM, s=45, color="#C05299", zorder=3)
    for k, (r, x, yv) in enumerate(zip(lean, xF, yM)):
        ax.annotate(TAG.get(r["trait"], r["trait"][:10]), (x, yv), fontsize=8,
                    xytext=(6, 6 if k % 2 == 0 else -11), textcoords="offset points")
    ax.set_xlabel("female β (per allele)", fontsize=fs)
    ax.set_ylabel("male β (per allele)", fontsize=fs)
    ax.tick_params(labelsize=fs)
    ax.margins(0.28)

    # c: specificity PheWAS (lollipop), diagonal small labels
    ax = fig.add_subplot(gs[1, 0]); panel(ax, "c")
    order = ["lean", "strength", "hormone", "size", "metabolic", "renal", "cv", "cognitive"]
    rws = [(r["trait"], float(r["zdiff"]), CAT.get(r["trait"], "size")) for r in ph]
    rws.sort(key=lambda x: (order.index(x[2]) if x[2] in order else 9, x[1]))
    y = np.arange(len(rws))
    ax.axvspan(-1.96, 1.96, color="0.93", zorder=0)
    ax.axvline(1.96, color="0.6", lw=0.6, ls=":")
    ax.hlines(y, 0, [r[1] for r in rws], color="0.82", lw=0.6, zorder=1)
    ax.scatter([r[1] for r in rws], y, c=[CATCOL[r[2]] for r in rws], s=22, zorder=2)
    for r, yy in zip(rws, y):
        if abs(r[1]) > 1.96:
            ax.text(r[1] + (0.12 if r[1] > 0 else -0.12), yy, f"{r[1]:+.1f}",
                    va="center", ha="left" if r[1] > 0 else "right", fontsize=6)
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rws], fontsize=6, rotation=20)
    ax.set_xlabel("sex-difference z", fontsize=fs)
    ax.text(-1.9, len(rws) + 0.4, "|z|<1.96 (n.s.)", fontsize=6.5, color="0.4")
    ax.set_ylim(-1, len(rws) + 1.6)
    ax.tick_params(axis="y", length=0)

    # d: strength (grip) F vs M, annotate significance
    ax = fig.add_subplot(gs[1, 1]); panel(ax, "d")
    g = load_row(ROOT + "/out/grip_mw/alm_sexdiff.tsv")
    if g:
        bF, sF, bM, sM = g[0] * BX, g[1] * abs(BX), g[2] * BX, g[3] * abs(BX)
        ax.errorbar(bF, 0.16, xerr=1.96 * sF, fmt="o", color="#C05299", ms=6, capsize=3, lw=1.0)
        ax.errorbar(bM, -0.16, xerr=1.96 * sM, fmt="s", color="#3E7CB1", ms=6, capsize=3, lw=1.0)
        zM = abs(bM / sM)
        pM = math.erfc(zM / math.sqrt(2))
        star = "***" if pM < 1e-3 else ("**" if pM < 1e-2 else "*")
        ax.text(bM + 1.96 * sM + 0.012, -0.16, star, ha="left", va="center", fontsize=13, color=EMPH)
        ax.text(bM, -0.46, f"β = {bM:+.3f}, P = {pM:.1e}", ha="center", fontsize=8, color=EMPH)
    ax.axvline(0, color="0.6", lw=0.7, ls="--")
    ax.set_yticks([0.16, -0.16]); ax.set_yticklabels(["Female", "Male"], fontsize=fs)
    ax.set_ylim(-0.7, 0.7)
    ax.set_xlabel("causal effect on muscle weakness (per allele A)", fontsize=fs)
    ax.tick_params(axis="x", labelsize=fs)
    save(fig, "figure2")


# ---------------------------------------------------------------- panel set 3
def fig3():
    """Fig3 V3: magenta/blue/teal palette (no red-green, no grey); a x-axis extended left; d legend top-right."""
    from sult2a1_crossplatform import read_somascan, SOMA, COMP
    from alm_locus_dump import dump
    fs = 9
    FC, MC, PC = "#C05299", "#3E7CB1", "#2E8B8B"   # female / male / protein
    LO, HI = 47_650_000, 47_930_000
    ss = read_somascan(SOMA)
    zss = {rs: ss[rs]["beta"] / ss[rs]["se"] for rs in ss
           if LO <= ss[rs]["pos"] <= HI and ss[rs]["se"] > 0}
    posss = {rs: ss[rs]["pos"] for rs in ss if LO <= ss[rs]["pos"] <= HI}
    ea_ss = {rs: ss[rs]["ea"] for rs in ss}

    def orient(raw, rs):
        pos, eo, oo, b, se, z = raw[rs]
        e = ea_ss.get(rs)
        if e is None:
            return None
        if eo == e:
            return b / se
        if oo == e or oo == COMP.get(eo) or eo == COMP.get(e):
            return -b / se
        return None

    fig = plt.figure(figsize=(9.2, 7.4))
    gs = fig.add_gridspec(2, 2, hspace=0.36, wspace=0.42)
    F = dump("GCST90000027", "19", LO, HI); M = dump("GCST90000026", "19", LO, HI)
    xs, zf, zm, zv = [], [], [], []
    for rs in (set(F) & set(M) & set(zss)):
        vf, vm = orient(F, rs), orient(M, rs)
        if vf is None or vm is None:
            continue
        xs.append(posss[rs] / 1e6); zf.append(vf); zm.append(vm); zv.append(zss[rs])
    o = np.argsort(xs)
    xs = np.array(xs)[o]; zf = np.array(zf)[o]; zm = np.array(zm)[o]; zv = np.array(zv)[o]

    # a
    ax = fig.add_subplot(gs[0, 0]); panel(ax, "a")
    ax.axvspan(47.8703, 47.8865, color="#F3E6CF", alpha=0.8, zorder=0)
    ax.axvline(47.871693, color="0.35", lw=0.9, ls=":")
    ax.axhline(0, color="0.7", lw=0.7, ls="--")
    ax.scatter(xs, zf, s=4, color=FC, alpha=0.35, label="ALM female", zorder=2)
    ax.scatter(xs, zm, s=4, color=MC, alpha=0.35, label="ALM male", zorder=2)
    ax.text(47.8695, 6.9, "SULT2A1", ha="right", va="top", fontsize=7.5, color="#8a6a2a")
    ax.text(47.8695, 5.7, "rs62129966", ha="right", va="top", fontsize=7.5, color="0.3")
    ax.set_ylabel("ALM z"); ax.set_xlabel("chr19 (Mb, GRCh38)", fontsize=fs)
    ax.set_xlim(47.65, 47.97); ax.set_ylim(-7.2, 7.2); ax.tick_params(labelsize=fs)
    ax.legend(loc="upper right", frameon=False, fontsize=8, handletextpad=0.3)

    # b
    ax = fig.add_subplot(gs[0, 1]); panel(ax, "b")
    w = [25, 50, 100, 200, 500, 1000]
    H4 = {"F": [0.839, 0.838, 0.835, 0.831, 0.799, 0.000],
          "M": [0.237, 0.235, 0.233, 0.110, 0.048, 0.000]}
    ax.axhline(0.75, color=EMPH, lw=0.9, ls="--", zorder=1)
    ax.text(55, 0.67, "H4 > 0.75: colocalised", fontsize=7, color=EMPH)
    for sex, col, mk, lab in (("F", FC, "o", "Female"), ("M", MC, "s", "Male")):
        ax.fill_between(w, 0, H4[sex], color=col, alpha=0.13, zorder=0)
        ax.plot(w, H4[sex], mk + "-", color=col, ms=6, lw=1.2,
                label=lab, zorder=2)
        ax.annotate(f"{H4[sex][0]:.2f}" + ("*" if sex == "F" else ""),
                    (w[0], H4[sex][0]), (6, 5), textcoords="offset points", fontsize=8, color=col)
    ax.set_xscale("log"); ax.set_xticks(w); ax.set_xticklabels([str(x) for x in w])
    ax.set_xlabel("coloc window (± kb)", fontsize=fs)
    ax.set_ylabel("PP.H4 (SULT2A1 × ALM)")
    ax.text(28, 1.02, "±1,000 kb: PP.H4 = 0, H3 = 1", fontsize=7, color="0.3", ha="left", va="top")
    ax.set_ylim(-0.02, 1.05); ax.tick_params(labelsize=fs)
    ax.legend(loc="upper right", frameon=False, fontsize=8, handletextpad=0.3)

    # c
    ax = fig.add_subplot(gs[1, 0]); panel(ax, "c")
    labs = ["SULT2A1 pQTL\n(SomaScan)", "ALM female", "ALM male"]
    vals = [-23.5, 4.0, -3.2]; cols = [PC, FC, MC]
    yy = np.arange(len(labs))[::-1]
    ax.axvline(0, color="0.7", lw=0.7, ls="--")
    ax.barh(yy, vals, color=cols, height=0.5, zorder=2)
    for yv, v in zip(yy, vals):
        ax.text(v + (0.9 if v > 0 else -0.9), yv, f"{v:+.1f}", va="center",
                ha="left" if v > 0 else "right", fontsize=8)
    ax.set_yticks(yy); ax.set_yticklabels(labs, fontsize=fs)
    ax.set_xlabel("z at rs62129966", fontsize=fs)
    ax.set_xlim(-30, 14); ax.tick_params(axis="y", length=0)

    # d (gene-region window 47.84-47.93, matching the reported r)
    ax = fig.add_subplot(gs[1, 1]); panel(ax, "d")
    msk = (xs >= 47.84) & (xs <= 47.93)
    zvg, zfg, zmg = zv[msk], zf[msk], zm[msk]
    for lab, col, zz in (("ALM female", FC, zfg), ("ALM male", MC, zmg)):
        ax.scatter(zvg, zz, s=7, color=col, alpha=0.4, label=lab, zorder=2)
        k, b = np.polyfit(zvg, zz, 1)
        xr = np.linspace(zvg.min(), zvg.max(), 50)
        ax.plot(xr, k * xr + b, color=col, lw=1.2, ls="--", zorder=3)
        r = np.corrcoef(zvg, zz)[0, 1]
        ax.text(0.04, 0.94 if col == FC else 0.85, f"r = {r:+.2f}",
                transform=ax.transAxes, fontsize=8, color=col, va="top")
    ax.axhline(0, color="0.7", lw=0.7, ls="--"); ax.axvline(0, color="0.7", lw=0.7, ls="--")
    ax.set_xlabel("SULT2A1 pQTL z", fontsize=fs); ax.set_ylabel("ALM z")
    ax.set_ylim(-7.2, 7.2); ax.tick_params(labelsize=fs)
    ax.legend(loc="upper right", frameon=False, fontsize=8, handletextpad=0.3)
    save(fig, "figure3")


# ---------------------------------------------------------------- panel set 1


# ---------------------------------------------------------------- panel set 4
def fig4():
    """Fig4 V1: Fig3-consistent palette/style; mediator triangle; forest 4c; 5-colour lollipop 4d."""
    fs = 9
    FC, MC = "#C05299", "#3E7CB1"
    fig = plt.figure(figsize=(9.2, 7.4))
    gs = fig.add_gridspec(2, 2, hspace=0.38, wspace=0.42)

    # a: IGF-1 colocalization (ribbons like Fig3b)
    ax = fig.add_subplot(gs[0, 0]); panel(ax, "a")
    data = {"F": [], "M": []}
    for r in csv.DictReader(open(ROOT + "/out/igf1_sult2a1_coloc.tsv"), delimiter="\t"):
        if r["H4_shared"]:
            data[r["sex"]].append((float(r["window_kb"]), float(r["H4_shared"])))
    for sex, col, mk, lab in (("F", FC, "o", "Female"), ("M", MC, "s", "Male")):
        d = sorted(data[sex])
        xs = [a for a, _ in d]; ys = [b for _, b in d]
        ax.fill_between(xs, 0, ys, color=col, alpha=0.13, zorder=0)
        ax.plot(xs, ys, mk + "-", color=col, ms=6, lw=1.2, label=lab, zorder=2)
        ax.annotate(f"{ys[0]:.2f}", (xs[0], ys[0]), (6, 5), textcoords="offset points", fontsize=7.5, color=col)
    ax.axhline(0.75, color=EMPH, lw=0.9, ls="--")
    ax.text(55, 0.78, "H4 > 0.75: colocalised", fontsize=7.5, color=EMPH)
    ax.set_xscale("log"); ax.set_xticks([50, 100, 250, 500, 1300])
    ax.set_xticklabels(["50", "100", "250", "500", "1300"])
    ax.set_xlabel("coloc window (± kb)", fontsize=fs); ax.set_ylabel("PP.H4 (SULT2A1 × IGF-1)", fontsize=fs)
    ax.set_ylim(-0.05, float(os.environ.get("FIG4_YMAX", "1.40")))
    ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.2))
    ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(0.1))
    ax.tick_params(axis="y", which="minor", length=2.0)
    ax.tick_params(labelsize=fs)
    ax.legend(loc="upper right", frameon=False, fontsize=7.5, handletextpad=0.3)

    # b: mediation TRIANGLE
    ax = fig.add_subplot(gs[0, 1]); ax.axis("off"); panel(ax, "b")
    def nb(cx, cy, t, fc):
        ax.add_patch(FancyBboxPatch((cx - 0.14, cy - 0.10), 0.28, 0.20, boxstyle="round,pad=0.02",
                                    fc=fc, ec="0.35", lw=0.8))
        ax.text(cx, cy, t, ha="center", va="center", fontsize=10)
    nb(0.17, 0.22, "SULT2A1", "#E4D7F5")
    nb(0.50, 0.70, "IGF-1", "#CFE2F3")
    nb(0.83, 0.22, "Lean mass", "#FFF2CC")
    # arrowheads, stopping at box edges
    ax.add_patch(FancyArrowPatch((0.30, 0.34), (0.38, 0.58), arrowstyle="-|>",
                                 mutation_scale=15, lw=1.1, shrinkA=0, shrinkB=0))
    ax.add_patch(FancyArrowPatch((0.62, 0.58), (0.70, 0.34), arrowstyle="-|>",
                                 mutation_scale=15, lw=1.1, shrinkA=0, shrinkB=0))
    ax.add_patch(FancyArrowPatch((0.35, 0.15), (0.65, 0.15), arrowstyle="-|>",
                                 mutation_scale=13, lw=0.9, ls="--", color="0.45", shrinkA=0, shrinkB=0))
    ax.text(0.285, 0.50, "β=0.034", fontsize=fs, ha="right")
    ax.text(0.715, 0.50, "β=0.17", fontsize=fs, ha="left")
    ax.text(0.50, 0.06, "direct", ha="center", fontsize=fs, color="0.45")
    ax.text(0.5, 0.94, "indirect (mediated) = 41% (female); male sign-discordant",
            ha="center", fontsize=fs)
    ax.set_xlim(0, 1); ax.set_ylim(0.02, 1.02)

    # c: reverse MR forest (points + 95% CI + null line)
    ax = fig.add_subplot(gs[1, 0]); panel(ax, "c")
    labs = ["Female", "Male"]; b = [-0.0219, -0.0120]; se = [0.0180, 0.0196]
    yy = [1, 0]
    ax.axvline(0, color="0.6", lw=0.7, ls="--")
    for lab, bb, ss, yv, col in zip(labs, b, se, yy, [FC, MC]):
        ax.errorbar(bb, yv, xerr=1.96 * ss, fmt="o", color=col, ms=6, capsize=3, lw=1.1)
        ax.text(bb, yv + 0.22, f"{bb:+.3f} ({bb-1.96*ss:+.3f}, {bb+1.96*ss:+.3f})",
                ha="center", fontsize=7.5)
    ax.set_yticks(yy); ax.set_yticklabels(labs, fontsize=fs)
    ax.set_xlabel("reverse MR: effect on SULT2A1 (β, 95% CI)", fontsize=fs)
    ax.set_xlim(-0.09, 0.06); ax.set_ylim(-0.5, 1.5); ax.tick_params(axis="y", length=0)

    # d: GTEx lollipop, 5 colours in Fig3 hue family, smaller value labels
    ax = fig.add_subplot(gs[1, 1]); panel(ax, "d")
    tis = ["Liver", "Adrenal gland", "Stomach", "Testis", "Skeletal muscle"]
    tpm = [365.8, 305.7, 4.74, 0.72, 0.0]
    cols = ["#2E8B8B", "#3E7CB1", "#6C6FB0", "#9C5FA6", "#C05299"]
    y = np.arange(len(tis))[::-1]
    ax.hlines(y, 0, tpm, color="0.87", lw=1.4, zorder=1)
    ax.scatter(tpm, y, s=70, color=cols, zorder=3)
    for yv, v in zip(y, tpm):
        ax.text(v + 12, yv, f"{v:g}", va="center", fontsize=7.5)
    ax.set_yticks(y); ax.set_yticklabels(tis, fontsize=fs)
    ax.set_xlabel("SULT2A1 expression (TPM)", fontsize=fs)
    ax.set_xlim(-25, 440); ax.tick_params(axis="y", length=0)
    save(fig, "figure4")


if __name__ == "__main__":
    fns = {"figure1": fig1, "figure2": fig2, "figure3": fig3, "figure4": fig4}
    want = sys.argv[1:]
    for k, fn in fns.items():
        if not want or k in want:
            fn()
