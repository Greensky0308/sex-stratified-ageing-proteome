#!/usr/bin/env python3
"""Build the phenotype-wise table of sex-differential causal effects at SULT2A1.

Convention: causal effect of the SULT2A1-lowering allele A (rs62129966, per allele).
The DXA ALM effect is converted from per-SD (Wald) to per-allele by multiplying by
the instrument effect on protein (bx = -0.378).
"""
import csv, os

ROOT = os.path.dirname(os.path.abspath(__file__))
BX = -0.377711                                   # instrument: A -> SULT2A1 (per allele)
OUT = os.path.join(ROOT, "tables")

# (phenotype, class, source, bF_perA, seF, bM_perA, seM, z, p)
def main():
    os.makedirs(OUT, exist_ok=True)
    rows = []


    def add(name, cls, src, bF, sF, bM, sM, z, p):
        rows.append([name, cls, src, bF, sF, bM, sM, z, p])


    # --- primary: DXA ALM, convert per-SD -> per-allele-A (multiply by BX) ---
    with open(os.path.join(ROOT, "out", "alm_ukbppp", "alm_sexdiff.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["protein"] == "SULT2A1":
                bF, sF = float(r["bF"]) * BX, float(r["seF"]) * abs(BX)
                bM, sM = float(r["bM"]) * BX, float(r["seM"]) * abs(BX)
                z = float(r["zdiff"])           # sign flips with allele flip
                add("Appendicular lean mass (DXA)", "Muscle mass", "UKB, 244730F/205513M",
                    bF, sF, bM, sM, -z, float(r["p"]))

    # --- replication + specificity + hormones: Neale panel (already per-allele-A) ---
    CLASS = {"Whole-body fat-free mass": "Muscle mass", "Whole-body water mass": "Muscle mass",
             "Arm fat-free mass (L)": "Muscle mass", "Leg fat-free mass (R)": "Muscle mass",
             "Whole-body water mass": "Muscle mass", "Basal metabolic rate": "Muscle mass",
             "Grip strength (L)": "Muscle strength",
             "Grip strength (R)": "Muscle strength", "IGF-1": "Hormone",
             "Testosterone": "Hormone", "SHBG": "Hormone", "Oestradiol": "Hormone",
             "Whole-body fat mass": "Control", "Arm fat mass (L)": "Control",
             "BMI": "Control", "Height": "Control", "Weight": "Control",
             "Heel BMD": "Control", "C-reactive protein": "Control",
             "Cystatin C": "Control", "Creatinine": "Control"}
    SRC = "UKB (Neale r2)"
    with open(os.path.join(ROOT, "out", "sult2a1_neale_phewas_all.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            t = r["trait"]
            if t not in CLASS:
                continue
            add(t, CLASS[t], SRC, float(r["bF"]), "", float(r["bM"]), "",
                float(r["zdiff"]), float(r["p"]))

    order = {"Muscle mass": 0, "Muscle strength": 1, "Hormone": 2, "Control": 3}
    rows.sort(key=lambda x: (order[x[1]], -abs(x[7])))

    # ---------------- Excel ----------------
    import openpyxl
    from openpyxl.styles import Font, Alignment
    wb = openpyxl.Workbook(); ws = wb.active; ws.title = "phenotypes"
    hdr = ["Phenotype", "Class", "Source", "beta_F (A)", "beta_M (A)", "z (F vs M)", "P"]
    ws.append(hdr)
    for c in ws[1]:
        c.font = Font(bold=True); c.alignment = Alignment(horizontal="center")
    for r in rows:
        ws.append([r[0], r[1], r[2], round(r[3], 4), round(r[5], 4),
                   round(r[7], 2), f"{r[8]:.2g}"])
    for col, w in zip("ABCDEFG", [30, 16, 20, 11, 11, 10, 10]):
        ws.column_dimensions[col].width = w
    wb.save(os.path.join(OUT, "phenotype_table.xlsx"))
    print("wrote phenotype_table.xlsx")

    # ---------------- Word three-line table ----------------
    import docx
    from docx.shared import Pt
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    doc = docx.Document()
    doc.add_paragraph("Sex-differential causal effect of SULT2A1 (per SULT2A1-lowering allele A) "
                      "across ageing phenotypes.")
    t = doc.add_table(rows=1, cols=7)
    t.style = "Table Grid"
    for i, h in enumerate(hdr):
        cell = t.rows[0].cells[i]; cell.text = h
        for p in cell.paragraphs:
            for run in p.runs:
                run.font.bold = True
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER


    def three_line(tbl):
        """Apply top/header/bottom borders only (three-line table)."""
        from docx.oxml.ns import qn
        from docx.oxml import OxmlElement
        def edge(cell, name, val):
            tcPr = cell._tc.get_or_add_tcPr()
            borders = tcPr.find(qn("w:tcBorders"))
            if borders is None:
                borders = OxmlElement("w:tcBorders"); tcPr.append(borders)
            e = OxmlElement("w:" + name); e.set(qn("w:val"), val)
            e.set(qn("w:sz"), "8"); borders.append(e)
        for r_i, row in enumerate(tbl.rows):
            for cell in row.cells:
                edge(cell, "top", "single" if r_i == 0 else "nil")
                edge(cell, "bottom", "single" if r_i in (0, len(tbl.rows) - 1) else "nil")
                edge(cell, "left", "nil"); edge(cell, "right", "nil")


    for r in rows:
        t.add_row().cells
        cells = t.rows[-1].cells
        for i, v in enumerate([r[0], r[1], r[2], f"{r[3]:+.4f}", f"{r[5]:+.4f}",
                               f"{r[7]:+.2f}", f"{r[8]:.2g}"]):
            cells[i].text = v
    three_line(t)
    doc.add_paragraph("beta_F, beta_M: causal effect of allele A on the phenotype in women/men. "
                      "z = (beta_F - beta_M)/SE, two-sided. A = SULT2A1-lowering allele (rs62129966). "
                      "Control rows are phenotypes expected to be null.")
    doc.save(os.path.join(OUT, "phenotype_table.docx"))
    print("wrote phenotype_table.docx")


if __name__ == "__main__":
    main()
