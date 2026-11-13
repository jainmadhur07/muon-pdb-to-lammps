"""CHARMM patches that plain residue (.rtp) lookups miss.
Values copied from the CHARMM36 toppar files:
  DISU  : top_all36_prot.rtf  (GROMACS port block CYS2)
  NGLB  : toppar_all36_carb_glycopeptide.str   (Asn-ND2 / GlcNAc-C1 link)
  14bb  : top_all36_carb.rtf   (Gal C1 -> GlcNAc O4 link)
Detected from geometry (heavy-atom distances), not from names."""
from collections import defaultdict
import math

NGLB_ASN = {"CG": ("CC2O1", 0.51), "OD1": ("OC2D1", -0.51), "ND2": ("NC2D1", -0.47),
            "HD22": ("HCP1", 0.31), "CB": ("CC321", -0.18),
            "HB1": ("HCA2", 0.09), "HB2": ("HCA2", 0.09)}
NGLB_SUGAR = {"C1": ("CC3162", 0.27), "H1": ("HCA1", 0.09), "C5": ("CC3163", 0.11),
              "H5": ("HCA1", 0.09), "O5": ("OC3C61", -0.40)}
P14BB_GLCNAC = {"C4": ("CC3161", 0.09), "O4": ("OC301", -0.36)}
P14BB_GAL = {"C1": ("CC3162", 0.29)}


def _d(a, b):
    return math.sqrt((a["x"]-b["x"])**2 + (a["y"]-b["y"])**2 + (a["z"]-b["z"])**2)


def apply_patches(atoms, residues, resdb, log):
    """Mutates atoms (fftype, charge). Returns counts."""
    by_res = {k: {a["name"]: a for a in v} for k, v in residues.items()}
    counts = defaultdict(int)

    def setvals(res_key, table):
        for name, (t, q) in table.items():
            a = by_res[res_key].get(name)
            if a is not None:
                a["fftype"], a["charge"] = t, q

    # --- disulfides: CYS SG-SG < 2.5 A between different residues
    sg = [(k, r["SG"]) for k, r in by_res.items()
          if "SG" in r and next(iter(r.values()))["resname"] in ("CYS", "CYS2")]
    cys2 = resdb.rtp["aminoacids"].get("CYS2")
    used = set()
    for i in range(len(sg) if cys2 else 0):
        for j in range(i+1, len(sg)):
            if sg[i][0][0] == sg[j][0][0] and sg[i][0] == sg[j][0]:
                continue
            if _d(sg[i][1], sg[j][1]) < 2.5:
                for k in (sg[i][0], sg[j][0]):
                    if k in used:
                        continue
                    used.add(k)
                    for name, a in by_res[k].items():
                        if name in cys2:
                            a["fftype"], a["charge"] = cys2[name]
                    counts["DISU_residues"] += 1

    # --- glycans
    sugar = {k: r for k, r in by_res.items()
             if next(iter(r.values()))["resname"] in ("BGLN", "BGLCNA", "BGAL")}
    glcnac = {k: r for k, r in sugar.items()
              if next(iter(r.values()))["resname"] in ("BGLN", "BGLCNA")}
    gal = {k: r for k, r in sugar.items() if next(iter(r.values()))["resname"] == "BGAL"}
    asn = {k: r for k, r in by_res.items()
           if next(iter(r.values()))["resname"] == "ASN" and "ND2" in r}
    for gk, gr in glcnac.items():
        c1 = gr.get("C1")
        for ak, ar in asn.items():
            if _d(ar["ND2"], c1) < 1.8:
                setvals(ak, NGLB_ASN)
                setvals(gk, NGLB_SUGAR)
                counts["NGLB"] += 1
                break
        o4 = gr.get("O4")
        for lk, lr in gal.items():
            if _d(o4, lr["C1"]) < 1.8:
                setvals(gk, P14BB_GLCNAC)
                setvals(lk, P14BB_GAL)
                counts["14BB"] += 1
                break
    return dict(counts)
