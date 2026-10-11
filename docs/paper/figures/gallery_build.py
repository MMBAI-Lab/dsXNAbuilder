"""Build the duplexes rendered in the figures (run from the repository root)."""
import json, os, sys
sys.path.insert(0, ".")
from dsxnabuilder.service import Session
from dsxnabuilder.reslib import ResidueLibrary

OUT = "docs/paper/figures/pdb"
GALLERY = [  # (file, label, top 5'->3', bottom 3'->5', model)
    ("b_dna", "B-DNA (#4)", "CGCGAATTCGCG", "GCGCTTAAGCGC", 4),
    ("a_dna", "A-DNA (#1)", "CGCGAATTCGCG", "GCGCTTAAGCGC", 1),
    ("c_dna", "C-DNA (#7)", "CGCGAATTCGCG", "GCGCTTAAGCGC", 7),
    ("d_dna", "D-DNA (#14)", "ATATATATATAT", "TATATATATATA", 14),
    ("z_dna", "Z-DNA (#15)", "CGCGCGCGCGCG", "GCGCGCGCGCGC", 15),
    ("s_dna", "S-DNA (#49)", "CGCGCGCGCGCG", "GCGCGCGCGCGC", 49),
    ("a_rna", "A-RNA (#20)", "GGCGCUUGCGUC", "CCGCGAACGCAG", 20),
    ("hybrid", "DNA:RNA hybrid (#22)", "r(GACUAGCUGACU)", "d(CTGATCGACTGA)", 22),
]
CASES = [
    ("gapmer", "LNA/PS gapmer : RNA (#22)", "r(GCAAUGUACUGAAUGC)",
     "+C*+G*+T*dT*dA*dC*dA*dT*dG*dA*dC*dT*dT*+A*+C*+G", 22),
    ("sirna", "siRNA 2'-OMe/2'-F/PS (#20)", "--mG*fC*mAfAmGfCmUfGmAfCmCfCmUfGmAfAmGfUmU*fU*mU",
     "fU*mU*fCmGfUmUfCmGfAmCfUmGfGmGfAmCfUmUfC*mA*fA--", 20),
    ("abasic", "abasic site (#4)", "CGCGA-TTCGCG", "GCGCTTAAGCGC", 4),
]
if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    s = Session(residue_library=ResidueLibrary.load_default())
    meta = {}
    for name, label, top, bottom, model in GALLERY + CASES:
        out = s.build({"top": top, "bottom": bottom, "model": model})
        assert out["ok"], (name, out["issues"])
        open(f"{OUT}/{name}.pdb", "w").write(out["pdb"])
        meta[name] = {"label": label, "top": top, "bottom": bottom, "model": model,
                      "summary": out["summary"], "issues": [i["code"] for i in out["issues"]]}
        print(f"{name:8s} {out['summary']['n_atoms']:5d} atoms  O3'-P {out['summary']['o3p']}  clashes {out['summary']['clashes']}")
    json.dump(meta, open(f"{OUT}/meta.json", "w"), indent=1)
