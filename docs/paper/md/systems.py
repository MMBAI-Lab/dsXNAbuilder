"""Use-case systems of the paper, written as a dsXNAbuilder user would type them.

Run from the repository root:  python docs/paper/md/systems.py <outdir>
Builds each duplex, writes the PDB and the advanced Amber/modXNA package
(OPC water, truncated octahedron, 12 A buffer) for the MD use case.
"""
import io, json, os, sys, zipfile
sys.path.insert(0, ".")
from dsxnabuilder.service import Session
from dsxnabuilder.reslib import ResidueLibrary


def alternating(core, overhang, first="m"):
    """2'-OMe/2'-F alternating strand, PS on the two terminal linkages at each end."""
    seq = core + overhang
    toks = [("m" if (i % 2 == 0) == (first == "m") else "f") + c for i, c in enumerate(seq)]
    out = toks[0]
    for i in range(1, len(toks)):
        out += ("*" if i <= 2 or i >= len(toks) - 2 else "") + toks[i]
    return out


def reverse_tokens(strand_5to3):
    """Write a 5'->3' strand 3'->5' (token order reversed, linkages kept between)."""
    import re
    toks = re.findall(r"\*[12RS]?|[a-z+]*[ACGTU](?:\^\w{3})?", strand_5to3)
    return "".join(reversed(toks))


CORE = "GCAAGCUGACCCUGAAGUU"                     # 19-bp siRNA core
ANTI_CORE = "".join({"A": "U", "U": "A", "G": "C", "C": "G"}[c] for c in reversed(CORE))
GAPMER = "+G*+C*+A*dT*dT*dC*dA*dG*dT*dA*dC*dA*dT*+T*+G*+C"   # 3-10-3 LNA/DNA, full PS

SYSTEMS = {
    "dickerson": {"top": "5'-CGCGAATTCGCG-3'", "bottom": "3'-GCGCTTAAGCGC-5'", "model": 4},
    "gapmer": {"top": "5'-r(GCAAUGUACUGAAUGC)-3'", "bottom": "3'-" + reverse_tokens(GAPMER) + "-5'",
               "model": 22},
    "sirna": {"top": "5'-" + "--" + alternating(CORE, "UU", "m") + "-3'",
              "bottom": "3'-" + reverse_tokens(alternating(ANTI_CORE, "UU", "f")) + "--" + "-5'",
              "model": 20},
}

if __name__ == "__main__":
    out = sys.argv[1]
    s = Session(residue_library=ResidueLibrary.load_default())
    for name, req in SYSTEMS.items():
        req = dict(req, amber={"water": "opc", "box": "octahedral", "buffer": 12.0})
        res = s.build(req)
        if not res["ok"]:
            print(name, "FAILED", [i["message"] for i in res["issues"] if i["level"] == "error"])
            continue
        d = os.path.join(out, name)
        os.makedirs(d, exist_ok=True)
        zipfile.ZipFile(io.BytesIO(res["zip"])).extractall(d)
        with open(os.path.join(d, "request.json"), "w") as fh:
            json.dump(req, fh, indent=1)
        sm = res["summary"]
        print(f"{name}: {sm['n_bp']} columns, {sm['n_residues']} nt, {sm['n_modified']} modXNA, "
              f"O3'-P {sm['o3p']}, clashes {sm['clashes']}; issues {[i['code'] for i in res['issues']]}")
        print("   top:   ", req["top"]); print("   bottom:", req["bottom"])
