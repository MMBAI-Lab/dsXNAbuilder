"""Two-strand sequence language."""

import pytest

from dsxnabuilder.seqlang import SequenceParseError, parse_duplex


def P(top, bottom, catalog, aliases, chem=("DNA", "DNA")):
    return parse_duplex(top, bottom, catalog, aliases, chem)


def test_plain_duplex_with_end_labels(catalog, aliases):
    d = P("5'-CGCGAATTCGCG-3'", "3'-GCGCTTAAGCGC-5'", catalog, aliases)
    assert d.n_columns == 12
    assert [r.letter for r in d.top.residues] == list("CGCGAATTCGCG")
    # bottom residues come out 5'->3', i.e. reversed from what is typed
    assert [r.letter for r in d.bottom.residues] == list("CGCGAATTCGCG")
    assert d.top.residues[0].cap == "5cap" and d.top.residues[-1].cap == "3cap"
    assert d.top.residues[0].backbone == "5PO"
    assert all(r.is_canonical for r in d.top.residues + d.bottom.residues)


def test_bottom_linkage_belongs_to_residue_on_its_left(catalog, aliases):
    # typed 3'->5': A*C G  ->  5'->3': G C*A? no: G C A with the '*' between C and A,
    # which is on the 3' side of C ... the linkage typed between A and C belongs
    # to A (left), the 3'-side residue in 5'->3' order.
    d = P("TGC", "A*CG", catalog, aliases)
    res = d.bottom.residues                     # 5'->3': G, C, A
    assert [r.letter for r in res] == ["G", "C", "A"]
    assert res[2].backbone == "PS1" and res[1].backbone == "DPO"


def test_top_linkage_belongs_to_following_residue(catalog, aliases):
    d = P("mA*mC*mG", "UGC", catalog, aliases, ("RNA", "RNA"))
    assert [r.backbone for r in d.top.residues] == ["5PO", "PS1", "PS1"]
    assert all(r.sugar == "OME" for r in d.top.residues)
    assert any(i.code == "PS_STEREO_UNSPECIFIED" for i in d.issues)


def test_explicit_stereo_and_groups(catalog, aliases):
    d = P("d(AC)*2r(G)", "TGC", catalog, aliases)
    top = d.top.residues
    assert [r.chem for r in top] == ["DNA", "DNA", "RNA"]
    assert top[2].backbone == "PS2" and top[2].sugar == "RC3"


def test_base_modifier_and_explicit_block(catalog, aliases):
    d = P("AC^5mCG", "TGC", catalog, aliases)
    assert d.top.residues[1].base == "M5C" and d.top.residues[1].letter == "C"
    assert not d.top.residues[1].is_canonical
    d = P("A{PS1/MOE/RUU}G", "TAC", catalog, aliases)
    r = d.top.residues[1]
    assert r.triplet() == ("PS1", "MOE", "RUU") and r.letter == "U"


def test_overhangs(catalog, aliases):
    d = P("--ACGT", "TTTGCA", catalog, aliases)
    assert [c is None for c in d.top.columns] == [True, True, False, False, False, False]
    assert d.top.residues[0].column == 2


@pytest.mark.parametrize("top,bottom,code", [
    ("AC-GT", "TGACA", "INTERNAL_GAP"),
    ("ACGT", "TGC", "LENGTH_MISMATCH"),
    ("-ACG", "-TGC", "EMPTY_COLUMN"),
    ("3'-ACGT", "TGCA", "WRONG_END_LABEL"),
    ("ACXT", "TGCA", "BAD_CHAR"),
    ("A(CGT", "TGCA", "UNBALANCED"),
    ("A", "T", "TOO_SHORT"),
])
def test_errors(catalog, aliases, top, bottom, code):
    with pytest.raises(SequenceParseError) as e:
        P(top, bottom, catalog, aliases)
    assert e.value.issue.code == code


def test_deoxy_u_and_ribo_t_are_not_amber_canonical(catalog, aliases):
    d = P("dU dA", "AT", catalog, aliases)
    assert not d.top.residues[0].is_canonical
    d = P("r(TA)", "r(AU)", catalog, aliases)
    assert not d.top.residues[0].is_canonical
