"""
Unit tests for the parsing/geometry helpers - the parts that were silently
wrong in the original scripts. Run with:  pytest -q

These tests need no network, no Open Babel and no smina.
"""

import os
import subprocess
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

import step2_build_ligands as s2  # noqa: E402
import step3_dock as s3  # noqa: E402


# --------------------------------------------------------------------------
# step 2 - file naming
# --------------------------------------------------------------------------

def test_safe_filename_strips_unsafe_characters():
    assert s2.safe_filename("Sterculia Gum (karaya gum)") == "Sterculia_Gum_karaya_gum"
    assert s2.safe_filename("cis-9,10-epoxy/stearate") == "cis-9_10-epoxy_stearate"
    assert s2.safe_filename("   ") == "unnamed"


def test_name_collisions_do_not_overwrite():
    """The original silently overwrote one compound with the other."""
    names = ["Vitamin C", "Vitamin/C", "Vitamin-C"]
    out = s2.unique_names(names)
    assert len(set(out)) == 3, out
    assert out[0] == "Vitamin_C" and out[1] == "Vitamin_C_2"


def test_pdbqt_validation_rejects_empty_and_stub_files(tmp_path):
    empty = tmp_path / "empty.pdbqt"
    empty.write_text("")
    assert not s2.pdbqt_is_valid(str(empty))

    stub = tmp_path / "stub.pdbqt"
    stub.write_text("REMARK  nothing here\nATOM      1  C   LIG\n")
    assert not s2.pdbqt_is_valid(str(stub))  # only 1 atom

    good = tmp_path / "good.pdbqt"
    good.write_text("".join(f"ATOM  {i:>5}  C   LIG A   1\n" for i in range(5)))
    assert s2.pdbqt_is_valid(str(good))


# --------------------------------------------------------------------------
# step 3 - score parsing (the most consequential original bug)
# --------------------------------------------------------------------------

SMINA_OUTPUT = """\
   _______  _______ _________ _        _______
  (  ____ \\(       )\\__   __/( (    /|(  ___  )

Using random seed: 42

mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
1       -7.5       0.000      0.000
Refine time 1.234
Loop time 12.5
"""

MULTI_MODE_OUTPUT = """\
mode |   affinity | dist from best mode
-----+------------+----------+----------
1       -9.2       0.000      0.000
2       -8.8       1.234      2.345
10      -6.1       3.210      5.432
1 molecule converted
"""


def test_parse_affinity_reads_mode_one():
    assert s3.parse_affinity(SMINA_OUTPUT) == -7.5


def test_parse_affinity_is_not_fooled_by_other_lines_starting_with_one():
    """
    The original matched any line beginning with '1' and, because it never
    broke out of the loop, the LAST match won - here mode 10, or the
    '1 molecule converted' line.
    """
    assert s3.parse_affinity(MULTI_MODE_OUTPUT) == -9.2


def test_parse_affinity_returns_none_instead_of_a_100_sentinel():
    assert s3.parse_affinity("smina: cannot open receptor file\n") is None
    assert s3.parse_affinity("") is None


# --------------------------------------------------------------------------
# step 3 - receptor cleaning and pocket geometry
# --------------------------------------------------------------------------

RAW_PDB = """\
MODEL        1
ATOM      1  N   ALA A   1      11.104  13.207  10.000  1.00 20.00           N
ATOM      2  CA AALA A   1      12.000  14.000  10.500  1.00 20.00           C
ATOM      3  CA BALA A   1      12.500  14.500  10.700  0.50 20.00           C
ATOM      4  H   ALA A   1      12.100  14.100  11.500  1.00 20.00           H
HETATM    5  O   HOH A 101      20.000  20.000  20.000  1.00 30.00           O
HETATM    6  C1  GOL A 102      21.000  21.000  21.000  1.00 30.00           C
HETATM    7  FE  HEM A 103      15.000  15.000  15.000  1.00 25.00          FE
ENDMDL
MODEL        2
ATOM      8  N   ALA A   1      99.000  99.000  99.000  1.00 20.00           N
ENDMDL
END
"""


def test_clean_pdb_removes_water_altlocs_hydrogens_and_model_two(tmp_path):
    src = tmp_path / "in.pdb"
    src.write_text(RAW_PDB)
    out = tmp_path / "clean.pdb"
    n = s3.clean_pdb(str(src), str(out))

    text = out.read_text()
    assert n == 2                      # N + CA(altloc A) only
    assert "HOH" not in text           # crystallographic water dropped
    assert "GOL" not in text           # cryoprotectant dropped
    assert "HEM" not in text           # HETATM dropped unless --keep-het
    assert "BALA" not in text          # alternate conformer B dropped
    assert "99.000" not in text        # second model dropped
    assert text.count("\n") == 3       # 2 atoms + END


def test_clean_pdb_can_keep_cofactors(tmp_path):
    src = tmp_path / "in.pdb"
    src.write_text(RAW_PDB)
    out = tmp_path / "clean_het.pdb"
    s3.clean_pdb(str(src), str(out), keep_het=True)
    text = out.read_text()
    assert "HEM" in text               # cofactor kept
    assert "HOH" not in text           # water still dropped
    assert "GOL" not in text


POCKET_PQR = """\
ATOM      1    C STP     1       0.000   0.000   0.000  0.00  3.00
ATOM      2    C STP     1      10.000   0.000   0.000  0.00  3.00
ATOM      3    C STP     1       0.000  10.000   0.000  0.00  3.00
ATOM      4    C STP     1       0.000   0.000  10.000  0.00  3.00
"""


def _make_pocket(tmp_path, filename, content):
    prot = tmp_path / "pdb1abc_Test.pdb"
    prot.write_text("ATOM      1  N   ALA A   1       1.000   1.000   1.000\n")
    pockets = tmp_path / "pdb1abc_Test_out" / "pockets"
    pockets.mkdir(parents=True)
    (pockets / filename).write_text(content)
    return str(prot)


def test_pocket_geometry_prefers_alpha_spheres(tmp_path):
    prot = _make_pocket(tmp_path, "pocket1_vert.pqr", POCKET_PQR)
    geom = s3.pocket_geometry(prot, rank=1, padding=4.0, box_min=18.0, box_max=30.0)
    assert geom["source"] == "alpha_spheres"
    assert geom["center"] == pytest.approx((2.5, 2.5, 2.5))
    # extent 10 A + 2*4 padding = 18, which equals box_min
    assert geom["size"] == pytest.approx((18.0, 18.0, 18.0))


def test_box_size_is_clamped(tmp_path):
    wide = POCKET_PQR + "ATOM      5    C STP     1      99.000   0.000   0.000  0.00  3.00\n"
    prot = _make_pocket(tmp_path, "pocket1_vert.pqr", wide)
    geom = s3.pocket_geometry(prot, rank=1, padding=4.0, box_min=18.0, box_max=30.0)
    assert geom["size"][0] == 30.0     # would be 107 without clamping


def test_pocket_geometry_missing_returns_none(tmp_path):
    prot = tmp_path / "nopocket.pdb"
    prot.write_text("ATOM      1  N   ALA A   1       1.000   1.000   1.000\n")
    assert s3.pocket_geometry(str(prot), 1, 4.0, 18.0, 30.0) is None


# --------------------------------------------------------------------------
# step 4 - end-to-end on synthetic data
# --------------------------------------------------------------------------

SYNTHETIC = """Chemical,Protein,Score(kcal/mol),Status
SpecificBinder,P1,-11.0,ok
SpecificBinder,P2,-6.0,ok
SpecificBinder,P3,-6.2,ok
SpecificBinder,P4,-5.9,ok
SpecificBinder,P5,-6.1,ok
SpecificBinder,P6,-6.3,ok
StickyPromiscuous,P1,-9.0,ok
StickyPromiscuous,P2,-9.1,ok
StickyPromiscuous,P3,-8.9,ok
StickyPromiscuous,P4,-9.0,ok
StickyPromiscuous,P5,-9.2,ok
StickyPromiscuous,P6,-8.8,ok
WeakBinder,P1,-4.0,ok
WeakBinder,P2,-4.2,ok
WeakBinder,P3,-3.9,ok
WeakBinder,P4,-4.1,ok
WeakBinder,P5,-4.0,ok
WeakBinder,P6,-4.3,ok
FailedRun,P1,,no_score: smina crashed
ClashingPose,P1,25.0,ok
"""


def test_step4_keeps_specific_drops_promiscuous(tmp_path):
    inp = tmp_path / "docking_results.csv"
    inp.write_text(SYNTHETIC)
    out = tmp_path / "druggable.csv"
    matrix = tmp_path / "matrix.csv"

    script = os.path.join(os.path.dirname(__file__), "..", "scripts",
                          "step4_analyze.py")
    proc = subprocess.run(
        [sys.executable, script, "-i", str(inp), "-o", str(out),
         "--matrix", str(matrix), "--ligand-index", "", "--min-proteins", "5"],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr

    import pandas as pd
    hits = pd.read_csv(out)
    assert set(hits["Chemical"]) == {"SpecificBinder"}
    assert hits.iloc[0]["Protein"] == "P1"
    # the flat, always-sticky profile must be rejected by the std filter
    assert "StickyPromiscuous" not in set(hits["Chemical"])
    # the failed run and the positive (clashing) score must never be scored
    assert "FailedRun" not in proc.stdout.split("Filtering audit")[0]
    assert matrix.exists()
