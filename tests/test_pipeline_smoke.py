"""
Smoke tests covering logic that does NOT require Vina/fpocket/PDBFixer/PyMOL
to be installed — these run anywhere (including CI without the full
scientific stack) and catch structural regressions in parsing/box-math.
Full end-to-end pipeline testing requires the real binaries and is left to
a separate integration-test job once deployed.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.services.pocket_detection import _box_from_coords, _parse_fpocket_info
from app.services.docking_engine import _parse_vina_stdout, _extract_model_coords


def test_box_from_coords_centers_and_pads():
    coords = [(0, 0, 0), (10, 10, 10)]
    box = _box_from_coords(coords, padding=5.0, source="test")
    assert box.center == (5.0, 5.0, 5.0)
    assert box.size == (20.0, 20.0, 20.0)  # 10 span + 2*5 padding


def test_box_from_coords_enforces_minimum_size():
    coords = [(0, 0, 0), (1, 1, 1)]  # tiny pocket
    box = _box_from_coords(coords, padding=1.0, source="test")
    assert all(s >= 15.0 for s in box.size)


def test_parse_vina_stdout_extracts_pose_table(tmp_path):
    stdout = """\
Detected 8 CPUs
Reading input ... done.
mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
   1       -7.2          0          0
   2       -6.8        1.9        3.2
   3       -6.1        2.4        4.0
"""
    poses = _parse_vina_stdout(stdout)
    assert len(poses) == 3
    assert poses[0].affinity == -7.2
    assert poses[0].rmsd_lb == 0
    assert poses[1].rmsd_ub == 3.2


def test_parse_fpocket_info(tmp_path):
    info_text = """\
Pocket 1 :
    Score :  0.5
    Druggability Score :  0.812
    Pocket Score :  22.1

Pocket 2 :
    Score :  0.3
    Druggability Score :  0.401
    Pocket Score :  15.0
"""
    info_file = tmp_path / "protein_info.txt"
    info_file.write_text(info_text)
    pockets = _parse_fpocket_info(info_file)
    assert len(pockets) == 2
    # ranked by druggability score, best first
    assert pockets[0]["pocket_id"] == "1"
    assert pockets[0]["druggability_score"] == 0.812


def test_extract_model_coords_from_pdbqt(tmp_path):
    pdbqt = tmp_path / "poses.pdbqt"
    pdbqt.write_text(
        "MODEL 1\n"
        "ATOM      1  C1  LIG A   1       1.000   2.000   3.000  1.00  0.00     0.000 C\n"
        "ENDMDL\n"
        "MODEL 2\n"
        "ATOM      1  C1  LIG A   1       4.000   5.000   6.000  1.00  0.00     0.000 C\n"
        "ENDMDL\n"
    )
    models = _extract_model_coords(pdbqt)
    assert len(models) == 2
    assert models[0] == [(1.0, 2.0, 3.0)]
    assert models[1] == [(4.0, 5.0, 6.0)]


if __name__ == "__main__":
    import inspect
    module = sys.modules[__name__]
    tests = [f for name, f in inspect.getmembers(module, inspect.isfunction) if name.startswith("test_")]
    passed, failed = 0, 0
    for t in tests:
        try:
            sig = inspect.signature(t)
            kwargs = {"tmp_path": Path("/tmp/docksmart_test")} if "tmp_path" in sig.parameters else {}
            if kwargs:
                kwargs["tmp_path"].mkdir(exist_ok=True)
            t(**kwargs)
            print(f"PASS  {t.__name__}")
            passed += 1
        except AssertionError as e:
            print(f"FAIL  {t.__name__}: {e}")
            failed += 1
    print(f"\n{passed} passed, {failed} failed")
    sys.exit(1 if failed else 0)
