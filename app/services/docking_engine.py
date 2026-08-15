"""
AutoDock Vina execution.

Uses the official `vina` Python bindings (pip install vina) when
available, and falls back to the Vina CLI binary otherwise — useful on
hosts where only the compiled binary can be installed and the Python
bindings aren't buildable.

IMPORTANT — RMSD semantics (read before trusting these numbers in a
manuscript):

  * The Vina CLI's stdout table reports "rmsd l.b." (symmetry-corrected,
    best-case) and "rmsd u.b." (naive, worst-case) RMSD of each pose
    *relative to the top-ranked pose of that same run* — NOT relative to
    a crystallographic reference. `_run_vina_cli` parses and reports
    these as-is.
  * The Python bindings (`vina.Vina.energies()`) do NOT expose this
    table at all — they return raw energy terms
    [total, inter, intra, torsions, best_intra]. There is no bundled
    symmetry-aware RMSD calculation available from the bindings.
    `_run_vina_python` therefore computes its own naive (non-symmetry-
    corrected, atom-order-based) Cartesian RMSD of each pose against the
    top pose, and reports it as `rmsd_ub` only — `rmsd_lb` is left as
    `None` to avoid implying a calculation that was not actually
    performed. State this distinction explicitly in the methods section
    of any paper reporting these numbers.
  * Neither number is RMSD-to-native-pose. For redocking/self-docking
    validation against a co-crystallized reference ligand (e.g. to
    benchmark against SeamDock as your abstract mentions), use
    `rmsd_to_reference()` below instead.
"""
from __future__ import annotations

import logging
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from app.config import settings
from app.services.pocket_detection import GridBox

logger = logging.getLogger(__name__)


class DockingError(RuntimeError):
    pass


@dataclass
class VinaPose:
    mode: int
    affinity: float                    # kcal/mol, predicted binding free energy
    rmsd_lb: Optional[float] = None    # symmetry-corrected RMSD vs. top pose (CLI path only)
    rmsd_ub: Optional[float] = None    # naive RMSD vs. top pose (CLI path, or our own calc)
    rmsd_note: str = ""


def run_vina(
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    box: GridBox,
    output_pdbqt: Path,
    log_file: Path,
    *,
    exhaustiveness: int = 8,
    num_modes: int = 9,
    energy_range: float = 3.0,
    seed: int | None = None,
) -> list[VinaPose]:
    try:
        return _run_vina_python(
            receptor_pdbqt, ligand_pdbqt, box, output_pdbqt,
            exhaustiveness, num_modes, energy_range, seed,
        )
    except ImportError:
        logger.info("`vina` python package unavailable, falling back to CLI binary")
        return _run_vina_cli(
            receptor_pdbqt, ligand_pdbqt, box, output_pdbqt, log_file,
            exhaustiveness, num_modes, energy_range, seed,
        )


def _run_vina_python(
    receptor_pdbqt, ligand_pdbqt, box, output_pdbqt,
    exhaustiveness, num_modes, energy_range, seed,
) -> list[VinaPose]:
    from vina import Vina  # raises ImportError if not installed

    v = Vina(sf_name="vina", cpu=settings.VINA_CPU, seed=seed if seed is not None else 0)
    v.set_receptor(str(receptor_pdbqt))
    v.set_ligand_from_file(str(ligand_pdbqt))
    v.compute_vina_maps(center=list(box.center), box_size=list(box.size))
    v.dock(exhaustiveness=exhaustiveness, n_poses=num_modes)
    v.write_poses(
        str(output_pdbqt), n_poses=num_modes, energy_range=energy_range, overwrite=True
    )

    # energies(): columns are [total, inter, intra, torsions, best_intra] — NOT rmsd.
    # `total` (column 0) is the predicted binding affinity in kcal/mol.
    energies = v.energies(n_poses=num_modes, energy_range=energy_range)

    naive_rmsds = _naive_rmsd_vs_top_pose(output_pdbqt)

    poses = []
    for i, row in enumerate(energies):
        poses.append(
            VinaPose(
                mode=i + 1,
                affinity=float(row[0]),
                rmsd_lb=None,
                rmsd_ub=naive_rmsds[i] if i < len(naive_rmsds) else None,
                rmsd_note=(
                    "naive atom-order Cartesian RMSD vs. top pose; not symmetry-corrected, "
                    "not RMSD-to-reference"
                ),
            )
        )
    logger.info("Vina (python bindings) produced %d poses", len(poses))
    return poses


def _run_vina_cli(
    receptor_pdbqt, ligand_pdbqt, box, output_pdbqt, log_file,
    exhaustiveness, num_modes, energy_range, seed,
) -> list[VinaPose]:
    output_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        settings.VINA_BINARY,
        "--receptor", str(receptor_pdbqt),
        "--ligand", str(ligand_pdbqt),
        "--center_x", str(box.center[0]), "--center_y", str(box.center[1]), "--center_z", str(box.center[2]),
        "--size_x", str(box.size[0]), "--size_y", str(box.size[1]), "--size_z", str(box.size[2]),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--energy_range", str(energy_range),
        "--cpu", str(settings.VINA_CPU),
        "--out", str(output_pdbqt),
    ]
    if seed is not None:
        cmd += ["--seed", str(seed)]

    try:
        proc = subprocess.run(
            cmd, check=True, capture_output=True, text=True,
            timeout=settings.JOB_TIMEOUT_SECONDS,
        )
    except FileNotFoundError as exc:
        raise DockingError(
            f"Vina binary '{settings.VINA_BINARY}' not found on PATH. "
            "Install AutoDock Vina or `pip install vina`."
        ) from exc
    except subprocess.CalledProcessError as exc:
        raise DockingError(f"Vina failed:\n{exc.stderr}") from exc
    except subprocess.TimeoutExpired as exc:
        raise DockingError(
            f"Vina exceeded the {settings.JOB_TIMEOUT_SECONDS}s job timeout."
        ) from exc

    log_file.parent.mkdir(parents=True, exist_ok=True)
    log_file.write_text(proc.stdout)
    poses = _parse_vina_stdout(proc.stdout)
    logger.info("Vina (CLI) produced %d poses", len(poses))
    return poses


def _parse_vina_stdout(stdout: str) -> list[VinaPose]:
    """Parse the Vina results table out of stdout: mode | affinity | rmsd l.b. | rmsd u.b."""
    poses: list[VinaPose] = []
    in_table = False
    for line in stdout.splitlines():
        if line.strip().startswith("-----+"):
            in_table = True
            continue
        if not in_table:
            continue
        parts = line.split()
        if len(parts) < 4 or not parts[0].isdigit():
            continue
        poses.append(
            VinaPose(
                mode=int(parts[0]),
                affinity=float(parts[1]),
                rmsd_lb=float(parts[2]),
                rmsd_ub=float(parts[3]),
                rmsd_note="Vina CLI symmetry-aware RMSD vs. top pose of this run",
            )
        )
    return poses


def _naive_rmsd_vs_top_pose(poses_pdbqt: Path) -> list[float]:
    """Cartesian RMSD (heavy atoms, atom-order matched — valid here since all
    poses share the same ligand connectivity/atom ordering from the same
    prepared PDBQT) of each pose against the first (best-affinity) pose."""
    coords_per_model = _extract_model_coords(poses_pdbqt)
    if not coords_per_model:
        return []
    ref = coords_per_model[0]
    rmsds = []
    for coords in coords_per_model:
        if len(coords) != len(ref):
            rmsds.append(float("nan"))
            continue
        sq_diffs = [
            (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2
            for a, b in zip(coords, ref)
        ]
        rmsds.append((sum(sq_diffs) / len(sq_diffs)) ** 0.5)
    return rmsds


def _extract_model_coords(pdbqt_file: Path) -> list[list[tuple[float, float, float]]]:
    models: list[list[tuple[float, float, float]]] = []
    current: list[tuple[float, float, float]] = []
    in_model = False
    for line in pdbqt_file.read_text().splitlines():
        if line.startswith("MODEL"):
            current, in_model = [], True
        elif line.startswith("ENDMDL"):
            models.append(current)
            in_model = False
        elif in_model and line.startswith(("ATOM", "HETATM")):
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            current.append((x, y, z))
    return models


def rmsd_to_reference(pose_pdb: Path, reference_ligand_pdb: Path) -> Optional[float]:
    """
    Heavy-atom RMSD of a docked pose against a reference (e.g. co-crystallized)
    ligand structure — the metric actually needed for redocking/self-docking
    validation studies (e.g. benchmarking against SeamDock).

    Uses RDKit's GetBestRMS, which performs an atom-mapping search so results
    are correct even when atom order differs between the two files, and is
    robust to trivial symmetry (e.g. a rotatable terminal -COOH or phenyl
    ring) — unlike the naive coordinate RMSD used elsewhere in this module.
    Both inputs must represent the *same* ligand (same connectivity).
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError as exc:
        raise DockingError("RDKit is required for rmsd_to_reference(). pip install rdkit") from exc

    pose_mol = Chem.MolFromPDBFile(str(pose_pdb), removeHs=True)
    ref_mol = Chem.MolFromPDBFile(str(reference_ligand_pdb), removeHs=True)
    if pose_mol is None or ref_mol is None:
        logger.warning("Could not parse one of the ligand files for RMSD-to-reference.")
        return None

    try:
        return AllChem.GetBestRMS(pose_mol, ref_mol)
    except (RuntimeError, ValueError) as exc:
        logger.warning("GetBestRMS failed (likely a connectivity/atom-mapping mismatch): %s", exc)
        return None
