"""
Ligand preparation.

DockSmart intentionally does NOT reimplement ChEMBL retrieval, ADMET
filtering, or descriptor computation — that is IntelligentLigPrep's job.
Three entry points are supported here:

  * SMILES string          -> embed 3D coords (ETKDGv3) + MMFF optimize
  * uploaded file (SDF/MOL2/PDB) -> read as-is
  * a file already produced by IntelligentLigPrep (SDF/PDBQT handoff)
    -> passed straight through, skipping re-preparation entirely
"""
from __future__ import annotations

import logging
from pathlib import Path

from app.config import settings

logger = logging.getLogger(__name__)


class LigandPrepError(RuntimeError):
    pass


def smiles_to_3d_sdf(smiles: str, output_sdf: Path) -> Path:
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError as exc:
        raise LigandPrepError("RDKit is not installed. Install with: pip install rdkit") from exc

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise LigandPrepError(f"Could not parse SMILES: '{smiles}'")

    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise LigandPrepError("3D embedding failed for the given SMILES.")
    AllChem.MMFFOptimizeMolecule(mol)

    output_sdf.parent.mkdir(parents=True, exist_ok=True)
    with Chem.SDWriter(str(output_sdf)) as writer:
        writer.write(mol)

    logger.info("Ligand embedded in 3D: %s", output_sdf)
    return output_sdf


def ligand_to_pdbqt(input_file: Path, output_pdbqt: Path) -> Path:
    """Convert a ligand (SDF/MOL2/PDB) to AutoDock PDBQT via Meeko."""
    import subprocess

    output_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["mk_prepare_ligand.py", "-i", str(input_file), "-o", str(output_pdbqt)]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=settings.MEEKO_TIMEOUT_SECONDS)
    except subprocess.CalledProcessError as exc:
        raise LigandPrepError(f"mk_prepare_ligand.py failed:\n{exc.stderr}") from exc
    except subprocess.TimeoutExpired as exc:
        raise LigandPrepError(
            f"mk_prepare_ligand.py exceeded the {settings.MEEKO_TIMEOUT_SECONDS}s timeout "
            "(usually a slow CPU allocation, not a stuck process) — try raising "
            "DOCKSMART_MEEKO_TIMEOUT."
        ) from exc
    except FileNotFoundError as exc:
        raise LigandPrepError(
            "mk_prepare_ligand.py not found on PATH (part of the 'meeko' package)."
        ) from exc

    if not output_pdbqt.exists():
        raise LigandPrepError("Ligand PDBQT was not produced — check Meeko output.")

    logger.info("Ligand converted to PDBQT: %s", output_pdbqt)
    return output_pdbqt
