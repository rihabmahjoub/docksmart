"""Pydantic schemas used by the API layer."""
from __future__ import annotations

from enum import Enum
from typing import Optional

from pydantic import BaseModel, Field


class JobStatus(str, Enum):
    QUEUED = "queued"
    RUNNING = "running"
    SUCCEEDED = "succeeded"
    FAILED = "failed"


class ReceptorInputType(str, Enum):
    PDB_ID = "pdb_id"
    ALPHAFOLD_UNIPROT = "alphafold_uniprot"
    UPLOAD = "upload"


class LigandInputType(str, Enum):
    SMILES = "smiles"
    CHEMBL_ID = "chembl_id"
    UPLOAD = "upload"
    FROM_INTELLIGENTLIGPREP = "from_intelligentligprep"  # prepared SDF/PDBQT handoff


class PocketSelectionMode(str, Enum):
    AUTO_FPOCKET = "auto_fpocket"          # rank pockets with fpocket, use top hit
    COCRYSTAL_LIGAND = "cocrystal_ligand"  # box built around existing HETATM ligand
    MANUAL = "manual"                      # user-supplied center/size


class DockingParams(BaseModel):
    exhaustiveness: int = Field(8, ge=1, le=32)
    num_modes: int = Field(9, ge=1, le=20)
    energy_range: float = Field(3.0, ge=1.0, le=10.0)
    seed: Optional[int] = None

    # Grid box — populated automatically unless pocket_mode == MANUAL
    pocket_mode: PocketSelectionMode = PocketSelectionMode.AUTO_FPOCKET
    center_x: Optional[float] = None
    center_y: Optional[float] = None
    center_z: Optional[float] = None
    size_x: float = 20.0
    size_y: float = 20.0
    size_z: float = 20.0

    remove_waters: bool = True
    add_hydrogens: bool = True
    ph: float = 7.4
    # Default False: for rigid-receptor docking, rebuilding unresolved loops
    # can introduce modeled geometry that was never observed experimentally
    # and is usually unnecessary unless the missing region is in/near the
    # binding site. Only missing atoms *within* resolved residues are always
    # fixed regardless of this flag. See receptor_prep.fix_receptor().
    rebuild_missing_loops: bool = False


class JobCreateResponse(BaseModel):
    job_id: str
    status: JobStatus


class JobStatusResponse(BaseModel):
    job_id: str
    status: JobStatus
    stage: str
    error_message: Optional[str] = None


class PoseResult(BaseModel):
    mode: int
    affinity_kcal_mol: float
    # See app/services/docking_engine.py module docstring: rmsd_lb is only
    # populated on the Vina-CLI execution path (symmetry-corrected RMSD vs.
    # the top pose of this run). rmsd_ub is always vs. the top pose of this
    # run, never vs. a crystallographic reference.
    rmsd_lb: Optional[float] = None
    rmsd_ub: Optional[float] = None
    rmsd_note: Optional[str] = None
    interaction_summary: Optional[dict] = None
    pose_pdb: Optional[str] = None
    pose_pdbqt: Optional[str] = None


class JobResult(BaseModel):
    poses: list[PoseResult]
    pocket_info: Optional[dict] = None
    best_pose_pdb: Optional[str] = None
    receptor_pdb: Optional[str] = None
    receptor_pdbqt: Optional[str] = None
    render_png: Optional[str] = None
    # Residues Meeko dropped from the receptor during PDBQT conversion due to
    # ambiguous inter-residue bond geometry (see receptor_prep.receptor_to_pdbqt).
    # Empty list in the normal case; surfaced in the UI as a caution banner
    # only when non-empty.
    receptor_dropped_residues: list[str] = []
