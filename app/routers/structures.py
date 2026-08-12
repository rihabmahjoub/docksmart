"""Endpoints for previewing/validating a receptor before launching a job."""
from __future__ import annotations

import tempfile
from pathlib import Path

from fastapi import APIRouter, HTTPException

from app.services.structure_fetch import StructureFetchError, fetch_alphafold_model, fetch_pdb_by_id

router = APIRouter(prefix="/api/structures", tags=["structures"])


@router.get("/pdb/{pdb_id}/preview")
def preview_pdb(pdb_id: str):
    """Fetch a PDB entry and return light metadata (used by the frontend to
    show 'found: 2 chains, resolution X Å, contains ligand Y' before the
    user commits to launching a full docking job)."""
    with tempfile.TemporaryDirectory() as tmp:
        try:
            pdb_path = fetch_pdb_by_id(pdb_id, Path(tmp))
        except StructureFetchError as exc:
            raise HTTPException(status_code=404, detail=str(exc)) from exc
        return _summarize_pdb(pdb_path, pdb_id)


@router.get("/alphafold/{uniprot_id}/preview")
def preview_alphafold(uniprot_id: str):
    with tempfile.TemporaryDirectory() as tmp:
        try:
            pdb_path = fetch_alphafold_model(uniprot_id, Path(tmp))
        except StructureFetchError as exc:
            raise HTTPException(status_code=404, detail=str(exc)) from exc
        summary = _summarize_pdb(pdb_path, uniprot_id)
        summary["source"] = "alphafold"
        summary["note"] = (
            "AlphaFold models have no bound ligand and no experimental resolution — "
            "per-residue pLDDT confidence (B-factor column) should be checked near "
            "the intended binding site before docking."
        )
        return summary


def _summarize_pdb(pdb_path: Path, identifier: str) -> dict:
    chains: set[str] = set()
    hetero_resnames: set[str] = set()
    n_atoms = 0
    for line in pdb_path.read_text(errors="ignore").splitlines():
        if line.startswith("ATOM"):
            n_atoms += 1
            chains.add(line[21])
        elif line.startswith("HETATM"):
            resname = line[17:20].strip()
            if resname not in ("HOH", "WAT"):
                hetero_resnames.add(resname)
            chains.add(line[21])
    return {
        "identifier": identifier,
        "n_protein_atoms": n_atoms,
        "chains": sorted(chains),
        "cocrystallized_ligands": sorted(hetero_resnames),
    }
