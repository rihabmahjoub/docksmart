"""Job lifecycle endpoints: create -> poll status -> fetch results/files."""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
import uuid
from pathlib import Path
from typing import Optional

from fastapi import APIRouter, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse

from app.config import settings
from app.database import create_job, get_job, list_recent_jobs
from app.models import DockingParams, JobCreateResponse, JobStatusResponse, PocketSelectionMode
from app.services.structure_fetch import StructureFetchError, fetch_alphafold_model, fetch_pdb_by_id

router = APIRouter(prefix="/api/jobs", tags=["jobs"])

MAX_UPLOAD_BYTES = settings.MAX_UPLOAD_MB * 1024 * 1024


@router.post("", response_model=JobCreateResponse)
async def create_docking_job(
    receptor_pdb_id: Optional[str] = Form(None),
    receptor_uniprot_id: Optional[str] = Form(None),
    receptor_file: Optional[UploadFile] = File(None),
    ligand_smiles: Optional[str] = Form(None),
    ligand_file: Optional[UploadFile] = File(None),
    exhaustiveness: int = Form(settings.DEFAULT_EXHAUSTIVENESS),
    num_modes: int = Form(settings.DEFAULT_NUM_MODES),
    pocket_mode: PocketSelectionMode = Form(PocketSelectionMode.AUTO_FPOCKET),
    center_x: Optional[float] = Form(None),
    center_y: Optional[float] = Form(None),
    center_z: Optional[float] = Form(None),
    size_x: float = Form(20.0),
    size_y: float = Form(20.0),
    size_z: float = Form(20.0),
    remove_waters: bool = Form(True),
    rebuild_missing_loops: bool = Form(False),
):
    receptor_inputs = [receptor_pdb_id, receptor_uniprot_id, receptor_file]
    if sum(x is not None for x in receptor_inputs) != 1:
        raise HTTPException(
            status_code=400,
            detail="Provide exactly one receptor source: PDB ID, UniProt ID, or a file upload.",
        )
    if (ligand_smiles is None) == (ligand_file is None):
        raise HTTPException(
            status_code=400, detail="Provide exactly one ligand source: SMILES or a file upload."
        )
    if pocket_mode == PocketSelectionMode.MANUAL and None in (center_x, center_y, center_z):
        raise HTTPException(
            status_code=400, detail="Manual pocket mode requires center_x, center_y, center_z."
        )

    job_id = str(uuid.uuid4())
    job_dir = settings.JOBS_DIR / job_id
    job_dir.mkdir(parents=True, exist_ok=True)

    receptor_pdb, receptor_source = await _resolve_receptor(
        job_dir, receptor_pdb_id, receptor_uniprot_id, receptor_file
    )
    ligand_path, ligand_is_smiles, ligand_source = await _resolve_ligand(
        job_dir, ligand_smiles, ligand_file
    )

    params = DockingParams(
        exhaustiveness=exhaustiveness,
        num_modes=num_modes,
        pocket_mode=pocket_mode,
        center_x=center_x, center_y=center_y, center_z=center_z,
        size_x=size_x, size_y=size_y, size_z=size_z,
        remove_waters=remove_waters,
        rebuild_missing_loops=rebuild_missing_loops,
    )

    create_job(job_id, receptor_source, ligand_source, params.model_dump())

    # Run the actual pipeline as a separate OS process rather than an
    # in-process FastAPI BackgroundTask. Confirmed on a real Render
    # deployment: running the heavy stages (PDBFixer/OpenMM, RDKit/Meeko,
    # Vina) inside the same long-lived uvicorn worker accumulates memory
    # across a job's stages — and across successive jobs, since nothing
    # ever recycles that process — and was OOM-killing a 512MB instance
    # even for a single job with Vina's own thread count already capped at
    # 1. A short-lived subprocess that exits when the job finishes hands
    # 100% of that memory back to the OS unconditionally. Parameters that
    # don't fit on a command line are written to params.json for
    # app/worker.py to read back; job status/results still flow through
    # the shared SQLite job database exactly as before, so nothing in the
    # frontend or polling logic needs to change.
    (job_dir / "params.json").write_text(params.model_dump_json())
    subprocess.Popen(
        [sys.executable, "-m", "app.worker", job_id, str(receptor_pdb), str(ligand_path),
         "1" if ligand_is_smiles else "0"],
        cwd=str(Path(__file__).resolve().parents[2]),  # repo root, so `python -m app...` resolves
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        start_new_session=True,  # detach fully — a client disconnect must not kill the job
    )

    return JobCreateResponse(job_id=job_id, status="queued")


@router.get("/{job_id}/status", response_model=JobStatusResponse)
def job_status(job_id: str):
    job = get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found.")
    return JobStatusResponse(
        job_id=job["job_id"], status=job["status"], stage=job["stage"],
        error_message=job.get("error_message"),
    )


@router.get("/{job_id}/result")
def job_result(job_id: str):
    job = get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found.")
    if job["status"] != "succeeded":
        raise HTTPException(status_code=409, detail=f"Job is '{job['status']}', not yet succeeded.")
    return job["result"]


@router.get("/{job_id}/file")
def job_file(job_id: str, path: str, download: bool = False):
    """Serve a specific output file (pose PDB, receptor PDBQT, render PNG,
    etc.) from within this job's directory only — path traversal is
    blocked by resolving against the job dir and refusing anything outside
    it.

    `download=true` forces a browser download (Content-Disposition:
    attachment) rather than inline rendering — needed for files like
    .pdbqt that have no browser-native viewer and would otherwise just
    print as unstyled text in the tab, which is what a plain FileResponse
    without an explicit disposition does for any content-type the browser
    doesn't recognize.
    """
    job_dir = (settings.JOBS_DIR / job_id).resolve()
    requested = (job_dir / path).resolve()
    if job_dir not in requested.parents and requested != job_dir:
        raise HTTPException(status_code=400, detail="Invalid file path.")
    if not requested.exists():
        raise HTTPException(status_code=404, detail="File not found.")

    media_type = _guess_media_type(requested.suffix)
    if download:
        return FileResponse(
            requested, media_type=media_type, filename=requested.name,
            headers={"Content-Disposition": f'attachment; filename="{requested.name}"'},
        )
    return FileResponse(requested, media_type=media_type)


def _guess_media_type(suffix: str) -> str:
    return {
        ".pdb": "chemical/x-pdb",
        ".pdbqt": "application/octet-stream",
        ".png": "image/png",
        ".sdf": "chemical/x-mdl-sdfile",
        ".log": "text/plain",
        ".csv": "text/csv",
    }.get(suffix.lower(), "application/octet-stream")


@router.get("/recent")
def recent_jobs(limit: int = 20):
    return list_recent_jobs(limit)


async def _resolve_receptor(job_dir, pdb_id, uniprot_id, upload) -> tuple[Path, str]:
    if pdb_id:
        try:
            return fetch_pdb_by_id(pdb_id, job_dir), f"pdb:{pdb_id.upper()}"
        except StructureFetchError as exc:
            raise HTTPException(status_code=404, detail=str(exc)) from exc
    if uniprot_id:
        try:
            return fetch_alphafold_model(uniprot_id, job_dir), f"alphafold:{uniprot_id.upper()}"
        except StructureFetchError as exc:
            raise HTTPException(status_code=404, detail=str(exc)) from exc
    return await _save_upload(upload, job_dir), f"upload:{upload.filename}"


async def _resolve_ligand(job_dir, smiles, upload) -> tuple[Path, bool, str]:
    if smiles:
        smiles_file = job_dir / "ligand_input.smi"
        smiles_file.write_text(smiles.strip())
        return smiles_file, True, f"smiles:{smiles}"
    return await _save_upload(upload, job_dir), False, f"upload:{upload.filename}"


async def _save_upload(upload: UploadFile, job_dir: Path) -> Path:
    dest = job_dir / Path(upload.filename).name
    size = 0
    with dest.open("wb") as out:
        while chunk := await upload.read(1024 * 1024):
            size += len(chunk)
            if size > MAX_UPLOAD_BYTES:
                dest.unlink(missing_ok=True)
                raise HTTPException(
                    status_code=413,
                    detail=f"File exceeds the {settings.MAX_UPLOAD_MB} MB upload limit.",
                )
            out.write(chunk)
    return dest
