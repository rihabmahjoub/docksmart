"""
Full docking pipeline orchestration.

This is the single place that sequences every stage described in the
project design (structure QC -> pocket detection -> PDBQT prep -> Vina
-> post-processing -> visualization) and persists progress/results to
the job database as it goes, so the frontend can poll job status.

Each stage is wrapped so a failure anywhere produces a clear, staged
error message rather than an opaque traceback — important for a tool
aimed at non-CLI-expert users.

Two reliability mechanisms live here, both added after real failures on a
free-tier deployment:

  * A process-wide semaphore (_JOB_SLOTS) limits how many jobs run their
    heavy stages concurrently. Confirmed necessary: running two docking
    jobs at once on Render's free tier (512MB RAM) OOM-killed the whole
    instance (exit 137). A second job now waits — reported to the user as
    the "queued" stage — instead of racing the first one for RAM.
  * PDBFixer is run through run_with_hard_timeout() rather than called
    directly, because it's the one stage with no subprocess boundary of
    its own and was confirmed to hang indefinitely (30+ minutes, no
    recovery short of a manual redeploy) on a real structure.
"""
from __future__ import annotations

import logging
import threading
import traceback
from pathlib import Path

from app.config import settings
from app.database import update_job
from app.models import DockingParams, PocketSelectionMode
from app.services import (
    docking_engine,
    ligand_prep,
    pocket_detection,
    postprocessing,
    receptor_prep,
    visualization,
)
from app.utils.proc_timeout import HardTimeoutError, run_with_hard_timeout

_KNOWN_STAGE_ERRORS = (
    receptor_prep.ReceptorPrepError,
    ligand_prep.LigandPrepError,
    pocket_detection.PocketDetectionError,
    docking_engine.DockingError,
    HardTimeoutError,
)

logger = logging.getLogger(__name__)

# Shared across all jobs handled by this process. FastAPI's BackgroundTasks
# run on worker threads, so a plain threading.Semaphore (not asyncio's) is
# the right primitive here — acquire() blocks the worker thread, which is
# fine since it isn't the event loop thread.
_JOB_SLOTS = threading.Semaphore(settings.MAX_CONCURRENT_JOBS)


class PipelineError(RuntimeError):
    def __init__(self, stage: str, message: str):
        self.stage = stage
        super().__init__(message)


def run_pipeline(
    job_id: str,
    receptor_pdb: Path,
    ligand_file: Path,
    ligand_is_smiles: bool,
    params: DockingParams,
) -> None:
    job_dir = settings.JOBS_DIR / job_id
    job_dir.mkdir(parents=True, exist_ok=True)

    # Fail fast on obviously-oversized receptors rather than let PDBFixer/
    # fpocket burn several minutes only to time out anyway. A plain text
    # scan for ATOM/HETATM lines is a cheap proxy for atom count and works
    # before we've even validated the file is well-formed.
    try:
        n_atoms = _count_atom_lines(receptor_pdb)
    except OSError as exc:
        update_job(job_id, status="failed", stage="structure_qc", error=f"Could not read receptor file: {exc}")
        return
    if n_atoms > settings.MAX_RECEPTOR_ATOMS:
        update_job(
            job_id, status="failed", stage="structure_qc",
            error=(
                f"Receptor has ~{n_atoms} atoms, over the {settings.MAX_RECEPTOR_ATOMS} "
                "limit for this deployment. Very large structures can take an "
                "impractically long time on constrained hosting — try a single "
                "chain or the relevant domain only, rather than a full complex."
            ),
        )
        return

    update_job(job_id, status="queued", stage="queued")
    acquired = _JOB_SLOTS.acquire(timeout=settings.JOB_TIMEOUT_SECONDS)
    if not acquired:
        update_job(
            job_id, status="failed", stage="queued",
            error="Timed out waiting for a free job slot — the server is busy with another job.",
        )
        return

    try:
        _run_pipeline_locked(job_id, job_dir, receptor_pdb, ligand_file, ligand_is_smiles, params)
    finally:
        _JOB_SLOTS.release()


def _run_pipeline_locked(
    job_id: str,
    job_dir: Path,
    receptor_pdb: Path,
    ligand_file: Path,
    ligand_is_smiles: bool,
    params: DockingParams,
) -> None:
    try:
        update_job(job_id, status="running", stage="structure_qc")
        fixed_receptor = run_with_hard_timeout(
            receptor_prep.fix_receptor,
            args=(receptor_pdb, job_dir / "receptor_fixed.pdb"),
            kwargs=dict(
                remove_waters=params.remove_waters,
                add_hydrogens=params.add_hydrogens,
                ph=params.ph,
                rebuild_missing_loops=params.rebuild_missing_loops,
            ),
            timeout_s=settings.STRUCTURE_QC_TIMEOUT_SECONDS,
        )

        update_job(job_id, stage="pocket_detection")
        box, pocket_info = _determine_box(fixed_receptor, job_dir, params)

        update_job(job_id, stage="receptor_pdbqt")
        receptor_pdbqt, receptor_dropped_residues = receptor_prep.receptor_to_pdbqt(
            fixed_receptor, job_dir / "receptor.pdbqt"
        )

        update_job(job_id, stage="ligand_pdbqt")
        ligand_prepared = ligand_file
        if ligand_is_smiles:
            ligand_prepared = ligand_prep.smiles_to_3d_sdf(
                ligand_file.read_text().strip(), job_dir / "ligand_3d.sdf"
            )
        ligand_pdbqt = ligand_prep.ligand_to_pdbqt(ligand_prepared, job_dir / "ligand.pdbqt")

        update_job(job_id, stage="docking")
        poses = docking_engine.run_vina(
            receptor_pdbqt,
            ligand_pdbqt,
            box,
            job_dir / "poses.pdbqt",
            job_dir / "vina.log",
            exhaustiveness=params.exhaustiveness,
            num_modes=params.num_modes,
            energy_range=params.energy_range,
            seed=params.seed,
        )
        if not poses:
            raise PipelineError("docking", "Vina completed but produced no poses.")

        update_job(job_id, stage="postprocessing")
        pose_pdb_files = postprocessing.split_poses_to_pdb(
            job_dir / "poses.pdbqt", job_dir / "poses_pdb"
        )
        best_pose_pdb = pose_pdb_files[0] if pose_pdb_files else None

        # Cap at the top 10 poses for reporting — Vina already returns them
        # best-affinity-first, and a 10-row table is what the UI presents.
        # (Docking itself still honors the user's requested num_modes; only
        # the reported/analyzed set is capped here.)
        poses = poses[:10]
        pose_pdb_files = pose_pdb_files[:10]

        pose_results = []
        for pose, pdb_path in zip(poses, pose_pdb_files):
            fp = postprocessing.interaction_fingerprint_prolif(fixed_receptor, pdb_path)
            if fp is None:
                fp = postprocessing.interaction_fingerprint(fixed_receptor, pdb_path)
            # split_poses_to_pdb() always writes pose_{i}.pdbqt alongside
            # pose_{i}.pdb in the same directory (the raw per-pose PDBQT it
            # converts from) — reconstruct that path rather than threading
            # another return value through, since the naming is guaranteed
            # deterministic by that function's own implementation.
            pose_pdbqt_path = pdb_path.with_suffix(".pdbqt")
            pose_results.append(
                {
                    "mode": pose.mode,
                    "affinity_kcal_mol": pose.affinity,
                    "rmsd_lb": pose.rmsd_lb,
                    "rmsd_ub": pose.rmsd_ub,
                    "rmsd_note": pose.rmsd_note,
                    "interaction_summary": fp,
                    "pose_pdb": str(pdb_path),
                    "pose_pdbqt": str(pose_pdbqt_path) if pose_pdbqt_path.exists() else None,
                }
            )

        render_png = None
        if best_pose_pdb is not None:
            update_job(job_id, stage="rendering")
            try:
                render_png = visualization.render_complex(
                    fixed_receptor, best_pose_pdb, box, job_dir / "render.png"
                )
            except visualization.VisualizationError as exc:
                logger.warning("Visualization skipped: %s", exc)

        result = {
            "poses": pose_results,
            "pocket_info": pocket_info,
            "best_pose_pdb": str(best_pose_pdb) if best_pose_pdb else None,
            "receptor_pdb": str(fixed_receptor),
            "receptor_pdbqt": str(receptor_pdbqt),
            "render_png": str(render_png) if render_png else None,
            "receptor_dropped_residues": receptor_dropped_residues,
        }
        update_job(job_id, status="succeeded", stage="done", result=result)
        logger.info("Job %s completed successfully", job_id)

    except PipelineError as exc:
        logger.error("Job %s failed at stage %s: %s", job_id, exc.stage, exc)
        update_job(job_id, status="failed", stage=exc.stage, error=str(exc))
    except _KNOWN_STAGE_ERRORS as exc:
        # A recognized service-level error (missing dependency, malformed
        # input, external tool failure, hard timeout) raised directly by a
        # stage function rather than wrapped in PipelineError — still a
        # clean, actionable message, not a crash. The job's `stage` column
        # already reflects the last update_job() call made before this
        # stage ran.
        logger.error("Job %s failed: %s", job_id, exc)
        update_job(job_id, status="failed", error=str(exc))
    except Exception as exc:  # noqa: BLE001 - genuinely unexpected: keep the traceback in logs
        logger.error("Job %s crashed unexpectedly:\n%s", job_id, traceback.format_exc())
        update_job(job_id, status="failed", error=f"Unexpected internal error: {exc}")


def _count_atom_lines(pdb_file: Path) -> int:
    count = 0
    with open(pdb_file, errors="ignore") as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                count += 1
    return count


def _determine_box(fixed_receptor: Path, job_dir: Path, params: DockingParams):
    if params.pocket_mode == PocketSelectionMode.MANUAL:
        if None in (params.center_x, params.center_y, params.center_z):
            raise PipelineError("pocket_detection", "Manual mode requires center_x/y/z.")
        box = pocket_detection.manual_box(
            params.center_x, params.center_y, params.center_z,
            params.size_x, params.size_y, params.size_z,
        )
        return box, {"mode": "manual"}

    if params.pocket_mode == PocketSelectionMode.COCRYSTAL_LIGAND:
        # NOTE: resname currently expected to be passed via params in a future
        # revision; for v1 this mode requires the caller to pre-resolve it.
        raise PipelineError(
            "pocket_detection",
            "Co-crystal ligand mode requires a HETATM residue name — not yet wired to the API.",
        )

    # default: AUTO_FPOCKET
    try:
        pockets, out_dir = pocket_detection.run_fpocket(fixed_receptor, job_dir / "fpocket")
    except pocket_detection.PocketDetectionError as exc:
        raise PipelineError("pocket_detection", str(exc)) from exc

    if not pockets:
        raise PipelineError("pocket_detection", "fpocket found no candidate pockets.")

    top = pockets[0]
    pocket_pdb = out_dir / "pockets" / f"pocket{top['pocket_id']}_atm.pdb"
    if not pocket_pdb.exists():
        raise PipelineError("pocket_detection", f"Expected pocket file missing: {pocket_pdb}")

    box = pocket_detection.box_from_fpocket_pocket(pocket_pdb)
    pocket_info = {
        "mode": "auto_fpocket",
        "selected_pocket": top,
        "all_pockets_ranked": pockets[:10],
    }
    return box, pocket_info
