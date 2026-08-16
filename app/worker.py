"""
Standalone entry point for running a single docking job in its own OS
process, invoked as:

    python -m app.worker <job_id> <receptor_pdb_path> <ligand_path> <ligand_is_smiles: 0|1>

Parameters that don't fit cleanly on a command line (the full DockingParams
object) are read from `<job_dir>/params.json`, written by the API handler
before this process is spawned.

Why this exists: see app/utils/proc_lock.py's module docstring. In short,
running the pipeline as an in-process FastAPI BackgroundTask meant every
job's memory (PDBFixer/OpenMM, then RDKit/Meeko, then Vina) accumulated
inside the single long-lived uvicorn worker and was never returned to the
OS — confirmed as the cause of a recurring Render OOM restart even after
capping Vina's own thread count. Running each job as a separate process
that exits when the job finishes means the OS reclaims 100% of that job's
memory unconditionally, regardless of what any library inside the pipeline
does internally.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

from app.models import DockingParams
from app.services.pipeline import run_pipeline


def main() -> int:
    if len(sys.argv) != 5:
        print(f"usage: python -m app.worker <job_id> <receptor_pdb> <ligand_path> <ligand_is_smiles:0|1>",
              file=sys.stderr)
        return 2

    job_id, receptor_pdb, ligand_path, ligand_is_smiles_flag = sys.argv[1:5]
    job_dir = Path(receptor_pdb).parent
    params_path = job_dir / "params.json"
    params = DockingParams(**json.loads(params_path.read_text()))

    run_pipeline(
        job_id=job_id,
        receptor_pdb=Path(receptor_pdb),
        ligand_file=Path(ligand_path),
        ligand_is_smiles=bool(int(ligand_is_smiles_flag)),
        params=params,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
