"""
DockSmart configuration.

All paths / tunables live here so the same codebase runs unchanged on a
local machine, PythonAnywhere, or Render — only environment variables
change between deployments.
"""
from __future__ import annotations

import os
from pathlib import Path


class Settings:
    # --- Identity -----------------------------------------------------
    APP_NAME: str = "DockSmart"
    APP_VERSION: str = "0.1.0"

    # --- Paths ----------------------------------------------------------
    BASE_DIR: Path = Path(__file__).resolve().parent.parent
    DATA_DIR: Path = Path(os.getenv("DOCKSMART_DATA_DIR", BASE_DIR / "data"))
    JOBS_DIR: Path = DATA_DIR / "jobs"
    DB_PATH: Path = DATA_DIR / "docksmart.db"

    # --- External binaries (override via env if not on PATH) -----------
    VINA_BINARY: str = os.getenv("VINA_BINARY", "vina")
    FPOCKET_BINARY: str = os.getenv("FPOCKET_BINARY", "fpocket")

    # --- Docking defaults -------------------------------------------------
    # Lowered from the AutoDock Vina "textbook" default of 8: on a
    # constrained free-tier CPU (e.g. Render's 0.1 CPU), exhaustiveness=8
    # measurably contributes to "this took too long" complaints. 4 is still
    # a reasonable search thoroughness for most targets; users who want the
    # more exhaustive search can raise it via the advanced parameters panel.
    DEFAULT_EXHAUSTIVENESS: int = int(os.getenv("DOCKSMART_EXHAUSTIVENESS", 4))
    DEFAULT_NUM_MODES: int = int(os.getenv("DOCKSMART_NUM_MODES", 10))
    DEFAULT_BOX_PADDING_A: float = 5.0  # angstroms around a co-crystal ligand / pocket

    # Vina's own thread count for a SINGLE job — separate from
    # MAX_CONCURRENT_JOBS below, which only bounds how many jobs run at
    # once. Left unset, Vina calls the host's CPU count (e.g. via
    # std::thread::hardware_concurrency / Python's os.cpu_count()) and
    # spawns that many search threads. On Render, the container reports the
    # *host machine's* full core count even though the container itself is
    # only entitled to a small memory budget (512MB on free tier) — so a
    # single job can OOM (exit 137) purely from Vina's own thread-local
    # search-state memory, with no other job running at all. Confirmed by a
    # real Render OOM restart on a single 1IEP docking job. Default of 1 is
    # deliberately conservative for free-tier memory safety; raise via env
    # var on a host with more RAM if faster wall-clock time is wanted.
    VINA_CPU: int = int(os.getenv("DOCKSMART_VINA_CPU", 1))

    # --- Per-stage subprocess timeouts -----------------------------------
    # These were originally hardcoded at 120s/180s, which is fine on a normal
    # CPU but genuinely too short on Render free tier (0.1 CPU): Meeko's
    # per-residue template matching and fpocket's cavity search are real
    # compute, not something hung/stuck — they just need more wall-clock
    # time when the CPU allocation is this small. Raise these via env vars
    # if you still hit timeouts on a larger receptor; there's no way to
    # know the right value in advance since it scales with both receptor
    # size and whatever CPU share the host happens to give you at that
    # moment.
    MEEKO_TIMEOUT_SECONDS: int = int(os.getenv("DOCKSMART_MEEKO_TIMEOUT", 600))
    FPOCKET_TIMEOUT_SECONDS: int = int(os.getenv("DOCKSMART_FPOCKET_TIMEOUT", 600))
    STRUCTURE_QC_TIMEOUT_SECONDS: int = int(os.getenv("DOCKSMART_STRUCTURE_QC_TIMEOUT", 300))
    # PDBFixer runs in-process (no subprocess of its own), so unlike the
    # tools above it had NO enforced timeout at all until this was added —
    # confirmed to hang indefinitely (30+ min, no recovery short of a
    # redeploy) on at least one real structure. See app/utils/proc_timeout.py.

    # Receptors above this atom count are rejected up front rather than
    # attempted — large complexes can make PDBFixer/fpocket genuinely take
    # a very long time even before hitting any timeout, which is a bad
    # experience (waiting minutes just to find out it times out anyway).
    # 15000 atoms is roughly a ~2000-residue single-chain protein or a
    # multi-chain complex of that combined size — generous for a typical
    # docking target, but well short of e.g. a full ribosome structure.
    MAX_RECEPTOR_ATOMS: int = int(os.getenv("DOCKSMART_MAX_RECEPTOR_ATOMS", 15000))

    # How many docking jobs may run their heavy stages at once. Confirmed by
    # a real Render free-tier crash (exit 137 = OOM kill) that running two
    # jobs concurrently exhausts the 512MB RAM budget. This is enforced by a
    # process-wide semaphore in pipeline.py, not just documented — a second
    # job will wait (and report "queued") rather than run concurrently and
    # risk taking the whole instance down with it.
    MAX_CONCURRENT_JOBS: int = int(os.getenv("DOCKSMART_MAX_CONCURRENT_JOBS", 1))

    # --- Limits (important on free-tier hosting: PythonAnywhere/Render) -
    MAX_UPLOAD_MB: int = int(os.getenv("DOCKSMART_MAX_UPLOAD_MB", 15))
    JOB_TIMEOUT_SECONDS: int = int(os.getenv("DOCKSMART_JOB_TIMEOUT", 900))  # 15 min hard cap

    # --- External services ----------------------------------------------
    RCSB_DOWNLOAD_URL: str = "https://files.rcsb.org/download/{pdb_id}.pdb"
    ALPHAFOLD_MODEL_URL: str = "https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"

    def ensure_dirs(self) -> None:
        self.DATA_DIR.mkdir(parents=True, exist_ok=True)
        self.JOBS_DIR.mkdir(parents=True, exist_ok=True)


settings = Settings()
settings.ensure_dirs()
