"""
Cross-process job-slot lock.

`pipeline.py` previously used a plain `threading.Semaphore` to cap how many
docking jobs run their heavy stages at once — this worked as long as every
job ran as a FastAPI BackgroundTask *inside the single long-lived uvicorn
worker process*. That in-process model turned out to be the real cause of a
recurring Render OOM restart: within one worker process, memory from
PDBFixer/OpenMM, then RDKit/Meeko, then Vina accumulates across a job's
stages (and across successive jobs, since nothing ever recycles the
process) and is never handed back to the OS, so a single 512MB-RAM instance
can be OOM-killed by one job's cumulative footprint alone — no concurrency
required. The fix is to run each job in its own short-lived OS subprocess
(see app/worker.py) so the OS fully reclaims all memory when a job finishes.

That move breaks `threading.Semaphore`, though: each subprocess gets its
own fresh Python interpreter and its own semaphore object, so it no longer
coordinates anything across processes. This module replaces it with an
`fcntl.flock`-based lock on a small set of on-disk lock files — flock is
enforced by the OS kernel, so it correctly serializes access across
however many separate subprocesses show up, including if the previous
subprocess crashed without cleanly releasing (the OS releases the lock
automatically when the holding process exits, for any reason).

Implemented as `settings.MAX_CONCURRENT_JOBS` separate lock files rather
than one, so a future increase of that setting (e.g. on a larger paid
instance) works without changing this module — each job takes whichever
slot file it can lock first.
"""
from __future__ import annotations

import fcntl
import time
from contextlib import contextmanager
from typing import Iterator

from app.config import settings


@contextmanager
def job_slot(poll_interval_s: float = 1.0, timeout_s: float | None = None) -> Iterator[None]:
    """Block until one of the MAX_CONCURRENT_JOBS slot files can be
    exclusively locked, hold it for the duration of the `with` block, then
    release it. Raises TimeoutError if timeout_s elapses first."""
    settings.DATA_DIR.mkdir(parents=True, exist_ok=True)
    n_slots = max(1, settings.MAX_CONCURRENT_JOBS)
    lock_paths = [settings.DATA_DIR / f"job_slot_{i}.lock" for i in range(n_slots)]

    start = time.monotonic()
    handle = None
    try:
        while handle is None:
            for path in lock_paths:
                f = open(path, "w")
                try:
                    fcntl.flock(f.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                    handle = f
                    break
                except BlockingIOError:
                    f.close()
            if handle is None:
                if timeout_s is not None and (time.monotonic() - start) > timeout_s:
                    raise TimeoutError("Timed out waiting for a free job slot.")
                time.sleep(poll_interval_s)
        yield
    finally:
        if handle is not None:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
            handle.close()
