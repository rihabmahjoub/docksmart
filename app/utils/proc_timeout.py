"""
Hard timeout wrapper for CPU-bound, in-process calls (PDBFixer/OpenMM) that
have no subprocess boundary of their own to attach a timeout to.

Why this exists: meeko, fpocket, and Vina's CLI fallback all run as
subprocesses, so `subprocess.run(..., timeout=N)` already gives them a real,
enforced wall-clock limit — Python kills the child process outright.
PDBFixer/OpenMM run as plain in-process Python/C calls with no such
boundary; a pathological structure (or just enough missing atoms on a very
large receptor) can block the worker thread indefinitely with nothing able
to interrupt it, which is exactly what produced an unrecoverable 30+ minute
hang requiring a manual redeploy to clear.

The fix: actually run the call in a separate OS process (not a thread —
threads can't be forcibly killed in Python, so a "timeout" on a thread only
stops *waiting* for it, the runaway computation keeps consuming CPU/RAM in
the background regardless) and hard-terminate that process if it exceeds
the deadline. This is real resource reclamation, not just a status update.
"""
from __future__ import annotations

import multiprocessing as mp
from typing import Any, Callable


class HardTimeoutError(RuntimeError):
    pass


def _worker(queue: mp.Queue, func: Callable, args: tuple, kwargs: dict) -> None:
    try:
        result = func(*args, **kwargs)
        queue.put(("ok", result))
    except Exception as exc:  # noqa: BLE001 - report anything back to the parent, don't swallow
        # Put the actual exception object on the queue (not a string
        # decomposition) so the parent can re-raise it with its original
        # type preserved — important because pipeline.py's error handling
        # branches on exception TYPE (ReceptorPrepError vs. everything
        # else) to report a clean stage-specific message. Simple custom
        # exceptions like ReceptorPrepError (a plain RuntimeError subclass
        # with just a message) pickle fine via the default Exception
        # __reduce__; only exotic exception types with unpicklable
        # attributes would need special handling, which none of this
        # codebase's error classes have.
        exc.__cause__ = None  # drop any unpicklable chained-exception context
        queue.put(("error", exc))


def run_with_hard_timeout(
    func: Callable,
    args: tuple = (),
    kwargs: dict | None = None,
    timeout_s: float = 300,
) -> Any:
    """
    Run func(*args, **kwargs) in a separate process. If it doesn't finish
    within timeout_s, the process is forcibly terminated (SIGTERM, then
    SIGKILL if it ignores that) and HardTimeoutError is raised. Any
    exception raised inside func is re-raised in the parent with its
    original type name and message preserved in the error text.

    func and its args/kwargs must be safe to use with the platform's
    multiprocessing start method (fork on Linux, which is what this ships
    on — see Dockerfile). Return value must be picklable (Path/str/dict/
    list of these are all fine; the receptor_prep functions this wraps
    return a Path).
    """
    kwargs = kwargs or {}
    ctx = mp.get_context("fork")
    queue: mp.Queue = ctx.Queue()
    process = ctx.Process(target=_worker, args=(queue, func, args, kwargs), daemon=True)
    process.start()
    process.join(timeout_s)

    if process.is_alive():
        process.terminate()
        process.join(5)
        if process.is_alive():
            process.kill()
            process.join()
        raise HardTimeoutError(
            f"Operation exceeded the {timeout_s:.0f}s limit and was terminated. "
            "This does not necessarily mean anything is broken — very large "
            "structures can genuinely take this long, especially on a "
            "constrained free-tier CPU. Try a smaller/single-chain receptor, "
            "or raise the relevant timeout setting if you control the deployment."
        )

    if queue.empty():
        # The process exited (or was killed by the OS, e.g. OOM) without
        # putting anything on the queue — most commonly an out-of-memory
        # kill on free-tier hosting, not a bug in the function itself.
        raise HardTimeoutError(
            "The worker process exited without returning a result — most likely "
            "it ran out of memory (common on free-tier hosting with a large "
            "structure) rather than a code error. Try a smaller receptor or a "
            "host with more RAM."
        )

    status, payload = queue.get()
    process.join()
    if status == "error":
        raise payload
    return payload
