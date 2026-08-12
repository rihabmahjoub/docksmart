"""
Binding-site / grid-box determination.

This is DockSmart's main methodological differentiator versus a "paste
in x/y/z coordinates" web form: pockets are detected and ranked
automatically, the user sees *why* a box was chosen, and can override it.

Three strategies, selected via DockingParams.pocket_mode:

  AUTO_FPOCKET       run fpocket, rank cavities by druggability score,
                     build a box around the top pocket's alpha spheres
  COCRYSTAL_LIGAND   if the receptor came from a PDB entry with a bound
                     HETATM ligand, center the box on that ligand
                     (classic co-crystallized re-docking / validation mode)
  MANUAL             user-supplied center + size, used as-is
"""
from __future__ import annotations

import logging
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from app.config import settings

logger = logging.getLogger(__name__)


class PocketDetectionError(RuntimeError):
    pass


@dataclass
class GridBox:
    center: tuple[float, float, float]
    size: tuple[float, float, float]
    source: str                      # "fpocket" | "cocrystal_ligand" | "manual"
    pocket_rank: Optional[int] = None
    druggability_score: Optional[float] = None
    pocket_score: Optional[float] = None


def run_fpocket(receptor_pdb: Path, workdir: Path) -> list[dict]:
    """Run fpocket and return ranked pocket summaries (best druggability first)."""
    workdir.mkdir(parents=True, exist_ok=True)
    local_copy = workdir / receptor_pdb.name
    if local_copy.resolve() != receptor_pdb.resolve():
        local_copy.write_bytes(receptor_pdb.read_bytes())

    cmd = [settings.FPOCKET_BINARY, "-f", str(local_copy)]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=settings.FPOCKET_TIMEOUT_SECONDS)
    except FileNotFoundError as exc:
        raise PocketDetectionError(
            f"'{settings.FPOCKET_BINARY}' not found. Install fpocket "
            "(https://github.com/Discngine/fpocket) or switch to manual/co-crystal mode."
        ) from exc
    except subprocess.CalledProcessError as exc:
        raise PocketDetectionError(f"fpocket failed:\n{exc.stderr}") from exc
    except subprocess.TimeoutExpired as exc:
        raise PocketDetectionError(
            f"fpocket exceeded the {settings.FPOCKET_TIMEOUT_SECONDS}s timeout "
            "(usually a slow CPU allocation, not a stuck process) — try raising "
            "DOCKSMART_FPOCKET_TIMEOUT, or a smaller/single-chain receptor."
        ) from exc

    out_dir = workdir / f"{local_copy.stem}_out"
    info_file = out_dir / f"{local_copy.stem}_info.txt"
    if not info_file.exists():
        raise PocketDetectionError("fpocket ran but no _info.txt output was found.")

    return _parse_fpocket_info(info_file), out_dir


import re

_POCKET_HEADER_RE = re.compile(r"^Pocket\s+(\d+)\s*:\s*$")


def _parse_fpocket_info(info_file: Path) -> list[dict]:
    """Parse fpocket's *_info.txt into a list of pocket dicts, ranked as reported.

    NOTE: a header line looks like 'Pocket 1 :' (just an index), while
    per-pocket fields include lines like 'Pocket Score : 22.1' — both start
    with the literal word 'Pocket', so a naive `startswith("Pocket")` check
    misclassifies the 'Pocket Score' field as a new pocket header. The
    regex below only matches the header form (word, integer, colon, nothing
    else) to disambiguate the two.
    """
    pockets: list[dict] = []
    current: dict = {}
    for raw_line in info_file.read_text().splitlines():
        line = raw_line.strip()
        header_match = _POCKET_HEADER_RE.match(line)
        if header_match:
            if current:
                pockets.append(current)
            current = {"pocket_id": header_match.group(1)}
        elif ":" in line and current:
            key, _, value = line.partition(":")
            key = key.strip().lower().replace(" ", "_")
            value = value.strip()
            try:
                current[key] = float(value)
            except ValueError:
                current[key] = value
    if current:
        pockets.append(current)

    pockets.sort(key=lambda p: p.get("druggability_score", 0.0), reverse=True)
    return pockets


def box_from_fpocket_pocket(pocket_pdb: Path, padding: float = None) -> GridBox:
    """Build a grid box from a pocket's alpha-sphere PDB (pocketN_atm.pdb / _vert.pqr)."""
    padding = padding if padding is not None else settings.DEFAULT_BOX_PADDING_A
    coords = _read_pdb_coords(pocket_pdb)
    if not coords:
        raise PocketDetectionError(f"No atoms found in pocket file {pocket_pdb}")
    return _box_from_coords(coords, padding, source="fpocket")


def box_from_cocrystal_ligand(receptor_pdb: Path, resname: str, padding: float = None) -> GridBox:
    """Build a grid box centered on a named HETATM residue in the original PDB file."""
    padding = padding if padding is not None else settings.DEFAULT_BOX_PADDING_A
    coords = [
        c for c in _read_pdb_coords(receptor_pdb, hetatm_only=True, resname=resname)
    ]
    if not coords:
        raise PocketDetectionError(
            f"No HETATM residue named '{resname}' found in {receptor_pdb.name}."
        )
    return _box_from_coords(coords, padding, source="cocrystal_ligand")


def _read_pdb_coords(
    pdb_file: Path, hetatm_only: bool = False, resname: Optional[str] = None
) -> list[tuple[float, float, float]]:
    coords = []
    record_types = ("HETATM",) if hetatm_only else ("ATOM", "HETATM")
    for line in pdb_file.read_text(errors="ignore").splitlines():
        if not line.startswith(record_types):
            continue
        if resname and line[17:20].strip() != resname:
            continue
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue
        coords.append((x, y, z))
    return coords


def _box_from_coords(
    coords: list[tuple[float, float, float]], padding: float, source: str
) -> GridBox:
    xs, ys, zs = zip(*coords)
    center = ((min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2, (min(zs) + max(zs)) / 2)
    size = (
        max(max(xs) - min(xs) + 2 * padding, 15.0),
        max(max(ys) - min(ys) + 2 * padding, 15.0),
        max(max(zs) - min(zs) + 2 * padding, 15.0),
    )
    return GridBox(center=center, size=size, source=source)


def manual_box(cx: float, cy: float, cz: float, sx: float, sy: float, sz: float) -> GridBox:
    return GridBox(center=(cx, cy, cz), size=(sx, sy, sz), source="manual")
