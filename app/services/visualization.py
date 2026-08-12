"""
Publication-style rendering via headless PyMOL (`pymol2` module, installed
with `pip install pymol-open-source` or `conda install pymol-open-source`).

Produces a static PNG (cartoon receptor + docking grid box + top pose
sticks + nearby contact residues) for the report/download. This
complements the in-browser interactive NGL.js viewer (app/static/js/
viewer.js), which renders directly from the raw PDB files and does not
depend on this service.
"""
from __future__ import annotations

import logging
from pathlib import Path

from app.services.pocket_detection import GridBox

logger = logging.getLogger(__name__)


class VisualizationError(RuntimeError):
    pass


def render_complex(
    receptor_pdb: Path,
    pose_pdb: Path,
    box: GridBox,
    output_png: Path,
    width: int = 1200,
    height: int = 900,
) -> Path:
    try:
        import pymol2
    except ImportError as exc:
        raise VisualizationError(
            "pymol2 is not installed. Install with: pip install pymol-open-source"
        ) from exc

    output_png.parent.mkdir(parents=True, exist_ok=True)

    with pymol2.PyMOL() as p:
        cmd = p.cmd
        cmd.load(str(receptor_pdb), "receptor")
        cmd.load(str(pose_pdb), "pose")

        cmd.hide("everything")
        cmd.show("cartoon", "receptor")
        cmd.color("gray80", "receptor")
        cmd.show("sticks", "pose")
        # NOTE: cmd.util.cbag("pose") was tried first but confirmed broken in
        # the pymol-open-source-whl build used here — util.cbag() internally
        # reaches for a different (uninitialized) global PyMOL instance than
        # the one pymol2.PyMOL() creates, raising "RuntimeError: Missing
        # PyMOL instance". Plain cmd.color() calls achieve the identical
        # visual result (green carbons, CPK-colored heteroatoms) and were
        # confirmed working directly against this same instance.
        cmd.color("green", "pose and elem C")
        cmd.color("atomic", "pose and not elem C")

        # Show side chains of residues actually contacting the pose (not the
        # whole box volume) — more informative for a methods/results figure
        # than an empty region of protein.
        cmd.select("pocket_env", "byres (receptor within 4.5 of pose)")
        cmd.show("sticks", "pocket_env and not (name C+N+O and not resn PRO)")
        cmd.color("skyblue", "pocket_env and elem C")

        _draw_grid_box(cmd, box)

        cmd.bg_color("white")
        cmd.set("ray_opaque_background", 0)
        cmd.set("cartoon_transparency", 0.15)
        cmd.orient("pose")
        # NOTE: `cmd.zoom("pose or box_outline", ...)` was tried first but
        # confirmed to raise a Selector-Error — PyMOL's boolean selection
        # algebra ("or") doesn't apply between an atom selection and a CGO
        # object name (box_outline isn't atoms). Non-fatal (PyMOL logs the
        # error and continues with a fallback view), but the camera framing
        # was wrong. Since the search box always encloses the pose region
        # by construction, a generous buffer on "pose" alone reliably keeps
        # the drawn box in frame without needing a mixed-type selection.
        cmd.zoom("pose", buffer=10)
        cmd.set("ray_trace_mode", 1)
        cmd.ray(width, height)
        cmd.png(str(output_png), dpi=300)

    logger.info("Rendered complex image: %s", output_png)
    return output_png


def _draw_grid_box(cmd, box: GridBox) -> None:
    """Draw the Vina search-space box as a wireframe CGO cube, so the figure
    documents which region of the receptor was actually searched — directly
    relevant for showing how the fpocket/co-crystal/manual box was chosen."""
    from pymol import cgo

    cx, cy, cz = box.center
    sx, sy, sz = (s / 2 for s in box.size)
    corners = [
        (cx + dx * sx, cy + dy * sy, cz + dz * sz)
        for dx in (-1, 1) for dy in (-1, 1) for dz in (-1, 1)
    ]
    edges = [
        (0, 1), (0, 2), (0, 4), (3, 1), (3, 2), (3, 7),
        (5, 1), (5, 4), (5, 7), (6, 2), (6, 4), (6, 7),
    ]
    obj = [cgo.COLOR, 1.0, 0.4, 0.0, cgo.LINEWIDTH, 2.0]
    for a, b in edges:
        obj += [cgo.BEGIN, cgo.LINES, cgo.VERTEX, *corners[a], cgo.VERTEX, *corners[b], cgo.END]
    cmd.load_cgo(obj, "box_outline")


def export_session(receptor_pdb: Path, pose_pdb: Path, output_pse: Path) -> Path:
    """Export a .pse PyMOL session so users can continue exploring locally."""
    try:
        import pymol2
    except ImportError as exc:
        raise VisualizationError(
            "pymol2 is not installed. Install with: pip install pymol-open-source"
        ) from exc

    output_pse.parent.mkdir(parents=True, exist_ok=True)
    with pymol2.PyMOL() as p:
        cmd = p.cmd
        cmd.load(str(receptor_pdb), "receptor")
        cmd.load(str(pose_pdb), "pose")
        cmd.show_as("cartoon", "receptor")
        cmd.show_as("sticks", "pose")
        cmd.save(str(output_pse))
    return output_pse
