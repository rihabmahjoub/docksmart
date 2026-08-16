"""
Receptor preparation.

Two steps, each swappable independently:

 1. Structural repair with PDBFixer (missing residues/atoms, non-standard
    residues, optional hydrogen addition at a given pH, water/heteroatom
    stripping) — plus alternate-location (altLoc) collapsing, which
    PDBFixer does NOT do on its own (see _collapse_altlocs).
 2. Conversion to AutoDock-flavoured PDBQT with Meeko, which has replaced
    the old (unmaintained) AutoDockTools `prepare_receptor4.py` scripts
    used in most legacy docking tutorials.

Both steps degrade gracefully with a clear error if the optional
dependency isn't installed, so the rest of the app can still be explored
/ demoed without a full scientific-stack install.
"""
from __future__ import annotations

import logging
import re
from pathlib import Path

from app.config import settings

logger = logging.getLogger(__name__)


class ReceptorPrepError(RuntimeError):
    pass


# Matches lines like:
#   matched with excess inter-residue bond(s): A:779
# emitted by Meeko's Polymer._build_padded_mols() when it perceives more
# inter-residue bonds for a residue pair than any known padding template
# accounts for (e.g. a distance-perceived contact that isn't really a
# disulfide). This is a DIFFERENT failure mode from the "did not match
# template" case that --allow_bad_res already covers: it happens one step
# later, during bond padding, so --allow_bad_res alone does not stop it
# from aborting the whole receptor.
_EXCESS_BOND_RE = re.compile(r"excess inter-residue bond\(s\):\s*([A-Za-z0-9]+:\d+)")


def _extract_offending_residues(stderr: str) -> list[str]:
    """Pull residue IDs (e.g. 'A:779') out of a Meeko excess-bond failure."""
    return sorted(set(_EXCESS_BOND_RE.findall(stderr or "")))


def _collapse_altlocs(pdb_text: str) -> str:
    """Keep exactly one alternate-location record per atom, blank its altLoc
    field, and drop the rest.

    Confirmed necessary by reproduction on a real structure (PDB 5K5X):
    PDBFixer does NOT resolve alternate locations on its own — despite this
    module's docstring having previously (incorrectly) implied it does —
    it passes duplicate ATOM records straight through when a residue has
    more than one refined conformation (altLoc 'A'/'B'/...). Meeko then
    tries to build a single RDKit molecule per residue from the raw PDB
    atom records; with both altLoc copies of an atom present, its
    altloc-handling path (_aux_altloc_mol_build) can end up perceiving an
    extra bond that pushes an atom's computed valence above what RDKit's
    sanitizer allows (observed: 'Explicit valence for atom # 6 C, 5, is
    greater than permitted'), aborting the whole receptor with an
    unhandled RDKit exception — not caught by --allow_bad_res, which only
    covers template-matching failures, not sanitization failures.

    For each (chain, residue, atom name) position, we keep the
    highest-occupancy record (ties broken by taking the first — usually
    'A' — alphabetically), matching standard PDB-cleaning practice used
    ahead of docking preparation.
    """
    best: dict[tuple, tuple] = {}  # (chain, resnum, icode, name) -> (occupancy, altloc, line_index)
    lines = pdb_text.splitlines()
    for idx, line in enumerate(lines):
        if not line.startswith(("ATOM", "HETATM")) or len(line) < 27:
            continue
        altloc = line[16]
        if altloc == " ":
            continue  # no alternate location on this atom; nothing to collapse
        key = (line[21], line[22:26], line[26], line[12:16])
        try:
            occupancy = float(line[54:60])
        except ValueError:
            occupancy = 1.0
        if key not in best or occupancy > best[key][0]:
            best[key] = (occupancy, altloc, idx)

    if not best:
        return pdb_text  # nothing had an altLoc — file unchanged

    winning_indices = {v[2] for v in best.values()}
    out_lines = []
    dropped = 0
    for idx, line in enumerate(lines):
        if not line.startswith(("ATOM", "HETATM")) or len(line) < 27 or line[16] == " ":
            out_lines.append(line)
            continue
        if idx in winning_indices:
            out_lines.append(line[:16] + " " + line[17:])  # blank the altLoc column
        else:
            dropped += 1  # a losing altLoc copy — drop it
    logger.info("Collapsed alternate locations: kept %d atoms, dropped %d duplicate altLoc records",
                len(best), dropped)
    return "\n".join(out_lines) + "\n"


def fix_receptor(
    input_pdb: Path,
    output_pdb: Path,
    *,
    remove_waters: bool = True,
    add_hydrogens: bool = True,
    ph: float = 7.4,
    rebuild_missing_loops: bool = False,
) -> Path:
    """
    Repair a receptor structure with PDBFixer and write a clean PDB.

    By default (rebuild_missing_loops=False), unresolved *loops* (gaps
    reported by findMissingResidues, typically flexible/disordered regions
    with no electron density) are left as-is rather than rebuilt: for rigid
    -receptor Vina docking, modeling in a loop the crystal structure never
    resolved risks introducing geometry that biases the grid box / clashes
    with the ligand without experimental support. Missing atoms *within*
    residues that ARE resolved (e.g. a truncated side chain) are always
    completed regardless of this flag — that is standard, low-risk repair.
    Set rebuild_missing_loops=True only if the disordered region is known to
    be structurally important (e.g. it lines the binding site) and you
    accept the modeling uncertainty this introduces.
    """
    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile
    except ImportError as exc:
        raise ReceptorPrepError(
            "PDBFixer/OpenMM is not installed. Install with: pip install pdbfixer openmm"
        ) from exc

    # Collapse alternate locations before anything else touches the
    # structure — see _collapse_altlocs' docstring. Written to a sibling
    # temp file rather than mutating input_pdb in place, since callers may
    # reasonably expect the original upload/download to be left untouched.
    raw_text = input_pdb.read_text()
    cleaned_text = _collapse_altlocs(raw_text)
    if cleaned_text is not raw_text:
        altloc_cleaned_pdb = input_pdb.with_name(input_pdb.stem + "_noaltloc.pdb")
        altloc_cleaned_pdb.write_text(cleaned_text)
        input_pdb = altloc_cleaned_pdb

    fixer = PDBFixer(filename=str(input_pdb))
    fixer.findMissingResidues()
    if not rebuild_missing_loops:
        # Keep the dict (required by addMissingAtoms' bookkeeping) but empty
        # it so no whole residues/loops get built in — only atom-level gaps
        # within existing residues will be filled below.
        fixer.missingResidues = {}

    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()

    if remove_waters:
        fixer.removeHeterogens(keepWater=False)
    else:
        fixer.removeHeterogens(keepWater=True)

    # Confirmed by reproduction: a receptor file with no standard ATOM
    # records (e.g. an entire structure deposited/exported as HETATM only,
    # or a ligand file mistakenly uploaded as the receptor) is stripped to
    # zero atoms by removeHeterogens() above. OpenMM then crashes deep
    # inside addMissingHydrogens() with an opaque "Cannot create a Context
    # for a System with no particles" error that gives the user no idea
    # what actually went wrong. Fail clearly here instead, before wasting
    # time on the rest of the pipeline.
    n_atoms = sum(1 for _ in fixer.topology.atoms())
    if n_atoms == 0:
        raise ReceptorPrepError(
            "The receptor file contains no usable protein atoms after cleanup "
            "(0 atoms remained after removing waters/heteroatoms). This usually "
            "means the uploaded file is not a standard protein PDB — check that "
            "it contains ATOM records (not just HETATM), and that you didn't "
            "accidentally upload a ligand or cofactor file as the receptor."
        )

    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    if add_hydrogens:
        # NOTE: addMissingHydrogens() uses a fixed-pH heuristic (not a true
        # pKa-prediction / titration method like PROPKA or H++). It is a
        # reasonable default for a physiological pH but is not a substitute
        # for dedicated protonation-state prediction on titratable
        # catalytic/binding-site residues if that precision matters for a
        # given target — flag this as a known limitation in the manuscript.
        fixer.addMissingHydrogens(ph)

    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    with open(output_pdb, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)

    logger.info("Receptor repaired: %s -> %s", input_pdb, output_pdb)
    return output_pdb


def receptor_to_pdbqt(input_pdb: Path, output_pdbqt: Path) -> tuple[Path, list[str]]:
    """Convert a prepared receptor PDB to PDBQT using Meeko.

    Returns (output_pdbqt, dropped_residues) — dropped_residues is normally
    empty; see the retry logic below for when it isn't. Returned rather than
    stashed on a module/function attribute so this stays safe if
    settings.MAX_CONCURRENT_JOBS > 1 (multiple jobs calling this
    concurrently on different threads).
    """
    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy
        from rdkit import Chem
    except ImportError as exc:
        raise ReceptorPrepError(
            "Meeko/RDKit is not installed. Install with: pip install meeko rdkit"
        ) from exc

    # Meeko's receptor prep path (mk_prepare_receptor.py) is normally invoked
    # as a CLI; we shell out to it for robustness across Meeko versions
    # rather than depending on internal APIs that change between releases.
    import subprocess

    output_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    base_cmd = [
        "mk_prepare_receptor.py",
        "--read_pdb", str(input_pdb),
        "-o", str(output_pdbqt.with_suffix("")),
        "-p",  # write pdbqt
        # NOTE: deliberately NOT passing -v/--write_vina_box here — that flag
        # requires --box_center/--box_size to already be known (it writes a
        # Vina config file), which we don't have yet at this stage: the grid
        # box is computed separately from fpocket/manual/co-crystal mode in
        # pocket_detection.py and passed directly to run_vina(), not through
        # a Meeko-generated config file. Confirmed by testing: passing -v
        # without box params makes mk_prepare_receptor.py exit non-zero.
        "--allow_bad_res",
        # Confirmed necessary in practice: Meeko infers inter-residue bonds
        # (including disulfides) from atomic distance, then tries to match
        # each bonded pair against a padding template. On real PDB structures
        # this sometimes finds a bond geometry that doesn't cleanly match any
        # known template (e.g. a distance-perceived S-S contact that isn't
        # really a disulfide, or a partially-resolved residue near one) and
        # raises RuntimeError, aborting the ENTIRE receptor rather than just
        # that residue. --allow_bad_res makes Meeko skip the offending
        # residue(s) with a warning instead of failing the whole conversion —
        # this is Meeko's own documented mechanism for this class of error,
        # not a DockSmart workaround. The skipped-residue warnings land in
        # stderr and are surfaced below if the command still fails for an
        # unrelated reason.
    ]

    dropped_residues: list[str] = []
    cmd = list(base_cmd)
    # Two attempts max: the first is the normal run; if it fails specifically
    # on the "excess inter-residue bond" padding error — which
    # --allow_bad_res does NOT cover, since it fires after template matching
    # already succeeded — retry once with those exact residues explicitly
    # deleted. This is safe for docking as long as they aren't residues that
    # line the binding site actually being used; they are almost always
    # spurious distance-perceived contacts (crystal packing, altloc
    # duplication) rather than real chemistry. If deletion still fails, the
    # original error is raised rather than silently masked.
    for attempt in range(2):
        try:
            subprocess.run(
                cmd, check=True, capture_output=True, text=True,
                timeout=settings.MEEKO_TIMEOUT_SECONDS,
            )
            break
        except subprocess.CalledProcessError as exc:
            offending = _extract_offending_residues(exc.stderr)
            if attempt == 0 and offending:
                dropped_residues = offending
                logger.warning(
                    "Meeko excess inter-residue bond(s) at %s — retrying with "
                    "these residues deleted from the receptor.",
                    ", ".join(offending),
                )
                cmd = base_cmd + ["--delete_residues", ",".join(offending)]
                continue
            raise ReceptorPrepError(
                f"mk_prepare_receptor.py failed:\n{exc.stderr}"
            ) from exc
        except subprocess.TimeoutExpired as exc:
            raise ReceptorPrepError(
                f"mk_prepare_receptor.py exceeded the {settings.MEEKO_TIMEOUT_SECONDS}s timeout. "
                "This is usually just a slow CPU allocation (common on free-tier hosting), not a "
                "stuck process — try raising DOCKSMART_MEEKO_TIMEOUT, or retry at a quieter time."
            ) from exc
        except FileNotFoundError as exc:
            raise ReceptorPrepError(
                "mk_prepare_receptor.py not found on PATH (part of the 'meeko' package)."
            ) from exc

    if not output_pdbqt.exists():
        raise ReceptorPrepError("Receptor PDBQT was not produced — check Meeko output.")

    if dropped_residues:
        logger.warning(
            "Receptor PDBQT produced with %d residue(s) dropped due to ambiguous "
            "inter-residue bond geometry: %s. If any of these line the binding "
            "site, treat this result with caution.",
            len(dropped_residues), ", ".join(dropped_residues),
        )

    logger.info("Receptor converted to PDBQT: %s", output_pdbqt)
    return output_pdbqt, dropped_residues
