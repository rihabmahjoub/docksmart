"""
Post-docking analysis.

 * split_poses_to_pdb  — pull individual MODEL blocks out of Vina's
   multi-pose output PDBQT and convert each to a plain PDB for
   visualization / downstream tools.

 * interaction_fingerprint — a lightweight, dependency-free distance-based
   contact analysis (H-bond donor/acceptor proximity, hydrophobic
   contacts) used to rerank/annotate poses beyond the raw Vina score.
   This is deliberately implemented without ProLIF/PLIP so the core
   pipeline has no hard dependency on them; if `prolif` is installed,
   `interaction_fingerprint_prolif` gives a richer, published-method
   fingerprint and should be preferred when available.

IMPORTANT — why this module shells out to the `obabel` CLI instead of
`from openbabel import pybel`: confirmed by direct testing that importing
`openbabel.pybel` in the same Python process as the `vina` package causes
an uncatchable C++ exception (`swig::stop_iteration`) that aborts the
whole process — a real ABI conflict between their bundled SWIG bindings,
not a hypothetical one. Since `docking_engine.py` imports `vina` earlier
in the pipeline, this module MUST keep using `subprocess` for OpenBabel,
never a direct Python import of `openbabel`/`pybel`, or the pipeline will
crash intermittently depending on import order.
"""
from __future__ import annotations

import logging
import subprocess
from pathlib import Path

logger = logging.getLogger(__name__)

# Very small default atom-type classification for the distance-based fallback.
HBOND_ELEMENTS = {"N", "O"}
HYDROPHOBIC_ELEMENTS = {"C"}
HBOND_CUTOFF_A = 3.5
HYDROPHOBIC_CUTOFF_A = 4.5


def split_poses_to_pdb(poses_pdbqt: Path, out_dir: Path) -> list[Path]:
    """Split a multi-MODEL Vina PDBQT into per-pose PDB files.

    Tries OpenBabel first (best general-purpose conversion), but ALWAYS
    falls back to a pure-Python PDBQT->PDB converter rather than silently
    handing back the raw .pdbqt file relabeled as .pdb — that fallback
    (the previous behavior) is confirmed to have caused "ligand doesn't
    show up in the viewer" reports: PDBQT looks superficially like PDB
    (same coordinate columns) but has AutoDock-specific keywords (ROOT,
    BRANCH, TORSDOF) and an atom-type column instead of a standard element
    column, which a strict PDB parser can choke on silently. The
    pure-Python converter below only ever emits well-formed PDB.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    models = _split_pdbqt_models(poses_pdbqt)
    pdb_paths = []
    for i, model_lines in enumerate(models, start=1):
        raw_pdbqt = out_dir / f"pose_{i}.pdbqt"
        raw_pdbqt.write_text("\n".join(model_lines))
        pdb_path = out_dir / f"pose_{i}.pdb"

        converted = False
        try:
            subprocess.run(
                ["obabel", str(raw_pdbqt), "-O", str(pdb_path)],
                check=True, capture_output=True, text=True, timeout=30,
            )
            # obabel can exit 0 but still write an empty/near-empty file on
            # malformed input — verify there's actually atom content before
            # trusting it, rather than just trusting the exit code.
            if pdb_path.exists() and _count_pdb_atoms(pdb_path) == _count_pdbqt_atoms(model_lines):
                converted = True
        except (FileNotFoundError, subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
            logger.warning("OpenBabel conversion failed for pose %d, using fallback converter: %s", i, exc)

        if not converted:
            pdb_path.write_text("\n".join(_pdbqt_lines_to_pdb(model_lines)))

        pdb_paths.append(pdb_path)
    return pdb_paths


def _count_pdb_atoms(pdb_path: Path) -> int:
    return sum(1 for line in pdb_path.read_text(errors="ignore").splitlines() if line.startswith(("ATOM", "HETATM")))


def _count_pdbqt_atoms(model_lines: list[str]) -> int:
    return sum(1 for line in model_lines if line.startswith(("ATOM", "HETATM")))


# AutoDock atom-type codes (PDBQT's equivalent of a PDB element column) ->
# standard element symbol, covering the types Meeko/AutoDock actually emit.
_AUTODOCK_TYPE_TO_ELEMENT = {
    "A": "C", "C": "C",              # aromatic / aliphatic carbon
    "N": "N", "NA": "N", "NS": "N",
    "OA": "O", "OS": "O", "O": "O",
    "SA": "S", "S": "S",
    "H": "H", "HD": "H", "HS": "H",
    "P": "P", "F": "F",
    "CL": "Cl", "Cl": "Cl", "BR": "Br", "Br": "Br", "I": "I",
    "MG": "Mg", "MN": "Mn", "ZN": "Zn", "CA": "Ca", "FE": "Fe",
    "NA+": "Na", "K": "K", "CU": "Cu",
}


def _pdbqt_lines_to_pdb(model_lines: list[str]) -> list[str]:
    """Convert one MODEL block's ATOM/HETATM lines from PDBQT to well-formed
    PDB, dropping AutoDock-specific keyword lines (ROOT/ENDROOT/BRANCH/
    ENDBRANCH/TORSDOF/REMARK) and mapping the AutoDock atom-type column to a
    standard element symbol so downstream tools (3Dmol.js, PyMOL, etc.) can
    parse bonds/elements correctly."""
    out = []
    serial = 0
    for line in model_lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        serial += 1
        # PDBQT shares PDB's column layout through the temperature factor
        # (cols 1-66); columns after that hold partial charge + AutoDock
        # type instead of PDB's element symbol.
        name = line[12:16]
        resname = line[17:20]
        chain = line[21] if len(line) > 21 else "A"
        resseq = line[22:26]
        x, y, z = line[30:38], line[38:46], line[46:54]
        occ = line[54:60] if len(line) >= 60 else "  1.00"
        temp = line[60:66] if len(line) >= 66 else "  0.00"
        adtype = line[77:].strip().split()[-1] if line[77:].strip() else name.strip()[:1]
        element = _AUTODOCK_TYPE_TO_ELEMENT.get(adtype, adtype[:1].upper())
        out.append(
            f"ATOM  {serial:>5} {name:<4} {resname:>3} {chain}{resseq:>4}    "
            f"{x:>8}{y:>8}{z:>8}{occ:>6}{temp:>6}          {element:>2}"
        )
    out.append("END")
    return out


def _split_pdbqt_models(pdbqt_file: Path) -> list[list[str]]:
    models: list[list[str]] = []
    current: list[str] = []
    for line in pdbqt_file.read_text().splitlines():
        if line.startswith("MODEL"):
            current = []
        elif line.startswith("ENDMDL"):
            models.append(current)
        else:
            current.append(line)
    return models


def interaction_fingerprint(receptor_pdb: Path, pose_pdb: Path) -> dict:
    """Cheap distance-based contact summary: counts of H-bond-capable and
    hydrophobic contacts within cutoff between receptor and pose atoms."""
    receptor_atoms = _read_atoms(receptor_pdb)
    ligand_atoms = _read_atoms(pose_pdb)

    hbond_contacts = 0
    hydrophobic_contacts = 0
    contact_residues: set[str] = set()

    for lx, ly, lz, lelem, *_ in ligand_atoms:
        for rx, ry, rz, relem, rresname, rresi, rchain in receptor_atoms:
            d2 = (lx - rx) ** 2 + (ly - ry) ** 2 + (lz - rz) ** 2
            if lelem in HBOND_ELEMENTS and relem in HBOND_ELEMENTS and d2 <= HBOND_CUTOFF_A ** 2:
                hbond_contacts += 1
                contact_residues.add(f"{rchain}:{rresname}{rresi}")
            elif (
                lelem in HYDROPHOBIC_ELEMENTS
                and relem in HYDROPHOBIC_ELEMENTS
                and d2 <= HYDROPHOBIC_CUTOFF_A ** 2
            ):
                hydrophobic_contacts += 1
                contact_residues.add(f"{rchain}:{rresname}{rresi}")

    return {
        "hbond_contacts": hbond_contacts,
        "hydrophobic_contacts": hydrophobic_contacts,
        "contact_residues": sorted(contact_residues),
        "method": "distance_cutoff_fallback",
    }


def interaction_fingerprint_prolif(receptor_pdb: Path, pose_pdb: Path) -> dict | None:
    """Preferred, richer fingerprint using ProLIF, if installed."""
    try:
        import MDAnalysis as mda
        import prolif as plf
    except ImportError:
        return None

    try:
        protein = mda.Universe(str(receptor_pdb))
        ligand = mda.Universe(str(pose_pdb))
        fp = plf.Fingerprint()
        fp.run(ligand.trajectory, ligand, protein)
        df = fp.to_dataframe()
        interactions = sorted({col[2] for col in df.columns})
        residues = sorted({str(col[1]) for col in df.columns})
        return {"interaction_types": interactions, "contact_residues": residues, "method": "prolif"}
    except Exception as exc:  # ProLIF/MDAnalysis edge cases vary a lot by input
        logger.warning("ProLIF fingerprinting failed, falling back: %s", exc)
        return None


def _read_atoms(pdb_file: Path):
    """Return (x, y, z, element, resname, resi, chain) tuples from a PDB file."""
    atoms = []
    for line in pdb_file.read_text(errors="ignore").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue
        element = (line[76:78].strip() or line[12:14].strip())[:1].upper()
        resname = line[17:20].strip()
        resi = line[22:26].strip()
        chain = line[21].strip() or "A"
        atoms.append((x, y, z, element, resname, resi, chain))
    return atoms
