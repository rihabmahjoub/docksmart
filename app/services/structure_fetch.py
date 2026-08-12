"""
Structure retrieval.

Gives the user three entry points for a receptor, matching how DockSmart
should sit relative to the rest of the platform:

  * a PDB identifier            -> fetched live from RCSB
  * a UniProt accession         -> fetched from the AlphaFold DB (handoff
                                    point for the platform's AlphaFold
                                    structure-analyser tool)
  * a local upload              -> user's own file (e.g. output of the
                                    protein-preparation tool)
"""
from __future__ import annotations

import logging
from pathlib import Path

import requests

from app.config import settings

logger = logging.getLogger(__name__)


class StructureFetchError(RuntimeError):
    pass


def fetch_pdb_by_id(pdb_id: str, dest_dir: Path) -> Path:
    pdb_id = pdb_id.strip().upper()
    if not (len(pdb_id) == 4 and pdb_id[0].isdigit()):
        raise StructureFetchError(f"'{pdb_id}' does not look like a valid 4-character PDB ID.")

    url = settings.RCSB_DOWNLOAD_URL.format(pdb_id=pdb_id)
    dest = dest_dir / f"{pdb_id}.pdb"
    _download(url, dest, not_found_hint=f"PDB entry {pdb_id} was not found on RCSB.")
    return dest


def fetch_alphafold_model(uniprot_id: str, dest_dir: Path) -> Path:
    uniprot_id = uniprot_id.strip().upper()
    url = settings.ALPHAFOLD_MODEL_URL.format(uniprot_id=uniprot_id)
    dest = dest_dir / f"AF-{uniprot_id}.pdb"
    _download(
        url,
        dest,
        not_found_hint=(
            f"No AlphaFold model found for UniProt '{uniprot_id}'. "
            "If you already have a predicted structure, upload it directly instead."
        ),
    )
    return dest


def _download(url: str, dest: Path, not_found_hint: str) -> None:
    try:
        resp = requests.get(url, timeout=30)
    except requests.RequestException as exc:
        raise StructureFetchError(f"Network error while fetching {url}: {exc}") from exc

    if resp.status_code == 404:
        raise StructureFetchError(not_found_hint)
    try:
        resp.raise_for_status()
    except requests.HTTPError as exc:
        # Any non-404 HTTP error (403/5xx/etc.) must also become a
        # StructureFetchError, not propagate as a raw requests exception —
        # otherwise the API router's `except StructureFetchError` clause
        # never catches it and the client gets an opaque 500.
        raise StructureFetchError(
            f"Fetching {url} failed with HTTP {resp.status_code}."
        ) from exc

    dest.write_bytes(resp.content)
    logger.info("Downloaded structure to %s", dest)


def save_upload(file_bytes: bytes, filename: str, dest_dir: Path) -> Path:
    dest = dest_dir / Path(filename).name
    dest.write_bytes(file_bytes)
    return dest
