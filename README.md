# DockSmart

An automated web tool for ligand–protein interaction studies. Automates
receptor/ligand preparation, binding-site detection, AutoDock Vina docking,
pose analysis, and visualization behind a single web interface.

Part of a larger structural bioinformatics platform: IntelligentLigPrep
(ligand curation/ADMET), a protein-preparation tool, DockSmart (this repo),
a molecular-dynamics platform, and an AlphaFold structure analyser.
DockSmart deliberately does **not** reimplement ligand ADMET filtering or
full protein preparation — it accepts prepared inputs from those tools and
focuses on the docking-specific stages (binding-site detection, PDBQT
conversion, Vina execution, pose interpretation).

## Pipeline

```
Receptor (PDB ID / UniProt→AlphaFold / upload)
Ligand   (SMILES / upload / IntelligentLigPrep handoff)
        │
        ▼
1. Structure QC — PDBFixer (missing atoms, protonation, water removal)
2. Binding-site detection — fpocket (ranked, druggability-scored) |
                              co-crystal ligand | manual box
3. PDBQT preparation — Meeko (receptor + ligand)
4. Docking — AutoDock Vina (Python bindings, CLI fallback)
5. Post-processing — pose splitting, interaction fingerprint,
                      redocking RMSD-to-reference (validation mode)
6. Visualization — headless PyMOL (static figure) + NGL.js (interactive)
```

## Running locally

```bash
pip install -r requirements.txt
# fpocket has no pip wheel — compile it from source (fast, small codebase,
# only needs gcc/g++/make — see Dockerfile for the exact recipe) or install
# it via your system package manager, then make sure the `fpocket` binary
# is on PATH.
uvicorn app.main:app --reload
```

Visit `http://localhost:8000`.

## Verified dependency stack (2026-08-10)

Every package version in `requirements.txt`, and every module the app
actually imports, was installed together in a clean venv and exercised
through a real end-to-end docking run (real PDB structure, real ligand,
real Vina docking, real PyMOL-rendered figure, through the actual HTTP
API) — not just resolved on paper. That run is what the current
`Dockerfile` is built around. Notable things that were NOT obvious going
in, and would have caused a second or third failed Render build if left
unpinned:

- `openbabel` on PyPI is a source distribution requiring SWIG + a
  pre-installed system OpenBabel — the actual installable wheel is the
  separate `openbabel-wheel` package.
- `meeko` has two runtime dependencies (`scipy`, `gemmi`) it does not
  declare in its own metadata — pip won't install them automatically.
- `pymol-open-source-whl` (not the ~185MB conda-forge `pymol-open-source`
  package) is a self-contained ~15MB wheel with no Qt/GStreamer/X11
  dependency chain, but hard-pins `numpy==1.26.4` — pin numpy explicitly
  in your own environment or a later `pip install` of something else can
  silently upgrade it and break PyMOL's import.
- **`import vina` followed by `from openbabel import pybel` in the same
  process crashes with an uncatchable C++ exception** (confirmed, not
  theoretical — a SWIG/ABI conflict between their bundled bindings).
  `app/services/postprocessing.py` therefore calls the `obabel` CLI via
  `subprocess`, never a direct Python import of `openbabel` — this is
  load-bearing, not a style choice; see the warning comment in that file.
- `cmd.util.cbag()` (PyMOL's standard "color by atom, green carbons"
  helper) raises `RuntimeError: Missing PyMOL instance` under
  `pymol-open-source-whl`'s multi-instance (`pymol2.PyMOL()`) API — it
  reaches for a different global instance internally. Worked around with
  two plain `cmd.color()` calls in `app/services/visualization.py`.
- `mk_prepare_receptor.py`'s `-v`/`--write_vina_box` flag requires box
  center/size to already be known (it writes a Vina config file) — don't
  pass it if, like this pipeline, you compute the box separately via
  fpocket and pass it directly to `run_vina()` instead.

## Deployment on Render (free tier)

The `Dockerfile` builds `fpocket` from source in a slim `debian:bookworm-
slim` stage (compiles in seconds, ~1.9MB binary, no BLAS/LAPACK actually
required despite some docs implying otherwise) and copies just the binary
into a `python:3.11-slim` final image — no conda/mamba anywhere in the
build. This was a deliberate rework after the first version (based on
`condaforge/mambaforge` + conda-installed PyMOL) pulled in ~185MB of
unnecessary Qt/GStreamer/X11 libraries, a bad trade against Render free
tier's 512MB RAM / 0.1 CPU ceiling.

Even with the leaner image, treat free tier as genuinely tight, not
comfortable:
- Free tier spins down after 15 min of inactivity; expect a 30–60s cold
  start on the next request.
- Keep `exhaustiveness`/`num_modes` low for your first live tests
  (e.g. 4 and 3) — confirm the pipeline completes at all before pushing
  compute up.
- If you hit OOM or a build/runtime timeout even after this rework,
  that's the signal to move to Render's paid Starter tier (512MB/0.5 CPU
  — same RAM, 5x the CPU) rather than trying to trim further; there isn't
  much fat left to cut in this image.
- Render's free-tier availability itself has been reported inconsistently
  across sources during 2026 — check Render's current pricing page before
  each deploy attempt rather than assuming it's still there.

## Scientific correctness notes (important for the manuscript's Methods section)

These are documented in code comments at the relevant module, summarized
here so they don't get lost:

1. **RMSD is not RMSD-to-native.** Vina's `rmsd_lb`/`rmsd_ub` (whether from
   the CLI table or our own fallback calculation on the Python-bindings
   path) describe how much each reported pose diverges from the *top
   pose of that same docking run* — never from a crystallographic
   reference. For redocking/self-docking benchmarking against a
   co-crystallized ligand (e.g. your planned SeamDock comparison), use
   `docking_engine.rmsd_to_reference()`, which computes true pose-vs-
   reference RMSD via RDKit's symmetry-aware `GetBestRMS`. Report exactly
   which of these two RMSD concepts appears in each table/figure —
   conflating them is a common and reviewer-visible error in docking
   papers.

2. **The Python Vina bindings don't expose the CLI's RMSD table at all.**
   `Vina.energies()` returns `[total, inter, intra, torsions,
   best_intra]` energy terms, not affinity+RMSD. `docking_engine.py`
   computes its own naive (non-symmetry-corrected) Cartesian RMSD in that
   code path and labels it explicitly (`rmsd_note`) so it's never confused
   with Vina CLI's symmetry-aware value.

3. **Missing-loop rebuilding is off by default.** PDBFixer can both (a)
   complete missing atoms within already-resolved residues (safe, always
   on) and (b) build in whole unresolved loops the crystal structure never
   observed (risky for rigid docking — can bias the grid box or clash
   with the ligand with no experimental support). DockSmart defaults to
   (a) only; loop rebuilding is an explicit opt-in
   (`rebuild_missing_loops`), and should only be used when the disordered
   region actually lines the binding site.

4. **Hydrogen addition is a fixed-pH heuristic, not true pKa prediction.**
   PDBFixer's `addMissingHydrogens(ph)` is not equivalent to PROPKA/H++
   -style titration. Fine as a default for a physiological pH; flag this
   as a limitation if a target has titratable catalytic/binding-site
   residues where protonation state materially affects docking.

5. **The built-in interaction fingerprint is a distance-cutoff heuristic,
   not a validated H-bond/π-stacking detector.** It classifies contacts
   purely by element and a distance cutoff (3.5 Å for N/O···N/O, 4.5 Å for
   C···C) with no donor–H···acceptor angle check and no aromatic-geometry
   check for π-stacking. `postprocessing.interaction_fingerprint_prolif()`
   (used automatically when ProLIF/MDAnalysis are installed) gives a
   published, geometry-aware fingerprint and should be the one reported in
   the manuscript if available; the built-in version exists purely so the
   pipeline still produces *some* interaction summary when those heavier
   dependencies aren't installed (e.g. on a constrained free host).

6. **fpocket ranking uses its own druggability score**, a machine-learned
   estimate from pocket geometry/physicochemical descriptors — it is a
   plausibility heuristic for "is this cavity drug-like," not a
   thermodynamic prediction of binding, and should be described as such.

7. **AlphaFold models carry no bound ligand and no crystallographic
   resolution.** Per-residue pLDDT (stored in the B-factor column) should
   be checked near the intended binding site before trusting a pocket
   detected there — DockSmart surfaces this as a warning in the
   AlphaFold-preview endpoint but does not currently gate docking on a
   pLDDT threshold; consider adding that gate before claiming AlphaFold
   support as validated in the paper.

8. **Protein/DNA-as-"ligand" docking is scientifically the least mature
   part of the pipeline as scaffolded.** Vina is designed and validated
   for small-molecule ligand docking against a rigid receptor; using it
   for protein–protein or protein–DNA docking is unconventional and
   under-validated compared with purpose-built tools (e.g. HDOCK,
   ClusPro, AutoDock CrankPep). If this claim goes in the abstract, budget
   real validation work (a benchmark set with known complexes) before
   presenting it as a supported feature rather than an experimental one.

## Reliability mechanisms (added after real free-tier failures)

Each of these was added to fix a specific, reproduced failure — not
preemptively:

- **Hard timeout on PDBFixer** (`app/utils/proc_timeout.py`): confirmed to
  hang indefinitely (30+ min, unrecoverable without a redeploy) on a real
  structure, because unlike the subprocess-based tools it had no wall-clock
  boundary at all. Now runs in a separate process that gets forcibly
  killed at `DOCKSMART_STRUCTURE_QC_TIMEOUT` (default 300s).
- **Receptor size pre-check**: rejects receptors over
  `DOCKSMART_MAX_RECEPTOR_ATOMS` (default 15000) immediately, before
  attempting PDBFixer/fpocket, rather than letting a huge structure run
  for minutes only to time out anyway.
- **Concurrency guard**: confirmed that two simultaneous docking jobs OOM-
  kill a 512MB Render instance (exit 137). `DOCKSMART_MAX_CONCURRENT_JOBS`
  (default 1) is now enforced by a real semaphore in `pipeline.py` — a
  second job waits (reported as the "queued" stage) instead of racing the
  first for RAM.
- **Zero-atom PDBFixer crash**: a receptor file with no ATOM records (only
  HETATM — e.g. a mislabeled export, or a ligand file uploaded as the
  receptor by mistake) used to crash deep inside OpenMM with an opaque
  "Cannot create a Context for a System with no particles" error.
  Reproduced and fixed with an explicit atom-count check and a clear
  message before that point.
- **Meeko/fpocket timeouts** (`DOCKSMART_MEEKO_TIMEOUT`,
  `DOCKSMART_FPOCKET_TIMEOUT`, both default 600s): originally hardcoded at
  120s/180s, confirmed too short even for a moderate receptor on Render's
  0.1 CPU free tier. Now configurable — raise further if a genuinely large
  receptor still times out.
- **Pose file conversion**: the previous PDBQT→PDB fallback silently
  handed back a raw `.pdbqt` file relabeled `.pdb` when OpenBabel failed,
  which produced poses that looked fine server-side but failed to render
  a ligand in the browser viewer. Replaced with a pure-Python converter
  (byte-verified against the PDB column spec) that never produces a
  malformed file.

## Known gaps to close before submission

- Co-crystal-ligand pocket mode (`pocket_detection.box_from_cocrystal_ligand`)
  is implemented but not yet wired to the `/api/jobs` endpoint (needs a
  `resname` parameter exposed in the form) — currently raises a clear
  "not yet wired" error if selected.
- No authentication/rate limiting — fine for a research demo, worth adding
  before wide public release given docking jobs are CPU-expensive.
- GNN rescoring is intentionally deferred to v2, as agreed — the
  `postprocessing.py` interaction fingerprint (and, when installed, the
  ProLIF-based one) is the v1 alternative method for going beyond raw
  Vina score.
- AlphaFold/UniProt receptor fetch was removed from the UI after repeated
  fetch failures in testing (the backend route still exists but is unused
  — worth investigating separately if you want it back, since P00533 does
  have an AlphaFold model and the failure was likely network-related on
  the deployment host, not a logic bug).

## Tests

```bash
python tests/test_pipeline_smoke.py
```

Covers pure parsing/geometry logic (fpocket info-file parsing, Vina stdout
table parsing, grid-box math) that runs without the heavy scientific
dependencies installed — useful for CI. Full pipeline integration testing
requires the real binaries and a test PDB/ligand pair, and should be run
manually against the Render deployment before submission.
