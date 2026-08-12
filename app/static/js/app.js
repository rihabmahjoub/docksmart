// ---- Shared helpers ---------------------------------------------------
function qs(sel) { return document.querySelector(sel); }
function qsa(sel) { return document.querySelectorAll(sel); }

function jobFileUrl(absolutePath, download) {
  // Result paths from the backend are absolute server-side paths (e.g.
  // ".../data/jobs/<job_id>/poses_pdb/pose_1.pdbqt") — pose files live in
  // a subdirectory, so stripping to just the basename (as an earlier
  // version of this function did) breaks the lookup. Preserve everything
  // after "/<job_id>/" instead, which the backend's job_file endpoint
  // resolves correctly (including subpaths) against that job's directory.
  const marker = `/${window.JOB_ID}/`;
  const idx = absolutePath.indexOf(marker);
  const relPath = idx >= 0 ? absolutePath.slice(idx + marker.length) : absolutePath.split("/").pop();
  const q = download ? "&download=true" : "";
  return `/api/jobs/${window.JOB_ID}/file?path=${encodeURIComponent(relPath)}${q}`;
}

// ---- index.html: tabs ---------------------------------------------------
qsa(".tab-btn").forEach(btn => {
  btn.addEventListener("click", () => {
    const group = btn.closest(".tabs");
    group.querySelectorAll(".tab-btn").forEach(b => b.classList.remove("active"));
    btn.classList.add("active");
    group.parentElement.querySelectorAll(".tab-panel").forEach(p => p.classList.remove("active"));
    document.getElementById(`tab-${btn.dataset.tab}`).classList.add("active");
  });
});

const advToggle = qs("#advanced-toggle");
if (advToggle) {
  advToggle.addEventListener("click", () => {
    const panel = qs("#advanced-panel");
    panel.style.display = panel.style.display === "block" ? "none" : "block";
  });
}

const pocketMode = qs("#pocket_mode");
if (pocketMode) {
  pocketMode.addEventListener("change", () => {
    qs("#manual-box").style.display = pocketMode.value === "manual" ? "block" : "none";
  });
}

// ---- index.html: form submission ---------------------------------------
const dockForm = qs("#dock-form");
if (dockForm) {
  dockForm.addEventListener("submit", async (e) => {
    e.preventDefault();
    const submitBtn = dockForm.querySelector("button[type=submit]");
    submitBtn.disabled = true;
    submitBtn.textContent = "Submitting…";

    const formData = new FormData();

    // Receptor tab: PDB ID or file upload (AlphaFold/UniProt option removed
    // after repeated fetch failures in testing — RCSB ID or local upload
    // cover the supported paths).
    const receptorTabId = document.querySelector(".card:nth-of-type(1) .tab-btn.active").dataset.tab;
    if (receptorTabId === "pdb") {
      formData.append("receptor_pdb_id", qs("#receptor_pdb_id").value.trim());
    } else {
      const f = qs("#receptor_file").files[0];
      if (f) formData.append("receptor_file", f);
    }

    const ligandTabId = document.querySelector(".card:nth-of-type(2) .tab-btn.active").dataset.tab;
    if (ligandTabId === "smiles") {
      formData.append("ligand_smiles", qs("#ligand_smiles").value.trim());
    } else {
      const f = qs("#ligand_file").files[0];
      if (f) formData.append("ligand_file", f);
    }

    formData.append("pocket_mode", pocketMode.value);
    if (pocketMode.value === "manual") {
      ["center_x", "center_y", "center_z", "size_x", "size_y", "size_z"].forEach(id => {
        formData.append(id, qs(`#${id}`).value);
      });
    }

    dockForm.querySelectorAll('[name="exhaustiveness"], [name="num_modes"]').forEach(el => {
      formData.append(el.name, el.value);
    });
    formData.append("remove_waters", dockForm.querySelector('[name="remove_waters"]').checked);
    formData.append("rebuild_missing_loops", dockForm.querySelector('[name="rebuild_missing_loops"]').checked);

    try {
      const resp = await fetch("/api/jobs", { method: "POST", body: formData });
      if (!resp.ok) {
        const err = await resp.json();
        alert("Could not start job: " + (err.detail || resp.statusText));
        submitBtn.disabled = false;
        submitBtn.textContent = "Run docking";
        return;
      }
      const data = await resp.json();
      window.location.href = `/jobs/${data.job_id}`;
    } catch (err) {
      alert("Network error: " + err.message);
      submitBtn.disabled = false;
      submitBtn.textContent = "Run docking";
    }
  });
}

// ---- job_status.html: polling + results rendering ------------------------
if (typeof window.JOB_ID !== "undefined") {
  const STAGE_ORDER = [
    "queued", "structure_qc", "pocket_detection", "receptor_pdbqt",
    "ligand_pdbqt", "docking", "postprocessing", "rendering", "done",
  ];

  let lastResult = null;
  let viewer = null;

  async function poll() {
    const resp = await fetch(`/api/jobs/${window.JOB_ID}/status`);
    const data = await resp.json();
    updateStageList(data.stage, data.status);

    if (data.status === "failed") {
      qs("#error-container").innerHTML =
        `<div class="error-box"><strong>Failed at stage: ${data.stage}</strong>\n${data.error_message || ""}</div>`;
      return;
    }
    if (data.status === "succeeded") {
      await loadResults();
      return;
    }
    setTimeout(poll, 2500);
  }

  function updateStageList(currentStage, status) {
    const idx = STAGE_ORDER.indexOf(currentStage);
    qsa("#stage-list li").forEach(li => {
      const liIdx = STAGE_ORDER.indexOf(li.dataset.stage);
      li.classList.remove("done", "active");
      if (status === "succeeded" || liIdx < idx) li.classList.add("done");
      else if (liIdx === idx) li.classList.add("active");
    });
  }

  async function loadResults() {
    const resp = await fetch(`/api/jobs/${window.JOB_ID}/result`);
    const result = await resp.json();
    lastResult = result;

    if (result.pocket_info) {
      qs("#pocket-card").style.display = "block";
      qs("#pocket-info").innerHTML = renderPocketInfo(result.pocket_info);
    }

    if (result.receptor_dropped_residues && result.receptor_dropped_residues.length > 0) {
      qs("#receptor-warning").innerHTML = `
        <div class="warning-box">
          <strong>Note:</strong> ${result.receptor_dropped_residues.length} receptor residue(s)
          (${result.receptor_dropped_residues.join(", ")}) were automatically excluded during
          receptor preparation because their inter-residue bond geometry was ambiguous
          (likely a crystal-packing artifact rather than real chemistry). If any of these sit
          in the binding site you used, treat this result with caution and inspect the input
          structure around those residues.
        </div>`;
    }

    renderPoseTable(result.poses);
    qs("#results-card").style.display = "block";

    if (result.receptor_pdb && result.poses.length > 0) {
      qs("#viewer-card").style.display = "block";
      setupViewer(result);
    }
  }

  function renderPoseTable(poses) {
    const tbody = qs("#pose-tbody");
    tbody.innerHTML = "";
    poses.forEach((pose, i) => {
      const tr = document.createElement("tr");
      if (i === 0) tr.classList.add("best");
      const fp = pose.interaction_summary || {};
      const downloadCell = pose.pose_pdbqt
        ? `<a href="${jobFileUrl(pose.pose_pdbqt, true)}" download>⬇ PDBQT</a>`
        : "—";
      tr.innerHTML = `
        <td>${pose.mode}</td>
        <td>${pose.affinity_kcal_mol.toFixed(2)}</td>
        <td>${pose.rmsd_ub === null || pose.rmsd_ub === undefined ? "—" : pose.rmsd_ub.toFixed(2)}</td>
        <td>${fp.hbond_contacts ?? "—"}</td>
        <td>${fp.hydrophobic_contacts ?? "—"}</td>
        <td class="pose-download">${downloadCell}</td>`;
      tbody.appendChild(tr);
    });
  }

  function renderPocketInfo(info) {
    if (info.mode === "manual") {
      return `<div class="hint">Manual binding-site box — center and size as specified in the job form.</div>`;
    }
    if (!info.selected_pocket) return "";
    const p = info.selected_pocket;
    const drug = p.druggability_score;
    const score = p.score;  // fpocket's general pocket "openness/quality" score — NOT a separate metric named "pocket score"
    return `
      <div class="score-display">
        <div class="score-block">
          <div class="score-value">${drug !== undefined ? drug.toFixed(2) : "n/a"}</div>
          <div class="score-label">Druggability score</div>
          <div class="score-def">fpocket's machine-learned estimate (0–1) of how likely this cavity is to bind a small molecule with reasonable affinity, based on pocket geometry and physicochemical properties. It is a plausibility heuristic, not a thermodynamic binding prediction.</div>
        </div>
        <div class="score-block">
          <div class="score-value">${score !== undefined ? score.toFixed(2) : "n/a"}</div>
          <div class="score-label">Pocket score</div>
          <div class="score-def">fpocket's general cavity-quality score, combining pocket shape and size. Higher generally means a larger, more well-defined cavity — it is a geometric descriptor, not a binding-affinity estimate.</div>
        </div>
      </div>
      <div class="hint" style="margin-top:10px;">Selected as the top-ranked cavity out of ${info.all_pockets_ranked ? info.all_pockets_ranked.length : "?"} detected by fpocket (pocket #${p.pocket_id}).</div>`;
  }

  // ---- CSV export (client-side, no backend endpoint needed) --------------
  qs("#download-csv-btn")?.addEventListener("click", () => {
    if (!lastResult) return;
    const rows = [["Mode", "Affinity (kcal/mol)", "RMSD u.b.", "H-bond contacts", "Hydrophobic contacts"]];
    lastResult.poses.forEach(p => {
      const fp = p.interaction_summary || {};
      rows.push([
        p.mode, p.affinity_kcal_mol.toFixed(2),
        p.rmsd_ub === null || p.rmsd_ub === undefined ? "" : p.rmsd_ub.toFixed(2),
        fp.hbond_contacts ?? "", fp.hydrophobic_contacts ?? "",
      ]);
    });
    const csv = rows.map(r => r.join(",")).join("\n");
    const blob = new Blob([csv], { type: "text/csv" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = `docksmart_${window.JOB_ID}_top_poses.csv`;
    a.click();
    URL.revokeObjectURL(a.href);
  });

  qs("#download-pdbqt-btn")?.addEventListener("click", () => {
    if (!lastResult || !lastResult.poses[0] || !lastResult.poses[0].pose_pdbqt) return;
    window.location.href = jobFileUrl(lastResult.poses[0].pose_pdbqt, true);
  });

  // ---- 3Dmol.js viewer ----------------------------------------------------

  // Wait for the next paint before touching layout-dependent APIs. The
  // #viewer-card is switched from display:none to display:block just
  // before setupViewer() is called; reading the container's size or
  // creating the 3Dmol canvas in that same synchronous tick can pick up
  // stale (zero / display:none-era) dimensions, which is what caused the
  // viewer box to render mis-sized and mis-positioned near the top-left of
  // the page instead of inside its card, and made zoomTo() frame the
  // camera incorrectly (so one of the two loaded models — usually the
  // ligand — fell outside the visible frame). Two rAFs reliably span a
  // full layout+paint cycle across browsers; one is not always enough.
  function nextPaint() {
    return new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve)));
  }

  async function setupViewer(result) {
    const poseSelect = qs("#pose-select");
    poseSelect.innerHTML = "";
    result.poses.forEach((p, i) => {
      const opt = document.createElement("option");
      opt.value = i;
      opt.textContent = `Pose ${p.mode} (${p.affinity_kcal_mol.toFixed(2)} kcal/mol)${i === 0 ? " — best" : ""}`;
      poseSelect.appendChild(opt);
    });

    await nextPaint();

    const container = qs("#viewer-container");

    function showViewerError(message) {
      console.error("[DockSmart viewer]", message);
      container.innerHTML = `<div style="padding:16px; color:#7a241d; font-size:0.85rem;">
        Could not display the 3D structure: ${message}. Open the browser console for details —
        the pose data itself is unaffected; you can still download each pose's PDBQT from the
        table below.
      </div>`;
    }

    // A structure fetch returning an HTTP error, or a job file path that no
    // longer resolves, previously failed *silently*: addModel() on an empty
    // or non-PDB response (e.g. a JSON 404 body) parses to zero atoms with
    // no exception, so nothing rendered and nothing was logged either. This
    // wrapper makes that failure visible instead — check the browser
    // console for the exact reason if the structure still doesn't appear.
    async function fetchStructureText(absolutePath, label) {
      const url = jobFileUrl(absolutePath);
      const resp = await fetch(url);
      if (!resp.ok) {
        throw new Error(`${label} fetch failed (HTTP ${resp.status} for ${url})`);
      }
      const text = await resp.text();
      if (!/^(ATOM|HETATM|MODEL|REMARK|CRYST1|HEADER)/m.test(text)) {
        throw new Error(`${label} response did not look like a PDB file (fetched from ${url})`);
      }
      return text;
    }

    let viewerReady = false;
    try {
      viewer = $3Dmol.createViewer("viewer-container", { backgroundColor: "white" });

      const receptorText = await fetchStructureText(result.receptor_pdb, "Receptor structure");
      viewer.addModel(receptorText, "pdb");  // model 0: receptor, stays loaded across pose switches
      viewerReady = true;
    } catch (err) {
      showViewerError(err.message || String(err));
      return;
    }

    async function loadPose(idx) {
      try {
        // Remove any previously-loaded ligand model (model 1) before adding
        // the newly selected one, rather than stacking them.
        const existing = viewer.getModel(1);
        if (existing) viewer.removeModel(existing);
        const pose = result.poses[idx];
        const poseText = await fetchStructureText(pose.pose_pdb, "Ligand pose");
        viewer.addModel(poseText, "pdb");
        applyStyles();
        viewer.zoomTo({ model: 1 });
        viewer.resize();
        viewer.render();
      } catch (err) {
        showViewerError(err.message || String(err));
      }
    }

    function applyStyles() {
      const bg = qs("#viz-bg").value;
      viewer.setBackgroundColor(bg === "dark" ? "#0a0e14" : "white");

      const style = qs("#viz-style").value;
      const colorMode = qs("#viz-color").value;

      let proteinStyle = {};
      const colorSpec = colorMode === "chain" ? { colorscheme: "chain" }
        : colorMode === "element" ? { colorscheme: "default" }
        : colorMode === "spectrum" ? { color: "spectrum" }
        : {};  // "interaction" mode colors are applied separately below

      if (style === "cartoon") proteinStyle = { cartoon: { ...colorSpec } };
      else if (style === "surface") proteinStyle = { cartoon: { ...colorSpec } };
      else if (style === "stick") proteinStyle = { stick: { ...colorSpec, radius: 0.15 } };
      else if (style === "line") proteinStyle = { line: { ...colorSpec } };

      viewer.setStyle({ model: 0 }, colorMode === "interaction" ? { cartoon: { color: "lightgrey" } } : proteinStyle);

      viewer.removeAllSurfaces();
      if (style === "surface") {
        // Opacity kept well under 1 (and a light grey rather than solid
        // white) specifically so the ligand sticks inside the pocket
        // remain visible through the surface rather than being visually
        // buried by it — a fully opaque/near-opaque surface was the other
        // reason the ligand looked "missing" in this mode specifically.
        viewer.addSurface(
          $3Dmol.SurfaceType.VDW,
          { opacity: 0.55, color: "#dfe3e6" },
          { model: 0 }
        );
      }

      if (colorMode === "interaction") {
        // Highlight residues within 4.5 Å of the ligand model (model 1) —
        // a simple, always-available visual proxy for "the binding site",
        // independent of whether the richer ProLIF fingerprint ran.
        viewer.setStyle(
          { model: 0, within: { distance: 4.5, sel: { model: 1 } } },
          { stick: { color: "orange", radius: 0.18 }, cartoon: { color: "lightgrey" } }
        );
      }

      // Ligand: always shown as colored sticks regardless of protein style,
      // so it never silently disappears depending on the chosen representation.
      viewer.setStyle({ model: 1 }, { stick: { colorscheme: "greenCarbon", radius: 0.25 } });
    }

    poseSelect.addEventListener("change", () => loadPose(parseInt(poseSelect.value, 10)));
    ["viz-bg", "viz-style", "viz-color"].forEach(id => {
      qs(`#${id}`).addEventListener("change", () => { applyStyles(); viewer.render(); });
    });

    // Keep the canvas matched to its container across browser resizes /
    // devtools panel toggles / responsive breakpoint changes, not just at
    // initial load.
    window.addEventListener("resize", () => { if (viewer) { viewer.resize(); viewer.render(); } });

    if (viewerReady) await loadPose(0);
  }

  poll();
}
