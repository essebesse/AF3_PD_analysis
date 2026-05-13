# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

A Streamlit web app for analyzing protein–protein structure prediction pulldown data.

Supported input formats (auto-detected per folder):

| Format | Detection | Metrics computed |
|---|---|---|
| **AF3 local pipeline** | `seed-N_sample-M/` subdirs + `ranking_scores.csv` | full: iPTM, PTM, ipSAE, interface pLDDT, PAE-binned contacts |
| **AF3 Server flat** | `fold_<name>_model_N.cif` + `fold_<name>_full_data_N.json` directly in folder | full (same as above) |
| **AlphaPulldown / AF2-multimer** | `ranking_debug.json` + `ranked_N.pdb` | partial: combined `iptm+ptm`, spatial Cβ contacts, interface pLDDT from PDB B-factors. No ipSAE/PAE (AF2 output has no per-residue PAE matrix) |
| **Legacy `format: 'af2'` caches** | pre-existing cache files with `iptm: None` + `iptm_ptm` | display-only — fields are normalized at load time |

Computes confidence metrics across **all 5 seed/sample models per prediction**, not just top-ranked.

## Running the App

```bash
# Using the launcher script (activates .venv automatically)
./run.sh

# Or manually
source .venv/bin/activate
streamlit run app.py
```

The app runs on `http://localhost:8501`.

## Dependencies

```bash
pip install -r requirements.txt
```

Key libraries: streamlit, pandas, numpy, matplotlib, gemmi (for CIF parsing), openpyxl, requests (UniProt API).

Python 3.11 venv is at `.venv/`.

## Architecture

**Single-page app with workflow steps** — `app.py` is the entry point and uses `st.radio` to switch between 4 steps. Each step imports a `show_*()` function from `pages/`:

- `app.py` — Main entry: page config, CSS, project folder browser (with recursive AF3 subproject scanner), step routing
- `pages/overview.py` — Step 1 "Load Data": scans AF3 directory, shows prediction summary
- `pages/batch_execution.py` — Step 2 "Analyze": runs analysis via local multiprocessing
- `pages/results.py` — Step 3 "Results": sortable table with filtering and CSV/XLSX export
- `pages/detailed_analysis.py` — Step 4 "Detailed Analysis": PAE plots, interface contacts, hub residues

**Core modules** (`core/`):

- `scanner.py` — `AF3Scanner` class + format detectors (`_is_af3_server_flat`, `_is_alphapulldown_format`). `load_prediction_data()` dispatches to AF3/server-flat/AlphaPulldown loaders. `find_af3_projects_recursive()` walks a directory tree to discover predictions in any of the supported formats.
- `analyzer.py` — Heavy analysis: `analyze_prediction_all_models()` dispatches by format (`_analyze_server_flat`, `_analyze_alphapulldown`). Computes ipSAE (Dunbrack 2025), interface contacts, interface pLDDT. `CIFParser` is **gemmi-backed** — preserves the legacy attribute-dict shape (`atoms[i] = {seq_id, chain_id, atom_name, res_name, x, y, z}`) so callers that walk `parser.atoms` in CIF order (e.g. for `atom_plddts` indexing) keep working.
- `interface_analyzer.py` — Contact classification (hydrophobic/polar/charged/aromatic), hub residue detection. `analyze_interface()` runs PAE + spatial filtering for AF3; `analyze_interface_spatial()` runs spatial-only for AlphaPulldown.
- `pae_plotter.py` — PAE matrix visualization: full matrix, interface zoom, model comparison. Uses green-white colormap.
- `viewer_3d.py` — 3Dmol.js-based 3D structure viewer with PAE interface coloring. **The 3Dmol-min.js bundle is inlined into the HTML via `_load_3dmol_js()` + `@lru_cache`** — no dependency on Streamlit's static serving. Supports both CIF (AF3) and PDB (AlphaPulldown) via `model_format` arg.
- `pymol_script.py` — PyMOL .py script generator for offline interface visualization.
- `slurm_manager.py` — SLURM job submission and status checking for cluster execution.
- `run_analysis_cli.py` — Standalone CLI entry point for SLURM jobs (bypasses Streamlit).
- `utils.py` — `format_score()`, `calculate_confidence_tier()`, `tier_color()`, `split_prediction_name()` (strips AF3 Server `fold_` prefix, splits on `_and_`), `is_valid_uniprot_accession()`, `normalize_cache_records()` (backfills `iptm`/`ptm`/`ranking_score` on legacy af2 caches), `fetch_gene_names_batch()` (UniProt API with pre-filter + per-accession fallback on HTTP 400).

## Expected Data Layout

Three formats, auto-detected per folder by the scanner:

```
# AF3 local pipeline
<project>/AF3/<bait>_and_<target>/
  ranking_scores.csv
  seed-N_sample-M/
    model.cif
    confidences.json
    summary_confidences.json

# AF3 Server flat (downloaded zips from alphafoldserver.com)
<project>/fold_<name>/
  fold_<name>_model_0.cif … fold_<name>_model_4.cif
  fold_<name>_full_data_0.json … fold_<name>_full_data_4.json
  fold_<name>_summary_confidences_0.json … _summary_confidences_4.json

# AlphaPulldown (AF2-multimer)
<project>/pulldown/<bait>_and_<prey>/
  ranking_debug.json
  ranked_0.pdb … ranked_4.pdb
  unrelaxed_model_*_multimer_v3_pred_0.pdb
```

Analysis results cache to `af3_app_all_models_analysis.json` at the predictions-folder level. Per-prediction `af3_app_analysis.json` files are also written inside each pair folder for fast Detailed-Analysis lookup.

## SLURM execution flow

1. User clicks Submit on the SLURM tab → `submit_slurm_jobs()` writes `_slurm_chunk_<i>.txt`, sbatches one job per chunk, sets `st.session_state['slurm_job_ids']`, then **`st.rerun()`** so the folder-state panel above the Submit button picks up the new chunk files.
2. Each SLURM job runs the analyzer over its chunk list and writes `_slurm_results_<i>.json`.
3. The user clicks **🔄 Refresh** on the SLURM tab (no auto-refresh — we used to have one, but `window.parent.location.reload()` was wiping session_state; a `@st.fragment(run_every="30s")` approach also had edge cases. Now it's a single explicit button.)
4. On the next render, `_auto_merge_if_needed()` (called from `app.py` top level AND from inside `_render_slurm_state_panel`) detects all chunks present + cache missing/stale, fires `merge_slurm_results()` inside a visible `st.status` block, then `st.rerun(scope="app")`.
5. `_slurm_merge_attempted_for: <folder>` marker prevents double-firing. Manual **📦 Merge chunks into cache** button is always available as a fallback / retry.

`_scan_slurm_state()` is **mtime-aware** — a result file is only considered "fresh" if its mtime is newer than its chunk file. This stops stale `[]` results from a previous failed run from being mistaken for completion.

## Key patterns / gotchas

- Pages use `sys.path.insert(0, ...)` to import from `core/` — keep this pattern when adding pages.
- `st.session_state` is heavily used. **Folder-scoped keys** (`selected_prediction`, `slurm_job_ids`, `slurm_num_chunks`, `_slurm_merge_attempted_for`, anything starting with `_detail_`/`_pae_plot_`/`_zoom_plot_`/`_iface_`/`_comp_plot_`) get auto-cleared in `app.py` when `project_path` changes — `_active_project` tracks the current folder.
- Scanner caching uses `@st.cache_data(ttl=300)`.
- Streamlit constraint: never programmatically write to a `session_state` key that is bound to a widget (`key=` parameter) — the widget owns that key.
- `fetch_gene_names_batch()` **must** pre-filter via `is_valid_uniprot_accession()` before sending the OR-query — UniProt's `/search` endpoint returns HTTP 400 if any single accession term has invalid syntax, killing the whole batch of 100. Symptom: gene names silently come back empty. Folders with `fold_` prefixes or non-UniProt IDs (gene-model IDs like `Cre07.g335750.t1.1`) hit this.
- Legacy `af2`-format caches store `iptm: None` + `iptm_ptm`. Always call `normalize_cache_records()` after `json.load()`-ing a cache file to backfill the standard numeric fields.
