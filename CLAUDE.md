# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

PulldownRanker: a Streamlit web app for analyzing protein–protein structure prediction pulldown data.

Supported input formats (auto-detected per folder):

| Format | Detection | Metrics computed |
|---|---|---|
| **AF3 local pipeline** (`format: 'af3'`) | `seed-N_sample-M/` subdirs + `ranking_scores.csv` | full: iPTM, PTM, ipSAE, interface pLDDT, PAE-binned contacts |
| **AF3 Server flat** (`format: 'server_flat'`) | `fold_<name>_model_N.cif` + `fold_<name>_full_data_N.json` directly in folder | full (same as above) |
| **AF2 with PAE CSVs** (`format: 'af2'`) | `ranking_debug.json` + `ranked_N.pdb` + `ranked_N_PAE.csv`, in the folder or in `seq/` | full, same definitions as AF3. `core/af2_pae.py` builds an AF3-shaped confidences dict so the AF3 code paths run unchanged. iPTM = combined `iptm+ptm` |
| **AlphaPulldown** (`format: 'alphapulldown'`) | same, without PAE CSVs | partial: combined `iptm+ptm`, spatial Cβ contacts, interface pLDDT from PDB B-factors. PAE only exists in `result_*.pkl` (~1 GB each), deliberately not read |
| **Legacy `format: 'af2'` caches** (from the retired af2_af3_analysis_app) | cache files with `iptm: None` + `iptm_ptm` | fields are normalized at load time |

Computes confidence metrics across **all 5 seed/sample models per prediction**, not just top-ranked. Models whose `seed-N_sample-M` files were deleted are skipped (`af3_models_on_disk()` in scanner, same check in `analyze_prediction_all_models`); otherwise `load_confidences` falls back to the top-level files and the top model would be reported under every label.

**Chain pair.** Every metric uses one pair: `interface_chain_pair(confidences)` in `analyzer.py` = the first two polymer chains (`polymer_chain_ids`: a chain whose tokens span more than one residue number; AF3 ligands have one token per atom on one residue number, ions one token). Fallback: first two chains of any kind. Structure-based code maps it with `_structure_pair()`. Detailed Analysis stores it as `_detail_pair_ids`; never use `chain_ids_resolved[0/1]` or `protein_names[0/1]` for "chain A/B" again. Chain labels: bait/prey names only when the folder name is `<bait>_and_<prey>`, else `Chain X`; non-polymer chains `RESNAME (X)` and drawn as sticks (`nonpolymer_chains` arg of `generate_viewer_html`).

**Record schema.** Every scan entry, prediction summary and result record has `format` (`af3`, `server_flat`, `af2`, `alphapulldown`). AF2 models whose PAE CSV exists but cannot be used carry `warning` (shown on Results and Detailed Analysis). AF3 local and Server records also carry `chain_pair` (e.g. `A:B`), `n_chains` (polymer chains) and `pair_iptm` (`chain_pair_iptm` for the scored pair, via `chain_pair_fields()`); Results shows Chains / Pair iPTM columns and a warning only when some record has more than two chains. Caches written before this lack the fields and show nothing until re-analysed.

**Planned next: all chain pairs.** Only the first polymer pair is scored now. Plan for automatic chain detection and evaluation of every cross-chain interface: `ROADMAP.md`.

`example_data/IFT27_pulldown/` (two small AF3 pairs) is for first runs and regression checks. `validation/ipsae_vs_reference.py` compares ipSAE with Dunbrack's `ipsae.py` (agreement < 5e-7).

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

- `af2_pae.py` — AF2 helpers: `af2_work_dir()` (folder or `seq/`), `load_af2_confidences()` (PAE CSV + PDB B-factors → AF3-shaped dict, heavy atoms only, aligned with `CIFParser.atoms`).
- `slurm_config.py` — reads optional `slurm_config.toml` (gitignored; template `slurm_config.example.toml`); `slurm_available()` hides the SLURM tab when `sbatch` is missing.
- `scanner.py` — `AF3Scanner` class + format detectors (`_is_af3_server_flat`, `_is_alphapulldown_format`). `load_prediction_data()` dispatches to AF3/server-flat/AlphaPulldown loaders. `find_af3_projects_recursive()` walks a directory tree to discover predictions in any of the supported formats.
- `analyzer.py` — Heavy analysis: `analyze_prediction_all_models()` dispatches by format (`_analyze_server_flat`, `_analyze_alphapulldown`). Computes ipSAE (Dunbrack 2025), interface contacts, interface pLDDT. `_pae_interface_metrics()` is the single contact/iPLDDT routine for AF3 local, AF3 Server and AF2. `CIFParser` is **gemmi-backed**, skips hydrogens, and preserves the legacy attribute-dict shape (`atoms[i] = {seq_id, chain_id, atom_name, res_name, x, y, z}`) so callers that walk `parser.atoms` in CIF order (e.g. for `atom_plddts` indexing) keep working.
- `interface_analyzer.py` — Interface contacts (residue pairs, Cβ–Cβ distance, PAE) and hub residue detection. No residue-type classification (removed 2026-11-18: Cβ-only distances cannot establish bond types). `analyze_interface()` runs PAE + spatial filtering for AF3; `analyze_interface_spatial()` runs spatial-only for AlphaPulldown.
- `pae_plotter.py` — PAE matrix visualization: full matrix, interface zoom, model comparison. Uses green-white colormap.
- `viewer_3d.py` — 3Dmol.js-based 3D structure viewer with PAE interface coloring. **The 3Dmol-min.js bundle is inlined into the HTML via `_load_3dmol_js()` + `@lru_cache`** — no dependency on Streamlit's static serving. Supports both CIF (AF3) and PDB (AlphaPulldown) via `model_format` arg. **Interactive residue selection is entirely client-side**: clicking a residue toggles side-chain sticks + a persistent label, hovering shows a transient label. `setBaseStyles()` (the chain/PAE cartoon coloring) is wrapped in a JS function and re-applied on every selection change, so sticks add/remove cleanly without losing interface colors. No Streamlit round-trip — picks never reach Python, so there's nothing in the server log when you click. Rendered with `st.iframe` when available (Streamlit >= 1.65; `components.html` is deprecated there), else `components.html`. Both run the JavaScript; tested in Chrome on 1.55 and 1.65.
- `pymol_script.py` — PyMOL .py script generator for offline interface visualization.
- `slurm_manager.py` — SLURM job submission and status checking for cluster execution.
- `run_analysis_cli.py` — Standalone CLI (no Streamlit). Finds predictions like the app (`find_predictions_folder` + `AF3Scanner`); exits 1 and lists prediction folders below `--input_dir` when no model was analysed. The SLURM path does not use it (the job-array script embeds its own worker).
- `utils.py` — `format_score()`, `calculate_confidence_tier()`, `tier_color()`, `split_prediction_name()` (strips AF3 Server `fold_` prefix, splits on `_and_`), `is_valid_uniprot_accession()`, `normalize_cache_records()` (backfills `iptm`/`ptm`/`ranking_score` on legacy af2 caches), `fetch_gene_names_batch()` (UniProt API with pre-filter + per-accession fallback on HTTP 400).

## Expected Data Layout

Four formats, auto-detected per folder by the scanner:

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

# AF2-multimer (AlphaPulldown, or AF2 with PAE CSVs); files may be in <bait>_and_<prey>/seq/
<project>/pulldown/<bait>_and_<prey>/
  ranking_debug.json
  ranked_0.pdb … ranked_4.pdb
  ranked_0_PAE.csv … ranked_4_PAE.csv   # optional
```

Analysis results cache to `af3_app_all_models_analysis.json` at the predictions-folder level. Per-prediction `af3_app_analysis.json` files are also written inside each pair folder for fast Detailed-Analysis lookup.

## SLURM execution flow

1. User clicks Submit on the SLURM tab → `submit_slurm_jobs()` writes `_slurm_chunk_<i>.txt`, sbatches ONE job array (`--array=0-N%max_concurrent`, options from `slurm_config.toml`), sets `st.session_state['slurm_job_ids']`, then **`st.rerun()`** so the folder-state panel above the Submit button picks up the new chunk files.
2. Each array task (`$SLURM_ARRAY_TASK_ID` = chunk index) runs the analyzer over its chunk list and writes `_slurm_results_<i>.json`.
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
