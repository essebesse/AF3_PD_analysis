# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

A Streamlit web app for analyzing AlphaFold 3 **and AlphaFold 2** pulldown
prediction data. It computes confidence metrics (iPTM / iPTM+pTM, ipSAE,
interface pLDDT) across all 5 models per prediction, not just top-ranked.

**Pairwise only**: the app scores the interface between the first two chains
of each prediction. Predictions with 3+ chains are still listed, but their
scoring is limited to Chain A vs Chain B (flagged with ⚠️ in the UI).

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
- `pages/detailed_analysis.py` — Step 4 "Detailed Analysis": PAE plots, interface contacts, hub residues, 3D viewer, PyMOL script export. Format-aware: routes data loading through `core/af2_detail.py` for AF2 and `core/analyzer.py` + `core/interface_analyzer.py` for AF3.

**Core modules** (`core/`):

- `scanner.py` — `AF3Scanner` class: lists prediction directories and tags each with `format` (`af3_local` / `af3_server_flat` / `af2`), `n_chains`, `pairwise_ok`, and (for AF2) `work_dir`. Also contains `detect_format()` and `find_predictions_recursive()` (with `find_af3_projects_recursive` kept as a backwards-compatible alias).
- `analyzer.py` — AF3 analysis: `analyze_prediction_all_models()` is the format-aware dispatcher that routes to AF3 local, AF3 Server flat, or AF2. Computes ipSAE (Dunbrack 2025), interface contacts (via the shared `count_interface_contacts()` helper), and pLDDT. Contains `CIFParser` and `extract_chain_indices`.
- `af2_analyzer.py` — AF2 analysis: `analyze_af2_prediction_all_models()`. Uses gemmi to parse `ranked_*.pdb` and reuses `calculate_ipsae` + `count_interface_contacts` from `analyzer.py`. AF2 results have `format='af2'`, `seed=0`, and `None` for AF3-only fields (iptm, ptm, ranking_score, fraction_disordered); iPTM+pTM is stored as `iptm_ptm`.
- `interface_analyzer.py` — Contact classification (hydrophobic/polar/charged/aromatic), hub residue detection. Imports from `analyzer.py`. Used by the AF3 Detailed Analysis path; AF2 uses `af2_detail.py` which reuses the same `classify_interaction` and `find_hub_residues` helpers via gemmi.
- `af2_detail.py` — Detailed-analysis adapters for AF2: `load_af2_model_data()` synthesizes an AF3-shaped `confidences` dict (PAE matrix, token_chain_ids, token_res_ids, atom_plddts) from `ranked_N_PAE.csv` + `ranked_N.pdb` + `ranked_N_pLDDT.csv`, so the existing PAE plot / 3D viewer code can render AF2 models without branching. Also provides `analyze_af2_interface()` and `compute_af2_interface_pae_per_residue()` (gemmi-based mirrors of the AF3 helpers).
- `pae_plotter.py` — PAE matrix visualization: full matrix, interface zoom, model comparison. Uses green-white colormap
- `viewer_3d.py` — 3Dmol.js-based 3D structure viewer with PAE interface coloring, embedded via `st.components.v1.html()`
- `pymol_script.py` — PyMOL .py script generator for offline interface visualization
- `slurm_manager.py` — SLURM job submission and status checking for cluster execution
- `run_analysis_cli.py` — Standalone CLI entry point for SLURM jobs (bypasses Streamlit)
- `utils.py` — `format_score()`, `calculate_confidence_tier()`, `tier_color()`, UniProt gene name batch lookup

## Expected Data Layout

Three formats are auto-detected per prediction folder:

1. **AF3 local pipeline** (default): `<project>/AF3/<bait>_and_<target>/` with `seed-N_sample-M/` subdirectories containing `model.cif`, `confidences.json`, `summary_confidences.json`.
2. **AF3 Server flat**: a single folder with `fold_<name>_model_N.cif` and `fold_<name>_summary_confidences_N.json` files directly inside (no seed subdirs).
3. **AF2 standard**: `ranked_N.pdb` + `ranked_N_PAE.csv` + `ranked_N_pLDDT.csv` + `ranking_debug.json`, either directly in the prediction folder *or* inside a `seq/` subdirectory.

Analysis results are cached to `af3_app_all_models_analysis.json` in the
project folder (name kept for backwards compatibility even when the cache
contains AF2 results). Each entry carries a `format` field for downstream
pages to key off. Old caches without `format` are treated as `af3_local` on
load.

## Known Limitations

- The 3D structure viewer (`viewer_3d.py`) **inlines** `static/3Dmol-min.js` into each iframe rather than loading it via `<script src>`. Streamlit's static server returns `Content-Type: text/plain` plus `X-Content-Type-Options: nosniff`, which makes browsers refuse to execute the file as JavaScript. Inlining sidesteps that and keeps the viewer fully offline. If the file is missing, the viewer tab renders a clear "source not found" message; other tabs still work.

## Key Patterns

- Pages use `sys.path.insert(0, ...)` to import from `core/` — keep this pattern when adding pages
- `st.session_state` is used heavily for cross-step state (e.g., `selected_prediction`, `navigate_to_detailed`, `analysis_results`)
- Scanner caching uses `@st.cache_data(ttl=300)`
- Streamlit constraint: never programmatically write to a `session_state` key that is bound to a widget (`key=` parameter) — the widget owns that key

## User-facing docs

- `Alpha_Pulldown_Analysis_User_Manual.txt` — end-user manual covering installation, workflow, metrics, troubleshooting. Update this when adding user-visible features (renamed from `AF3_Analysis_App_User_Manual.txt` when AF2 support landed).
- `paper_text_Alpha_Pulldown_methods.txt` — methods/results paragraphs for the upcoming paper. Keep the metric formulas and tier thresholds in sync with the implementation in `core/analyzer.py`.
- `README.md` — short setup + usage. Quick-start, not detail.
