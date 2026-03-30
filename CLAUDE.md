# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

A Streamlit web app for analyzing AlphaFold 3 pulldown (AF3-PD) prediction data. It computes confidence metrics (iPTM, ipSAE, interface pLDDT) across all 5 seed/sample models per prediction, not just top-ranked.

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

- `scanner.py` — `AF3Scanner` class: reads `summary_confidences.json` and `ranking_scores.csv` from each prediction directory for lightweight overview. Also contains `find_af3_projects_recursive()` for recursive subproject discovery
- `analyzer.py` — Heavy analysis: `analyze_prediction_all_models()` computes ipSAE (Dunbrack 2025), interface contacts, pLDDT using custom CIF parsing. Contains `CIFParser` and `extract_chain_indices`
- `interface_analyzer.py` — Contact classification (hydrophobic/polar/charged/aromatic), hub residue detection. Imports from `analyzer.py`
- `pae_plotter.py` — PAE matrix visualization: full matrix, interface zoom, model comparison. Uses green-white colormap
- `viewer_3d.py` — 3Dmol.js-based 3D structure viewer with PAE interface coloring, embedded via `st.components.v1.html()`
- `pymol_script.py` — PyMOL .py script generator for offline interface visualization
- `slurm_manager.py` — SLURM job submission and status checking for cluster execution
- `run_analysis_cli.py` — Standalone CLI entry point for SLURM jobs (bypasses Streamlit)
- `utils.py` — `format_score()`, `calculate_confidence_tier()`, `tier_color()`, UniProt gene name batch lookup

## Expected Data Layout

The app expects: `<project>/AF3/<bait>_and_<target>/` with `seed-N_sample-M/` subdirectories containing `model.cif`, `confidences.json`, `summary_confidences.json`.

Analysis results are cached to `af3_app_all_models_analysis.json` in the project folder.

## Known Limitations

- The 3D structure viewer (`viewer_3d.py`) loads 3Dmol.js from `cdn.jsdelivr.net`. On air-gapped/isolated networks, the viewer tab shows a fallback error message. All other tabs (PAE plots, contacts, exports, PyMOL scripts) work fully offline. Fix if needed: bundle `3Dmol-min.js` as a Streamlit static file and serve locally.

## Key Patterns

- Pages use `sys.path.insert(0, ...)` to import from `core/` — keep this pattern when adding pages
- `st.session_state` is used heavily for cross-step state (e.g., `selected_prediction`, `navigate_to_detailed`, `analysis_results`)
- Scanner caching uses `@st.cache_data(ttl=300)`
- Streamlit constraint: never programmatically write to a `session_state` key that is bound to a widget (`key=` parameter) — the widget owns that key
