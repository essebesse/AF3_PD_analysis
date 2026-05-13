# AF3 / AlphaPulldown Analysis Web App

A Streamlit-based GUI for analyzing protein–protein structure prediction pulldown data.
Originally written for AlphaFold 3 (AF3-PD) output; now also handles AlphaFold Server
flat-format folders and AlphaPulldown (AF2-multimer) output.

## Features

- **Quick Overview**: Browse all predictions with iPTM, PTM, and ranking scores
- **All-Models Analysis**: Analyze all 5 seed/sample models per prediction (not just top-ranked)
- **Sortable Results**: Filter and sort by confidence metrics (iPTM, ipSAE, interface pLDDT)
- **Detailed Analysis**: PAE matrices, interface contacts, hub residues, interactive 3D viewer (3Dmol.js, fully offline)
- **Batch Execution**: Local multiprocessing or SLURM cluster submission with auto-merge of chunked results
- **Export**: Download results as CSV / XLSX; per-model PAE PNGs; PyMOL `.py` scripts for offline visualization
- **AF3 Server support**: Analyze predictions downloaded directly from the AlphaFold Server (flat `fold_<name>_*` files)
- **AlphaPulldown / AF2-multimer support**: Partial analysis (combined iptm+ptm, spatial Cβ contacts, pLDDT-from-B-factor interface) for folders containing `ranked_*.pdb` + `ranking_debug.json`
- **UniProt gene-name lookup**: Batch UniProt API lookup with a per-accession fallback so one bad ID doesn't kill 99 good ones; AF3 Server `fold_` prefixes are stripped automatically

## Local Setup / Installation

It is recommended to use an isolated Python environment. Either `venv` or
`conda` works; pick whichever your site supports.

### Option 1 — venv (Python 3.11 recommended)

```bash
# Clone or copy the repository
cd /path/to/af3_analysis_app/

# Ubuntu/Debian only: install the venv stdlib module if missing
sudo apt install python3-venv     # skip on macOS and on clusters where it's already present

# Create and activate the virtual environment
python3 -m venv .venv
source .venv/bin/activate         # Linux / macOS

# Install dependencies
pip install -r requirements.txt
```

### Option 2 — conda (no sudo, useful on shared clusters)

```bash
conda create -n af3_app python=3.11 -y
conda activate af3_app
pip install -r requirements.txt
```

### Install behind a corporate/institutional proxy

Some clusters intercept TLS and pip fails with SSL verification errors. If you
hit that, add the following flags:

```bash
pip install --trusted-host pypi.org --trusted-host files.pythonhosted.org -r requirements.txt
```

### macOS note

On some macOS Python builds you may see a `urllib3 LibreSSL/OpenSSL` warning
the first time a UniProt lookup runs. It is harmless and can be ignored.

### Dependencies

- streamlit >= 1.30
- pandas
- numpy
- matplotlib
- gemmi
- openpyxl
- requests

## Usage

With the virtual environment active, launch the app with the provided script:

```bash
bash run.sh
```

Or manually:

```bash
source .venv/bin/activate
streamlit run app.py

# If `streamlit` is not on your PATH (common after pip install on clusters),
# run it as a module instead:
python -m streamlit run app.py
```

Then open your browser to `http://localhost:8501`

## Expected Data Structure

The scanner auto-detects three formats per folder. You can mix subfolders of
different formats under one project — each gets handled appropriately.

### Local AF3 pipeline output

```
ProjectFolder/
  AF3/
    <bait>_and_<target>/
      ranking_scores.csv
      seed-N_sample-M/
        model.cif
        confidences.json
        summary_confidences.json
```

### AlphaFold Server output (flat format)

Predictions downloaded from the AlphaFold Server can be loaded directly after
unzipping. Point the app at the folder containing the flat files:

```
fold_<name>/
  fold_<name>_model_0.cif
  fold_<name>_model_1.cif
  ...
  fold_<name>_summary_confidences_0.json
  fold_<name>_summary_confidences_1.json
  ...
  fold_<name>_full_data_0.json
  fold_<name>_full_data_1.json
  ...
  fold_<name>_job_request.json
```

### AlphaPulldown / AlphaFold 2 multimer output

Folders produced by AlphaPulldown are recognised by the presence of
`ranking_debug.json` next to `ranked_*.pdb`:

```
ProjectFolder/
  pulldown/
    <bait>_and_<prey>/
      ranking_debug.json
      ranked_0.pdb … ranked_4.pdb
      unrelaxed_model_*_multimer_v3_pred_0.pdb
      timings.json
```

> **Note:** AF2-multimer output does not include a per-residue PAE matrix, so
> only a partial analysis is computed for AlphaPulldown folders:
>
> | Metric | Available | Source |
> |---|---|---|
> | `iPTM` column | combined `iptm+ptm` | `ranking_debug.json` (true iPTM/PTM aren't separable) |
> | Spatial interface contacts | yes | Cβ–Cβ ≤ 8 Å from the PDB |
> | Interface pLDDT | yes | mean B-factor of interface Cβ atoms |
> | ipSAE | **no** | needs PAE |
> | PAE matrix plots | **no** | needs PAE |
> | PAE-binned contact counts | **no** | needs PAE |
>
> A clear banner is shown on the Detailed Analysis page for AlphaPulldown predictions so users aren't confused by the missing tabs.

## Workflow

1. **Load Data** — select a project folder. The scanner lists subfolders containing predictions.
2. **Analyze** — run the per-model analysis. Three options:
   * **Local** — multiprocessing on the current machine; live estimate of runtime based on prediction count, CPUs, and format.
   * **SLURM** — split predictions into N chunks, sbatch one job per chunk to a cluster. After the jobs finish on the cluster, click 🔄 **Refresh** in the SLURM panel: the app detects all chunk files are present and automatically merges them into the analysis cache (with a visible progress box). A manual **📦 Merge chunks into cache** button is always available as a fallback.
3. **Results** — sortable / filterable table of all models with CSV / XLSX export.
4. **Detailed Analysis** — pick one prediction to view PAE plots, interactive 3D structure (3Dmol.js), interface contacts table, hub residues, side-by-side model comparison, and PyMOL script export.

## Project Structure

```
af3_analysis_app/
├── app.py                       # Main Streamlit entry point (folder picker, step routing, folder-scoped state cleanup)
├── run.sh                       # Launcher script (sources .venv, runs streamlit)
├── pages/
│   ├── overview.py              # Step 1: Load Data
│   ├── batch_execution.py       # Step 2: Analyze (local multiprocessing + SLURM with auto-merge)
│   ├── results.py               # Step 3: Results table (filter, sort, CSV/XLSX export)
│   └── detailed_analysis.py     # Step 4: Detailed Analysis (PAE / 3D / contacts / hubs / PyMOL)
├── core/
│   ├── scanner.py               # Directory scanning + format detection (AF3, AF3 Server flat, AlphaPulldown)
│   ├── analyzer.py              # Analysis dispatch by format (ipSAE, contacts, pLDDT); gemmi-backed CIFParser
│   ├── interface_analyzer.py    # Contact classification + hub residues (PAE-aware and spatial-only paths)
│   ├── pae_plotter.py           # PAE matrix visualization
│   ├── viewer_3d.py             # 3Dmol.js 3D viewer (JS bundle inlined into HTML, fully offline)
│   ├── pymol_script.py          # PyMOL .py script generator
│   ├── slurm_manager.py         # SLURM job submission and status checking
│   ├── run_analysis_cli.py      # CLI entry point used by SLURM jobs (bypasses Streamlit)
│   └── utils.py                 # Helpers: score formatting, UniProt batch lookup, prediction-name parsing, legacy-cache normalization
├── static/
│   └── 3Dmol-min.js             # Bundled 3Dmol.js (~526 KB), inlined into the viewer HTML at render time
├── .streamlit/config.toml       # Streamlit settings (sidebar nav off, toolbar in viewer mode)
├── requirements.txt
└── README.md
```

## Analysis Metrics

- **iPTM**: Interface predicted template modeling score. For AlphaPulldown / AF2-multimer the displayed value is the combined `iptm+ptm` from `ranking_debug.json`, since AF2-multimer doesn't separate the two.
- **ipSAE**: Interface structure-based alignment score (Dunbrack 2025). AF3 only.
- **iPLDDT**: Mean pLDDT of interface residues. From AF3 `atom_plddts` or from the PDB B-factor column for AlphaPulldown.
- **PAE contacts (PAE ≤ 3/5/8)**: Number of Cβ–Cβ ≤ 8 Å residue pairs whose PAE is also below 3 / 5 / 8 Å. AF3 only.
- **Spatial contacts**: Cβ–Cβ ≤ 8 Å pairs without a PAE filter. AlphaPulldown.
- **Confidence tiers**: High (ipSAE > 0.7), Medium (0.5 – 0.7), Low (0.3 – 0.5), Very Low (≤ 0.3). Falls back to iPTM-based thresholds (> 0.8 / > 0.6 / > 0.4) when ipSAE is unavailable.

## Result caches and reproducibility

- `af3_app_all_models_analysis.json` — full per-model result list at the predictions-folder level. The Results and Detailed Analysis pages read this.
- `af3_app_analysis.json` — per-prediction subset written inside each pair folder, used as a fast-path by Detailed Analysis.
- `af3_app_analysis_summary.txt` — human-readable summary grouped by confidence tier, written after local analysis.
- Legacy caches (`format: 'af2'` with `iptm: None` + `iptm_ptm`) from older analyzers are normalized at load time — no need to re-run.

## License

MIT License — see [LICENSE](LICENSE).
