# AlphaFold Pulldown Analysis Web App

A Streamlit-based GUI for analyzing AlphaFold **3** and AlphaFold **2**
pulldown prediction data.

> **Pairwise predictions only.** This app analyzes the interface between the
> first two chains of each prediction. Predictions with three or more chains
> are still listed, but their interface scoring is limited to Chain A vs
> Chain B — a ⚠️ marker flags those rows.

## Features

- **Quick Overview**: Browse all predictions with iPTM / iPTM+pTM and ranking scores
- **All-Models Analysis**: Analyze all 5 seed/sample (AF3) or ranked_* (AF2) models per prediction
- **Sortable Results**: Filter and sort by confidence metrics (iPTM, ipSAE, interface pLDDT)
- **Detailed Analysis (AF3 only for now)**: PAE matrices, interface contacts, hub residues
- **Batch Execution**: Local multiprocessing (all formats) or SLURM cluster submission (AF3 only)
- **Export**: Download results as CSV or XLSX
- **Multiple formats supported**:
  - AF3 local pipeline (`seed-*_sample-*` subdirs + `ranking_scores.csv`)
  - AF3 Server flat format (unzipped folder of `fold_*_model_N.cif` etc.)
  - AF2 standard output (`ranked_*.pdb` + `ranked_*_PAE.csv` + `ranking_debug.json`),
    either flat in the prediction folder or inside a `seq/` subdirectory

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
python3 -m venv venv
source venv/bin/activate          # Linux / macOS

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
source venv/bin/activate
streamlit run app.py

# If `streamlit` is not on your PATH (common after pip install on clusters),
# run it as a module instead:
python -m streamlit run app.py
```

Then open your browser to `http://localhost:8501`

## Expected Data Structure

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

### AlphaFold 2 standard output

Two layouts are supported. The app auto-detects which one is in use.

**Flat layout** — ranked files directly inside the prediction folder:

```
ProjectFolder/
  <pred_name>/
    ranked_0.pdb
    ranked_0_PAE.csv
    ranked_0_pLDDT.csv
    ranked_1.pdb
    ...
    ranking_debug.json
```

**`seq/` layout** — ranked files inside a `seq/` subdirectory:

```
ProjectFolder/
  <pred_name>/
    seq/
      ranked_0.pdb
      ranked_0_PAE.csv
      ranked_0_pLDDT.csv
      ...
      ranking_debug.json
```

Detailed Analysis (PAE plot, 3D viewer, interface table) is **not yet**
wired up for AF2 — screening metrics (ipSAE, PAE contacts, interface
pLDDT) are fully supported on the **Results** page.

## Project Structure

```
af3_analysis_app/
├── app.py                       # Main Streamlit entry point
├── run.sh                       # Launcher script
├── pages/
│   ├── overview.py              # Step 1: Load Data
│   ├── batch_execution.py       # Step 2: Analyze
│   ├── results.py               # Step 3: Results table
│   └── detailed_analysis.py     # Step 4: Detailed Analysis
├── core/
│   ├── scanner.py               # Directory scanning (AF3 local + server + AF2)
│   ├── analyzer.py              # AF3 analysis logic (ipSAE, contacts, pLDDT) + shared helpers
│   ├── af2_analyzer.py          # AF2 analysis logic (uses shared helpers)
│   ├── interface_analyzer.py    # Contact classification and hub residues
│   ├── pae_plotter.py           # PAE matrix visualization
│   ├── viewer_3d.py             # 3Dmol.js 3D structure viewer
│   ├── pymol_script.py          # PyMOL script generator
│   ├── slurm_manager.py         # SLURM job submission
│   └── utils.py                 # Shared utilities
├── requirements.txt
└── README.md
```

## Analysis Metrics

- **iPTM**: Interface predicted template modeling score
- **ipSAE**: Interface structure-based alignment score (Dunbrack 2025)
- **ipLDDT**: Interface local distance difference score
- **PAE contacts**: Number of residue pairs below PAE threshold
- **Confidence tiers**: High (>0.7), Medium (0.5-0.7), Low (0.3-0.5)

## License

MIT License — see [LICENSE](LICENSE).
