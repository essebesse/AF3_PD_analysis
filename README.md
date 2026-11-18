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
- **AlphaFold 2 multimer support**: full analysis (ipSAE, PAE contacts, PAE plots) when `ranked_N_PAE.csv` files are present; partial analysis (combined iptm+ptm, spatial Cβ contacts, B-factor pLDDT) for plain AlphaPulldown folders (`ranked_*.pdb` + `ranking_debug.json`)
- **ipSAE checked against the reference implementation**: see [Validation](#validation)
- **UniProt gene-name lookup**: Batch UniProt API lookup with a per-accession fallback so one bad ID doesn't kill 99 good ones; AF3 Server `fold_` prefixes are stripped automatically

## Local Setup / Installation

It is recommended to use an isolated Python environment. Either `venv` or
`conda` works; pick whichever your site supports.

### Option 1: venv (Python 3.11 recommended)

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

### Option 2: conda (no sudo, useful on shared clusters)

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

### Try it on the example data

`example_data/IFT27_pulldown/` holds two AF3 predictions from a human IFT27
pulldown (12 MB):

| Folder | Pair | Models on disk | Top-model iPTM | Top-model ipSAE |
|---|---|---|---|---|
| `q9bw83_and_q9y547` | IFT27–IFT25, a known heterodimer | all 5 | 0.84 | 0.65 |
| `q9bw83_and_o75952` | IFT27–CABYR, no known interaction | top model only | 0.39 | 0.02 |

In step 1, enter the full path to `example_data/IFT27_pulldown` as the project
folder. These files are AlphaFold 3 output and are covered by the AlphaFold 3
Output Terms of Use (`TERMS_OF_USE.md` in each folder), not by the MIT license
of the code.

## Expected Data Structure

The scanner auto-detects three formats per folder. You can mix subfolders of
different formats under one project; each gets handled appropriately.

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

### AlphaFold 2 multimer output (AlphaPulldown or AF2)

A folder is recognised as AF2-multimer output when it holds
`ranking_debug.json` and `ranked_0.pdb`, either directly or in a `seq/`
subfolder:

```
ProjectFolder/
  pulldown/
    <bait>_and_<prey>/
      ranking_debug.json
      ranked_0.pdb … ranked_4.pdb
      ranked_0_PAE.csv … ranked_4_PAE.csv     # optional, enables full analysis
```

What can be computed depends on whether a PAE matrix is available as
`ranked_N_PAE.csv` (an L×L matrix with a header row and an index column, one
row per residue in PDB order):

| Metric | With `ranked_N_PAE.csv` | Without (plain AlphaPulldown) |
|---|---|---|
| `iPTM` column | combined `iptm+ptm` from `ranking_debug.json` | same |
| ipSAE | yes | no |
| PAE-binned contacts (PAE ≤ 3/5/8 Å) | yes | no |
| Spatial contacts (Cβ–Cβ ≤ 8 Å, no PAE filter) | no | yes |
| Interface pLDDT | mean heavy-atom pLDDT (B-factor) of residues in contacts with PAE ≤ 6 Å | mean B-factor of interface Cβ atoms |
| PAE plots, PAE-coloured 3D view, PAE export | yes | no |

AF2-multimer's `ranking_debug.json` only holds the combined ranking confidence
(0.8·iPTM + 0.2·pTM), so the app shows that value in the iPTM and PTM columns.
AlphaPulldown stores the PAE matrix only inside the `result_*.pkl` files
(about 1 GB each); the app does not read those.

## Workflow

1. **Load Data**: select a project folder. The scanner lists subfolders containing predictions.
2. **Analyze**: run the per-model analysis. Two options:
   * **Local**: multiprocessing on the current machine; live estimate of runtime based on prediction count, CPUs, and format.
   * **SLURM**: only shown when `sbatch` is on PATH. Splits the predictions into N chunks and submits them as one job array with a limit on how many tasks run at once. After the tasks finish, click 🔄 **Refresh** in the SLURM panel: the app detects all chunk files are present and automatically merges them into the analysis cache (with a visible progress box). A manual **📦 Merge chunks into cache** button is always available as a fallback. See [SLURM settings](#slurm-settings).
3. **Results**: sortable / filterable table of all models with CSV / XLSX export.
4. **Detailed Analysis**: pick one prediction to view PAE plots, interactive 3D structure (3Dmol.js), interface contacts table, hub residues, side-by-side model comparison, and PyMOL script export.

## SLURM settings

Cluster-specific options live in `slurm_config.toml` in the app folder. Copy
`slurm_config.example.toml` and set your partition, QOS, account, time limit,
CPU/memory choices and the line that activates the app's Python environment on
a compute node. Every key is optional; empty values leave the `#SBATCH` option
out so the cluster default applies. Without the file, the app's own `.venv` is
used and no partition/QOS/account is requested. `slurm_config.toml` is not
tracked by git.

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
│   ├── scanner.py               # Directory scanning + format detection (AF3, AF3 Server flat, AF2)
│   ├── analyzer.py              # Analysis dispatch by format (ipSAE, contacts, pLDDT); gemmi-backed CIFParser
│   ├── af2_pae.py               # AF2 input: seq/ layout, ranked_N_PAE.csv → AF3-shaped confidences
│   ├── slurm_config.py          # Reads slurm_config.toml (site-specific SLURM options)
│   ├── interface_analyzer.py    # Contact classification + hub residues (PAE-aware and spatial-only paths)
│   ├── pae_plotter.py           # PAE matrix visualization
│   ├── viewer_3d.py             # 3Dmol.js 3D viewer (JS bundle inlined into HTML, fully offline)
│   ├── pymol_script.py          # PyMOL .py script generator
│   ├── slurm_manager.py         # SLURM job submission and status checking
│   ├── run_analysis_cli.py      # CLI entry point used by SLURM jobs (bypasses Streamlit)
│   └── utils.py                 # Helpers: score formatting, UniProt batch lookup, prediction-name parsing, legacy-cache normalization
├── static/
│   └── 3Dmol-min.js             # Bundled 3Dmol.js (~526 KB), inlined into the viewer HTML at render time
├── example_data/                # Two AF3 predictions for a first run (see above)
├── validation/
│   └── ipsae_vs_reference.py    # ipSAE: this app vs Dunbrack's ipsae.py
├── slurm_config.example.toml    # Template for site-specific SLURM options
├── .streamlit/config.toml       # Streamlit settings (sidebar nav off, toolbar in viewer mode)
├── requirements.txt
└── README.md
```

## Analysis Metrics

- **iPTM**: Interface predicted template modeling score. For AF2-multimer the displayed value is the combined ranking confidence (0.8·iPTM + 0.2·pTM) from `ranking_debug.json`, since that file doesn't store the two separately.
- **ipSAE**: Interface predicted structure alignment error score (Dunbrack 2025), PAE cutoff 10 Å by default, the larger of the two directions (A→B, B→A). AF3, and AF2 with PAE CSVs.
- **iPLDDT**: Mean per-atom pLDDT over residues in contacts with PAE ≤ 6 Å. From AF3 `atom_plddts`, or from the PDB B-factor column for AF2 (hydrogens excluded).
- **PAE contacts (PAE ≤ 3/5/8)**: Number of Cβ–Cβ ≤ 8 Å residue pairs (Cα for glycine) whose PAE is also at or below 3 / 5 / 8 Å. AF3, and AF2 with PAE CSVs.
- **Spatial contacts**: Cβ–Cβ ≤ 8 Å pairs without a PAE filter. AF2 without PAE CSVs.
- **Models analysed**: all models listed in `ranking_scores.csv` whose `seed-N_sample-M` files are on disk. When those folders were deleted to save space, only the top model (the top-level files) is analysed.
- **Confidence tiers**: High (ipSAE > 0.7), Medium (0.5 – 0.7), Low (0.3 – 0.5), Very Low (< 0.3). Falls back to iPTM-based thresholds (> 0.8 / > 0.6 / > 0.4) when ipSAE is unavailable.

## Result caches and reproducibility

- `af3_app_all_models_analysis.json`: full per-model result list at the predictions-folder level. The Results and Detailed Analysis pages read this.
- `af3_app_analysis.json`: per-prediction subset written inside each pair folder, used as a fast-path by Detailed Analysis.
- `af3_app_analysis_summary.txt`: human-readable summary grouped by confidence tier, written after local analysis.
- Legacy caches (`format: 'af2'` with `iptm: None` + `iptm_ptm`) from older analyzers are normalized at load time; no need to re-run.

## Validation

`validation/ipsae_vs_reference.py` scores every model in one or more AF3
folders with both this app and the reference `ipsae.py` (version 4) from
[DunbrackLab/IPSAE](https://github.com/DunbrackLab/IPSAE), PAE cutoff 10 Å:

```bash
curl -LO https://raw.githubusercontent.com/DunbrackLab/IPSAE/main/ipsae.py
python validation/ipsae_vs_reference.py ipsae.py                # example data
python validation/ipsae_vs_reference.py ipsae.py /path/to/AF3   # your own folder
```

On the 7 models of the example data the largest difference is 4.7×10⁻⁷. The
same agreement (below 5×10⁻⁷) was found for the 5 models of an AF2-multimer
prediction with PAE CSVs (ARL13B–INPP5E, not bundled).

## License

MIT License, see [LICENSE](LICENSE).
