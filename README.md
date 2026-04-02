# AF3 Pulldown Analysis Web App

A Streamlit-based GUI for analyzing AlphaFold 3 pulldown (AF3-PD) prediction data.

## Features

- **Quick Overview**: Browse all predictions with iPTM, PTM, and ranking scores
- **All-Models Analysis**: Analyze all 5 seed/sample models per prediction (not just top-ranked)
- **Sortable Results**: Filter and sort by confidence metrics (iPTM, ipSAE, interface pLDDT)
- **Detailed Analysis**: PAE matrices, interface contacts, hub residues
- **Batch Execution**: Local multiprocessing or SLURM cluster submission
- **Export**: Download results as CSV or XLSX
- **AF3 Server support**: Analyze predictions downloaded directly from the AlphaFold Server

## Local Setup / Installation

It is recommended to use a Python virtual environment to avoid polluting your global packages.

```bash
# Clone or copy the repository
cd /path/to/af3_analysis_app/

# Create a virtual environment (Python 3.11 recommended)
python3 -m venv venv

# Activate the virtual environment
source venv/bin/activate   # Linux / macOS

# Install dependencies
pip install -r requirements.txt
```

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
│   ├── scanner.py               # AF3 directory scanning (local + server format)
│   ├── analyzer.py              # Core analysis logic (ipSAE, contacts, pLDDT)
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
