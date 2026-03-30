# AF3 Pulldown Analysis Web App

A Streamlit-based GUI for analyzing AlphaFold 3 pulldown (AF3-PD) prediction data.

## Features

- **Quick Overview**: Browse all predictions with iPTM, PTM, and ranking scores
- **All-Models Analysis**: Analyze all 5 seed/sample models per prediction (not just top-ranked)
- **Sortable Results**: Filter and sort by confidence metrics (iPTM, ipSAE, interface pLDDT)
- **Detailed Analysis**: PAE matrices, interface contacts, hub residues
- **Batch Execution**: Local multiprocessing or SLURM cluster submission
- **Export**: Download results as CSV or XLSX

## Installation

```bash
cd /path/to/SCRIPTS/af3_analysis_app/
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

```bash
cd /path/to/SCRIPTS/af3_analysis_app/
streamlit run app.py
```

Then open your browser to `http://localhost:8501`

## Project Structure

```
SCRIPTS/af3_analysis_app/
├── app.py                          # Main Streamlit entry point
├── pages/
│   ├── 1_Overview.py               # Project selection + summary
│   ├── 2_Results.py                # Sortable results table
│   ├── 3_Detailed_Analysis.py      # Per-prediction deep dive
│   └── 4_Batch_Execution.py        # Local/SLURM execution
├── core/
│   ├── __init__.py
│   ├── scanner.py                  # AF3 directory scanning
│   ├── analyzer.py                 # Core analysis logic
│   └── utils.py                    # Shared utilities
├── requirements.txt
└── README.md
```

## Expected Data Structure

The app expects AF3 prediction directories in this format:

```
ProjectFolder/
  AF3/
    <bait>_and_<target>/
      <name>_model.cif
      <name>_confidences.json
      <name>_summary_confidences.json
      <name>_data.json
      ranking_scores.csv
      seed-N_sample-M/
        confidences.json
        model.cif
        summary_confidences.json
```

## Analysis Metrics

- **iPTM**: Interface predicted template modeling score
- **ipSAE**: Interface structure-based alignment score (Dunbrack 2025)
- **ipLDDT**: Interface local distance difference score
- **PAE contacts**: Number of residue pairs below PAE threshold
- **Confidence tiers**: High (>0.7), Medium (0.5-0.7), Low (0.3-0.5)

## Based On

This app extracts and adapts core analysis functions from:
- `AF3_PD_analysis_v4.py` - Main batch analysis with ipSAE
- `AF3_complex_analysis_v4.py` - Multi-subunit analysis
- `interface_analysis/` - Interface contact analysis

## License

Same as the parent AF3_PD_analysis project.