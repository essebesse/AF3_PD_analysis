# AF3 Analysis App — Implementation Tasks for Phases 2-6

## Important Context Before You Start

- **Working directory**: `/home/au14762/AF3_PD_analysis/SCRIPTS/af3_analysis_app/`
- **Virtual environment**: Always activate with `source .venv/bin/activate` before testing
- **Test dataset**: `/home/au14762/AF3_PD_analysis/Q9BW83_IFT27/AF3/` — 1533 predictions, each with 5 seed/sample models
- **Phase 1 is done and working**: `app.py`, `core/scanner.py`, `core/analyzer.py`, `core/utils.py`, all pages
- **Page files are in `pages/`**: `overview.py`, `results.py`, `detailed_analysis.py`, `batch_execution.py` (NOT numbered — they were renamed from `1_Overview.py` etc.)
- **The app runs**: `streamlit run app.py` launches on port 8501

### Critical Rules (Qwen made these mistakes last time — do NOT repeat them)

1. **data.json structure is NESTED**: `seq_info['protein']['id']` NOT `seq_info.get('protein.id')`. The `protein` key contains a dict with `id`, `sequence`, etc.
2. **CIF column headers start with `_atom_site.`** NOT `#`. The parser reads lines starting with `_atom_site.` to find column indices, then data lines start with `ATOM`.
3. **CB atoms**: Always prefer `CB` atom, fall back to `CA` for glycine. Use `parser.get_representative_atom(chain, seq_id)`.
4. **`st.pyplot()` requires a figure object**: Create with `fig, ax = plt.subplots()`, then call `st.pyplot(fig)` and `plt.close(fig)`.
5. **Never use `st.switch_page()`** — navigation is via sidebar radio buttons in `app.py`.
6. **Preserve chain order**: Use list-based ordering (first appearance in `token_chain_ids`), NOT `set()` which randomizes.
7. **Import paths**: Pages import from `core.*` (e.g., `from core.analyzer import ...`). The `sys.path.insert` at top of each page file handles this.
8. **No emojis in metric labels** — Streamlit renders them inconsistently.

---

## TASK 1: Wire Up the All-Models Analysis Pipeline

**Goal**: Make the "Start Analysis" button in `pages/batch_execution.py` actually run `analyze_prediction_all_models()` for every prediction, with multiprocessing, progress tracking, and JSON cache output.

### What exists now
- `core/analyzer.py` has working `analyze_prediction_all_models(pred_dir)` that returns a list of dicts (one per seed/sample model)
- `pages/batch_execution.py` has `run_local_analysis()` that iterates dirs but never calls the analyzer — it just updates a progress bar with empty results

### What to do

1. **Edit `pages/batch_execution.py` function `run_local_analysis()`**:
   - Collect all prediction dirs from `af3_folder` (skip dirs starting with `seed-`)
   - Use `multiprocessing.Pool(num_cpus)` to call `analyze_prediction_all_models` in parallel
   - Since Streamlit reruns on interaction, use `st.status()` container to show progress
   - After completion, flatten all results into one list and save to `{project_path}/AF3/af3_app_all_models_analysis.json`
   - Handle the fact that `analyze_prediction_all_models` loads large confidences.json files (~20-45MB each) — this is memory-intensive. Process in batches if needed.
   - The `sequences` field in results contains full protein sequences — strip these before saving to JSON cache (they're huge and already in data.json)

2. **Important**: `multiprocessing` + Streamlit is tricky. The simplest working approach:
   ```python
   from multiprocessing import Pool
   from core.analyzer import analyze_prediction_all_models

   with Pool(num_cpus) as pool:
       all_results = []
       for i, result_list in enumerate(pool.imap_unordered(analyze_prediction_all_models, pred_dirs)):
           all_results.extend(result_list)
           progress_bar.progress((i + 1) / total)
           status_text.text(f"Analyzed {i+1}/{total} predictions ({len(all_results)} models)")
   ```

3. **Remove the `while True` loop in `monitor_job()`** — it blocks the Streamlit event loop. Replace with a single `squeue` check that the user can re-run by clicking a "Refresh Status" button.

### Test
```bash
source .venv/bin/activate
streamlit run app.py
# Navigate to Batch Execution, set CPUs to 4, click Start
# Should process predictions and save JSON cache
# Then navigate to Results — should show full metrics
```

### Verify
- `af3_app_all_models_analysis.json` is created in `Q9BW83_IFT27/AF3/`
- Each entry has: `prediction_name`, `seed`, `sample`, `iptm`, `ptm`, `ipsae`, `interface_plddt`, `contacts_pae3/5/6/8`, `is_top_ranked`
- Results page loads from cache and shows ipSAE column populated

---

## TASK 2: Fix the Results Page to Use Real Per-Model Data

**Goal**: When cached analysis exists, Results page should show per-model metrics (each seed/sample has its OWN iPTM, ipSAE, etc.), not repeat the top-level prediction iPTM for every row.

### What exists now
- `pages/results.py` function `show_results()` falls back to scanner data when no cache exists
- When using scanner data, it copies the prediction-level iPTM to every seed/sample row (wrong — each model has its own iPTM from its own summary_confidences.json)
- The sort-by dropdown lists "ipSAE" but sorting isn't actually implemented

### What to do

1. **When cache exists**: Load from `af3_app_all_models_analysis.json` — each entry already has per-model `iptm`, `ptm`, `ipsae`, etc. Display directly.

2. **When no cache exists**: Keep current behavior but add a warning banner: "Run full analysis (Batch Execution page) for per-model metrics including ipSAE."

3. **Implement actual sorting**: After building the `filtered_results` list, sort it:
   ```python
   sort_column_map = {
       "ipSAE": "ipsae", "iPTM": "iptm", "ranking_score": "ranking_score",
       "interface_plddt": "interface_plddt", "contacts_pae6": "contacts_pae6"
   }
   sort_key = sort_column_map.get(sort_by, "iptm")
   filtered_results.sort(key=lambda r: r.get(sort_key) or 0, reverse=True)
   ```

4. **Fix XLSX export** — replace `df.to_excel(index=False)` (needs a file path) with:
   ```python
   import io
   buffer = io.BytesIO()
   df.to_excel(buffer, index=False, engine='openpyxl')
   buffer.seek(0)
   st.download_button(label="Download XLSX", data=buffer, ...)
   ```

5. **Add prediction selection**: When user clicks a prediction name, store it in `st.session_state['selected_prediction']` so the Detailed Analysis page can pick it up.

### Test
- With no cache: shows scanner-based overview with warning
- After running analysis: shows per-model metrics, sorting works, export works

---

## TASK 3: Create `core/pae_plotter.py` — Generalized PAE Plot Module

**Goal**: Reusable PAE plot generation that `pages/detailed_analysis.py` calls instead of inline matplotlib code. Based on `SCRIPTS/interface_analysis/create_pae_plots.py` but generalized (no hardcoded protein names).

### What to create: `core/pae_plotter.py`

```python
"""Generalized PAE plot generation for AF3 predictions."""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from typing import List, Optional, Tuple, Dict


def plot_pae_matrix(pae_matrix: np.ndarray,
                    chain_names: List[str],
                    chain_lengths: List[int],
                    title: str = "PAE Matrix",
                    vmax: float = 30.0) -> plt.Figure:
    """
    Plot full PAE matrix with chain boundary lines and labels.

    Args:
        pae_matrix: NxN PAE matrix from confidences.json
        chain_names: List of chain names (e.g., ["IFT27", "RAB23"])
        chain_lengths: List of residue counts per chain
        title: Plot title
        vmax: Maximum PAE value for colorscale

    Returns:
        matplotlib Figure object (caller must do st.pyplot(fig) and plt.close(fig))
    """
    # Implementation:
    # - Use imshow with 'RdYlBu_r' colormap, origin='lower', vmin=0, vmax=vmax
    # - Draw white dashed lines at chain boundaries
    # - Label axes with chain names at midpoints
    # - Add colorbar with label "Expected Position Error (Angstrom)"
    pass


def plot_pae_interface_zoom(pae_matrix: np.ndarray,
                            chain_a_indices: List[int],
                            chain_b_indices: List[int],
                            chain_a_name: str,
                            chain_b_name: str,
                            pae_cutoff: float = 10.0) -> plt.Figure:
    """
    Plot zoomed-in PAE for the interface region only (chain A vs chain B block).

    Extract the off-diagonal block pae_matrix[chain_a_indices][:, chain_b_indices]
    and plot it.
    """
    pass


def plot_all_models_comparison(pae_matrices: List[np.ndarray],
                               model_labels: List[str],
                               chain_lengths: List[int],
                               chain_names: List[str]) -> plt.Figure:
    """
    Side-by-side PAE plots for all 5 models of a prediction.

    Use plt.subplots(1, N, figsize=(4*N, 4)) where N = len(pae_matrices).
    Share colorbar across all subplots.
    """
    pass
```

### Key details from the existing `create_pae_plots.py` to replicate:
- Chain boundaries: cumulative sum of chain_lengths gives boundary positions
- Color scheme: `RdYlBu_r` with vmin=0, vmax=30
- Interface zoom: slice `pae_matrix[a_start:a_end, b_start:b_end]`
- The existing script uses hardcoded names like "ZMYND19" — replace with the `chain_names` parameter
- Get chain names from data.json `sequences[i]['protein']['id']` or use UniProt gene name if available

### Then update `pages/detailed_analysis.py`:
- Import and call `plot_pae_matrix()` instead of inline matplotlib code
- Add the interface zoom plot as a second figure
- In the "Model Comparison" tab, load all 5 models' PAE matrices and call `plot_all_models_comparison()`

---

## TASK 4: Create `core/interface_analyzer.py` — Generalized Interface Analysis

**Goal**: Extract the spatial interface analysis from `SCRIPTS/interface_analysis/spatial_interface_analysis.py` into a generic module. No hardcoded protein names/sequences.

### What to create: `core/interface_analyzer.py`

```python
"""Generalized interface analysis for AF3 predictions."""

import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from core.analyzer import CIFParser, extract_chain_indices


def analyze_interface(cif_path: Path,
                      confidences: dict,
                      pae_cutoff: float = 10.0,
                      distance_cutoff: float = 8.0) -> dict:
    """
    Full interface analysis combining PAE and spatial filtering.

    Returns dict with:
        'contacts': list of contact dicts (chain_a_res, chain_b_res, distance, pae,
                    chain_a_aa, chain_b_aa, interaction_type)
        'hub_residues_a': list of (res_id, res_name, contact_count) sorted by count
        'hub_residues_b': same for chain B
        'chemical_breakdown': dict of interaction_type -> count
        'summary': dict with total_contacts, mean_pae, mean_distance, interface_area_estimate
    """
    pass


def classify_interaction(res_name_a: str, atom_name_a: str,
                         res_name_b: str, atom_name_b: str,
                         distance: float) -> str:
    """
    Classify chemical interaction type.

    Returns one of: 'hydrophobic', 'hydrogen_bond', 'salt_bridge',
                    'aromatic', 'polar', 'other'

    Rules (simplified):
    - salt_bridge: charged pairs (R/K/H with D/E) and dist < 4.0
    - hydrogen_bond: polar atoms (N, O, S) and dist < 3.5
    - hydrophobic: nonpolar residues (A, V, L, I, M, F, W, P) and dist < 5.0
    - aromatic: aromatic residues (F, Y, W, H) both sides and dist < 5.5
    - polar: everything else with dist < 4.0
    - other: fallback
    """
    pass


def find_hub_residues(contacts: List[dict], chain_key: str,
                      min_contacts: int = 3) -> List[Tuple[int, str, int]]:
    """
    Find hub residues (residues with >= min_contacts interface contacts).

    Args:
        contacts: list of contact dicts from analyze_interface
        chain_key: 'chain_a_res' or 'chain_b_res'
        min_contacts: minimum contacts to qualify as hub

    Returns:
        List of (res_id, res_name, contact_count) sorted by count descending
    """
    pass


def map_interface_regions(contacts: List[dict],
                          chain_a_length: int,
                          chain_b_length: int,
                          window: int = 10) -> dict:
    """
    Map interface to sequence regions using a sliding window.

    Returns dict with:
        'chain_a_regions': list of (start, end) tuples of interface regions
        'chain_b_regions': same
        'chain_a_coverage': fraction of chain A at interface
        'chain_b_coverage': fraction of chain B at interface
    """
    pass
```

### Then update `pages/detailed_analysis.py` tab2 (Interface Contacts):
- Replace the inline hub residue calculation with `analyze_interface()`
- Add chemical interaction breakdown as a bar chart (`st.bar_chart`)
- Add interface region mapping display
- Show amino acid types in the contacts table (need res_name from CIFParser)

### Reference: Look at these files for the logic to generalize:
- `SCRIPTS/interface_analysis/spatial_interface_analysis.py` — the PAE+CB spatial filtering logic
- `SCRIPTS/interface_analysis/interface_analysis.py` — chemical interaction classification
- Both have hardcoded protein names/sequences — replace with dynamic detection from data.json

---

## TASK 5: Create `core/slurm_manager.py` — Clean SLURM Module

**Goal**: Extract SLURM logic from `pages/batch_execution.py` into a proper module. Fix the blocking monitor loop.

### What to create: `core/slurm_manager.py`

```python
"""SLURM job management for AF3 analysis."""

import subprocess
from pathlib import Path
from typing import Optional, Dict


def generate_slurm_script(project_path: str,
                          app_dir: str,
                          partition: str = "vader",
                          num_cpus: int = 32,
                          memory: str = "128G",
                          qos: str = "normal",
                          time_limit: str = "23:59:59",
                          analyze_all_models: bool = True) -> str:
    """Generate SLURM batch script content."""
    # The script should:
    # 1. Activate the app's venv
    # 2. Run a standalone analysis script (not streamlit)
    # 3. Save results to the JSON cache
    pass


def submit_job(script_content: str, working_dir: str) -> Dict:
    """
    Submit SLURM job.
    Returns {'success': bool, 'job_id': str or None, 'error': str or None}
    """
    pass


def check_job_status(job_id: str) -> Dict:
    """
    Check SLURM job status (non-blocking, single check).
    Returns {'status': str, 'running': bool, 'completed': bool}
    """
    pass


def get_job_output(job_id: str, working_dir: str) -> Optional[str]:
    """Read job stdout file if it exists."""
    pass
```

### Then update `pages/batch_execution.py`:
- Import from `core.slurm_manager` instead of inline code
- Replace `monitor_job()` `while True` loop with a "Check Status" button that calls `check_job_status()` once per click
- Show job output in an expander when available

### Also create: `core/run_analysis_cli.py`
A standalone CLI script that the SLURM job calls (since Streamlit can't run inside SLURM):
```python
"""Standalone CLI for running AF3 analysis (used by SLURM jobs)."""
import argparse
import json
from pathlib import Path
from multiprocessing import Pool
from analyzer import analyze_prediction_all_models

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', required=True)
    parser.add_argument('--cpus', type=int, default=32)
    parser.add_argument('--pae_cutoff', type=float, default=10.0)
    args = parser.parse_args()

    # Find prediction dirs, run analysis, save JSON cache
    ...
```

---

## TASK 6: Polish — Session State, Export, Error Handling

**Goal**: Make the app feel production-ready.

### 6a. Session state for navigation
- Store `selected_prediction` and `selected_model` in `st.session_state` so switching between pages preserves selection
- In `pages/results.py`, add row selection: when user picks a prediction from the table, store it and show a "View Details" button
- In `pages/detailed_analysis.py`, pre-select from session state if available

### 6b. Fix XLSX export in results page
Replace the broken `df.to_excel()` call:
```python
import io
buffer = io.BytesIO()
df.to_excel(buffer, index=False, engine='openpyxl')
buffer.seek(0)
st.download_button("Download XLSX", data=buffer,
                   file_name="af3_results.xlsx",
                   mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
```

### 6c. Caching with @st.cache_data
- Wrap the scanner.scan() call with `@st.cache_data` to avoid re-scanning on every page interaction:
```python
@st.cache_data(ttl=300)  # 5 min cache
def cached_scan(project_path: str):
    scanner = AF3Scanner(Path(project_path))
    predictions = scanner.scan()
    return predictions, scanner.get_summary_stats()
```

### 6d. Error handling
- Wrap all JSON loads in try/except with user-friendly error messages
- Handle missing files gracefully (some predictions may lack ranking_scores.csv)
- Add `st.toast()` for success/failure notifications

### 6e. Add iPTM distribution plot to Overview page
```python
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(10, 4))
iptm_values = [p['iptm'] for p in predictions]
ax.hist(iptm_values, bins=50, color='steelblue', edgecolor='white')
ax.axvline(x=0.5, color='red', linestyle='--', label='iPTM=0.5 threshold')
ax.set_xlabel('iPTM')
ax.set_ylabel('Count')
ax.set_title('iPTM Distribution')
ax.legend()
st.pyplot(fig)
plt.close(fig)
```

---

## Order of Implementation

**Do these in order — each builds on the previous:**

1. **TASK 1** (analysis pipeline) — most critical, enables everything else
2. **TASK 2** (results page) — needs TASK 1 output
3. **TASK 6b,c,d** (quick polish fixes) — small, independent
4. **TASK 3** (PAE plotter) — standalone module
5. **TASK 4** (interface analyzer) — standalone module
6. **TASK 5** (SLURM manager) — only needed for cluster use
7. **TASK 6a,e** (session state, histogram) — nice-to-have

## How to Test After Each Task

```bash
cd ~/AF3_PD_analysis/SCRIPTS/af3_analysis_app
source .venv/bin/activate

# Syntax check all files:
python3 -c "import py_compile; [py_compile.compile(f, doraise=True) for f in ['app.py', 'core/scanner.py', 'core/analyzer.py', 'core/utils.py', 'pages/overview.py', 'pages/results.py', 'pages/detailed_analysis.py', 'pages/batch_execution.py']]"

# Test imports:
python3 -c "from core.scanner import AF3Scanner; from core.analyzer import *; from core.utils import *; print('OK')"

# Launch app:
streamlit run app.py
```
