# AF3 Pulldown Analysis App: User Guide

Step-by-step instructions for analysing pulldown-style structure predictions
(one bait against many preys) with this web app. Installation is described in
`README.md`.

---

## 1. Getting started

### 1.1. Prerequisites
- Terminal or SSH access to the machine where the app is installed.
- A Python environment with the packages in `requirements.txt` (see `README.md`;
  the `.venv` folder is created during installation, it is not part of the repository).

### 1.2. Launch the app

```bash
cd /path/to/af3_analysis_app
./run.sh
```

Or manually:
```bash
source .venv/bin/activate
streamlit run app.py
```

Then open `http://localhost:8501` in a web browser.

> **Remote server:** forward the port over SSH, then browse to
> `http://localhost:8501` on your own machine:
> ```bash
> ssh -L 8501:localhost:8501 user@server
> ```

### 1.3. Example data

`example_data/IFT27_pulldown/` contains two AF3 predictions: IFT27 with IFT25
(a known heterodimer) and IFT27 with CABYR (no known interaction). Use the full
path to `example_data/IFT27_pulldown` as the project folder in Step 1.

---

## 2. Prepare your data

The app detects the format of each prediction folder automatically.

### Option A: AlphaFold 3, local pipeline output
```
MyProject/
  AF3/
    BAIT_and_PREY/
      ranking_scores.csv
      seed-N_sample-M/
        model.cif
        confidences.json
        summary_confidences.json
```

If the `seed-N_sample-M/` folders were deleted to save space, only the top
model (the files directly in `BAIT_and_PREY/`) is analysed.

### Option B: AlphaFold Server (downloaded zips)
```
MyProject/
  fold_IFT27_and_TUBG1/
    fold_IFT27_and_TUBG1_model_0.cif
    fold_IFT27_and_TUBG1_model_1.cif
    ...
    fold_IFT27_and_TUBG1_full_data_0.json
    ...
```

### Option C: AlphaFold 2 multimer (AlphaPulldown or AF2)
```
MyProject/
  pulldown/
    BAIT_and_PREY/
      ranking_debug.json
      ranked_0.pdb ... ranked_4.pdb
      ranked_0_PAE.csv ... ranked_4_PAE.csv    (optional)
```

The files may also sit in a `seq/` subfolder of the prediction folder.

- **With `ranked_N_PAE.csv` files:** the same analysis as for AF3 (ipSAE, PAE
  contacts, PAE plots).
- **Without them (plain AlphaPulldown output):** a partial analysis. AlphaPulldown
  keeps the PAE matrix only inside its `result_*.pkl` files, which the app does
  not read. ipSAE, PAE contacts and PAE plots are not available; spatial contacts,
  interface pLDDT and the 3D view are. The Detailed Analysis page shows a notice.

---

## 3. The 4-step workflow

### Step 1: Load Data

1. Type or browse to the **project folder** that contains your prediction
   subfolders (for example the folder that holds `AF3/` or `pulldown/`).
2. Click **Select this folder**.
3. The app shows the number of predictions found, a preview table (if an
   analysis from an earlier run exists) and an iPTM histogram.

> **Nested projects:** **Scan for AF3 subprojects** searches the folder tree
> for prediction folders.

---

### Step 2: Analyze

#### A) Local analysis (on the current machine)
1. Choose:
   - **Number of CPUs** (default 8).
   - **PAE cutoff** in Å for ipSAE (default 10).
   - **Analyze all 5 models per prediction:** metrics for every model, not only
     the top-ranked one.
   - **Skip predictions already in cache:** only analyse folders added since the
     last run.
2. Click **Start Local Analysis**. A progress bar shows the run.
3. When it finishes, `af3_app_analysis_summary.txt` is written next to the predictions.

#### B) SLURM cluster
The **SLURM Cluster** tab only appears when the `sbatch` command is available.
Partition, QOS, account, time limit and the Python environment for the compute
nodes are read from `slurm_config.toml` (template: `slurm_config.example.toml`).

1. Choose:
   - **Number of chunks:** the predictions are split into this many array tasks.
   - **Max tasks running at once:** limits load on the cluster and shared storage.
   - **CPUs / memory per task.**
   - The same PAE cutoff and all-models option as for local analysis.
2. Click **Submit to SLURM**. One job array is submitted.
3. Check progress with `squeue -u $USER | grep AF3app`.
4. When the tasks have finished, click **🔄 Refresh** in the SLURM panel. The app
   merges the chunk results into one cache file. **📦 Merge chunks into cache**
   does the same by hand if the automatic merge did not run.

---

### Step 3: Results

1. A sortable table lists every prediction (and every model if all 5 were analysed).
2. Filters at the top of the page:
   - **Confidence tier:** High / Medium / Low / Very Low
   - **iPTM range**
   - **Models:** all models or top-ranked only
   - **Sort by:** ipSAE, iPTM, ranking score, interface pLDDT or contacts
3. Gene names are fetched from UniProt when the folder names contain UniProt accessions.
4. Click a row to select a prediction, then **Go to Detailed Analysis**.
5. The filtered table can be downloaded as **CSV** or **XLSX**.

---

### Step 4: Detailed Analysis

1. Select a prediction (pre-selected if you clicked a row in Results). You can
   search by gene name or UniProt accession.
2. Pick a **model**. The top-ranked model is listed first.
3. Tabs:

#### 3D Structure
- Interactive 3D view (3Dmol.js) with interface residues coloured by PAE.
- Click a residue to show or hide its side chain and a label; hover to identify it.

#### PAE Plots
- Full PAE matrix and a zoom on the chain A / chain B blocks.
- Not available for AF2 predictions without PAE CSV files.

#### Interface Contacts
- Residue pairs with PAE ≤ 10 Å and Cβ–Cβ distance ≤ 8 Å (Cα for glycine). For
  AF2 without PAE: Cβ–Cβ ≤ 8 Å only.
- **Hub residues:** residues with the most contacts on each chain.
- **Interface regions:** contiguous stretches of interface residues on each protein.

#### Model Comparison
- Table of all models on disk (iPTM, PTM, ipSAE, contacts, ranking score).
- PAE matrices side by side (when PAE is available).

#### PyMOL Scripts
- A `.py` script to run in PyMOL (`File → Run Script`), with the commands
  `show_interface`, `show_contacts` and `show_pae_only`.
- Colours: forest green PAE < 3 Å, green 3–5 Å, yellow 5–8 Å, orange 8–12 Å.

#### Export
- PAE matrix as CSV (when PAE is available) and interface contacts as CSV.

---

## 4. Metrics

| Metric | Definition | Available for |
|---|---|---|
| **iPTM** | Interface predicted TM-score from the model. For AF2-multimer: the combined ranking confidence 0.8·iPTM + 0.2·pTM, because `ranking_debug.json` stores only that value. | All formats |
| **ipSAE** | Interface predicted structure alignment error score (Dunbrack 2025), PAE cutoff 10 Å, larger of the two directions. Matches the reference `ipsae.py` (see README, Validation). | AF3; AF2 with PAE CSVs |
| **iPLDDT** | Mean per-atom pLDDT over residues in contacts with PAE ≤ 6 Å. For AF2 without PAE: mean B-factor of interface Cβ atoms. | All formats |
| **PAE ≤ 3/5/8** | Number of residue pairs with Cβ–Cβ ≤ 8 Å and PAE at or below 3, 5 or 8 Å. | AF3; AF2 with PAE CSVs |
| **Cβ ≤ 8** | Residue pairs with Cβ–Cβ < 8 Å, no PAE filter. | AF2 without PAE CSVs |
| **Ranking score** | AF3's model ranking score; for AF2 the combined ranking confidence. | All formats |

### Confidence tiers
Each model is assigned a tier from its ipSAE:

- **High:** ipSAE > 0.7
- **Medium:** ipSAE 0.5–0.7
- **Low:** ipSAE 0.3–0.5
- **Very Low:** ipSAE < 0.3

Without ipSAE (AF2 without PAE CSVs) the tiers use iPTM: > 0.8, > 0.6, > 0.4.

---

## 5. Files written by the app

| File | Content |
|---|---|
| `af3_app_all_models_analysis.json` | All computed metrics, read by Results and Detailed Analysis. |
| `<prediction folder>/af3_app_analysis.json` | Metrics for one prediction (faster loading in Detailed Analysis). |
| `af3_app_analysis_summary.txt` | Text summary grouped by confidence tier. |

---

## 6. Troubleshooting

### No predictions found
- Point the app at the folder that **contains** the prediction folders, not at a
  single prediction folder. (A single AF3 Server or AF2 prediction folder also works.)
- For AlphaFold Server output, the `.cif` and `.json` files must be in the same folder.

### SLURM tasks failed
- Read `AF3_app_chunk*.log` and `*.err` next to the predictions.
- Check `env_setup` in `slurm_config.toml`: the compute nodes must be able to
  activate the app's Python environment.

### Gene names missing
- Folder names must contain valid UniProt accessions (gene-model IDs such as
  `Cre07.g335750.t1.1` cannot be looked up).
- UniProt lookups need internet access; behind a proxy they can fail. This does
  not affect the analysis.

### Adding predictions later
- In **Analyze**, tick **Skip predictions already in cache** so only new folders are processed.

### Comparing all models
- In **Analyze**, tick **Analyze all 5 models per prediction**; in **Results**, set
  the model filter to **All models**.

---

## 7. Command reference

```bash
./run.sh                                  # launch the app
squeue -u $USER | grep AF3app             # SLURM job status
ssh -L 8501:localhost:8501 user@server    # port forwarding
```
