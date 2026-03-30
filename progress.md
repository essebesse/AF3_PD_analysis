# AF3 Analysis App - Progress Report

## Overall Status: Phase 1 ~80% complete, Phases 2-6 scaffolded but non-functional

The Qwen model created the full file structure and a reasonable skeleton of the app, but most of the code is **UI scaffolding only** — the actual analysis logic has significant bugs and the app likely won't run correctly as-is.

---

## What Was Done Well

1. **Project structure matches the plan exactly** — `app.py`, `pages/`, `core/`, `requirements.txt`, `README.md` all created as specified
2. **Streamlit app skeleton** — multi-page layout with sidebar navigation, settings controls, and filter UI all present
3. **`core/scanner.py`** — reasonable implementation of directory scanning with both naming conventions handled
4. **`core/utils.py`** — utility functions (confidence tiers, formatting, UniProt lookup) are clean and correct
5. **`core/analyzer.py`** — ipSAE calculation (`calculate_ipsae`, `calculate_asymmetric_ipsae`) faithfully extracted from the v4 script with correct Dunbrack 2025 math
6. **`requirements.txt` and `README.md`** — complete and accurate

---

## Issues and Bugs Found

### Critical Bugs (will crash or produce wrong results)

| File | Issue | Severity |
|------|-------|----------|
| `core/analyzer.py:241` | **Syntax error**: `extract_chain_indices` return type annotation has mismatched brackets — `Tuple[..., Dict[str, List[int]]:` missing closing `]` | CRASH |
| `core/analyzer.py:66-184` | **CIFParser is broken**: Column header parsing logic is wrong — it looks for lines starting with `#` for headers, but CIF `_atom_site.*` headers start with `_`, not `#`. The parser will fail to read any CIF file correctly | CRASH |
| `core/analyzer.py:156-184` | **`get_cbeta_coords` doesn't actually select CB/CA atoms** — it just stores the last atom coord per residue regardless of atom name. Should filter for `CB` (or `CA` for GLY) | WRONG RESULTS |
| `core/scanner.py:105-106` | **data.json parsing wrong**: Uses `seq_info.get('protein.id')` but the actual JSON structure is nested: `seq_info['protein']['id']`. Dot notation doesn't work on dicts | SILENT FAIL |
| `pages/3_Detailed_Analysis.py:109` | **`st.pyplot(figsize=...)` is not valid Streamlit API** — should create figure with `plt.subplots()` first, then pass to `st.pyplot(fig)` | CRASH |
| `pages/1_Overview.py:125` | **References undefined variables** `partition` and `memory` — these are only defined in the SLURM branch of `app.py` but the page receives them via function args that may not include them | CRASH |
| `app.py:134-148` | **Page imports use wrong module names** — imports `pages._1_Overview` but files are named `1_Overview.py` (leading digit = invalid Python module name). The `_` prefix convention in imports doesn't match the actual filenames | CRASH |
| `core/utils.py:101-123` | **`get_gene_name` uses `requests` but never imports it** — missing `import requests` | CRASH |

### Non-Critical Issues

| File | Issue |
|------|-------|
| `pages/4_Batch_Execution.py:81-137` | **`run_local_analysis` is a stub** — iterates dirs but never actually calls `analyze_prediction_all_models`. Just updates a progress bar with empty results, then saves an empty `[]` to cache |
| `pages/4_Batch_Execution.py:267-291` | **`monitor_job` blocks forever** — uses `while True` + `time.sleep(30)` which will freeze the Streamlit app. Should use `st.status` or async polling |
| `pages/2_Results.py:157-177` | **XLSX export broken** — `df.to_excel()` requires a file path or buffer, not called correctly for download |
| `pages/2_Results.py:54` | **iPTM from summary used for ALL models** — uses the top-level prediction iPTM for every seed/sample row instead of loading per-sample summary_confidences.json |
| `pages/3_Detailed_Analysis.py:94-104` | **Always loads top-level confidences.json** regardless of which seed/sample model is selected. Ignores the `seed`/`sample` selection |
| `core/analyzer.py:331-425` | **`analyze_sample` path logic confused** — takes `sample_path` but calls `sample_path.parent.parent` for data.json, suggesting wrong directory depth assumptions |

---

## Progress by Plan Phase

### Phase 1: Core Infrastructure + Quick Overview (MVP)
| Item | Status | Notes |
|------|--------|-------|
| Streamlit app skeleton with multi-page layout | Done (with import bug) | File naming vs import mismatch |
| `scanner.py`: scan AF3 dirs, load summary JSONs | 80% done | data.json parsing has nested-key bug |
| Page 1 (Overview): folder selection + preview table | Done | Works if scanner bug is fixed |
| Page 2 (Results): sortable table with iPTM/PTM | 70% done | Shows table but uses wrong iPTM for sub-models; export broken |
| `utils.py` | 90% done | Missing `import requests` |

**Phase 1 verdict: ~75% complete. Needs bug fixes to be functional.**

### Phase 2: All-Models Analysis Engine
| Item | Status | Notes |
|------|--------|-------|
| Extract CIFParser from v4 script | Attempted | Parser is broken (wrong column detection, no CB filtering) |
| `calculate_ipsae()` | Done | Math looks correct |
| `all_models_analyzer.py` | NOT CREATED | Plan called for this as a separate module; logic partially in `analyzer.py` instead |
| Contact counting at PAE thresholds | Done | In `analyze_sample()` — logic is correct but uses broken CIFParser |
| Multiprocessing support | Stub only | Code is commented out in batch page |
| JSON cache save/load | Partially done | Save path exists, but analysis never actually runs |

**Phase 2 verdict: ~30% complete. Core math done, but CIF parsing broken and no working pipeline.**

### Phase 3: Detailed Interface Analysis
| Item | Status | Notes |
|------|--------|-------|
| `pae_plotter.py` (generalized PAE plots) | NOT CREATED | Inline matplotlib in page 3 instead; has API bugs |
| `interface_analyzer.py` (generalized contacts) | NOT CREATED | Contact logic is in `analyzer.py` but tied to broken CIF parser |
| Page 3 with tabs | Scaffolded | PAE plot tab crashes; contacts tab may work if parser fixed; model comparison shows only ranking_score |
| Model comparison view (all 5 side-by-side) | Minimal | Shows table only, no PAE plot comparison |

**Phase 3 verdict: ~20% complete. UI structure exists but no working analysis.**

### Phase 4: Complex (Multi-Chain) Support
| Item | Status | Notes |
|------|--------|-------|
| Bait/prey auto-detection | NOT STARTED | |
| Chain assignment UI | NOT STARTED | Radio button exists in sidebar but does nothing |
| Multi-chain analysis | NOT STARTED | |

**Phase 4 verdict: 0% complete.**

### Phase 5: SLURM Integration
| Item | Status | Notes |
|------|--------|-------|
| `slurm_manager.py` | NOT CREATED | SLURM logic is inline in page 4 |
| SLURM script generation | Done | Reasonable template |
| Job submission via `sbatch` | Done | Basic implementation |
| Job monitoring | Done but broken | `while True` loop freezes Streamlit |

**Phase 5 verdict: ~40% complete. Script generation works; monitoring needs rewrite.**

### Phase 6: Polish & Export
| Item | Status | Notes |
|------|--------|-------|
| CSV export | Partially done | CSV works; XLSX broken |
| PyMOL script generation | NOT STARTED | |
| Session persistence | NOT STARTED | |
| Error handling / logging | Minimal | Some try/except blocks |

**Phase 6 verdict: ~15% complete.**

---

## Summary Table

| Phase | Plan Description | Completion | Blocking Issues |
|-------|-----------------|------------|-----------------|
| 1 | Core Infrastructure + Overview | ~75% | Import bug, scanner data.json bug |
| 2 | All-Models Analysis Engine | ~30% | CIFParser broken, no `all_models_analyzer.py` |
| 3 | Detailed Interface Analysis | ~20% | No `pae_plotter.py` or `interface_analyzer.py` |
| 4 | Complex (Multi-Chain) | 0% | Not started |
| 5 | SLURM Integration | ~40% | No `slurm_manager.py`, monitor loop blocks |
| 6 | Polish & Export | ~15% | XLSX broken, no PyMOL, no persistence |

## Recommended Next Steps (Priority Order)

1. **Fix critical bugs** — syntax error in analyzer.py, CIFParser, page imports, data.json parsing
2. **Get Phase 1 actually running** — verify app launches and scans a real AF3 directory
3. **Fix CIFParser** — either properly port from v4 script or switch to using `gemmi` (already in requirements)
4. **Wire up the analysis pipeline** — make "Run Analysis" actually call `analyze_prediction_all_models` with multiprocessing
5. **Fix per-model data loading** — detailed analysis page should load the selected seed/sample, not always the top-level files
6. **Create missing modules** — `all_models_analyzer.py`, `pae_plotter.py`, `interface_analyzer.py`, `slurm_manager.py`
