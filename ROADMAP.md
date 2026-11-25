# Roadmap

## Multi-chain complexes (all chain pairs)

**Status:** steps 1–9 implemented and checked 2026-11-24 (not yet committed). Step 10 (paper) open. Results of the checks: see "Checks done" at the end.

**Current limitation.** Every metric (ipSAE, iPLDDT, PAE contacts) is computed for one pair only: the first two polymer chains (`interface_chain_pair()` in `core/analyzer.py`). With a 2:2 job, or a bait with several partners, the other interfaces are not scored, and A:B may be two copies of the same protein. Since `2a16b30`, Results shows which pair was scored (Chains column), the pair's own iPTM (Pair iPTM, from `chain_pair_iptm`), and a warning when a record has more than two polymer chains. That marks the problem but does not solve it.

**Goal.** Detect all polymer chains automatically, evaluate every cross-chain interface, and rank each model by its best bait–prey interface.

### Decisions (agreed 2026-11-24)

1. **Bait chains:** by default the entity of chain A and all its copies; every other polymer chain is prey. Folder names with `_with_` (AlphaPulldown convention `bait1_bait2_with_prey`): every chain before the last is bait. The user can override per project ("bait chains: A,B") on the Analyze page; the choice is stored with the results.
2. **Ranking score of a complex:** the best single bait–prey chain pair by ipSAE (Dunbrack's per-pair definition, checkable against `ipsae.py`). No pooling of chains into groups.
3. **Results default:** one row per model (scored by its best bait–prey pair); a toggle shows one row per chain pair.

### The existing script `elo_lab/SCRIPTS/AF3_complex_analysis_v6.py`

Reviewed 2026-11-24. Keep the idea (bait/prey roles, the `_with_` rule, per-chain interface pLDDT). Do not copy its code:

- `contacts_pae8` and `contacts_pae12` are wrong: residue pairs are only considered when PAE ≤ 6 Å (`count_bait_prey_contacts`), so both are in effect the PAE < 6 count.
- ipSAE pools all bait chains against all prey chains. That is not the per-chain-pair definition: d0 depends on the pooled residue count, and the value cannot be compared with `ipsae.py`.
- iPTM is the whole-complex value; `chain_pair_iptm` is not used.
- Ligand and ion chains count as chains (all chains of `token_chain_ids`). With the `_with_` rule a trailing ion chain becomes the prey.
- 3+ chains need `--bait_chains/--prey_chains` or an interactive prompt; in `--auto` mode those folders are skipped.
- Top-ranked model only; residue numbers are taken as index within chain + 1 (assumes numbering from 1 without gaps); interface metrics are skipped when iPTM < 0.3; pure-Python O(N_bait × N_prey) contact loops.

Use it for comparison in the checks (step 9), with the differences above expected.

### Steps

1. **Chains and copies.** `polymer_chain_ids()` gives the polymer chains (ligands and ions excluded: a chain is a polymer when it has more atoms than tokens, AF3 only). Group copies of one protein by mmCIF entity ID (`_struct_asym.entity_id`, read with gemmi in the CIF parse that already happens per model). Note: the `sequences` field that the earlier version of this plan relied on no longer exists (the per-model `data.json` read was removed, see below). AF3 Server mmCIFs carry entity IDs too. AF2 (PDB, no entities): group chains with identical sequences read from the PDB.
2. **Bait assignment.** `bait_chains(confidences, structure, folder_name, override)` implementing decision 1. Roles per pair: `bait-prey`, `bait-bait`, `prey-prey`.
3. **Per-pair metrics.** For every unordered pair (i, j) of polymer chains, per model:
   - ipSAE i→j, j→i and the larger (`calculate_ipsae` already works on any two token-index lists);
   - pair iPTM from `chain_pair_iptm[i][j]` (AF3 local and Server);
   - contacts at PAE ≤ 3, 5, 8 Å and interface pLDDT via `_pae_interface_metrics()`, which already takes two chains' token indices; map to structure chains with `_structure_pair()`.
   Parse the confidences JSON and the CIF once per model and reuse them for all pairs. Expected extra time is small: JSON parsing is most of the per-model time (speed tests 2026-11-24: 5.0 of 10.2 s for the largest ARL3 model, 4,866 tokens).
4. **Record schema.** One record per model, as now. The existing fields (`ipsae`, `contacts_pae3/5/8`, `interface_plddt`, `chain_pair`, `pair_iptm`) hold the summary pair = best bait–prey pair by ipSAE, so two-chain results, Results columns and exports stay the same. New fields:
   - `pairs`: list of `{chains: "A:C", role, entities, ipsae, ipsae_ab, ipsae_ba, pair_iptm, contacts_pae3, contacts_pae5, contacts_pae8, interface_plddt}`;
   - `bait_chains` (list) and `chain_entities` (chain → entity ID).
   Old caches lack `pairs`: in `normalize_cache_records()` build a one-entry `pairs` from the summary fields.
5. **Analyze page.** Optional per-project bait override (text field, e.g. `A,B`), shown only when some prediction has more than two polymer chains. Store it in the result file so a reopened project shows how it was scored.
6. **Results page.** Default one row per model (summary pair, with its chains and role). Toggle: one row per chain pair, with a role filter (bait–prey / bait–bait / prey–prey). Excel export: second sheet with all pairs. Replace the current ">2 chains" warning with a note on how the summary pair was chosen.
7. **Detailed Analysis.** Chain-pair selector replacing the fixed `_detail_pair_ids` (default: summary pair). It drives the PAE zoom, contact table, hub residues, interface regions, 3D colouring and PyMOL script. Add a chain × chain heatmap of ipSAE (and pair iPTM) as an overview of the complex.
8. **AF2, SLURM, CLI.** AF2 runs through the same code (`core/af2_pae.py` already builds all chains); Pair iPTM stays empty (`ranking_debug.json` has no per-pair values). AlphaPulldown without PAE: spatial contacts per pair only. The SLURM job script and `core/run_analysis_cli.py` call `analyze_prediction_all_models`, so they follow; check that `merge_slurm_results()` and the CLI keep `pairs`, `bait_chains`, `chain_entities`.
9. **Checks.**
   - Two-chain results unchanged: all fields on `example_data/IFT27_pulldown/` before and after.
   - Per-pair ipSAE against Dunbrack's `ipsae.py` (reports every chain pair) with `validation/ipsae_vs_reference.py`, extended to all pairs, on a real 3+ chain case: CCDC92/CCDC198 against the IFT-B subcomplexes (`/emcc/au14762/AF/Lotte_Pedersen/CCDC92_CCDC198/IFTB_AF3/`) and one 2:2 job.
   - Side by side with `AF3_complex_analysis_v6.py` on the same folders; differences expected from pooled ipSAE and the PAE 8/12 contact counts.
   - Folder layouts: `validation/folder_layouts.py` still passes.
   - Speed: rerun `manuscript/speed_tests/speed_bench.py sizes` on the same predictions and compare.
10. **Paper.** The PulldownRanker manuscript (local draft `manuscript/manuscript.tex`, outside git; also `essebesse/Publications`, folder `PulldownRanker`) lists the one-pair restriction as the main limitation. Update the Metrics section (chain detection now atoms per token; bait assignment; summary pair), the Limitations section and the Results/Detailed Analysis descriptions.

### Deviations from the plan

- Step 1: copies are grouped by identical residue sequence, labelled with the mmCIF entity ID. Entity IDs alone fail when a protein is entered twice as separate AF3 inputs: in the KIF7 dimer screen (`AF3_APD/Q2M1P5_KIF7_dimer/`) the two KIF7 chains have entity IDs 1 and 2.
- Decision 1, `_with_` rule: the last polymer chain *and its copies* are prey (so `X_with_Y` run as X, X, Y, Y gives bait X, X and prey Y, Y).
- `AlphaPulldown` without PAE: the summary pair is the bait–prey pair with the most spatial contacts (there is no ipSAE).

### Checks done (2026-11-24)

- Two-chain results unchanged, every field: example data (6 models), AF2 with PAE CSVs (IFT27–O00471 from `AF/HsIFT25_27_AF2_screen/`, CSVs made from the npz, 5 models), AlphaPulldown without PAE (IFT25–A0A1B0GUW6, 5 models).
- Per-pair ipSAE against `ipsae.py`: 16 chain pairs (CCDC198_IFT46_56, 3 chains; CCDC198_IFT88_70_52_46, 5 chains; KIF7 dimer with A2A3K4, 2:1), largest difference 4.3e-7. The alphajudge `interfaces.csv` in the CCDC198_IFT46_56 folder has the same A:B and A:C ipSAE.
- `AF3_complex_analysis_v6.py` on CCDC198_IFT46_56 (top model): same per-pair ipSAE (A:B 0.1625, A:C 0.5424); its pooled A vs B+C ipSAE is 0.574; contacts at PAE ≤ 8 are 50 against 64 here (v6 only counts pairs with PAE ≤ 6), PAE ≤ 3 13 against 16, interface pLDDT 85.9 against 86.6.
- `validation/folder_layouts.py --af2 …`: 7 layouts pass.
- Headless Streamlit test of Analyze, Results (both views, role filter) and Detailed Analysis (pair selector, switching pairs) on a scratch copy of the three multi-chain top models plus one example pair.
- Speed (`speed_bench.py one`, same predictions as `sizes`): time and memory unchanged (largest model, 4,866 tokens: 10.2 → 9.5 s, 2.17 GB both).
- Not done: a 2:2 job (none found; the KIF7 2:1 job covers copies and `_with_`). Done 2026-11-25, see below.

### Checks done (2026-11-25)

Test project of file symlinks (deleted afterwards): `AF/HsKCTD10_IFT25_27/runs/KCTD10x2_IFT25-27x2_GTP_Mg` (AF3, KCTD10/IFT25/IFT27 2:2:2 + 2 GTP + 2 Mg, 25 models), `runs/KCTD10x4` (AF3 homotetramer, 25 models), the KCTD10x2 seed-11 models relinked as AF3 Server files (`fold_*`, 5 models), `AF/Vito/CROCC_1_CEP250_N500tetramer` (AF2 2:2 with PAE CSVs, named `crocc_x2_with_cep250_x2`) and the same without PAE CSVs (AlphaPulldown).

- SLURM: `submit_slurm_jobs()` (2 chunks, `%1`, vader9, bait override `A,B`) then `_auto_merge_if_needed()`: 65 models, identical in every field to `run_analysis_cli.py --bait_chains A,B`. Without the override only the homotetramer changes (bait A–D, summary A:B bait-bait; with `A,B`: B:C bait-prey).
- Ligands: 6 polymer chains, 15 pairs, GTP/Mg chains excluded.
- AF3 Server layout vs local layout, same 5 models: identical in every field and pair.
- ipSAE vs `ipsae.py`, top models: KCTD10x2 (15 pairs) and KCTD10x4 (6 pairs), largest difference 4.4e-7; AF2 2:2 CROCC/CEP250, all 5 models, 30 pairs, within the 4-decimal rounding of stored values (≤ 4.9e-5). `validation/ipsae_vs_reference.py` now uses `prediction_file_prefix()` (failed on renamed folders).
- Headless Streamlit test (Analyze, Results both views, Detailed Analysis with pair switching) on all five: no errors. 3D viewer rendered (headless browser screenshots): colouring follows the selected pair, GTP drawn as sticks/spheres.
- Local analysis with the `spawn` start method (macOS/Windows default): runs.
- Scale: synthetic result file of 7,800 models (1,540 preys × 5 + 20 six-chain predictions), 7.2 MB. Results page redraw 1.4 s for a two-chain-only project (same as before); 4.3 s when multi-chain predictions are present, because the XLSX (now with the "Chain pairs" sheet) is rebuilt on every redraw. First load ~30 s, of which 17 s is the UniProt gene-name lookup (unchanged). Open: build the XLSX only on click.

### Fixes from P. Meyke's notes (NCS1 predictions, 2026-11-25)

- 3D viewer coloured every chain pair, including protein–ion pairs, when a model had only one polymer pair: EF-hand residues binding Ca²⁺ were coloured as interface (NCS1–QARS1: 14 residues with no NCS1–QARS1 interface). Now always only the selected pair.
- The same overlay restyled ion chains as cartoon, so ions with confident PAE vanished (NCS1–VPS4A). Ligand chains are now never restyled.
- Model comparison with many models (51, plus a Ca chain): colourbar 70% of a 60-inch figure, empty top band, Ca tick/name labels on top of the prey's. Fixed margins and colourbar height in inches; chains below 3% of tokens (≥5 tokens) get no tick or name label (also in the full PAE plot).
- Not bugs (checked on the top models): WFS1 632–654 lies on NCS1 but has PAE ≥ 17.5 Å, so all metrics are 0; Ric8a 130–176 is 8–22 Å from NCS1 in the RIC8A_apo model. Not reproduced: hub list (note 24), zoom jumping back (note 4).
- Open feature requests: notes 1, 2, 3, 6, 7, 8, 10–14, 17, 19, 20, 22, 23, 25, 26.


### Notes from X. (2026-11-25, tested TTC21A + IFT122 + WDR35 + IFT43 at home)

- Fixed: the "Chain pairs" maps (ipSAE, pair iPTM) had no download button.
- Fixed (also P. Meyke's note 23): Interface Zoom starred every residue pair below the PAE cutoff, with no distance filter (TTC21A–IFT122: 7,288 starred pairs, 25 within 8 Å), so the stars read as contacts. Stars are now only pairs with Cβ–Cβ ≤ 8 Å (`pair_contact_mask()` in `core/analyzer.py`, same definition as the contact counts); legend says so. NCS1–CEP89 test model: 2,805 → 93 stars.
- Fixed: upload a .zip or folder instead of typing a path
- Open: AlphaJudge `interface_meta_score` as a Results column (all AlphaJudge metrics elsewhere). X. offers to do the AlphaJudge part and asks for GitHub access.

### Notes from J. Wang (2026-11-25): declined

- Protein names per chain in Results: not added. Chains stay A, B, C...; the user knows which protein is which.
- One-click interface sticks / polar contacts in the 3D viewer (as in ProtView, Mol*): not added.
