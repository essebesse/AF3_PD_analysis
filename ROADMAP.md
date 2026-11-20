# Roadmap

## Multi-chain complexes (all chain pairs)

**Status:** planned, not started. To be done on ghost, where real multi-chain predictions are available.

**Current limitation.** Every metric (ipSAE, iPLDDT, PAE contacts) is computed for one pair only: the first two polymer chains (`interface_chain_pair()` in `core/analyzer.py`). With a 2:2 job, or a bait with several partners, the other interfaces are not scored, and A:B may be two copies of the same protein. Since `2a16b30`, Results shows which pair was scored (Chains column), the pair's own iPTM (Pair iPTM, from `chain_pair_iptm`), and a warning when a record has more than two polymer chains. That marks the problem but does not solve it.

**Goal.** Detect all polymer chains automatically and evaluate every cross-chain interface.

Steps:

1. **Chains.** `polymer_chain_ids()` already finds the polymer chains (ligands and ions are excluded). Group chains with identical sequences (copies of one protein) using the `sequences` field, so that a 2:2 job shows bait–prey interfaces and bait–bait / prey–prey interfaces separately.
2. **Per-pair metrics.** For every unordered pair (i, j) of polymer chains: ipSAE (max of both directions; `calculate_ipsae` already works on any two blocks of the PAE matrix), pair iPTM from `chain_pair_iptm`, PAE contacts and iPLDDT via `_pae_interface_metrics()`, which already takes the two chains' token indices. Map to structure chains with `_structure_pair()`.
3. **Record schema.** Keep one record per model with a summary (best pair by ipSAE, plus `chain_pair` / `pair_iptm` as now), and add a `pairs` list with the metrics of every pair. Old caches lack `pairs`; treat them as one pair, as `normalize_cache_records()` does for other legacy fields.
4. **Results.** Default view: one row per model, scored by its best pair (or by the bait–prey pair when the bait is known from the folder name). Option: one row per pair.
5. **Detailed Analysis.** Chain-pair selector, replacing the fixed `_detail_pair_ids`. Feed it to the PAE zoom, contact table, hub residues, 3D view colouring and PyMOL script.
6. **AF2.** Same code through `core/af2_pae.py` (the AF3-shaped confidences dict already carries all chains). `ranking_debug.json` has no per-pair iPTM, so Pair iPTM stays empty for AF2.
7. **SLURM.** The job-array script (written by `pages/batch_execution.py`) imports `analyze_prediction_all_models`, so it follows the analyzer. Check that `merge_slurm_results()` keeps the new `pairs` field, and that `core/run_analysis_cli.py` does too.
8. **Checks.** Two-chain results must not change (compare with `example_data/IFT27_pulldown/` as for `2a16b30`). Add a real 2:2 or 1:1:1 test case. Re-run `validation/ipsae_vs_reference.py`: Dunbrack's `ipsae.py` reports all chain pairs, so every pair can be checked against it.

**Paper.** The PulldownRanker manuscript (`essebesse/Publications`, folder `PulldownRanker`) lists the one-pair restriction as the main limitation. Update the Limitations section and the metrics description when this is done.
