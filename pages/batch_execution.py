"""
Step 2: Analyze - Local Multiprocessing Execution

Runs full analysis on all predictions using local CPUs.
"""

import streamlit as st
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _auto_merge_if_needed(af3_folder: str) -> bool:
    """If a SLURM run's chunks are all written and the cache is missing or
    stale, merge them in place. Returns True if a merge was performed (so
    callers may want to ``st.rerun()``).
    """
    state = _scan_slurm_state(af3_folder)
    expected = state['expected_chunks']
    if expected == 0:
        return False
    if state['completed_chunks'] < expected:
        return False  # still waiting on jobs
    if state['empty_chunks'] == expected:
        return False  # nothing useful to merge
    if state['cache_is_current']:
        return False  # already merged

    # Guard against running twice in the same session for this folder
    if st.session_state.get('_slurm_merge_attempted_for') == af3_folder:
        return False
    st.session_state['_slurm_merge_attempted_for'] = af3_folder

    non_empty = expected - state['empty_chunks']
    with st.status(
        f"SLURM run finished — merging {non_empty} chunk(s) into analysis cache...",
        expanded=True,
    ) as status:
        try:
            merge_slurm_results(af3_folder, expected)
            status.update(label=f"✅ Merge complete — {non_empty} chunks combined.",
                          state="complete")
            return True
        except Exception as e:
            status.update(label=f"❌ Merge failed: {e} — go to Analyze tab and click manual merge.",
                          state="error")
            # Marker stays set so we don't loop; user can manually retry.
            return False


def show_analyze(project_path: str, af3_folder: str):
    """Display the Analyze step."""

    st.header("🚀 Run Analysis")

    st.markdown("""
    Run full analysis on all predictions. This will:
    - Calculate ipSAE scores for all models
    - Compute interface contacts with spatial filtering
    - Generate PAE matrices and interface pLDDT metrics
    - Save results to cache file for quick loading
    """)

    if not os.path.isdir(af3_folder):
        st.error(f"Predictions folder not found: {af3_folder}")
        return

    # If a SLURM run finished while the user was elsewhere, merge it now
    # (regardless of which sub-tab is active). Re-runs the page after merge
    # so the cache-aware UI below picks up the new state.
    if _auto_merge_if_needed(af3_folder):
        st.rerun()

    # Use the scanner — handles AF3 local pipeline, AF3 Server flat,
    # AlphaPulldown, and the "folder is itself one prediction" case.
    from core.scanner import AF3Scanner
    try:
        _preds = AF3Scanner(Path(af3_folder)).scan()
    except Exception:
        _preds = []
    pred_count = len(_preds)
    # Detect dominant format for time-estimate heuristics
    _fmt_counts = {}
    for p in _preds:
        _fmt_counts[p.get('format', 'af3')] = _fmt_counts.get(p.get('format', 'af3'), 0) + 1
    dominant_format = max(_fmt_counts, key=_fmt_counts.get) if _fmt_counts else 'af3'

    st.info(f"Found {pred_count} predictions to analyze ({dominant_format})")

    # Show prominent notice when analysis cache already exists
    cache_file = Path(af3_folder) / "af3_app_all_models_analysis.json"
    already_analyzed = cache_file.exists()
    if already_analyzed:
        import json
        try:
            with open(cache_file, 'r') as f:
                cached_count = len(json.load(f))
        except Exception:
            cached_count = 0
        st.success(
            f"Results already available ({cached_count} models in cache). "
            "You can go directly to Results, or re-analyze below if needed."
        )
        if st.button("Next: View Results →", type="primary"):
            st.session_state['_navigate_to'] = "3. Results"
            st.rerun()

    st.divider()

    tab_local, tab_slurm = st.tabs(["🖥️ Local (this machine)", "🖧 SLURM Cluster"])

    with tab_local:
        show_local_execution(project_path, af3_folder, pred_count, dominant_format)

    with tab_slurm:
        show_slurm_execution(project_path, af3_folder, pred_count)


def _format_duration(seconds: float) -> str:
    """Pretty-print a duration as 's', 'm:ss', or 'h:mm:ss'."""
    seconds = max(0, int(seconds))
    if seconds < 60:
        return f"{seconds}s"
    if seconds < 3600:
        return f"{seconds // 60}m {seconds % 60:02d}s"
    return f"{seconds // 3600}h {(seconds % 3600) // 60:02d}m"


# Rough wall-clock estimate (seconds/prediction, single-CPU equivalent)
# calibrated from real runs: AF3 top ~0.5s, AF3 all-models ~3s,
# AlphaPulldown ~0.1s (no PAE math), AF2-with-PAE ~0.4s.
_PER_PRED_SECS = {
    ('af3', False):                0.6,
    ('af3', True):                 3.0,
    ('server_flat', False):        0.6,
    ('server_flat', True):         3.0,
    ('alphapulldown', False):      0.15,
    ('alphapulldown', True):       0.4,
    ('af2', False):                0.4,
    ('af2', True):                 1.5,
}


def _estimate_local_runtime(pred_count: int, num_cpus: int, analyze_all: bool,
                             dominant_format: str) -> str:
    if pred_count <= 0:
        return "—"
    unit = _PER_PRED_SECS.get((dominant_format, analyze_all),
                              _PER_PRED_SECS[('af3', analyze_all)])
    # Multiprocessing overhead: ~20% loss at high CPU counts
    effective_cpus = max(1, num_cpus * 0.8)
    total = pred_count * unit / effective_cpus
    # Add fixed ~5s overhead for pool startup + UniProt batch + summary write
    total += 5
    return _format_duration(total)


def show_local_execution(project_path: str, af3_folder: str,
                          pred_count: int, dominant_format: str):
    """Display local execution options."""

    st.subheader("🖥️ Local Execution")

    cache_file = Path(af3_folder) / "af3_app_all_models_analysis.json"

    col1, col2 = st.columns(2)

    with col1:
        num_cpus = st.slider("Number of CPUs:", 1, 32, 8)
        pae_cutoff = st.slider("PAE cutoff (Å):", 5.0, 15.0, 10.0, 0.5)

    with col2:
        analyze_all = st.checkbox("Analyze all 5 models per prediction", value=False,
                                   help="Default: top-ranked model only. Check to analyze all seed/sample models.")

        skip_cached = False
        if cache_file.exists():
            skip_cached = st.checkbox(
                "Skip predictions already in cache",
                value=False,
                help=(
                    "Analyze only predictions not yet in the cache — useful when "
                    "new predictions were added to a folder that was already run. "
                    "Leave unchecked if you changed the PAE cutoff or the all-models "
                    "option, otherwise the cache will mix results from different settings."
                ),
            )

    # Pre-run summary: how many predictions, how long it'll take
    n_to_run = pred_count
    if skip_cached and cache_file.exists():
        try:
            import json as _json
            with open(cache_file) as _f:
                _cached_names = {r.get('prediction_name') for r in _json.load(_f) if r.get('prediction_name')}
            from core.scanner import AF3Scanner
            _all = AF3Scanner(Path(af3_folder)).scan()
            n_to_run = sum(1 for p in _all if p['name'] not in _cached_names)
        except Exception:
            n_to_run = pred_count

    eta = _estimate_local_runtime(n_to_run, num_cpus, analyze_all, dominant_format)

    sum_col1, sum_col2, sum_col3 = st.columns(3)
    sum_col1.metric("Predictions to analyze", n_to_run)
    sum_col2.metric("Models per prediction", "all 5" if analyze_all else "top only")
    sum_col3.metric("Estimated runtime", eta,
                    help=(
                        f"Rough estimate: ~{_PER_PRED_SECS.get((dominant_format, analyze_all), 1.0):.1f}s "
                        f"per prediction (single-CPU), divided by {num_cpus} CPUs × 0.8 efficiency, "
                        "plus ~5s fixed overhead. Real time varies with chain size and disk speed."
                    ))

    st.divider()

    # Execution button — label reflects what's about to happen
    if skip_cached:
        btn_label = f"Analyze {n_to_run} Missing Prediction(s)"
    elif cache_file.exists():
        btn_label = f"Re-analyze All {n_to_run} Predictions"
    else:
        btn_label = f"Start Local Analysis ({n_to_run} predictions)"
    if st.button(btn_label, type="primary", disabled=(n_to_run == 0)):
        run_local_analysis(project_path, af3_folder, num_cpus, pae_cutoff,
                           analyze_all, skip_cached)


def write_summary_txt(results: list, out_path: Path, pae_cutoff: float, analyze_all: bool,
                      gene_cache: dict = None):
    """Write a human-readable summary .txt in the style of AF3_PD_analysis_v4_summary.txt.

    gene_cache: pre-built dict mapping ACCESSION -> gene name (or accession if unknown).
                If None, no gene name annotation is added.
    """
    from datetime import datetime

    gene_cache = gene_cache or {}

    def gene_label(acc):
        acc = acc.upper()
        gene = gene_cache.get(acc, acc)
        return f"{gene} ({acc})" if gene != acc else acc

    # Collapse to best model per prediction (highest ipSAE)
    best = {}
    for r in results:
        name = r['prediction_name']
        ipsae = r.get('ipsae') or 0
        if name not in best or ipsae > (best[name].get('ipsae') or 0):
            best[name] = r
    per_pred = list(best.values())
    per_pred.sort(key=lambda r: r.get('ipsae') or 0, reverse=True)

    high   = [r for r in per_pred if (r.get('ipsae') or 0) >  0.7]
    medium = [r for r in per_pred if 0.5 < (r.get('ipsae') or 0) <= 0.7]
    low    = [r for r in per_pred if 0.3 < (r.get('ipsae') or 0) <= 0.5]
    vlow   = [r for r in per_pred if (r.get('ipsae') or 0) <= 0.3]

    total = len(per_pred)
    models_label = "all models" if analyze_all else "top-ranked model only"

    lines = []
    lines.append("AF3 Pulldown Analysis - App Summary")
    lines.append("=" * 60)
    lines.append(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"ipSAE PAE Cutoff: {pae_cutoff}Å")
    lines.append(f"Models analyzed: {models_label}")
    lines.append(f"Total models in cache: {len(results)}")
    lines.append("")
    lines.append("ipSAE Confidence Classification:")
    lines.append("  High      (ipSAE > 0.7)  : strong evidence of interaction")
    lines.append("  Medium    (ipSAE 0.5-0.7): promising, likely genuine")
    lines.append("  Low       (ipSAE 0.3-0.5): weak signal, inspect manually")
    lines.append("  Very Low  (ipSAE < 0.3)  : likely not an interaction")
    lines.append("")
    lines.append(f"Total Predictions: {total}")
    lines.append(f"  High:      {len(high):4d}  ({100*len(high)/total:.1f}%)" if total else "  High:         0")
    lines.append(f"  Medium:    {len(medium):4d}  ({100*len(medium)/total:.1f}%)" if total else "  Medium:       0")
    lines.append(f"  Low:       {len(low):4d}  ({100*len(low)/total:.1f}%)" if total else "  Low:          0")
    lines.append(f"  Very Low:  {len(vlow):4d}  ({100*len(vlow)/total:.1f}%)" if total else "  Very Low:     0")
    lines.append("")

    name_w = 44
    hdr = f"{'Bait':<22} {'Prey':<22} {'ipSAE':>6}  {'iPTM':>5}  {'PAE<=3':>6}  {'PAE<=5':>6}  {'PAE<=8':>6}  {'iPLDDT':>6}"
    sep = "-" * len(hdr)

    for tier_label, tier_results in [
        (f"HIGH CONFIDENCE (ipSAE > 0.7) [N={len(high)}]", high),
        (f"MEDIUM CONFIDENCE (ipSAE 0.5-0.7) [N={len(medium)}]", medium),
        (f"LOW / AMBIGUOUS (ipSAE 0.3-0.5) [N={len(low)}]", low),
    ]:
        lines.append("=" * 80)
        lines.append(tier_label)
        lines.append("=" * 80)
        lines.append(hdr)
        lines.append(sep)
        if not tier_results:
            lines.append("  (none)")
        for r in tier_results:
            pred_name = r['prediction_name']
            from core.utils import split_prediction_name as _split
            bait_acc, prey_acc = _split(pred_name)
            bait_col = gene_label(bait_acc)[:22]
            prey_col = gene_label(prey_acc)[:22]
            ipsae  = f"{r.get('ipsae') or 0:.3f}"
            iptm   = f"{r.get('iptm') or 0:.3f}"
            pae3   = str(r.get('contacts_pae3') or 0)
            pae5   = str(r.get('contacts_pae5') or 0)
            pae8   = str(r.get('contacts_pae8') or 0)
            iplddt = f"{r.get('interface_plddt') or 0:.1f}"
            lines.append(f"{bait_col:<22} {prey_col:<22} {ipsae:>6}  {iptm:>5}  {pae3:>6}  {pae5:>6}  {pae8:>6}  {iplddt:>6}")
        lines.append("")

    with open(out_path, 'w') as f:
        f.write("\n".join(lines))


def run_local_analysis(project_path: str, af3_folder: str, num_cpus: int, pae_cutoff: float,
                       analyze_all: bool = True, skip_cached: bool = False):
    """Run analysis locally with multiprocessing."""

    progress_bar = st.progress(0)
    status_text = st.empty()
    log_area = st.empty()
    status_container = st.status("Analysis in progress...", expanded=True)

    log_lines = []
    results = []

    try:
        # Import analyzer
        from core.scanner import AF3Scanner, resolve_prediction_dir
        from core.analyzer import analyze_prediction_all_models
        import json

        # Get all prediction directories (resolve AF3 Server nesting)
        all_pred_dirs = []
        for d in os.listdir(af3_folder):
            d_path = Path(af3_folder) / d
            if d_path.is_dir() and not d.startswith('seed-') and not d.startswith('.'):
                all_pred_dirs.append(resolve_prediction_dir(d_path))

        if not all_pred_dirs:
            st.warning("No prediction directories found.")
            return

        # Optionally filter out predictions that are already cached
        cache_file = Path(af3_folder) / "af3_app_all_models_analysis.json"
        cached_results = []
        pred_dirs = all_pred_dirs
        if skip_cached and cache_file.exists():
            try:
                with open(cache_file) as f:
                    cached_results = json.load(f)
                cached_names = {r.get('prediction_name') for r in cached_results if r.get('prediction_name')}
                pred_dirs = [d for d in all_pred_dirs if d.name not in cached_names]
                skipped = len(all_pred_dirs) - len(pred_dirs)
                st.info(f"Skipping {skipped} predictions already in cache; {len(pred_dirs)} to analyze.")
            except Exception as e:
                st.warning(f"Could not read existing cache ({e}); running full analysis instead.")
                cached_results = []
                pred_dirs = all_pred_dirs

        total = len(pred_dirs)
        if total == 0:
            status_container.update(label="Nothing to do", state="complete")
            st.success("All predictions are already in the cache. Go to **3. Results** to view them.")
            return

        status_text.text(f"Starting analysis with {num_cpus} CPUs on {total} predictions...")

        # Collect all unique accessions from ALL directory names (so gene lookup covers the whole cache)
        from core.utils import fetch_gene_names_batch, split_prediction_name

        unique_accs = sorted({
            acc
            for d in all_pred_dirs
            for acc in split_prediction_name(d.name)
            if acc
        })

        # Use multiprocessing pool with imap_unordered for streaming results
        from multiprocessing import Pool, cpu_count
        from functools import partial

        # Cap CPUs at available cores
        actual_cpus = min(num_cpus, cpu_count())
        status_text.text(f"Using {actual_cpus} CPUs...")

        analyze_fn = partial(analyze_prediction_all_models,
                             ipsae_pae_cutoff=pae_cutoff,
                             top_only=not analyze_all)

        all_results = []
        try:
            with Pool(actual_cpus) as pool:
                for i, result_list in enumerate(pool.imap_unordered(
                    analyze_fn,
                    pred_dirs
                )):
                    if result_list:
                        all_results.extend(result_list)

                    progress = (i + 1) / total
                    progress_bar.progress(progress)
                    pred_label = result_list[0]['prediction_name'] if result_list else '(no results)'
                    status_text.text(f"Analyzed {i+1}/{total} predictions ({len(all_results)} models)")

                    with status_container:
                        st.info(f"Completed: {pred_label}")
        except (BrokenPipeError, ConnectionResetError, EOFError, OSError):
            # Worker pipe broke (e.g. Streamlit rerun) — use whatever results came in
            st.warning(f"Pool interrupted — saving {len(all_results)} models collected so far.")

        # Strip sequences from new results before saving (huge and already in data.json)
        for r in all_results:
            if 'sequences' in r:
                del r['sequences']

        # Per-prediction JSON files: only rewrite for predictions we just re-analyzed
        from collections import defaultdict
        by_pred = defaultdict(list)
        for r in all_results:
            by_pred[r.get('prediction_name', '')].append(r)
        pred_dir_map = {d.name: d for d in pred_dirs}
        for pred_name, entries in by_pred.items():
            if pred_name:
                # Fallback path also resolves AF3 Server <name>/<name>/ nesting
                target_dir = pred_dir_map.get(
                    pred_name,
                    resolve_prediction_dir(Path(af3_folder) / pred_name),
                )
                pred_json = target_dir / "af3_app_analysis.json"
                try:
                    with open(pred_json, 'w') as pf:
                        json.dump(entries, pf)  # no indent — 3× faster on networked FS
                except OSError as e:
                    st.warning(f"Could not save {pred_json.name}: {e}")

        # Merge freshly analyzed results with any untouched cached entries
        results = cached_results + all_results

        with open(cache_file, 'w') as f:
            json.dump(results, f)

        # Fetch gene names via UniProt batch API with live progress
        n_batches = max(1, (len(unique_accs) + 99) // 100)
        status_text.text(f"Fetching gene names for {len(unique_accs)} proteins (0/{n_batches} batches)...")
        progress_bar.progress(0.0)

        gene_cache = {}

        def _gene_progress(done, total):
            batch_num = (done + 99) // 100
            status_text.text(f"Fetching gene names: {done}/{total} proteins ({batch_num}/{n_batches} batches)...")
            progress_bar.progress(min(done / total, 1.0))

        gene_cache = fetch_gene_names_batch(unique_accs, progress_callback=_gene_progress)

        # Push gene names into session_state cache for immediate use in tables
        st.session_state.setdefault('gene_name_cache', {}).update(gene_cache)

        # Write human-readable summary
        txt_file = Path(af3_folder) / "af3_app_analysis_summary.txt"
        write_summary_txt(results, txt_file, pae_cutoff, analyze_all, gene_cache)

        status_container.update(label="Analysis complete!", state="complete")
        status_text.text("Analysis complete!")
        progress_bar.progress(1.0)

        st.success(f"Results saved to {cache_file}")
        st.success(f"Summary written to {txt_file}")
        if cached_results:
            st.success(
                f"Newly analyzed: {len(all_results)} models. "
                f"Cache now contains {len(results)} models total "
                f"({len(cached_results)} kept from previous run)."
            )
        else:
            st.success(f"Total models analyzed: {len(results)}")

    except Exception as e:
        st.error(f"Analysis failed: {e}")
        import traceback
        log_area.text(traceback.format_exc())


def _scan_slurm_state(af3_folder: str) -> dict:
    """Folder-driven view of SLURM run state. Survives browser refresh.

    For each chunk file ``_slurm_chunk_<i>.txt`` we check the matching
    ``_slurm_results_<i>.json``. A result file with mtime ≥ its chunk file
    is considered ``fresh`` (i.e. produced by this run, not left over from
    a previous submission). An empty fresh result (``[]``, 2 bytes) means
    the job ran but the analyzer found no models — usually wrong input format.
    """
    import time
    af3_path = Path(af3_folder)
    chunks = sorted(af3_path.glob('_slurm_chunk_*.txt'))
    cache_file = af3_path / 'af3_app_all_models_analysis.json'

    chunk_info = []
    for chunk_file in chunks:
        idx = int(chunk_file.stem.rsplit('_', 1)[1])
        result_file = af3_path / f"_slurm_results_{idx}.json"
        chunk_mtime = chunk_file.stat().st_mtime
        if result_file.exists():
            r_mtime = result_file.stat().st_mtime
            r_size = result_file.stat().st_size
            fresh = r_mtime >= chunk_mtime
            empty = r_size <= 2  # "[]"
        else:
            r_mtime = None
            r_size = None
            fresh = False
            empty = False
        chunk_info.append({
            'idx': idx,
            'chunk_file': chunk_file,
            'result_file': result_file if result_file.exists() else None,
            'chunk_mtime': chunk_mtime,
            'result_mtime': r_mtime,
            'result_size': r_size,
            'fresh': fresh,
            'empty': empty,
        })

    expected = len(chunks)
    completed = sum(1 for c in chunk_info if c['fresh'])
    empty_count = sum(1 for c in chunk_info if c['fresh'] and c['empty'])

    cache_mtime = cache_file.stat().st_mtime if cache_file.exists() else 0
    latest_result_mtime = max((c['result_mtime'] for c in chunk_info if c['result_mtime']), default=0)
    cache_is_current = cache_file.exists() and cache_mtime >= latest_result_mtime

    return {
        'expected_chunks': expected,
        'completed_chunks': completed,
        'empty_chunks': empty_count,
        'chunk_info': chunk_info,
        'cache_file': cache_file,
        'cache_exists': cache_file.exists(),
        'cache_size': cache_file.stat().st_size if cache_file.exists() else 0,
        'cache_mtime': cache_mtime,
        'cache_is_current': cache_is_current,
        'now': time.time(),
    }


def _humanize_ago(seconds: float) -> str:
    if seconds < 60:
        return f"{int(seconds)} s ago"
    if seconds < 3600:
        return f"{int(seconds/60)} min ago"
    return f"{seconds/3600:.1f} h ago"


def _render_slurm_state_panel(af3_folder: str):
    """The SLURM-run status panel (folder-driven, no auto-refresh).

    User clicks 🔄 Refresh once their jobs finish (typically 1–5 min).
    The auto-merge fires the moment a render shows all chunks complete.
    """
    state = _scan_slurm_state(af3_folder)
    expected = state['expected_chunks']

    if expected == 0 and not state['cache_exists']:
        return  # nothing to show

    st.subheader("SLURM run state")

    if expected > 0:
        col1, col2, col3 = st.columns(3)
        col1.metric("Chunks submitted", expected)
        col2.metric("Chunks finished",
                    f"{state['completed_chunks']} / {expected}")
        col3.metric("Empty / failed",
                    state['empty_chunks'])

        # Per-chunk detail
        with st.expander("Per-chunk detail", expanded=False):
            for c in state['chunk_info']:
                if c['fresh']:
                    age = _humanize_ago(state['now'] - c['result_mtime'])
                    if c['empty']:
                        st.text(f"  ⚠ chunk {c['idx']}: empty result ({age})")
                    else:
                        kb = c['result_size'] / 1024
                        st.text(f"  ✅ chunk {c['idx']}: {kb:.1f} KB ({age})")
                elif c['result_file']:
                    age = _humanize_ago(state['now'] - c['result_mtime'])
                    st.text(f"  🟡 chunk {c['idx']}: stale result from previous run ({age}) — current job still pending")
                else:
                    age = _humanize_ago(state['now'] - c['chunk_mtime'])
                    st.text(f"  ⏳ chunk {c['idx']}: submitted {age}, no result yet")

    pending = expected - state['completed_chunks']

    # Decide which prompt to show
    if pending > 0:
        st.info(
            f"⏳ **{pending} of {expected} chunks still running on the cluster.**  \n"
            f"Click 🔄 Refresh below once your jobs are done (typically 1–5 minutes). "
            f"When all chunks are complete, the next render will auto-merge into the analysis cache "
            f"and show the progress here."
        )
        st.code("squeue -u $USER | grep AF3app", language="bash")
        if st.button("🔄 Refresh", type="primary"):
            st.rerun()
    elif expected > 0 and not state['cache_is_current']:
        # All chunks present but cache not yet built (or stale).
        if state['empty_chunks'] == expected:
            st.error(
                f"❌ All {expected} chunks finished, but every result file is empty. "
                "The analyzer found no models — usually this means the input folder isn't "
                "AF3 or AlphaPulldown format. Nothing to merge."
            )
            # Even so, offer a cleanup button so scratch files don't linger
            if st.button("🗑 Clean up empty chunk files"):
                _cleanup_slurm_scratch(af3_folder)
                st.rerun()
        else:
            non_empty = expected - state['empty_chunks']
            st.success(f"✅ All {expected} chunks complete — {non_empty} with results.")
            if state['empty_chunks'] > 0:
                st.warning(f"{state['empty_chunks']} chunk(s) returned no models — partial merge.")

            merge_marker = st.session_state.get('_slurm_merge_attempted_for')
            already_attempted = (merge_marker == af3_folder)

            # Manual merge button — always visible as a fallback. The auto-merge
            # path below will also fire on first render, but if anything went
            # wrong (e.g. Streamlit didn't pick up code changes, browser didn't
            # reload, auto-merge failed silently), this button still works.
            col_btn, col_status = st.columns([1, 2])
            with col_btn:
                manual_clicked = st.button(
                    "📦 Merge chunks into cache",
                    type="primary",
                    help="Combine the chunk result files into af3_app_all_models_analysis.json. "
                         "Safe to click even if auto-merge already ran.",
                )
            with col_status:
                if already_attempted:
                    st.caption("Auto-merge was attempted in this session — click to retry.")
                else:
                    st.caption("Auto-merge will fire below; you can also trigger it manually here.")

            if manual_clicked:
                st.session_state['_slurm_merge_attempted_for'] = af3_folder
                with st.status(
                    f"Merging {non_empty} chunk(s) into analysis cache...",
                    expanded=True,
                ) as status:
                    try:
                        merge_slurm_results(af3_folder, expected)
                        status.update(label=f"✅ Merge complete — {non_empty} chunks combined.",
                                      state="complete")
                    except Exception as e:
                        status.update(label=f"❌ Merge failed: {e}", state="error")
                st.rerun(scope="app")  # whole-app rerun so post-merge state shows everywhere

            # Auto-merge — only fires once per folder per session. On error
            # the marker stays set so we don't infinite-loop; user can click
            # the manual button above to retry explicitly.
            elif not already_attempted:
                st.session_state['_slurm_merge_attempted_for'] = af3_folder
                with st.status(
                    f"Auto-merging {non_empty} chunk(s) into analysis cache...",
                    expanded=True,
                ) as status:
                    try:
                        merge_slurm_results(af3_folder, expected)
                        status.update(label=f"✅ Merge complete — {non_empty} chunks combined.",
                                      state="complete")
                    except Exception as e:
                        status.update(label=f"❌ Auto-merge failed: {e} — try the manual button above.",
                                      state="error")
                st.rerun(scope="app")  # whole-app rerun so post-merge state shows everywhere
    elif state['cache_is_current'] and expected > 0:
        # Cache built from this run — happy path
        age = _humanize_ago(state['now'] - state['cache_mtime'])
        kb = state['cache_size'] / 1024
        st.success(
            f"✅ Analysis cache built ({kb:.1f} KB, {age}). "
            "Go to **3. Results** to view."
        )
        col_a, col_b = st.columns([1, 1])
        with col_a:
            if st.button("🗑 Clean up SLURM scratch files"):
                _cleanup_slurm_scratch(af3_folder)
                st.rerun()

    elif state['cache_exists']:
        # Cache exists, no current SLURM run
        age = _humanize_ago(state['now'] - state['cache_mtime'])
        kb = state['cache_size'] / 1024
        st.caption(f"Existing analysis cache: {kb:.1f} KB ({age}). "
                   "Submit a new SLURM run below to recompute, or skip to **3. Results**.")

    st.divider()


def _cleanup_slurm_scratch(af3_folder: str):
    """Delete _slurm_chunk_*.txt, _slurm_results_*.json, and SLURM .log/.err files."""
    af3_path = Path(af3_folder)
    n = 0
    for pattern in ['_slurm_chunk_*.txt', '_slurm_results_*.json',
                    'AF3_app_chunk*.log', 'AF3_app_chunk*.err']:
        for f in af3_path.glob(pattern):
            try:
                f.unlink()
                n += 1
            except OSError:
                pass
    st.success(f"Removed {n} SLURM scratch file(s).")


def show_slurm_execution(project_path: str, af3_folder: str, pred_count: int):
    """Display SLURM cluster submission options."""

    # Always-visible folder-driven run-state panel (works without session_state)
    _render_slurm_state_panel(af3_folder)

    st.markdown("""
    Submit analysis to the CPU cluster (vader nodes). The predictions will be
    split into chunks and distributed across multiple SLURM jobs.
    """)

    col1, col2 = st.columns(2)

    with col1:
        num_jobs = st.number_input("Number of SLURM jobs:", min_value=1, max_value=50,
                                   value=6,
                                   help="Predictions will be split evenly across this many jobs (6 = one per vader node)")
        cpus_per_job = st.selectbox("CPUs per job:", [32, 48, 64, 72], index=2)
        memory = st.selectbox("Memory per job:", ["64G", "128G", "256G"], index=1)

    with col2:
        pae_cutoff = st.slider("PAE cutoff (Å):", 5.0, 15.0, 10.0, 0.5, key="slurm_pae")
        analyze_all = st.checkbox("Analyze all 5 models per prediction", value=False, key="slurm_all_models",
                                   help="Default: top-ranked model only.")
        qos = "normal"

    preds_per_job = (pred_count + num_jobs - 1) // num_jobs if num_jobs > 0 else pred_count
    st.caption(f"~{preds_per_job} predictions per job, {cpus_per_job} CPUs each")

    st.divider()

    # Per-job squeue/sacct detail — only when we have job IDs from this session
    if st.session_state.get('slurm_job_ids'):
        from core.slurm_manager import check_job_status

        with st.expander("Per-job cluster status (squeue / sacct)", expanded=False):
            for job_id in st.session_state['slurm_job_ids']:
                status = check_job_status(job_id)
                s = status['status']
                icon = {"RUNNING": "🟢", "PENDING": "🟡", "COMPLETED": "✅",
                        "FAILED": "❌"}.get(s, "⚪")
                st.text(f"  {icon} Job {job_id}: {s}")

    # Submit button
    if st.button("Submit to SLURM", type="primary"):
        submit_slurm_jobs(af3_folder, num_jobs, cpus_per_job, memory, qos,
                          pae_cutoff, analyze_all)


def submit_slurm_jobs(af3_folder: str, num_jobs: int, cpus_per_job: int,
                       memory: str, qos: str, pae_cutoff: float, analyze_all: bool):
    """Split predictions into chunks and submit SLURM jobs."""
    import json

    from core.scanner import resolve_prediction_dir

    app_dir = str(Path(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    af3_path = Path(af3_folder)

    # List all prediction directories (resolve AF3 Server nesting)
    pred_dirs = sorted([
        str(resolve_prediction_dir(d)) for d in af3_path.iterdir()
        if d.is_dir() and not d.name.startswith('seed-') and not d.name.startswith('.')
    ])
    total = len(pred_dirs)

    if total == 0:
        st.error("No prediction directories found.")
        return

    # Split into chunks
    chunk_size = (total + num_jobs - 1) // num_jobs
    chunks = [pred_dirs[i:i + chunk_size] for i in range(0, total, chunk_size)]

    # Remove orphan SLURM scratch from any previous run so the new run's
    # mtime-based completion detection doesn't see stale [] results as "done".
    n_removed = 0
    for pattern in ['_slurm_chunk_*.txt', '_slurm_results_*.json',
                    'AF3_app_chunk*.log', 'AF3_app_chunk*.err']:
        for old in af3_path.glob(pattern):
            try:
                old.unlink()
                n_removed += 1
            except OSError:
                pass
    if n_removed > 0:
        st.info(f"Cleared {n_removed} SLURM scratch file(s) from a previous run.")

    # Write chunk lists and generate SLURM scripts
    from core.slurm_manager import submit_job

    partition = "cpu"
    venv_activate = f"{app_dir}/.venv/bin/activate"

    job_ids = []
    for i, chunk in enumerate(chunks):
        # Write chunk file listing the prediction folder names
        chunk_file = af3_path / f"_slurm_chunk_{i}.txt"
        with open(chunk_file, 'w') as f:
            f.write('\n'.join(chunk))

        output_file = af3_path / f"_slurm_results_{i}.json"
        job_name = f"AF3app_{Path(af3_folder).name}_chunk{i}"

        # Escape paths for safe embedding in shell single-quoted strings
        esc_chunk = str(chunk_file).replace("'", "'\\''")
        esc_af3 = str(af3_folder).replace("'", "'\\''")
        esc_output = str(output_file).replace("'", "'\\''")
        esc_venv = venv_activate.replace("'", "'\\''")
        esc_app = app_dir.replace("'", "'\\''")

        script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpus_per_job}
#SBATCH --mem={memory}
#SBATCH --time=23:59:59
#SBATCH --qos={qos}
#SBATCH --partition={partition}
#SBATCH --output='{esc_af3}/AF3_app_chunk{i}_%j.log'
#SBATCH --error='{esc_af3}/AF3_app_chunk{i}_%j.err'

echo "AF3 Analysis App - Chunk {i+1}/{len(chunks)}"
echo "Job ID: $SLURM_JOB_ID | Node: $SLURMD_NODENAME | CPUs: $SLURM_CPUS_PER_TASK"
echo "Predictions in chunk: {len(chunk)}"
echo "Start: $(date)"

source '{esc_venv}'
cd '{esc_app}'

# Pass paths via environment variables to avoid shell/Python quoting issues
export AF3_CHUNK_FILE='{esc_chunk}'
export AF3_FOLDER='{esc_af3}'
export AF3_OUTPUT_FILE='{esc_output}'

python << 'PYEOF'
import json, sys, os
from pathlib import Path
from multiprocessing import Pool
from functools import partial
from core.analyzer import analyze_prediction_all_models

chunk_file = os.environ["AF3_CHUNK_FILE"]
af3_folder = os.environ["AF3_FOLDER"]
output_file = os.environ["AF3_OUTPUT_FILE"]

with open(chunk_file) as f:
    pred_paths = [line.strip() for line in f if line.strip()]

pred_dirs = [Path(p) for p in pred_paths]
print(f"Analyzing {{len(pred_dirs)}} predictions with {cpus_per_job} CPUs...")

analyze_fn = partial(analyze_prediction_all_models,
                     ipsae_pae_cutoff={pae_cutoff},
                     top_only={not analyze_all})

all_results = []
with Pool({cpus_per_job}) as pool:
    for i, result_list in enumerate(pool.imap_unordered(analyze_fn, pred_dirs)):
        if result_list:
            for r in result_list:
                r.pop("sequences", None)
            all_results.extend(result_list)
            pn = result_list[0].get("prediction_name", "")
            if pn:
                pdir = next((d for d in pred_dirs if d.name == pn), None)
                if pdir:
                    try:
                        with open(pdir / "af3_app_analysis.json", "w") as pf:
                            json.dump(result_list, pf, indent=2)
                    except Exception:
                        pass
        if (i+1) % 100 == 0:
            print(f"  {{i+1}}/{{len(pred_dirs)}} done ({{len(all_results)}} models)")

with open(output_file, "w") as f:
    json.dump(all_results, f)

print(f"Done: {{len(all_results)}} models saved to {{output_file}}")
PYEOF

echo "Finished: $(date)"
"""
        result = submit_job(script, af3_folder)
        if result['success']:
            job_ids.append(result['job_id'])
            st.success(f"Chunk {i+1}/{len(chunks)}: submitted job {result['job_id']} ({len(chunk)} predictions)")
        else:
            st.error(f"Chunk {i+1} failed: {result['error']}")

    if job_ids:
        st.session_state['slurm_job_ids'] = job_ids
        st.session_state['slurm_num_chunks'] = len(chunks)
        # Clear any stale merge marker from previous SLURM runs in this session
        st.session_state.pop('_slurm_merge_attempted_for', None)
        st.info(f"Submitted {len(job_ids)} jobs. Check progress with the command below.")
        st.code(f"squeue -u $USER | grep AF3app", language="bash")
        # Re-render the page so the folder-state panel above picks up the
        # newly-written chunk files and starts the 30 s auto-refresh.
        # Without this rerun the panel stays on its pre-submit state and
        # never injects the auto-refresh JS, leaving the user stuck.
        st.rerun()


def merge_slurm_results(af3_folder: str, num_chunks: int):
    """Merge chunk result files into the final analysis cache."""
    import json

    af3_path = Path(af3_folder)
    all_results = []

    for i in range(num_chunks):
        chunk_file = af3_path / f"_slurm_results_{i}.json"
        if chunk_file.exists():
            with open(chunk_file, 'r') as f:
                all_results.extend(json.load(f))
            st.text(f"  Loaded chunk {i}: {chunk_file.name}")

    if not all_results:
        st.error("No chunk results found. Check job logs for errors.")
        return

    # Save merged results (no indent — 3× smaller, faster on networked FS)
    cache_file = af3_path / "af3_app_all_models_analysis.json"
    with open(cache_file, 'w') as f:
        json.dump(all_results, f)

    st.text(f"  Saved {len(all_results)} models to {cache_file.name}")

    # Save per-prediction JSON files for fast lookup in Detailed Analysis.
    # On networked storage these many-tiny-file writes are I/O-bound; we
    # skip the slow resolve_prediction_dir() call when no AF3 Server-style
    # nesting exists in this folder, and we batch-update progress.
    from collections import defaultdict
    from core.scanner import resolve_prediction_dir
    by_pred = defaultdict(list)
    for r in all_results:
        by_pred[r.get('prediction_name', '')].append(r)

    pred_names = [k for k in by_pred if k]
    total = len(pred_names)

    # Detect AF3 Server <name>/<name>/ nesting once; if absent we can skip
    # the per-iteration iterdir() call inside resolve_prediction_dir.
    needs_resolve = any(
        (af3_path / p).is_dir() and (af3_path / p / p).is_dir()
        for p in pred_names[:5]   # cheap sample
    )

    progress = st.progress(0, text=f"Writing per-prediction files... 0/{total}")
    n_pred_files = 0
    update_every = max(1, total // 50)  # ~50 progress updates total
    for idx, pred_name in enumerate(pred_names):
        if needs_resolve:
            target_dir = resolve_prediction_dir(af3_path / pred_name)
        else:
            target_dir = af3_path / pred_name
        pred_json = target_dir / "af3_app_analysis.json"
        try:
            with open(pred_json, 'w') as pf:
                json.dump(by_pred[pred_name], pf)  # no indent — faster
            n_pred_files += 1
        except Exception:
            pass
        if idx % update_every == 0 or idx == total - 1:
            progress.progress((idx + 1) / total,
                              text=f"Writing per-prediction files... {idx + 1}/{total}")
    progress.empty()

    st.text(f"  Saved {n_pred_files} per-prediction JSON files")

    # Clean up chunk files
    for i in range(num_chunks):
        for pattern in [f"_slurm_chunk_{i}.txt", f"_slurm_results_{i}.json"]:
            f = af3_path / pattern
            f.unlink(missing_ok=True)
