"""
Step 2: Analyze - Local Multiprocessing Execution

Runs full analysis on all predictions using local CPUs.
"""

import streamlit as st
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


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

    # Count predictions
    pred_count = len([d for d in os.listdir(af3_folder)
                     if os.path.isdir(os.path.join(af3_folder, d))
                     and not d.startswith('seed-') and not d.startswith('.')])

    st.info(f"Found {pred_count} predictions to analyze")

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
        show_local_execution(project_path, af3_folder)

    with tab_slurm:
        show_slurm_execution(project_path, af3_folder, pred_count)


def show_local_execution(project_path: str, af3_folder: str):
    """Display local execution options."""

    st.subheader("🖥️ Local Execution")

    col1, col2 = st.columns(2)

    with col1:
        num_cpus = st.slider("Number of CPUs:", 1, 32, 8)
        pae_cutoff = st.slider("PAE cutoff (Å):", 5.0, 15.0, 10.0, 0.5)

    with col2:
        analyze_all = st.checkbox("Analyze all 5 models per prediction", value=False,
                                   help="Default: top-ranked model only. Check to analyze all seed/sample models.")

    st.divider()

    # Execution button — label changes if results already exist
    cache_file = Path(af3_folder) / "af3_app_all_models_analysis.json"
    btn_label = "Re-analyze All Predictions" if cache_file.exists() else "Start Local Analysis"
    if st.button(btn_label, type="primary"):
        run_local_analysis(project_path, af3_folder, num_cpus, pae_cutoff, analyze_all)


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
            if '_and_' in pred_name:
                bait_acc, prey_acc = pred_name.split('_and_', 1)
            else:
                bait_acc, prey_acc = pred_name, ''
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
                       analyze_all: bool = True):
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

        # Get all prediction directories (resolve AF3 Server nesting)
        pred_dirs = []
        for d in os.listdir(af3_folder):
            d_path = Path(af3_folder) / d
            if d_path.is_dir() and not d.startswith('seed-') and not d.startswith('.'):
                pred_dirs.append(resolve_prediction_dir(d_path))

        total = len(pred_dirs)

        if total == 0:
            st.warning("No prediction directories found.")
            return

        status_text.text(f"Starting analysis with {num_cpus} CPUs on {total} predictions...")

        # Collect all unique accessions from directory names
        from core.utils import fetch_gene_names_batch

        unique_accs = sorted({
            acc.upper()
            for d in pred_dirs
            for part in (d.name.split('_and_', 1) if '_and_' in d.name else [d.name])
            for acc in [part.upper()]
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
        except BrokenPipeError:
            # Worker pipe broke (e.g. Streamlit rerun) — use whatever results came in
            st.warning(f"Pool interrupted — saving {len(all_results)} models collected so far.")

        results = all_results

        # Strip sequences from results before saving (they're huge and already in data.json)
        for r in results:
            if 'sequences' in r:
                del r['sequences']

        # Save per-prediction JSON files for fast lookup in detailed analysis
        import json
        from collections import defaultdict
        by_pred = defaultdict(list)
        for r in results:
            by_pred[r.get('prediction_name', '')].append(r)
        # Build name -> resolved dir lookup
        pred_dir_map = {d.name: d for d in pred_dirs}
        for pred_name, entries in by_pred.items():
            if pred_name:
                target_dir = pred_dir_map.get(pred_name, Path(af3_folder) / pred_name)
                pred_json = target_dir / "af3_app_analysis.json"
                try:
                    with open(pred_json, 'w') as pf:
                        json.dump(entries, pf, indent=2)
                except Exception:
                    pass

        # Save combined results (big JSON)
        cache_file = Path(af3_folder) / "af3_app_all_models_analysis.json"

        with open(cache_file, 'w') as f:
            json.dump(results, f, indent=2)

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
        st.success(f"Total models analyzed: {len(results)}")

    except Exception as e:
        st.error(f"Analysis failed: {e}")
        import traceback
        log_area.text(traceback.format_exc())


def show_slurm_execution(project_path: str, af3_folder: str, pred_count: int):
    """Display SLURM cluster submission options."""

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

    # ── SLURM Job Monitor (always visible when jobs exist) ──
    has_jobs = 'slurm_job_ids' in st.session_state and st.session_state['slurm_job_ids']

    if has_jobs:
        from core.slurm_manager import check_job_status

        st.subheader("SLURM Job Monitor")

        n_running = 0
        n_pending = 0
        n_done = 0
        n_failed = 0
        job_statuses = []
        for job_id in st.session_state['slurm_job_ids']:
            status = check_job_status(job_id)
            s = status['status']
            job_statuses.append((job_id, s))
            if s == 'RUNNING':
                n_running += 1
            elif s == 'PENDING':
                n_pending += 1
            elif s == 'FAILED':
                n_failed += 1
            else:
                n_done += 1

        total_jobs = len(st.session_state['slurm_job_ids'])
        all_done = (n_running == 0 and n_pending == 0)

        # Summary bar
        col1, col2, col3, col4 = st.columns(4)
        col1.metric("Running", n_running)
        col2.metric("Pending", n_pending)
        col3.metric("Completed", n_done)
        col4.metric("Failed", n_failed)

        # Per-job status
        with st.expander("Job details", expanded=False):
            for job_id, s in job_statuses:
                icon = {"RUNNING": "🟢", "PENDING": "🟡", "COMPLETED": "✅",
                        "FAILED": "❌"}.get(s, "⚪")
                st.text(f"  {icon} Job {job_id}: {s}")

        if not all_done:
            st.info(f"{n_running + n_pending} job(s) still running/pending.")
            if st.button("Refresh Status", type="primary"):
                st.rerun()
        else:
            # All done — auto-merge
            if n_failed > 0:
                st.warning(f"{n_failed} job(s) failed. Partial results will be merged.")

            merge_key = f"merged_{','.join(st.session_state['slurm_job_ids'])}"
            if merge_key not in st.session_state:
                with st.spinner("Merging results..."):
                    merge_slurm_results(af3_folder, st.session_state.get('slurm_num_chunks', num_jobs))
                    st.session_state[merge_key] = True

            st.success(f"All {total_jobs} jobs finished. Results merged and ready.")
            st.info("Go to **3. Results** to view the analysis.")

            if st.button("Clear job monitor"):
                st.session_state.pop('slurm_job_ids', None)
                st.session_state.pop('slurm_num_chunks', None)
                st.rerun()

        st.divider()

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

        # Escape paths for safe embedding in shell and Python string literals
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
#SBATCH --output="{af3_folder}/AF3_app_chunk{i}_%j.log"
#SBATCH --error="{af3_folder}/AF3_app_chunk{i}_%j.err"

echo "AF3 Analysis App - Chunk {i+1}/{len(chunks)}"
echo "Job ID: $SLURM_JOB_ID | Node: $SLURMD_NODENAME | CPUs: $SLURM_CPUS_PER_TASK"
echo "Predictions in chunk: {len(chunk)}"
echo "Start: $(date)"

source '{esc_venv}'
cd '{esc_app}'

python -c '
import json, sys
from pathlib import Path
from multiprocessing import Pool
from functools import partial
from core.analyzer import analyze_prediction_all_models

chunk_file = "{esc_chunk}"
af3_folder = "{esc_af3}"
output_file = "{esc_output}"

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
            # Save per-prediction JSON in the prediction dir
            pn = result_list[0].get("prediction_name", "")
            if pn:
                # Find matching resolved dir by name
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
'

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
        st.info(f"Submitted {len(job_ids)} jobs. Check progress with the command below.")
        st.code(f"squeue -u $USER | grep AF3app", language="bash")


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

    # Save merged results
    cache_file = af3_path / "af3_app_all_models_analysis.json"
    with open(cache_file, 'w') as f:
        json.dump(all_results, f, indent=2)

    st.text(f"  Saved {len(all_results)} models to {cache_file.name}")

    # Save per-prediction JSON files for fast lookup
    from collections import defaultdict
    from core.scanner import resolve_prediction_dir
    by_pred = defaultdict(list)
    for r in all_results:
        by_pred[r.get('prediction_name', '')].append(r)

    progress = st.progress(0, text="Writing per-prediction files...")
    n_pred_files = 0
    pred_names = [k for k in by_pred if k]
    for idx, pred_name in enumerate(pred_names):
        # Resolve AF3 Server nesting for correct save location
        target_dir = resolve_prediction_dir(af3_path / pred_name)
        pred_json = target_dir / "af3_app_analysis.json"
        try:
            with open(pred_json, 'w') as pf:
                json.dump(by_pred[pred_name], pf, indent=2)
            n_pred_files += 1
        except Exception:
            pass
        if idx % 50 == 0 or idx == len(pred_names) - 1:
            progress.progress((idx + 1) / len(pred_names), text=f"Writing per-prediction files... {idx + 1}/{len(pred_names)}")
    progress.empty()

    st.text(f"  Saved {n_pred_files} per-prediction JSON files")

    # Clean up chunk files
    for i in range(num_chunks):
        for pattern in [f"_slurm_chunk_{i}.txt", f"_slurm_results_{i}.json"]:
            f = af3_path / pattern
            f.unlink(missing_ok=True)
