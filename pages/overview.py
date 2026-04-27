"""
Page 1: Overview - Project Selection and Data Loading

Shows prediction list instantly (directory listing only).
Full data loaded from analysis cache if available.
"""

import streamlit as st
import os
from pathlib import Path

import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.scanner import AF3Scanner
from core.utils import format_score


@st.cache_data(ttl=300)
def cached_fast_scan(af3_folder: str):
    """Fast scan: just lists directories, no file I/O."""
    scanner = AF3Scanner(Path(af3_folder))
    return scanner.scan()


def show_load_data(project_path: str, af3_folder: str):
    """Display the Load Data step."""

    st.header("📊 Load Data")

    if not os.path.isdir(af3_folder):
        st.error(f"Predictions folder not found: {af3_folder}")
        return

    # Instant: just list directories
    predictions = cached_fast_scan(str(af3_folder))
    total = len(predictions)

    if total == 0:
        st.warning("No prediction directories found.")
        return

    import pandas as pd

    # Check if full analysis cache exists
    cache_file = Path(af3_folder) / "af3_app_all_models_analysis.json"
    has_cache = cache_file.exists()

    if has_cache:
        # Show rich view from cache
        import json
        with open(cache_file, 'r') as f:
            cached = json.load(f)

        # Which formats are in the cache? Drives column labels and which
        # metric we show (iPTM for AF3 vs iPTM+pTM for AF2).
        cache_formats = {r.get('format', 'af3_local') for r in cached}
        is_af2_cache = cache_formats == {'af2'}
        is_mixed = len(cache_formats) > 1

        # Unified "primary score" per entry:
        #   AF3 -> iptm, AF2 -> iptm_ptm. Both live on the unified 'iptm' key
        #   for AF3; for AF2 we prefer iptm_ptm.
        def _primary_score(r):
            if r.get('format') == 'af2':
                v = r.get('iptm_ptm')
            else:
                v = r.get('iptm')
            return v if v is not None else 0

        # Collapse to best model per prediction (by ipSAE if present, else primary score)
        best = {}
        for r in cached:
            name = r['prediction_name']
            key = r.get('ipsae') if r.get('ipsae') is not None else _primary_score(r)
            prev = best.get(name)
            prev_key = (prev.get('ipsae') if prev and prev.get('ipsae') is not None
                        else (_primary_score(prev) if prev else None))
            if prev is None or (key is not None and (prev_key is None or key > prev_key)):
                best[name] = r
        per_pred = list(best.values())
        primary_scores = [_primary_score(r) for r in per_pred]
        ipsae_scores = [r.get('ipsae') for r in per_pred if r.get('ipsae') is not None]

        primary_label = "iPTM+pTM" if is_af2_cache else "iPTM"
        primary_thresh = 0.5 if is_af2_cache else 0.5

        col1, col2, col3, col4 = st.columns(4)
        col1.metric("Total Predictions", total)
        col2.metric("Analyzed", len(per_pred))
        col3.metric(
            f"Avg {primary_label}",
            f"{sum(primary_scores)/len(primary_scores):.3f}" if primary_scores else "-",
        )
        # Use ipSAE-based high-confidence count when available (unified metric)
        if ipsae_scores:
            col4.metric("High Confidence (ipSAE>0.5)", sum(1 for i in ipsae_scores if i > 0.5))
        else:
            col4.metric(
                f"High Confidence ({primary_label}>{primary_thresh})",
                sum(1 for i in primary_scores if i > primary_thresh),
            )

        # Score histogram (primary score for the cache)
        import matplotlib.pyplot as plt
        fig_hist, ax_hist = plt.subplots(figsize=(10, 4))
        ax_hist.hist(primary_scores, bins=50, color='steelblue', edgecolor='white')
        ax_hist.axvline(x=0.5, color='red', linestyle='--', label=f'{primary_label}=0.5')
        ax_hist.axvline(x=0.8, color='green', linestyle='--', label=f'{primary_label}=0.8')
        ax_hist.set_xlabel(primary_label)
        ax_hist.set_ylabel('Count')
        ax_hist.set_title(f'{primary_label} Distribution ({len(per_pred)} predictions)')
        ax_hist.legend()
        st.pyplot(fig_hist)
        plt.close(fig_hist)

        st.divider()

        tab_preview, tab_full = st.tabs(["📋 Top Predictions", "🔬 All Models"])

        with tab_preview:
            if 'gene_name_cache' not in st.session_state:
                st.session_state['gene_name_cache'] = {}
            gene_cache = st.session_state['gene_name_cache']
            top_20 = sorted(per_pred, key=_primary_score, reverse=True)[:20]

            display_data = []
            for r in top_20:
                pred_name = r['prediction_name']
                fmt = r.get('format', 'af3_local')
                if fmt != 'af2' and '_and_' in pred_name:
                    bait, prey = pred_name.split('_and_', 1)
                    bait_acc, prey_acc = bait.upper(), prey.upper()
                    bait_gene = gene_cache.get(bait_acc, '')
                    prey_gene = gene_cache.get(prey_acc, '')
                    bait_label = f"{bait_gene} ({bait_acc})" if bait_gene else bait_acc
                    prey_label = f"{prey_gene} ({prey_acc})" if prey_gene else prey_acc
                else:
                    bait_label, prey_label = pred_name, ''
                # ⚠️ marker for multi-chain predictions
                if r.get('n_chains') and r['n_chains'] > 2:
                    bait_label = f"⚠️ {bait_label}"
                row = {
                    'Bait': bait_label,
                    'Prey': prey_label,
                    primary_label: format_score(_primary_score(r)),
                    'ipSAE': format_score(r.get('ipsae')) if r.get('ipsae') is not None else "-",
                }
                if is_mixed:
                    row['Format'] = {'af2': 'AF2', 'af3_local': 'AF3', 'af3_server_flat': 'AF3 Server'}.get(fmt, fmt)
                if not is_af2_cache:
                    row['Ranking Score'] = format_score(r.get('ranking_score')) if r.get('ranking_score') is not None else "-"
                display_data.append(row)
            df_preview = pd.DataFrame(display_data)
            st.dataframe(df_preview, width='stretch', hide_index=True)
            st.caption(f"Top 20 of {len(per_pred)} predictions by {primary_label}.")

        with tab_full:
            def _contact(r, key):
                v = r.get(key)
                return int(v) if v is not None else "-"

            if 'gene_name_cache' not in st.session_state:
                st.session_state['gene_name_cache'] = {}
            gene_cache = st.session_state['gene_name_cache']

            # Only AF3 names follow the <acc>_and_<acc> convention
            unique_accs = set()
            for r in cached:
                if r.get('format') == 'af2':
                    continue
                pred_name = r['prediction_name']
                if '_and_' in pred_name:
                    bait, prey = pred_name.split('_and_', 1)
                    unique_accs.add(bait.upper())
                    unique_accs.add(prey.upper())

            missing_accs = [acc for acc in unique_accs if acc not in gene_cache]
            if missing_accs:
                st.info(f"Fetching gene names for {len(missing_accs)} proteins...")
                from core.utils import fetch_gene_names_batch
                try:
                    new_genes = fetch_gene_names_batch(missing_accs)
                    gene_cache.update(new_genes)
                except Exception:
                    st.warning("Could not fetch gene names from UniProt.")

            rows = []
            for r in cached:
                pred_name = r['prediction_name']
                fmt = r.get('format', 'af3_local')
                if fmt != 'af2' and '_and_' in pred_name:
                    bait, prey = pred_name.split('_and_', 1)
                    bait_acc, prey_acc = bait.upper(), prey.upper()
                    bait_gene = gene_cache.get(bait_acc, '')
                    prey_gene = gene_cache.get(prey_acc, '')
                    bait_label = f"{bait_gene} ({bait_acc})" if bait_gene else bait_acc
                    prey_label = f"{prey_gene} ({prey_acc})" if prey_gene else prey_acc
                else:
                    bait_label, prey_label = pred_name, ''
                if r.get('n_chains') and r['n_chains'] > 2:
                    bait_label = f"⚠️ {bait_label}"
                model_label = "Top" if r.get('is_top_ranked') else f"s{r['seed']}-m{r['sample']}"
                row = {
                    'Bait': bait_label,
                    'Prey': prey_label,
                    'Model': model_label,
                    primary_label: format_score(_primary_score(r)),
                    'ipSAE': format_score(r.get('ipsae')) if r.get('ipsae') is not None else "-",
                    'iPLDDT': format_score(r.get('interface_plddt')) if r.get('interface_plddt') is not None else "-",
                    'PAE≤3': _contact(r, 'contacts_pae3'),
                    'PAE≤5': _contact(r, 'contacts_pae5'),
                    'PAE≤8': _contact(r, 'contacts_pae8'),
                }
                if is_mixed:
                    row['Format'] = {'af2': 'AF2', 'af3_local': 'AF3', 'af3_server_flat': 'AF3 Server'}.get(fmt, fmt)
                if not is_af2_cache:
                    row['Ranking Score'] = format_score(r.get('ranking_score')) if r.get('ranking_score') is not None else "-"
                rows.append(row)

            df_analysis = pd.DataFrame(rows)
            df_analysis = df_analysis.sort_values(
                by='ipSAE', ascending=False,
                key=lambda col: pd.to_numeric(col, errors='coerce').fillna(-1)
            )
            st.dataframe(
                df_analysis,
                width='stretch',
                hide_index=True,
                column_config={
                    "Model": st.column_config.TextColumn(width="small"),
                    "PAE≤3": st.column_config.NumberColumn(width="small"),
                    "PAE≤5": st.column_config.NumberColumn(width="small"),
                    "PAE≤8": st.column_config.NumberColumn(width="small"),
                }
            )
            st.caption(f"{len(cached)} models from cache · sorted by ipSAE")

    else:
        # No cache — just show the folder list (with format + chain warnings)
        st.metric("Total Predictions", total)
        st.info("No analysis results yet. Go to **Analyze** step to run full analysis.")

        fmts_here = {p.get('format', 'unknown') for p in predictions}
        show_format_col = len(fmts_here) > 1

        display_data = []
        for pred in predictions[:100]:
            bait_label = pred['bait'].upper() if pred.get('bait') else pred['name']
            prey_label = pred['prey'].upper() if pred.get('prey') else ''
            if pred.get('n_chains') and pred['n_chains'] > 2:
                bait_label = f"⚠️ {bait_label}"
            row = {'Bait': bait_label, 'Prey': prey_label}
            if show_format_col:
                row['Format'] = {'af2': 'AF2', 'af3_local': 'AF3', 'af3_server_flat': 'AF3 Server'}.get(
                    pred.get('format', 'unknown'), pred.get('format', 'unknown')
                )
            display_data.append(row)

        df = pd.DataFrame(display_data)
        st.dataframe(df, width='stretch', hide_index=True)
        if total > 100:
            st.caption(f"Showing first 100 of {total} predictions.")

    # Next step button
    st.divider()
    if st.button("Next: Run Analysis →", type="primary"):
        st.session_state['_navigate_to'] = "2. Analyze"
        st.rerun()
