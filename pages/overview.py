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

        # Summary stats from cache
        best = {}
        for r in cached:
            name = r['prediction_name']
            iptm = r.get('iptm', 0)
            if name not in best or iptm > best[name].get('iptm', 0):
                best[name] = r
        per_pred = list(best.values())
        iptms = [r.get('iptm', 0) for r in per_pred]

        col1, col2, col3, col4 = st.columns(4)
        col1.metric("Total Predictions", total)
        col2.metric("Analyzed", len(per_pred))
        col3.metric("Avg iPTM", f"{sum(iptms)/len(iptms):.3f}" if iptms else "-")
        col4.metric("High Confidence (iPTM>0.5)", sum(1 for i in iptms if i > 0.5))

        # iPTM histogram
        import matplotlib.pyplot as plt
        fig_hist, ax_hist = plt.subplots(figsize=(10, 4))
        ax_hist.hist(iptms, bins=50, color='steelblue', edgecolor='white')
        ax_hist.axvline(x=0.5, color='red', linestyle='--', label='iPTM=0.5')
        ax_hist.axvline(x=0.8, color='green', linestyle='--', label='iPTM=0.8')
        ax_hist.set_xlabel('iPTM')
        ax_hist.set_ylabel('Count')
        ax_hist.set_title(f'iPTM Distribution ({len(per_pred)} predictions)')
        ax_hist.legend()
        st.pyplot(fig_hist)
        plt.close(fig_hist)

        st.divider()

        tab_preview, tab_full = st.tabs(["📋 Top Predictions", "🔬 All Models"])

        with tab_preview:
            gene_cache = st.session_state.get('gene_name_cache', {})
            top_20 = sorted(per_pred, key=lambda x: x.get('iptm', 0), reverse=True)[:20]

            display_data = []
            for r in top_20:
                pred_name = r['prediction_name']
                if '_and_' in pred_name:
                    bait, prey = pred_name.split('_and_', 1)
                else:
                    bait, prey = pred_name, ''
                bait_acc = bait.upper()
                prey_acc = prey.upper()
                bait_gene = gene_cache.get(bait_acc, '')
                prey_gene = gene_cache.get(prey_acc, '')
                bait_label = f"{bait_gene} ({bait_acc})" if bait_gene else bait_acc
                prey_label = f"{prey_gene} ({prey_acc})" if prey_gene else prey_acc
                display_data.append({
                    'Bait': bait_label,
                    'Prey': prey_label,
                    'iPTM': format_score(r.get('iptm', 0)),
                    'ipSAE': format_score(r.get('ipsae')) if r.get('ipsae') is not None else "-",
                    'Ranking Score': format_score(r.get('ranking_score', 0)),
                })
            df_preview = pd.DataFrame(display_data)
            st.dataframe(df_preview, width='stretch', hide_index=True)
            st.caption(f"Top 20 of {len(per_pred)} predictions by iPTM.")

        with tab_full:
            def _contact(r, key):
                v = r.get(key)
                return int(v) if v is not None else "-"

            gene_cache = st.session_state.get('gene_name_cache', {})

            unique_accs = set()
            for r in cached:
                pred_name = r['prediction_name']
                if '_and_' in pred_name:
                    bait, prey = pred_name.split('_and_', 1)
                    unique_accs.add(bait.upper())
                    unique_accs.add(prey.upper())

            missing_accs = [acc for acc in unique_accs if acc not in gene_cache]
            if missing_accs:
                st.info(f"Fetching gene names for {len(missing_accs)} proteins...")
                from core.utils import fetch_gene_names_batch
                new_genes = fetch_gene_names_batch(missing_accs)
                st.session_state.setdefault('gene_name_cache', {}).update(new_genes)
                gene_cache.update(new_genes)

            rows = []
            for r in cached:
                pred_name = r['prediction_name']
                if '_and_' in pred_name:
                    bait, prey = pred_name.split('_and_', 1)
                else:
                    bait, prey = pred_name, ''
                bait_acc = bait.upper()
                prey_acc = prey.upper()
                bait_gene = gene_cache.get(bait_acc, '')
                prey_gene = gene_cache.get(prey_acc, '')
                bait_label = f"{bait_gene} ({bait_acc})" if bait_gene else bait_acc
                prey_label = f"{prey_gene} ({prey_acc})" if prey_gene else prey_acc
                model_label = "Top" if r.get('is_top_ranked') else f"s{r['seed']}-m{r['sample']}"
                rows.append({
                    'Bait': bait_label,
                    'Prey': prey_label,
                    'Model': model_label,
                    'iPTM': format_score(r['iptm']),
                    'ipSAE': format_score(r.get('ipsae')) if r.get('ipsae') is not None else "-",
                    'Ranking Score': format_score(r['ranking_score']),
                    'iPLDDT': format_score(r.get('interface_plddt')) if r.get('interface_plddt') is not None else "-",
                    'PAE≤3': _contact(r, 'contacts_pae3'),
                    'PAE≤5': _contact(r, 'contacts_pae5'),
                    'PAE≤8': _contact(r, 'contacts_pae8'),
                })

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
        # No cache — just show the folder list
        st.metric("Total Predictions", total)
        st.info("No analysis results yet. Go to **Analyze** step to run full analysis.")

        # Show prediction names as a simple table
        display_data = []
        for pred in predictions[:100]:
            bait_label = pred['bait'].upper()
            prey_label = pred['prey'].upper() if pred['prey'] else ''
            display_data.append({'Bait': bait_label, 'Prey': prey_label})

        df = pd.DataFrame(display_data)
        st.dataframe(df, width='stretch', hide_index=True)
        if total > 100:
            st.caption(f"Showing first 100 of {total} predictions.")

    # Next step button
    st.divider()
    if st.button("Next: Run Analysis →", type="primary"):
        st.session_state['_navigate_to'] = "2. Analyze"
        st.rerun()
