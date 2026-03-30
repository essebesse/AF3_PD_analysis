"""
Page 2: Results - Sortable Table with All Models

Displays full analysis results with sorting, filtering, and export capabilities.
"""

import streamlit as st
import os
import sys
from pathlib import Path
from typing import Optional

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.scanner import AF3Scanner
from core.utils import format_score, calculate_confidence_tier, tier_color


def show_results(project_path: str, af3_folder: str):
    """Display the Results page."""

    st.header("📊 Analysis Results")

    if st.button("← Back: Load Data"):
        st.session_state['_navigate_to'] = "1. Load Data"
        st.rerun()

    if not os.path.isdir(af3_folder):
        st.error(f"AF3 folder not found: {af3_folder}")
        return

    # Load cached analysis or scan
    cache_file = Path(af3_folder) / "af3_app_all_models_analysis.json"

    if cache_file.exists():
        st.success(f"✓ Found cached analysis: {cache_file}")
        import json
        with open(cache_file, 'r') as f:
            analysis_results = json.load(f)
    else:
        st.warning("⚠️ No cached analysis found. Run full analysis on the Analyze step for per-model metrics including ipSAE.")

        # Quick scan with summary data for a basic overview
        try:
            from core.scanner import load_prediction_data
            scanner = AF3Scanner(Path(af3_folder))
            predictions = scanner.scan()
            scanner.load_all()  # reads summary_confidences.json + ranking_scores.csv

            # Build results from summary data (one row per prediction, top-ranked only)
            analysis_results = []
            for pred in scanner.predictions:
                if 'iptm' not in pred:
                    continue
                analysis_results.append({
                    'prediction_name': pred['name'],
                    'bait': pred['bait'],
                    'prey': pred['prey'],
                    'seed': 0,
                    'sample': 0,
                    'is_top_ranked': True,
                    'iptm': pred.get('iptm', 0),
                    'ptm': pred.get('ptm', 0),
                    'ranking_score': pred.get('ranking_score', 0),
                    'ipsae': None,  # Not calculated yet
                    'interface_plddt': None,
                    'contacts_pae3': None,
                    'contacts_pae5': None,
                    'contacts_pae8': None,
                    'fraction_disordered': pred.get('fraction_disordered', 0),
                    'has_full_analysis': False
                })
        except Exception as e:
            st.error(f"Error loading results: {e}")
            return

    # Filters
    st.subheader("🔍 Filters")
    col1, col2, col3, col4 = st.columns(4)

    with col1:
        confidence_filter = st.selectbox(
            "Confidence:",
            ["All", "High", "Medium", "Low", "Very Low"]
        )

    with col2:
        iptm_min, iptm_max = 0.0, 1.0
        iptm_range = st.slider(
            "iPTM Range:",
            min_value=0.0, max_value=1.0,
            value=(0.0, 1.0)
        )

    with col3:
        model_filter = st.radio("Models:", ["All models", "Top-ranked only"],
                                horizontal=True, label_visibility="collapsed")

    with col4:
        sort_by = st.selectbox(
            "Sort by:",
            ["ipSAE", "iPTM", "ranking_score", "interface_plddt", "contacts_pae5"]
        )

    # Apply filters
    filtered_results = analysis_results

    if confidence_filter != "All":
        def get_tier(r):
            if r.get('ipsae') is not None:
                return calculate_confidence_tier(r['iptm'], r['ipsae'])
            return calculate_confidence_tier(r['iptm'])

        filtered_results = [r for r in filtered_results if get_tier(r) == confidence_filter]

    if iptm_range != (0.0, 1.0):
        filtered_results = [r for r in filtered_results if iptm_range[0] <= r['iptm'] <= iptm_range[1]]

    if model_filter == "Top-ranked only":
        filtered_results = [r for r in filtered_results if r.get('is_top_ranked', False)]

    # Resolve gene names FIRST so we can use them everywhere
    import pandas as pd
    gene_cache = st.session_state.get('gene_name_cache', {})

    unique_accs = set()
    for r in filtered_results:
        pred_name = r['prediction_name']
        if '_and_' in pred_name:
            bait, prey = pred_name.split('_and_', 1)
            unique_accs.add(bait.upper())
            unique_accs.add(prey.upper())

    missing_accs = [acc for acc in unique_accs if acc not in gene_cache]
    if missing_accs:
        from core.utils import fetch_gene_names_batch
        new_genes = fetch_gene_names_batch(missing_accs)
        st.session_state.setdefault('gene_name_cache', {}).update(new_genes)
        gene_cache.update(new_genes)

    def _pred_label(pred_name):
        """Turn 'q9bw83_and_q9ujt0' into 'IFT27 × TUBE1 (Q9BW83 × Q9UJT0)'."""
        if '_and_' in pred_name:
            bait, prey = pred_name.split('_and_', 1)
            bait_acc = bait.upper()
            prey_acc = prey.upper()
            bait_gene = gene_cache.get(bait_acc, '')
            prey_gene = gene_cache.get(prey_acc, '')
            if bait_gene and prey_gene:
                return f"{bait_gene} × {prey_gene}"
            elif bait_gene:
                return f"{bait_gene} × {prey_acc}"
            else:
                return f"{bait_acc} × {prey_acc}"
        return pred_name

    # Display results count
    st.caption(f"Showing {len(filtered_results)} models from {len(set(r['prediction_name'] for r in filtered_results))} predictions")

    # Implement actual sorting
    sort_column_map = {
        "ipSAE": "ipsae",
        "iPTM": "iptm",
        "ranking_score": "ranking_score",
        "interface_plddt": "interface_plddt",
        "contacts_pae5": "contacts_pae5"
    }
    sort_key = sort_column_map.get(sort_by, "iptm")
    filtered_results.sort(key=lambda r: r.get(sort_key) or 0, reverse=True)

    display_data = []
    # Keep mapping from display row index to prediction name
    row_to_pred = []
    for r in filtered_results:
        model_label = "Top" if r.get('is_top_ranked') else f"s{r['seed']}-m{r['sample']}"

        pred_name = r['prediction_name']
        if '_and_' in pred_name:
            bait, prey = pred_name.split('_and_', 1)
        else:
            bait, prey = pred_name, ''
        bait_acc = bait.upper()
        prey_acc = prey.upper()
        bait_gene = gene_cache.get(bait_acc, '')
        prey_gene = gene_cache.get(prey_acc, '')
        # Always show gene name (if available) followed by UniProt ID
        bait_label = f"{bait_gene} ({bait_acc})" if bait_gene else bait_acc
        prey_label = f"{prey_gene} ({prey_acc})" if prey_gene else prey_acc

        def _contact(key):
            v = r.get(key)
            return int(v) if v is not None else "-"

        row_to_pred.append(pred_name)
        display_data.append({
            'Bait': bait_label,
            'Prey': prey_label,
            'Model': model_label,
            'iPTM': format_score(r['iptm']),
            'ipSAE': format_score(r.get('ipsae')) if r.get('ipsae') is not None else "-",
            'Ranking Score': format_score(r['ranking_score']),
            'iPLDDT': format_score(r.get('interface_plddt')) if r.get('interface_plddt') is not None else "-",
            'PAE≤3': _contact('contacts_pae3'),
            'PAE≤5': _contact('contacts_pae5'),
            'PAE≤8': _contact('contacts_pae8'),
        })

    df = pd.DataFrame(display_data)

    # Display with row selection — click a row to go to Detailed Analysis
    st.caption("Click a row to select it for detailed analysis")
    event = st.dataframe(
        df,
        width='stretch',
        hide_index=True,
        on_select="rerun",
        selection_mode="single-row",
        column_config={
            "Model": st.column_config.TextColumn(width="small"),
            "PAE≤3": st.column_config.NumberColumn(width="small"),
            "PAE≤5": st.column_config.NumberColumn(width="small"),
            "PAE≤8": st.column_config.NumberColumn(width="small"),
        },
        key="results_table"
    )

    # Handle row selection
    selected_rows = event.selection.rows if event and event.selection else []
    if selected_rows:
        row_idx = selected_rows[0]
        if row_idx < len(row_to_pred):
            selected_pred_name = row_to_pred[row_idx]
            st.session_state['selected_prediction'] = selected_pred_name
            pred_label = _pred_label(selected_pred_name)
            col1, col2 = st.columns([3, 1])
            with col1:
                st.info(f"Selected: **{pred_label}**")
            with col2:
                if st.button("Go to Detailed Analysis", type="primary"):
                    st.session_state['_navigate_to'] = "4. Detailed Analysis"
                    st.rerun()

    # Export options
    st.divider()
    col1, col2 = st.columns(2)

    with col1:
        csv = df.to_csv(index=False)
        st.download_button(
            label="Download CSV",
            data=csv,
            file_name="af3_analysis_results.csv",
            mime="text/csv"
        )

    with col2:
        try:
            import io
            buffer = io.BytesIO()
            df.to_excel(buffer, index=False, engine='openpyxl')
            buffer.seek(0)
            st.download_button(
                label="Download XLSX",
                data=buffer,
                file_name="af3_analysis_results.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )
        except Exception as e:
            st.error(f"XLSX export failed: {e}")