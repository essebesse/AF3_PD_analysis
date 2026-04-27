"""
Page 3: Detailed Analysis - Per-Prediction Deep Dive

PAE plots, interface contacts, hub residues, and model comparison.
"""

import streamlit as st
import os
import sys
from pathlib import Path
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.scanner import AF3Scanner, load_prediction_data
from core.utils import format_score


def show_detailed(project_path: str, af3_folder: str):
    """Display the Detailed Analysis page."""

    st.header("🔬 Detailed Analysis")

    if st.button("← Back to Results"):
        st.session_state['_navigate_to'] = "3. Results"
        st.rerun()

    if not os.path.isdir(af3_folder):
        st.error(f"AF3 folder not found: {af3_folder}")
        return

    # Get list of predictions for selection (cached to avoid re-listing 44k dirs on every rerun)
    @st.cache_data(ttl=300)
    def _cached_scan(af3_dir: str):
        s = AF3Scanner(Path(af3_dir))
        return s.scan()

    try:
        predictions = _cached_scan(af3_folder)
    except Exception as e:
        st.error(f"Error loading predictions: {e}")
        return

    if not predictions:
        st.warning("No predictions found.")
        return

    # Prediction selector with search
    st.subheader("📁 Select Prediction")

    # Build gene name lookup for search
    gene_cache = st.session_state.get('gene_name_cache', {})
    pred_names = [p['name'] for p in predictions]

    # Search box: filter by gene name or UniProt ID
    search_query = st.text_input(
        "Search by gene name or UniProt ID:",
        placeholder="e.g. TUBG1, IFT27, q9bw83",
        key="detailed_search"
    )

    if search_query:
        q = search_query.strip().lower()
        filtered_names = []
        for name in pred_names:
            # Match against folder name
            if q in name.lower():
                filtered_names.append(name)
                continue
            # Match against gene names from cache
            if '_and_' in name:
                bait, prey = name.split('_and_', 1)
                bait_gene = gene_cache.get(bait.upper(), '').lower()
                prey_gene = gene_cache.get(prey.upper(), '').lower()
                if q in bait_gene or q in prey_gene:
                    filtered_names.append(name)
        if filtered_names:
            pred_names = filtered_names
        else:
            st.warning(f"No predictions matching '{search_query}'")

    col1, col2 = st.columns([3, 1])

    with col1:
        # Pre-select from session state if navigating from Results page
        default_idx = 0
        if 'selected_prediction' in st.session_state:
            try:
                default_idx = pred_names.index(st.session_state['selected_prediction'])
            except ValueError:
                default_idx = 0

        # Build display labels with gene names
        def _label(name):
            if '_and_' in name:
                bait, prey = name.split('_and_', 1)
                bg = gene_cache.get(bait.upper(), '')
                pg = gene_cache.get(prey.upper(), '')
                if bg and pg:
                    return f"{bg} × {pg}  ({name})"
                elif bg:
                    return f"{bg} × {prey.upper()}  ({name})"
            return name

        display_labels = [_label(n) for n in pred_names]
        label_to_name = dict(zip(display_labels, pred_names))

        selected_label = st.selectbox("Prediction:", display_labels, index=default_idx)
        selected_pred_name = label_to_name.get(selected_label, pred_names[0])

    with col2:
        # Lazily load full data for the selected prediction (cached)
        @st.cache_data(ttl=300)
        def _cached_load_pred(pred_dir_str):
            return load_prediction_data(Path(pred_dir_str))

        selected_pred = None
        for p in predictions:
            if p['name'] == selected_pred_name:
                if 'iptm' not in p:
                    full = _cached_load_pred(p['prediction_dir'])
                    if full:
                        p.update(full)
                selected_pred = p
                break

        if selected_pred:
            fmt = selected_pred.get('format', 'af3_local')
            if fmt == 'af2':
                iptm_ptm_val = selected_pred.get('iptm_ptm')
                if iptm_ptm_val is None:
                    iptm_ptm_val = selected_pred.get('iptm')
                if iptm_ptm_val is not None:
                    st.metric("iPTM+pTM", format_score(iptm_ptm_val))
            else:
                iptm_val = selected_pred.get('iptm')
                ptm_val = selected_pred.get('ptm')
                if iptm_val is not None:
                    st.metric("iPTM", format_score(iptm_val))
                if ptm_val is not None:
                    st.metric("PTM", format_score(ptm_val))

    if not selected_pred:
        return

    is_af2 = selected_pred.get('format') == 'af2'

    n = selected_pred.get('n_chains')
    if n and n > 2:
        st.info(
            f"This prediction has **{n} chains**. Pairwise scoring is limited "
            "to chains A and B; the third chain is shown in the structure but "
            "not in interface metrics."
        )

    st.divider()

    # Model selector within prediction
    st.subheader("🎲 Select Model")
    seed_samples = selected_pred.get('seed_samples', [])

    if seed_samples:
        # Find top-ranked model. AF2 ranks by iptm+ptm (stored as ranking_score
        # in the seed_samples list); AF3 by ranking_score directly.
        best_ss = max(seed_samples, key=lambda ss: ss.get('ranking_score') or 0)
        top_seed, top_sample = best_ss['seed'], best_ss['sample']

        # Build model options with top-ranked first (default)
        def _model_label(ss, is_top):
            score = ss.get('ranking_score') or 0
            if is_af2:
                core = f"ranked_{ss['sample']} (iPTM+pTM={score:.4f})"
            else:
                core = f"seed-{ss['seed']}_sample-{ss['sample']} (score={score:.4f})"
            return f"Top-ranked — {core}" if is_top else core

        top_option = None
        other_options = []
        for ss in seed_samples:
            is_top = ss['seed'] == top_seed and ss['sample'] == top_sample
            opt = (ss['seed'], ss['sample'], _model_label(ss, is_top))
            if is_top:
                top_option = opt
            else:
                other_options.append(opt)

        model_options = []
        if top_option:
            model_options.append(top_option)
        model_options.extend(other_options)

        selected_model = st.selectbox(
            "Model:",
            model_options,
            format_func=lambda x: x[2]
        )

        if selected_model:
            seed, sample = selected_model[0], selected_model[1]
            model_label = f"ranked_{sample}" if is_af2 else f"seed-{seed}_sample-{sample}"

            # Load detailed data for selected model
            st.divider()
            st.subheader(f"📊 Model: {model_label}")

            import json
            import numpy as np
            import pandas as pd
            import streamlit.components.v1 as components
            from core.analyzer import load_confidences
            from core.pae_plotter import (
                plot_pae_matrix,
                plot_pae_interface_zoom,
                get_chain_info_from_confidences,
                get_chain_indices_from_confidences
            )
            from core.viewer_3d import generate_viewer_html, compute_interface_pae_per_residue

            pred_dir = Path(selected_pred['prediction_dir'])
            af2_work_dir = Path(selected_pred['work_dir']) if is_af2 and selected_pred.get('work_dir') else None

            # ── Compute everything ONCE per model, store in session state ──
            model_key = f"{selected_pred_name}_{seed}_{sample}"
            if st.session_state.get('_detail_model_key') != model_key:
                # Model changed — clear old cached plots/data
                old_key = st.session_state.get('_detail_model_key', '')
                old_pred = old_key.rsplit('_', 2)[0] if old_key else ''
                for k in list(st.session_state.keys()):
                    if old_key and k.startswith(('_pae_plot_', '_zoom_plot_', '_iface_')) and old_key in k:
                        del st.session_state[k]
                    elif old_pred and k.startswith('_comp_plot_') and old_pred in k:
                        del st.session_state[k]

                confidences = None
                chain_ids_resolved: list = []
                chain_lengths_resolved: list = []
                protein_names = ['Chain A', 'Chain B']
                structure_file: Path = Path('')
                structure_format = 'cif'
                structure_content: str = ''
                pae_data = None
                viewer_html = None

                if is_af2 and af2_work_dir is not None:
                    # AF2 path: synthesize confidences + load PDB
                    from core.af2_detail import (
                        load_af2_model_data,
                        compute_af2_interface_pae_per_residue,
                    )
                    bundle = load_af2_model_data(af2_work_dir, sample)
                    if bundle is not None:
                        confidences = bundle['confidences']
                        chain_ids_resolved = bundle['chain_ids']
                        chain_lengths_resolved = bundle['chain_lengths']
                        structure_file = bundle['pdb_path']
                        structure_format = 'pdb'
                        structure_content = bundle['pdb_content']

                        # AF2 folder names don't follow <bait>_and_<prey>; show
                        # gene names from the cache when the folder name happens
                        # to match a known accession, else fall back to chain IDs
                        protein_names = []
                        for i, cid in enumerate(chain_ids_resolved):
                            protein_names.append(cid if cid else f"Chain {i+1}")

                        if chain_ids_resolved and chain_lengths_resolved:
                            pae_data = compute_af2_interface_pae_per_residue(
                                structure_file, confidences,
                                chain_ids_resolved, chain_lengths_resolved,
                            )
                else:
                    # AF3 path (unchanged)
                    confidences = load_confidences(pred_dir, selected_pred_name, seed, sample)
                    if confidences is not None:
                        chain_ids_resolved, chain_lengths_resolved = get_chain_info_from_confidences(confidences)
                        if '_and_' in selected_pred_name:
                            bait_acc, prey_acc = selected_pred_name.split('_and_', 1)
                            accs = [bait_acc.upper(), prey_acc.upper()]
                        else:
                            accs = [selected_pred_name.upper()]
                        protein_names = []
                        for i, cid in enumerate(chain_ids_resolved):
                            if i < len(accs):
                                acc = accs[i]
                                protein_names.append(gene_cache.get(acc, acc))
                            else:
                                protein_names.append(cid)

                    sample_dir = pred_dir / f"seed-{seed}_sample-{sample}"
                    cif_file = sample_dir / f"{selected_pred_name}_seed-{seed}_sample-{sample}_model.cif"
                    if not cif_file.exists():
                        cif_file = sample_dir / "model.cif"
                    if not cif_file.exists():
                        cif_file = pred_dir / f"{selected_pred_name}_model.cif"
                    if not cif_file.exists():
                        cif_file = pred_dir / "model.cif"
                    structure_file = cif_file
                    structure_format = 'cif'
                    if cif_file.exists():
                        structure_content = cif_file.read_text()

                    if confidences is not None and chain_ids_resolved and chain_lengths_resolved and cif_file.exists():
                        pae_data = compute_interface_pae_per_residue(
                            cif_file, confidences, chain_ids_resolved, chain_lengths_resolved
                        )

                if structure_content:
                    viewer_html = generate_viewer_html(
                        structure_content,
                        pae_residue_data=pae_data,
                        chain_names=protein_names,
                        chain_ids_all=chain_ids_resolved,
                        height=600,
                        structure_format=structure_format,
                    )

                # Store everything in session state
                st.session_state['_detail_model_key'] = model_key
                st.session_state['_detail_confidences'] = confidences
                st.session_state['_detail_chain_ids'] = chain_ids_resolved
                st.session_state['_detail_chain_lengths'] = chain_lengths_resolved
                st.session_state['_detail_protein_names'] = protein_names
                st.session_state['_detail_cif_file'] = str(structure_file)
                st.session_state['_detail_structure_format'] = structure_format
                st.session_state['_detail_pae_data'] = pae_data
                st.session_state['_detail_viewer_html'] = viewer_html

            # Retrieve from session state (instant — no recomputation)
            confidences = st.session_state.get('_detail_confidences')
            chain_ids_resolved = st.session_state.get('_detail_chain_ids', [])
            chain_lengths_resolved = st.session_state.get('_detail_chain_lengths', [])
            protein_names = st.session_state.get('_detail_protein_names', ['Chain A', 'Chain B'])
            cif_file = Path(st.session_state.get('_detail_cif_file', ''))
            structure_format = st.session_state.get('_detail_structure_format', 'cif')
            pae_data = st.session_state.get('_detail_pae_data')
            viewer_html = st.session_state.get('_detail_viewer_html')

            # Tabs for different analysis views
            tab0, tab1, tab2, tab3, tab5, tab4 = st.tabs(["3D Structure", "PAE Plots", "Interface Contacts", "Model Comparison", "PyMOL Scripts", "Export"])

            with tab0:
                st.markdown("### 3D Structure Viewer")

                if viewer_html:
                    components.html(viewer_html, height=620, scrolling=False)
                else:
                    st.warning(f"Structure file not found for {model_label}.")

            with tab1:
                st.markdown("### PAE Matrix Visualization")

                if confidences is not None:
                    pae_matrix = np.array(confidences.get('pae', []))

                    if pae_matrix.size > 0:
                        chain_indices = get_chain_indices_from_confidences(confidences)

                        # Use cached PNG from session state
                        pae_plot_key = f"_pae_plot_{model_key}"
                        if pae_plot_key not in st.session_state:
                            fig = plot_pae_matrix(
                                pae_matrix, protein_names, chain_lengths_resolved,
                                title=f'{selected_pred_name} (seed-{seed}_sample-{sample})'
                            )
                            import io
                            buf = io.BytesIO()
                            fig.savefig(buf, format='png', dpi=150, bbox_inches='tight')
                            plt.close(fig)
                            st.session_state[pae_plot_key] = buf.getvalue()
                        st.image(st.session_state[pae_plot_key])
                        pae_fname = f"{'_'.join(protein_names)}_s{seed}_m{sample}_PAE.png"
                        st.download_button("Download PAE plot", st.session_state[pae_plot_key],
                                           file_name=pae_fname, mime="image/png")

                        if len(chain_ids_resolved) >= 2 and len(chain_indices) >= 2:
                            st.subheader("Interface Zoom")
                            chain_a_name = protein_names[0]
                            chain_b_name = protein_names[1]
                            chain_a_idx = chain_indices[chain_ids_resolved[0]]
                            chain_b_idx = chain_indices[chain_ids_resolved[1]]

                            zoom_plot_key = f"_zoom_plot_{model_key}"
                            if zoom_plot_key not in st.session_state:
                                fig_zoom = plot_pae_interface_zoom(
                                    pae_matrix, chain_a_idx, chain_b_idx,
                                    chain_a_name, chain_b_name
                                )
                                buf = io.BytesIO()
                                fig_zoom.savefig(buf, format='png', dpi=150, bbox_inches='tight')
                                plt.close(fig_zoom)
                                st.session_state[zoom_plot_key] = buf.getvalue()
                            st.image(st.session_state[zoom_plot_key])
                            zoom_fname = f"{'_'.join(protein_names)}_s{seed}_m{sample}_PAE_interface.png"
                            st.download_button("Download interface zoom", st.session_state[zoom_plot_key],
                                               file_name=zoom_fname, mime="image/png")
                    else:
                        st.warning("PAE matrix not available in confidences file.")
                else:
                    st.warning(f"Could not load confidences for seed-{seed}_sample-{sample}. Check that the files exist.")

            with tab2:
                st.markdown("### Interface Contacts")

                from core.interface_analyzer import analyze_interface, map_interface_regions
                from core.af2_detail import analyze_af2_interface

                if cif_file.exists() and confidences is not None:
                    # Run interface analysis (cached in session state)
                    iface_key = f"_iface_{model_key}"
                    if iface_key not in st.session_state:
                        if is_af2:
                            st.session_state[iface_key] = analyze_af2_interface(
                                cif_file, confidences, pae_cutoff=10.0, distance_cutoff=8.0
                            )
                        else:
                            st.session_state[iface_key] = analyze_interface(
                                cif_file, confidences, pae_cutoff=10.0, distance_cutoff=8.0
                            )
                    interface_result = st.session_state[iface_key]

                    contacts = interface_result.get('contacts', [])
                    summary = interface_result.get('summary', {})
                    chemical_breakdown = interface_result.get('chemical_breakdown', {})
                    hub_residues_a = interface_result.get('hub_residues_a', [])
                    hub_residues_b = interface_result.get('hub_residues_b', [])

                    if contacts:
                        st.success(f"Found {len(contacts)} interface contacts (PAE≤10Å, Cβ≤8Å)")

                        # Show summary metrics
                        col1, col2, col3, col4 = st.columns(4)
                        col1.metric("Total Contacts", summary.get('total_contacts', 0))
                        col2.metric("Mean PAE", f"{summary.get('mean_pae', 0):.2f} Å")
                        col3.metric("Mean Distance", f"{summary.get('mean_distance', 0):.2f} Å")
                        col4.metric("Interface Area", f"{summary.get('interface_area_estimate', 0)} residues")

                        # Show chemical interaction breakdown
                        st.subheader("Chemical Interaction Types")
                        import pandas as pd
                        breakdown_df = pd.DataFrame([
                            {'Type': k, 'Count': v}
                            for k, v in sorted(chemical_breakdown.items(), key=lambda x: x[1], reverse=True)
                        ])
                        st.bar_chart(breakdown_df.set_index('Type'))

                        # Show top contacts table with protein names as column headers
                        contact_df = pd.DataFrame(contacts[:50])  # Show top 50
                        name_a = protein_names[0] if len(protein_names) > 0 else "Chain A"
                        name_b = protein_names[1] if len(protein_names) > 1 else "Chain B"
                        col_rename = {
                            'chain_a_res': f'{name_a} res',
                            'chain_b_res': f'{name_b} res',
                            'chain_a_aa': f'{name_a} aa',
                            'chain_b_aa': f'{name_b} aa',
                        }
                        contact_df = contact_df.rename(columns=col_rename)
                        st.dataframe(contact_df, width='stretch')

                        # Hub residues analysis
                        st.subheader("Hub Residues")
                        name_a = protein_names[0] if len(protein_names) > 0 else "Chain A"
                        name_b = protein_names[1] if len(protein_names) > 1 else "Chain B"

                        col1, col2 = st.columns(2)
                        with col1:
                            st.markdown(f"### {name_a} Hubs")
                            if hub_residues_a:
                                for res_id, res_name, count in hub_residues_a[:10]:
                                    st.write(f"Res {res_id} ({res_name}): {count} contacts")
                            else:
                                st.write("No hub residues found")

                        with col2:
                            st.markdown(f"### {name_b} Hubs")
                            if hub_residues_b:
                                for res_id, res_name, count in hub_residues_b[:10]:
                                    st.write(f"Res {res_id} ({res_name}): {count} contacts")
                            else:
                                st.write("No hub residues found")

                        # Interface region mapping
                        if len(chain_lengths_resolved) >= 2:
                            st.subheader("Interface Sequence Regions")
                            region_info = map_interface_regions(contacts, chain_lengths_resolved[0], chain_lengths_resolved[1])

                            col1, col2 = st.columns(2)
                            with col1:
                                st.markdown(f"### {name_a} Regions")
                                for start, end in region_info['chain_a_regions']:
                                    st.write(f"Residues {start}-{end}")
                                st.write(f"Coverage: {region_info['chain_a_coverage']*100:.1f}%")

                            with col2:
                                st.markdown(f"### {name_b} Regions")
                                for start, end in region_info['chain_b_regions']:
                                    st.write(f"Residues {start}-{end}")
                                st.write(f"Coverage: {region_info['chain_b_coverage']*100:.1f}%")
                    else:
                        st.warning("No interface contacts found.")
                else:
                    st.warning("Required files not found. Run full analysis first.")

            with tab3:
                st.markdown("### All Models Comparison")

                from core.analyzer import load_summary
                from core.pae_plotter import plot_all_models_comparison

                # Load analysis results: prefer per-prediction JSON (fast), fall back to big JSON
                @st.cache_data(ttl=300)
                def _load_pred_analysis(pred_dir_str, pred_name):
                    """Load per-prediction af3_app_analysis.json (small, fast)."""
                    pred_json = Path(pred_dir_str) / "af3_app_analysis.json"
                    if pred_json.exists():
                        with open(pred_json, 'r') as f:
                            entries = json.load(f)
                        return {(e['seed'], e['sample']): e for e in entries}
                    return None

                @st.cache_data(ttl=300)
                def _load_big_analysis_index(af3_dir):
                    """Fallback: load big JSON and index by prediction_name."""
                    cache_file = Path(af3_dir) / "af3_app_all_models_analysis.json"
                    index = {}
                    if cache_file.exists():
                        with open(cache_file, 'r') as f:
                            for entry in json.load(f):
                                pn = entry['prediction_name']
                                if pn not in index:
                                    index[pn] = {}
                                index[pn][(entry['seed'], entry['sample'])] = entry
                    return index

                # Try per-prediction file first (tiny, instant)
                cached_models = _load_pred_analysis(str(pred_dir), selected_pred_name)
                if cached_models is None:
                    # Fall back to big JSON (cached after first load)
                    all_analysis = _load_big_analysis_index(af3_folder)
                    cached_models = all_analysis.get(selected_pred_name, {})

                # Compare all 5 models in a table with per-model metrics
                comparison_data = []
                for ss in seed_samples:
                    is_top = ss['seed'] == top_seed and ss['sample'] == top_sample
                    cache_key = (ss['seed'], ss['sample'])
                    cached = cached_models.get(cache_key, {})

                    if is_af2:
                        # AF2 has no separate iPTM/PTM — only iPTM+pTM
                        iptm_ptm_val = ss.get('iptm_ptm') or ss.get('ranking_score') or 0
                        comparison_data.append({
                            'Model': f"ranked_{ss['sample']}",
                            'iPTM+pTM': format_score(iptm_ptm_val),
                            'ipSAE': format_score(cached.get('ipsae')) if cached.get('ipsae') is not None else "-",
                            'iPLDDT': format_score(cached.get('interface_plddt')) if cached.get('interface_plddt') is not None else "-",
                            'PAE≤5': str(cached.get('contacts_pae5', '-')),
                            'Top-ranked': "Yes" if is_top else "No",
                        })
                    else:
                        # Get per-model iPTM/PTM from summary or fall back to top-level
                        model_summary = load_summary(pred_dir, selected_pred_name, ss['seed'], ss['sample'])
                        model_iptm = model_summary.get('iptm', 0) if model_summary else selected_pred.get('iptm', 0)
                        model_ptm = model_summary.get('ptm', 0) if model_summary else selected_pred.get('ptm', 0)
                        comparison_data.append({
                            'Model': f"seed-{ss['seed']}_sample-{ss['sample']}",
                            'iPTM': format_score(model_iptm),
                            'PTM': format_score(model_ptm),
                            'ipSAE': format_score(cached.get('ipsae')) if cached.get('ipsae') is not None else "-",
                            'iPLDDT': format_score(cached.get('interface_plddt')) if cached.get('interface_plddt') is not None else "-",
                            'PAE≤5': str(cached.get('contacts_pae5', '-')),
                            'Ranking Score': format_score(ss['ranking_score']),
                            'Top-ranked': "Yes" if is_top else "No"
                        })

                comp_df = pd.DataFrame(comparison_data)
                st.dataframe(comp_df, width='stretch')

                # Load PAE matrices for all models and show side-by-side comparison
                st.subheader("PAE Matrix Comparison")

                comp_plot_key = f"_comp_plot_{selected_pred_name}"
                if comp_plot_key not in st.session_state:
                    pae_matrices = []
                    model_labels = []
                    if is_af2 and af2_work_dir is not None:
                        from core.af2_detail import load_af2_model_data as _load_af2
                        for ss in seed_samples:
                            is_top_ss = ss['seed'] == top_seed and ss['sample'] == top_sample
                            label = "Top" if is_top_ss else f"ranked_{ss['sample']}"
                            bundle = _load_af2(af2_work_dir, ss['sample'])
                            if bundle and bundle['confidences'].get('pae'):
                                pae_matrices.append(np.array(bundle['confidences']['pae']))
                                model_labels.append(label)
                    else:
                        for ss in seed_samples:
                            is_top_ss = ss['seed'] == top_seed and ss['sample'] == top_sample
                            label = "Top" if is_top_ss else f"s{ss['seed']}-m{ss['sample']}"
                            conf = load_confidences(pred_dir, selected_pred_name, ss['seed'], ss['sample'])
                            if conf and conf.get('pae'):
                                pae_matrices.append(np.array(conf['pae']))
                                model_labels.append(label)

                    if pae_matrices:
                        import io
                        fig = plot_all_models_comparison(pae_matrices, model_labels, chain_lengths_resolved, protein_names)
                        buf = io.BytesIO()
                        fig.savefig(buf, format='png', dpi=150, bbox_inches='tight')
                        plt.close(fig)
                        st.session_state[comp_plot_key] = buf.getvalue()
                    else:
                        st.session_state[comp_plot_key] = None

                if st.session_state[comp_plot_key]:
                    st.image(st.session_state[comp_plot_key])
                else:
                    st.warning("Could not load PAE matrices for comparison.")

            with tab5:
                st.markdown("### PyMOL Visualization Scripts")
                st.caption("Generate .py scripts to load in PyMOL for offline analysis")

                from core.pymol_script import generate_pymol_script

                if cif_file.exists() and confidences is not None and len(chain_ids_resolved) >= 2:
                    # Reuse PAE data and contacts from session state (already computed)
                    pae_data_pymol = pae_data

                    # Reuse interface analysis from tab2 (or compute if not yet done)
                    iface_key = f"_iface_{model_key}"
                    if iface_key not in st.session_state:
                        if is_af2:
                            from core.af2_detail import analyze_af2_interface as _analyze_iface_af2
                            st.session_state[iface_key] = _analyze_iface_af2(
                                cif_file, confidences, pae_cutoff=10.0, distance_cutoff=8.0
                            )
                        else:
                            from core.interface_analyzer import analyze_interface as _analyze_iface
                            st.session_state[iface_key] = _analyze_iface(
                                cif_file, confidences, pae_cutoff=10.0, distance_cutoff=8.0
                            )
                    iface_result = st.session_state[iface_key]
                    contacts_pymol = iface_result.get('contacts', []) if iface_result else None

                    name_a = protein_names[0] if len(protein_names) > 0 else chain_ids_resolved[0]
                    name_b = protein_names[1] if len(protein_names) > 1 else chain_ids_resolved[1]

                    pymol_script = generate_pymol_script(
                        cif_path=str(cif_file.resolve()),
                        prediction_name=selected_pred_name,
                        chain_a_id=chain_ids_resolved[0],
                        chain_b_id=chain_ids_resolved[1],
                        protein_a_name=name_a,
                        protein_b_name=name_b,
                        interface_residues=pae_data_pymol,
                        contacts=contacts_pymol,
                        seed=seed,
                        sample=sample,
                    )

                    script_filename = f"{selected_pred_name}_s{seed}_m{sample}_pymol.py"

                    st.download_button(
                        label="Download PyMOL Interface Script",
                        data=pymol_script,
                        file_name=script_filename,
                        mime="text/x-python",
                    )

                    with st.expander("Preview script"):
                        st.code(pymol_script, language="python")

                    st.markdown("""
**Usage in PyMOL:**
1. Open PyMOL
2. `File` → `Run Script` → select the downloaded `.py` file
3. Use commands: `show_interface`, `show_contacts`, `show_pae_only`

**Color scheme:**
- Forest green: PAE < 3 A (very high confidence)
- Green: PAE 3-5 A (high confidence)
- Yellow: PAE 5-8 A (moderate)
- Orange: PAE 8-12 A (low)
- Light blue / light orange: non-interface chain coloring
""")
                else:
                    st.warning("CIF file and confidences required for PyMOL script generation.")

            with tab4:
                st.markdown("### Export Options")

                import numpy as np
                import io

                if confidences is not None:
                    pae_matrix = np.array(confidences.get('pae', []))

                    if pae_matrix.size > 0:
                        csv_buffer = io.StringIO()
                        np.savetxt(csv_buffer, pae_matrix, delimiter=',', fmt='%.4f')

                        st.download_button(
                            label="Download PAE Matrix (CSV)",
                            data=csv_buffer.getvalue(),
                            file_name=f"{selected_pred_name}_s{seed}_m{sample}_pae_matrix.csv",
                            mime="text/csv"
                        )

                    # Export contacts (reuse from interface analysis if available)
                    if cif_file.exists():
                        iface_key = f"_iface_{model_key}"
                        if iface_key in st.session_state and st.session_state[iface_key]:
                            contacts = st.session_state[iface_key].get('contacts', [])
                        elif is_af2:
                            from core.af2_detail import analyze_af2_interface as _af2_iface
                            contacts = _af2_iface(cif_file, confidences).get('contacts', [])
                        else:
                            from core.analyzer import calculate_interface_contacts
                            contacts = calculate_interface_contacts(cif_file, confidences)

                        if contacts:
                            import pandas as pd
                            contact_df = pd.DataFrame(contacts)
                            csv = contact_df.to_csv(index=False)

                            st.download_button(
                                label="Download Interface Contacts (CSV)",
                                data=csv,
                                file_name=f"{selected_pred_name}_s{seed}_m{sample}_contacts.csv",
                                mime="text/csv"
                            )
                else:
                    st.warning("No confidences loaded for this model.")