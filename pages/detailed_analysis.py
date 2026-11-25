"""
Page 3: Detailed Analysis - Per-Prediction Deep Dive

PAE plots, interface contacts, hub residues, and model comparison.
"""

import io
import streamlit as st
import os
import sys
from pathlib import Path
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.scanner import AF3Scanner, load_prediction_data
from core.utils import TABLE_H, format_score, split_prediction_name


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
            bait_acc, prey_acc = split_prediction_name(name)
            bait_gene = gene_cache.get(bait_acc, '').lower()
            prey_gene = gene_cache.get(prey_acc, '').lower()
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
            bait_acc, prey_acc = split_prediction_name(name)
            if prey_acc:
                bg = gene_cache.get(bait_acc, '')
                pg = gene_cache.get(prey_acc, '')
                if bg and pg:
                    return f"{bg} × {pg}  ({name})"
                elif bg:
                    return f"{bg} × {prey_acc}  ({name})"
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

        if selected_pred and 'iptm' in selected_pred:
            st.metric("iPTM", format_score(selected_pred['iptm']))
            st.metric("PTM", format_score(selected_pred['ptm']))

    if not selected_pred:
        return

    # Format-specific notice — AlphaPulldown / AF2-multimer has no PAE matrix,
    # so PAE plots, ipSAE, and PAE-filtered contacts are not available.
    if selected_pred.get('format') == 'alphapulldown':
        st.warning(
            "**AF2-multimer prediction without PAE files.** "
            "No `ranked_N_PAE.csv` files were found (AlphaPulldown keeps the PAE "
            "matrix only inside its `result_*.pkl` files, which the app does not read). "
            "Available: 3D structure, spatial Cβ–Cβ interface contacts, "
            "interface pLDDT (from B-factors), hub residues, PyMOL/contacts export. "
            "Not available: ipSAE, PAE plots, PAE-filtered contacts, PAE matrix export."
        )

    st.divider()

    # Model selector within prediction
    st.subheader("🎲 Select Model")
    seed_samples = selected_pred.get('seed_samples', [])
    if selected_pred.get('format', 'af3') == 'af3' and selected_pred.get('prediction_dir'):
        # Drop models whose seed-N_sample-M files were deleted (see scanner)
        from core.scanner import af3_models_on_disk
        n_listed = len(seed_samples)
        seed_samples = af3_models_on_disk(Path(selected_pred['prediction_dir']), seed_samples)
        if len(seed_samples) < n_listed:
            st.caption(f"{n_listed - len(seed_samples)} of {n_listed} models are listed in "
                       "ranking_scores.csv but their seed/sample files are not on disk.")

    if seed_samples:
        # Find top-ranked model (highest ranking score)
        best_ss = max(seed_samples, key=lambda ss: ss['ranking_score'])
        top_seed, top_sample = best_ss['seed'], best_ss['sample']

        # Build model options with top-ranked first (default)
        top_option = None
        other_options = []
        for ss in seed_samples:
            is_top = ss['seed'] == top_seed and ss['sample'] == top_sample
            if is_top:
                label = f"Top-ranked (seed={ss['seed']}, sample={ss['sample']}, score={ss['ranking_score']:.4f})"
                top_option = (ss['seed'], ss['sample'], label)
            else:
                label = f"seed-{ss['seed']}_sample-{ss['sample']} (score={ss['ranking_score']:.4f})"
                other_options.append((ss['seed'], ss['sample'], label))

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

            # Load detailed data for selected model
            st.divider()
            st.subheader(f"📊 Model: seed-{seed}_sample-{sample}")

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
            # Prefix of the AF3 file names; differs from the folder name if renamed
            from core.scanner import prediction_file_prefix
            file_prefix = prediction_file_prefix(pred_dir)
            # AF2 ranked_N.pdb may sit in pred_dir/seq/
            _is_af2 = selected_pred.get('format') in ('alphapulldown', 'af2')
            struct_dir = pred_dir
            if _is_af2:
                from core.af2_pae import af2_work_dir
                struct_dir = af2_work_dir(pred_dir) or pred_dir

            # ── Compute everything ONCE per model, store in session state ──
            model_key = f"{selected_pred_name}_{seed}_{sample}"
            if st.session_state.get('_detail_model_key') != model_key:
                # Model changed. Exact-name keys (_pae_plot_, _comp_plot_) are
                # popped directly; the per-pair caches (interface zoom, interface
                # contacts, per-pair viewer HTML and interface PAE) carry the
                # model key as a prefix, so every key starting with it is dropped
                # (prefix, not substring, so predictions whose names share a
                # prefix do not cross-clear).
                old_key = st.session_state.get('_detail_model_key', '')
                old_pred = old_key.rsplit('_', 2)[0] if old_key else ''
                to_drop = []
                if old_key:
                    to_drop.append(f'_pae_plot_{old_key}')
                if old_pred:
                    to_drop.append(f'_comp_plot_{old_pred}')
                for k in to_drop:
                    st.session_state.pop(k, None)
                if old_key:
                    for _prefix in (f'_zoom_plot_{old_key}',
                                    f'_iface_{old_key}',
                                    f'_detail_pae_data_{old_key}',
                                    f'_detail_viewer_html_{old_key}'):
                        for k in list(st.session_state.keys()):
                            if k.startswith(_prefix):
                                st.session_state.pop(k, None)

                # Recompute all derived data
                confidences = load_confidences(pred_dir, file_prefix, seed, sample)
                pae_problem = None
                if confidences is None and _is_af2:
                    from core.af2_pae import load_af2_confidences_checked
                    pae_problem = load_af2_confidences_checked(struct_dir, sample)[1]

                chain_ids_resolved = []
                chain_lengths_resolved = []
                protein_names = ['Chain A', 'Chain B']
                pred_format = selected_pred.get('format', 'af3')

                if confidences is not None:
                    chain_ids_resolved, chain_lengths_resolved = get_chain_info_from_confidences(confidences)
                elif _is_af2:
                    # No PAE/confidences — derive chain info from the PDB itself.
                    try:
                        import gemmi as _gemmi
                        _structure = _gemmi.read_structure(str(
                            struct_dir / f"ranked_{sample}.pdb"
                        ))
                        if len(_structure) > 0:
                            _model = _structure[0]
                            for _ch in _model:
                                _residues = [r for r in _ch]
                                if _residues:
                                    chain_ids_resolved.append(_ch.name)
                                    chain_lengths_resolved.append(len(_residues))
                    except Exception:
                        pass

                # Find structure file — local-pipeline CIF, AF3 Server flat,
                # or AlphaPulldown (AF2-multimer) PDB
                sample_dir = pred_dir / f"seed-{seed}_sample-{sample}"
                cif_candidates = [
                    sample_dir / f"{file_prefix}_seed-{seed}_sample-{sample}_model.cif",
                    sample_dir / "model.cif",
                    pred_dir / f"{file_prefix}_model_{sample}.cif",  # AF3 Server flat
                    struct_dir / f"ranked_{sample}.pdb",                       # AlphaPulldown
                    pred_dir / f"{file_prefix}_model.cif",
                    pred_dir / "model.cif",
                ]
                cif_file = next((p for p in cif_candidates if p.exists()),
                                cif_candidates[0])

                # Chain labels. When there are exactly two polymer chains they
                # get the bait / prey names from the folder name if it has the
                # form <bait>_and_<prey>; with more than two polymer chains the
                # folder name no longer describes a single pair, so every
                # polymer chain is labelled "Chain X". Ligand and ion chains
                # get their residue name, e.g. "ATP (C)"; anything else "Chain X".
                nonpolymer_ids = []
                poly_ids = []
                if chain_ids_resolved:
                    from core.analyzer import polymer_chain_ids
                    if confidences is not None:
                        poly_ids = list(polymer_chain_ids(confidences))
                        _polymers = set(poly_ids)
                        nonpolymer_ids = [c for c in chain_ids_resolved if c not in _polymers]
                    else:
                        # AF2 without PAE/confidences: every parsed chain counts
                        # (polymer vs. ligand cannot be told apart).
                        poly_ids = list(chain_ids_resolved)
                    resnames = {}
                    if nonpolymer_ids and cif_file.exists():
                        try:
                            import gemmi as _gemmi
                            for _ch in _gemmi.read_structure(str(cif_file))[0]:
                                if _ch.name in nonpolymer_ids and len(_ch) > 0:
                                    resnames[_ch.name] = _ch[0].name
                        except Exception:
                            pass
                    pair_labels = {}
                    bait_acc, prey_acc = split_prediction_name(selected_pred_name)
                    if prey_acc and len(poly_ids) == 2:
                        # gene_cache stores '' for unresolved accessions, so
                        # fall back to the name from the folder.
                        pair_labels = {poly_ids[0]: gene_cache.get(bait_acc) or bait_acc,
                                       poly_ids[1]: gene_cache.get(prey_acc) or prey_acc}
                    protein_names = []
                    for cid in chain_ids_resolved:
                        if cid in pair_labels:
                            protein_names.append(pair_labels[cid])
                        elif cid in resnames:
                            protein_names.append(f"{resnames[cid]} ({cid})")
                        else:
                            protein_names.append(f"Chain {cid}")

                # Store everything in session state
                st.session_state['_detail_model_key'] = model_key
                st.session_state['_detail_confidences'] = confidences
                st.session_state['_detail_chain_ids'] = chain_ids_resolved
                st.session_state['_detail_pae_problem'] = pae_problem
                st.session_state['_detail_chain_lengths'] = chain_lengths_resolved
                st.session_state['_detail_protein_names'] = protein_names
                st.session_state['_detail_nonpolymer_ids'] = nonpolymer_ids
                st.session_state['_detail_polymer_ids'] = poly_ids
                st.session_state['_detail_cif_file'] = str(cif_file)

            # Retrieve from session state (instant — no recomputation)
            confidences = st.session_state.get('_detail_confidences')
            chain_ids_resolved = st.session_state.get('_detail_chain_ids', [])
            chain_lengths_resolved = st.session_state.get('_detail_chain_lengths', [])
            protein_names = st.session_state.get('_detail_protein_names', ['Chain A', 'Chain B'])
            nonpolymer_ids = st.session_state.get('_detail_nonpolymer_ids', [])
            cif_file = Path(st.session_state.get('_detail_cif_file', ''))
            poly_ids = st.session_state.get('_detail_polymer_ids')
            if poly_ids is None:
                # Session left over from an older app version: recompute.
                if confidences is not None:
                    from core.analyzer import polymer_chain_ids
                    poly_ids = list(polymer_chain_ids(confidences))
                else:
                    poly_ids = list(chain_ids_resolved)

            # ── Chain pair selection ─────────────────────────────────────
            # Load the per-model analysis record: it supplies the summary chain
            # pair (default selection) plus per-pair ipSAE / pair iPTM values
            # used for the option labels and the chain pair overview heatmaps.
            # Prefer per-prediction JSON (fast), fall back to the big JSON.
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

            cached_models = _load_pred_analysis(str(pred_dir), selected_pred_name)
            if cached_models is None:
                # Fall back to big JSON (cached after first load)
                all_analysis = _load_big_analysis_index(af3_folder)
                cached_models = all_analysis.get(selected_pred_name, {})

            record = (cached_models or {}).get((seed, sample), {})

            def _split_pair(chains):
                """Parse an 'A:C' chain pair into its two chain ids."""
                if isinstance(chains, (list, tuple)):
                    return [str(c) for c in chains]
                parts = str(chains or '').split(':')
                return parts if len(parts) == 2 else []

            # All unordered pairs of polymer chains, in chain order.
            candidate_pairs = []
            for _i in range(len(poly_ids)):
                for _j in range(_i + 1, len(poly_ids)):
                    candidate_pairs.append((poly_ids[_i], poly_ids[_j]))

            # Per-pair record entries, keyed by unordered chain set.
            _pairs_by_set = {}
            for _p in (record.get('pairs') or []):
                if not isinstance(_p, dict):
                    continue
                _ids = _split_pair(_p.get('chains'))
                if len(_ids) == 2:
                    _pairs_by_set[frozenset(_ids)] = _p

            def _pair_label(pair):
                """Option label, e.g. 'A:C  bait-prey  ipSAE 0.542  pair iPTM 0.64'."""
                info = _pairs_by_set.get(frozenset(pair))
                label = f"{pair[0]}:{pair[1]}"
                if info:
                    if info.get('role'):
                        label += f"  {info['role']}"
                    if info.get('ipsae') is not None:
                        label += f"  ipSAE {info['ipsae']:.3f}"
                    if info.get('pair_iptm') is not None:
                        label += f"  pair iPTM {info['pair_iptm']:.2f}"
                return label

            # Default pair: the summary pair from the record, else the first
            # pair in chain order.
            default_pair = None
            _cp_ids = _split_pair(record.get('chain_pair'))
            if len(_cp_ids) == 2:
                for _cand in candidate_pairs:
                    if set(_cand) == set(_cp_ids):
                        default_pair = list(_cand)
                        break
            if default_pair is None and candidate_pairs:
                default_pair = list(candidate_pairs[0])
            if default_pair is None:
                default_pair = list(chain_ids_resolved[:2])

            if len(candidate_pairs) > 1:
                # More than two polymer chains: let the user pick which pair
                # is analysed (defaults to the record's summary pair).
                _select_key = f"_detail_pair_select_{model_key}"
                _default_idx = 0
                for _i, _cand in enumerate(candidate_pairs):
                    if list(_cand) == default_pair:
                        _default_idx = _i
                        break
                # Reset the widget when its last value is no longer a valid
                # option (Streamlit would raise on an unknown value).
                if _select_key in st.session_state and st.session_state[_select_key] not in candidate_pairs:
                    del st.session_state[_select_key]
                _chosen = st.selectbox(
                    "Chain pair:",
                    candidate_pairs,
                    index=_default_idx,
                    format_func=_pair_label,
                    key=_select_key,
                )
                pair_ids = list(_chosen)
            else:
                pair_ids = default_pair

            _name_of = dict(zip(chain_ids_resolved, protein_names))
            pair_names = [_name_of.get(c, f"Chain {c}") for c in pair_ids]
            pair_suffix = "_".join(pair_ids) if pair_ids else "none"

            # ── Per-(model, pair) derived data ────────────────────────────
            # The 3D viewer colours only the selected pair's interface once
            # there is a choice, so viewer HTML and interface PAE are cached
            # per (model, pair); everything above stays per model.
            _pair_pae_key = f"_detail_pae_data_{model_key}_{pair_suffix}"
            _pair_viewer_key = f"_detail_viewer_html_{model_key}_{pair_suffix}"
            pae_data = st.session_state.get(_pair_pae_key)
            viewer_html = st.session_state.get(_pair_viewer_key)
            if _pair_viewer_key not in st.session_state:
                # Compute interface PAE (uses confidences already in memory — no re-read)
                pae_data = None
                if confidences is not None and chain_ids_resolved and chain_lengths_resolved and cif_file.exists():
                    # Always only the selected (protein) pair: with all chain
                    # pairs, ions and ligands (e.g. Ca2+ in EF-hands) coloured
                    # their binding residues as if they were an interface.
                    only_pair = tuple(pair_ids) if len(pair_ids) == 2 else None
                    pae_data = compute_interface_pae_per_residue(
                        cif_file, confidences, chain_ids_resolved, chain_lengths_resolved,
                        only_pair=only_pair
                    )

                # Build 3D viewer HTML ONCE (the only expensive part)
                viewer_html = None
                if cif_file.exists():
                    cif_content = cif_file.read_text()
                    model_format = "pdb" if cif_file.suffix.lower() == ".pdb" else "cif"
                    if model_format == "pdb":
                        # AF2 relaxed PDBs carry explicit hydrogens (AF3 CIFs do
                        # not); drop them so the inlined viewer payload stays small.
                        from core.af2_pae import HYDROGEN_ELEMENTS
                        cif_content = "\n".join(
                            ln for ln in cif_content.splitlines()
                            if not (ln.startswith(("ATOM", "HETATM"))
                                    and ln[76:78].strip() in HYDROGEN_ELEMENTS)
                        )
                    viewer_html = generate_viewer_html(
                        cif_content,
                        pae_residue_data=pae_data,
                        chain_names=protein_names,
                        chain_ids_all=chain_ids_resolved,
                        height=600,
                        model_format=model_format,
                        nonpolymer_chains=nonpolymer_ids,
                    )
                st.session_state[_pair_pae_key] = pae_data
                st.session_state[_pair_viewer_key] = viewer_html

            # Tabs for different analysis views
            if st.session_state.get('_detail_pae_problem'):
                st.warning(st.session_state['_detail_pae_problem'])

            tab0, tab1, tab2, tab3, tab5, tab4 = st.tabs(["3D Structure", "PAE Plots", "Interface Contacts", "Model Comparison", "PyMOL Scripts", "Export"])

            with tab0:
                st.markdown("### 3D Structure Viewer")

                if viewer_html:
                    # st.iframe (Streamlit >= 1.65) replaces the deprecated
                    # components.html; keep the old call for older Streamlit.
                    # viewer_html is generated by the app (names are escaped).
                    if hasattr(st, 'iframe'):
                        st.iframe(viewer_html, height=620)
                    else:
                        components.html(viewer_html, height=620, scrolling=False)
                else:
                    st.warning(f"CIF file not found for seed-{seed}_sample-{sample}.")

            with tab1:
                st.markdown("### PAE Matrix Visualization")

                # Chain pairs: with more than two polymer chains every pair was
                # scored per model; show ipSAE and pair iPTM as two chain x
                # chain heatmaps for a quick overview.
                _pairs_list = record.get('pairs') or []
                if len(poly_ids) > 2 and any(
                        isinstance(_p, dict) and _p.get('ipsae') is not None
                        for _p in _pairs_list):
                    st.subheader("Chain pairs")
                    _entity_of = record.get('chain_entities') or {}
                    _n = len(poly_ids)
                    _mat_ipsae = np.full((_n, _n), np.nan)
                    _mat_iptm = np.full((_n, _n), np.nan)
                    for _p in _pairs_list:
                        if not isinstance(_p, dict):
                            continue
                        _ids = _split_pair(_p.get('chains'))
                        if len(_ids) != 2:
                            continue
                        try:
                            _ia = poly_ids.index(_ids[0])
                            _ib = poly_ids.index(_ids[1])
                        except ValueError:
                            continue
                        if _p.get('ipsae') is not None:
                            _mat_ipsae[_ia, _ib] = _p['ipsae']
                            _mat_ipsae[_ib, _ia] = _p['ipsae']
                        if _p.get('pair_iptm') is not None:
                            _mat_iptm[_ia, _ib] = _p['pair_iptm']
                            _mat_iptm[_ib, _ia] = _p['pair_iptm']

                    def _axis_label(cid):
                        # Chain ID plus molecule ID from the record, e.g. 'A (1)'.
                        _ent = _entity_of.get(cid)
                        if _ent in (None, ''):
                            return cid
                        return f"{cid} ({_ent})"

                    _fig, (_ax_ipsae, _ax_iptm) = plt.subplots(
                        1, 2, figsize=(8.0, max(3.5, 0.45 * _n + 2.5)))
                    for _ax, _mat, _title in (
                            (_ax_ipsae, _mat_ipsae, 'ipSAE'),
                            (_ax_iptm, _mat_iptm, 'pair iPTM')):
                        _ax.imshow(_mat, cmap='Greens', vmin=0, vmax=1)
                        _ax.set_xticks(range(_n))
                        _ax.set_yticks(range(_n))
                        _ax.set_xticklabels([_axis_label(c) for c in poly_ids],
                                            rotation=90, fontsize=7)
                        _ax.set_yticklabels([_axis_label(c) for c in poly_ids],
                                            fontsize=7)
                        _ax.set_title(_title, fontsize=9)
                        for _i in range(_n):
                            for _j in range(_n):
                                _v = _mat[_i, _j]
                                if not np.isnan(_v):
                                    _ax.text(_j, _i, f"{_v:.2f}",
                                             ha='center', va='center', fontsize=6)
                    _fig.tight_layout()
                    st.pyplot(_fig)
                    plt.close(_fig)

                if confidences is None and _is_af2:
                    st.info(
                        "**No PAE matrix for this AF2-multimer prediction.** "
                        "PAE plots need `ranked_N_PAE.csv` files next to `ranked_N.pdb`."
                    )
                elif confidences is not None:
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
                            buf = io.BytesIO()
                            fig.savefig(buf, format='png', dpi=150, bbox_inches='tight')
                            plt.close(fig)
                            st.session_state[pae_plot_key] = buf.getvalue()
                        st.image(st.session_state[pae_plot_key])
                        pae_fname = f"{'_'.join(protein_names)}_s{seed}_m{sample}_PAE.png"
                        st.download_button("Download PAE plot", st.session_state[pae_plot_key],
                                           file_name=pae_fname, mime="image/png")

                        if len(pair_ids) >= 2 and len(chain_indices) >= 2:
                            st.subheader("Interface Zoom")
                            chain_a_name, chain_b_name = pair_names[0], pair_names[1]
                            chain_a_idx = chain_indices[pair_ids[0]]
                            chain_b_idx = chain_indices[pair_ids[1]]

                            zoom_plot_key = f"_zoom_plot_{model_key}_{pair_suffix}"
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

                from core.interface_analyzer import (
                    analyze_interface, analyze_interface_spatial, map_interface_regions,
                )

                _is_ap = _is_af2 and confidences is None   # AF2 model without PAE
                _can_full = cif_file.exists() and confidences is not None
                _can_spatial = cif_file.exists() and _is_ap

                if _is_ap and _can_spatial:
                    st.caption(
                        "AF2-multimer: spatial Cβ–Cβ contacts only (no PAE filter)."
                    )

                if _can_full or _can_spatial:
                    iface_key = f"_iface_{model_key}_{pair_suffix}"
                    if iface_key not in st.session_state:
                        if _can_full:
                            st.session_state[iface_key] = analyze_interface(
                                cif_file, confidences, pae_cutoff=10.0, distance_cutoff=8.0,
                                chain_pair=tuple(pair_ids)
                            )
                        else:
                            st.session_state[iface_key] = analyze_interface_spatial(
                                cif_file, distance_cutoff=8.0,
                                chain_pair=tuple(pair_ids)
                            )
                    interface_result = st.session_state[iface_key]

                    contacts = interface_result.get('contacts', [])
                    summary = interface_result.get('summary', {})
                    hub_residues_a = interface_result.get('hub_residues_a', [])
                    hub_residues_b = interface_result.get('hub_residues_b', [])

                    if contacts:
                        if _is_ap:
                            st.success(f"Found {len(contacts)} interface contacts (Cβ≤8Å, no PAE filter)")
                        else:
                            st.success(f"Found {len(contacts)} interface contacts (PAE≤10Å, Cβ≤8Å)")

                        # Show summary metrics — Mean PAE is unavailable for AlphaPulldown
                        mean_pae = summary.get('mean_pae')
                        mean_dist = summary.get('mean_distance')
                        if mean_pae is None:
                            col1, col2, col3 = st.columns(3)
                            col1.metric("Total Contacts", summary.get('total_contacts', 0))
                            col2.metric("Mean Distance",
                                        f"{mean_dist:.2f} Å" if mean_dist is not None else "—")
                            col3.metric("Interface Area",
                                        f"{summary.get('interface_area_estimate', 0)} residues")
                        else:
                            col1, col2, col3, col4 = st.columns(4)
                            col1.metric("Total Contacts", summary.get('total_contacts', 0))
                            col2.metric("Mean PAE", f"{mean_pae:.2f} Å")
                            col3.metric("Mean Distance",
                                        f"{mean_dist:.2f} Å" if mean_dist is not None else "—")
                            col4.metric("Interface Area",
                                        f"{summary.get('interface_area_estimate', 0)} residues")

                        import pandas as pd

                        # Show top contacts table with protein names as column headers
                        contact_df = pd.DataFrame(contacts[:50])  # Show top 50
                        # Header stems must be distinct: empty names fall back to
                        # "Chain A/B", and identical labels (e.g. homodimers) get
                        # a chain-ID suffix, so the renamed DataFrame never has
                        # duplicate column names (which st.dataframe rejects).
                        name_a = (pair_names[0] if len(pair_names) > 0 else "") or "Chain A"
                        name_b = (pair_names[1] if len(pair_names) > 1 else "") or "Chain B"
                        if name_a == name_b:
                            cid_a = pair_ids[0] if len(pair_ids) > 0 else "A"
                            cid_b = pair_ids[1] if len(pair_ids) > 1 else "B"
                            name_a = f"{name_a} ({cid_a})"
                            name_b = f"{name_b} ({cid_b})"
                        col_rename = {
                            'chain_a_res': f'{name_a} res',
                            'chain_b_res': f'{name_b} res',
                            'chain_a_aa': f'{name_a} aa',
                            'chain_b_aa': f'{name_b} aa',
                            'distance': 'Cβ–Cβ (Å)',
                            'pae': 'PAE (Å)',
                        }
                        # Residue/aa pairs side by side, then the measurements
                        order = [c for c in ('chain_a_res', 'chain_a_aa', 'chain_b_res', 'chain_b_aa',
                                             'distance', 'pae')
                                 if c in contact_df.columns]
                        contact_df = contact_df[order + [c for c in contact_df.columns if c not in order]]
                        contact_df = contact_df.rename(columns=col_rename)
                        st.dataframe(contact_df, width='stretch', hide_index=True, height=TABLE_H)

                        # Hub residues analysis
                        st.subheader("Hub Residues")
                        name_a = pair_names[0] if len(pair_names) > 0 else "Chain A"
                        name_b = pair_names[1] if len(pair_names) > 1 else "Chain B"

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

                        # Interface region mapping. Lengths of the two chosen
                        # chains, looked up by their index in chain_ids_resolved
                        # (the pair is no longer the first two chains).
                        _chain_len_by_id = dict(zip(chain_ids_resolved, chain_lengths_resolved))
                        _len_a = _chain_len_by_id.get(pair_ids[0]) if pair_ids else None
                        _len_b = _chain_len_by_id.get(pair_ids[1]) if len(pair_ids) > 1 else None
                        if _len_a is not None and _len_b is not None:
                            st.subheader("Interface Sequence Regions")
                            region_info = map_interface_regions(contacts, _len_a, _len_b)

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
                    if _is_ap:
                        st.warning("Structure file (ranked_*.pdb) not found for this model.")
                    else:
                        st.warning("Required files not found. Run full analysis first.")

            with tab3:
                st.markdown("### All Models Comparison")

                from core.analyzer import load_summary
                from core.pae_plotter import plot_all_models_comparison

                # cached_models (and the selected model's record) were loaded
                # above the tabs for the chain pair selector; reuse them here.

                # Compare all 5 models in a table with per-model metrics
                comparison_data = []
                for ss in seed_samples:
                    is_top = ss['seed'] == top_seed and ss['sample'] == top_sample
                    cache_key = (ss['seed'], ss['sample'])
                    cached = cached_models.get(cache_key, {})

                    # Per-model iPTM/PTM. AlphaPulldown stores combined iptm+ptm
                    # directly on seed_samples; AF3 needs load_summary.
                    if ss.get('iptm') is not None and _is_af2:
                        model_iptm = ss['iptm']
                        model_ptm = ss.get('ptm', model_iptm)
                    else:
                        model_summary = load_summary(pred_dir, file_prefix, ss['seed'], ss['sample'])
                        model_iptm = model_summary.get('iptm', 0) if model_summary else selected_pred['iptm']
                        model_ptm = model_summary.get('ptm', 0) if model_summary else selected_pred['ptm']

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
                st.dataframe(comp_df, width='stretch', hide_index=True)

                # Side-by-side PAE matrix comparison — AF3 only
                _show_pae_comparison = selected_pred.get('format') != 'alphapulldown'
                if not _show_pae_comparison:
                    st.caption(
                        "PAE matrix side-by-side comparison is not available for "
                        "AF2-multimer predictions without ranked_N_PAE.csv files."
                    )

                comp_plot_key = f"_comp_plot_{selected_pred_name}"
                if _show_pae_comparison and comp_plot_key not in st.session_state:
                    pae_matrices = []
                    model_labels = []
                    for ss in seed_samples:
                        is_top = ss['seed'] == top_seed and ss['sample'] == top_sample
                        label = "Top" if is_top else f"s{ss['seed']}-m{ss['sample']}"
                        conf = load_confidences(pred_dir, file_prefix, ss['seed'], ss['sample'])
                        if conf and conf.get('pae'):
                            pae_matrices.append(np.array(conf['pae']))
                            model_labels.append(label)

                    if pae_matrices:
                        fig = plot_all_models_comparison(pae_matrices, model_labels, chain_lengths_resolved, protein_names)
                        buf = io.BytesIO()
                        fig.savefig(buf, format='png', dpi=150, bbox_inches='tight')
                        plt.close(fig)
                        st.session_state[comp_plot_key] = buf.getvalue()
                    else:
                        st.session_state[comp_plot_key] = None

                if _show_pae_comparison:
                    if st.session_state.get(comp_plot_key):
                        st.image(st.session_state[comp_plot_key])
                    else:
                        st.warning("Could not load PAE matrices for comparison.")

            with tab5:
                st.markdown("### PyMOL Visualization Scripts")
                st.caption("Generate .py scripts to load in PyMOL for offline analysis")

                from core.pymol_script import generate_pymol_script

                if cif_file.exists() and len(pair_ids) >= 2:
                    # Reuse PAE data (None for AlphaPulldown — no PAE tiers in script)
                    pae_data_pymol = pae_data

                    # Reuse interface analysis from tab2 (compute spatially if AP)
                    iface_key = f"_iface_{model_key}_{pair_suffix}"
                    if iface_key not in st.session_state:
                        if confidences is not None:
                            from core.interface_analyzer import analyze_interface as _analyze_iface
                            st.session_state[iface_key] = _analyze_iface(
                                cif_file, confidences, pae_cutoff=10.0, distance_cutoff=8.0,
                                chain_pair=tuple(pair_ids)
                            )
                        else:
                            from core.interface_analyzer import analyze_interface_spatial as _analyze_iface_spatial
                            st.session_state[iface_key] = _analyze_iface_spatial(
                                cif_file, distance_cutoff=8.0,
                                chain_pair=tuple(pair_ids)
                            )
                    iface_result = st.session_state[iface_key]
                    contacts_pymol = iface_result.get('contacts', []) if iface_result else None

                    name_a = pair_names[0] if len(pair_names) > 0 else pair_ids[0]
                    name_b = pair_names[1] if len(pair_names) > 1 else pair_ids[1]

                    pymol_script = generate_pymol_script(
                        cif_path=str(cif_file.resolve()),
                        prediction_name=selected_pred_name,
                        chain_a_id=pair_ids[0],
                        chain_b_id=pair_ids[1],
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
                    st.warning("Structure file and at least two chains required for PyMOL script generation.")

            with tab4:
                st.markdown("### Export Options")

                import numpy as np

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
                elif _is_af2:
                    st.caption("PAE matrix export needs ranked_N_PAE.csv files for AF2-multimer predictions.")

                if cif_file.exists():
                    iface_key = f"_iface_{model_key}_{pair_suffix}"
                    if iface_key in st.session_state and st.session_state[iface_key]:
                        contacts = st.session_state[iface_key].get('contacts', [])
                    elif confidences is not None:
                        from core.interface_analyzer import analyze_interface as _analyze_iface_export
                        contacts = _analyze_iface_export(
                            cif_file, confidences, pae_cutoff=10.0, distance_cutoff=8.0,
                            chain_pair=tuple(pair_ids)
                        ).get('contacts', [])
                    else:
                        from core.interface_analyzer import analyze_interface_spatial
                        contacts = analyze_interface_spatial(
                            cif_file, distance_cutoff=8.0, chain_pair=tuple(pair_ids)
                        ).get('contacts', [])

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
                    st.warning("Structure file not found for this model.")
