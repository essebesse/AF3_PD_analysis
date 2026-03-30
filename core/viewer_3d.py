"""
3D structure viewer for Streamlit using 3Dmol.js (served as static file).

Embeds an interactive 3D viewer via st.components.v1.html() with
PAE-based interface coloring. Interface residues are identified using both
inter-chain PAE and spatial proximity.

3Dmol.js is loaded from CDN (jsdelivr),
keeping the HTML small (~5KB + CIF) for fast model switching.
"""

import json as _json
import numpy as np
import gemmi
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def compute_interface_pae_per_residue(
    cif_path: Path,
    confidences: Dict,
    chain_ids: List[str],
    chain_lengths: List[int],
    pae_cutoff: float = 12.0,
    distance_cutoff: float = 8.0,
) -> Dict[str, List[Tuple[int, float]]]:
    """
    Identify interface residues using both inter-chain PAE and spatial proximity.

    Strategy (PAE-first for speed):
    1. Extract inter-chain PAE blocks and find residue pairs below pae_cutoff
    2. Only for those pairs, check CB-CB distance in the 3D structure
    3. Color residues that pass both filters

    Returns dict mapping chain_id -> [(residue_number, best_pae), ...]
    Only interface residues are included.
    """
    pae_matrix = np.array(confidences.get('pae', []))
    if pae_matrix.size == 0 or len(chain_ids) < 2:
        return {}

    # Build PAE index boundaries
    boundaries = [0]
    for length in chain_lengths:
        boundaries.append(boundaries[-1] + length)

    # Step 1: Find inter-chain residue pairs with low PAE (vectorized)
    candidates = []
    for i in range(len(chain_ids)):
        for j in range(i + 1, len(chain_ids)):
            si, ei = boundaries[i], boundaries[i + 1]
            sj, ej = boundaries[j], boundaries[j + 1]

            pae_ij = pae_matrix[si:ei, sj:ej]
            pae_ji = pae_matrix[sj:ej, si:ei]
            bidir = (pae_ij + pae_ji.T) / 2.0

            ri_arr, rj_arr = np.where(bidir < pae_cutoff)
            if len(ri_arr) == 0:
                continue

            pae_vals = bidir[ri_arr, rj_arr]
            for k in range(len(ri_arr)):
                candidates.append((i, int(ri_arr[k]), j, int(rj_arr[k]), float(pae_vals[k])))

    if not candidates:
        return {}

    # Step 2: Parse CIF and build coordinate arrays per chain
    try:
        structure = gemmi.read_structure(str(cif_path))
    except Exception:
        return {}

    model = structure[0]

    chain_coords: Dict[str, np.ndarray] = {}
    for ci, chain_id in enumerate(chain_ids):
        n_res = chain_lengths[ci]
        arr = np.full((n_res, 3), np.nan)
        for chain in model:
            if chain.name != chain_id:
                continue
            for residue in chain:
                seq_idx = residue.seqid.num - 1
                if seq_idx < 0 or seq_idx >= n_res:
                    continue
                atom = residue.find_atom('CB', '\0') or residue.find_atom('CA', '\0')
                if atom:
                    arr[seq_idx] = [atom.pos.x, atom.pos.y, atom.pos.z]
            break
        chain_coords[chain_id] = arr

    # Step 3: Filter candidates by spatial distance
    dist_sq_cutoff = distance_cutoff * distance_cutoff
    interface_best_pae: Dict[Tuple[str, int], float] = {}

    for ci, ri, cj, rj, pae_val in candidates:
        chain_i = chain_ids[ci]
        chain_j = chain_ids[cj]
        coord_i = chain_coords[chain_i][ri]
        coord_j = chain_coords[chain_j][rj]

        if np.isnan(coord_i[0]) or np.isnan(coord_j[0]):
            continue

        dx = coord_i[0] - coord_j[0]
        dy = coord_i[1] - coord_j[1]
        dz = coord_i[2] - coord_j[2]
        dist_sq = dx * dx + dy * dy + dz * dz

        if dist_sq > dist_sq_cutoff:
            continue

        res_i = ri + 1
        res_j = rj + 1
        key_i = (chain_i, res_i)
        key_j = (chain_j, res_j)
        if key_i not in interface_best_pae or pae_val < interface_best_pae[key_i]:
            interface_best_pae[key_i] = pae_val
        if key_j not in interface_best_pae or pae_val < interface_best_pae[key_j]:
            interface_best_pae[key_j] = pae_val

    result: Dict[str, List[Tuple[int, float]]] = {}
    for (chain_id, res_num), best_pae in sorted(interface_best_pae.items()):
        if chain_id not in result:
            result[chain_id] = []
        result[chain_id].append((res_num, best_pae))

    return result


def _pae_to_hex(pae_value: float) -> str:
    """Map interface PAE to hex color."""
    if pae_value < 3.0:
        return '#228B22'
    elif pae_value < 5.0:
        return '#00CC00'
    elif pae_value < 8.0:
        return '#FFDC00'
    return '#FF7800'


CHAIN_COLORS = [
    '#8CB4DC',  # Muted blue   - chain A
    '#DC9696',  # Muted salmon - chain B
    '#A0C8A0',  # Muted green  - chain C
    '#C8B4D2',  # Muted purple - chain D
    '#D2C396',  # Muted tan    - chain E
    '#B4D2D2',  # Muted teal   - chain F
]


def generate_viewer_html(
    cif_content: str,
    pae_residue_data: Optional[Dict[str, List[Tuple[int, float]]]] = None,
    chain_names: Optional[List[str]] = None,
    chain_ids_all: Optional[List[str]] = None,
    height: int = 600,
) -> str:
    """
    Generate HTML with 3Dmol.js viewer.

    3Dmol.js is loaded from Streamlit static file server (not inline).
    PAE interface coloring is always applied when data is available.
    """
    all_chains = chain_ids_all or (list(pae_residue_data.keys()) if pae_residue_data else [])

    # Build style commands
    style_lines = []
    for idx, chain_id in enumerate(all_chains):
        base_color = CHAIN_COLORS[idx % len(CHAIN_COLORS)]
        style_lines.append(
            f'    viewer.setStyle({_json.dumps({"chain": chain_id})},'
            f'{_json.dumps({"cartoon": {"color": base_color}})});'
        )
    if not all_chains:
        style_lines.append(
            '    viewer.setStyle({},{"cartoon":{"colorscheme":"chain"}});'
        )

    # Add PAE overlay
    has_pae = pae_residue_data and any(pae_residue_data.values())
    n_interface = 0
    if has_pae:
        for chain_id, residue_paes in pae_residue_data.items():
            for res_num, pae_val in residue_paes:
                color = _pae_to_hex(pae_val)
                style_lines.append(
                    f'    viewer.setStyle({_json.dumps({"chain": chain_id, "resi": res_num})},'
                    f'{_json.dumps({"cartoon": {"color": color}})});'
                )
                n_interface += 1

    style_js = '\n'.join(style_lines)
    cif_escaped = _json.dumps(cif_content)

    # Build legend
    chain_legend_items = ''
    for i, cid in enumerate(all_chains):
        c = CHAIN_COLORS[i % len(CHAIN_COLORS)]
        label = chain_names[i] if chain_names and i < len(chain_names) else f'Chain {cid}'
        chain_legend_items += (
            f'<div style="display:flex;align-items:center;gap:6px;margin:2px 0;">'
            f'<div style="width:14px;height:14px;border-radius:3px;border:1px solid #aaa;'
            f'background:{c};flex-shrink:0;"></div>'
            f'<span>{label}</span></div>'
        )

    legend_html = ''
    if has_pae:
        legend_html = f'''
    <div style="position:absolute;bottom:8px;left:8px;background:rgba(255,255,255,0.92);
        border-radius:6px;padding:6px 10px;font-family:Arial,sans-serif;font-size:12px;
        z-index:1000;box-shadow:0 1px 4px rgba(0,0,0,0.15);">
        <div style="font-weight:bold;margin-bottom:4px;">Interface PAE ({n_interface} residues)</div>
        <div style="display:flex;align-items:center;gap:6px;margin:2px 0;">
            <div style="width:14px;height:14px;border-radius:3px;border:1px solid #aaa;
                 background:#228B22;flex-shrink:0;"></div>
            <span>&lt; 3 &#197; (very high)</span></div>
        <div style="display:flex;align-items:center;gap:6px;margin:2px 0;">
            <div style="width:14px;height:14px;border-radius:3px;border:1px solid #aaa;
                 background:#00CC00;flex-shrink:0;"></div>
            <span>3 - 5 &#197; (high)</span></div>
        <div style="display:flex;align-items:center;gap:6px;margin:2px 0;">
            <div style="width:14px;height:14px;border-radius:3px;border:1px solid #aaa;
                 background:#FFDC00;flex-shrink:0;"></div>
            <span>5 - 8 &#197; (moderate)</span></div>
        <div style="display:flex;align-items:center;gap:6px;margin:2px 0;">
            <div style="width:14px;height:14px;border-radius:3px;border:1px solid #aaa;
                 background:#FF7800;flex-shrink:0;"></div>
            <span>8 - 12 &#197; (low)</span></div>
        <div style="font-weight:bold;margin:4px 0 2px 0;">Chains</div>
        {chain_legend_items}
    </div>'''
    elif chain_legend_items:
        legend_html = f'''
    <div style="position:absolute;bottom:8px;left:8px;background:rgba(255,255,255,0.92);
        border-radius:6px;padding:6px 10px;font-family:Arial,sans-serif;font-size:12px;
        z-index:1000;box-shadow:0 1px 4px rgba(0,0,0,0.15);">
        <div style="font-weight:bold;margin-bottom:4px;">Chains</div>
        {chain_legend_items}
    </div>'''

    # No-interface message overlay
    no_iface_msg = ''
    if not has_pae and pae_residue_data is not None:
        no_iface_msg = '''
    <div style="position:absolute;top:8px;left:8px;background:rgba(255,240,220,0.95);
        border-radius:6px;padding:8px 14px;font-family:Arial,sans-serif;font-size:13px;
        z-index:1000;box-shadow:0 1px 4px rgba(0,0,0,0.15);border:1px solid #e0c080;">
        No confident interface residues detected (PAE &lt; 12 &#197; + distance &lt; 8 &#197;)
    </div>'''

    html = f'''<!DOCTYPE html>
<html>
<head>
<style>* {{ margin:0; padding:0; box-sizing:border-box; }}</style>
</head>
<body>
<div id="container" style="width:100%;height:{height}px;position:relative;">
    <div id="viewer" style="width:100%;height:100%;">
        <div id="viewer_fallback" style="padding:40px;text-align:center;color:#b91c1c;font-family:Arial,sans-serif;">
            3Dmol.js failed to load.<br>
            This usually means you are on an isolated network without CDN access.<br>
            The PAE plots and other analysis tabs still work.
        </div>
    </div>
    {legend_html}
    {no_iface_msg}
</div>
<script src="https://cdn.jsdelivr.net/npm/3dmol@2.5.4/build/3Dmol-min.js"></script>
<script>
    if (typeof $3Dmol !== "undefined") {{
        document.getElementById("viewer_fallback").remove();
        var viewer = $3Dmol.createViewer("viewer", {{backgroundColor: "white"}});
        viewer.addModel({cif_escaped}, "cif");
{style_js}
        viewer.zoomTo();
        viewer.render();
    }}
</script>
</body>
</html>'''

    return html
