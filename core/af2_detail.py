"""
AF2 detailed-analysis adapters.

Bridges AF2 inputs (ranked_N.pdb + ranked_N_PAE.csv + ranked_N_pLDDT.csv) to the
AF3-shaped data the existing detail tabs (PAE plots, interface contacts, 3D
viewer, PyMOL export) expect. The strategy is to synthesize a "confidences"
dict identical in shape to AF3's confidences.json so pae_plotter.py and the
PAE-overlay logic can be reused unchanged.

Pairwise-only: when more than two chains are present, only chains A and B
(first two in PDB order) are scored.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

try:
    import gemmi
except ImportError as e:  # pragma: no cover
    raise ImportError(
        "AF2 detail support requires gemmi. Install it with `pip install gemmi`."
    ) from e

from core.interface_analyzer import (
    classify_interaction,
    find_hub_residues,
)


# Three-letter to one-letter conversion shared with interface_analyzer
_THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}


def _representative_atom(residue: "gemmi.Residue") -> Optional["gemmi.Atom"]:
    """Return the CB atom (or CA fallback) for a residue."""
    cb = None
    ca = None
    for atom in residue:
        if atom.name == "CB" and cb is None:
            cb = atom
        elif atom.name == "CA" and ca is None:
            ca = atom
    return cb or ca


def _walk_residues(structure: "gemmi.Structure"):
    """Yield (chain_name, residue, flat_idx) over the first model in PDB order,
    skipping water/HOH entries to keep the residue count aligned with the PAE
    matrix (which AF2 produces over real residues only)."""
    flat_idx = 0
    model = structure[0]
    for chain in model:
        for residue in chain:
            if residue.het_flag == "H" and residue.name in {"HOH", "WAT"}:
                continue
            yield chain.name, residue, flat_idx
            flat_idx += 1


def load_af2_model_data(
    work_dir: Path,
    model_idx: int,
) -> Optional[Dict]:
    """
    Load all data needed to populate the AF3-shaped detail-tab pipeline for a
    single AF2 ranked_N model.

    Returns:
        {
            'confidences': {                # shaped like AF3 confidences.json
                'pae': list[list[float]],   # L x L
                'token_chain_ids': list[str],
                'token_res_ids': list[int],
                'atom_plddts': list[float], # one per residue (broadcast for AF2)
            },
            'pdb_path': Path,
            'pdb_content': str,
            'chain_ids': list[str],         # PDB chain order
            'chain_lengths': list[int],     # residues per chain
        }
        or None if required files are missing.
    """
    work_dir = Path(work_dir)
    pdb_file = work_dir / f"ranked_{model_idx}.pdb"
    pae_file = work_dir / f"ranked_{model_idx}_PAE.csv"
    plddt_file = work_dir / f"ranked_{model_idx}_pLDDT.csv"

    if not (pdb_file.exists() and pae_file.exists()):
        return None

    try:
        pae_matrix = pd.read_csv(pae_file, index_col=0).values.astype(float)
    except Exception as exc:
        print(f"AF2 detail: failed to read {pae_file.name}: {exc}")
        return None
    if pae_matrix.size == 0:
        return None

    try:
        structure = gemmi.read_structure(str(pdb_file))
    except Exception as exc:
        print(f"AF2 detail: failed to read {pdb_file.name}: {exc}")
        return None

    plddt_per_residue: List[float] = []
    if plddt_file.exists():
        try:
            df = pd.read_csv(plddt_file)
            if "pLDDT" in df.columns:
                plddt_per_residue = [float(v) for v in df["pLDDT"].values]
            else:
                plddt_per_residue = [float(v) for v in df.iloc[:, -1].values]
        except Exception as exc:
            print(f"AF2 detail: failed to read {plddt_file.name}: {exc}")

    pae_size = pae_matrix.shape[0]
    token_chain_ids: List[str] = []
    token_res_ids: List[int] = []
    chain_lengths_map: Dict[str, int] = {}
    chain_order: List[str] = []

    for chain_name, residue, flat_idx in _walk_residues(structure):
        if flat_idx >= pae_size:
            break
        token_chain_ids.append(chain_name)
        token_res_ids.append(int(residue.seqid.num))
        if chain_name not in chain_lengths_map:
            chain_lengths_map[chain_name] = 0
            chain_order.append(chain_name)
        chain_lengths_map[chain_name] += 1

    # AF2 pLDDT is per-residue; if length mismatch, pad/truncate so downstream
    # consumers still work without IndexError. Preferring the residue ordering
    # that matches token_chain_ids.
    if plddt_per_residue and len(plddt_per_residue) != len(token_chain_ids):
        # Truncate or pad to match
        if len(plddt_per_residue) > len(token_chain_ids):
            plddt_per_residue = plddt_per_residue[: len(token_chain_ids)]
        else:
            plddt_per_residue = plddt_per_residue + [0.0] * (
                len(token_chain_ids) - len(plddt_per_residue)
            )

    confidences = {
        "pae": pae_matrix.tolist(),
        "token_chain_ids": token_chain_ids,
        "token_res_ids": token_res_ids,
        "atom_plddts": plddt_per_residue,  # one entry per residue (AF2 has no per-atom pLDDT)
    }

    chain_ids = chain_order
    chain_lengths = [chain_lengths_map[c] for c in chain_order]

    # 3Dmol renders cartoons + heavy-atom sticks; embedded H atoms inflate the
    # iframe srcdoc by 2-3x (AF2 ColabFold PDBs have explicit H, AF3 CIFs do
    # not). Strip them so the viewer payload stays comparable to AF3.
    pdb_content = "\n".join(
        line for line in pdb_file.read_text().splitlines()
        if not (line.startswith(("ATOM", "HETATM")) and len(line) >= 78
                and line[76:78].strip() == "H")
    )

    return {
        "confidences": confidences,
        "pdb_path": pdb_file,
        "pdb_content": pdb_content,
        "chain_ids": chain_ids,
        "chain_lengths": chain_lengths,
    }


def _build_residue_index(
    pdb_path: Path,
) -> Tuple[
    List[str],                                  # chain_ids in PDB order
    Dict[Tuple[str, int], Tuple[float, float, float]],  # (chain, resnum) -> (x,y,z)
    Dict[Tuple[str, int], str],                 # (chain, resnum) -> three-letter aa
    Dict[str, List[int]],                       # chain -> ordered residue numbers
]:
    """Single gemmi pass that builds everything the interface analysis needs."""
    structure = gemmi.read_structure(str(pdb_path))
    chain_ids: List[str] = []
    coords: Dict[Tuple[str, int], Tuple[float, float, float]] = {}
    res_names: Dict[Tuple[str, int], str] = {}
    chain_residues: Dict[str, List[int]] = {}

    model = structure[0]
    for chain in model:
        cid = chain.name
        if cid not in chain_residues:
            chain_ids.append(cid)
            chain_residues[cid] = []
        for residue in chain:
            if residue.het_flag == "H" and residue.name in {"HOH", "WAT"}:
                continue
            resnum = int(residue.seqid.num)
            atom = _representative_atom(residue)
            if atom is None:
                continue
            key = (cid, resnum)
            coords[key] = (float(atom.pos.x), float(atom.pos.y), float(atom.pos.z))
            res_names[key] = residue.name
            chain_residues[cid].append(resnum)

    return chain_ids, coords, res_names, chain_residues


def analyze_af2_interface(
    pdb_path: Path,
    confidences: Dict,
    pae_cutoff: float = 10.0,
    distance_cutoff: float = 8.0,
) -> Dict:
    """
    Full interface analysis for an AF2 model — mirror of
    interface_analyzer.analyze_interface but PDB-based via gemmi.

    Returns the same dict shape as the AF3 helper:
      contacts, hub_residues_a, hub_residues_b, chemical_breakdown, summary
    """
    empty = {
        "contacts": [],
        "hub_residues_a": [],
        "hub_residues_b": [],
        "chemical_breakdown": {},
        "summary": {},
    }

    pdb_path = Path(pdb_path)
    if not pdb_path.exists():
        return empty

    pae_matrix = np.array(confidences.get("pae", []))
    if pae_matrix.size == 0:
        return empty

    token_chain_ids = confidences.get("token_chain_ids", [])
    token_res_ids = confidences.get("token_res_ids", [])
    if not token_chain_ids or not token_res_ids:
        return empty

    # Build (chain,res) -> PAE-token-index lookup from the synthesized confidences
    res_to_idx: Dict[Tuple[str, int], int] = {}
    for i, (cid, rnum) in enumerate(zip(token_chain_ids, token_res_ids)):
        res_to_idx[(cid, int(rnum))] = i

    chain_ids, coords, res_names, chain_residues = _build_residue_index(pdb_path)
    if len(chain_ids) < 2:
        return empty

    chain_a_id, chain_b_id = chain_ids[0], chain_ids[1]
    a_residues = chain_residues[chain_a_id]
    b_residues = chain_residues[chain_b_id]

    pae_rows, pae_cols = pae_matrix.shape
    contacts: List[Dict] = []

    for res_a in a_residues:
        coord_a = coords.get((chain_a_id, res_a))
        a_idx = res_to_idx.get((chain_a_id, res_a), -1)
        if coord_a is None or a_idx < 0 or a_idx >= pae_rows:
            continue
        for res_b in b_residues:
            coord_b = coords.get((chain_b_id, res_b))
            b_idx = res_to_idx.get((chain_b_id, res_b), -1)
            if coord_b is None or b_idx < 0 or b_idx >= pae_cols:
                continue

            pae_val = float(pae_matrix[a_idx, b_idx])
            if pae_val > pae_cutoff:
                continue

            dx = coord_a[0] - coord_b[0]
            dy = coord_a[1] - coord_b[1]
            dz = coord_a[2] - coord_b[2]
            dist = (dx * dx + dy * dy + dz * dz) ** 0.5
            if dist > distance_cutoff:
                continue

            res_name_a = res_names.get((chain_a_id, res_a), "UNK")
            res_name_b = res_names.get((chain_b_id, res_b), "UNK")

            contacts.append({
                "chain_a_res": res_a,
                "chain_b_res": res_b,
                "distance": round(dist, 2),
                "pae": round(pae_val, 2),
                "chain_a_aa": res_name_a,
                "chain_b_aa": res_name_b,
                "interaction_type": classify_interaction(res_name_a, res_name_b, dist),
            })

    hub_residues_a = find_hub_residues(contacts, "chain_a_res")
    hub_residues_b = find_hub_residues(contacts, "chain_b_res")

    chemical_breakdown: Dict[str, int] = {}
    for c in contacts:
        itype = c.get("interaction_type", "other")
        chemical_breakdown[itype] = chemical_breakdown.get(itype, 0) + 1

    if contacts:
        mean_pae = float(np.mean([c["pae"] for c in contacts]))
        mean_distance = float(np.mean([c["distance"] for c in contacts]))
    else:
        mean_pae = 0.0
        mean_distance = 0.0

    summary = {
        "total_contacts": len(contacts),
        "mean_pae": round(mean_pae, 2),
        "mean_distance": round(mean_distance, 2),
        "interface_area_estimate": (
            len({c["chain_a_res"] for c in contacts})
            + len({c["chain_b_res"] for c in contacts})
        ),
    }

    return {
        "contacts": contacts,
        "hub_residues_a": hub_residues_a,
        "hub_residues_b": hub_residues_b,
        "chemical_breakdown": chemical_breakdown,
        "summary": summary,
    }


def compute_af2_interface_pae_per_residue(
    pdb_path: Path,
    confidences: Dict,
    chain_ids: List[str],
    chain_lengths: List[int],
    pae_cutoff: float = 12.0,
    distance_cutoff: float = 8.0,
) -> Dict[str, List[Tuple[int, float]]]:
    """
    PAE-first interface residue identification for AF2 — mirror of
    viewer_3d.compute_interface_pae_per_residue but using residue numbers
    drawn from the PDB rather than 1..N indexing.

    Returns dict mapping chain_id -> [(residue_number, best_pae), ...].
    """
    pae_matrix = np.array(confidences.get("pae", []))
    if pae_matrix.size == 0 or len(chain_ids) < 2:
        return {}

    token_chain_ids = confidences.get("token_chain_ids", [])
    token_res_ids = confidences.get("token_res_ids", [])
    if not token_chain_ids or not token_res_ids:
        return {}

    # Group token indices by chain, preserving order
    chain_token_indices: Dict[str, List[int]] = {}
    chain_token_resnums: Dict[str, List[int]] = {}
    for i, (cid, rnum) in enumerate(zip(token_chain_ids, token_res_ids)):
        chain_token_indices.setdefault(cid, []).append(i)
        chain_token_resnums.setdefault(cid, []).append(int(rnum))

    # Find low-PAE inter-chain candidate pairs (only the chains we care about)
    candidates: List[Tuple[str, int, str, int, float]] = []
    for ci in range(len(chain_ids)):
        for cj in range(ci + 1, len(chain_ids)):
            chain_i, chain_j = chain_ids[ci], chain_ids[cj]
            idx_i = chain_token_indices.get(chain_i, [])
            idx_j = chain_token_indices.get(chain_j, [])
            res_i = chain_token_resnums.get(chain_i, [])
            res_j = chain_token_resnums.get(chain_j, [])
            if not idx_i or not idx_j:
                continue
            sub = pae_matrix[np.ix_(idx_i, idx_j)]
            sub_t = pae_matrix[np.ix_(idx_j, idx_i)]
            bidir = (sub + sub_t.T) / 2.0
            ri_arr, rj_arr = np.where(bidir < pae_cutoff)
            if len(ri_arr) == 0:
                continue
            pae_vals = bidir[ri_arr, rj_arr]
            for k in range(len(ri_arr)):
                rnum_i = res_i[int(ri_arr[k])]
                rnum_j = res_j[int(rj_arr[k])]
                candidates.append((chain_i, rnum_i, chain_j, rnum_j, float(pae_vals[k])))

    if not candidates:
        return {}

    # Distance filter — load coords once
    _, coords, _, _ = _build_residue_index(Path(pdb_path))

    dist_sq_cutoff = distance_cutoff * distance_cutoff
    interface_best_pae: Dict[Tuple[str, int], float] = {}
    for chain_i, rnum_i, chain_j, rnum_j, pae_val in candidates:
        coord_i = coords.get((chain_i, rnum_i))
        coord_j = coords.get((chain_j, rnum_j))
        if coord_i is None or coord_j is None:
            continue
        dx = coord_i[0] - coord_j[0]
        dy = coord_i[1] - coord_j[1]
        dz = coord_i[2] - coord_j[2]
        if dx * dx + dy * dy + dz * dz > dist_sq_cutoff:
            continue
        key_i = (chain_i, rnum_i)
        key_j = (chain_j, rnum_j)
        if key_i not in interface_best_pae or pae_val < interface_best_pae[key_i]:
            interface_best_pae[key_i] = pae_val
        if key_j not in interface_best_pae or pae_val < interface_best_pae[key_j]:
            interface_best_pae[key_j] = pae_val

    result: Dict[str, List[Tuple[int, float]]] = {}
    for (chain_id, res_num), best_pae in sorted(interface_best_pae.items()):
        result.setdefault(chain_id, []).append((res_num, best_pae))
    return result


def map_af2_interface_regions(
    contacts: List[Dict],
    chain_a_length: int,
    chain_b_length: int,
    window: int = 10,
) -> Dict:
    """Re-export of map_interface_regions for symmetry with the AF3 path —
    uses the existing chain-length-based logic (which is residue-count, not
    PDB-numbering, so works identically for AF2 once you pass the right L)."""
    from core.interface_analyzer import map_interface_regions
    return map_interface_regions(contacts, chain_a_length, chain_b_length, window)
