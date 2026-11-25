"""
Generalized interface analysis module for AF3 predictions.
Based on spatial_interface_analysis.py and interface_analysis.py but generalized
(no hardcoded protein names/sequences).
"""

import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from core.analyzer import CIFParser, chain_token_indices, _structure_pair


def analyze_interface(cif_path: Path,
                      confidences: dict,
                      pae_cutoff: float = 10.0,
                      distance_cutoff: float = 8.0,
                      chain_pair: Optional[Tuple[str, str]] = None) -> dict:
    """
    Full interface analysis combining PAE and spatial filtering.

    Args:
        cif_path: Path to CIF file
        confidences: Dictionary loaded from confidences.json
        pae_cutoff: PAE cutoff for interface definition
        distance_cutoff: Cbeta-Cbeta distance cutoff
        chain_pair: (chain A, chain B) to analyse; default the first two
            polymer chains (interface_chain_pair)

    Returns dict with:
        'contacts': list of contact dicts
        'hub_residues_a': list of (res_id, res_name, contact_count)
        'hub_residues_b': same for chain B
        'summary': dict with total_contacts, mean_pae, mean_distance, interface_area_estimate
    """
    if not cif_path.exists():
        return {
            'contacts': [],
            'hub_residues_a': [],
            'hub_residues_b': [],
            'summary': {}
        }

    parser = CIFParser(str(cif_path), verbose=False)
    if not parser.parse_atoms():
        return {
            'contacts': [],
            'hub_residues_a': [],
            'hub_residues_b': [],
            'summary': {}
        }

    chain_ids = parser.get_chain_ids()
    if len(chain_ids) < 2:
        return {
            'contacts': [],
            'hub_residues_a': [],
            'hub_residues_b': [],
            'summary': {}
        }

    pae_matrix = np.array(confidences.get('pae', []))
    if pae_matrix.size == 0:
        return {
            'contacts': [],
            'hub_residues_a': [],
            'hub_residues_b': [],
            'summary': {}
        }

    if chain_pair and all(c in chain_ids for c in chain_pair):
        chain_a_id, chain_b_id = chain_pair
    else:
        chain_a_id, chain_b_id = _structure_pair(chain_ids, confidences)
    tok = chain_token_indices(confidences)
    chain_a_indices = tok.get(chain_a_id, [])
    chain_b_indices = tok.get(chain_b_id, [])
    if not chain_a_indices or not chain_b_indices:
        return {
            'contacts': [],
            'hub_residues_a': [],
            'hub_residues_b': [],
            'summary': {}
        }

    chain_a_residues = parser.get_chain_residues(chain_a_id)
    chain_b_residues = parser.get_chain_residues(chain_b_id)

    # Pre-build CB coordinate dictionary for O(1) lookup
    cb_coords = {}     # (chain_id, seq_id) -> (x, y, z)
    cb_is_cb = {}      # (chain_id, seq_id) -> True if CB, False if CA fallback
    res_names = {}
    for atom in parser.atoms:
        key = (atom['chain_id'], atom['seq_id'])
        name = atom['atom_name']
        # Prefer CB, fall back to CA for glycine
        if name == 'CB':
            cb_coords[key] = (atom['x'], atom['y'], atom['z'])
            cb_is_cb[key] = True
        elif name == 'CA' and key not in cb_coords:
            cb_coords[key] = (atom['x'], atom['y'], atom['z'])
            cb_is_cb[key] = False
        if key not in res_names:
            res_names[key] = atom['res_name']

    contacts = []

    for i, res_a in enumerate(chain_a_residues):
        # Get CB coordinate for residue A
        coord_a = cb_coords.get((chain_a_id, res_a))
        if coord_a is None:
            continue

        # Map residue index to PAE matrix index
        a_idx = chain_a_indices[i] if i < len(chain_a_indices) else -1
        if a_idx < 0 or a_idx >= len(pae_matrix):
            continue

        for j, res_b in enumerate(chain_b_residues):
            # Get CB coordinate for residue B
            coord_b = cb_coords.get((chain_b_id, res_b))
            if coord_b is None:
                continue

            # Map residue index to PAE matrix index
            b_idx = chain_b_indices[j] if j < len(chain_b_indices) else -1
            if b_idx < 0 or b_idx >= len(pae_matrix[0]):
                continue

            # PAE check first (cheap array lookup) before distance calculation
            pae_val = pae_matrix[a_idx, b_idx]
            if pae_val > pae_cutoff:
                continue

            # Distance check only for PAE-passing pairs
            dx = coord_a[0] - coord_b[0]
            dy = coord_a[1] - coord_b[1]
            dz = coord_a[2] - coord_b[2]
            dist = (dx*dx + dy*dy + dz*dz) ** 0.5
            if dist > distance_cutoff:
                continue

            # Get residue names from pre-built dict
            res_name_a = res_names.get((chain_a_id, res_a), 'UNK')
            res_name_b = res_names.get((chain_b_id, res_b), 'UNK')

            contacts.append({
                'chain_a_res': res_a,
                'chain_b_res': res_b,
                'distance': round(dist, 2),
                'pae': round(float(pae_val), 2),
                'chain_a_aa': res_name_a,
                'chain_b_aa': res_name_b,
            })

    # Find hub residues
    hub_residues_a = find_hub_residues(contacts, 'chain_a_res')
    hub_residues_b = find_hub_residues(contacts, 'chain_b_res')

    # Summary statistics
    if contacts:
        mean_pae = float(np.mean([c['pae'] for c in contacts]))
        mean_distance = float(np.mean([c['distance'] for c in contacts]))
    else:
        mean_pae = 0
        mean_distance = 0

    summary = {
        'total_contacts': len(contacts),
        'mean_pae': round(mean_pae, 2),
        'mean_distance': round(mean_distance, 2),
        'interface_area_estimate': len(set(c['chain_a_res'] for c in contacts)) +
                                   len(set(c['chain_b_res'] for c in contacts))
    }

    return {
        'contacts': contacts,
        'hub_residues_a': hub_residues_a,
        'hub_residues_b': hub_residues_b,
        'summary': summary
    }


def analyze_interface_spatial(cif_path: Path,
                               distance_cutoff: float = 8.0,
                               chain_pair: Optional[Tuple[str, str]] = None) -> dict:
    """Spatial-only interface analysis (no PAE filter).

    Used for AlphaPulldown / AF2-multimer output where the per-residue PAE
    matrix is not available. Returns the same dict shape as
    ``analyze_interface``; the ``pae`` field of each contact is ``None``.
    """
    if not cif_path.exists():
        return {'contacts': [], 'hub_residues_a': [], 'hub_residues_b': [],
                'summary': {}}

    parser = CIFParser(str(cif_path), verbose=False)
    if not parser.parse_atoms():
        return {'contacts': [], 'hub_residues_a': [], 'hub_residues_b': [],
                'summary': {}}

    chain_ids = parser.get_chain_ids()
    if len(chain_ids) < 2:
        return {'contacts': [], 'hub_residues_a': [], 'hub_residues_b': [],
                'summary': {}}

    if chain_pair and all(c in chain_ids for c in chain_pair):
        chain_a_id, chain_b_id = chain_pair
    else:
        chain_a_id, chain_b_id = chain_ids[0], chain_ids[1]
    chain_a_residues = parser.get_chain_residues(chain_a_id)
    chain_b_residues = parser.get_chain_residues(chain_b_id)

    cb_coords = {}
    res_names = {}
    for atom in parser.atoms:
        key = (atom['chain_id'], atom['seq_id'])
        if atom['atom_name'] == 'CB':
            cb_coords[key] = (atom['x'], atom['y'], atom['z'])
        elif atom['atom_name'] == 'CA' and key not in cb_coords:
            cb_coords[key] = (atom['x'], atom['y'], atom['z'])
        if key not in res_names:
            res_names[key] = atom['res_name']

    contacts = []
    for res_a in chain_a_residues:
        coord_a = cb_coords.get((chain_a_id, res_a))
        if coord_a is None:
            continue
        for res_b in chain_b_residues:
            coord_b = cb_coords.get((chain_b_id, res_b))
            if coord_b is None:
                continue
            dx = coord_a[0] - coord_b[0]
            dy = coord_a[1] - coord_b[1]
            dz = coord_a[2] - coord_b[2]
            dist = (dx*dx + dy*dy + dz*dz) ** 0.5
            if dist > distance_cutoff:
                continue
            res_name_a = res_names.get((chain_a_id, res_a), 'UNK')
            res_name_b = res_names.get((chain_b_id, res_b), 'UNK')
            contacts.append({
                'chain_a_res': res_a,
                'chain_b_res': res_b,
                'distance': round(dist, 2),
                'pae': None,
                'chain_a_aa': res_name_a,
                'chain_b_aa': res_name_b,
            })

    hub_residues_a = find_hub_residues(contacts, 'chain_a_res')
    hub_residues_b = find_hub_residues(contacts, 'chain_b_res')

    if contacts:
        mean_distance = float(np.mean([c['distance'] for c in contacts]))
    else:
        mean_distance = 0

    summary = {
        'total_contacts': len(contacts),
        'mean_pae': None,
        'mean_distance': round(mean_distance, 2),
        'interface_area_estimate': len({c['chain_a_res'] for c in contacts}) +
                                   len({c['chain_b_res'] for c in contacts}),
    }

    return {
        'contacts': contacts,
        'hub_residues_a': hub_residues_a,
        'hub_residues_b': hub_residues_b,
        'summary': summary,
    }


def find_hub_residues(contacts: List[dict], chain_key: str,
                      min_contacts: int = 3) -> List[Tuple[int, str, int]]:
    """
    Find hub residues (residues with >= min_contacts interface contacts).

    Args:
        contacts: list of contact dicts from analyze_interface
        chain_key: 'chain_a_res' or 'chain_b_res'
        min_contacts: minimum contacts to qualify as hub

    Returns:
        List of (res_id, res_name, contact_count) sorted by count descending
    """
    contact_counts = {}
    residue_names = {}

    for contact in contacts:
        res_id = contact.get(chain_key)
        if res_id is None:
            continue

        contact_counts[res_id] = contact_counts.get(res_id, 0) + 1

        # Store residue name
        if chain_key == 'chain_a_res':
            residue_names[res_id] = contact.get('chain_a_aa', 'UNK')
        else:
            residue_names[res_id] = contact.get('chain_b_aa', 'UNK')

    # Filter by minimum contacts and sort
    hubs = [
        (res_id, residue_names.get(res_id, 'UNK'), count)
        for res_id, count in contact_counts.items()
        if count >= min_contacts
    ]
    hubs.sort(key=lambda x: x[2], reverse=True)

    return hubs


def map_interface_regions(contacts: List[dict],
                          chain_a_length: int,
                          chain_b_length: int,
                          window: int = 10) -> dict:
    """
    Map interface to sequence regions using a sliding window.

    Args:
        contacts: list of contact dicts from analyze_interface
        chain_a_length: Number of residues in chain A
        chain_b_length: Number of residues in chain B
        window: Sliding window size for region detection

    Returns dict with:
        'chain_a_regions': list of (start, end) tuples
        'chain_b_regions': list of (start, end) tuples
        'chain_a_coverage': fraction of chain A at interface
        'chain_b_coverage': fraction of chain B at interface
    """
    # Get unique interface residues for each chain
    chain_a_interface = set(c['chain_a_res'] for c in contacts)
    chain_b_interface = set(c['chain_b_res'] for c in contacts)

    # Find contiguous regions using sliding window
    def find_regions(interface_residues, chain_length):
        if not interface_residues:
            return []

        sorted_residues = sorted(interface_residues)
        regions = []
        start = sorted_residues[0]
        end = sorted_residues[0]

        for res in sorted_residues[1:]:
            if res - end <= window:
                end = res
            else:
                regions.append((start, end))
                start = res
                end = res
        regions.append((start, end))

        return regions

    chain_a_regions = find_regions(chain_a_interface, chain_a_length)
    chain_b_regions = find_regions(chain_b_interface, chain_b_length)

    # Calculate coverage
    chain_a_coverage = len(chain_a_interface) / chain_a_length if chain_a_length > 0 else 0
    chain_b_coverage = len(chain_b_interface) / chain_b_length if chain_b_length > 0 else 0

    return {
        'chain_a_regions': chain_a_regions,
        'chain_b_regions': chain_b_regions,
        'chain_a_coverage': round(chain_a_coverage, 3),
        'chain_b_coverage': round(chain_b_coverage, 3)
    }