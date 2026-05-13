"""
Core Analysis Module for AF3 Pulldown Analysis App

Extracts and adapts core analysis functions from AF3_PD_analysis_v4.py
to analyze ALL models (not just top-ranked) in AF3 prediction directories.
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
import re
import gemmi


def calculate_asymmetric_ipsae(pae_matrix: np.ndarray,
                               aligned_chain_indices: List[int],
                               scored_chain_indices: List[int],
                               pae_cutoff: float = 10.0) -> float:
    """
    Calculate asymmetric ipSAE score (A -> B direction).
    From Dunbrack 2025 (Equation 14).
    """
    per_residue_ipsae_scores = []

    for i in aligned_chain_indices:
        pae_vector = pae_matrix[i, scored_chain_indices]
        filtered_paes = pae_vector[pae_vector < pae_cutoff]

        if len(filtered_paes) == 0:
            continue

        L = len(filtered_paes)
        if L >= 27:
            d_o = 1.24 * np.power(L - 15, 1.0/3.0) - 1.8
        else:
            d_o = 1.0

        terms = 1.0 / (1.0 + np.power(filtered_paes / d_o, 2))
        mean_score = np.mean(terms)
        per_residue_ipsae_scores.append(mean_score)

    if len(per_residue_ipsae_scores) == 0:
        return 0.0

    return float(np.max(per_residue_ipsae_scores))


def calculate_ipsae(pae_matrix: np.ndarray,
                   chain_a_indices: List[int],
                   chain_b_indices: List[int],
                   pae_cutoff: float = 10.0) -> float:
    """
    Calculate symmetric ipSAE score between two chains.
    ipSAE(A, B) = max[ipSAE(A -> B), ipSAE(B -> A)]
    """
    ipsae_A_to_B = calculate_asymmetric_ipsae(
        pae_matrix, chain_a_indices, chain_b_indices, pae_cutoff
    )
    ipsae_B_to_A = calculate_asymmetric_ipsae(
        pae_matrix, chain_b_indices, chain_a_indices, pae_cutoff
    )
    return max(ipsae_A_to_B, ipsae_B_to_A)


class CIFParser:
    """mmCIF parser used by the analysis pipeline.

    Backed by gemmi (handles edge cases the legacy hand-rolled tokenizer
    missed: multi-character atom/residue names, quoted strings, alternate
    locations, varying column orderings). The ``atoms`` list preserves the
    same shape as before — ordered by CIF appearance — so downstream code
    that relies on positional ``atom_plddts`` indexing keeps working.
    """

    def __init__(self, cif_file: str, verbose: bool = False):
        self.cif_file = cif_file
        self.atoms = []
        self.chain_mapping = {}
        self.verbose = verbose

    def parse_atoms(self) -> bool:
        """Parse atomic coordinates from CIF file."""
        try:
            structure = gemmi.read_structure(self.cif_file)
        except Exception as e:
            if self.verbose:
                print(f"Error parsing CIF file {self.cif_file}: {e}")
            return False

        if len(structure) == 0:
            return False

        # Use the first model; AF3 always emits a single model per CIF.
        model = structure[0]
        for chain in model:
            cid = chain.name
            for residue in chain:
                seq_id = residue.seqid.num
                res_name = residue.name
                for atom in residue:
                    self.atoms.append({
                        'seq_id': seq_id,
                        'chain_id': cid,
                        'atom_name': atom.name,
                        'res_name': res_name,
                        'x': atom.pos.x,
                        'y': atom.pos.y,
                        'z': atom.pos.z,
                    })

        if self.verbose:
            print(f"Parsed {len(self.atoms)} atoms from {self.cif_file}")
        return len(self.atoms) > 0

    def get_residue_atoms(self, chain_id: str, seq_id: int) -> List[Dict]:
        """Get all atoms for a specific residue."""
        return [atom for atom in self.atoms
                if atom['chain_id'] == chain_id and atom['seq_id'] == seq_id]

    def get_representative_atom(self, chain_id: str, seq_id: int) -> Optional[Dict]:
        """Get CB atom (or CA for glycine) as representative of residue position."""
        residue_atoms = self.get_residue_atoms(chain_id, seq_id)
        if not residue_atoms:
            return None

        # Try CB first (best representative including side chain direction)
        for atom in residue_atoms:
            if atom['atom_name'] == 'CB':
                return atom

        # Fallback to CA (for glycine or if CB missing)
        for atom in residue_atoms:
            if atom['atom_name'] == 'CA':
                return atom

        # Last resort: first atom
        return residue_atoms[0]

    def calculate_cb_distance(self, chain1: str, seq1: int, chain2: str, seq2: int) -> float:
        """Calculate CB-CB distance (CA for glycine)."""
        atom1 = self.get_representative_atom(chain1, seq1)
        atom2 = self.get_representative_atom(chain2, seq2)
        if atom1 is None or atom2 is None:
            return float('inf')
        return float(np.sqrt(
            (atom1['x'] - atom2['x'])**2 +
            (atom1['y'] - atom2['y'])**2 +
            (atom1['z'] - atom2['z'])**2
        ))

    def get_chain_ids(self) -> List[str]:
        """Get unique chain IDs in order."""
        seen = []
        for atom in self.atoms:
            if atom['chain_id'] not in seen:
                seen.append(atom['chain_id'])
        return seen

    def get_chain_residues(self, chain_id: str) -> List[int]:
        """Get sorted list of residue seq_ids for a chain."""
        res_ids = set()
        for atom in self.atoms:
            if atom['chain_id'] == chain_id:
                res_ids.add(atom['seq_id'])
        return sorted(res_ids)

    def get_chain_lengths(self) -> Dict[str, int]:
        """Get the number of residues per chain."""
        return {chain: len(self.get_chain_residues(chain)) for chain in self.get_chain_ids()}


def load_confidences(pred_dir: Path, prefix: str, seed: Optional[int] = None,
                     sample: Optional[int] = None) -> Optional[Dict]:
    """
    Load confidences.json from a prediction directory.
    Handles local-pipeline naming (with seed-N_sample-M subdirs), the
    AF3 Server flat layout ({prefix}_full_data_{sample}.json with no
    seed dirs), and the top-level fallback when only the top model is stored.
    """
    candidates = []

    if seed is not None and sample is not None:
        # Local pipeline — sample-specific files
        candidates.append(pred_dir / f"{prefix}_seed-{seed}_sample-{sample}_confidences.json")
        sample_dir = pred_dir / f"seed-{seed}_sample-{sample}"
        candidates.append(sample_dir / "confidences.json")
        candidates.append(sample_dir / f"{prefix}_seed-{seed}_sample-{sample}_confidences.json")
        # AF3 Server flat — full_data_<sample>.json holds the same PAE/atom_plddts payload
        candidates.append(pred_dir / f"{prefix}_full_data_{sample}.json")

    # Always fall back to top-level files
    candidates.append(pred_dir / f"{prefix}_confidences.json")
    candidates.append(pred_dir / "confidences.json")

    for conf_file in candidates:
        if conf_file.exists():
            with open(conf_file, 'r') as f:
                return json.load(f)
    return None


def load_summary(pred_dir: Path, prefix: str, seed: Optional[int] = None,
                 sample: Optional[int] = None) -> Optional[Dict]:
    """
    Load summary_confidences.json from a prediction directory.
    Handles local-pipeline + AF3 Server flat layouts and falls back to
    top-level files when only the top model is stored.
    """
    candidates = []

    if seed is not None and sample is not None:
        candidates.append(pred_dir / f"{prefix}_seed-{seed}_sample-{sample}_summary_confidences.json")
        sample_dir = pred_dir / f"seed-{seed}_sample-{sample}"
        candidates.append(sample_dir / "summary_confidences.json")
        candidates.append(sample_dir / f"{prefix}_seed-{seed}_sample-{sample}_summary_confidences.json")
        # AF3 Server flat
        candidates.append(pred_dir / f"{prefix}_summary_confidences_{sample}.json")

    candidates.append(pred_dir / f"{prefix}_summary_confidences.json")
    candidates.append(pred_dir / "summary_confidences.json")

    for summary_file in candidates:
        if summary_file.exists():
            with open(summary_file, 'r') as f:
                return json.load(f)
    return None


def extract_chain_indices(confidences: Dict) -> Tuple[List[int], List[int], List[int], Dict[str, List[int]]]:
    """
    Extract chain indices from confidences.json.

    Returns:
        Tuple of (all_indices, chain_a_indices, chain_b_indices, chain_indices_dict)
    """
    token_chain_ids = confidences.get('token_chain_ids', [])

    # Preserve order of first appearance (not set() which randomizes)
    unique_chains = []
    for c in token_chain_ids:
        if c not in unique_chains:
            unique_chains.append(c)

    chain_indices = {chain: [] for chain in unique_chains}
    all_indices = []

    for i, chain in enumerate(token_chain_ids):
        all_indices.append(i)
        chain_indices[chain].append(i)

    chain_a = chain_indices.get(unique_chains[0], []) if len(unique_chains) > 0 else []
    chain_b = chain_indices.get(unique_chains[1], []) if len(unique_chains) > 1 else []

    return all_indices, chain_a, chain_b, chain_indices


def calculate_interface_contacts(cif_path: Path, confidences: Dict,
                                 pae_cutoff: float = 10.0,
                                 distance_cutoff: float = 8.0) -> List[Dict]:
    """
    Calculate interface contacts using both PAE and spatial (Cbeta-Cbeta) filtering.

    Returns list of contact dictionaries with details.
    """
    if not cif_path.exists():
        return []

    parser = CIFParser(str(cif_path))
    if not parser.parse_atoms():
        return []

    chain_ids = parser.get_chain_ids()
    if len(chain_ids) < 2:
        return []

    pae_matrix = np.array(confidences.get('pae', []))
    if pae_matrix.size == 0:
        return []

    _, chain_a_indices, chain_b_indices, _ = extract_chain_indices(confidences)
    if not chain_a_indices or not chain_b_indices:
        return []

    chain_a_id = chain_ids[0]
    chain_b_id = chain_ids[1]
    chain_a_residues = parser.get_chain_residues(chain_a_id)
    chain_b_residues = parser.get_chain_residues(chain_b_id)

    contacts = []

    for i, res_a in enumerate(chain_a_residues):
        for j, res_b in enumerate(chain_b_residues):
            # Check spatial distance (CB-CB)
            dist = parser.calculate_cb_distance(chain_a_id, res_a, chain_b_id, res_b)
            if dist > distance_cutoff:
                continue

            # Map residue indices to PAE matrix indices
            a_idx = chain_a_indices[i] if i < len(chain_a_indices) else -1
            b_idx = chain_b_indices[j] if j < len(chain_b_indices) else -1

            if a_idx < 0 or b_idx < 0 or a_idx >= len(pae_matrix) or b_idx >= len(pae_matrix[0]):
                continue

            pae_val = pae_matrix[a_idx, b_idx]
            if pae_val > pae_cutoff:
                continue

            contacts.append({
                'chain_a_res': res_a,
                'chain_b_res': res_b,
                'distance': round(dist, 2),
                'pae': round(float(pae_val), 2)
            })

    return contacts


def analyze_sample(pred_dir: Path, prefix: str, seed: int, sample_num: int,
                   ipsae_pae_cutoff: float = 10.0) -> Optional[Dict]:
    """
    Analyze a single seed/sample model.

    Args:
        pred_dir: The prediction directory (e.g., .../AF3/q9bw83_and_p12345/)
        prefix: Prediction name prefix
        seed: Seed number
        sample_num: Sample number
        ipsae_pae_cutoff: PAE cutoff for ipSAE calculation

    Returns dictionary with analysis results.
    """
    try:
        # Load confidences
        confidences = load_confidences(pred_dir, prefix, seed, sample_num)
        if not confidences:
            return None

        # Load summary
        summary = load_summary(pred_dir, prefix, seed, sample_num)
        if not summary:
            return None

        # Load data.json for sequences
        data_file = pred_dir / f"{prefix}_data.json"
        if not data_file.exists():
            data_file = pred_dir / "data.json"

        sequences = {}
        if data_file.exists():
            with open(data_file, 'r') as f:
                data = json.load(f)
            if 'sequences' in data:
                for seq_info in data['sequences']:
                    protein = seq_info.get('protein', {})
                    if isinstance(protein, dict):
                        chain_id = protein.get('id', 'unknown')
                        sequences[chain_id] = protein.get('sequence', '')
                    else:
                        # Flat format fallback
                        chain_id = seq_info.get('protein.id', 'unknown')
                        sequences[chain_id] = seq_info.get('protein.sequence', '')

        # Get PAE matrix and calculate ipSAE
        pae_matrix = np.array(confidences.get('pae', []))
        if pae_matrix.size == 0:
            return None

        _, chain_a_indices, chain_b_indices, _ = extract_chain_indices(confidences)

        if not chain_a_indices or not chain_b_indices:
            return None

        ipsae = calculate_ipsae(pae_matrix, chain_a_indices, chain_b_indices, ipsae_pae_cutoff)

        # Find CIF file for this sample (spatial filtering requires structure)
        sample_dir = pred_dir / f"seed-{seed}_sample-{sample_num}"
        cif_candidates = [
            sample_dir / f"{prefix}_seed-{seed}_sample-{sample_num}_model.cif",
            sample_dir / "model.cif",
            pred_dir / f"{prefix}_seed-{seed}_sample-{sample_num}_model.cif",
            pred_dir / f"{prefix}_model.cif",   # top-ranked fallback
            pred_dir / "model.cif",
        ]
        cif_path = next((p for p in cif_candidates if p.exists()), None)

        # Count interface contacts and collect interface residues for pLDDT
        contacts_pae3 = contacts_pae5 = contacts_pae8 = 0
        interface_residue_indices = set()  # token indices of interface residues
        if cif_path:
            parser = CIFParser(str(cif_path))
            if parser.parse_atoms():
                chain_ids = parser.get_chain_ids()
                if len(chain_ids) >= 2:
                    chain_a_id, chain_b_id = chain_ids[0], chain_ids[1]
                    chain_a_res = parser.get_chain_residues(chain_a_id)
                    chain_b_res = parser.get_chain_residues(chain_b_id)

                    # Pre-build Cb coordinate dict: (chain_id, seq_id) -> (x, y, z)
                    cb_coords = {}
                    for atom in parser.atoms:
                        key = (atom['chain_id'], atom['seq_id'])
                        name = atom['atom_name']
                        if name == 'CB' or (name == 'CA' and key not in cb_coords):
                            cb_coords[key] = (atom['x'], atom['y'], atom['z'])

                    for i, res_a in enumerate(chain_a_res):
                        a_idx = chain_a_indices[i] if i < len(chain_a_indices) else -1
                        if a_idx < 0:
                            continue
                        coord_a = cb_coords.get((chain_a_id, res_a))
                        if coord_a is None:
                            continue

                        for j, res_b in enumerate(chain_b_res):
                            b_idx = chain_b_indices[j] if j < len(chain_b_indices) else -1
                            if b_idx < 0 or a_idx >= len(pae_matrix) or b_idx >= pae_matrix.shape[1]:
                                continue

                            # PAE check first — cheap array lookup
                            pae_val = pae_matrix[a_idx, b_idx]
                            if pae_val > 8.0:
                                continue

                            # Distance check only for PAE-passing pairs
                            coord_b = cb_coords.get((chain_b_id, res_b))
                            if coord_b is None:
                                continue
                            dx = coord_a[0] - coord_b[0]
                            dy = coord_a[1] - coord_b[1]
                            dz = coord_a[2] - coord_b[2]
                            dist = (dx*dx + dy*dy + dz*dz) ** 0.5
                            if dist > 8.0:
                                continue

                            # Both filters passed — bin by PAE threshold
                            if pae_val <= 3:
                                contacts_pae3 += 1
                            if pae_val <= 5:
                                contacts_pae5 += 1
                            contacts_pae8 += 1

                            # Track interface residues (PAE ≤ 6Å) for pLDDT calc
                            if pae_val <= 6.0:
                                interface_residue_indices.add(a_idx)
                                interface_residue_indices.add(b_idx)

        # Calculate interface pLDDT from spatially-validated interface residues only
        atom_plddts = confidences.get('atom_plddts', [])
        interface_plddt = 0.0
        if atom_plddts and interface_residue_indices and cif_path:
            # Map interface token indices back to (chain_id, res_id)
            token_chain_ids = confidences.get('token_chain_ids', [])
            token_res_ids = confidences.get('token_res_ids', [])
            interface_residues = set()
            for idx in interface_residue_indices:
                if idx < len(token_chain_ids) and idx < len(token_res_ids):
                    interface_residues.add((token_chain_ids[idx], token_res_ids[idx]))

            # Get per-atom pLDDT for interface residues from CIF atom ordering
            interface_plddt_values = []
            for atom_idx, atom in enumerate(parser.atoms):
                if atom_idx >= len(atom_plddts):
                    break
                if (atom['chain_id'], atom['seq_id']) in interface_residues:
                    interface_plddt_values.append(float(atom_plddts[atom_idx]))

            if interface_plddt_values:
                interface_plddt = float(np.mean(interface_plddt_values))
        elif atom_plddts and not interface_residue_indices:
            # No spatially-validated interface — fall back to 0
            interface_plddt = 0.0

        # Calculate fraction disordered
        fraction_disordered = summary.get('fraction_disordered', 0)

        return {
            'seed': seed,
            'sample': sample_num,
            'iptm': float(summary.get('iptm', 0)),
            'ptm': float(summary.get('ptm', 0)),
            'ranking_score': float(summary.get('ranking_score', 0)),
            'ipsae': round(ipsae, 4),
            'interface_plddt': round(float(interface_plddt), 2),
            'contacts_pae3': contacts_pae3,
            'contacts_pae5': contacts_pae5,
            'contacts_pae8': contacts_pae8,
            'fraction_disordered': fraction_disordered,
            'pae_mean': float(np.mean(pae_matrix)) if pae_matrix.size > 0 else 0,
            'sequences': sequences
        }

    except Exception as e:
        print(f"Error analyzing sample {seed}-{sample_num}: {e}")
        return None


def analyze_prediction_all_models(pred_dir: Path, ipsae_pae_cutoff: float = 10.0,
                                   top_only: bool = False) -> List[Dict]:
    """
    Analyze models in a prediction directory.

    Supports three input formats:
      * AF3 local pipeline (seed-N_sample-M/ subdirs + ranking_scores.csv)
      * AF3 Server flat (fold_*_model_N.cif + fold_*_full_data_N.json)
      * AlphaPulldown / AF2-multimer (ranked_N.pdb + ranking_debug.json)

    The AlphaPulldown path is partial: no PAE matrix is available, so ipSAE
    and PAE-binned contacts are ``None``; iPTM is the combined iptm+ptm
    score, and ``interface_plddt`` is mean B-factor pLDDT over CB atoms in
    spatial CB-CB contact (< 8 Å).

    Args:
        top_only: If True, only analyze the top-ranked model.
    Returns list of analysis results, one per model.
    """
    from core.scanner import resolve_prediction_dir, _is_af3_server_flat, _is_alphapulldown_format
    pred_dir = resolve_prediction_dir(pred_dir)

    if _is_af3_server_flat(pred_dir):
        return _analyze_server_flat(pred_dir, ipsae_pae_cutoff, top_only)

    if _is_alphapulldown_format(pred_dir):
        return _analyze_alphapulldown(pred_dir, top_only)

    results = []
    pred_name = pred_dir.name
    prefix = pred_name

    ranking_file = pred_dir / f"{pred_name}_ranking_scores.csv"
    if not ranking_file.exists():
        ranking_file = pred_dir / "ranking_scores.csv"
    if not ranking_file.exists():
        return []

    df = pd.read_csv(ranking_file)
    top_ranked = df.iloc[0] if len(df) > 0 else None

    rows = df.iloc[:1] if top_only else df
    for _, row in rows.iterrows():
        seed = int(row['seed'])
        sample = int(row['sample'])

        result = analyze_sample(pred_dir, prefix, seed, sample, ipsae_pae_cutoff)
        if result:
            is_top = (top_ranked is not None and
                     int(top_ranked['seed']) == seed and
                     int(top_ranked['sample']) == sample)
            result['is_top_ranked'] = is_top
            result['prediction_name'] = pred_name
            results.append(result)

    return results


def _analyze_server_flat(pred_dir: Path, ipsae_pae_cutoff: float = 10.0,
                          top_only: bool = False) -> List[Dict]:
    """
    Analyze an AF3 Server flat-format prediction folder.

    Files are named fold_<name>_model_N.cif, fold_<name>_full_data_N.json,
    fold_<name>_summary_confidences_N.json where N = 0..4.
    full_data_N.json contains the PAE matrix (same structure as confidences.json).
    """
    import re as _re

    pred_name = pred_dir.name
    # Infer prefix from summary files
    summary_files = sorted(pred_dir.glob('*_summary_confidences_*.json'))
    if not summary_files:
        return []
    prefix_match = _re.match(r'^(.+)_summary_confidences_\d+\.json$', summary_files[0].name)
    if not prefix_match:
        return []
    prefix = prefix_match.group(1)

    # Collect model indices from available model files
    model_indices = sorted(
        int(m.group(1))
        for f in pred_dir.glob(f'{prefix}_model_*.cif')
        if (m := _re.search(r'_model_(\d+)\.cif$', f.name))
    )
    if not model_indices:
        return []

    # Determine top-ranked by ranking_score
    scores = {}
    for idx in model_indices:
        sf = pred_dir / f"{prefix}_summary_confidences_{idx}.json"
        if sf.exists():
            with open(sf) as f:
                s = json.load(f)
            scores[idx] = float(s.get('ranking_score', 0))
    top_idx = max(scores, key=scores.get) if scores else model_indices[0]

    indices_to_run = [top_idx] if top_only else model_indices

    results = []
    for idx in indices_to_run:
        result = _analyze_server_model(pred_dir, prefix, idx, ipsae_pae_cutoff)
        if result:
            result['is_top_ranked'] = (idx == top_idx)
            result['prediction_name'] = pred_name
            results.append(result)

    return results


def _analyze_server_model(pred_dir: Path, prefix: str, model_idx: int,
                           ipsae_pae_cutoff: float = 10.0) -> Optional[Dict]:
    """
    Analyze a single model from an AF3 Server flat-format folder.
    """
    try:
        # Load summary
        summary_file = pred_dir / f"{prefix}_summary_confidences_{model_idx}.json"
        if not summary_file.exists():
            return None
        with open(summary_file) as f:
            summary = json.load(f)

        # Load full_data (contains PAE, atom_plddts, token_chain_ids, etc.)
        full_data_file = pred_dir / f"{prefix}_full_data_{model_idx}.json"
        if not full_data_file.exists():
            return None
        with open(full_data_file) as f:
            confidences = json.load(f)

        pae_matrix = np.array(confidences.get('pae', []))
        if pae_matrix.size == 0:
            return None

        _, chain_a_indices, chain_b_indices, _ = extract_chain_indices(confidences)
        if not chain_a_indices or not chain_b_indices:
            return None

        ipsae = calculate_ipsae(pae_matrix, chain_a_indices, chain_b_indices, ipsae_pae_cutoff)

        # CIF file for spatial contact calculation
        cif_path = pred_dir / f"{prefix}_model_{model_idx}.cif"

        contacts_pae3 = contacts_pae5 = contacts_pae8 = 0
        interface_residue_indices: Set[int] = set()

        if cif_path.exists():
            parser = CIFParser(str(cif_path))
            if parser.parse_atoms():
                chain_ids = parser.get_chain_ids()
                if len(chain_ids) >= 2:
                    chain_a_id, chain_b_id = chain_ids[0], chain_ids[1]
                    chain_a_res = parser.get_chain_residues(chain_a_id)
                    chain_b_res = parser.get_chain_residues(chain_b_id)

                    cb_coords: Dict[tuple, tuple] = {}
                    for atom in parser.atoms:
                        key = (atom['chain_id'], atom['seq_id'])
                        name = atom['atom_name']
                        if name == 'CB' or (name == 'CA' and key not in cb_coords):
                            cb_coords[key] = (atom['x'], atom['y'], atom['z'])

                    for i, res_a in enumerate(chain_a_res):
                        a_idx = chain_a_indices[i] if i < len(chain_a_indices) else -1
                        if a_idx < 0:
                            continue
                        coord_a = cb_coords.get((chain_a_id, res_a))
                        if coord_a is None:
                            continue
                        for j, res_b in enumerate(chain_b_res):
                            b_idx = chain_b_indices[j] if j < len(chain_b_indices) else -1
                            if b_idx < 0 or a_idx >= len(pae_matrix) or b_idx >= pae_matrix.shape[1]:
                                continue
                            pae_val = pae_matrix[a_idx, b_idx]
                            if pae_val > 8.0:
                                continue
                            coord_b = cb_coords.get((chain_b_id, res_b))
                            if coord_b is None:
                                continue
                            dx = coord_a[0] - coord_b[0]
                            dy = coord_a[1] - coord_b[1]
                            dz = coord_a[2] - coord_b[2]
                            dist = (dx*dx + dy*dy + dz*dz) ** 0.5
                            if dist > 8.0:
                                continue
                            if pae_val <= 3:
                                contacts_pae3 += 1
                            if pae_val <= 5:
                                contacts_pae5 += 1
                            contacts_pae8 += 1
                            if pae_val <= 6.0:
                                interface_residue_indices.add(a_idx)
                                interface_residue_indices.add(b_idx)

        # Interface pLDDT
        atom_plddts = confidences.get('atom_plddts', [])
        interface_plddt = 0.0
        if atom_plddts and interface_residue_indices and cif_path.exists():
            token_chain_ids = confidences.get('token_chain_ids', [])
            token_res_ids = confidences.get('token_res_ids', [])
            interface_residues = set()
            for idx in interface_residue_indices:
                if idx < len(token_chain_ids) and idx < len(token_res_ids):
                    interface_residues.add((token_chain_ids[idx], token_res_ids[idx]))
            interface_plddt_values = []
            for atom_idx, atom in enumerate(parser.atoms):
                if atom_idx >= len(atom_plddts):
                    break
                if (atom['chain_id'], atom['seq_id']) in interface_residues:
                    interface_plddt_values.append(float(atom_plddts[atom_idx]))
            if interface_plddt_values:
                interface_plddt = float(np.mean(interface_plddt_values))

        return {
            'seed': 0,
            'sample': model_idx,
            'iptm': float(summary.get('iptm', 0)),
            'ptm': float(summary.get('ptm', 0)),
            'ranking_score': float(summary.get('ranking_score', 0)),
            'ipsae': round(ipsae, 4),
            'interface_plddt': round(float(interface_plddt), 2),
            'contacts_pae3': contacts_pae3,
            'contacts_pae5': contacts_pae5,
            'contacts_pae8': contacts_pae8,
            'fraction_disordered': float(summary.get('fraction_disordered', 0)),
            'pae_mean': float(np.mean(pae_matrix)) if pae_matrix.size > 0 else 0,
            'cif_path': str(cif_path) if cif_path.exists() else None,
            'format': 'server_flat',
        }

    except Exception as e:
        print(f"Error analyzing server model {model_idx}: {e}")
        return None


def _analyze_alphapulldown(pred_dir: Path, top_only: bool = False) -> List[Dict]:
    """Analyze an AlphaPulldown (AF2-multimer) prediction folder.

    Reads ``ranking_debug.json`` for per-model iptm+ptm scores, then
    computes spatial interface contacts (CB-CB < 8 Å) and interface pLDDT
    (mean B-factor of CB atoms at the interface) from each ``ranked_N.pdb``.
    No PAE information is available in this output format.
    """
    rd_file = pred_dir / "ranking_debug.json"
    if not rd_file.is_file():
        return []
    try:
        with open(rd_file) as f:
            rd = json.load(f)
    except Exception:
        return []

    scores = rd.get("iptm+ptm") or rd.get("iptm_ptm") or {}
    order = rd.get("order") or list(scores.keys())
    if not scores or not order:
        return []

    pred_name = pred_dir.name
    indices = [0] if top_only else list(range(len(order)))

    results = []
    for idx in indices:
        model_key = order[idx] if idx < len(order) else None
        score = float(scores.get(model_key, 0)) if model_key else 0.0
        pdb_path = pred_dir / f"ranked_{idx}.pdb"
        if not pdb_path.is_file():
            continue

        contacts_count, iface_plddt = _interface_metrics_from_pdb(pdb_path)

        results.append({
            'seed': 0,
            'sample': idx,
            'iptm': score,          # combined iptm+ptm — only confidence score available
            'ptm': score,
            'ranking_score': score,
            'ipsae': None,
            'interface_plddt': round(iface_plddt, 2) if iface_plddt is not None else None,
            'contacts_pae3': None,
            'contacts_pae5': None,
            'contacts_pae8': None,
            'contacts_spatial': contacts_count,   # CB-CB pairs < 8 Å, no PAE filter
            'fraction_disordered': 0,
            'pae_mean': None,
            'is_top_ranked': (idx == 0),
            'prediction_name': pred_name,
            'format': 'alphapulldown',
            'model_key': model_key,
        })

    return results


def _interface_metrics_from_pdb(pdb_path: Path,
                                  distance_cutoff: float = 8.0
                                  ) -> Tuple[int, Optional[float]]:
    """Spatial interface contact count + mean CB pLDDT (B-factor) for an
    AF2-multimer PDB. No PAE filtering. Returns ``(n_contacts, mean_plddt)``;
    mean_plddt is ``None`` if no interface was found.
    """
    import gemmi
    try:
        structure = gemmi.read_structure(str(pdb_path))
    except Exception:
        return 0, None
    if len(structure) == 0:
        return 0, None
    model = structure[0]

    # Collect CB coords per chain (CA fallback for glycine), with B-factor pLDDT
    cb_by_chain: Dict[str, List[Tuple[int, float, float, float, float]]] = {}
    for chain in model:
        entries = []
        for residue in chain:
            cb = residue.find_atom('CB', '\0') or residue.find_atom('CA', '\0')
            if cb is not None:
                entries.append((
                    residue.seqid.num,
                    cb.pos.x, cb.pos.y, cb.pos.z,
                    cb.b_iso,
                ))
        if entries:
            cb_by_chain[chain.name] = entries

    chain_ids = list(cb_by_chain.keys())
    if len(chain_ids) < 2:
        return 0, None
    a_id, b_id = chain_ids[0], chain_ids[1]
    a_atoms = cb_by_chain[a_id]
    b_atoms = cb_by_chain[b_id]

    dist_sq = distance_cutoff * distance_cutoff
    n_contacts = 0
    iface_a, iface_b = set(), set()
    for (ra, ax, ay, az, _) in a_atoms:
        for (rb, bx, by, bz, _) in b_atoms:
            dx, dy, dz = ax - bx, ay - by, az - bz
            if dx*dx + dy*dy + dz*dz <= dist_sq:
                n_contacts += 1
                iface_a.add(ra)
                iface_b.add(rb)

    if not iface_a and not iface_b:
        return 0, None

    plddts = [bf for (r, _, _, _, bf) in a_atoms if r in iface_a]
    plddts += [bf for (r, _, _, _, bf) in b_atoms if r in iface_b]
    mean_plddt = sum(plddts) / len(plddts) if plddts else None
    return n_contacts, mean_plddt