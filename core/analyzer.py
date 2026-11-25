"""
Core Analysis Module for PulldownRanker

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

from core.af2_pae import is_hydrogen


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
        self.structure = None   # gemmi.Structure, kept for entity IDs

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
        self.structure = structure

        # Use the first model; AF3 always emits a single model per CIF.
        model = structure[0]
        for chain in model:
            cid = chain.name
            for residue in chain:
                seq_id = residue.seqid.num
                res_name = residue.name
                for atom in residue:
                    # AF3 CIFs have no hydrogens; AF2 relaxed PDBs do. Skip them
                    # so atom lists (and atom_plddts alignment) are heavy-atom only.
                    if is_hydrogen(atom):
                        continue
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

    # AF2 with ranked_N_PAE.csv: synthesize an AF3-shaped dict (sample = N)
    if sample is not None:
        from core.af2_pae import af2_work_dir, load_af2_confidences
        work_dir = af2_work_dir(pred_dir)
        if work_dir is not None:
            return load_af2_confidences(work_dir, sample)
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

    a_id, b_id = interface_chain_pair(confidences)
    chain_a = chain_indices.get(a_id, []) if a_id else []
    chain_b = chain_indices.get(b_id, []) if b_id else []

    return all_indices, chain_a, chain_b, chain_indices


def polymer_chain_ids(confidences: Dict) -> List[str]:
    """Chain IDs (in token order) of polymer chains (protein, DNA, RNA).

    AF3 uses one token per residue for polymers and one token per atom for
    ligands and ions. With atom_chain_ids (AF3 confidences), a chain counts
    as a polymer when it has more atoms than tokens; ligand and ion chains
    have exactly one atom per token, also when several ions share one chain
    (3 Ca2+ entered as one ligand: one chain, residue numbers 1-3).
    Otherwise (the AF2 path, proteins only) a chain counts as a polymer when
    its tokens span more than one residue number. Without either field every
    chain is treated as a polymer.
    """
    chains = confidences.get('token_chain_ids', [])
    res_ids = confidences.get('token_res_ids', [])
    atom_chains = confidences.get('atom_chain_ids', [])
    order: List[str] = []
    residues: Dict[str, set] = {}
    n_tokens: Dict[str, int] = {}
    for i, c in enumerate(chains):
        if c not in residues:
            order.append(c)
            residues[c] = set()
            n_tokens[c] = 0
        n_tokens[c] += 1
        if i < len(res_ids):
            residues[c].add(res_ids[i])
    if atom_chains:
        n_atoms: Dict[str, int] = {}
        for c in atom_chains:
            n_atoms[c] = n_atoms.get(c, 0) + 1
        return [c for c in order if n_atoms.get(c, 0) > n_tokens[c]]
    if not res_ids:
        return order
    return [c for c in order if len(residues[c]) > 1]


def interface_chain_pair(confidences: Dict) -> Tuple[Optional[str], Optional[str]]:
    """The two chains whose interface is scored: the first two polymer chains.

    Ligand and ion chains are skipped, so a ligand entered before the
    proteins is never scored as the bait. Falls back to the first two chains
    of any kind when there are fewer than two polymer chains.
    """
    pair = polymer_chain_ids(confidences)
    if len(pair) < 2:
        pair = list(dict.fromkeys(confidences.get('token_chain_ids', [])))
    a = pair[0] if len(pair) > 0 else None
    b = pair[1] if len(pair) > 1 else None
    return a, b


def _structure_pair(parser_chain_ids: List[str], confidences: Dict) -> Tuple[Optional[str], Optional[str]]:
    """Chain pair to analyse in a parsed structure: the scored pair from the
    confidences when both chains are present, else the first two chains."""
    a, b = interface_chain_pair(confidences)
    if a in parser_chain_ids and b in parser_chain_ids:
        return a, b
    if len(parser_chain_ids) >= 2:
        return parser_chain_ids[0], parser_chain_ids[1]
    return None, None


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

    chain_a_id, chain_b_id = _structure_pair(chain_ids, confidences)
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


def chain_token_indices(confidences: Dict) -> Dict[str, List[int]]:
    """Token indices of every chain, chains in order of first appearance."""
    out: Dict[str, List[int]] = {}
    for i, c in enumerate(confidences.get('token_chain_ids', [])):
        out.setdefault(c, []).append(i)
    return out


def chain_entities(structure, chain_ids: List[str]) -> Dict[str, str]:
    """Chain ID -> molecule ID; copies of one molecule share the ID.

    Copies are chains with identical residue sequences. mmCIF entity IDs
    (``_struct_asym.entity_id``) alone are not enough: when the same protein
    is entered twice as separate AF3 inputs, its copies get two entity IDs.
    The ID shown is the mmCIF entity ID of the first chain with that sequence;
    without entity IDs (AF2 PDB files) sequences are numbered 1, 2, ...
    """
    if structure is None or len(structure) == 0:
        return {}
    sub_to_ent = {}
    for ent in structure.entities:
        for sub in ent.subchains:
            sub_to_ent[sub] = ent.name
    by_seq: Dict[tuple, str] = {}
    out: Dict[str, str] = {}
    for chain in structure[0]:
        c = chain.name
        if c not in chain_ids or c in out:
            continue
        seq = tuple(r.name for r in chain)
        if seq not in by_seq:
            ent = next((sub_to_ent[r.subchain] for r in chain if r.subchain in sub_to_ent), None)
            if ent is None or ent in by_seq.values():
                ent = str(len(by_seq) + 1)
                while ent in by_seq.values():
                    ent += "'"
            by_seq[seq] = ent
        out[c] = by_seq[seq]
    return out


def parse_bait_override(value) -> List[str]:
    """'A,B' / 'A B' / ['A', 'B'] -> ['A', 'B']; empty -> []."""
    if not value:
        return []
    if isinstance(value, str):
        value = re.split(r'[\s,;]+', value)
    return [str(v).strip() for v in value if str(v).strip()]


def assign_bait_chains(polymers: List[str], entities: Dict[str, str],
                       folder_name: str = '', override=None) -> List[str]:
    """Bait chains among the polymer chains (the rest is prey).

    1. ``override`` (e.g. 'A,B'): those chains, if any of them is a polymer chain.
    2. Folder name with ``_with_`` (AlphaPulldown ``bait1_bait2_with_prey``):
       the last polymer chain is the prey, together with its copies (same
       entity); all other polymer chains are bait.
    3. Otherwise the first polymer chain and all its copies.
    """
    if not polymers:
        return []
    chosen = [c for c in parse_bait_override(override) if c in polymers]
    if chosen:
        return chosen

    def ent(c):
        return entities.get(c, f"chain:{c}")

    if '_with_' in (folder_name or '').lower() and len(polymers) > 1:
        prey_ent = ent(polymers[-1])
        bait = [c for c in polymers if ent(c) != prey_ent]
        return bait or polymers[:-1]
    first = ent(polymers[0])
    return [c for c in polymers if ent(c) == first]


def pair_role(a: str, b: str, bait: List[str]) -> str:
    n = (a in bait) + (b in bait)
    return {2: 'bait-bait', 1: 'bait-prey', 0: 'prey-prey'}[n]


def _scored_chains(confidences: Dict) -> List[str]:
    """Chains whose pairs are scored: the polymer chains, or all chains when
    there are fewer than two polymer chains (as interface_chain_pair)."""
    chains = polymer_chain_ids(confidences)
    if len(chains) < 2:
        chains = list(dict.fromkeys(confidences.get('token_chain_ids', [])))
    return chains


class _StructureIndex:
    """One parse of a model's structure file, reused for every chain pair:
    per chain the sorted residue numbers and their Cβ (Cα fallback)
    coordinates, and per atom a residue key aligned with ``atom_plddts``."""

    def __init__(self, struct_path: Optional[Path]):
        self.ok = False
        self.structure = None
        self.chain_ids: List[str] = []
        self.residues: Dict[str, List[int]] = {}
        self.coords: Dict[str, np.ndarray] = {}
        self.atom_keys: List[tuple] = []
        if struct_path is None or not Path(struct_path).exists():
            return
        parser = CIFParser(str(struct_path))
        if not parser.parse_atoms():
            return
        self.structure = parser.structure
        cb: Dict[tuple, tuple] = {}
        res_sets: Dict[str, set] = {}
        for atom in parser.atoms:
            key = (atom['chain_id'], atom['seq_id'])
            self.atom_keys.append(key)
            if key[0] not in res_sets:
                res_sets[key[0]] = set()
                self.chain_ids.append(key[0])
            res_sets[key[0]].add(key[1])
            name = atom['atom_name']
            if name == 'CB' or (name == 'CA' and key not in cb):
                cb[key] = (atom['x'], atom['y'], atom['z'])
        nan3 = (np.nan, np.nan, np.nan)
        for c in self.chain_ids:
            res = sorted(res_sets[c])
            self.residues[c] = res
            self.coords[c] = np.array([cb.get((c, r), nan3) for r in res], dtype=float).reshape(-1, 3)
        self.ok = True


def pair_contact_mask(struct_path: Path, chain_a: str, chain_b: str, n_a: int, n_b: int,
                      distance_cutoff: float = 8.0) -> Optional[np.ndarray]:
    """Boolean (n_a x n_b) array: True where residue i of chain A and residue j
    of chain B have Cβ–Cβ (Cα for Gly) within ``distance_cutoff``. Rows and
    columns follow the chains' token order (i-th residue = i-th token), as in
    the contact counts. None if the structure cannot be read."""
    sidx = _StructureIndex(struct_path)
    if not sidx.ok or chain_a not in sidx.coords or chain_b not in sidx.coords:
        return None
    xa = sidx.coords[chain_a][:n_a]
    xb = sidx.coords[chain_b][:n_b]
    mask = np.zeros((n_a, n_b), dtype=bool)
    with np.errstate(invalid='ignore'):
        d = np.sqrt(((xa[:, None, :] - xb[None, :, :]) ** 2).sum(-1))
        mask[:len(xa), :len(xb)] = d <= distance_cutoff   # NaN (no Cβ/Cα) -> False
    return mask


def _pair_interface_metrics(sidx: _StructureIndex, confidences: Dict, pae_matrix: np.ndarray,
                            chain_a: str, chain_b: str, idx_a: List[int], idx_b: List[int],
                            atom_plddts: Optional[np.ndarray] = None,
                            atom_res_ids: Optional[np.ndarray] = None,
                            res_id_of: Optional[Dict[tuple, int]] = None) -> Tuple[int, int, int, float]:
    """PAE-binned interface contacts and interface pLDDT for one chain pair.

    A contact is an inter-chain residue pair with Cβ–Cβ (Cα for Gly) ≤ 8 Å and
    PAE(a→b) ≤ 3 / 5 / 8 Å. Interface pLDDT is the mean per-atom pLDDT over
    residues in contacts with PAE ≤ 6 Å. The i-th residue (sorted residue
    number) of a chain in the structure is the chain's i-th token.

    Returns ``(pae3, pae5, pae8, interface_plddt)``; zeros if either chain is
    missing from the structure.
    """
    if not sidx.ok or chain_a not in sidx.residues or chain_b not in sidx.residues:
        return 0, 0, 0, 0.0
    n_rows, n_cols = pae_matrix.shape

    def side(chain, idx, limit):
        m = min(len(sidx.residues[chain]), len(idx))
        tok = np.asarray(idx[:m], dtype=int)
        xyz = sidx.coords[chain][:m]
        keep = (tok >= 0) & (tok < limit)
        return tok[keep], xyz[keep]

    ta, xa = side(chain_a, idx_a, n_rows)
    tb, xb = side(chain_b, idx_b, n_cols)
    if len(ta) == 0 or len(tb) == 0:
        return 0, 0, 0, 0.0
    pae = pae_matrix[np.ix_(ta, tb)]
    d2 = ((xa[:, 0, None] - xb[None, :, 0]) ** 2
          + (xa[:, 1, None] - xb[None, :, 1]) ** 2
          + (xa[:, 2, None] - xb[None, :, 2]) ** 2)
    with np.errstate(invalid='ignore'):
        contact = (pae <= 8.0) & (np.sqrt(d2) <= 8.0)   # NaN (no Cβ/Cα) -> False
    pae3 = int(np.count_nonzero(contact & (pae <= 3)))
    pae5 = int(np.count_nonzero(contact & (pae <= 5)))
    pae8 = int(np.count_nonzero(contact))

    interface_plddt = 0.0
    strong = contact & (pae <= 6.0)
    if atom_plddts is not None and len(atom_plddts) and strong.any():
        rows, cols = np.nonzero(strong)
        tokens = set(ta[np.unique(rows)].tolist()) | set(tb[np.unique(cols)].tolist())
        tc = confidences.get('token_chain_ids', [])
        tr = confidences.get('token_res_ids', [])
        ids = {res_id_of[(tc[k], tr[k])] for k in tokens
               if k < len(tc) and k < len(tr) and (tc[k], tr[k]) in res_id_of}
        if ids:
            mask = np.isin(atom_res_ids, list(ids))
            if mask.any():
                interface_plddt = float(np.mean(atom_plddts[mask]))
    return pae3, pae5, pae8, interface_plddt


def _pair_iptm(summary: Optional[Dict], all_chains: List[str], a: str, b: str) -> Optional[float]:
    cpi = (summary or {}).get('chain_pair_iptm')
    if not cpi or a not in all_chains or b not in all_chains:
        return None
    try:
        v = cpi[all_chains.index(a)][all_chains.index(b)]
        return float(v) if v is not None else None
    except (IndexError, TypeError, ValueError):
        return None


def _summary_pair(pairs: List[Dict], key: str = 'ipsae') -> Optional[Dict]:
    """Best bait–prey pair by ``key`` (first one on ties); if there is no
    bait–prey pair (e.g. a homodimer), the best pair of any role."""
    pool = [p for p in pairs if p['role'] == 'bait-prey'] or pairs
    if not pool:
        return None
    return max(pool, key=lambda p: p.get(key) or 0)   # max keeps the first on ties


SUMMARY_PAIR_FIELDS = ('ipsae', 'interface_plddt', 'contacts_pae3', 'contacts_pae5',
                       'contacts_pae8', 'pair_iptm')


def analyze_model_pairs(struct_path: Optional[Path], confidences: Dict,
                        summary: Optional[Dict], pae_matrix: np.ndarray,
                        folder_name: str = '', bait_override=None,
                        ipsae_pae_cutoff: float = 10.0) -> Optional[Dict]:
    """Interface metrics for every pair of polymer chains of one model.

    Returns the record fields: the summary pair (best bait–prey pair by
    ipSAE) in ``ipsae``, ``interface_plddt``, ``contacts_pae3/5/8``,
    ``chain_pair``, ``pair_iptm``; plus ``n_chains``, ``pairs`` (one entry per
    chain pair), ``bait_chains`` and ``chain_entities``. ``None`` when fewer
    than two chains are present.
    """
    tok = chain_token_indices(confidences)
    chains = [c for c in _scored_chains(confidences) if tok.get(c)]
    if len(chains) < 2:
        return None
    all_chains = list(tok.keys())
    sidx = _StructureIndex(struct_path)
    entities = chain_entities(sidx.structure, chains)
    bait = assign_bait_chains(chains, entities, folder_name, bait_override)

    atom_plddts = np.asarray(confidences.get('atom_plddts', []), dtype=float)
    res_id_of: Dict[tuple, int] = {}
    atom_res_ids = None
    if sidx.ok and len(atom_plddts):
        keys = sidx.atom_keys[:len(atom_plddts)]   # atom_plddts is in structure atom order
        atom_res_ids = np.array([res_id_of.setdefault(k, len(res_id_of)) for k in keys], dtype=int)
        atom_plddts = atom_plddts[:len(keys)]

    pairs = []
    for i in range(len(chains)):
        for j in range(i + 1, len(chains)):
            a, b = chains[i], chains[j]
            ab = calculate_asymmetric_ipsae(pae_matrix, tok[a], tok[b], ipsae_pae_cutoff)
            ba = calculate_asymmetric_ipsae(pae_matrix, tok[b], tok[a], ipsae_pae_cutoff)
            p3, p5, p8, ipl = _pair_interface_metrics(
                sidx, confidences, pae_matrix, a, b, tok[a], tok[b],
                atom_plddts, atom_res_ids, res_id_of)
            pairs.append({
                'chains': f"{a}:{b}",
                'role': pair_role(a, b, bait),
                'entities': f"{entities.get(a, '?')}:{entities.get(b, '?')}",
                'ipsae': round(max(ab, ba), 4),
                'ipsae_ab': round(ab, 4),
                'ipsae_ba': round(ba, 4),
                'pair_iptm': _pair_iptm(summary, all_chains, a, b),
                'contacts_pae3': p3,
                'contacts_pae5': p5,
                'contacts_pae8': p8,
                'interface_plddt': round(float(ipl), 2),
            })
    best = _summary_pair(pairs)
    out = {k: best[k] for k in SUMMARY_PAIR_FIELDS}
    out.update({
        'chain_pair': best['chains'],
        'pair_role': best['role'],
        'n_chains': len(polymer_chain_ids(confidences)),
        'pairs': pairs,
        'bait_chains': bait,
        'chain_entities': entities,
    })
    return out


def _spatial_pairs_from_pdb(pdb_path: Path, folder_name: str = '', bait_override=None,
                            distance_cutoff: float = 8.0) -> Optional[Dict]:
    """AF2 model without PAE: spatial Cβ–Cβ contacts and B-factor interface
    pLDDT for every chain pair; summary pair = bait–prey pair with most contacts."""
    try:
        structure = gemmi.read_structure(str(pdb_path))
    except Exception:
        return None
    if len(structure) == 0:
        return None
    cb_by_chain: Dict[str, List[Tuple[int, float, float, float, float]]] = {}
    for chain in structure[0]:
        entries = []
        for residue in chain:
            cb = residue.find_atom('CB', '\0') or residue.find_atom('CA', '\0')
            if cb is not None:
                entries.append((residue.seqid.num, cb.pos.x, cb.pos.y, cb.pos.z, cb.b_iso))
        if entries:
            cb_by_chain[chain.name] = entries
    chains = list(cb_by_chain)
    if len(chains) < 2:
        return None
    entities = chain_entities(structure, chains)
    bait = assign_bait_chains(chains, entities, folder_name, bait_override)
    arr = {c: np.array(v, dtype=float) for c, v in cb_by_chain.items()}
    pairs = []
    for i in range(len(chains)):
        for j in range(i + 1, len(chains)):
            a, b = chains[i], chains[j]
            A, B = arr[a], arr[b]
            d2 = ((A[:, 1, None] - B[None, :, 1]) ** 2 + (A[:, 2, None] - B[None, :, 2]) ** 2
                  + (A[:, 3, None] - B[None, :, 3]) ** 2)
            hit = d2 <= distance_cutoff * distance_cutoff
            n = int(np.count_nonzero(hit))
            ia, ib = hit.any(axis=1), hit.any(axis=0)
            vals = list(A[ia, 4]) + list(B[ib, 4])
            pairs.append({
                'chains': f"{a}:{b}",
                'role': pair_role(a, b, bait),
                'entities': f"{entities.get(a, '?')}:{entities.get(b, '?')}",
                'ipsae': None, 'ipsae_ab': None, 'ipsae_ba': None, 'pair_iptm': None,
                'contacts_pae3': None, 'contacts_pae5': None, 'contacts_pae8': None,
                'contacts_spatial': n,
                'interface_plddt': round(float(sum(vals)) / len(vals), 2) if vals else None,
            })
    best = _summary_pair(pairs, key='contacts_spatial')
    return {
        'contacts_spatial': best['contacts_spatial'],
        'interface_plddt': best['interface_plddt'],
        'chain_pair': best['chains'],
        'pair_role': best['role'],
        'pair_iptm': None,
        'n_chains': len(chains),
        'pairs': pairs,
        'bait_chains': bait,
        'chain_entities': entities,
    }


def analyze_sample(pred_dir: Path, prefix: str, seed: int, sample_num: int,
                   ipsae_pae_cutoff: float = 10.0, bait_override=None) -> Optional[Dict]:
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

        # Get PAE matrix and calculate ipSAE
        pae_matrix = np.array(confidences.get('pae', []))
        if pae_matrix.size == 0:
            return None

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

        pair_fields = analyze_model_pairs(cif_path, confidences, summary, pae_matrix,
                                          pred_dir.name, bait_override, ipsae_pae_cutoff)
        if pair_fields is None:
            return None

        return {
            'seed': seed,
            'sample': sample_num,
            'iptm': float(summary.get('iptm', 0)),
            'ptm': float(summary.get('ptm', 0)),
            'ranking_score': float(summary.get('ranking_score', 0)),
            'fraction_disordered': summary.get('fraction_disordered', 0),
            'pae_mean': float(np.mean(pae_matrix)) if pae_matrix.size > 0 else 0,
            'format': 'af3',
            **pair_fields,
        }

    except Exception as e:
        print(f"Error analyzing sample {seed}-{sample_num}: {e}")
        return None


def analyze_prediction_all_models(pred_dir: Path, ipsae_pae_cutoff: float = 10.0,
                                   top_only: bool = False, bait_override=None) -> List[Dict]:
    """
    Analyze models in a prediction directory.

    Supports three input formats:
      * AF3 local pipeline (seed-N_sample-M/ subdirs + ranking_scores.csv)
      * AF3 Server flat (fold_*_model_N.cif + fold_*_full_data_N.json)
      * AF2-multimer (ranked_N.pdb + ranking_debug.json, in the folder or
        in a seq/ subfolder). With ranked_N_PAE.csv files the analysis is
        full (format 'af2'); without them (plain AlphaPulldown output) it is
        partial (format 'alphapulldown'): no ipSAE or PAE-binned contacts,
        spatial Cβ contacts and B-factor pLDDT only. In both, iPTM is AF2's
        combined ranking confidence (0.8·iPTM + 0.2·pTM).

    Args:
        top_only: If True, only analyze the top-ranked model.
        bait_override: bait chain IDs ('A,B' or a list) for complexes with
            more than two chains; default rule in assign_bait_chains().
    Returns list of analysis results, one per model.
    """
    from core.scanner import resolve_prediction_dir
    pred_dir = resolve_prediction_dir(Path(pred_dir))
    override = parse_bait_override(bait_override)
    results = _analyze_prediction(pred_dir, ipsae_pae_cutoff, top_only, override)
    for r in results:
        r['prediction_dir'] = str(pred_dir)   # where to write af3_app_analysis.json
        r['bait_override'] = ','.join(override) if override else None
    return results


def _analyze_prediction(pred_dir: Path, ipsae_pae_cutoff: float, top_only: bool,
                        bait_override=None) -> List[Dict]:
    from core.scanner import (_is_af3_server_flat, _is_alphapulldown_format,
                              af3_models_on_disk, prediction_file_prefix)

    if _is_af3_server_flat(pred_dir):
        return _analyze_server_flat(pred_dir, ipsae_pae_cutoff, top_only, bait_override)

    if _is_alphapulldown_format(pred_dir):
        return _analyze_alphapulldown(pred_dir, top_only, ipsae_pae_cutoff, bait_override)

    try:
        results = []
        pred_name = pred_dir.name
        prefix = prediction_file_prefix(pred_dir)   # differs from pred_name if renamed

        ranking_file = pred_dir / f"{prefix}_ranking_scores.csv"
        if not ranking_file.exists():
            ranking_file = pred_dir / "ranking_scores.csv"
        if not ranking_file.exists():
            return []

        df = pd.read_csv(ranking_file)
        if len(df) == 0:
            return []

        # Top-ranked = highest ranking_score. AF3's local ranking_scores.csv is
        # written in (seed, sample) generation order, NOT sorted by score, so we
        # must not assume row 0 is best. The rest of the app (Detailed Analysis,
        # server-flat, AlphaPulldown) also defines "top" by max ranking_score —
        # using idxmax here keeps them consistent and fixes top_only analysing
        # the wrong model.
        if 'ranking_score' in df.columns:
            top_pos = int(df['ranking_score'].astype(float).idxmax())
        else:
            top_pos = df.index[0]
        top_ranked = df.loc[top_pos]

        top_seed, top_sample = int(top_ranked['seed']), int(top_ranked['sample'])
        rows = df.loc[[top_pos]] if top_only else df
        models = [{'seed': int(r['seed']), 'sample': int(r['sample']),
                   'ranking_score': float(r['ranking_score']) if 'ranking_score' in df.columns else 0.0}
                  for _, r in rows.iterrows()]
        if not top_only:
            # Skip models whose seed-N_sample-M files were deleted; otherwise
            # load_confidences falls back to the top-level (top model) files
            # and the top model is reported under every label.
            models = af3_models_on_disk(pred_dir, models)
        for m in models:
            seed, sample = m['seed'], m['sample']
            is_top = (seed == top_seed and sample == top_sample)

            result = analyze_sample(pred_dir, prefix, seed, sample, ipsae_pae_cutoff, bait_override)
            if result:
                result['is_top_ranked'] = is_top
                result['prediction_name'] = pred_name
                results.append(result)

        return results
    except Exception as e:
        # One malformed folder must not abort an entire batch run (the local and
        # SLURM pools re-raise worker exceptions). Skip it and keep going.
        print(f"Error analyzing prediction {pred_dir}: {e}")
        return []


def _analyze_server_flat(pred_dir: Path, ipsae_pae_cutoff: float = 10.0,
                          top_only: bool = False, bait_override=None) -> List[Dict]:
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
        result = _analyze_server_model(pred_dir, prefix, idx, ipsae_pae_cutoff, bait_override)
        if result:
            result['is_top_ranked'] = (idx == top_idx)
            result['prediction_name'] = pred_name
            results.append(result)

    return results


def _analyze_server_model(pred_dir: Path, prefix: str, model_idx: int,
                           ipsae_pae_cutoff: float = 10.0, bait_override=None) -> Optional[Dict]:
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

        # CIF file for spatial contact calculation
        cif_path = pred_dir / f"{prefix}_model_{model_idx}.cif"

        pair_fields = analyze_model_pairs(cif_path, confidences, summary, pae_matrix,
                                          pred_dir.name, bait_override, ipsae_pae_cutoff)
        if pair_fields is None:
            return None

        return {
            'seed': 0,
            'sample': model_idx,
            'iptm': float(summary.get('iptm', 0)),
            'ptm': float(summary.get('ptm', 0)),
            'ranking_score': float(summary.get('ranking_score', 0)),
            'fraction_disordered': float(summary.get('fraction_disordered', 0)),
            'pae_mean': float(np.mean(pae_matrix)) if pae_matrix.size > 0 else 0,
            'cif_path': str(cif_path) if cif_path.exists() else None,
            'format': 'server_flat',
            **pair_fields,
        }

    except Exception as e:
        print(f"Error analyzing server model {model_idx}: {e}")
        return None


def _analyze_alphapulldown(pred_dir: Path, top_only: bool = False,
                           ipsae_pae_cutoff: float = 10.0, bait_override=None) -> List[Dict]:
    """Analyze an AF2-multimer prediction folder (AlphaPulldown or AF2).

    Reads ``ranking_debug.json`` for per-model iptm+ptm scores. For each
    ``ranked_N.pdb``:

    * with ``ranked_N_PAE.csv``: full analysis like AF3 (ipSAE, PAE-binned
      contacts, interface pLDDT), ``format='af2'``;
    * without it: spatial contacts (Cβ–Cβ < 8 Å) and interface pLDDT (mean
      B-factor of interface Cβ atoms) only, ``format='alphapulldown'``.
    """
    from core.af2_pae import af2_work_dir, load_af2_confidences_checked

    work_dir = af2_work_dir(pred_dir)
    if work_dir is None:
        return []
    rd_file = work_dir / "ranking_debug.json"
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
        pdb_path = work_dir / f"ranked_{idx}.pdb"
        if not pdb_path.is_file():
            continue

        confidences, pae_problem = load_af2_confidences_checked(work_dir, idx)
        if confidences is not None:
            r = _analyze_af2_with_pae(pdb_path, confidences, ipsae_pae_cutoff,
                                      pred_name, bait_override)
            if r is not None:
                r.update({
                    'seed': 0,
                    'sample': idx,
                    'iptm': score,          # combined 0.8*iptm + 0.2*ptm
                    'ptm': score,
                    'ranking_score': score,
                    'iptm_ptm': score,
                    'fraction_disordered': 0,
                    'is_top_ranked': (idx == 0),
                    'prediction_name': pred_name,
                    'format': 'af2',
                    'model_key': model_key,
                })
                results.append(r)
                continue

        spatial = _spatial_pairs_from_pdb(pdb_path, pred_name, bait_override) or {
            'contacts_spatial': 0, 'interface_plddt': None}

        results.append({
            'seed': 0,
            'sample': idx,
            'iptm': score,          # combined iptm+ptm — only confidence score available
            'ptm': score,
            'ranking_score': score,
            'ipsae': None,
            'contacts_pae3': None,
            'contacts_pae5': None,
            'contacts_pae8': None,
            **spatial,   # contacts_spatial: Cβ–Cβ pairs ≤ 8 Å, no PAE filter
            'fraction_disordered': 0,
            'pae_mean': None,
            'is_top_ranked': (idx == 0),
            'prediction_name': pred_name,
            'format': 'alphapulldown',
            'model_key': model_key,
            'warning': pae_problem,   # set when a PAE CSV exists but could not be used
        })

    return results


def _analyze_af2_with_pae(pdb_path: Path, confidences: Dict,
                          ipsae_pae_cutoff: float = 10.0, folder_name: str = '',
                          bait_override=None) -> Optional[Dict]:
    """ipSAE, PAE-binned contacts and interface pLDDT for one AF2 model with a
    PAE matrix. Same definitions as the AF3 path."""
    try:
        pae_matrix = np.asarray(confidences['pae'], dtype=float)
        pair_fields = analyze_model_pairs(pdb_path, confidences, None, pae_matrix,
                                          folder_name, bait_override, ipsae_pae_cutoff)
        if pair_fields is None:
            return None
        return {**pair_fields, 'pae_mean': float(np.mean(pae_matrix))}
    except Exception as e:
        print(f"Error analyzing AF2 model {pdb_path}: {e}")
        return None
