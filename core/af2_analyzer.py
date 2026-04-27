"""
AF2 analyzer — produces per-model metrics for AlphaFold 2 predictions in a
schema aligned with the AF3 analyzer so the downstream pages can consume
both. Pairwise-only: if the model has three or more chains, only chains A
and B are scored (pairwise_ok=False is flagged so the UI can warn the user).

Inputs per model N (N = 0..4):
  ranked_N.pdb          — structure (parsed with gemmi)
  ranked_N_PAE.csv      — PAE matrix (square L×L; first row and col are indices)
  ranked_N_pLDDT.csv    — per-residue pLDDT (column 'pLDDT')
  ranking_debug.json    — iptm+ptm scores per model, plus 'order'

ipSAE, interface contacts, and interface pLDDT are computed using the same
shared helpers as the AF3 analyzer in core/analyzer.py.
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
        "The AF2 analyzer requires gemmi. Install it with `pip install gemmi`."
    ) from e

from core.analyzer import calculate_ipsae, count_interface_contacts


def _load_pae_matrix(pae_csv: Path) -> np.ndarray:
    """Read a square PAE matrix from ranked_N_PAE.csv.

    The file has an index column (row labels 0..L-1) and a header row of the
    same shape. pandas.read_csv(index_col=0) drops both, leaving the L×L float
    matrix.
    """
    return pd.read_csv(pae_csv, index_col=0).values


def _load_plddt(plddt_csv: Path) -> np.ndarray:
    """Read per-residue pLDDT values from ranked_N_pLDDT.csv."""
    df = pd.read_csv(plddt_csv)
    if "pLDDT" in df.columns:
        return df["pLDDT"].values
    # Fallback: assume the last column is the pLDDT
    return df.iloc[:, -1].values


def _representative_coord(residue: "gemmi.Residue") -> Optional[Tuple[float, float, float]]:
    """Return CB coordinate (CA for glycine / missing CB) as (x, y, z)."""
    cb = None
    ca = None
    for atom in residue:
        if atom.name == "CB" and cb is None:
            cb = atom
        elif atom.name == "CA" and ca is None:
            ca = atom
    chosen = cb or ca
    if chosen is None:
        return None
    p = chosen.pos
    return (float(p.x), float(p.y), float(p.z))


def _chain_coords_and_indices(
    structure: "gemmi.Structure",
    pae_size: int,
) -> Tuple[List[str], Dict[str, List[Optional[Tuple[float, float, float]]]], Dict[str, List[int]]]:
    """
    Walk the structure's first model in chain order. For each chain, collect
    per-residue CB/CA coordinates (in PDB order) AND the corresponding
    position in the concatenated PAE matrix (which mirrors the residue order
    in the PDB file, chain by chain).

    Returns (chain_names, chain_coords, chain_pae_indices).
    """
    chain_names: List[str] = []
    chain_coords: Dict[str, List[Optional[Tuple[float, float, float]]]] = {}
    chain_pae_indices: Dict[str, List[int]] = {}

    flat_idx = 0
    model = structure[0]
    for chain in model:
        name = chain.name
        chain_names.append(name)
        chain_coords[name] = []
        chain_pae_indices[name] = []
        for residue in chain:
            # Skip water / ligand so the residue count matches the PAE matrix
            if residue.het_flag == "H" and residue.name in {"HOH", "WAT"}:
                continue
            if flat_idx >= pae_size:
                # PAE matrix shorter than structure — stop tracking
                break
            coord = _representative_coord(residue)
            chain_coords[name].append(coord)
            chain_pae_indices[name].append(flat_idx)
            flat_idx += 1

    return chain_names, chain_coords, chain_pae_indices


def _analyze_single_model(
    work_dir: Path,
    model_idx: int,
    iptm_ptm_by_model: Dict[str, float],
    ordered_model_names: List[str],
    ipsae_pae_cutoff: float,
) -> Optional[Dict]:
    """Analyze one ranked_N model. Returns the per-model result dict or None
    if required files are missing or malformed."""
    pdb_file = work_dir / f"ranked_{model_idx}.pdb"
    pae_file = work_dir / f"ranked_{model_idx}_PAE.csv"
    plddt_file = work_dir / f"ranked_{model_idx}_pLDDT.csv"

    if not (pdb_file.exists() and pae_file.exists()):
        return None

    try:
        pae_matrix = _load_pae_matrix(pae_file)
    except Exception as exc:
        print(f"AF2: failed to read {pae_file.name}: {exc}")
        return None
    if pae_matrix.size == 0:
        return None

    try:
        structure = gemmi.read_structure(str(pdb_file))
    except Exception as exc:
        print(f"AF2: failed to read {pdb_file.name}: {exc}")
        return None

    chain_names, chain_coords, chain_pae_indices = _chain_coords_and_indices(
        structure, pae_size=pae_matrix.shape[0]
    )

    if len(chain_names) < 2:
        # Need at least two chains for an interface
        return None

    # Pairwise analysis uses chains A and B = first two PDB chains
    a_name, b_name = chain_names[0], chain_names[1]
    a_coords = chain_coords[a_name]
    b_coords = chain_coords[b_name]
    a_idxs = chain_pae_indices[a_name]
    b_idxs = chain_pae_indices[b_name]

    ipsae = calculate_ipsae(pae_matrix, a_idxs, b_idxs, ipsae_pae_cutoff)

    counts, interface_residue_indices = count_interface_contacts(
        pae_matrix, a_coords, b_coords, a_idxs, b_idxs,
    )

    # Interface pLDDT — mean over residues that were classified as interface
    interface_plddt = 0.0
    if plddt_file.exists() and interface_residue_indices:
        try:
            plddt = _load_plddt(plddt_file)
            picked = [float(plddt[i]) for i in interface_residue_indices if i < len(plddt)]
            if picked:
                interface_plddt = float(np.mean(picked))
        except Exception as exc:
            print(f"AF2: failed to read {plddt_file.name}: {exc}")

    # iPTM+pTM for this model: ranked_N corresponds to ordered_model_names[N]
    iptm_ptm: Optional[float] = None
    if model_idx < len(ordered_model_names):
        iptm_ptm = float(iptm_ptm_by_model.get(ordered_model_names[model_idx], 0.0))

    n_chains = len(chain_names)

    return {
        "format": "af2",
        "seed": 0,
        "sample": model_idx,
        "is_top_ranked": model_idx == 0,
        "iptm_ptm": iptm_ptm,
        "iptm": None,
        "ptm": None,
        "ranking_score": None,
        "ipsae": round(float(ipsae), 4),
        "interface_plddt": round(float(interface_plddt), 2),
        "contacts_pae3": counts["pae3"],
        "contacts_pae5": counts["pae5"],
        "contacts_pae8": counts["pae8"],
        "fraction_disordered": None,
        "pae_mean": float(np.mean(pae_matrix)),
        "n_chains": n_chains,
        "pairwise_ok": n_chains == 2,
    }


def analyze_af2_prediction_all_models(
    pred_dir: Path,
    ipsae_pae_cutoff: float = 10.0,
    top_only: bool = False,
) -> List[Dict]:
    """
    Analyze all ranked_N models (0..4) for an AF2 prediction.

    pred_dir may be the prediction folder itself or its seq/ subdir — this
    function resolves the actual work directory by looking for
    ranking_debug.json. If top_only, only ranked_0 is analyzed.
    """
    from core.scanner import _af2_work_dir  # local import to avoid cycle

    work_dir = _af2_work_dir(pred_dir)
    if work_dir is None:
        return []

    ranking_file = work_dir / "ranking_debug.json"
    try:
        with open(ranking_file) as f:
            ranking = json.load(f)
    except Exception as exc:
        print(f"AF2: failed to read {ranking_file}: {exc}")
        return []

    iptm_ptm_by_model = ranking.get("iptm+ptm", {}) or {}
    ordered = ranking.get("order", []) or []

    # Determine which ranked_* indices actually exist
    indices = sorted(
        int(p.stem.split("_")[1])
        for p in work_dir.glob("ranked_*.pdb")
        if p.stem.split("_")[1].isdigit()
    )
    if not indices:
        return []
    if top_only:
        indices = [0] if 0 in indices else indices[:1]

    pred_name = pred_dir.name

    results: List[Dict] = []
    for idx in indices:
        r = _analyze_single_model(
            work_dir=work_dir,
            model_idx=idx,
            iptm_ptm_by_model=iptm_ptm_by_model,
            ordered_model_names=ordered,
            ipsae_pae_cutoff=ipsae_pae_cutoff,
        )
        if r is None:
            continue
        r["prediction_name"] = pred_name
        results.append(r)

    return results
