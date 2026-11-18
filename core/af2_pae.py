"""
AF2 / AlphaPulldown input helpers.

AF2-multimer output comes in two flavours that this app reads:

* **AlphaPulldown** — ``ranked_N.pdb`` + ``ranking_debug.json``. The PAE matrix
  is only stored inside the multi-GB ``result_*.pkl`` files, which the app does
  not read. Analysis is limited to spatial contacts and pLDDT.
* **AF2 with PAE CSVs** — the same files plus ``ranked_N_PAE.csv`` (square L×L
  matrix with an index row and column). Full analysis: ipSAE, PAE-binned
  contacts, PAE plots.

Either layout may sit directly in the prediction folder or in a ``seq/``
subfolder of it.

``load_af2_confidences()`` builds a dict shaped like AF3's ``confidences.json``
(``pae``, ``token_chain_ids``, ``token_res_ids``, ``atom_plddts``) so the AF3
code paths (ipSAE, contacts, PAE plots, 3D PAE colouring) run unchanged.
``atom_plddts`` follows the heavy-atom order of ``CIFParser.atoms``.
"""

from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np

try:
    import gemmi
except ImportError as e:  # pragma: no cover
    raise ImportError("AF2 support requires gemmi (`pip install gemmi`).") from e


def af2_work_dir(d: Path) -> Optional[Path]:
    """Return the folder holding ``ranking_debug.json`` + ``ranked_0.pdb``
    (``d`` itself or ``d/seq``), or ``None`` if ``d`` is not an AF2 prediction."""
    try:
        for candidate in (d, d / "seq"):
            if (candidate / "ranking_debug.json").is_file() and \
                    (candidate / "ranked_0.pdb").is_file():
                return candidate
    except (PermissionError, OSError):
        pass
    return None


def af2_pae_file(work_dir: Path, model_idx: int) -> Optional[Path]:
    """Path to ``ranked_<idx>_PAE.csv`` if it exists."""
    p = work_dir / f"ranked_{model_idx}_PAE.csv"
    return p if p.is_file() else None


def af2_has_pae(work_dir: Path) -> bool:
    return af2_pae_file(work_dir, 0) is not None


# Single definition of "hydrogen" for the whole app. CIFParser.atoms and the
# AF2 atom_plddts list are aligned by position, so both must skip the same atoms.
HYDROGEN_ELEMENTS = ("H", "D")


def is_hydrogen(atom: "gemmi.Atom") -> bool:
    return atom.element.name in HYDROGEN_ELEMENTS


def load_af2_pae(pae_csv: Path) -> np.ndarray:
    """Read a square PAE matrix from ``ranked_N_PAE.csv`` (first row and first
    column are 0-based residue indices)."""
    m = np.loadtxt(pae_csv, delimiter=",", skiprows=1)
    if m.ndim != 2 or m.shape[1] < 2:
        raise ValueError(f"unexpected PAE CSV shape {m.shape}")
    m = m[:, 1:]
    if m.shape[0] != m.shape[1]:
        raise ValueError(f"PAE matrix is {m.shape[0]}x{m.shape[1]}, not square")
    return m


def load_af2_confidences(work_dir: Path, model_idx: int) -> Optional[Dict]:
    """AF3-shaped confidences dict for AF2 model ``ranked_<idx>``, or ``None``
    when the PAE CSV or PDB is missing / unreadable. See
    ``load_af2_confidences_checked`` for the reason."""
    return load_af2_confidences_checked(work_dir, model_idx)[0]


def load_af2_confidences_checked(work_dir: Path, model_idx: int) -> Tuple[Optional[Dict], Optional[str]]:
    """``(confidences, problem)``. ``problem`` is a sentence for the user when
    a ``ranked_<idx>_PAE.csv`` exists but cannot be used (the model then gets
    the analysis without PAE); ``None`` when there is no PAE file or it loaded.

    One token per residue in PDB order (AF2 PAE is per residue). pLDDT is
    taken from the PDB B-factors, one value per heavy atom.
    """
    pae_csv = af2_pae_file(work_dir, model_idx)
    pdb = work_dir / f"ranked_{model_idx}.pdb"
    if pae_csv is None or not pdb.is_file():
        return None, None
    try:
        pae = load_af2_pae(pae_csv)
        structure = gemmi.read_structure(str(pdb))
    except Exception as e:
        return None, f"{pae_csv.name} could not be used ({e}); analysed without PAE."
    if len(structure) == 0:
        return None, f"{pdb.name} contains no model; analysed without PAE."

    token_chain_ids, token_res_ids, atom_plddts = [], [], []
    for chain in structure[0]:
        for residue in chain:
            if residue.is_water():
                continue
            token_chain_ids.append(chain.name)
            token_res_ids.append(int(residue.seqid.num))
            for atom in residue:
                if not is_hydrogen(atom):
                    atom_plddts.append(float(atom.b_iso))

    if len(token_chain_ids) != pae.shape[0]:
        return None, (f"{pae_csv.name} is {pae.shape[0]}x{pae.shape[1]} but {pdb.name} has "
                      f"{len(token_chain_ids)} residues; analysed without PAE.")

    return {
        "pae": pae.tolist(),   # list of lists, same as AF3 confidences.json
        "token_chain_ids": token_chain_ids,
        "token_res_ids": token_res_ids,
        "atom_plddts": atom_plddts,
    }, None
