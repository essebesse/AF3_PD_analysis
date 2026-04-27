"""
AF3 / AF2 Directory Scanner Module

Fast directory listing + lazy per-prediction loading.
Scan = list folders. Data loaded only when needed.

Supports:
  - AF3 local pipeline (seed-N_sample-M/ subdirs + ranking_scores.csv)
  - AF3 Server flat format (fold_*_model_N.cif directly in folder)
  - AF2 standard (ranked_*.pdb + ranking_debug.json + PAE/pLDDT CSVs),
    either flat in the prediction folder or inside a seq/ subdir.
"""

import csv
import json
from pathlib import Path
from typing import Dict, List, Literal, Optional, Tuple
import re


FormatName = Literal["af3_server_flat", "af3_local", "af2", "unknown"]


def resolve_prediction_dir(d: Path) -> Path:
    """
    Resolve extra nesting from AF3 Server / Universal jobs.

    AF3 Server sometimes creates: AF3/name/name/ where the inner folder
    holds the actual data (seed dirs, ranking_scores.csv, etc.).
    This detects that pattern and returns the inner directory.
    """
    inner = d / d.name
    if inner.is_dir():
        # Check if the inner dir has AF3 data files (seed dirs or ranking_scores)
        has_seeds = any(
            c.is_dir() and c.name.startswith('seed-') for c in inner.iterdir()
        )
        has_ranking = (inner / f"{d.name}_ranking_scores.csv").exists() or (inner / "ranking_scores.csv").exists()
        has_summary = (inner / f"{d.name}_summary_confidences.json").exists() or (inner / "summary_confidences.json").exists()
        if has_seeds or has_ranking or has_summary:
            return inner
    return d


def _is_af3_server_flat(d: Path) -> bool:
    """
    Detect AF3 Server flat-file format.

    The AlphaFold Server produces a folder with files named:
      fold_<name>_model_N.cif
      fold_<name>_summary_confidences_N.json
      fold_<name>_full_data_N.json
    directly inside the folder (no seed subdirectories).
    """
    try:
        files = list(d.iterdir())
        has_model = any(
            f.is_file() and re.search(r'_model_\d+\.cif$', f.name)
            for f in files
        )
        has_summary = any(
            f.is_file() and re.search(r'_summary_confidences_\d+\.json$', f.name)
            for f in files
        )
        return has_model and has_summary
    except PermissionError:
        return False


def _is_af3_local(d: Path) -> bool:
    """Detect an AF3 local-pipeline prediction folder (the one with ranking_scores.csv
    or summary_confidences.json, plus typical seed-*_sample-* subdirs)."""
    try:
        if any(d.glob('*summary_confidences.json')):
            return True
        for child in d.iterdir():
            if child.is_dir() and child.name.startswith('seed-'):
                return True
    except (PermissionError, OSError):
        pass
    return False


def _af2_work_dir(d: Path) -> Optional[Path]:
    """If d contains an AF2 prediction, return the directory that actually holds
    the ranked_*.pdb files (either d itself or d/seq). Otherwise return None."""
    try:
        for candidate in (d, d / "seq"):
            if not candidate.is_dir():
                continue
            if not (candidate / "ranking_debug.json").exists():
                continue
            has_pdb = any(candidate.glob("ranked_*.pdb"))
            has_pae = any(candidate.glob("ranked_*_PAE.csv"))
            if has_pdb and has_pae:
                return candidate
    except (PermissionError, OSError):
        pass
    return None


def detect_format(d: Path) -> FormatName:
    """Identify what kind of prediction folder d is. First match wins:

    1. AF3 Server flat (fold_*_model_N.cif directly in d)
    2. AF3 local pipeline (summary_confidences.json / seed-* subdirs)
    3. AF2 (ranking_debug.json + ranked_*.pdb + ranked_*_PAE.csv), possibly in d/seq
    4. unknown
    """
    if _is_af3_server_flat(d):
        return "af3_server_flat"
    if _is_af3_local(d):
        return "af3_local"
    if _af2_work_dir(d) is not None:
        return "af2"
    return "unknown"


def _count_chains_af3_local(pred_dir: Path) -> int:
    """Count unique chains in an AF3 local prediction by reading
    token_chain_ids from the top-level confidences.json. Cheap, no CIF parse."""
    pred_name = pred_dir.name
    candidates = [
        pred_dir / f"{pred_name}_confidences.json",
        pred_dir / "confidences.json",
    ]
    for cf in candidates:
        if cf.exists():
            try:
                with open(cf) as f:
                    data = json.load(f)
                ids = data.get('token_chain_ids', [])
                return len(set(ids)) if ids else 0
            except (OSError, ValueError):
                continue
    return 0


def _count_chains_af3_server_flat(pred_dir: Path) -> int:
    """Count chains for an AF3 server flat prediction by reading the first
    full_data_*.json. (That file mirrors confidences.json and has token_chain_ids.)"""
    full_data = sorted(pred_dir.glob('*_full_data_*.json'))
    if not full_data:
        return 0
    try:
        with open(full_data[0]) as f:
            data = json.load(f)
        ids = data.get('token_chain_ids', [])
        return len(set(ids)) if ids else 0
    except (OSError, ValueError):
        return 0


def _count_chains_af2(work_dir: Path) -> int:
    """Count unique chain IDs in the top-ranked PDB for an AF2 prediction."""
    pdb = work_dir / "ranked_0.pdb"
    if not pdb.exists():
        return 0
    chains = set()
    try:
        with open(pdb) as f:
            for line in f:
                if line.startswith('ATOM') and len(line) >= 22:
                    chains.add(line[21])
    except OSError:
        return 0
    return len(chains)


_AF2_SPECIES_PREFIX_RE = re.compile(r'^([A-Z][a-z])_(?=[A-Z])')


def _parse_af2_bait_prey(folder_name: str) -> Tuple[str, str]:
    """Best-effort bait/prey parse for AF2 folder names like
    'Arl13B_INPP5E', 'Arl13B_INPP5E_Tulp3', or 'Hs_ADCY3_ADCY6_ANKMY2'.

    Strips a leading two-letter species prefix (Hs_, Mm_, Cr_, …) only if it
    is followed by an uppercase letter. Returns the first two remaining tokens
    as (bait, prey). For a multi-gene name with three+ tokens, prey is the
    second token; the full folder name is always preserved on the scan entry
    so the user can see it.
    """
    name = _AF2_SPECIES_PREFIX_RE.sub('', folder_name)
    tokens = name.split('_') if name else []
    if len(tokens) >= 2:
        return tokens[0], tokens[1]
    if len(tokens) == 1:
        return tokens[0], ""
    return folder_name, ""


def _dir_has_predictions(d: Path) -> bool:
    """
    Check if a directory directly contains AF3 or AF2 prediction data — either
    because d itself is a prediction folder, or because one of its children is.
    """
    # Flat server format: the files are directly inside d (no subdirs needed)
    if _is_af3_server_flat(d):
        return True
    # AF2 can have ranked_*.pdb directly in d or in d/seq/
    if _af2_work_dir(d) is not None:
        return True
    try:
        for child in d.iterdir():
            if not child.is_dir() or child.name.startswith('.'):
                continue
            resolved = resolve_prediction_dir(child)
            # AF3 Server flat format in a child dir
            if _is_af3_server_flat(resolved):
                return True
            # AF2 in a child dir (flat or via seq/)
            if _af2_work_dir(resolved) is not None:
                return True
            # Local pipeline: summary_confidences.json present
            if any(resolved.glob('*summary_confidences.json')):
                return True
            # Local pipeline: seed subdirectories
            try:
                if any(c.is_dir() and c.name.startswith('seed-') for c in resolved.iterdir()):
                    return True
            except (PermissionError, OSError):
                pass
    except PermissionError:
        pass
    return False


def _folder_format_label(d: Path) -> str:
    """Inspect a folder's children and return a short badge describing the
    mix of prediction formats found inside. Used for the recursive-scan UI."""
    formats_seen = set()
    # The project folder itself may be a single flat server prediction
    if _is_af3_server_flat(d):
        formats_seen.add('AF3 Server')
    if _af2_work_dir(d) is not None:
        formats_seen.add('AF2')
    try:
        for child in d.iterdir():
            if not child.is_dir() or child.name.startswith('.'):
                continue
            resolved = resolve_prediction_dir(child)
            if _is_af3_server_flat(resolved):
                formats_seen.add('AF3 Server')
            elif _is_af3_local(resolved):
                formats_seen.add('AF3')
            elif _af2_work_dir(resolved) is not None:
                formats_seen.add('AF2')
    except (PermissionError, OSError):
        pass
    if not formats_seen:
        return 'unknown'
    return ', '.join(sorted(formats_seen))


def find_predictions_recursive(root: str, max_depth: int = 6) -> List[Dict]:
    """
    Recursively find all directories under `root` that contain AlphaFold
    prediction data (AF3 local, AF3 Server flat, or AF2).

    Returns a list of dicts with keys:
      - project_path: the parent folder a user would select as the "project"
      - af3_folder:   the folder with the actual predictions (kept as
                      'af3_folder' for backwards compatibility with callers,
                      even when it contains AF2 data)
      - label:        a human-readable relative-path label
      - format:       'AF3', 'AF3 Server', 'AF2', or a mix-separated list

    Prunes descent once a match is found (predictions dirs are leaf-level).
    """
    results = []
    root_path = Path(root)

    # Check if root itself is a flat server prediction (user browsed into it)
    if _is_af3_server_flat(root_path):
        return [{
            'project_path': str(root_path),
            'af3_folder': str(root_path),
            'label': root_path.name,
            'format': 'AF3 Server',
        }]
    # Check if root itself is an AF2 prediction folder
    if _af2_work_dir(root_path) is not None:
        return [{
            'project_path': str(root_path),
            'af3_folder': str(root_path),
            'label': root_path.name,
            'format': 'AF2',
        }]

    def _walk(current: Path, depth: int):
        if depth > max_depth:
            return
        try:
            subdirs = [d for d in current.iterdir()
                       if d.is_dir() and not d.name.startswith('.')]
        except PermissionError:
            return

        for d in sorted(subdirs):
            # d is an AF3 Server flat prediction folder itself
            if _is_af3_server_flat(d):
                results.append({
                    'project_path': str(d),
                    'af3_folder': str(d),
                    'label': str(d.relative_to(root_path)),
                    'format': 'AF3 Server',
                })
                continue  # don't descend further

            # d is an AF2 prediction folder (flat or via seq/)
            if _af2_work_dir(d) is not None:
                results.append({
                    'project_path': str(d),
                    'af3_folder': str(d),
                    'label': str(d.relative_to(root_path)),
                    'format': 'AF2',
                })
                continue

            # d/AF3/ is a predictions folder
            af3_sub = d / "AF3"
            if af3_sub.is_dir() and _dir_has_predictions(af3_sub):
                results.append({
                    'project_path': str(d),
                    'af3_folder': str(af3_sub),
                    'label': str(d.relative_to(root_path)),
                    'format': _folder_format_label(af3_sub),
                })
                continue

            # d itself is a predictions folder
            if _dir_has_predictions(d):
                results.append({
                    'project_path': str(d),
                    'af3_folder': str(d),
                    'label': str(d.relative_to(root_path)),
                    'format': _folder_format_label(d),
                })
                continue

            # Not a match — keep descending
            _walk(d, depth + 1)

    _walk(root_path, 0)
    return results


# Backwards-compatible alias (callers in app.py still import this name).
find_af3_projects_recursive = find_predictions_recursive


def find_predictions_folder(project_path: str) -> Optional[str]:
    """
    Auto-detect where AF3 predictions live.
    If <project_path>/AF3/ exists, use that. Otherwise use <project_path>/ itself.
    """
    project = Path(project_path)
    af3_sub = project / "AF3"
    if af3_sub.is_dir():
        return str(af3_sub)
    if project.is_dir():
        return str(project)
    return None


def _load_server_flat_prediction(pred_dir: Path) -> Optional[Dict]:
    """
    Load summary data for an AF3 Server flat-format prediction folder.

    Reads all fold_<name>_summary_confidences_N.json files (N=0-4) and returns
    the entry with the highest ranking_score as the representative summary, plus
    a seed_samples list ordered by ranking_score descending.
    """
    try:
        # Discover all summary files: fold_<name>_summary_confidences_N.json
        summary_files = sorted(pred_dir.glob('*_summary_confidences_*.json'))
        if not summary_files:
            return None

        # Infer prefix from first file: strip _summary_confidences_N.json
        first = summary_files[0].name
        prefix_match = re.match(r'^(.+)_summary_confidences_\d+\.json$', first)
        if not prefix_match:
            return None
        prefix = prefix_match.group(1)

        # Parse the prediction name — AF3 Server uses fold_<bait>_<target>
        pred_name = pred_dir.name
        # Try to split on last underscore-separated word after "fold_"
        name_part = pred_name
        if name_part.startswith('fold_'):
            name_part = name_part[5:]  # strip "fold_"
        # Try _and_ split first; fall back to using the full name
        and_match = re.match(r'(.+)_and_(.+)', name_part)
        if and_match:
            bait_name = and_match.group(1)
            prey_name = and_match.group(2)
        else:
            # No _and_ — store whole name as bait
            bait_name = name_part
            prey_name = ""

        # Load all per-model summaries
        seed_samples = []
        for sf in sorted(summary_files):
            idx_match = re.search(r'_summary_confidences_(\d+)\.json$', sf.name)
            if not idx_match:
                continue
            model_idx = int(idx_match.group(1))
            with open(sf, 'r') as f:
                s = json.load(f)
            seed_samples.append({
                'seed': 0,
                'sample': model_idx,
                'ranking_score': float(s.get('ranking_score', 0)),
                'iptm': float(s.get('iptm', 0)),
                'ptm': float(s.get('ptm', 0)),
            })

        if not seed_samples:
            return None

        # Sort by ranking_score descending; top-ranked model is representative
        seed_samples.sort(key=lambda x: x['ranking_score'], reverse=True)
        top = seed_samples[0]

        # Load the top-ranked summary file for detailed fields
        top_summary_file = pred_dir / f"{prefix}_summary_confidences_{top['sample']}.json"
        with open(top_summary_file, 'r') as f:
            top_summary = json.load(f)

        n_chains = _count_chains_af3_server_flat(pred_dir)
        return {
            'name': pred_name,
            'bait': bait_name,
            'prey': prey_name,
            'iptm': float(top_summary.get('iptm', 0)),
            'ptm': float(top_summary.get('ptm', 0)),
            'ranking_score': float(top_summary.get('ranking_score', 0)),
            'fraction_disordered': float(top_summary.get('fraction_disordered', 0)),
            'has_clash': float(top_summary.get('has_clash', 0)),
            'chain_iptm': top_summary.get('chain_iptm', []),
            'chain_pair_iptm': top_summary.get('chain_pair_iptm', []),
            'chain_ptm': top_summary.get('chain_ptm', []),
            'seed_samples': seed_samples,
            'prediction_dir': str(pred_dir),
            'format': 'af3_server_flat',
            'server_prefix': prefix,
            'n_chains': n_chains,
            'pairwise_ok': n_chains == 2,
        }
    except Exception as e:
        print(f"Error loading server flat prediction {pred_dir}: {e}")
        return None


def _load_af2_prediction(pred_dir: Path) -> Optional[Dict]:
    """
    Load summary data for an AF2 prediction folder (either flat layout or
    using a seq/ subdir for the ranked_*.pdb files). Returns a summary dict
    compatible with the overview page, with format='af2'.
    """
    try:
        work_dir = _af2_work_dir(pred_dir)
        if work_dir is None:
            return None

        ranking_file = work_dir / "ranking_debug.json"
        with open(ranking_file) as f:
            ranking = json.load(f)

        iptm_ptm_scores = ranking.get("iptm+ptm", {})
        order = ranking.get("order", [])

        # Build per-model list in ranked order: ranked_N corresponds to order[N]
        seed_samples = []
        for n, model_name in enumerate(order):
            score = float(iptm_ptm_scores.get(model_name, 0.0))
            seed_samples.append({
                'seed': 0,
                'sample': n,
                'ranking_score': score,   # reuse key so overview sort still works
                'iptm_ptm': score,
                'model_name': model_name,
            })
        if not seed_samples and iptm_ptm_scores:
            # Fall back to any available scores if 'order' is missing
            for n, (model_name, score) in enumerate(iptm_ptm_scores.items()):
                seed_samples.append({
                    'seed': 0,
                    'sample': n,
                    'ranking_score': float(score),
                    'iptm_ptm': float(score),
                    'model_name': model_name,
                })

        pred_name = pred_dir.name
        bait_name, prey_name = _parse_af2_bait_prey(pred_name)
        top_score = seed_samples[0]['iptm_ptm'] if seed_samples else 0.0
        n_chains = _count_chains_af2(work_dir)

        return {
            'name': pred_name,
            'bait': bait_name,
            'prey': prey_name,
            'iptm_ptm': top_score,
            'iptm': top_score,               # unified field the overview reads
            'ranking_score': top_score,      # unified field the overview sorts on
            'seed_samples': seed_samples,
            'prediction_dir': str(pred_dir),
            'work_dir': str(work_dir),
            'format': 'af2',
            'n_chains': n_chains,
            'pairwise_ok': n_chains == 2,
        }
    except Exception as e:
        print(f"Error loading AF2 prediction {pred_dir}: {e}")
        return None


def load_prediction_data(pred_dir: Path) -> Optional[Dict]:
    """
    Load summary data for a single prediction directory.
    Reads only summary_confidences.json (~330B) and ranking_scores.csv (~150B).
    Handles AF3 Server nesting (name/name/) and flat server format automatically.
    AF2 predictions (ranked_*.pdb + ranking_debug.json) are also supported.
    """
    try:
        # Resolve AF3 Server extra nesting
        pred_dir = resolve_prediction_dir(pred_dir)
        pred_name = pred_dir.name

        # --- AF2 ---
        if _af2_work_dir(pred_dir) is not None:
            return _load_af2_prediction(pred_dir)

        # --- AF3 Server flat format ---
        if _is_af3_server_flat(pred_dir):
            return _load_server_flat_prediction(pred_dir)

        # Find summary_confidences.json
        summary_file = pred_dir / f"{pred_name}_summary_confidences.json"
        if not summary_file.exists():
            summary_file = pred_dir / "summary_confidences.json"
        if not summary_file.exists():
            return None

        match = re.match(r'(.+)_and_(.+)', pred_name)
        if match:
            bait_name = match.group(1)
            prey_name = match.group(2)
        else:
            bait_name = pred_name
            prey_name = ""

        with open(summary_file, 'r') as f:
            summary = json.load(f)

        ranking_file = pred_dir / f"{pred_name}_ranking_scores.csv"
        if not ranking_file.exists():
            ranking_file = pred_dir / "ranking_scores.csv"
        seed_samples = []
        if ranking_file.exists():
            with open(ranking_file, 'r') as rf:
                reader = csv.DictReader(rf)
                for row in reader:
                    seed_samples.append({
                        'seed': int(row['seed']),
                        'sample': int(row['sample']),
                        'ranking_score': float(row['ranking_score'])
                    })

        n_chains = len({c for c in summary.get('token_chain_ids', [])}) if summary.get('token_chain_ids') else _count_chains_af3_local(pred_dir)
        return {
            'name': pred_name,
            'bait': bait_name,
            'prey': prey_name,
            'iptm': float(summary.get('iptm', 0)),
            'ptm': float(summary.get('ptm', 0)),
            'ranking_score': float(summary.get('ranking_score', 0)),
            'fraction_disordered': float(summary.get('fraction_disordered', 0)),
            'has_clash': float(summary.get('has_clash', 0)),
            'chain_iptm': summary.get('chain_iptm', []),
            'chain_pair_iptm': summary.get('chain_pair_iptm', []),
            'chain_ptm': summary.get('chain_ptm', []),
            'seed_samples': seed_samples,
            'prediction_dir': str(pred_dir),
            'format': 'af3_local',
            'n_chains': n_chains,
            'pairwise_ok': n_chains == 2,
        }
    except Exception as e:
        print(f"Error loading {pred_dir}: {e}")
        return None


class AF3Scanner:
    """Fast scanner: lists directories instantly, loads data on demand."""

    def __init__(self, predictions_folder: Path):
        self.af3_folder = Path(predictions_folder)
        self.predictions = []

    def scan(self) -> List[Dict]:
        """
        Fast scan: lists subdirectories (and detects AF3 Server flat format).
        Returns list of dicts with name, bait, prey, prediction_dir.
        """
        if not self.af3_folder.exists():
            raise ValueError(f"Predictions folder not found: {self.af3_folder}")

        self.predictions = []

        # If the folder itself is an AF3 Server flat prediction (user pointed
        # directly at the unzipped server folder), treat it as one prediction.
        if _is_af3_server_flat(self.af3_folder):
            pred_name = self.af3_folder.name
            name_part = pred_name[5:] if pred_name.startswith('fold_') else pred_name
            and_match = re.match(r'(.+)_and_(.+)', name_part)
            if and_match:
                bait_name, prey_name = and_match.group(1), and_match.group(2)
            else:
                bait_name, prey_name = name_part, ""
            n_chains = _count_chains_af3_server_flat(self.af3_folder)
            self.predictions.append({
                'name': pred_name,
                'bait': bait_name,
                'prey': prey_name,
                'prediction_dir': str(self.af3_folder),
                'format': 'af3_server_flat',
                'n_chains': n_chains,
                'pairwise_ok': n_chains == 2,
            })
            return self.predictions

        # Folder itself is an AF2 prediction (flat or via seq/)
        af2_wd = _af2_work_dir(self.af3_folder)
        if af2_wd is not None:
            pred_name = self.af3_folder.name
            bait_name, prey_name = _parse_af2_bait_prey(pred_name)
            n_chains = _count_chains_af2(af2_wd)
            self.predictions.append({
                'name': pred_name,
                'bait': bait_name,
                'prey': prey_name,
                'prediction_dir': str(self.af3_folder),
                'work_dir': str(af2_wd),
                'format': 'af2',
                'n_chains': n_chains,
                'pairwise_ok': n_chains == 2,
            })
            return self.predictions

        for d in sorted(self.af3_folder.iterdir()):
            if not d.is_dir() or d.name.startswith('seed-') or d.name.startswith('.'):
                continue

            # Resolve AF3 Server extra nesting (name/name/)
            resolved = resolve_prediction_dir(d)

            # AF3 Server flat prediction folder inside the scanned directory
            if _is_af3_server_flat(resolved):
                pred_name = resolved.name
                name_part = pred_name[5:] if pred_name.startswith('fold_') else pred_name
                and_match = re.match(r'(.+)_and_(.+)', name_part)
                if and_match:
                    bait_name, prey_name = and_match.group(1), and_match.group(2)
                else:
                    bait_name, prey_name = name_part, ""
                n_chains = _count_chains_af3_server_flat(resolved)
                self.predictions.append({
                    'name': pred_name,
                    'bait': bait_name,
                    'prey': prey_name,
                    'prediction_dir': str(resolved),
                    'format': 'af3_server_flat',
                    'n_chains': n_chains,
                    'pairwise_ok': n_chains == 2,
                })
                continue

            # AF2 prediction (flat or seq/)
            af2_child_wd = _af2_work_dir(resolved)
            if af2_child_wd is not None:
                pred_name = resolved.name
                bait_name, prey_name = _parse_af2_bait_prey(pred_name)
                n_chains = _count_chains_af2(af2_child_wd)
                self.predictions.append({
                    'name': pred_name,
                    'bait': bait_name,
                    'prey': prey_name,
                    'prediction_dir': str(resolved),
                    'work_dir': str(af2_child_wd),
                    'format': 'af2',
                    'n_chains': n_chains,
                    'pairwise_ok': n_chains == 2,
                })
                continue

            pred_name = resolved.name
            match = re.match(r'(.+)_and_(.+)', pred_name)
            if match:
                bait_name = match.group(1)
                prey_name = match.group(2)
            else:
                bait_name = pred_name
                prey_name = ""

            n_chains = _count_chains_af3_local(resolved)
            self.predictions.append({
                'name': pred_name,
                'bait': bait_name,
                'prey': prey_name,
                'prediction_dir': str(resolved),
                'format': 'af3_local',
                'n_chains': n_chains,
                'pairwise_ok': n_chains == 2,
            })

        return self.predictions

    def load_all(self) -> List[Dict]:
        """Load full summary data for all predictions (reads files)."""
        loaded = []
        for pred in self.predictions:
            info = load_prediction_data(Path(pred['prediction_dir']))
            if info:
                loaded.append(info)
        self.predictions = loaded
        return self.predictions

    def get_summary_stats(self) -> Dict:
        """Get summary statistics from loaded predictions."""
        total = len(self.predictions)
        if not total:
            return {}

        iptms = [p['iptm'] for p in self.predictions if 'iptm' in p]
        ranking_scores = [p['ranking_score'] for p in self.predictions if 'ranking_score' in p]

        return {
            'total_predictions': total,
            'avg_iptm': sum(iptms) / len(iptms) if iptms else 0,
            'max_iptm': max(iptms) if iptms else 0,
            'min_iptm': min(iptms) if iptms else 0,
            'avg_ranking_score': sum(ranking_scores) / len(ranking_scores) if ranking_scores else 0,
            'high_confidence_count': sum(1 for i in iptms if i > 0.5),
            'analyzed_count': 0,
        }

    def get_prediction_by_name(self, name: str) -> Optional[Dict]:
        """Get prediction info by name. Loads data lazily if needed."""
        for pred in self.predictions:
            if pred['name'] == name:
                if 'iptm' not in pred:
                    full = load_prediction_data(Path(pred['prediction_dir']))
                    if full:
                        pred.update(full)
                return pred
        return None


def scan_af3_directory(predictions_folder: str) -> Tuple[List[Dict], Dict]:
    """Convenience function to scan a predictions directory and load summary data."""
    scanner = AF3Scanner(Path(predictions_folder))
    scanner.scan()
    predictions = scanner.load_all()
    stats = scanner.get_summary_stats()
    return predictions, stats
