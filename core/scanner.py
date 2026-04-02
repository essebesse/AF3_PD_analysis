"""
AF3 Directory Scanner Module

Fast directory listing + lazy per-prediction loading.
Scan = list folders. Data loaded only when needed.
"""

import csv
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import re


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


def _dir_has_predictions(d: Path) -> bool:
    """Check if a directory directly contains AF3 prediction subdirectories."""
    try:
        for child in d.iterdir():
            if not child.is_dir() or child.name.startswith('.'):
                continue
            resolved = resolve_prediction_dir(child)
            # AF3 Server flat format
            if _is_af3_server_flat(resolved):
                return True
            # Check for summary_confidences.json (with or without name prefix)
            if any(resolved.glob('*summary_confidences.json')):
                return True
            # Check for seed directories
            if any(c.is_dir() and c.name.startswith('seed-') for c in resolved.iterdir()):
                return True
    except PermissionError:
        pass
    return False


def find_af3_projects_recursive(root: str, max_depth: int = 6) -> List[Dict]:
    """
    Recursively find all directories under `root` that contain AF3 prediction data.

    Returns a list of dicts with 'project_path' and 'af3_folder' (the folder
    with the actual predictions, which may be project/AF3/ or project/ itself).

    Prunes descent once a match is found (predictions dirs are leaf-level).
    """
    results = []
    root_path = Path(root)

    def _walk(current: Path, depth: int):
        if depth > max_depth:
            return
        try:
            subdirs = [d for d in current.iterdir()
                       if d.is_dir() and not d.name.startswith('.')]
        except PermissionError:
            return

        for d in sorted(subdirs):
            # Check if d is an AF3 Server flat prediction folder itself
            if _is_af3_server_flat(d):
                results.append({
                    'project_path': str(d),
                    'af3_folder': str(d),
                    'label': str(d.relative_to(root_path)),
                    'format': 'server_flat',
                })
                continue  # don't descend further

            # Check if d/AF3/ is a predictions folder
            af3_sub = d / "AF3"
            if af3_sub.is_dir() and _dir_has_predictions(af3_sub):
                results.append({
                    'project_path': str(d),
                    'af3_folder': str(af3_sub),
                    'label': str(d.relative_to(root_path)),
                })
                continue  # don't descend further into this project

            # Check if d itself is a predictions folder
            if _dir_has_predictions(d):
                results.append({
                    'project_path': str(d),
                    'af3_folder': str(d),
                    'label': str(d.relative_to(root_path)),
                })
                continue  # don't descend further

            # Not a match — keep descending
            _walk(d, depth + 1)

    _walk(root_path, 0)
    return results


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
            'format': 'server_flat',
            'server_prefix': prefix,
        }
    except Exception as e:
        print(f"Error loading server flat prediction {pred_dir}: {e}")
        return None


def load_prediction_data(pred_dir: Path) -> Optional[Dict]:
    """
    Load summary data for a single prediction directory.
    Reads only summary_confidences.json (~330B) and ranking_scores.csv (~150B).
    Handles AF3 Server nesting (name/name/) and flat server format automatically.
    """
    try:
        # Resolve AF3 Server extra nesting
        pred_dir = resolve_prediction_dir(pred_dir)
        pred_name = pred_dir.name

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
            self.predictions.append({
                'name': pred_name,
                'bait': bait_name,
                'prey': prey_name,
                'prediction_dir': str(self.af3_folder),
                'format': 'server_flat',
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
                self.predictions.append({
                    'name': pred_name,
                    'bait': bait_name,
                    'prey': prey_name,
                    'prediction_dir': str(resolved),
                    'format': 'server_flat',
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

            self.predictions.append({
                'name': pred_name,
                'bait': bait_name,
                'prey': prey_name,
                'prediction_dir': str(resolved),
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
