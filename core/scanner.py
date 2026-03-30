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


def _dir_has_predictions(d: Path) -> bool:
    """Check if a directory directly contains AF3 prediction subdirectories."""
    try:
        for child in d.iterdir():
            if not child.is_dir() or child.name.startswith('.'):
                continue
            resolved = resolve_prediction_dir(child)
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


def load_prediction_data(pred_dir: Path) -> Optional[Dict]:
    """
    Load summary data for a single prediction directory.
    Reads only summary_confidences.json (~330B) and ranking_scores.csv (~150B).
    Handles AF3 Server nesting (name/name/) automatically.
    """
    try:
        # Resolve AF3 Server extra nesting
        pred_dir = resolve_prediction_dir(pred_dir)
        pred_name = pred_dir.name

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
        Fast scan: just lists subdirectories. No file I/O.
        Returns list of dicts with name, bait, prey, prediction_dir.
        """
        if not self.af3_folder.exists():
            raise ValueError(f"Predictions folder not found: {self.af3_folder}")

        self.predictions = []
        for d in sorted(self.af3_folder.iterdir()):
            if not d.is_dir() or d.name.startswith('seed-') or d.name.startswith('.'):
                continue

            # Resolve AF3 Server extra nesting (name/name/)
            resolved = resolve_prediction_dir(d)
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
