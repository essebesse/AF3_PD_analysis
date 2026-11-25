"""
AF3 Directory Scanner Module

Fast directory listing + lazy per-prediction loading.
Scan = list folders. Data loaded only when needed.
"""

import csv
import os
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import re


def _has_af3_output(d: Path) -> bool:
    """d holds AF3 local-pipeline output: a ranking_scores.csv or
    summary_confidences.json (with or without the job-name prefix) or
    seed-N_sample-M/ folders."""
    try:
        with os.scandir(d) as it:
            for e in it:
                n = e.name
                if n.endswith('ranking_scores.csv') or n.endswith('summary_confidences.json'):
                    if e.is_file():
                        return True
                elif n.startswith('seed-') and '_sample-' in n and e.is_dir():
                    return True
    except OSError:
        return False
    return False


def prediction_kind(d: Path) -> Optional[str]:
    """Format of the prediction output held directly in d ('server_flat',
    'af2', 'alphapulldown' or 'af3'), or None when d holds none. AF2 counts
    when its files are in d or in d/seq/."""
    if _is_af3_server_flat(d):
        return 'server_flat'
    wd = _af2_work_dir(d)
    if wd is not None:
        return _af2_format(wd)
    if _has_af3_output(d):
        return 'af3'
    return None


def _subdirs(d: Path) -> List[Path]:
    """Sub-folders worth searching: not hidden, not AF3 seed-N_sample-M/."""
    try:
        with os.scandir(d) as it:
            return sorted(Path(e.path) for e in it
                          if e.is_dir() and not e.name.startswith(('.', 'seed-', '__pycache__')))
    except OSError:
        return []


def find_prediction_dirs(root: Path, max_depth: int = 4) -> List[Tuple[Path, str]]:
    """Every prediction folder at or below root, at most max_depth levels down,
    as (folder, format). The search does not go below a prediction folder,
    so a normal pulldown costs one listing per prediction."""
    found: List[Tuple[Path, str]] = []

    def walk(d: Path, depth: int):
        kind = prediction_kind(d)
        if kind:
            found.append((d, kind))
            return
        if depth < max_depth:
            for sub in _subdirs(d):
                walk(sub, depth + 1)

    walk(Path(root), 0)
    return found


def resolve_prediction_dir(d: Path) -> Path:
    """
    The folder below d that holds the output of the prediction named d.name.

    d itself when it holds prediction output. Otherwise a folder with the same
    name up to 3 levels down, searched only through small folders (at most 5
    sub-folders), so a run folder with extra nesting is resolved:
      name/name/            (AF3 Server / Universal jobs)
      name/AF3/name/        (run folders with msa/ beside the AF3 output)
      name/output/run1/name/
    Requiring the same name keeps a folder of separate projects from being
    taken for one prediction.
    """
    d = Path(d)
    if prediction_kind(d):
        return d

    def walk(cur: Path, depth: int) -> Optional[Path]:
        subs = _subdirs(cur)
        if len(subs) > 5:
            return None
        for sub in subs:
            if sub.name == d.name and prediction_kind(sub):
                return sub
        if depth < 3:
            for sub in subs:
                hit = walk(sub, depth + 1)
                if hit is not None:
                    return hit
        return None

    return walk(d, 1) or d


def prediction_file_prefix(pred_dir: Path) -> str:
    """Job-name prefix of the AF3 files in pred_dir (<prefix>_ranking_scores.csv,
    <prefix>_summary_confidences.json, <prefix>_summary_confidences_0.json).
    Usually the folder name; differs when the folder was renamed."""
    pred_dir = Path(pred_dir)
    name = pred_dir.name
    suffixes = ('_ranking_scores.csv', '_summary_confidences.json', '_summary_confidences_0.json')
    if any((pred_dir / f"{name}{suf}").exists() for suf in suffixes):
        return name
    for suf in suffixes:
        hits = sorted(pred_dir.glob(f"*{suf}"))
        if hits:
            return hits[0].name[:-len(suf)]
    return name


def af3_models_on_disk(pred_dir: Path, seed_samples: List[Dict]) -> List[Dict]:
    """Keep only the (seed, sample) rows whose per-model files exist, plus the
    top-ranked row (whose files are the top-level ones).

    AF3 lists all models in ranking_scores.csv, but the seed-N_sample-M/
    folders are often deleted to save space. Without this filter the
    top-level (top model) files would be reused for every row. One directory
    listing per call.
    """
    if not seed_samples or 'ranking_score' not in seed_samples[0]:
        return seed_samples
    try:
        names = set(os.listdir(pred_dir))
    except OSError:
        return seed_samples
    prefix = prediction_file_prefix(pred_dir)
    top = max(seed_samples, key=lambda ss: float(ss['ranking_score']))
    kept = []
    for ss in seed_samples:
        tag = f"seed-{ss['seed']}_sample-{ss['sample']}"
        if ss is top or tag in names or f"{prefix}_{tag}_confidences.json" in names:
            kept.append(ss)
    return kept


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
    except OSError:   # missing folder, not a folder, no permission
        return False


def _is_alphapulldown_format(d: Path) -> bool:
    """
    Detect AF2-multimer output (AlphaPulldown or AF2 with PAE CSVs).

    ``ranked_<N>.pdb`` + ``ranking_debug.json`` either directly in ``d`` or in
    ``d/seq``. See ``core/af2_pae.py`` for the two flavours.
    """
    from core.af2_pae import af2_work_dir
    return af2_work_dir(d) is not None


def _af2_work_dir(d: Path) -> Optional[Path]:
    from core.af2_pae import af2_work_dir
    return af2_work_dir(d)


def _af2_format(work_dir: Path) -> str:
    """'af2' when ranked_0_PAE.csv is present (full analysis), else
    'alphapulldown' (no PAE). ``work_dir`` is the folder holding ranked_*.pdb."""
    from core.af2_pae import af2_has_pae
    return 'af2' if af2_has_pae(work_dir) else 'alphapulldown'


def _af2_entry(d: Path, work_dir: Path) -> Dict:
    """Scan entry for an AF2 prediction folder ``d`` (files in ``work_dir``,
    which is ``d`` or ``d/seq``)."""
    pred_name = d.name
    match = re.match(r'(.+)_and_(.+)', pred_name)
    if match:
        bait_name, prey_name = match.group(1), match.group(2)
    else:
        bait_name, prey_name = pred_name, ""
    return {
        'name': pred_name,
        'bait': bait_name,
        'prey': prey_name,
        'prediction_dir': str(d),
        'format': _af2_format(work_dir),
    }


def _dir_has_predictions(d: Path) -> bool:
    """
    Check if a directory directly contains prediction data: AF3 Server flat
    files in d itself, or a sub-folder that holds prediction output (also
    behind same-name nesting, see resolve_prediction_dir).
    """
    if _is_af3_server_flat(d):
        return True
    for child in _subdirs(d):
        if prediction_kind(resolve_prediction_dir(child)):
            return True
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

    # Check if root itself is a flat server prediction (user browsed into it)
    if _is_af3_server_flat(root_path):
        pred_name = root_path.name
        name_part = pred_name[5:] if pred_name.startswith('fold_') else pred_name
        and_match = re.match(r'(.+)_and_(.+)', name_part)
        if and_match:
            bait_name, prey_name = and_match.group(1), and_match.group(2)
        else:
            bait_name, prey_name = name_part, ""
        return [{
            'project_path': str(root_path),
            'af3_folder': str(root_path),
            'label': root_path.name,
            'format': 'server_flat',
        }]
    # Root itself is one AF2 prediction (ranked_*.pdb here or in seq/)
    root_wd = _af2_work_dir(root_path)
    if root_wd is not None:
        return [{
            'project_path': str(root_path),
            'af3_folder': str(root_path),
            'label': root_path.name,
            'format': _af2_format(root_wd),
        }]

    # Root itself holds prediction folders (user browsed into a pulldown):
    # offer root, not every prediction as its own project
    if _dir_has_predictions(root_path):
        return [{
            'project_path': str(root_path),
            'af3_folder': str(root_path),
            'label': root_path.name,
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
            # Check if d is an AF3 Server flat prediction folder itself
            if _is_af3_server_flat(d):
                results.append({
                    'project_path': str(d),
                    'af3_folder': str(d),
                    'label': str(d.relative_to(root_path)),
                    'format': 'server_flat',
                })
                continue  # don't descend further

            # d is a single AF2 prediction whose files sit in d/seq/. (A flat
            # AlphaPulldown pair folder is found via its parent instead, so a
            # pulldown/ folder is listed once, not once per pair.)
            d_wd = _af2_work_dir(d)
            if d_wd is not None and d_wd != d:
                results.append({
                    'project_path': str(d),
                    'af3_folder': str(d),
                    'label': str(d.relative_to(root_path)),
                    'format': _af2_format(d_wd),
                })
                continue

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


def _load_alphapulldown_prediction(pred_dir: Path, work_dir: Optional[Path] = None) -> Optional[Dict]:
    """
    Load summary data for an AlphaPulldown (AF2-multimer) prediction.

    Reads ``ranking_debug.json`` to get the per-model combined iptm+ptm
    scores. iPTM and PTM cannot be separated in this format, so we report
    the combined value as ``iptm`` (the only ranking score available) and
    leave PAE-derived fields unset.
    """
    try:
        work_dir = work_dir or _af2_work_dir(pred_dir)
        if work_dir is None:
            return None
        rd_file = work_dir / "ranking_debug.json"
        with open(rd_file, 'r') as f:
            rd = json.load(f)

        scores = rd.get("iptm+ptm") or rd.get("iptm_ptm") or {}
        order = rd.get("order") or list(scores.keys())
        if not scores or not order:
            return None

        # Top-ranked model first, then the rest in scoring order
        seed_samples = []
        for idx, model_key in enumerate(order):
            score = float(scores.get(model_key, 0))
            seed_samples.append({
                'seed': 0,
                'sample': idx,
                'ranking_score': score,
                'iptm': score,       # combined iptm+ptm, kept here for per-model display
                'ptm': score,
                'model_key': model_key,  # original AF2 model name (model_1_multimer_v3_pred_0)
            })

        top = seed_samples[0]
        pred_name = pred_dir.name
        match = re.match(r'(.+)_and_(.+)', pred_name)
        if match:
            bait_name, prey_name = match.group(1), match.group(2)
        else:
            bait_name, prey_name = pred_name, ""

        return {
            'name': pred_name,
            'bait': bait_name,
            'prey': prey_name,
            'iptm': top['ranking_score'],
            'ptm': top['ranking_score'],
            'ranking_score': top['ranking_score'],
            'fraction_disordered': 0,
            'has_clash': 0,
            'chain_iptm': [],
            'chain_pair_iptm': [],
            'chain_ptm': [],
            'seed_samples': seed_samples,
            'prediction_dir': str(pred_dir),
            'format': _af2_format(work_dir),
        }
    except Exception as e:
        print(f"Error loading AlphaPulldown prediction {pred_dir}: {e}")
        return None


def load_prediction_data(pred_dir: Path) -> Optional[Dict]:
    """
    Load summary data for a single prediction directory.
    Handles AF3 local pipeline, AF3 Server flat, AF3 Server nesting,
    and AlphaPulldown (AF2-multimer) formats automatically.
    """
    try:
        # Resolve AF3 Server extra nesting
        pred_dir = resolve_prediction_dir(pred_dir)
        pred_name = pred_dir.name

        # --- AF3 Server flat format ---
        if _is_af3_server_flat(pred_dir):
            return _load_server_flat_prediction(pred_dir)

        # --- AlphaPulldown / AF2-multimer ---
        af2_wd = _af2_work_dir(pred_dir)
        if af2_wd is not None:
            return _load_alphapulldown_prediction(pred_dir, af2_wd)

        # Find summary_confidences.json (file prefix differs from the
        # folder name when the folder was renamed)
        prefix = prediction_file_prefix(pred_dir)
        summary_file = pred_dir / f"{prefix}_summary_confidences.json"
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

        ranking_file = pred_dir / f"{prefix}_ranking_scores.csv"
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
            'format': 'af3',
        }
    except Exception as e:
        print(f"Error loading {pred_dir}: {e}")
        return None


def _scan_entry(d: Path, kind: str) -> Dict:
    """Scan entry for prediction folder d of format kind."""
    if kind in ('af2', 'alphapulldown'):
        return _af2_entry(d, _af2_work_dir(d))
    name = d.name
    name_part = name[5:] if kind == 'server_flat' and name.startswith('fold_') else name
    m = re.match(r'(.+)_and_(.+)', name_part)
    bait, prey = (m.group(1), m.group(2)) if m else (name_part, "")
    return {'name': name, 'bait': bait, 'prey': prey,
            'prediction_dir': str(d), 'format': kind}


class AF3Scanner:
    """Fast scanner: lists directories instantly, loads data on demand."""

    def __init__(self, predictions_folder: Path):
        self.af3_folder = Path(predictions_folder)
        self.predictions = []
        self.skipped = []
        self.duplicates = []

    def scan(self) -> List[Dict]:
        """
        Find every prediction folder at or below the predictions folder (up to
        4 levels down; the folder itself may be one prediction). Returns dicts
        with name, bait, prey, prediction_dir, format. Folders directly below
        the predictions folder that hold no prediction output (logs, MSAs,
        unfinished runs) are listed in self.skipped.
        """
        if not self.af3_folder.exists():
            raise ValueError(f"Predictions folder not found: {self.af3_folder}")

        self.predictions = []
        self.skipped = []
        found = find_prediction_dirs(self.af3_folder)
        for d, kind in found:
            self.predictions.append(_scan_entry(d, kind))

        # Results are keyed by prediction name, so two folders with the same
        # name (e.g. v1/pairA/ and v2/pairA/) would be mixed up: report them.
        seen: Dict[str, str] = {}
        self.duplicates = []
        for p in self.predictions:
            if p['name'] in seen:
                self.duplicates.append((seen[p['name']], p['prediction_dir']))
            else:
                seen[p['name']] = p['prediction_dir']

        if not (found and found[0][0] == self.af3_folder):
            covered = {d.relative_to(self.af3_folder).parts[0] for d, _ in found}
            self.skipped = [c.name for c in _subdirs(self.af3_folder) if c.name not in covered]
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
