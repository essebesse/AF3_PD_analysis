"""Check that the app finds and analyses predictions in different folder layouts.

Builds each layout in a temporary folder on local disk (/tmp) from the
example data, runs the same discovery and analysis as the app
(find_af3_projects_recursive -> find_predictions_folder -> AF3Scanner.scan ->
analyze_prediction_all_models) and compares the number of predictions and
models found with what the layout holds.

    python validation/folder_layouts.py [--af2 /path/to/af2_prediction]

--af2: an AF2 prediction folder (ranked_N.pdb, ranking_debug.json, optional
ranked_N_PAE.csv, in the folder or in seq/) for the AF2 layouts; skipped
without it.
"""
import argparse
import shutil
import sys
import tempfile
from pathlib import Path

APP = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP))
from core.scanner import AF3Scanner, find_af3_projects_recursive, find_predictions_folder  # noqa: E402
from core.analyzer import analyze_prediction_all_models  # noqa: E402

EX = APP / 'example_data/IFT27_pulldown/AF3'
A = 'q9bw83_and_q9y547'   # 5 models (seed-1_sample-0..4)
B = 'q9bw83_and_o75952'   # top model only
MODELS = {A: 5, B: 1}


def put(src_name, dest: Path, rename=None):
    """Copy one example prediction to dest (folder renamed to `rename`, files keep their prefix)."""
    out = dest / (rename or src_name)
    shutil.copytree(EX / src_name, out)
    return out


def af2_copy(src: Path, dest: Path, flat: bool):
    wd = src / 'seq' if (src / 'seq' / 'ranking_debug.json').exists() else src
    out = dest / src.name
    target = out if flat else out / 'seq'
    target.mkdir(parents=True)
    for f in wd.iterdir():
        if f.name == 'ranking_debug.json' or (f.name.startswith('ranked_') and f.suffix in ('.pdb', '.csv')):
            shutil.copy2(f, target / f.name)
    return out


def layouts(root: Path, af2: Path = None):
    """Yield (name, folder the user selects, expected predictions, expected models)."""
    def new(name):
        p = root / name
        p.mkdir()
        return p

    p = new('standard')                       # project/AF3/<pred>/
    put(A, p / 'AF3'); put(B, p / 'AF3')
    yield 'project/AF3/<pred>/', p, 2, 6

    p = new('flat')                           # project/<pred>/
    put(A, p); put(B, p)
    yield 'project/<pred>/', p, 2, 6

    p = new('server_nesting')                 # project/<pred>/<pred>/
    put(A, p / A); put(B, p / B)
    yield 'project/<pred>/<pred>/', p, 2, 6

    p = new('run_folders')                    # project/<pred>/AF3/<pred>/ + msa/
    for n in (A, B):
        put(n, p / n / 'AF3')
        (p / n / 'msa' / n).mkdir(parents=True)
        shutil.copy2(EX / n / f'{n}_data.json', p / n / 'msa' / n)
    yield 'project/<pred>/AF3/<pred>/ (+ msa/)', p, 2, 6

    p = new('deeper')                         # project/<pred>/output/run1/<pred>/
    for n in (A, B):
        put(n, p / n / 'output' / 'run1')
    yield 'project/<pred>/output/run1/<pred>/', p, 2, 6

    p = new('renamed')                        # folder renamed, files keep the AF3 job name
    put(A, p, rename='job_1'); put(B, p, rename='job_2')
    yield 'project/<renamed>/ (files keep job name)', p, 2, 6

    p = new('junk_siblings')                  # other folders next to the predictions
    put(A, p / 'AF3'); put(B, p / 'AF3')
    for j in ('logs', 'msas', 'slurm_out'):
        (p / 'AF3' / j).mkdir()
        (p / 'AF3' / j / 'x.txt').write_text('x')
    yield 'project/AF3/{<pred>, logs/, msas/}', p, 2, 6

    p = new('batches')                        # predictions grouped in sub-folders
    put(A, p / 'batch1'); put(B, p / 'batch2')
    yield 'project/{batch1,batch2}/<pred>/', p, 2, 6

    p = new('single')                         # user selects one prediction folder
    q = put(A, p)
    yield 'one prediction folder selected', q, 1, 5

    p = new('unfinished')                     # a prediction folder with no output yet
    put(A, p / 'AF3'); put(B, p / 'AF3')
    (p / 'AF3' / 'q9bw83_and_p00000').mkdir()
    (p / 'AF3' / 'q9bw83_and_p00000' / 'q9bw83_and_p00000_data.json').write_text('{}')
    yield 'project/AF3/ with an unfinished prediction', p, 2, 6

    if af2:
        p = new('af2_seq')
        af2_copy(af2, p, flat=False)
        yield 'AF2 project/<pair>/seq/', p, 1, 5
        p = new('af2_flat')
        af2_copy(af2, p, flat=True)
        yield 'AF2 project/<pair>/', p, 1, 5
        p = new('mixed')
        af2_copy(af2, p, flat=False); put(A, p)
        yield 'AF3 and AF2 in one folder', p, 2, 10


def run(sel: Path):
    """The app's path: project browser -> predictions folder -> scan -> analysis."""
    projects = find_af3_projects_recursive(str(sel.parent), max_depth=6)
    # listed: the "Scan for AF3 subprojects" button offers this folder or a
    # project inside it (batch1/, batch2/ are offered separately)
    proj = [x for x in projects if Path(x['project_path']) == sel or sel in Path(x['project_path']).parents
            or Path(x['project_path']) in sel.parents]
    folder = find_predictions_folder(str(sel))   # folder typed or selected by the user
    preds = AF3Scanner(Path(folder)).scan()
    n_models = 0
    for pr in preds:
        n_models += len(analyze_prediction_all_models(Path(pr['prediction_dir'])))
    return bool(proj), len(preds), n_models


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--af2', type=Path)
    args = ap.parse_args()
    fails = 0
    with tempfile.TemporaryDirectory(prefix='pdr_layouts_') as tmp:
        for name, sel, exp_p, exp_m in layouts(Path(tmp), args.af2):
            # each layout in its own parent so the project browser sees only it
            listed, n_p, n_m = run(sel)
            ok = listed and n_p == exp_p and n_m == exp_m
            fails += not ok
            print(f"{'ok  ' if ok else 'FAIL'} {name:45s} listed={listed!s:5s} "
                  f"predictions {n_p}/{exp_p}  models {n_m}/{exp_m}")
    print(f"{fails} failed")
    sys.exit(1 if fails else 0)


if __name__ == '__main__':
    main()
