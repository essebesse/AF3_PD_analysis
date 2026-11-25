"""Compare the app's ipSAE with the reference implementation of Dunbrack (2025).

Get the reference script first:
    curl -LO https://raw.githubusercontent.com/DunbrackLab/IPSAE/main/ipsae.py

Usage (from the app folder):
    python validation/ipsae_vs_reference.py /path/to/ipsae.py [--top] [folder ...]

A folder is an AF3 local prediction or a folder of them. Without folders, the
bundled example_data/IFT27_pulldown/AF3 is used. Every model with
confidences + structure files is scored by both implementations (PAE cutoff
10 Å), for every pair of polymer chains (ipsae.py reports all chain pairs).
--top: top-ranked model only. Inputs are copied to a temporary folder
because ipsae.py writes its output next to its inputs.
"""
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

APP = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(APP))

import numpy as np  # noqa: E402
from core.analyzer import analyze_model_pairs, calculate_ipsae, chain_token_indices  # noqa: E402
from core.scanner import prediction_file_prefix  # noqa: E402

CUTOFF = 10


def reference_ipsae(ref_script: Path, pae_file: Path, struct_file: Path) -> dict:
    """{(chain1, chain2): ipSAE} for every chain pair, symmetric ('max') rows,
    both orders as keys."""
    subprocess.run([sys.executable, str(ref_script), str(pae_file), str(struct_file),
                    str(CUTOFF), str(CUTOFF)], check=True, capture_output=True)
    out = Path(f"{struct_file.with_suffix('')}_{CUTOFF:02d}_{CUTOFF:02d}.txt")
    rows = [line.split() for line in out.read_text().splitlines() if line.strip()]
    head = next(r for r in rows if 'Type' in r)
    i_c1, i_c2 = head.index('Chn1'), head.index('Chn2')
    i_type, i_val = head.index('Type'), head.index('ipSAE')
    vals = {}
    for r in rows:
        if r is not head and len(r) > i_val and r[i_type] == 'max':
            vals[(r[i_c1], r[i_c2])] = vals[(r[i_c2], r[i_c1])] = float(r[i_val])
    return vals


def af3_models(pair: Path, top_only: bool = False):
    """(label, confidences.json, model.cif) for each model of an AF3 local prediction."""
    name = prediction_file_prefix(pair)   # differs from the folder name after a rename
    for d in ([] if top_only else sorted(pair.glob('seed-*_sample-*'))):
        conf = d / f"{name}_{d.name}_confidences.json"
        cif = d / f"{name}_{d.name}_model.cif"
        if not conf.exists():
            conf, cif = d / 'confidences.json', d / 'model.cif'
        if conf.exists() and cif.exists():
            yield d.name, conf, cif
    conf, cif = pair / f"{name}_confidences.json", pair / f"{name}_model.cif"
    if conf.exists() and cif.exists():
        yield 'top', conf, cif


def main():
    args = [a for a in sys.argv[1:] if a != '--top']
    top_only = '--top' in sys.argv
    if not args:
        sys.exit(__doc__)
    ref = Path(args[0]).resolve()
    folders = [Path(p) for p in args[1:]] or [APP / 'example_data/IFT27_pulldown/AF3']

    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        for folder in folders:
            is_pred = any(folder.glob('*ranking_scores.csv')) or any(folder.glob('seed-*_sample-*'))
            preds = [folder] if is_pred else sorted(p for p in folder.iterdir() if p.is_dir())
            for pair in preds:
                for label, conf, cif in af3_models(pair, top_only):
                    work = Path(tmp) / pair.name / label
                    work.mkdir(parents=True)
                    conf_c = Path(shutil.copy(conf, work))
                    cif_c = Path(shutil.copy(cif, work))
                    c = json.load(open(conf_c))
                    app = analyze_model_pairs(cif_c, c, None, np.array(c['pae']),
                                              pair.name, None, float(CUTOFF))
                    refv = reference_ipsae(ref, conf_c, cif_c)
                    tok = chain_token_indices(c)
                    for p in app['pairs']:
                        a, b = p['chains'].split(':')
                        # stored values are rounded to 4 decimals; compare unrounded
                        v = calculate_ipsae(np.array(c['pae']), tok[a], tok[b], float(CUTOFF))
                        assert round(v, 4) == p['ipsae'], (p['chains'], v, p['ipsae'])
                        rows.append((pair.name, label, p['chains'], v, refv.get((a, b))))

    print(f"{'prediction':32} {'model':22} {'pair':5} {'app':>8} {'reference':>10} {'diff':>9}")
    for pair, label, ch, av, rv in rows:
        rs = f"{rv:10.4f}" if rv is not None else f"{'missing':>10}"
        ds = f"{av - rv:9.2e}" if rv is not None else ""
        print(f"{pair:32} {label:22} {ch:5} {av:8.4f} {rs} {ds}")
    found = [r for r in rows if r[4] is not None]
    if found:
        print(f"\n{len(rows)} chain pairs ({len(rows) - len(found)} missing in reference), "
              f"max |difference| = {max(abs(r[3] - r[4]) for r in found):.2e}")


if __name__ == '__main__':
    main()
