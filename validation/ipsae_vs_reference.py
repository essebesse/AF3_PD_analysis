"""Compare the app's ipSAE with the reference implementation of Dunbrack (2025).

Get the reference script first:
    curl -LO https://raw.githubusercontent.com/DunbrackLab/IPSAE/main/ipsae.py

Usage (from the app folder):
    python validation/ipsae_vs_reference.py /path/to/ipsae.py [AF3_folder ...]

Without AF3 folders, the bundled example_data/IFT27_pulldown/AF3 is used.
Every model with confidences + structure files is scored by both
implementations (PAE cutoff 10 Å). Inputs are copied to a temporary folder
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
from core.analyzer import calculate_ipsae, extract_chain_indices  # noqa: E402

CUTOFF = 10


def reference_ipsae(ref_script: Path, pae_file: Path, struct_file: Path) -> float:
    subprocess.run([sys.executable, str(ref_script), str(pae_file), str(struct_file),
                    str(CUTOFF), str(CUTOFF)], check=True, capture_output=True)
    out = Path(f"{struct_file.with_suffix('')}_{CUTOFF:02d}_{CUTOFF:02d}.txt")
    rows = [line.split() for line in out.read_text().splitlines() if line.strip()]
    i_type, i_val = rows[0].index('Type'), rows[0].index('ipSAE')
    return next(float(r[i_val]) for r in rows[1:] if r[i_type] == 'max')


def af3_models(pair: Path):
    """(label, confidences.json, model.cif) for each model of an AF3 local prediction."""
    name = pair.name
    for d in sorted(pair.glob('seed-*_sample-*')):
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
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    ref = Path(sys.argv[1]).resolve()
    folders = [Path(p) for p in sys.argv[2:]] or [APP / 'example_data/IFT27_pulldown/AF3']

    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        for folder in folders:
            for pair in sorted(p for p in folder.iterdir() if p.is_dir()):
                for label, conf, cif in af3_models(pair):
                    work = Path(tmp) / pair.name / label
                    work.mkdir(parents=True)
                    conf_c = Path(shutil.copy(conf, work))
                    cif_c = Path(shutil.copy(cif, work))
                    c = json.load(open(conf_c))
                    _, a, b, _ = extract_chain_indices(c)
                    app_v = calculate_ipsae(np.array(c['pae']), a, b, float(CUTOFF))
                    rows.append((pair.name, label, app_v, reference_ipsae(ref, conf_c, cif_c)))

    print(f"{'prediction':28} {'model':16} {'app':>8} {'reference':>10} {'diff':>9}")
    for pair, label, av, rv in rows:
        print(f"{pair:28} {label:16} {av:8.4f} {rv:10.4f} {av - rv:9.2e}")
    if rows:
        print(f"\n{len(rows)} models, max |difference| = {max(abs(r[2] - r[3]) for r in rows):.2e}")


if __name__ == '__main__':
    main()
