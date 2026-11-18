"""
Site-specific SLURM settings.

Read from ``slurm_config.toml`` in the app folder (not tracked by git; copy
``slurm_config.example.toml`` to start). Every key is optional. Empty strings
mean "leave out this #SBATCH line and use the cluster default".
"""

import shutil
import tomllib
from pathlib import Path
from typing import Dict

APP_DIR = Path(__file__).resolve().parent.parent
CONFIG_FILE = APP_DIR / "slurm_config.toml"

DEFAULTS: Dict = {
    "partition": "",
    "qos": "",
    "account": "",
    "time": "23:59:59",
    "default_jobs": 4,          # number of chunks (array tasks)
    "max_concurrent": 2,        # array throttle: %N
    "cpus_options": [8, 16, 32, 64],
    "default_cpus": 16,
    "memory_options": ["16G", "32G", "64G", "128G"],
    "default_memory": "32G",
    # Shell line that makes the app's Python environment available on the
    # compute node. Default: the app's own .venv.
    "env_setup": f"source '{APP_DIR}/.venv/bin/activate'",
    # Extra raw lines placed after the #SBATCH header, e.g.
    # ["#SBATCH --constraint=avx2", "module load python/3.11"]
    "extra_lines": [],
}


def slurm_available() -> bool:
    """True when ``sbatch`` is on PATH."""
    return shutil.which("sbatch") is not None


def load_slurm_config() -> Dict:
    cfg = dict(DEFAULTS)
    if CONFIG_FILE.is_file():
        try:
            with open(CONFIG_FILE, "rb") as f:
                user = tomllib.load(f)
            cfg.update({k: v for k, v in user.items() if k in DEFAULTS})
        except (OSError, tomllib.TOMLDecodeError) as e:
            print(f"Ignoring {CONFIG_FILE.name}: {e}")
    return cfg
