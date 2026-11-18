"""
Site-specific SLURM settings.

Read from ``slurm_config.toml`` in the app folder (not tracked by git; copy
``slurm_config.example.toml`` to start). Every key is optional. Empty strings
mean "leave out this #SBATCH line and use the cluster default". Invalid
values are replaced by the defaults below.
"""

import shlex
import shutil
from pathlib import Path
from typing import Dict

try:                                # Python 3.11+
    import tomllib as _toml
except ImportError:                 # Python <= 3.10: optional `pip install tomli`
    try:
        import tomli as _toml
    except ImportError:
        _toml = None

APP_DIR = Path(__file__).resolve().parent.parent
CONFIG_FILE = APP_DIR / "slurm_config.toml"

MAX_JOBS = 50                       # upper bound of the chunk / concurrency inputs

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
    "env_setup": "source " + shlex.quote(str(APP_DIR / ".venv" / "bin" / "activate")),
    # Extra raw lines placed after the #SBATCH header, e.g.
    # ["#SBATCH --constraint=avx2", "module load python/3.11"]
    "extra_lines": [],
}


def slurm_available() -> bool:
    """True when ``sbatch`` is on PATH."""
    return shutil.which("sbatch") is not None


def _valid(key: str, value) -> bool:
    if key in ("default_jobs", "max_concurrent"):
        return isinstance(value, int) and not isinstance(value, bool) and 1 <= value <= MAX_JOBS
    if key == "default_cpus":
        return isinstance(value, int) and not isinstance(value, bool) and value >= 1
    if key == "cpus_options":
        return (isinstance(value, list) and len(value) > 0
                and all(isinstance(v, int) and not isinstance(v, bool) and v >= 1 for v in value))
    if key == "memory_options":
        return isinstance(value, list) and len(value) > 0 and all(isinstance(v, str) and v for v in value)
    if key == "extra_lines":
        return isinstance(value, list) and all(isinstance(v, str) for v in value)
    if key == "env_setup":
        return isinstance(value, str) and value.strip() != ""
    return isinstance(value, str)       # partition, qos, account, time, default_memory


def load_slurm_config() -> Dict:
    cfg = dict(DEFAULTS)
    if not CONFIG_FILE.is_file():
        return cfg
    if _toml is None:
        print(f"Ignoring {CONFIG_FILE.name}: reading it needs Python 3.11+ or `pip install tomli`")
        return cfg
    try:
        with open(CONFIG_FILE, "rb") as f:
            user = _toml.load(f)
    except (OSError, ValueError) as e:      # TOMLDecodeError is a ValueError
        print(f"Ignoring {CONFIG_FILE.name}: {e}")
        return cfg
    for k, v in user.items():
        if k not in DEFAULTS:
            print(f"{CONFIG_FILE.name}: unknown key '{k}' ignored")
        elif _valid(k, v):
            cfg[k] = v
        else:
            print(f"{CONFIG_FILE.name}: invalid value for '{k}' ({v!r}); using default {DEFAULTS[k]!r}")
    # Defaults must be one of the offered options
    if cfg["default_cpus"] not in cfg["cpus_options"]:
        cfg["default_cpus"] = cfg["cpus_options"][0]
    if cfg["default_memory"] not in cfg["memory_options"]:
        cfg["default_memory"] = cfg["memory_options"][0]
    return cfg
