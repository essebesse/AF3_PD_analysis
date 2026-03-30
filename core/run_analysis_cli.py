"""
Standalone CLI for running AF3 analysis (used by SLURM jobs).
This script can be run independently of Streamlit for batch execution.
"""

import argparse
import json
import sys
from pathlib import Path
from multiprocessing import Pool, cpu_count

# Ensure the project root is on sys.path so 'core.*' imports work
# regardless of working directory (e.g., when spawned by SLURM)
_project_root = str(Path(__file__).resolve().parent.parent)
if _project_root not in sys.path:
    sys.path.insert(0, _project_root)


def find_prediction_dirs(af3_folder: Path) -> list:
    """Find all prediction directories (skip seed-* and hidden dirs)."""
    pred_dirs = []
    for d in sorted(af3_folder.iterdir()):
        if d.is_dir() and not d.name.startswith('seed-') and not d.name.startswith('.'):
            pred_dirs.append(d)
    return pred_dirs


def analyze_single_prediction(pred_dir: Path, pae_cutoff: float = 10.0,
                               all_models: bool = False) -> list:
    """
    Analyze a single prediction directory.
    Default: top-ranked model only. --all_models: all 5 seed/sample models.
    """
    from core.analyzer import analyze_prediction_all_models
    return analyze_prediction_all_models(
        pred_dir, ipsae_pae_cutoff=pae_cutoff, top_only=(not all_models)
    )


def main():
    parser = argparse.ArgumentParser(
        description='Run AF3 pulldown analysis on all predictions in a directory'
    )
    parser.add_argument('--input_dir', required=True,
                        help='Path to AF3 directory (e.g., Q9BW83_IFT27/AF3)')
    parser.add_argument('--cpus', type=int, default=32,
                        help='Number of CPUs to use (default: 32)')
    parser.add_argument('--pae_cutoff', type=float, default=10.0,
                        help='PAE cutoff for ipSAE calculation (default: 10.0)')
    parser.add_argument('--all_models', action='store_true',
                        help='Analyze all 5 models per prediction (default: top-ranked only)')
    parser.add_argument('--output_file', type=str, default=None,
                        help='Output JSON file (default: af3_app_all_models_analysis.json in input_dir)')

    args = parser.parse_args()

    # Resolve paths
    input_path = Path(args.input_dir)
    if not input_path.is_absolute():
        input_path = Path.cwd() / input_path

    if not input_path.exists():
        print(f"Error: Input directory not found: {input_path}")
        sys.exit(1)

    # Determine output file
    if args.output_file:
        output_file = Path(args.output_file)
    else:
        output_file = input_path / "af3_app_all_models_analysis.json"

    # Find all prediction directories
    pred_dirs = find_prediction_dirs(input_path)
    total = len(pred_dirs)

    if total == 0:
        print("No prediction directories found.")
        sys.exit(0)

    mode = "all 5 models" if args.all_models else "top-ranked model only"
    print(f"Found {total} predictions to analyze ({mode})")
    print(f"Using {args.cpus} CPUs")
    print(f"PAE cutoff: {args.pae_cutoff}")

    # Cap CPUs at available cores
    actual_cpus = min(args.cpus, cpu_count())
    print(f"Using {actual_cpus} CPUs (capped at available cores)")

    # Run analysis with multiprocessing
    all_results = []
    total_models = 0

    from functools import partial
    analyze_fn = partial(analyze_single_prediction,
                         pae_cutoff=args.pae_cutoff,
                         all_models=args.all_models)

    with Pool(actual_cpus) as pool:
        for i, result_list in enumerate(pool.imap_unordered(
            analyze_fn,
            pred_dirs
        )):
            if result_list:
                # Strip sequences before storing (they're huge)
                for r in result_list:
                    if 'sequences' in r:
                        del r['sequences']
                all_results.extend(result_list)
                total_models += len(result_list)

                # Save per-prediction JSON for fast lookup in detailed analysis
                # Use resolve_prediction_dir to handle AF3 Server nesting (name/name/)
                pred_name = result_list[0].get('prediction_name', '')
                if pred_name:
                    from core.scanner import resolve_prediction_dir
                    target_dir = resolve_prediction_dir(input_path / pred_name)
                    pred_json = target_dir / "af3_app_analysis.json"
                    try:
                        with open(pred_json, 'w') as pf:
                            json.dump(result_list, pf, indent=2)
                    except Exception:
                        pass  # Non-critical - big JSON is the fallback

            # Progress update
            print(f"Progress: {i+1}/{total} predictions analyzed ({total_models} models)")

    # Save combined results (big JSON)
    with open(output_file, 'w') as f:
        json.dump(all_results, f, indent=2)

    print(f"\nAnalysis complete!")
    print(f"Total models analyzed: {total_models}")
    print(f"Results saved to: {output_file}")
    print(f"Per-prediction JSON files saved to each prediction directory")


if __name__ == '__main__':
    main()
