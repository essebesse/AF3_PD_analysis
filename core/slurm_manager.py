"""
SLURM job management module for AF3 analysis.
Extracted from pages/batch_execution.py for cleaner organization.
"""

import subprocess
from pathlib import Path
from typing import Optional, Dict


def submit_job(script_content: str, working_dir: str) -> Dict:
    """
    Submit SLURM job.

    Args:
        script_content: Content of the SLURM script
        working_dir: Working directory for the job submission

    Returns:
        Dict with {'success': bool, 'job_id': str or None, 'error': str or None}
    """
    import tempfile

    try:
        # Save script to a temporary file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.slurm', delete=False) as f:
            f.write(script_content)
            temp_script = f.name

        # Submit job
        result = subprocess.run(
            ['sbatch', temp_script],
            capture_output=True,
            text=True,
            cwd=working_dir,
            timeout=60
        )

        # Clean up temp file
        Path(temp_script).unlink(missing_ok=True)

        if result.returncode == 0:
            # Parse job ID from output (format: "Submitted batch job <id>")
            job_id = result.stdout.strip().split()[-1]
            return {'success': True, 'job_id': job_id, 'error': None}
        else:
            return {'success': False, 'job_id': None, 'error': result.stderr}

    except Exception as e:
        return {'success': False, 'job_id': None, 'error': str(e)}


def check_job_status(job_id: str) -> Dict:
    """
    Check SLURM job status (non-blocking, single check).

    Args:
        job_id: SLURM job ID

    Returns:
        Dict with {'status': str, 'running': bool, 'completed': bool}
    """
    try:
        result = subprocess.run(
            ['squeue', '-j', job_id, '--noheader', '-o', '%T'],
            capture_output=True,
            text=True,
            timeout=30
        )

        status = result.stdout.strip()

        if status == 'RUNNING':
            return {'status': 'RUNNING', 'running': True, 'completed': False}
        elif status == 'PENDING':
            return {'status': 'PENDING', 'running': False, 'completed': False}
        elif status in ('FAILED', 'CANCELLED', 'TIMEOUT', 'OUT_OF_MEMORY', 'NODE_FAIL'):
            return {'status': status, 'running': False, 'completed': True}
        elif not status:
            # Job no longer in queue — likely completed, failed, or cancelled.
            # Use sacct for definitive status if available.
            try:
                sacct = subprocess.run(
                    ['sacct', '-j', job_id, '--noheader', '-o', 'State', '-P'],
                    capture_output=True, text=True, timeout=15
                )
                sacct_status = sacct.stdout.strip().split('\n')[0].strip() if sacct.stdout.strip() else ''
                if sacct_status:
                    done = sacct_status not in ('RUNNING', 'PENDING')
                    return {'status': sacct_status, 'running': not done, 'completed': done}
            except Exception:
                pass
            return {'status': 'COMPLETED', 'running': False, 'completed': True}
        else:
            return {'status': status, 'running': status in ['RUNNING', 'PENDING', 'CONFIGURING'], 'completed': False}

    except Exception as e:
        return {'status': 'ERROR', 'running': False, 'completed': False, 'error': str(e)}


def get_job_output(job_id: str, working_dir: str) -> Optional[str]:
    """
    Read job stdout file if it exists.

    Args:
        job_id: SLURM job ID
        working_dir: Directory where job output files are located

    Returns:
        Job output content as string, or None if not found
    """
    # Try both naming patterns (chunk logs and legacy)
    import glob
    working = Path(working_dir)
    candidates = list(working.glob(f"AF3_app_chunk*_{job_id}.log")) + \
                 [working / f"AF3_analysis_{job_id}.out"]
    for f in candidates:
        if f.exists():
            try:
                return f.read_text()
            except Exception:
                return None
    return None


def cancel_job(job_id: str) -> Dict:
    """
    Cancel a running SLURM job.

    Args:
        job_id: SLURM job ID

    Returns:
        Dict with {'success': bool, 'error': str or None}
    """
    try:
        result = subprocess.run(
            ['scancel', job_id],
            capture_output=True,
            text=True,
            timeout=30
        )

        if result.returncode == 0:
            return {'success': True, 'error': None}
        else:
            return {'success': False, 'error': result.stderr}

    except Exception as e:
        return {'success': False, 'error': str(e)}