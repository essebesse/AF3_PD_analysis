"""
PulldownRanker - Main Entry Point

A Streamlit-based GUI for analyzing AlphaFold3 pulldown prediction data.
Provides a modern interface for running analysis scripts, viewing results
for ALL models (not just top-ranked), and performing detailed interface analysis.
"""

import streamlit as st
import os
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

# Page configuration (must be first)
st.set_page_config(
    page_title="PulldownRanker",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# Custom CSS
st.markdown("""
<style>
    .main-header {
        font-size: 2rem;
        font-weight: bold;
        color: #1e40af;
        margin-bottom: 1rem;
    }
</style>
""", unsafe_allow_html=True)

# App title
st.markdown('<div class="main-header">🔬 PulldownRanker</div>', unsafe_allow_html=True)
st.markdown("### Web-based analysis for AlphaFold 3 protein interaction predictions")

# Check if core modules exist
try:
    from core.scanner import AF3Scanner
except ImportError as e:
    st.error(f"Failed to import core modules: {e}")
    st.stop()

# ── Native file dialogs ───────────────────────────────────────────────────────
# Each dialog runs in a separate process: tkinter must own the main thread, and
# Streamlit runs this script in a worker thread, so calling it in-process can
# hang or crash the app. The dialog opens on the machine running Streamlit,
# which is the user's own computer when the app is started locally (./run.sh).
def pick_folder_dialog(start_dir):
    """Open the operating system's native folder picker.

    Linux: zenity (GNOME) or kdialog (KDE). macOS: Finder via osascript.
    Windows, or Linux without either tool: tkinter (native on Windows).
    Returns (path or None, error or None); cancelling gives (None, None).
    """
    title = "Select project folder"
    # Fallback dialog: tkinter is in the standard library, but on Linux it is a
    # system package (python3-tk), not a pip one, so it may be missing.
    tk_script = (
        "import sys, tkinter as tk\n"
        "from tkinter import filedialog\n"
        "root = tk.Tk(); root.withdraw(); root.attributes('-topmost', True)\n"
        "print(filedialog.askdirectory(initialdir=sys.argv[1], title=sys.argv[2]))\n"
    )
    tk_cmd = [sys.executable, "-c", tk_script, start_dir, title]

    # Pick the dialog that looks native on this system
    if sys.platform == "darwin":
        apple = (f'POSIX path of (choose folder with prompt "{title}" '
                 f'default location (POSIX file "{start_dir}"))')
        cmd = ["osascript", "-e", apple]
    elif sys.platform.startswith("linux") and shutil.which("zenity"):
        # The GTK4 folder chooser (zenity 4) lists files greyed out; a filter
        # that no file name can match ('*/') hides them, folders stay visible.
        cmd = ["zenity", "--file-selection", "--directory",
               f"--title={title}", f"--filename={start_dir.rstrip('/')}/",
               "--file-filter=Folders | */"]
    elif sys.platform.startswith("linux") and shutil.which("kdialog"):
        cmd = ["kdialog", "--getexistingdirectory", start_dir, "--title", title]
    else:
        cmd = tk_cmd

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    except (subprocess.TimeoutExpired, OSError) as e:
        return None, f"Could not open folder dialog: {e}"
    if result.returncode != 0:
        # zenity, kdialog and osascript exit with 1 when the user cancels
        if cmd is not tk_cmd and result.returncode == 1:
            return None, None
        lines = result.stderr.strip().splitlines()
        detail = lines[-1] if lines else ""
        return None, f"Could not open folder dialog. {detail}"
    # osascript ends folder paths with '/'; drop it so the path matches typed ones
    path = result.stdout.strip()
    if len(path) > 1:
        path = path.rstrip("/")
    return (path or None), None


def pick_zip_files_dialog(start_dir):
    """Open the native file picker for one or more .zip files.

    Same backends as pick_folder_dialog(). Returns (list of paths, error or None);
    cancelling gives ([], None).
    """
    title = "Select .zip files"
    # Fallback dialog, as in pick_folder_dialog(); prints one path per line
    tk_script = (
        "import sys, tkinter as tk\n"
        "from tkinter import filedialog\n"
        "root = tk.Tk(); root.withdraw(); root.attributes('-topmost', True)\n"
        "files = filedialog.askopenfilenames(initialdir=sys.argv[1], title=sys.argv[2],\n"
        "    filetypes=[('ZIP archives', '*.zip *.ZIP')])\n"
        "print('\\n'.join(files))\n"
    )
    tk_cmd = [sys.executable, "-c", tk_script, start_dir, title]

    # Pick the dialog that looks native on this system. Only .zip files are
    # listed; folders stay visible so the user can navigate to them.
    if sys.platform == "darwin":
        # AppleScript returns a list of aliases; join their POSIX paths by line
        apple = (
            f'set fs to choose file with prompt "{title}" of type {{"zip"}} '
            f'default location (POSIX file "{start_dir}") with multiple selections allowed\n'
            'set out to ""\n'
            'repeat with f in fs\n'
            'set out to out & POSIX path of f & linefeed\n'
            'end repeat\n'
            'return out'
        )
        cmd = ["osascript", "-e", apple]
    elif sys.platform.startswith("linux") and shutil.which("zenity"):
        cmd = ["zenity", "--file-selection", "--multiple", "--separator=\n",
               f"--title={title}", f"--filename={start_dir.rstrip('/')}/",
               "--file-filter=ZIP archives | *.zip *.ZIP"]
    elif sys.platform.startswith("linux") and shutil.which("kdialog"):
        cmd = ["kdialog", "--getopenfilename", start_dir, "*.zip *.ZIP|ZIP archives",
               "--multiple", "--separate-output", "--title", title]
    else:
        cmd = tk_cmd

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    except (subprocess.TimeoutExpired, OSError) as e:
        return [], f"Could not open file dialog: {e}"
    if result.returncode != 0:
        # zenity, kdialog and osascript exit with 1 when the user cancels
        if cmd is not tk_cmd and result.returncode == 1:
            return [], None
        lines = result.stderr.strip().splitlines()
        detail = lines[-1] if lines else ""
        return [], f"Could not open file dialog. {detail}"
    return [p for p in result.stdout.splitlines() if p.strip()], None


def extract_zips(zip_paths):
    """Extract .zip files and return (project folder, messages, errors).

    One zip: extracted to <zip folder>/<zip name>/, which becomes the project.
    Several: each to <folder of the first zip>/unzipped_predictions/<zip name>/,
    and unzipped_predictions/ becomes the project. A target folder that already
    has content is reused, not overwritten.
    """
    # Several zips get a shared parent folder instead of being extracted next
    # to themselves: the zips usually sit in Downloads, and using that as the
    # project would make the scanner walk everything else stored there.
    zips = [Path(p) for p in zip_paths]
    if len(zips) == 1:
        project = zips[0].parent / zips[0].stem
        targets = [project]
    else:
        project = zips[0].parent / "unzipped_predictions"
        targets = [project / z.stem for z in zips]

    messages, errors = [], []
    for z, target in zip(zips, targets):
        # Reuse an earlier extraction: it may already hold the analysis cache
        if target.is_dir() and any(target.iterdir()):
            messages.append(f"`{z.name}`: already extracted, reusing `{target}`")
            continue
        try:
            with zipfile.ZipFile(z) as zf:
                target.mkdir(parents=True, exist_ok=True)
                zf.extractall(target)  # strips absolute paths and '..'
            messages.append(f"`{z.name}`: extracted to `{target}`")
        except (zipfile.BadZipFile, OSError) as e:
            errors.append(f"`{z.name}`: {e}")

    # No project when every zip failed before its folder was created
    if not project.is_dir():
        return None, messages, errors
    return str(project), messages, errors


# ── Project folder selection ──────────────────────────────────────────────────
# Initialize browser state — start from the directory the app was launched in
if 'browser_dir' not in st.session_state:
    st.session_state['browser_dir'] = os.getcwd()
if 'project_path' not in st.session_state:
    st.session_state['project_path'] = ""

project_path_current = st.session_state.get('project_path', "")

if project_path_current and os.path.isdir(project_path_current):
    # Project already selected — show compact path with option to change
    col_path, col_change = st.columns([4, 1])
    with col_path:
        st.text_input("Project folder:", value=project_path_current, disabled=True)
    with col_change:
        st.markdown("<br>", unsafe_allow_html=True)
        if st.button("Change"):
            st.session_state['project_path'] = ""
            st.session_state['browser_dir'] = project_path_current
            st.rerun()
else:
    # No project selected — show text input + browser
    typed_path = st.text_input(
        "Project folder:",
        value=project_path_current,
        placeholder="/path/to/your/project  (folder with AF3 prediction subdirectories)",
    )
    if typed_path:
        typed_path = os.path.expanduser(typed_path)
    if typed_path and os.path.isdir(typed_path):
        st.session_state['project_path'] = typed_path
        st.rerun()

    # Native dialogs: pick a folder, or one or more .zip files that are
    # extracted and then used as the project folder
    col_folder, col_zip, _ = st.columns([1, 1, 3])
    with col_folder:
        if st.button("📁 Select folder…"):
            chosen, err = pick_folder_dialog(st.session_state['browser_dir'])
            if err:
                st.error(err)
            elif chosen:
                st.session_state['project_path'] = chosen
                st.rerun()
    with col_zip:
        pick_zips = st.button("📦 Select .zip files…")
    # Handled outside the column so messages use the full page width
    if pick_zips:
        zip_paths, err = pick_zip_files_dialog(st.session_state['browser_dir'])
        if err:
            st.error(err)
        elif zip_paths:
            with st.spinner(f"Extracting {len(zip_paths)} .zip file(s)..."):
                project, messages, errors = extract_zips(zip_paths)
            if project and not errors:
                st.session_state.pop('_zip_result', None)
                st.session_state['project_path'] = project
                st.rerun()
            # Some zips failed: keep the report across reruns so the user can
            # read it and still continue with the zips that did extract.
            st.session_state['_zip_result'] = (project, messages, errors)

    # Report from a partly failed extraction, kept until the user continues
    zip_result = st.session_state.get('_zip_result')
    if zip_result:
        project, messages, errors = zip_result
        for e in errors:
            st.error(e)
        for m in messages:
            st.caption(m)
        if project and st.button(f"Continue with {project}", key="zip_continue"):
            st.session_state.pop('_zip_result', None)
            st.session_state['project_path'] = project
            st.rerun()

    st.markdown("##### Browse folders")
    browse_dir = Path(st.session_state['browser_dir'])
    if not browse_dir.is_dir():
        browse_dir = Path.home()
        st.session_state['browser_dir'] = str(browse_dir)

    nav_col1, nav_col2, nav_col3 = st.columns([1, 3, 1])
    with nav_col1:
        if st.button("⬆ Parent"):
            st.session_state['browser_dir'] = str(browse_dir.parent)
            st.rerun()
    with nav_col2:
        st.caption(f"📂 `{browse_dir}`")
    with nav_col3:
        if st.button("✓ Select this folder", type="primary"):
            st.session_state['project_path'] = str(browse_dir)
            st.rerun()

    try:
        subdirs = sorted(
            [d for d in browse_dir.iterdir() if d.is_dir() and not d.name.startswith('.')],
            key=lambda d: d.name.lower()
        )
    except PermissionError:
        subdirs = []
        st.warning("Permission denied.")

    if subdirs:
        cols = st.columns(4)
        for i, d in enumerate(subdirs[:60]):
            if cols[i % 4].button(d.name, key=f"dir_{i}", width="stretch"):
                st.session_state['browser_dir'] = str(d)
                st.rerun()
        if len(subdirs) > 60:
            st.caption(f"...and {len(subdirs) - 60} more folders")
    else:
        st.caption("(no subdirectories)")

    # ── Recursive AF3 subproject scanner ──────────────────────────────
    st.divider()
    st.markdown("##### Scan for AF3 subprojects")
    st.caption("Recursively search for all folders containing AF3 prediction data under the current directory.")

    if st.button("Scan for AF3 subprojects", key="scan_recursive"):
        from core.scanner import find_af3_projects_recursive
        with st.spinner(f"Scanning `{browse_dir}` ..."):
            found = find_af3_projects_recursive(str(browse_dir))
        st.session_state['recursive_scan_results'] = found
        st.session_state['recursive_scan_root'] = str(browse_dir)

    # Show results if we have them and we're still in the same browse_dir
    scan_results = st.session_state.get('recursive_scan_results')
    scan_root = st.session_state.get('recursive_scan_root', '')
    if scan_results is not None and scan_root == str(browse_dir):
        if scan_results:
            st.success(f"Found {len(scan_results)} AF3 project(s) under `{browse_dir.name}/`")
            for i, proj in enumerate(scan_results):
                col_label, col_btn = st.columns([4, 1])
                with col_label:
                    st.code(proj['label'], language=None)
                with col_btn:
                    if st.button("Select", key=f"scanpick_{i}", type="primary"):
                        st.session_state['project_path'] = proj['project_path']
                        st.session_state.pop('recursive_scan_results', None)
                        st.rerun()
        else:
            st.warning("No AF3 prediction data found under this directory.")

    st.divider()

project_path = st.session_state.get('project_path', "")

if not project_path:
    st.info("Enter a project folder path above, or use Browse to navigate your filesystem. "
            "The folder should contain AF3 prediction subdirectories (with .cif and .json files).")
    st.stop()

# Auto-detect where predictions live:
# 1) <project>/AF3/ if it exists and has predictions
# 2) <project>/ itself if it directly contains predictions
from core.scanner import find_predictions_folder
af3_folder = find_predictions_folder(project_path)
if not af3_folder:
    st.error(f"Folder not found: {project_path}")
    st.stop()

# Project folder changed since last render — wipe folder-scoped session state
# so we don't carry forward stale SLURM job IDs, selected predictions, or
# cached PAE/iface/viewer entries from the previous folder.
if st.session_state.get('_active_project') != project_path:
    folder_scoped_keys = {
        'selected_prediction', 'slurm_job_ids', 'slurm_num_chunks',
        '_slurm_merge_attempted_for',
    }
    folder_scoped_prefixes = (
        '_detail_', '_pae_plot_', '_zoom_plot_', '_iface_', '_comp_plot_',
    )
    for k in list(st.session_state.keys()):
        if k in folder_scoped_keys or k.startswith(folder_scoped_prefixes):
            del st.session_state[k]
    # Drop any st.cache_data caches keyed by the old folder (per-folder
    # scanners cache by path arg, so old entries are inert, but clearing
    # also gets rid of the legacy cache reads).
    try:
        st.cache_data.clear()
    except Exception:
        pass
    st.session_state['_active_project'] = project_path

# If a SLURM run finished while the user was on any other workflow step
# (Results, Detailed Analysis, Load Data), merge it now so they see fresh
# data the moment they navigate. The merge runs at top level so it fires
# regardless of which step is currently selected.
from pages.batch_execution import _auto_merge_if_needed as _auto_merge
if _auto_merge(af3_folder):
    st.rerun()

# Sidebar: app info
with st.sidebar:
    st.header("ℹ️ About")
    st.caption("Analysis settings are on the **Analyze** step.")


# Workflow step navigation (horizontal radio in main content area)
STEPS = ["1. Load Data", "2. Analyze", "3. Results", "4. Detailed Analysis"]

# Handle pending navigation from page buttons
if '_navigate_to' in st.session_state:
    target = st.session_state.pop('_navigate_to')
    if target in STEPS:
        st.session_state['workflow_step'] = target

step = st.radio("Workflow Step", STEPS, horizontal=True, key="workflow_step", label_visibility="collapsed")

# Main content based on selected step
if step == "1. Load Data":
    from pages.overview import show_load_data
    show_load_data(project_path, af3_folder)

elif step == "2. Analyze":
    from pages.batch_execution import show_analyze
    show_analyze(project_path, af3_folder)

elif step == "3. Results":
    from pages.results import show_results
    show_results(project_path, af3_folder)

elif step == "4. Detailed Analysis":
    from pages.detailed_analysis import show_detailed
    show_detailed(project_path, af3_folder)

# Footer
st.divider()
st.caption("Built with Streamlit | Based on AF3_PD_analysis_v4.py with ipSAE scoring")