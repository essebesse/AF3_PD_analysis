"""
AF3 Pulldown Analysis Web App - Main Entry Point

A Streamlit-based GUI for analyzing AlphaFold3 pulldown prediction data.
Provides a modern interface for running analysis scripts, viewing results
for ALL models (not just top-ranked), and performing detailed interface analysis.
"""

import streamlit as st
import os
from pathlib import Path

# Page configuration (must be first)
st.set_page_config(
    page_title="AF3 Pulldown Analysis",
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
st.markdown('<div class="main-header">🔬 AlphaFold Pulldown Analysis</div>', unsafe_allow_html=True)
st.markdown("### Web-based analysis for AlphaFold 3 and AlphaFold 2 protein interaction predictions")

# Pairwise-only notice: the app only scores the first two chains. Keep this
# compact but persistent so users with 3-chain predictions see it every time.
with st.expander("ℹ️ Scope: pairwise predictions only", expanded=False):
    st.markdown(
        "This app analyzes the interface between the **first two chains** of each "
        "prediction. Predictions with three or more chains are still listed, but "
        "their interface scoring is only Chain A vs Chain B — a ⚠️ marker flags "
        "those rows."
    )

# Check if core modules exist
try:
    from core.scanner import AF3Scanner
except ImportError as e:
    st.error(f"Failed to import core modules: {e}")
    st.stop()

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
    if typed_path and os.path.isdir(typed_path):
        st.session_state['project_path'] = typed_path
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

    # ── Recursive AlphaFold-predictions scanner ───────────────────────
    st.divider()
    st.markdown("##### Scan for AlphaFold predictions")
    st.caption("Recursively search for all folders containing AF3 or AF2 prediction data under the current directory.")

    if st.button("Scan for AlphaFold predictions", key="scan_recursive"):
        from core.scanner import find_predictions_recursive
        with st.spinner(f"Scanning `{browse_dir}` ..."):
            found = find_predictions_recursive(str(browse_dir))
        st.session_state['recursive_scan_results'] = found
        st.session_state['recursive_scan_root'] = str(browse_dir)

    # Show results if we have them and we're still in the same browse_dir
    scan_results = st.session_state.get('recursive_scan_results')
    scan_root = st.session_state.get('recursive_scan_root', '')
    if scan_results is not None and scan_root == str(browse_dir):
        if scan_results:
            st.success(f"Found {len(scan_results)} AlphaFold project(s) under `{browse_dir.name}/`")
            for i, proj in enumerate(scan_results):
                col_label, col_fmt, col_btn = st.columns([3, 1, 1])
                with col_label:
                    st.code(proj['label'], language=None)
                with col_fmt:
                    badge = proj.get('format', 'unknown')
                    st.caption(f"📄 {badge}")
                with col_btn:
                    if st.button("Select", key=f"scanpick_{i}", type="primary"):
                        st.session_state['project_path'] = proj['project_path']
                        st.session_state.pop('recursive_scan_results', None)
                        st.rerun()
        else:
            st.warning("No AlphaFold prediction data found under this directory.")

    st.divider()

project_path = st.session_state.get('project_path', "")

if not project_path:
    st.info("Enter a project folder path above, or use Browse to navigate your filesystem. "
            "The folder should contain AF3 or AF2 prediction subdirectories.")
    st.stop()

# Auto-detect where predictions live:
# 1) <project>/AF3/ if it exists and has predictions
# 2) <project>/ itself if it directly contains predictions
from core.scanner import find_predictions_folder
af3_folder = find_predictions_folder(project_path)
if not af3_folder:
    st.error(f"Folder not found: {project_path}")
    st.stop()

# Format override + detection chip
from core.scanner import AF3Scanner
FORMAT_CHOICES = ["Auto-detect", "AF3", "AF2"]
col_fmt, col_chip = st.columns([1, 3])
with col_fmt:
    format_override = st.radio(
        "Format",
        FORMAT_CHOICES,
        horizontal=True,
        index=FORMAT_CHOICES.index(st.session_state.get("format_override", "Auto-detect")),
        help=(
            "Auto-detect inspects each prediction folder. Override if you "
            "want to force one format."
        ),
        key="format_override",
    )
with col_chip:
    # Summarize detected formats across all predictions in af3_folder.
    @st.cache_data(ttl=300)
    def _summarize_formats(folder: str):
        preds = AF3Scanner(Path(folder)).scan()
        counts: dict[str, int] = {}
        pairwise_warn = 0
        for p in preds:
            counts[p.get('format', 'unknown')] = counts.get(p.get('format', 'unknown'), 0) + 1
            if p.get('n_chains', 0) > 2:
                pairwise_warn += 1
        return counts, pairwise_warn, len(preds)

    try:
        fmt_counts, pairwise_warn, n_total = _summarize_formats(str(af3_folder))
    except Exception as e:
        fmt_counts, pairwise_warn, n_total = {}, 0, 0
        st.caption(f"Could not summarize folder: {e}")

    if n_total == 0:
        st.caption("❌ No AF2 or AF3 predictions found in this folder.")
    else:
        label_map = {
            'af3_local': 'AF3',
            'af3_server_flat': 'AF3 Server',
            'af2': 'AF2',
            'unknown': 'unknown',
        }
        parts = [f"{n} {label_map.get(k, k)}" for k, n in sorted(fmt_counts.items(), key=lambda kv: -kv[1])]
        summary = ", ".join(parts)
        if len(fmt_counts) == 1:
            chip = f"✅ Detected {summary} — {n_total} predictions"
        else:
            chip = f"⚠️ Mixed formats: {summary} ({n_total} predictions)"
        if format_override != "Auto-detect":
            chip += f"  ·  🔧 Format locked to **{format_override}** (ignoring auto-detect)"
        if pairwise_warn:
            chip += f"  ·  ⚠️ {pairwise_warn} prediction(s) with 3+ chains — pairwise A/B only"
        st.caption(chip)

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