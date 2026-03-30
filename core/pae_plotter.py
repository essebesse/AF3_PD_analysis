"""
Generalized PAE plot generation module for AF3 predictions.
Based on create_pae_plots.py — uses green-white colormap matching the
original interface_analysis scripts.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
from typing import List, Optional, Tuple, Dict


def create_green_white_cmap():
    """Green (low PAE = high confidence) to white (high PAE = low confidence)."""
    colors = ['#006400', '#228B22', '#32CD32', '#90EE90', '#F0FFF0', '#FFFFFF']
    return LinearSegmentedColormap.from_list('green_white', colors, N=100)


def plot_pae_matrix(pae_matrix: np.ndarray,
                    chain_names: List[str],
                    chain_lengths: List[int],
                    title: str = "PAE Matrix",
                    vmax: float = 31.75) -> plt.Figure:
    """
    Plot full PAE matrix with chain boundary lines, labels, and
    red dashed rectangles marking inter-chain (interface) quadrants.

    Args:
        pae_matrix: NxN PAE matrix from confidences.json
        chain_names: List of chain names (e.g., ["IFT27", "RAB23"])
        chain_lengths: List of residue counts per chain
        title: Plot title
        vmax: Maximum PAE value for colorscale

    Returns:
        matplotlib Figure object (caller must do st.pyplot(fig) and plt.close(fig))
    """
    fig, ax = plt.subplots(figsize=(7, 6))
    cmap = create_green_white_cmap()

    im = ax.imshow(pae_matrix, cmap=cmap, vmin=0, vmax=vmax,
                   aspect='auto', origin='lower')

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Predicted Aligned Error (Angstrom)', rotation=270,
                   labelpad=20, fontsize=12)

    # Calculate chain boundaries
    boundaries = [0]
    for length in chain_lengths:
        boundaries.append(boundaries[-1] + length)
    total = boundaries[-1]

    # Black solid chain boundary lines
    for boundary in boundaries[1:-1]:
        ax.axhline(y=boundary - 0.5, color='black', linewidth=2)
        ax.axvline(x=boundary - 0.5, color='black', linewidth=2)

    # Red dashed rectangles on inter-chain quadrants
    if len(chain_lengths) == 2:
        len_a, len_b = chain_lengths[0], chain_lengths[1]
        # Upper-left off-diagonal: chain B rows, chain A cols
        rect1 = patches.Rectangle((0, len_a), len_a, len_b,
                                   linewidth=3, edgecolor='red',
                                   facecolor='none', linestyle='--')
        ax.add_patch(rect1)
        # Lower-right off-diagonal: chain A rows, chain B cols
        rect2 = patches.Rectangle((len_a, 0), len_b, len_a,
                                   linewidth=3, edgecolor='red',
                                   facecolor='none', linestyle='--')
        ax.add_patch(rect2)

    # Residue tick labels — per-chain numbering, avoid overlap at boundaries
    ticks = []
    labels = []
    # Scale min gap with total size to prevent label overlap on x-axis
    min_tick_gap = max(20, total // 30)

    for chain_idx, (name, length) in enumerate(zip(chain_names, chain_lengths)):
        offset = boundaries[chain_idx]
        raw_interval = max(1, length // 4)
        nice_values = [10, 20, 50, 100, 200, 500, 1000]
        interval = nice_values[0]
        for nv in nice_values:
            if nv >= raw_interval:
                interval = nv
                break

        # First residue of chain
        ticks.append(offset)
        labels.append('1')

        # Interior ticks — also check distance to previous tick
        for i in range(interval, length, interval):
            pos = offset + i
            too_close = False
            for b in boundaries[1:-1]:
                if abs(pos - b) < min_tick_gap:
                    too_close = True
                    break
            if not too_close and (not ticks or abs(pos - ticks[-1]) >= min_tick_gap):
                ticks.append(pos)
                labels.append(str(i + 1))

        # Last residue — skip if too close to previous tick or boundary
        last_pos = offset + length - 1
        if not ticks or abs(last_pos - ticks[-1]) >= min_tick_gap:
            too_close = False
            for b in boundaries[1:-1]:
                if abs(last_pos - (b - 1)) < min_tick_gap and last_pos != offset:
                    too_close = True
                    break
            if not too_close:
                ticks.append(last_pos)
                labels.append(str(length))

    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, fontsize=8, rotation=45, ha='right')
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels, fontsize=8)

    # Chain name labels — use xlabel area, no separate "Residue Number" label
    # to avoid overlap. Place chain names with enough room below ticks.
    ax.set_xlabel('')  # no generic label — chain names serve as x-labels
    ax.set_ylabel('Residue Number', fontsize=12, fontweight='bold')

    for i, (name, length) in enumerate(zip(chain_names, chain_lengths)):
        center = boundaries[i] + length / 2
        # Use annotate with axes fraction for consistent placement below plot
        ax.annotate(name, xy=(center, 0), xycoords=('data', 'axes fraction'),
                    xytext=(0, -35), textcoords='offset points',
                    ha='center', fontsize=13, fontweight='bold',
                    annotation_clip=False)

    ax.set_title(f'{title}\nGreen = Low PAE (high confidence), '
                 'White = High PAE (low confidence)',
                 fontsize=13, fontweight='bold', pad=20)

    plt.tight_layout()
    fig.subplots_adjust(bottom=0.15)  # extra room for rotated ticks + chain name labels
    return fig


def plot_pae_interface_zoom(pae_matrix: np.ndarray,
                            chain_a_indices: List[int],
                            chain_b_indices: List[int],
                            chain_a_name: str,
                            chain_b_name: str,
                            pae_cutoff: float = 6.0) -> plt.Figure:
    """
    Plot zoomed-in PAE for the interface region (chain A vs chain B block).
    Red dots mark contacts below pae_cutoff.

    Args:
        pae_matrix: NxN PAE matrix from confidences.json
        chain_a_indices: List of residue indices for chain A
        chain_b_indices: List of residue indices for chain B
        chain_a_name: Name of chain A
        chain_b_name: Name of chain B
        pae_cutoff: PAE cutoff for highlighting contacts with red dots

    Returns:
        matplotlib Figure object
    """
    fig, ax = plt.subplots(figsize=(7, 6))
    cmap = create_green_white_cmap()

    # Extract the off-diagonal block: chain A rows, chain B columns
    a_start, a_end = chain_a_indices[0], chain_a_indices[-1] + 1
    b_start, b_end = chain_b_indices[0], chain_b_indices[-1] + 1
    interface_block = pae_matrix[a_start:a_end, b_start:b_end]

    im = ax.imshow(interface_block, cmap=cmap, vmin=0, vmax=31.75,
                   aspect='auto', origin='lower')

    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Predicted Aligned Error (Angstrom)', rotation=270,
                   labelpad=20, fontsize=12)


    # Tick labels — scale interval with protein length to avoid overlap
    len_a = a_end - a_start
    len_b = b_end - b_start

    def _tick_interval(length):
        if length > 800:
            return 200
        elif length > 400:
            return 100
        elif length > 200:
            return 50
        else:
            return 20

    x_interval = _tick_interval(len_b)
    xticks = list(range(0, len_b + 1, x_interval))
    ax.set_xticks(xticks)
    ax.set_xticklabels([str(t + 1) for t in xticks], fontsize=8)

    y_interval = _tick_interval(len_a)
    yticks = list(range(0, len_a + 1, y_interval))
    ax.set_yticks(yticks)
    ax.set_yticklabels([str(t + 1) for t in yticks], fontsize=8)

    # Mark residue pairs with PAE < 6 and < 3 using stars
    rows_6, cols_6 = np.where((interface_block < pae_cutoff) & (interface_block >= 3.0))
    rows_3, cols_3 = np.where(interface_block < 3.0)

    if len(rows_6) > 0:
        ax.scatter(cols_6, rows_6, marker='*', s=15, c='orange', alpha=0.8,
                   edgecolors='none', label=f'PAE < {pae_cutoff:.0f}Å ({len(rows_6)})')
    if len(rows_3) > 0:
        ax.scatter(cols_3, rows_3, marker='*', s=25, c='red', alpha=0.9,
                   edgecolors='none', label=f'PAE < 3Å ({len(rows_3)})')

    if len(rows_6) > 0 or len(rows_3) > 0:
        ax.legend(loc='upper right', fontsize=9, framealpha=0.9)

    ax.set_xlabel(f'{chain_b_name} Residue Number', fontsize=12, fontweight='bold')
    ax.set_ylabel(f'{chain_a_name} Residue Number', fontsize=12, fontweight='bold')
    ax.set_title(f'Interface PAE: {chain_a_name} vs {chain_b_name}',
                 fontsize=13, fontweight='bold', pad=15)

    plt.tight_layout()
    return fig


def plot_all_models_comparison(pae_matrices: List[np.ndarray],
                               model_labels: List[str],
                               chain_lengths: List[int],
                               chain_names: List[str]) -> plt.Figure:
    """
    Side-by-side PAE plots for all models of a prediction.

    Args:
        pae_matrices: List of NxN PAE matrices (one per model)
        model_labels: List of model labels (e.g., ["Top", "s1-m0", ...])
        chain_lengths: List of residue counts per chain
        chain_names: List of chain names

    Returns:
        matplotlib Figure object with subplots
    """
    n_models = len(pae_matrices)
    if n_models == 0:
        raise ValueError("At least one PAE matrix required")

    # Layout: up to 3 per row
    ncols = min(n_models, 3)
    nrows = (n_models + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.5 * nrows))
    if n_models == 1:
        axes = [axes]
    else:
        axes = axes.flatten()

    cmap = create_green_white_cmap()
    vmax = 31.75

    # Calculate boundaries once
    boundaries = [0]
    for length in chain_lengths:
        boundaries.append(boundaries[-1] + length)

    im = None
    for i in range(n_models):
        ax = axes[i]
        pae_matrix = pae_matrices[i]
        label = model_labels[i]

        im = ax.imshow(pae_matrix, cmap=cmap, vmin=0, vmax=vmax,
                       aspect='auto', origin='lower')

        # Black chain boundary lines
        for boundary in boundaries[1:-1]:
            ax.axhline(y=boundary - 0.5, color='black', linewidth=1.5)
            ax.axvline(x=boundary - 0.5, color='black', linewidth=1.5)

        # Red dashed interface rectangles for 2-chain
        if len(chain_lengths) == 2:
            len_a, len_b = chain_lengths[0], chain_lengths[1]
            rect = patches.Rectangle((len_a, 0), len_b, len_a,
                                      linewidth=2, edgecolor='red',
                                      facecolor='none', linestyle='--')
            ax.add_patch(rect)

        # Contact counts for subtitle
        if len(chain_lengths) == 2:
            pae_inter = pae_matrix[:chain_lengths[0], chain_lengths[0]:]
            n6 = int(np.sum(pae_inter < 6.0))
            n3 = int(np.sum(pae_inter < 3.0))
            ax.set_title(f'{label}\n<3: {n3}, <6: {n6}',
                         fontsize=11, fontweight='bold')
        else:
            ax.set_title(label, fontsize=11, fontweight='bold')

        # Minimal ticks
        tick_positions = [0, boundaries[-1] - 1]
        tick_labels_list = ['1', str(boundaries[-1])]
        for b in boundaries[1:-1]:
            tick_positions.extend([b - 1, b])
            tick_labels_list.extend([str(b), '1'])
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels_list, fontsize=8)
        ax.set_yticks(tick_positions)
        ax.set_yticklabels(tick_labels_list, fontsize=8)

        # Chain name labels
        for j, (name, length) in enumerate(zip(chain_names, chain_lengths)):
            center = boundaries[j] + length / 2
            ax.text(center, -boundaries[-1] * 0.06, name, ha='center',
                    fontsize=8, fontweight='bold')

    # Hide unused subplots
    for i in range(n_models, len(axes)):
        fig.delaxes(axes[i])

    # Shared colorbar
    if im is not None:
        cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
        cbar = fig.colorbar(im, cax=cbar_ax)
        cbar.set_label('PAE (Angstrom)', rotation=270, labelpad=20, fontsize=12)

    plt.tight_layout(rect=[0, 0, 0.9, 0.96])
    return fig


def get_chain_info_from_confidences(confidences: Dict) -> Tuple[List[str], List[int]]:
    """
    Extract chain names and lengths from confidences.json.

    Args:
        confidences: Dictionary loaded from confidences.json

    Returns:
        Tuple of (chain_names, chain_lengths)
    """
    token_chain_ids = confidences.get('token_chain_ids', [])

    # Preserve order of first appearance
    unique_chains = []
    for c in token_chain_ids:
        if c not in unique_chains:
            unique_chains.append(c)

    # Calculate chain lengths
    chain_lengths = []
    for chain in unique_chains:
        count = sum(1 for c in token_chain_ids if c == chain)
        chain_lengths.append(count)

    # Use chain IDs as names (can be overridden later with UniProt gene names)
    chain_names = unique_chains

    return chain_names, chain_lengths


def get_chain_indices_from_confidences(confidences: Dict) -> Dict[str, List[int]]:
    """
    Get residue indices for each chain from confidences.json.

    Args:
        confidences: Dictionary loaded from confidences.json

    Returns:
        Dictionary mapping chain name to list of residue indices
    """
    token_chain_ids = confidences.get('token_chain_ids', [])

    chain_indices = {}
    for i, chain in enumerate(token_chain_ids):
        if chain not in chain_indices:
            chain_indices[chain] = []
        chain_indices[chain].append(i)

    return chain_indices
