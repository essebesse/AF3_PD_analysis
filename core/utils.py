"""
Shared utility functions for the AF3 Analysis App.
"""

from typing import Optional
import requests


def format_score(value: float, decimals: int = 3) -> str:
    """Format a score value for display."""
    if value is None:
        return "-"
    return f"{value:.{decimals}f}"


def calculate_confidence_tier(iptm: float, ipsae: Optional[float] = None) -> str:
    """
    Calculate confidence tier based on iPTM and optionally ipSAE.

    Args:
        iptm: iPTM score
        ipsae: Optional ipSAE score for more refined classification

    Returns:
        Confidence tier string: "High", "Medium", "Low", or "Very Low"
    """
    # Use ipSAE if available, otherwise fall back to iPTM
    if ipsae is not None:
        # ipSAE-based classification (Dunbrack 2025)
        if ipsae > 0.7:
            return "High"
        elif ipsae >= 0.5:
            return "Medium"
        elif ipsae >= 0.3:
            return "Low"
        else:
            return "Very Low"
    else:
        # iPTM-based classification
        if iptm > 0.8:
            return "High"
        elif iptm > 0.6:
            return "Medium"
        elif iptm > 0.4:
            return "Low"
        else:
            return "Very Low"


def tier_color(tier: str) -> str:
    """Get color code for confidence tier."""
    colors = {
        "High": "#22c55e",      # Green
        "Medium": "#eab308",    # Yellow
        "Low": "#f97316",       # Orange
        "Very Low": "#ef4444",  # Red
    }
    return colors.get(tier, "#6b7280")  # Gray for unknown


def parse_prediction_name(name: str) -> tuple:
    """
    Parse prediction name into bait and prey components.

    Args:
        name: Prediction name (e.g., "q9bw83_and_p12345")

    Returns:
        Tuple of (bait, prey) strings
    """
    if "_and_" in name:
        parts = name.split("_and_")
        return parts[0], parts[1]
    return name, ""


def format_model_label(seed: int, sample: int, is_top_ranked: bool = False) -> str:
    """
    Format seed/sample as a display label.

    Args:
        seed: Seed number
        sample: Sample number
        is_top_ranked: Whether this is the top-ranked model

    Returns:
        Formatted label string
    """
    if is_top_ranked:
        return "Top"
    return f"s{seed}-m{sample}"


def fetch_gene_names_batch(accessions: list, progress_callback=None) -> dict:
    """
    Fetch gene names for many accessions using UniProt batch search API.
    Processes in batches of 100; calls progress_callback(done, total) after each batch.

    Returns dict mapping ACCESSION -> gene name (falls back to accession if not found).
    """
    result = {}
    total = len(accessions)
    done = 0
    batch_size = 100

    for i in range(0, total, batch_size):
        batch = accessions[i:i + batch_size]
        query = " OR ".join(f"accession:{acc}" for acc in batch)
        url = (
            "https://rest.uniprot.org/uniprotkb/search"
            f"?query={requests.utils.quote(query)}"
            "&fields=accession,gene_names&format=tsv&size=500"
        )
        try:
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                lines = response.text.strip().split('\n')
                for line in lines[1:]:  # skip header
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        acc = parts[0].strip().upper()
                        gene = parts[1].strip().split()[0] if parts[1].strip() else ''
                        if acc and gene:
                            result[acc] = gene
        except Exception:
            pass

        done += len(batch)
        if progress_callback:
            progress_callback(done, total)

    # Fill in accessions not returned by API (e.g. obsolete IDs) with empty string
    # so display code can distinguish "no gene name found" from "has gene name"
    for acc in accessions:
        if acc.upper() not in result:
            result[acc.upper()] = ''

    return result


