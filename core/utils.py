"""
Shared utility functions for the AF3 Analysis App.
"""

import re
from typing import Optional, Tuple
import requests


# UniProt accession format (covers both 6- and 10-character canonical forms).
# See https://www.uniprot.org/help/accession_numbers
_UNIPROT_ACC_RE = re.compile(
    r'^[OPQ][0-9][A-Z0-9]{3}[0-9]$'
    r'|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$'
)


def is_valid_uniprot_accession(acc: str) -> bool:
    """True if the string looks like a canonical UniProtKB accession."""
    if not acc:
        return False
    return bool(_UNIPROT_ACC_RE.match(acc.strip().upper()))


def normalize_cache_records(records: list) -> list:
    """Backfill missing iPTM/PTM/ranking_score on legacy cache records.

    Some pre-existing caches (``format: 'af2'`` from earlier analyzers) store
    only the combined ``iptm_ptm`` and leave ``iptm`` / ``ptm`` / ``ranking_score``
    as ``None``. The display code expects numeric values everywhere; this
    populates the standard fields in place from whatever combined score is
    available (or 0 as last resort). Returns the same list for chaining.
    """
    for r in records:
        combined = r.get('iptm_ptm')
        if r.get('iptm') is None:
            r['iptm'] = float(combined) if combined is not None else 0.0
        if r.get('ptm') is None:
            r['ptm'] = float(combined) if combined is not None else 0.0
        if r.get('ranking_score') is None:
            r['ranking_score'] = float(combined) if combined is not None else float(r.get('iptm') or 0.0)
    return records


def split_prediction_name(name: str) -> Tuple[str, str]:
    """Split an AF3 prediction folder name into ``(bait_acc, prey_acc)``.

    Strips the AF3 Server ``fold_`` prefix and uppercases both sides. If the
    name has no ``_and_`` separator, returns ``(name.upper(), '')``.
    """
    if not name:
        return '', ''
    s = name
    if s.startswith('fold_'):
        s = s[5:]
    if '_and_' in s:
        bait, prey = s.split('_and_', 1)
        return bait.upper(), prey.upper()
    return s.upper(), ''


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
    iptm = iptm or 0.0
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
        parts = name.split("_and_", 1)
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


def _parse_uniprot_tsv(text: str, result: dict) -> None:
    """Parse a TSV response body from UniProt and merge into ``result``."""
    for line in text.strip().split('\n')[1:]:  # skip header
        parts = line.split('\t')
        if len(parts) >= 2:
            acc = parts[0].strip().upper()
            gene = parts[1].strip().split()[0] if parts[1].strip() else ''
            if acc and gene:
                result[acc] = gene


def fetch_gene_names_batch(accessions: list, progress_callback=None) -> dict:
    """
    Fetch gene names for many accessions using UniProt batch search API.

    Processes valid-shape accessions in batches of 100. UniProt rejects the
    whole batch with HTTP 400 if *any* term has invalid syntax (e.g. an AF3
    Server ``fold_xxx`` ID), so we pre-filter to canonical-form accessions
    before querying; on 400 we still fall back to per-accession requests.

    Calls ``progress_callback(done, total)`` after each batch.

    Returns dict mapping ``ACCESSION -> gene name`` (empty string if not
    found). Every input accession appears as a key.
    """
    import time

    # Normalise + de-duplicate while preserving original keys for the return map
    normalised = {acc.strip().upper() for acc in accessions if acc}
    queryable = sorted(a for a in normalised if is_valid_uniprot_accession(a))

    result = {}
    total = len(queryable)
    done = 0
    batch_size = 100

    def _do_query(query_accs):
        """Return (parsed_result_dict, status_code) for one batch."""
        if not query_accs:
            return {}, 200
        query = " OR ".join(f"accession:{a}" for a in query_accs)
        url = (
            "https://rest.uniprot.org/uniprotkb/search"
            f"?query={requests.utils.quote(query)}"
            "&fields=accession,gene_names&format=tsv&size=500"
        )
        try:
            resp = requests.get(url, timeout=15)
        except requests.RequestException:
            return {}, None
        if resp.status_code == 429:
            time.sleep(2)
            try:
                resp = requests.get(url, timeout=15)
            except requests.RequestException:
                return {}, None
        batch_result = {}
        if resp.status_code == 200:
            _parse_uniprot_tsv(resp.text, batch_result)
        return batch_result, resp.status_code

    for i in range(0, total, batch_size):
        batch = queryable[i:i + batch_size]
        batch_result, status = _do_query(batch)

        # If the batch itself was rejected (typically 400 because one term has
        # bad syntax we didn't catch), fall back to per-accession queries so
        # one bad ID doesn't take down 99 good ones.
        if status == 400:
            for single in batch:
                single_result, _ = _do_query([single])
                batch_result.update(single_result)
                time.sleep(0.1)

        result.update(batch_result)

        done += len(batch)
        if progress_callback:
            progress_callback(done, total)

        # Rate limit: ~1 request/second to respect UniProt API limits
        if i + batch_size < total:
            time.sleep(1.0)

    # Fill in every input accession (incl. non-UniProt-shaped ones) with
    # the empty string when no gene name was found, so callers can
    # distinguish "looked up" from "not yet looked up".
    for acc in normalised:
        result.setdefault(acc, '')

    return result


