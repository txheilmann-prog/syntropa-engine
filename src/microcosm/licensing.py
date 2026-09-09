"""
SOURCE LICENSING + REDISTRIBUTION POLICY -- the single machine-readable answer to "may we hand this model out,
and under what terms". Any download or export feature built on this substrate MUST call policy_for() and refuse
(or link out instead of serving bytes) whenever redistribute is False.

Every entry is grounded in the AUTHORITATIVE license text, not recalled: a license is a legal fact and must be
backed by its source. Policy for the sources represented in the open demo library:
  - BiGG / BioModels: each model retains its source publication's terms; BiGG requests citation. Case-by-case;
    default to link + cite, not bulk re-serve. Curated BioModels confirmed CC0 at ingest are freely reusable.
  - Hand-authored components (first-principles / synthetic-lab-table / biomass / biochemistry): ours.

The map is keyed by a PREFIX of a component's source string (the strings carry version suffixes, e.g.
"BioModels-CC0 (BMDB)"), longest-prefix-wins, so a new minor version inherits the policy automatically.
Deployments that ingest additional model sources supply those sources' policies through an optional companion
module, merged in below when present.
"""
from __future__ import annotations

# Licensing is a per-model ANNOTATION + a query-time FILTER, not an ingest gate. Two orthogonal axes drive the
# subset filter:
#   redistribute -- may we hand out the model BYTES ourselves.
#   commercial   -- may the model be used for a COMMERCIAL purpose (False for any NonCommercial license;
#                   None = unknown/depends, treated conservatively).
# `tier` is a short human label for display/filtering. Every entry is grounded in the authoritative license text,
# ordered longest/most-specific first; policy_for takes the first source.startswith(prefix) hit.
_POLICIES = [
    ("BiGG", {
        "license": "BiGG (UCSD): cite the model's own source publication; academic / non-commercial use -- "
                   "commercial use may require a separate UCSD license (bigg.ucsd.edu). Treated conservatively as "
                   "non-commercial.",
        "license_url": "http://bigg.ucsd.edu/",
        "redistribute": False, "commercial": False, "tier": "per-publication (non-commercial)", "cite_required": True,
        "attribution": "BiGG Models (King et al. 2016) + the model's own source publication",
        "cite": "King ZA, et al. Nucleic Acids Res. 2016;44(D1):D515-522.",
        "notes": "Each model keeps its source publication's terms; link + cite, resolve redistribution per model."}),
    ("BioModels-CC0", {  # curated BioModels confirmed CC0 at ingest (ebi.ac.uk terms)
        "license": "CC0 1.0", "license_url": "https://creativecommons.org/publicdomain/zero/1.0/",
        "redistribute": True, "commercial": True, "tier": "CC0", "cite_required": False,
        "attribution": "BioModels (Malik-Sheriff et al. 2020)",
        "cite": "Malik-Sheriff RS, et al. Nucleic Acids Res. 2020;48(D1):D407-D415.",
        "notes": "BioModels hosted models are CC0 (ebi.ac.uk terms); fully reusable."}),
    ("BioModels", {
        "license": "CC0 / per-submission (BioModels)",
        "license_url": "https://www.ebi.ac.uk/biomodels/",
        "redistribute": False, "commercial": None, "tier": "per-submission", "cite_required": True,
        "attribution": "BioModels (Malik-Sheriff et al. 2020)",
        "cite": "Malik-Sheriff RS, et al. Nucleic Acids Res. 2020;48(D1):D407-D415.",
        "notes": "Curated models are CC0 (use the BioModels-CC0 source tag when confirmed); otherwise link + cite."}),
]
try:                                                    # additional data-source policies from a site-local module
    from ._corpus_licensing import CORPUS_POLICIES      # (present in full deployments; not part of the open substrate)
    _POLICIES = CORPUS_POLICIES + _POLICIES
except ImportError:
    pass

# Hand-authored components we built from first principles -- ours to license.
_OURS_PREFIXES = ("first-principles", "synthetic-lab-table", "biomass", "biochemistry")
_OURS = {
    "license": "CC-BY 4.0 (microcosm original)", "license_url": "https://creativecommons.org/licenses/by/4.0/",
    "redistribute": True, "commercial": True, "tier": "CC-BY", "cite_required": True,
    "attribution": "microcosm project (original component)",
    "cite": "microcosm project.", "notes": "Original hand-authored component; ours to license (CC-BY 4.0)."}

_UNKNOWN = {
    "license": "Unknown", "license_url": None, "redistribute": False, "commercial": None, "tier": "unknown",
    "cite_required": True, "attribution": None, "cite": None,
    "notes": "Source not recognized; treated as reference only; DENY redistribution + commercial use until the "
             "license is established (fail-closed). The report tells the user to establish terms with the source."}


def _prefix_hit(source: str, prefix: str, strict: bool) -> bool:
    """Prefix match with a BOUNDARY for permissive (redistribute=True) policies -- closes the fail-open where a
    spoofed/mistyped source like 'BioModels-CC0_fork' would inherit redistribute=True by a bare startswith. A
    permissive prefix matches only on an EXACT hit or a real boundary (space or '(' -- the shape of every real
    source string, e.g. 'BioModels-CC0 (BMDB)'). Deny-by-default policies keep loose startswith: over-matching
    there DENIES more, which is fail-closed and safe on a legal question."""
    if not source.startswith(prefix):
        return False
    if not strict:
        return True
    return len(source) == len(prefix) or source[len(prefix)] in " ("


def policy_for(source: str | None) -> dict:
    """The redistribution policy for a component's source string. Longest/most-specific prefix wins; an
    unrecognized source DENIES redistribution by default (fail-closed, never fail-open on a legal question)."""
    if not source:
        return dict(_UNKNOWN)
    for prefix in _OURS_PREFIXES:
        if _prefix_hit(source, prefix, strict=True):                 # _OURS = redistribute=True -> boundary-strict
            return dict(_OURS)
    for prefix, pol in _POLICIES:
        if _prefix_hit(source, prefix, strict=bool(pol["redistribute"])):
            return dict(pol)
    return dict(_UNKNOWN)


def may_redistribute(source: str | None) -> bool:
    """Fast gate for a download/export path: True only if we may hand out the model file itself."""
    return bool(policy_for(source)["redistribute"])


def may_use_commercially(source: str | None) -> bool:
    """True only if the source's license clearly permits COMMERCIAL use. None (unknown) -> False (fail-closed)."""
    return policy_for(source).get("commercial") is True


# The named licensing SUBSETS a consumer can restrict a model set to: the collection holds everything reachable;
# the consumer picks the slice their intended use allows. Each is a predicate over a source string.
_SUBSETS = {
    "all":            lambda p: True,                                   # everything (no license filter)
    "redistributable": lambda p: p["redistribute"] is True,            # we may re-serve the bytes
    "commercial":     lambda p: p.get("commercial") is True,           # usable commercially (excludes NC + unknown)
    "cc0":            lambda p: p.get("tier") == "CC0",                # public domain only
    "open":           lambda p: p.get("tier") in ("CC0", "CC-BY"),     # permissive, commercial-friendly, re-servable
}


def subsets() -> list[str]:
    """The names of the licensing subsets a consumer may restrict a model set to."""
    return list(_SUBSETS)


def in_subset(source: str | None, subset: str = "all") -> bool:
    """Does `source` belong to the named licensing `subset`? Used to filter a model set at query time so a consumer
    only uses models their intended use permits (e.g. subset='commercial' excludes NonCommercial models)."""
    pred = _SUBSETS.get(subset)
    if pred is None:
        raise ValueError(f"unknown licensing subset {subset!r}; choose from {list(_SUBSETS)}")
    return pred(policy_for(source))


def usage_note(source: str | None) -> str:
    """One-line, report-facing disclosure for a model that was USED in a result: what it is licensed under and,
    when we cannot hand out the bytes, where the consumer must obtain it. This is how a non-redistributable or
    NonCommercial model can PARTICIPATE in a result while the user is told exactly how to license it."""
    p = policy_for(source)
    lic = p.get("license", "Unknown")
    if p.get("redistribute") and p.get("commercial") is True:
        return f"{lic} -- reusable (incl. commercial) with attribution: {p.get('attribution') or source}"
    where = p.get("license_url") or p.get("attribution") or source
    constraint = "non-commercial reuse only" if (p.get("redistribute") and not p.get("commercial")) else \
                 "obtain a license from the source before reuse"
    return f"{lic} -- {constraint}; source: {where}"
