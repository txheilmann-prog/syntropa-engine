"""
Component schema (v0.1) -- the black box. DELIBERATELY MINIMAL: we build a few real, diverse
components against this, run them, and let what breaks drive the schema before we freeze v1. Do not
add fields speculatively.

A Component is a declarative manifest (JSON, versioned in git) describing:
  - identity + provenance (MULTI-SOURCE: >=1 source, cross-confirmation is the librarian's job)
  - PORTS: its interface to the shared medium, each a species with a GROUNDED chemical identity
    (formula) and a role. Identity is how boxes free-associate: two boxes that name the same medium
    species MUST agree on its formula, or it is an error (the legibility guard, enforced at wiring).
  - a TRANSFER spec (how the box computes outputs from inputs): function | lookup | sbml (v0.1).
  - a VALIDITY domain: conditions under which the box is meaningful. Composing/driving a box outside
    its domain is a first-class warning -- honest composition, not plausible nonsense.

The schema is versioned (SCHEMA_VERSION); every manifest carries schema_version; migrations are
explicit. Extra keys are preserved (forward-compatible) but validated keys are strict.
"""
from __future__ import annotations

import json
import os
import re
import sys
from dataclasses import dataclass, field

# Units + dimensions are owned by units.py (the SSoT for physical-quantity discipline). Re-exported here
# so existing callers keep importing them from component.
from .units import (ATOMIC_WEIGHTS, CANONICAL_UNIT, canonical_formula, dimension_of,  # noqa: F401
                    from_canonical, molar_mass, parse_formula, to_canonical, validate_unit)

SCHEMA_VERSION = "1.0"          # frozen at v1; v2 extends additively via migrations
KINDS = {"chemistry", "reaction_network", "organism", "lookup", "logical"}
TRANSFER_TYPES = {"function", "lookup", "sbml", "fba"}
ROLES = {"input", "output", "catalyst", "modifier", "environment"}


@dataclass(frozen=True)
class Port:
    species: str                 # medium species id (the free-association key)
    formula: str | None          # grounded chemical identity; None = a grounding GAP (chemical, not yet grounded)
    role: str                    # one of ROLES
    chemical: bool = True        # False = protein/pool/abstract state: conserved by identity, NOT mass-checkable.
                                 # This DISAMBIGUATES "legitimately not a chemical" from "chemical but ungrounded".
    entity_id: str | None = None  # canonical identity: InChIKey/ChEBI (chemical) or UniProt/name (non-chemical).
                                 # IDENTITY, not atoms: it distinguishes HCO3- from CO3-- where formula/charge cannot.
    unit: str = CANONICAL_UNIT   # the unit THIS box speaks for the species; the engine reconciles to canonical
    charge: float | None = None  # net charge (signed) of the species AS TRACKED, for charge-balance checking.
                                 # Formula carries the ACTUAL atoms of that protonation state (HCO3- = C H O3).

    def atoms(self) -> dict[str, float]:
        # only chemical species with a grounded formula contribute atoms; proteins/pools contribute none
        return parse_formula(self.formula) if (self.chemical and self.formula) else {}


@dataclass
class Component:
    id: str
    name: str
    kind: str
    provenance: list[dict]
    transfer: dict               # {"type": ..., ...spec...}
    ports: dict[str, Port]       # species id -> Port
    params: dict = field(default_factory=dict)
    validity: dict = field(default_factory=dict)
    notes: str = ""
    schema_version: str = SCHEMA_VERSION
    source_path: str | None = None
    payload: bytes | None = None      # SBML/table content, stored in the DB (not a loose file)
    open_boundary: bool = False       # True = exchanges matter with an unmodelled compartment (air, feed,
                                      # biomass); its net atom-flux is a DECLARED boundary, not a leak.
    context: dict = field(default_factory=dict)  # v2 additive: the environment the box ASSUMES, for the
                                      # composition-validity gate. {phase, requires:[species], bridges:[phases],
                                      # ranges:{intensive_var:[lo,hi]}}. Empty = unconstrained (v1-compatible).

    def inputs(self) -> list[Port]:
        return [p for p in self.ports.values() if p.role in ("input", "catalyst", "modifier")]

    def outputs(self) -> list[Port]:
        return [p for p in self.ports.values() if p.role == "output"]

    def validity_warnings(self, state: dict[str, float]) -> list[str]:
        """Honest composition: driving a box outside its declared domain is a
        first-class WARNING, not a silent plausible-nonsense result. A validity entry {species: [lo, hi]}
        bounds a medium condition IN THAT PORT'S UNIT (e.g. temperature [0,40] in degC), so the canonical
        medium state (K) is converted to the port unit before the range check."""
        warns = []
        for key, dom in self.validity.items():
            if isinstance(dom, list) and len(dom) == 2 and key in state:
                p = self.ports.get(key)
                v = from_canonical(state[key], p.unit, p.formula) if p else state[key]
                if v < dom[0] or v > dom[1]:
                    unit = f" {p.unit}" if p and p.unit else ""
                    warns.append(f"{self.id}: {key}={v:g}{unit} outside validity domain {dom}")
        return warns

    def validate(self) -> list[str]:
        """Return a list of problems; empty == valid. Loud, mechanical, no silent acceptance."""
        errs = []
        if self.schema_version != SCHEMA_VERSION:
            errs.append(f"schema_version {self.schema_version!r} != {SCHEMA_VERSION!r} (needs migration)")
        if not self.id or not re.match(r"^[a-z0-9_]+$", self.id):
            errs.append(f"id {self.id!r} must be lower_snake_case")
        if self.kind not in KINDS:
            errs.append(f"kind {self.kind!r} not in {sorted(KINDS)}")
        if not self.provenance:
            errs.append("provenance is empty (every component needs >=1 source)")
        t = self.transfer.get("type")
        if t not in TRANSFER_TYPES:
            errs.append(f"transfer.type {t!r} not in {sorted(TRANSFER_TYPES)}")
        if not self.ports:
            errs.append("no ports declared")
        for sid, p in self.ports.items():
            if p.species != sid:
                errs.append(f"port key {sid!r} != port.species {p.species!r}")
            if p.role not in ROLES:
                errs.append(f"port {sid!r} role {p.role!r} not in {sorted(ROLES)}")
            if not p.chemical and p.formula:
                errs.append(f"port {sid!r} is non-chemical but carries a formula {p.formula!r} "
                            f"(contradiction: a non-chemical species cannot be mass-checked)")
            unit_err = validate_unit(p.unit, chemical=p.chemical, formula=p.formula)   # units discipline
            if unit_err:
                errs.append(f"port {sid!r}: {unit_err}")
        return errs


def resolve_model_path(ref: str, anchor_path: str) -> str:
    """Portably resolve a stored payload/model path. Manifests may carry ABSOLUTE paths from the machine that
    created them (e.g. a Windows `C:\\...\\library\\...` path); on another OS those don't exist. Use
    the path as-is if it exists, else graft its `library/...` tail onto the library root of the anchor (a locally
    correct path -- the manifest's own source_path). A no-op where the path already resolves; makes models OS-portable."""
    if ref and os.path.exists(ref):
        return ref
    q = (ref or "").replace("\\", "/")
    a = (anchor_path or "").replace("\\", "/")
    j = a.rfind("/library/")
    i = q.rfind("library/")
    if j >= 0 and i >= 0:
        cand = os.path.join(a[:j], q[i:])
        if os.path.exists(cand):
            return cand
    if ref and os.path.isabs(ref):
        return ref
    return os.path.join(os.path.dirname(anchor_path or ""), ref or "")


def _fetch_bigg_model(provenance, dest: str) -> None:
    """A demo model whose upstream license does not permit redistribution ships LINK-ONLY (see MODELS.md): the
    manifest is present, the model file is not. Fetch it from BiGG -- the source named in the manifest's
    provenance -- into `dest`, so the bundled demos run from a clean checkout. Fetched once, then reused.

    Hardened against a hostile manifest: the BiGG id is charset-validated (no injected path/URL segments), the
    download follows only same-host redirects, and the response is size-capped."""
    import re
    import time
    import urllib.error
    import urllib.request
    from urllib.parse import urlparse

    bigg_id = None
    for e in provenance or []:
        if (e.get("source") or "").lower() == "bigg":
            parts = (e.get("ref") or "").split("(")[0].split()
            if parts:
                bigg_id = parts[0]
                break
    if not bigg_id:
        raise FileNotFoundError(f"{dest} is not bundled and its manifest names no fetchable BiGG model; see MODELS.md.")
    if not re.fullmatch(r"[A-Za-z0-9_]+", bigg_id):        # never let a manifest inject path or URL segments
        raise ValueError(f"refusing to fetch: {bigg_id!r} is not a valid BiGG model id.")

    allowed_hosts = {"bigg.ucsd.edu", "bigg.bio"}          # BiGG is migrating ucsd.edu -> bigg.bio; allow both

    class _SameHost(urllib.request.HTTPRedirectHandler):   # refuse a redirect off the BiGG hosts (SSRF guard)
        def redirect_request(self, req, fp, code, msg, headers, newurl):
            if urlparse(newurl).hostname not in allowed_hosts:
                raise urllib.error.HTTPError(newurl, code, "cross-host redirect refused", headers, fp)
            return super().redirect_request(req, fp, code, msg, headers, newurl)

    url = f"https://bigg.ucsd.edu/static/models/{bigg_id}.xml"
    max_bytes = 200 * 1024 * 1024
    opener = urllib.request.build_opener(_SameHost())
    os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)
    sys.stderr.write(f"[microcosm] fetching link-only model {bigg_id} from BiGG ...\n")
    last = None
    for attempt in range(3):                       # a transient DNS/network hiccup should not fail the demo
        try:
            with opener.open(url, timeout=90) as resp:
                if urlparse(resp.url).hostname not in allowed_hosts:
                    raise ValueError(f"refusing content from unexpected host: {resp.url!r}")
                data = resp.read(max_bytes + 1)
            if len(data) > max_bytes:
                raise ValueError(f"refusing {bigg_id}: response exceeds the {max_bytes}-byte cap.")
            with open(dest, "wb") as fh:
                fh.write(data)
            return
        except Exception as ex:
            last = ex
            time.sleep(2 * (attempt + 1))
    raise FileNotFoundError(
        f"{dest} is not bundled (link-only model; see MODELS.md) and the fetch from {url} failed: {last}. "
        f"Obtain the model from BiGG and place it at that path.") from last


def _library_root(anchor):
    """The `library` directory at or above the manifest anchor -- the safe root for fetched model files. Falls
    back to the anchor's own directory when there is no `library` ancestor."""
    a = os.path.realpath(os.path.dirname(anchor or "."))
    parts = a.replace("\\", "/").split("/")
    if "library" in parts:
        idx = len(parts) - 1 - parts[::-1].index("library")           # last `library` segment in the path
        return os.path.realpath("/".join(parts[:idx + 1]) or "/")
    return a


def _within(path, root):
    """True iff `path` resolves inside `root` (mismatched drives / relative-vs-absolute -> treated as outside)."""
    try:
        return os.path.commonpath([os.path.realpath(path), root]) == root
    except ValueError:
        return False


def ensure_model_file(comp) -> str:
    """Resolve a component's model-payload path, fetching a link-only demo model from its cited source if the
    file is not present (see MODELS.md). One place both the byte loader and the FBA path resolver go through."""
    ref = comp.transfer.get("path")
    rp = resolve_model_path(ref, getattr(comp, "source_path", None))
    if ref and not os.path.exists(rp):
        root = _library_root(getattr(comp, "source_path", None))      # untrusted-manifest guard: only ever fetch
        if root and not _within(rp, root):                            # into the library tree, never via traversal
            raise ValueError(f"refusing to fetch a link-only model outside the library tree "
                             f"(resolved path {rp!r} escapes {root!r}); see MODELS.md.")
        _fetch_bigg_model(getattr(comp, "provenance", None), rp)
    return rp


def load_component(path: str) -> Component:
    """Load + validate a manifest. Raises ValueError with ALL problems (fail loud, fail complete)."""
    with open(path, encoding="utf-8") as f:
        m = json.load(f)
    ports = {sid: Port(species=sid, formula=canonical_formula(pv.get("formula")), role=pv.get("role", "input"),
                       chemical=pv.get("chemical", True), entity_id=pv.get("entity_id"),
                       unit=pv.get("unit", CANONICAL_UNIT), charge=pv.get("charge"))
             for sid, pv in (m.get("ports") or {}).items()}
    # (resolve_model_path defined below is used here + by community._model_path)
    c = Component(id=m.get("id", ""), name=m.get("name", ""), kind=m.get("kind", ""),
                  provenance=m.get("provenance", []), transfer=m.get("transfer", {}),
                  ports=ports, params=m.get("params", {}), validity=m.get("validity", {}),
                  notes=m.get("notes", ""), schema_version=m.get("schema_version", "?"),
                  open_boundary=bool(m.get("open_boundary", False)), context=m.get("context", {}),
                  source_path=os.path.abspath(path))
    # payload (SBML/table) travels WITH the component so the DB, not the filesystem, is the store
    ref = c.transfer.get("path")
    if ref:
        rp = ensure_model_file(c)
        with open(rp, "rb") as pf:
            c.payload = pf.read()
    problems = c.validate()
    if problems:
        raise ValueError(f"invalid component {path}:\n  - " + "\n  - ".join(problems))
    return c
