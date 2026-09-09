"""
Engine: stitch components into a netlist over a SHARED GROUNDED MEDIUM and drive it forward.

Netlist (declarative, versioned): which components are present, the initial medium, environmental
drivers, duration/dt. Boxes FREE-ASSOCIATE through the medium: each binds its ports to medium species
by id. GUARD (enforced, not assumed): if two components declare the same medium species with
DIFFERENT chemical formulas, that is an identity conflict and the netlist is rejected -- you cannot
free-associate species you disagree about. This is the legibility guarantee at the wiring boundary.

Runs on Vivarium (each box = a Process); the ConservationMonitor is applied to the whole-medium
trajectory. The monitor is VENDORED (microcosm.conservation) so the package stands alone.
"""
from __future__ import annotations

import math                                               # noqa: E402
try:                                                      # vivarium drives the netlist SIMULATION only; the
    from vivarium.core.process import Process             # noqa: E402  FBA and validity paths do not need
    from vivarium.core.engine import Engine               # noqa: E402  it, so it stays optional here -- a
    _HAS_VIVARIUM = True                                  # caller that does not simulate still imports
except ModuleNotFoundError:                               # cleanly.
    _HAS_VIVARIUM = False

    class Process:                                        # stub base so BlackBoxProcess parses; run_netlist guards
        pass
    Engine = None
from .conservation import ConservationMonitor             # noqa: E402
from .transfer import build_transfer                      # noqa: E402
from .units import (NONACCUMULABLE_DIMS, UNITS, amount_to_concentration,  # noqa: E402
                    concentration_to_amount, dimension_of, from_canonical, parse_formula, to_canonical)


def unit_reconciled(transfer, ports, volume: float = 1.0):
    """Wrap a box's transfer so it speaks its OWN declared units while the medium stays canonical (amount
    in mmol; conditions in K, bar). Reading in: convert each medium value (ABSOLUTE) to the box's unit --
    an amount unit converts within the amount dimension; a CONCENTRATION unit is bridged amount->conc via
    volume, then to the box unit; a condition unit converts within its dimension. Writing out: reverse,
    as DELTAS (affine offsets are not applied to a rate), and GUARD -- a box may READ an intensive
    condition (temperature/pressure/pH/potential) but may not emit an accumulating delta for one (no
    combine rule yet). Boxes already entirely in canonical mmol are returned unwrapped (no overhead)."""
    plan = {s: (dimension_of(p.unit), p.unit, p.formula) for s, p in ports.items()}
    trivial = all(dim == "amount" and UNITS.get(u, (None, ("", 0)))[1] == ("linear", 1.0)
                  for dim, u, f in plan.values())
    if trivial:
        return transfer

    def to_box(s, medium_val):
        dim, u, f = plan[s]
        if dim == "concentration":
            return from_canonical(amount_to_concentration(medium_val, volume), u, f)
        return from_canonical(medium_val, u, f)                # amount / condition: within-dimension

    def to_medium(s, box_delta):
        dim, u, f = plan[s]
        if dim in NONACCUMULABLE_DIMS:
            raise ValueError(f"box emitted a delta for '{s}' ({dim}): intensive conditions cannot be "
                             f"accumulated (a box may read them, not write an additive change) in v1")
        if dim == "concentration":
            return concentration_to_amount(to_canonical(box_delta, u, f, delta=True), volume)
        return to_canonical(box_delta, u, f, delta=True)

    def wrapped(avail, dt):
        box_avail = dict(avail)
        for s in plan:
            if s in avail:
                box_avail[s] = to_box(s, avail[s])
        deltas = transfer(box_avail, dt)
        return {s: (to_medium(s, v) if s in plan else v) for s, v in deltas.items()}
    return wrapped


KNOWN_PHASES = {"aqueous", "gas", "membrane", "surface", "solid"}


def composition_validity(components) -> list[str]:
    """Is this netlist even POSSIBLE as one physical system? A validity axis PRIOR to conservation and
    CHEAPER (static, no simulation): conservation assumes the boxes can coexist; often they cannot.
    Checks the boxes' DECLARED environmental preconditions (component.context) for mutual satisfiability
    -- design-by-contract / assume-guarantee + the Gibbs phase rule. Returns a list of problems ([] = a
    valid, coexistable composition). Components without a declared context are unconstrained (v1-compatible).

      1. PHASE consistency: all declared phases equal, OR an explicit interphase-BRIDGE component connects
         them (the aeration box is the gas<->aqueous template). Unbridged multi-phase -> impossible.
      2. INTENSIVE-VARIABLE range intersection: intersect each box's required range for a shared variable
         (pH, temperature, ...); an EMPTY intersection means no single medium state satisfies all boxes.
      3. IMPLICATION closure: a required constituent must be present in the medium; a charged species
         implies an aqueous phase.
    """
    problems: list[str] = []
    all_species = {s for c in components for s in c.ports}

    phases, bridged = {}, set()
    for c in components:
        ph = c.context.get("phase")
        if ph:
            if ph not in KNOWN_PHASES:
                problems.append(f"{c.id}: unknown phase {ph!r} (not in {sorted(KNOWN_PHASES)})")
            phases[c.id] = ph
        bridged.update(c.context.get("bridges", []))
    distinct = set(phases.values())
    unbridged = distinct - bridged
    if len(unbridged) > 1:
        problems.append(f"phase conflict: components require phases {sorted(distinct)} but no interphase "
                        f"bridge connects {sorted(unbridged)} -- add a bridge component (e.g. a gas<->aqueous "
                        f"transfer like aeration), or they cannot share a medium")

    ranges: dict[str, list] = {}
    for c in components:
        for var, rng in c.context.get("ranges", {}).items():
            if isinstance(rng, list) and len(rng) == 2:
                ranges.setdefault(var, []).append((c.id, float(rng[0]), float(rng[1])))
    for var, lst in ranges.items():
        lo, hi = max(x[1] for x in lst), min(x[2] for x in lst)
        if lo > hi:
            who = ", ".join(f"{i}[{l:g},{h:g}]" for i, l, h in lst)
            problems.append(f"no common {var}: required ranges {who} do not overlap -- the ecosystem cannot "
                            f"satisfy every box at one {var}")

    for c in components:
        for req in c.context.get("requires", []):
            if req not in all_species:
                problems.append(f"{c.id} requires '{req}' present in the medium, but no component provides it")

    charged = sorted({s for c in components for s, p in c.ports.items() if p.chemical and p.charge})
    if charged and distinct and "aqueous" not in distinct and "aqueous" not in bridged:
        problems.append(f"charged species {charged} imply an aqueous phase, but the declared phases are "
                        f"{sorted(distinct)} (ions require solution)")
    return problems


def medium_grounding(components) -> dict[str, str]:
    """Union of all ports' species -> formula, REJECTING identity conflicts (free-association guard).
    Identity is layered: a canonical `entity_id` (ChEBI/InChIKey) is the true
    identity -- two boxes that give the same medium species DIFFERENT entity_ids mean different things and
    cannot free-associate. The `formula` is ATOMS for conservation (same pool -> same atoms), and
    chemical-ness is global. Only chemical, grounded species are returned (proteins are conserved by id)."""
    grounding: dict[str, str | None] = {}
    chemical: dict[str, bool] = {}
    identity: dict[str, str] = {}
    for c in components:
        for sid, p in c.ports.items():
            if p.entity_id:
                if sid in identity and identity[sid] != p.entity_id:
                    raise ValueError(f"identity conflict on '{sid}': {c.id} says entity_id "
                                     f"{p.entity_id!r} but a prior component says {identity[sid]!r} -- a "
                                     f"different canonical identity means a different species.")
                identity[sid] = p.entity_id
            if sid in chemical and chemical[sid] != p.chemical:
                raise ValueError(f"identity conflict on '{sid}': {c.id} says chemical={p.chemical} but a "
                                 f"prior component says chemical={chemical[sid]}. A species is chemical "
                                 f"(mass-checkable) or not, globally -- boxes cannot disagree.")
            chemical[sid] = p.chemical
            if p.chemical:
                prev, cur = grounding.get(sid), p.formula   # formulas are canonicalised at load
                if prev and cur and prev != cur:            # conflict ONLY if both declared and differ
                    raise ValueError(f"identity conflict on '{sid}': {c.id} says atoms {cur!r} but a prior "
                                     f"component says {prev!r}. The shared medium pool has one atom count -- "
                                     f"boxes cannot disagree.")
                grounding[sid] = prev or cur                # a null formula is a GAP -> fill it, not conflict
    return {s: f for s, f in grounding.items() if f}      # only chemically-grounded species


class BlackBoxProcess(Process):
    defaults = {"transfer": lambda a, dt: {}, "species": []}

    def ports_schema(self):
        return {"medium": {s: {"_default": 0.0, "_emit": True, "_updater": "accumulate"}
                           for s in self.parameters["species"]}}

    def next_update(self, dt, states):
        avail = dict(states["medium"])
        return {"medium": self.parameters["transfer"](avail, dt)}


def _atom_flux(deltas: dict, grounding: dict) -> dict:
    """Net atoms a box moved into the medium this step = sum(delta * formula) over its grounded species."""
    flux: dict[str, float] = {}
    for sp, d in deltas.items():
        for el, n in parse_formula(grounding.get(sp)).items():
            flux[el] = flux.get(el, 0.0) + n * d
    return flux


def _recording(transfer, sink: list):
    """Wrap a transfer so every applied delta dict is appended to `sink` -- the actual per-step boundary
    flux of an open box (recorded, not reconstructed, because FBA optima are degenerate)."""
    def wrapped(avail, dt):
        deltas = transfer(avail, dt)
        sink.append(dict(deltas))
        return deltas
    return wrapped


def run_netlist(components, netlist: dict):
    """components: list[Component]. netlist: {initial_medium, duration, dt?}. Returns (series, report).

    OPEN-BOUNDARY handling: a box flagged `open_boundary` exchanges matter with an
    UNMODELLED compartment (atmosphere, feed, biomass) -- e.g. FBA drains atoms into biomass, aeration
    trades O2 with the air. Its net port atom-flux each step is a DECLARED boundary the monitor subtracts,
    so a legitimately-open box does not read as a leak, while any leak BEYOND the declaration (a bug in a
    closed box, a mis-composition) is still caught. Limitation: an open box's own internal imbalance is
    trusted (a biomass-composition cross-check is future work)."""
    if not _HAS_VIVARIUM:
        raise RuntimeError("run_netlist needs vivarium (the netlist simulator); not installed in this environment")
    composition_problems = composition_validity(components)  # is this even POSSIBLE? (static, pre-sim gate)
    grounding = medium_grounding(components)               # runs the identity guard
    all_species = sorted({s for c in components for s in c.ports})
    procs, topo, records = {}, {}, {}
    for c in components:
        wt = unit_reconciled(build_transfer(c), c.ports)
        if getattr(c, "open_boundary", False):
            # RECORD the actual applied deltas (FBA has degenerate optima -- re-solving would give a
            # different flux distribution; only the deltas actually applied are the true boundary flux)
            records[c.id] = []
            wt = _recording(wt, records[c.id])
        procs[c.id] = BlackBoxProcess({"transfer": wt, "species": list(c.ports.keys())})
        topo[c.id] = {"medium": ("medium",)}
    init = {"medium": {s: float(netlist.get("initial_medium", {}).get(s, 0.0)) for s in all_species}}
    eng = Engine(processes=procs, topology=topo, initial_state=init, display_info=False, progress_bar=False)
    eng.update(float(netlist.get("duration", 20.0)))
    data = eng.emitter.get_data()
    times = sorted(data)
    series = [{s: data[t]["medium"].get(s, 0.0) for s in all_species} for t in times]

    # NUMERICAL-DIVERGENCE guard: the fixed-step accumulate integrator can blow up on a STIFF composition
    # (a real integrator-stability limit -- NOT a modelling truth, and NOT something to paper over with a
    # fabricated value). A non-finite / None medium value means the trajectory diverged; keep the finite
    # prefix (physically meaningful up to divergence) and report the divergence honestly, rather than
    # crashing a downstream unit conversion or the conservation check with an opaque error. The proper cure
    # is an adaptive/stiff integrator (future work); until then a diverged netlist is a truthful outcome.
    def _finite(st):
        return all(v is not None and math.isfinite(v) for v in st.values())
    diverged = None
    for i, st in enumerate(series):
        if not _finite(st):
            diverged = times[i]
            series, times = series[:i], times[:i]
            break

    # declare each open box's per-step boundary flux from the ACTUAL applied deltas (recorded above)
    exchanges = None
    if records and len(series) > 1:
        n = len(series) - 1
        exchanges = []
        for t in range(n):
            net: dict[str, float] = {}
            for rec in records.values():
                if t < len(rec):
                    for el, v in _atom_flux(rec[t], grounding).items():
                        net[el] = net.get(el, 0.0) + v
            exchanges.append(net)

    # charge-balance too: an acid/base box (carbonate) must conserve charge, not only atoms
    charges = {sid: p.charge for c in components for sid, p in c.ports.items()
               if p.chemical and p.charge is not None}
    report = ConservationMonitor(grounding, charges).check(
        [{s: st[s] for s in grounding} for st in series], exchanges=exchanges)
    # validity domains checked against the driving conditions (initial state); honest, not enforced-away
    warnings = sorted({w for c in components for w in c.validity_warnings(series[0] if series else {})})
    return {"times": times, "species": all_species, "series": series,
            "validity_warnings": warnings, "declared_exchanges": exchanges is not None,
            "diverged": diverged, "composition_problems": composition_problems}, report
