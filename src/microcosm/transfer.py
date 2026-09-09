"""
Transfer functions: turn a component's `transfer` spec into a callable
    transfer(available: dict[species->amount], dt: float) -> dict[species->delta]
so every black box -- a python function, an experimental lookup table, or a whole SBML model -- is
driven identically by the engine. New kinds register via TRANSFER_BUILDERS (the extensibility seam).

  function : a registered named callable (manifests are data; code is referenced by name).
  lookup   : an experimental table {input-condition -> output rates}, nearest-neighbour for v0.1.
  sbml     : a real kinetic model, stepped by roadrunner via state-exchange each dt (co-simulation).
             v0.1 assumes unit compartment volume (amount ~ concentration); a known schema limitation.
"""
from __future__ import annotations

# registry of named python transfer functions, referenced by manifests
FUNCTION_REGISTRY: dict[str, callable] = {}


def register(name):
    def deco(fn):
        FUNCTION_REGISTRY[name] = fn
        return fn
    return deco


def _build_function(spec, comp):
    name = spec.get("name")
    if name not in FUNCTION_REGISTRY:
        raise ValueError(f"function transfer references unregistered '{name}'")
    fn = FUNCTION_REGISTRY[name]
    params = comp.params
    return lambda avail, dt: fn(avail, dt, params)


def _build_lookup(spec, comp):
    # spec: {"inputs": [sp,...], "table": [{"when": {sp:val,...}, "rates": {sp:rate,...}}, ...]}
    table = spec["table"]
    keys = spec["inputs"]

    def transfer(avail, dt):
        # nearest row by Euclidean distance in the declared input space
        pt = [avail.get(k, 0.0) for k in keys]
        best = min(table, key=lambda r: sum((r["when"].get(k, 0.0) - v) ** 2 for k, v in zip(keys, pt)))
        return {sp: rate * dt for sp, rate in best["rates"].items()}
    return transfer


def _build_sbml(spec, comp):
    import roadrunner
    roadrunner.Logger.setLevel(roadrunner.Logger.LOG_ERROR)
    if comp.payload is None:
        raise ValueError(f"sbml component {comp.id} has no payload (SBML content) loaded")
    sbml = comp.payload.decode("utf-8") if isinstance(comp.payload, (bytes, bytearray)) else comp.payload
    port_map = spec.get("port_map", {})            # model-species -> medium-species
    rr = roadrunner.RoadRunner(sbml)
    model_ids = set(rr.model.getFloatingSpeciesIds())
    mapped = {m: s for m, s in port_map.items() if m in model_ids}
    if not mapped:
        raise ValueError(f"sbml transfer: no port_map species found in model {comp.id}")

    def transfer(avail, dt):
        # push current medium levels into the model's mapped species, advance dt, read back deltas
        before = {}
        for m, s in mapped.items():
            if s in avail:
                try:
                    rr[m] = avail[s]
                except Exception:
                    pass
            before[m] = rr[m]
        rr.simulate(0, dt, 2)
        return {s: rr[m] - before[m] for m, s in mapped.items()}
    return transfer


def _load_fba_model(comp):
    """Load a cobra model from an fba component's SBML payload (cobra stays behind this seam)."""
    import cobra
    import os
    import tempfile
    if comp.payload is None:
        raise ValueError(f"fba component {comp.id} has no payload (SBML content) loaded")
    with tempfile.NamedTemporaryFile("wb", suffix=".xml", delete=False) as tf:
        tf.write(comp.payload)
        tmp = tf.name
    try:
        cobra.io.sbml.LOGGER.disabled = True
        return cobra.io.read_sbml_model(tmp)
    finally:
        os.unlink(tmp)


def _fba_input_reactions(model, comp):
    """[(exchange-reaction id, medium species)] for the component's INPUT ports present in the model."""
    port_map = comp.transfer.get("port_map", {})
    return [(rid, port_map[rid]) for rid in port_map
            if model.reactions.has_id(rid) and comp.ports.get(port_map[rid])
            and comp.ports[port_map[rid]].role in ("input", "catalyst")]


def _build_fba(spec, comp):
    """Genome-scale metabolism as a black box via DYNAMIC FBA. A genome-scale model has no rate laws
    (roadrunner cannot drive it); instead, each dt we cap each INPUT exchange's uptake by what the
    medium currently holds, maximize the model's objective (growth), and return the exchange fluxes as
    medium deltas (uptake negative, secretion positive). This box is intrinsically OPEN: biomass drains
    atoms out of the tracked medium -- the exchanges ARE its declared boundary."""
    model = _load_fba_model(comp)
    port_map = spec.get("port_map", {})            # exchange-reaction id -> medium species
    vmax = float(spec.get("vmax", 10.0))           # max uptake rate cap (mmol/gDW/h scale)
    rxn_ids = [rid for rid in port_map if model.reactions.has_id(rid)]
    if not rxn_ids:
        raise ValueError(f"fba transfer: no port_map exchange reactions found in model {comp.id}")
    inputs = dict(_fba_input_reactions(model, comp))

    def transfer(avail, dt):
        with model:                                # scratch context: bound edits revert after solve
            for rid in rxn_ids:
                if rid in inputs:                                  # cap uptake by availability (dFBA)
                    amt = max(0.0, avail.get(inputs[rid], 0.0))
                    rxn = model.reactions.get_by_id(rid)
                    rxn.upper_bound = max(rxn.upper_bound, 0.0)     # never FORCE uptake -- the medium controls
                    rxn.lower_bound = -min(vmax, amt / dt) if amt > 0 else 0.0
            sol = model.optimize()
            if sol.status != "optimal":
                return {}
            return {port_map[rid]: float(sol.fluxes[rid]) * dt for rid in rxn_ids}
    return transfer


def fba_analyzer(comp):
    """Return growth_under(availability) -> the REAL FBA growth rate (objective value) for the component
    under a given nutrient availability. Standard flux-balance essentiality analysis (Orth-Thiele-Palsson
    2010, Nat Biotechnol): a nutrient present -> its uptake is allowed (up to vmax); absent -> its
    exchange is closed (a knockout). Loads the model once; cobra stays behind this seam."""
    model = _load_fba_model(comp)
    vmax = float(comp.transfer.get("vmax", 10.0))
    inputs = _fba_input_reactions(model, comp)

    def growth_under(availability):
        with model:
            for rid, sp in inputs:
                present = availability.get(sp, 0.0) > 0
                rxn = model.reactions.get_by_id(rid)
                rxn.upper_bound = max(rxn.upper_bound, 0.0)         # never FORCE uptake; availability controls
                rxn.lower_bound = -vmax if present else 0.0
            sol = model.optimize()
            return float(sol.objective_value) if sol.status == "optimal" else 0.0
    return growth_under


TRANSFER_BUILDERS = {"function": _build_function, "lookup": _build_lookup,
                     "sbml": _build_sbml, "fba": _build_fba}


def build_transfer(comp):
    t = comp.transfer.get("type")
    if t not in TRANSFER_BUILDERS:
        raise ValueError(f"unknown transfer type {t!r}")
    return TRANSFER_BUILDERS[t](comp.transfer, comp)
