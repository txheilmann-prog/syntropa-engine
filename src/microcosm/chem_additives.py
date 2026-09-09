"""
CHEMISTRY-ADDITIVE LAYER -- the abiotic organic/inorganic boxes that BRIDGE, CATALYZE, and MODIFY the
shared medium. This is the connective tissue between organism black boxes: it is where under-explored,
possibly-novel A->B routes live (most metabolic modeling ignores the abiotic chemistry that couples cells).

Design: TWO generic, DATA-DRIVEN transfer primitives so the layer scales by adding MANIFESTS, not code.
Every additive box is then just a manifest carrying REAL equilibrium/kinetic constants: the
constant RATIOS are physically real -- an acid-dissociation Ka, a stability constant K, a Henry constant;
only the absolute relaxation timescale is nominal and is labelled as such, exactly as carbonate_speciation
already does). New chemistry = new manifest data, reviewed for its constants.

  mass_action        : reversible law-of-mass-action over a list of elementary reactions. kf/kr ratio = Keq
                       (encode the real pKa / logK / Henry constant as that ratio in the manifest). BRIDGE
                       (speciation, phase transfer) and MODIFY (chelation, acid/base) boxes use this.
  catalyzed_reaction : the same, but each reaction rate is multiplied by the amount of a CATALYST species
                       that is NOT consumed (enzyme, mineral surface, metal). This is the CATALYZE role --
                       presence enables/accelerates a step. kcat is nominal (labelled); the driven
                       reaction's kf/kr ratio is still the real Keq so it cannot violate thermodynamics.

Law of mass action: Guldberg & Waage (authoritative). Reaction network + constants come from the manifest.
"""
from __future__ import annotations

from .transfer import register


def _reaction_rate(rxn, avail, scale=1.0):
    """Net forward extent-rate of one elementary reversible reaction: (kf*[reactants] - kr*[products])."""
    kf, kr = rxn["kf"], rxn["kr"]
    fwd = kf
    for s, n in rxn["reactants"].items():
        fwd *= avail.get(s, 0.0) ** n
    rev = kr
    for s, n in rxn["products"].items():
        rev *= avail.get(s, 0.0) ** n
    return (fwd - rev) * scale


def _apply(rxn, v, deltas):
    for s, n in rxn["reactants"].items():
        deltas[s] = deltas.get(s, 0.0) - n * v
    for s, n in rxn["products"].items():
        deltas[s] = deltas.get(s, 0.0) + n * v


@register("mass_action")
def mass_action(avail, dt, params):
    """Reversible mass-action over params['reactions'] = [{reactants:{sp:stoich}, products:{sp:stoich},
    kf, kr}, ...]. kf/kr encodes the real equilibrium constant; the network is atom/charge-balanced BY
    CONSTRUCTION in the manifest and the conservation monitor verifies it."""
    deltas: dict[str, float] = {}
    for rxn in params["reactions"]:
        _apply(rxn, _reaction_rate(rxn, avail) * dt, deltas)
    return deltas


@register("catalyzed_reaction")
def catalyzed_reaction(avail, dt, params):
    """Catalyst-enabled mass-action: each reaction's rate is multiplied by kcat * amount(catalyst). The
    catalyst species is NOT consumed (it never appears in reactants/products), so it only ENABLES/speeds
    the transform -- the CATALYZE role. The driven reaction still carries the real Keq via kf/kr, so a
    catalyst cannot push a reaction past its thermodynamic equilibrium (no free lunch)."""
    cat = avail.get(params["catalyst"], 0.0)
    kcat = params.get("kcat", 1.0)
    deltas: dict[str, float] = {}
    if cat <= 0.0:
        return deltas                                  # no catalyst present -> the step stays dead
    for rxn in params["reactions"]:
        _apply(rxn, _reaction_rate(rxn, avail, scale=kcat * cat) * dt, deltas)
    return deltas
