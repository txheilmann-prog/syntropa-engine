"""Built-in transfer functions referenced by name from `function` manifests.

SOURCING (every number/formula is authoritative or explicitly labelled nominal):
  - STOICHIOMETRY of every box is an authoritative textbook net equation (atom-balanced; the conservation
    check verifies it).
  - photosynthesis / aerobic_respiration / ethanol_fermentation: the first-order rate law is a
    PHENOMENOLOGICAL, NOMINAL placeholder -- a tunable parameter, NOT a sourced kinetic model. These are
    illustrative net-reaction boxes; a REAL kinetic model is wrapped as an SBML or FBA box instead
    (Teusink/Pritchard/Kholodenko/e_coli_core).
  - carbonate_speciation: mass-action relaxation whose EQUILIBRIUM is real -- the rate-constant ratios are
    the measured pK1/pK2 supplied by the manifest; the relaxation timescale is nominal.
  - aeration: O2 saturation is REAL (Benson-Krause / USGS solubility table in the manifest) x Henry's law.
"""
from .transfer import register
from . import chem_additives  # noqa: F401  (registers the chemistry-additive primitives: bridge/catalyze/modify)


@register("photosynthesis")
def photosynthesis(avail, dt, params):
    """Stoichiometry (authoritative, textbook): 6 CO2 + 6 H2O -> C6H12O6 + 6 O2, atom-balanced.
    Kinetics: PHENOMENOLOGICAL first-order min-limiting law; the rate constant is a nominal tunable
    parameter, NOT a sourced measurement -- an illustrative net-reaction box, not a validated model."""
    k = params.get("rate", {}).get("default", 0.3) if isinstance(params.get("rate"), dict) else params.get("rate", 0.3)
    r = k * min(avail.get("light", 0.0), avail.get("CO2", 0.0) / 6.0, avail.get("H2O", 0.0) / 6.0) * dt
    return {"light": -6 * r, "CO2": -6 * r, "H2O": -6 * r, "glucose": +r, "O2": +6 * r}


@register("ethanol_fermentation")
def ethanol_fermentation(avail, dt, params):
    """Stoichiometry (authoritative, textbook): C6H12O6 -> 2 C2H6O + 2 CO2, atom-balanced.
    Kinetics: PHENOMENOLOGICAL first-order in glucose; rate constant is a nominal tunable parameter, NOT a
    sourced kinetic model."""
    k = params.get("rate", {}).get("default", 0.5) if isinstance(params.get("rate"), dict) else params.get("rate", 0.5)
    r = k * avail.get("glucose", 0.0) * dt
    return {"glucose": -r, "ethanol": +2 * r, "CO2": +2 * r}


@register("aerobic_respiration")
def aerobic_respiration(avail, dt, params):
    """Stoichiometry (authoritative, textbook): C6H12O6 + 6 O2 -> 6 CO2 + 6 H2O, atom-balanced (the mirror
    of photosynthesis). Kinetics: PHENOMENOLOGICAL first-order, glucose/O2-limited; rate constant is a
    nominal tunable parameter, NOT a sourced kinetic model. Reuses only already-grounded species."""
    k = params.get("rate", {}).get("default", 0.2) if isinstance(params.get("rate"), dict) else params.get("rate", 0.2)
    r = k * min(avail.get("glucose", 0.0), avail.get("O2", 0.0) / 6.0) * dt
    return {"glucose": -r, "O2": -6 * r, "CO2": +6 * r, "H2O": +6 * r}


@register("carbonate_speciation")
def carbonate_speciation(avail, dt, params):
    """Carbonate/bicarbonate equilibrium, mass-action relaxation (atom- AND charge-balanced BY
    CONSTRUCTION -- the conservation monitor verifies both):
        R1: CO2 + H2O <-> HCO3- + H+       R2: HCO3- <-> CO3-- + H+
    Species carry their ACTUAL protonation (HCO3- = C H O3, charge -1; CO3-- = C O3, charge -2; H+ = H,
    charge +1) -- NOT a neutralized formula, so H and charge accounting stay honest and pH is real.
    The rate constants come FROM THE MANIFEST (required, no fabricated default): their ratios reproduce
    the measured pK1/pK2, so the equilibrium and pH are physically real; the relaxation timescale is
    nominal. The proton pool is what a pH readout reads (via the acidity/log bridge)."""
    co2, h2o = avail.get("CO2", 0.0), avail.get("H2O", 0.0)
    hco3, co3, h = avail.get("bicarbonate", 0.0), avail.get("carbonate", 0.0), avail.get("proton", 0.0)
    v1 = (params["kf1"] * co2 - params["kr1"] * hco3 * h) * dt      # kf1/kr1 = Ka1 (real pK1)
    v2 = (params["kf2"] * hco3 - params["kr2"] * co3 * h) * dt      # kf2/kr2 = Ka2 (real pK2)
    return {"CO2": -v1, "H2O": -v1, "bicarbonate": v1 - v2, "carbonate": v2, "proton": v1 + v2}


def _interp(x, table):
    """Linear interpolation over a sorted [[x0,y0],[x1,y1],...] table; clamps outside the range."""
    if x <= table[0][0]:
        return table[0][1]
    if x >= table[-1][0]:
        return table[-1][1]
    for (x0, y0), (x1, y1) in zip(table, table[1:]):
        if x0 <= x <= x1:
            return y0 + (y1 - y0) * (x - x0) / (x1 - x0)
    return table[-1][1]


@register("aeration")
def aeration(avail, dt, params):
    """Surface O2 exchange, Henry's-law relaxation to saturation: dO2 = kLa * (O2_sat(T,P) - O2) * dt.
    O2_sat is REAL: interpolated from the manifest's Benson-Krause / USGS dissolved-oxygen solubility
    table (mg/L vs degC at 1 atm), scaled by pressure per Henry's law (O2_sat proportional to P). No
    fabricated temperature term -- the temperature dependence is the measured solubility curve. kLa is a
    genuine, system-specific surface-transfer coefficient (a real modelling parameter, tunable per
    system). Units are the box's own -- O2 mg/L, T degC, P atm -- reconciled by the engine."""
    o2 = avail.get("O2", 0.0)                        # mg/L
    t_c = avail.get("temperature", 20.0)            # degC
    p_atm = avail.get("pressure", 1.0)              # atm
    sat = _interp(t_c, params["O2_solubility_mgL"]) * p_atm      # real solubility(T) * Henry P-proportionality
    kla = params.get("kLa", 0.5)                    # 1/time surface-transfer coefficient (system-specific)
    return {"O2": kla * (sat - o2) * dt}            # mg/L delta
