"""
Units + dimensions -- the SSoT for physical-quantity discipline. Every
medium quantity has a DIMENSION with ONE canonical unit; every box
port declares the unit IT speaks; the engine reconciles at the boundary.

Two kinds of conversion, kept distinct on purpose:
  1. WITHIN-dimension unit conversion  (mg/L <-> mM, degC <-> K, atm <-> bar) -- `to_canonical`.
  2. CROSS-dimension derived relationships, each needing OTHER state, as named functions:
       concentration <-> amount   via VOLUME     (concentration = amount / volume)
       pH            <-> [H+]      via LOG        ([H+] = 10^-pH)          (reserved)
       gas amount    <-> volume    via PV=nRT      (ideal gas)             (reserved)
So `to_canonical` never silently crosses dimensions; crossing is always an explicit, state-aware call.

Discipline: `validate_unit` is called at component load -- an unknown unit, a chemical species carrying
a non-substance unit, or a mass-concentration unit with no molar mass is a LOUD error, not a runtime
surprise. This is the load-time discipline for physical quantities.
"""
from __future__ import annotations

import math
import re

_ELEMENT_RE = re.compile(r"([A-Z][a-z]?)(\d*)")

# IUPAC standard atomic weights (g/mol), common biological elements; extend as needed.
ATOMIC_WEIGHTS = {  # IUPAC standard atomic weights (2021); extended for biological cofactors
    "H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "P": 30.974, "S": 32.06,
    "Na": 22.990, "K": 39.098, "Cl": 35.45, "Ca": 40.078, "Mg": 24.305, "Fe": 55.845,
    # metals + metalloids + halogens that carry real cofactor chemistry (B12=Co, molybdopterin=Mo, SeCys=Se, ...)
    "Zn": 65.38, "Mn": 54.938, "Cu": 63.546, "Co": 58.933, "Ni": 58.693, "Mo": 95.95, "Se": 78.971,
    "B": 10.81, "F": 18.998, "I": 126.904, "Br": 79.904, "Si": 28.085, "Al": 26.982, "As": 74.922,
    "W": 183.84, "Cd": 112.414, "Hg": 200.592, "V": 50.942, "Sr": 87.62, "Ba": 137.327}

# dimension -> (canonical unit, extensivity). EXTENSIVE quantities add when systems combine (the engine
# accumulates their deltas); INTENSIVE ones do not (mixing two 20C volumes gives 20C, not 40C).
DIMENSIONS = {
    "amount":        ("mmol", "extensive"),
    "concentration": ("mM",   "intensive"),    # backed by an extensive amount via volume
    "volume":        ("L",    "extensive"),
    "mass":          ("g",    "extensive"),
    "temperature":   ("K",    "intensive"),
    "pressure":      ("bar",  "intensive"),
    "potential":     ("mV",   "intensive"),    # membrane / redox potential (reserved)
    "acidity":       ("pH",   "intensive"),    # pH, log-related to [H+] (reserved)
    "dimensionless": ("",     "intensive"),
}

# Substance quantities (a chemical species must use one of these); the rest are environment conditions.
SUBSTANCE_DIMS = {"amount", "concentration"}
# Intensive environment conditions a box may READ but not emit as an accumulating delta (no combine rule
# yet). concentration is intensive but proxies an extensive amount, so it IS accumulable.
NONACCUMULABLE_DIMS = {"temperature", "pressure", "potential", "acidity"}

# unit -> (dimension, transform). transforms (WITHIN-dimension):
#   ("linear", scale)          canonical = value * scale
#   ("affine", scale, offset)  canonical = value * scale + offset   (offset skipped when converting a delta)
#   ("mass_conc", factor)      mass concentration -> mM via molar mass: canonical = value*factor / molar_mass
UNITS = {
    # amount (extensive)
    "mmol": ("amount", ("linear", 1.0)), "umol": ("amount", ("linear", 1e-3)),
    "mol": ("amount", ("linear", 1000.0)),
    "mmol/gDW/h": ("amount", ("linear", 1.0)),          # FBA flux*dt with gDW=1 (documented shortcut)
    # concentration -- molar
    "mM": ("concentration", ("linear", 1.0)), "uM": ("concentration", ("linear", 1e-3)),
    "M": ("concentration", ("linear", 1000.0)),
    # concentration -- mass (needs molar mass)
    "g/L": ("concentration", ("mass_conc", 1e3)), "mg/L": ("concentration", ("mass_conc", 1.0)),
    "ug/L": ("concentration", ("mass_conc", 1e-3)),
    # volume (extensive)
    "L": ("volume", ("linear", 1.0)), "mL": ("volume", ("linear", 1e-3)),
    "uL": ("volume", ("linear", 1e-6)), "fL": ("volume", ("linear", 1e-15)),
    # mass (extensive)
    "g": ("mass", ("linear", 1.0)), "mg": ("mass", ("linear", 1e-3)),
    "gDW": ("mass", ("linear", 1.0)),                   # grams dry weight (biomass basis)
    # temperature (affine)
    "K": ("temperature", ("affine", 1.0, 0.0)), "degC": ("temperature", ("affine", 1.0, 273.15)),
    "degF": ("temperature", ("affine", 5.0 / 9.0, 273.15 - 32.0 * 5.0 / 9.0)),
    # pressure
    "bar": ("pressure", ("linear", 1.0)), "atm": ("pressure", ("linear", 1.01325)),
    "kPa": ("pressure", ("linear", 0.01)), "Pa": ("pressure", ("linear", 1e-5)),
    "mmHg": ("pressure", ("linear", 0.00133322)), "psi": ("pressure", ("linear", 0.0689476)),
    # potential (reserved)
    "mV": ("potential", ("linear", 1.0)), "V": ("potential", ("linear", 1000.0)),
    # acidity (reserved; pH is canonical of its dimension -- the log link to [H+] is a cross-dim bridge)
    "pH": ("acidity", ("linear", 1.0)),
    # dimensionless
    "dimensionless": ("dimensionless", ("linear", 1.0)), "": ("dimensionless", ("linear", 1.0)),
}

CANONICAL_UNIT = "mmol"                                  # the default (amount) canonical unit
R_GAS = 0.0831446                                        # L*bar / (mmol...) gas constant, reserved for PV=nRT
STANDARD_GRAVITY = 9.80665                               # m/s^2; a CONSTANT (its altitude variation ~0.3%/km
#   is negligible). Gravity never enters cellular biochemistry directly (molecules are too small; thermal
#   motion dominates). It matters only through PRESSURE: barometric pressure falls with altitude (-> gas
#   solubility, Henry's law) and hydrostatic pressure rises with water depth (P = rho*g*h). Both SET the
#   pressure dimension; g is a parameter of the hydrostatic relationship, not a medium quantity. Cell
#   sedimentation/buoyancy (also ~g) only appears once the medium is SPATIAL -- reserved with geometry.
WATER_DENSITY = 1000.0                                   # kg/m^3, for hydrostatic pressure


def parse_formula(formula: str | None) -> dict[str, float]:
    """Non-strict formula -> {element: count}. Empty/None -> {}."""
    atoms: dict[str, float] = {}
    for el, n in _ELEMENT_RE.findall(formula or ""):
        if el:
            atoms[el] = atoms.get(el, 0.0) + (float(n) if n else 1.0)
    return atoms


def canonical_formula(formula: str | None) -> str | None:
    """Normalise a formula to a canonical Hill string (C, H, then alphabetical; omit count 1) so that
    different WRITINGS of the same molecule compare equal -- 'CH4O1' and 'CH4O' both -> 'CH4O', 'CH1O2' and
    'CHO2' both -> 'CHO2'. This popped out of mass-ingesting BiGG GEMs (explicit-1 vs implicit-1)."""
    atoms = parse_formula(formula)
    if not atoms:
        return None
    order = (["C"] if "C" in atoms else []) + (["H"] if "H" in atoms else []) \
        + sorted(e for e in atoms if e not in ("C", "H"))

    def part(e):
        n = atoms[e]
        n = int(n) if float(n).is_integer() else n
        return e + ("" if n == 1 else str(n))
    return "".join(part(e) for e in order)


def molar_mass(formula: str | None) -> float | None:
    """g/mol from a grounded formula, or None if ungrounded / has an element we don't tabulate. Returns None for
    VARIABLE-LENGTH polymers / R-groups: '(C6H10O5)n' (glycogen), 'C20H36O(C5H8)n', anything with
    parens / '*' / a wildcard R-group has NO fixed molecular mass -- the old code silently returned the monomer
    mass (a wrong number that then poisoned canonical_formula + any mass check)."""
    if formula and re.search(r"[()*]|R(?![abefghnu])", formula):     # polymer/complex/R-group -> no fixed mass
        return None
    atoms = parse_formula(formula)
    if not atoms:
        return None
    total = 0.0
    for el, n in atoms.items():
        if el not in ATOMIC_WEIGHTS:
            return None
        total += ATOMIC_WEIGHTS[el] * n
    return total


def dimension_of(unit: str | None) -> str | None:
    e = UNITS.get(unit or "")
    return e[0] if e else None


def canonical_unit(dimension: str) -> str:
    return DIMENSIONS[dimension][0]


def is_extensive(unit: str | None) -> bool:
    dim = dimension_of(unit)
    return dim is not None and DIMENSIONS[dim][1] == "extensive"


def to_canonical(value: float, unit: str | None, formula: str | None = None, *, delta: bool = False) -> float:
    """A value in `unit` -> its dimension's canonical unit (WITHIN-dimension only). `delta=True` means a
    difference/rate, so an affine OFFSET is skipped (5 degC of heating is 5 K, not 278 K)."""
    e = UNITS.get(unit or "")
    if e is None:
        if not unit:
            return value
        raise ValueError(f"unknown unit {unit!r}")
    kind = e[1][0]
    if kind == "linear":
        return value * e[1][1]
    if kind == "affine":
        return value * e[1][1] + (0.0 if delta else e[1][2])
    if kind == "mass_conc":
        m = molar_mass(formula)
        if not m:
            raise ValueError(f"cannot convert {unit} without a molar mass (formula={formula!r})")
        return value * e[1][1] / m                       # -> mM
    raise ValueError(f"unknown transform for unit {unit!r}")


def from_canonical(value: float, unit: str | None, formula: str | None = None, *, delta: bool = False) -> float:
    """Canonical unit -> a value in `unit` (inverse of to_canonical)."""
    e = UNITS.get(unit or "")
    if e is None:
        if not unit:
            return value
        raise ValueError(f"unknown unit {unit!r}")
    kind = e[1][0]
    if kind == "linear":
        return value / e[1][1]
    if kind == "affine":
        return (value - (0.0 if delta else e[1][2])) / e[1][1]
    if kind == "mass_conc":
        m = molar_mass(formula)
        if not m:
            raise ValueError(f"cannot convert {unit} without a molar mass (formula={formula!r})")
        return value * m / e[1][1]
    raise ValueError(f"unknown transform for unit {unit!r}")


# --- cross-dimension derived relationships (each needs OTHER state; never hidden inside to_canonical) --

def concentration_to_amount(conc_mM: float, volume_L: float = 1.0) -> float:
    """concentration (mM, canonical) * volume (L) -> amount (mmol). V=1 is the documented v1 shortcut;
    volume is a reserved first-class quantity so this becomes real once compartments carry volumes."""
    return conc_mM * volume_L


def amount_to_concentration(amount_mmol: float, volume_L: float = 1.0) -> float:
    return amount_mmol / volume_L if volume_L else amount_mmol


def ph_to_concentration(ph: float) -> float:
    """pH -> [H+] in mM (canonical concentration). LOG transform (the reserved acidity<->concentration
    bridge). Uses concentration, an approximation of activity (ideal-solution assumption)."""
    return (10.0 ** (-ph)) * 1000.0


def concentration_to_ph(cH_mM: float) -> float:
    if cH_mM <= 0:
        raise ValueError("pH undefined for non-positive [H+]")
    return -math.log10(cH_mM / 1000.0)


def barometric_pressure(altitude_m: float, sea_level_bar: float = 1.01325) -> float:
    """Atmospheric pressure (bar) at an altitude, ISA troposphere model (valid to ~11 km). This is how
    ALTITUDE enters the sim: lower pressure -> lower gas partial pressure -> lower dissolved-gas
    saturation (Henry's law). Sets the pressure dimension; g is folded into the ISA constants."""
    return sea_level_bar * (1.0 - 2.25577e-5 * altitude_m) ** 5.25588


def hydrostatic_pressure(depth_m: float, surface_bar: float = 1.01325,
                         density_kg_m3: float = WATER_DENSITY, g: float = STANDARD_GRAVITY) -> float:
    """Pressure (bar) at a water DEPTH: P = P_surface + rho*g*h (~1 atm per 10 m). This is where gravity
    literally appears -- as a constant in the hydrostatic term. Sets the pressure dimension for deep-water
    organisms."""
    return surface_bar + density_kg_m3 * g * depth_m / 1e5   # Pa -> bar


# --- discipline: validate a declared unit at component load (fail fast, not mid-simulation) -----------

def validate_unit(unit: str | None, *, chemical: bool, formula: str | None) -> str | None:
    """Return an error string if the unit is undisciplined, else None. Catches: unknown/mistyped units;
    a chemical species carrying a non-substance unit (e.g. O2 in degC); a mass-concentration unit with
    no molar mass to convert it."""
    if unit and unit not in UNITS:
        return f"unknown unit {unit!r} (not in the units registry)"
    dim = dimension_of(unit) or "dimensionless"
    if chemical and dim not in SUBSTANCE_DIMS:
        return (f"chemical species must use an amount/concentration unit, not a {dim} unit ({unit!r})")
    e = UNITS.get(unit or "")
    if e and e[1][0] == "mass_conc" and not molar_mass(formula):
        return f"unit {unit!r} needs a grounded formula (molar mass) to be convertible"
    return None
