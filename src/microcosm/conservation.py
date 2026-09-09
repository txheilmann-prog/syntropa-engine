"""
Composition-boundary conservation monitor -- part of microcosm's physics core. Standard atom (and, if
charges are given, charge) bookkeeping using microcosm's own units.parse_formula. Simulator-agnostic: given a grounding (species id -> chemical formula) and a time
series of species amounts, it checks whole-system atom (and, if charges are given, charge) conservation
over time, subtracting DECLARED open-boundary exchanges.

For a CLOSED composite (no declared influx/efflux) the atom totals must be constant; any drift localizes a
conservation bug introduced by the composition, named per element. Open-boundary fluxes (the EXCHANGE
analogue) are declared and subtracted (see `exchanges`).
"""
from __future__ import annotations

import math

from dataclasses import dataclass

from .units import parse_formula


@dataclass
class ConservationReport:
    elements: dict[str, float]          # element -> max |drift from t0| over the series
    leaking: list[str]                  # elements whose drift exceeds tolerance
    n_steps: int

    @property
    def conserved(self) -> bool:
        return not self.leaking

    def summary(self) -> str:
        if self.conserved:
            return f"atom-conserved across {self.n_steps} steps (max drift negligible)"
        worst = ", ".join(f"{e} {self.elements[e]:+.3g}" for e in self.leaking)
        return f"CONSERVATION VIOLATION across {self.n_steps} steps -- leaking: {worst}"


class ConservationMonitor:
    """Checks atom (and, if charges given, charge) conservation of a composed system over time."""

    def __init__(self, grounding: dict[str, str], charges: dict[str, float] | None = None):
        self.formulas = {sp: parse_formula(f) for sp, f in grounding.items() if f}
        self.charges = dict(charges or {})

    def atom_totals(self, state: dict[str, float]) -> dict[str, float]:
        """Total atoms (and pseudo-element 'charge') in a single state = sum species amount * formula."""
        totals: dict[str, float] = {}
        for sp, amt in state.items():
            f = self.formulas.get(sp)
            if f:
                for el, n in f.items():
                    totals[el] = totals.get(el, 0.0) + n * amt
            if sp in self.charges:
                totals["charge"] = totals.get("charge", 0.0) + self.charges[sp] * amt
        return totals

    def check(self, series: list[dict[str, float]], *,
              exchanges: list[dict[str, float]] | None = None,
              rtol: float = 1e-6, atol: float = 1e-9) -> ConservationReport:
        """`series`: species-amount state at each timestep (index 0 = initial). `exchanges`: optional
        per-step declared net atom influx (positive = added to the system) to subtract before checking,
        for open composites. Returns per-element max drift from the (exchange-adjusted) initial totals."""
        if not series:
            return ConservationReport({}, [], 0)
        base = self.atom_totals(series[0])
        cum: dict[str, float] = {}       # cumulative declared exchange up to step t
        drift: dict[str, float] = {e: 0.0 for e in base}
        for t, state in enumerate(series):
            cur = self.atom_totals(state)
            for e in set(cur) | set(base):
                expected = base.get(e, 0.0) + cum.get(e, 0.0)
                d = abs(cur.get(e, 0.0) - expected)
                if not math.isfinite(d):               # a NaN/inf amount is a broken trajectory, never "conserved"
                    d = float("inf")
                drift[e] = max(drift.get(e, 0.0), d)
            if exchanges and t < len(exchanges):
                for e, v in exchanges[t].items():
                    cum[e] = cum.get(e, 0.0) + v
        leaking = [e for e, d in drift.items() if d > atol + rtol * abs(base.get(e, 0.0))]
        leaking.sort(key=lambda e: -drift[e])
        return ConservationReport(drift, leaking, len(series))
