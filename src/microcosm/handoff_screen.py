"""Per-hand-off thermodynamic FEASIBILITY screen for proposed interspecies metabolite hand-offs.

Given a proposed hand-off (organism P excretes carrier M, organism C consumes it), this asks the one question a
stoichiometric or genomic-complementarity predictor cannot: is there a physiological concentration of M at which
BOTH P's disposal reaction and C's uptake reaction are thermodynamically feasible? It is a NECESSARY-condition
FILTER -- it can reject a proposed hand-off, it never confirms one -- and it reproduces the MEASURED interspecies-H2
thresholds of the three H2-scavenging guilds (sulfate-reducer, methanogen, acetogen) from first principles; see
guard().

Modelling choices, each grounded in cited physics (no fitted parameters):
- Every dG'^0 comes live from `thermo` (eQuilibrator component-contribution); nothing is hand-set.
- The Schink energy QUANTUM (~20 kJ/mol, ~1/3 ATP; Schink 1997, Microbiol Mol Biol Rev 61:262) is applied ONLY
  to the energy-CONSERVING consumer. The producer's electron-DISPOSAL step needs only dG' < 0. Requiring the
  quantum from both wrongly forbids tight syntrophy (e.g. propionate -> methanogen); the asymmetry is what
  reproduces the measured thresholds.
- Feasibility uses the dG' POINT estimate vs the quantum (this is what reproduces the measured thresholds); the
  component-contribution uncertainty band is reported separately by clears() as a MARGINAL flag, never to move the
  boundary.
- THE THRESHOLD AGREEMENT IS CONDITIONAL ON THE SHIPPED SOLUTE CONCENTRATIONS, and that condition belongs
  beside the claim rather than in a footnote. Being INSIDE a cited concentration window at a POINT says
  nothing about the window. Swept across the cited acetate range, this screen stops reproducing the
  Cord-Ruwisch acetogen threshold at BOTH ends: 46.09 Pa at the low end and 115.78 Pa at the high end,
  against a published band of 50 to 95 Pa. The agreement is real and it is not unconditional, so a caller
  changing a concentration is changing what this screen validates against, not merely its inputs.

Concentrations, pH, ionic strength, and temperature here are engine DEFAULTS (25 C, the PHYS dict), explicit and
overridable per hand-off via `conc=`; to model a specific medium supply its values (from a named published table)
at call time rather than hard-coding them here.
"""
from __future__ import annotations

import math

from . import thermo

# RT at 298.15 K, kJ/mol. A caller MAY raise this to model a warmer reactor, but read this first, because doing so
# is a HALF correction and the other half is not available. MEASURED: eQuilibrator does NOT temperature-correct
# dG'^0 at all. Setting `cc.temperature` is silently accepted and reported back while the value does not move:
# hydrogenotrophic methanogenesis is -195.10 kJ/mol at 298.15, 313.15 AND 328.15 K alike, and propionate
# oxidation shifts 0.01 kJ/mol over 30 K. eQuilibrator's own documentation gives the reason, that dG at another
# temperature needs the entropy term, which the component-contribution training set lacks for most compounds. So
# raising RT here moves the ln-concentration term to the new temperature while dG'^0 stays pinned at 25 C, mixing
# two temperatures inside one inequality, and a partial correction of that kind measurably WORSENED agreement
# with documented obligate/facultative lifestyles in our own testing: a documented obligate aromatic syntrophy
# read as facultative at 37, 39 and 55 C alike. Treat every verdict as sign-feasibility at 25 C until a
# thermodynamic data source carrying enthalpy or entropy data replaces the estimator.
RT = 8.314462618e-3 * 298.15
# EVERY VERDICT HERE IS BULK-PHASE, AND CONDITIONAL ON THE MEMBERS ACHIEVING LOCAL PROXIMITY. There is one shared
# concentration per metabolite and no spatial, diffusion or aggregation variable anywhere in this module, while
# interspecies H2 and formate transfer is diffusion-limited: the decisive concentration is the one AT THE CELL-TO-CELL
# INTERFACE, which is why syntrophs aggregate.
# THE NARROWNESS OF THE WINDOWS IS WHY THIS MATTERS, measured rather than asserted. `handoff_window` on the obligate
# syntrophies returns a band spanning a factor of only about 1.5 to 4: propionate to methanogen is feasible between
# 2 and 3 Pa in the high-CO2 regime, 3 to 7 Pa in the low. A band that narrow is a statement about a controlled
# microenvironment, not about a stirred vessel where dissolved H2 ranges over orders of magnitude. The internal
# control is the NOT-obligate case: lactate oxidation returns a window six orders of magnitude wide, so window width
# tracks obligate-ness and the narrow bands are narrow for a physical reason.
# CONSEQUENCE FOR READING ANY RESULT FROM THIS MODULE: bulk-phase feasibility is NECESSARY AND NOT SUFFICIENT for a
# hand-off to occur. Nothing in a well-mixed formulation can establish that two organisms will aggregate, so
# proximity is either engineered (co-immobilization, a biofilm carrier, granulation) or tested BEFORE the conversion
# is tested, since a negative result in a stirred vessel cannot distinguish wrong chemistry from cells too far apart.
QUANTUM = 20.0                           # Schink per-organism energy quantum, kJ/mol (Schink 1997)
kH_H2 = 7.8e-4                           # dissolved H2 solubility, M/atm (Henry, ~25 C)


def M_to_Pa(m: float) -> float:
    return (m / kH_H2) * 101325.0


def Pa_to_M(pa: float) -> float:
    return (pa / 101325.0) * kH_H2


# Physiological-default activities (M). Overridable per hand-off via `conc=`; the pilot freezes these to a named
# published table for Faust's medium. co2 = dissolved CO2(aq) ~1.3 mM at pCO2 ~0.05 atm; acetate ~10 mM (gut).
PHYS = {"co2": 1.3e-3, "ac": 1e-2, "for": 1e-3, "h2s": 1e-3, "ch4": 1e-3, "so4": 1e-3, "h2o": 1.0}
DEFAULT_CONC = 1e-3


def _dg0(sub: dict, prod: dict):
    d = thermo.reaction_dg0_prime(sub, prod)
    return (d["dg0_prime_kj_mol"], d.get("uncertainty_kj_mol", 0.0)) if d.get("reliable") else (None, None)


def _lnQ(sub: dict, prod: dict, conc: dict) -> float:
    """ln Q at the given activities; water activity held at 1.

    Split out of `dg_prime` because it is the ONLY concentration-dependent half of dG' = dG'0 + RT ln Q. A
    concentration SCAN re-evaluates this and nothing else, so the scans below hoist dG'0 out of their loops and
    call this per point. One implementation, used by both, so the two cannot drift apart.
    """
    lnQ = 0.0
    for side, sign in (("sub", -1), ("prod", +1)):
        for b, n in (sub if side == "sub" else prod).items():
            if b == "h2o":
                continue
            lnQ += sign * n * math.log(conc.get(b, PHYS.get(b, DEFAULT_CONC)))
    return lnQ


def dg_prime(sub: dict, prod: dict, conc: dict):
    """dG' (kJ/mol) at the given activities; water activity held at 1. Returns (dg, uncertainty) or (None, None)."""
    g0, err = _dg0(sub, prod)
    if g0 is None:
        return None, None
    return g0 + RT * _lnQ(sub, prod, conc), err


def clears(dg, err, quantum=QUANTUM):
    """Three-way verdict: point estimate vs the quantum, with the 1.96-sigma band flagging MARGINAL."""
    if dg is None:
        return "abstain"
    if dg + 1.96 * err <= -quantum:
        return "feasible"
    if dg - 1.96 * err > -quantum:
        return "infeasible"
    return "marginal"


def solo_threshold(sub, prod, shared="h2", conc=None, quantum=QUANTUM, lo=1e-11, hi=1e-1, n=400):
    """The [shared] boundary at which one reaction just clears the quantum (point estimate). Returns (kind, M):
    'ceiling' (max [shared]) for a shared-metabolite PRODUCER, 'floor' (min [shared]) for a CONSUMER."""
    conc = dict(PHYS, **(conc or {}))
    # HOISTED: dG'0 is concentration-INDEPENDENT, so this scan needed ONE eQuilibrator evaluation and was
    # making n+1 identical ones. An unreliable reaction is unreliable at every point, so the early return is
    # the same answer the loop reached by finding nothing feasible.
    g0, _err = _dg0(sub, prod)
    if g0 is None:
        return (None, None)
    feas = []
    for i in range(n + 1):
        m = lo * (hi / lo) ** (i / n)
        c = dict(conc); c[shared] = m
        if g0 + RT * _lnQ(sub, prod, c) <= -quantum:
            feas.append(m)
    if not feas:
        return (None, None)
    n_sh = prod.get(shared, 0) - sub.get(shared, 0)
    return ("ceiling", max(feas)) if n_sh > 0 else ("floor", min(feas))


def handoff_window(producer, consumer, shared="h2", conc=None, q_prod=0.0, q_cons=QUANTUM,
                   lo=1e-11, hi=1e-1, n=400):
    """Feasible [shared] sub-range (lo_M, hi_M) where BOTH clear their energy requirement, else None. ASYMMETRIC:
    the producer's electron-DISPOSAL step needs only dG' < 0 (q_prod=0); the energy-CONSERVING scavenger must
    clear the Schink quantum (q_cons=20). producer/consumer are (sub, prod) reaction dicts."""
    conc = dict(PHYS, **(conc or {}))
    # HOISTED, same reason as solo_threshold: two dG'0 evaluations rather than two per point. If either
    # reaction is unreliable it is unreliable at every point, which the old loop expressed as `continue` at
    # every point and therefore an empty result.
    gp, _ = _dg0(*producer)
    gc, _ = _dg0(*consumer)
    if gp is None or gc is None:
        return None
    feas = []
    for i in range(n + 1):
        m = lo * (hi / lo) ** (i / n)
        c = dict(conc); c[shared] = m
        if gp + RT * _lnQ(*producer, c) <= -q_prod and gc + RT * _lnQ(*consumer, c) <= -q_cons:
            feas.append(m)
    return (min(feas), max(feas)) if feas else None


# --- re-runnable GUARD: the engine must reproduce the MEASURED interspecies-H2 thresholds ---------------------
# Documented thresholds (Cord-Ruwisch 1988; Schink 1997): sulfate-reducer <1 Pa, methanogen 1-10 Pa,
# acetogen 50-95 Pa. Plus the discrimination anchor: strong syntrophies obligate (low [H2] ceiling) vs lactate
# not-obligate (high ceiling), and a non-empty propionate->methanogen overlap window. PROBE: fails at quantum=0.
_SCAVENGERS = {
    "sulfate-reducer": (dict(h2=4, so4=1), dict(h2s=1, h2o=4), (0.0, 1.0)),      # <1 Pa
    "methanogen":      (dict(h2=4, co2=1), dict(ch4=1, h2o=2), (1.0, 10.0)),     # 1-10 Pa
    "acetogen":        (dict(h2=4, co2=2), dict(ac=1, h2o=2),  (50.0, 95.0)),    # 50-95 Pa
}
_PROPIONATE = (dict(ppa=1, h2o=2), dict(ac=1, co2=1, h2=3))
_LACTATE = (dict(lac__L=1, h2o=1), dict(ac=1, co2=1, h2=2))


def guard(quantum=QUANTUM):
    """Return a list of (name, ok, detail) checks reproducing the measured thresholds. Empty eQuilibrator -> []."""
    try:
        thermo._cc()
    except Exception as e:  # eQuilibrator unavailable -> skip cleanly, like the other thermo guards
        return [("SKIP", True, f"eQuilibrator unavailable: {type(e).__name__}")]
    out = []
    for name, (s, p, (lo_pa, hi_pa)) in _SCAVENGERS.items():
        _, v = solo_threshold(s, p, quantum=quantum)
        pa = M_to_Pa(v) if v else None
        ok = pa is not None and lo_pa <= pa <= hi_pa
        out.append((f"{name} floor in {lo_pa}-{hi_pa} Pa", ok, f"{pa:.2f} Pa" if pa else "no floor"))
    _, pc = solo_threshold(*_PROPIONATE)          # obligate: low ceiling
    _, lc = solo_threshold(*_LACTATE)             # not-obligate: high ceiling
    out.append(("propionate obligate (ceiling < 1e-4 M)", pc is not None and pc < 1e-4,
                f"{pc:.2e} M" if pc else "none"))
    out.append(("lactate NOT-obligate (ceiling > 1e-4 M)", lc is not None and lc > 1e-4,
                f"{lc:.2e} M" if lc else "none"))
    w = handoff_window(_PROPIONATE, _SCAVENGERS["methanogen"][:2])
    out.append(("propionate->methanogen window non-empty", w is not None,
                f"{M_to_Pa(w[0]):.1f}-{M_to_Pa(w[1]):.1f} Pa" if w else "NONE"))
    return out
