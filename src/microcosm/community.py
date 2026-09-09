"""
RIGOROUS multi-organism community metabolism via COMMUNITY FBA (micom; Diener, Gibbons & Resendis-Antonio,
mSystems 2020, DOI 10.1128/mSystems.00606-19). The authoritative Python implementation -- no hand-rolling.

WHY THIS EXISTS: independent per-organism dFBA lets simultaneous boxes each size their uptake as if they
owned the whole shared pool, so N consumers over-draw a limiting substrate (the medium goes negative and
diverges) -- competition was DISCARDED, not modelled. Community FBA solves ALL organisms together in ONE
optimization where the shared-medium mass balance is a HARD constraint: they compete for the same pool
inside a single LP, so over-draw is impossible and coexistence + cross-feeding are resolved rigorously
(exactly why SteadyCom/MICOM exist -- independent FBA gets competition wrong).

SOLVER + UNIQUENESS: a plain community LP max-growth solution is a NON-UNIQUE vertex of a
degenerate optimum -- the reported cross-feeding flips with solver pivot order (a real reproducibility bug).
We therefore solve with PARSIMONIOUS FBA (pFBA; Lewis et al., Mol Syst Biol 2010): minimize total flux at
max community growth. That is LP-only (glpk), deterministic run-to-run (verified), and picks a defensible
representative (minimal enzyme investment). pFBA NARROWS but does not fully remove degeneracy -- FVA ranges
on a reported exchange are the honest full disclosure (a deeper check, future). micom's QP cooperative-
tradeoff (fairness de-bias of fast growers) needs a robust QP solver; osqp -- the only open one -- returns
garbage on GEM-scale QPs, so it is gated behind `fair=True` for when Gurobi/CPLEX is present.
LIMITATION (recorded, not hidden): community FBA is stoichiometric and THERMODYNAMICS-BLIND -- a route it
scores by yield can be second-law-infeasible. A thermodynamic gate (TMFA + loopless + MDF ranking) must sit
UPSTREAM of any route-optimality claim; until it does, treat yields as necessary-not-sufficient.

micom stays BEHIND this module: it is a sanctioned community-FBA simulation seam.
"""
from __future__ import annotations

import os
import re

# LP solver for the pFBA community solve. GLPK can wedge forever on some GEM-scale community LPs (observed as a
# hang with idle CPU during a long batch run). HiGHS (optlang 'hybrid' interface: HiGHS for LP,
# OSQP only for QP) is a robust, fast, deterministic LP solver that does NOT wedge. Because the fair=False path
# is pFBA (LP-only), it uses HiGHS -- NOT OSQP -- so the author's "osqp returns garbage on GEM-scale QPs" concern
# (which is about the fair=True QP cooperative-tradeoff) does NOT apply here. Env-overridable; glpk stays as a
# fallback. A per-solve timeout (below) is a second guard so nothing can wedge even if a solver misbehaves.
_SOLVER = os.environ.get("MICROCOSM_SOLVER", "hybrid")
# seconds; a longer solve ABORTS rather than hangs. The cap exists to stop a misbehaving solver WEDGING (the
# GLPK hang noted above), and 300s is generous for that: a real solve is ~0.03s (slim_optimize on a
# ~1,000-reaction community), so anything past 300s is pathological, not merely slow.
# WHY A GENEROUS CAP MATTERS: a timeout returns nan, and nan is INDISTINGUISHABLE from an infeasible LP at the
# return value -- so under load a too-tight cap makes the gate report "community cannot grow" for what is really
# a solver abort: infrastructure impersonating biology, and load-dependent, so it need not reproduce on a quiet
# box. A generous cap spends TIME (recoverable) instead of CORRECTNESS (not), and the solver_failed flag below
# still separates "could not answer" from "the answer is no".
_SOLVE_TIMEOUT = float(os.environ.get("MICROCOSM_SOLVE_TIMEOUT", "300"))

# Uptake rate for offered auxotrophy-support metabolites (amino acids / vitamins / nucleosides) -- i.e. every
# species offered via `available` that is NOT named in `feed`. NOMINAL + TUNABLE, calibrated against the
# validation set, NOT chosen by taste. It must be high enough that auxotrophic curated GEMs still
# grow, and low enough that the FEEDSTOCK -- not the amino-acid floor -- is the carbon source, or
# super-additivity becomes unreachable and no syntrophy can be detected. Applies to the directed-production
# path only; community_transform keeps its own simpler default-medium semantics.
GROWTH_FACTOR_RATE = float(os.environ.get("MICROCOSM_GROWTH_FACTOR_RATE", "0.1"))
# Max CARBON flux a single NON-WHITELIST organic requirement may contribute (an AGORA2 model gap-filled to
# "require" a polysaccharide/sugar/dipeptide). NOMINAL + TUNABLE. It exists because a fixed MOLAR rate on a
# large polymer dumps enormous carbon: for example a medium carrying amylopectin (5400 C),
# pullulan (7200 C) and arabinogalactan (3862 C) each at 0.1 molar = 1646 carbon vs a lactose feedstock's
# 120, so lactate ran to 64.97 (>the ~40 lactose can yield) and NOTHING was feedstock-dependent. Capping the
# CARBON (rate = MAX_REQ_CARBON_FLUX / carbons) forces the feedstock to dominate while keeping such an
# organism feasible on a trace. Small whitelisted growth factors (GROWTH_FACTORS) are exempt -- their carbon
# is already small and they are the intended auxotrophy support.
MAX_REQ_CARBON_FLUX = float(os.environ.get("MICROCOSM_MAX_REQ_CARBON_FLUX", "0.1"))
# Uptake ceiling for TRACE inorganic NUTRIENTS (N, P, S, K, Mg, metals -- see TRACE below) when they are
# offered but not part of the feedstock. NOMINAL, TUNABLE, calibrated against the validation set --
# not chosen by taste. It is deliberately NOT the model's own default: those run
# to 1000 and that unconstrained soup is what let a community import free sulfide and free copper and
# "produce" CuS=500.673 while its biology made 0.372. Generous enough that N/P/S never limit growth (a cell
# needs them stoichiometrically with carbon), bounded enough that nothing is free. Carbon-free by
# construction, so a generous rate cannot manufacture a cross-feeding signal.
TRACE_RATE = float(os.environ.get("MICROCOSM_TRACE_RATE", "10"))

# SINGLE-CARBON SOURCES. The defined-carbon closure only closed sources with >=2 carbons, so every C1 organic
# compound in a model's default medium stayed open -- and methanogens/acetogens/methylotrophs are C1 organisms.
# For example, a solo methanogen whose default medium leaves methanol (a direct methylotrophic-methanogenesis
# substrate) open makes methane with no partner and no feedstock: the closure shuts the sugars while leaving the
# methanogen's actual substrate wide open, so the consortium and the feedstock are both irrelevant to the answer
# -- the unit-of-analysis error where a free handout, not the BLEND + FEEDSTOCK, drives the system.
# `for` (formate) was ALSO explicitly exempt here and is removed for the same reason: formate is an organic C1
# substrate AND the currency of interspecies formate transfer -- gifting it from the medium short-circuits the
# very syntrophy we are trying to detect. Closing a community exchange only blocks IMPORT from outside; members
# still cross-feed each other through the shared medium compartment, so real syntrophy remains reachable.
# Anything genuinely offered as a feedstock (methanol/formate/syngas scenarios exist) is exempt via
# `offered_bases`, so this does not break those runs.
_NON_CATABOLIC_C = ("co2", "hco3", "h2", "h2s")   # inorganic carbon (carbon cycle) + non-carbon. NOT `for`/`meoh`.

# ABIOTIC CHEMISTRY of the shared medium (see microcosm.abiotic). OFF by default so no existing result shifts
# silently; callers that need it turn it on. Chemistry is ENVIRONMENT, not a member: it is injected into the
# community's medium compartment and is always on for that community, so the tuple space stays C(organisms,n).
ABIOTIC = os.environ.get("MICROCOSM_ABIOTIC", "0") not in ("0", "", "false", "False")
_MIN_CARBON_TO_CLOSE = 1                          # close C1 sources too; a >=2 threshold leaves every C1 open

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # the microcosm project root
LIBRARY = os.path.join(ROOT, "library")

# inorganic trace nutrients kept freely available so GEMs can actually grow (mirrors the dFBA "map the
# carbon/energy sources, leave trace minerals free" policy). Uptake of these is not the A->B question.
# NOTE: O2 and CO2 are DELIBERATELY NOT trace-open -- oxygen is a defining environmental/redox variable
# (leaving it open silently makes every community aerobic and kills fermentation/overflow routes), and CO2
# is a carbon-cycle variable. Both are set explicitly per scenario via `available`/`feed`.
TRACE = {"h2o", "h", "pi", "so4", "k", "na1", "cl", "ca2", "mg2", "fe2", "fe3", "mn2", "zn2", "cu2",
         "cobalt2", "mobd", "ni2", "sel", "slnt", "tungs", "cd2", "hg2", "nh4"}


def _base(rxn_or_met):
    """EX_glc__D_m -> glc__D ; glc__D_e -> glc__D ; AGORA/VMH EX_12ppd_S(e) -> 12ppd_S (strip EX_ prefix and
    the compartment tag in either BiGG '_e' or AGORA/VMH '(e)'/'[e]' notation, so the medium exchange name
    matches what micom builds for the community)."""
    s = re.sub(r"^EX_", "", rxn_or_met)
    s = re.sub(r"[\[(][a-z][\])]$", "", s)                  # AGORA/VMH: 12ppd_S(e) or _S[e] -> 12ppd_S
    return re.sub(r"_[a-z0-9]$", "", s)


# Feedstock scenario aliases -> BiGG metabolite base (the scenarios speak human names, the models speak bases).
_FEED_ALIAS = {"glucose": "glc__D", "lactate": "lac__L", "sulfate": "so4", "methanol": "meoh", "acetate": "ac",
               "hydrogen": "h2", "co": "co", "fructose": "fru", "xylose": "xyl__D", "ethanol": "etoh",
               "glycerol": "glyc", "succinate": "succ", "pyruvate": "pyr", "formate": "for",
               "cellobiose": "cellb", "sucrose": "sucr", "arabinose": "arab__L",
               "propionate": "ppa", "butyrate": "but", "valerate": "va",
               "caproate": "hxa", "hexanoate": "hxa",
               "benzoate": "bz", "4-hydroxybenzoate": "4hbz"}


def _ex_ids(base):
    """Community exchange ids for a base in BOTH namespace conventions -- BiGG double-underscore (glc__D, ala__L)
    AND AGORA2 single-underscore (glc_D, ala_L). AGORA2 metabolite ids drop the second underscore of the BiGG
    stereo suffix, so EX_glc__D_m never matches an AGORA2 community while EX_glc_D_m does. Trying both makes the
    medium/target resolution namespace-agnostic, so an AGORA2 community is not silently starved."""
    variants = {base, base.replace("__", "_")}
    return ["EX_" + b + "_m" for b in variants]


# GROWTH FACTORS: biosynthetic building blocks (amino acids, vitamins/cofactors, nucleosides/bases, polyamines)
# that curated GEMs are commonly auxotrophic for. These stay OPEN from the default medium even though they contain
# carbon -- they are NOT catabolic carbon/energy sources. Everything else multi-carbon in the default medium (sugars,
# fibers, organic acids) is CLOSED so the defined feedstock is the SOLE multi-carbon source -> cross-feeding
# emergence (a product built from the feedstock via a relay no single member completes) becomes DETECTABLE instead
# of masked by a rich diet where every member already makes the whole palette alone.
GROWTH_FACTORS = frozenset({
    "ala__L", "arg__L", "asn__L", "asp__L", "cys__L", "gln__L", "glu__L", "gly", "his__L", "ile__L", "leu__L",
    "lys__L", "met__L", "phe__L", "pro__L", "ser__L", "thr__L", "trp__L", "tyr__L", "val__L",
    "btn", "fol", "ribflv", "thm", "thmpp", "pydxn", "pydx5p", "pnto__R", "nac", "4abz", "4hbz", "cbl1", "cbl2",
    "cbi", "adocbl", "pheme", "sheme", "hemeO", "q8", "2dmmq8", "mqn8", "mqn7", "5mthf", "10fthf", "amet",
    "ade", "gua", "ura", "csn", "hxan", "xan", "ins", "adn", "gsn", "cytd", "uri", "duri", "thymd", "dad_2",
    "dcyt", "dgsn", "din", "xtsn", "spmd", "ptrc", "23camp", "chol", "inost",
})


def _carbon_count(formula):
    """Number of carbon atoms in a metabolite formula (C6H12O6 -> 6). Avoids matching Ca/Cl/Co/Cu/Cd (C + a
    lowercase letter). No formula -> 0 (treated as non-carbon)."""
    if not formula:
        return 0
    m = re.search(r"C(?![a-z])(\d*)", formula)
    if not m:
        return 0
    return int(m.group(1)) if m.group(1) else 1


def _exchange_formula(com, rxn_id):
    """The chemical formula of an exchange reaction's single metabolite (for carbon counting the requirement
    medium). Returns '' if unavailable -- treated as non-carbon by _carbon_count, i.e. left at the molar rate,
    which is the safe direction (a metabolite whose carbon we cannot read is not throttled, not over-thrown)."""
    try:
        mets = list(com.reactions.get_by_id(rxn_id).metabolites)
        return mets[0].formula if mets else ""
    except Exception:
        return ""


def _fba_components(components):
    """The genome-scale FBA members -- only these join the community LP. Others (kinetic/chemistry) do not."""
    return [c for c in components if c.transfer.get("type") == "fba"]


def _model_path(comp):
    from .component import resolve_model_path, _fetch_bigg_model
    p = comp.transfer.get("path")
    anchor = getattr(comp, "source_path", None) or os.path.join(LIBRARY, "targets", "x.json")
    rp = resolve_model_path(p, anchor)
    if p and not os.path.exists(rp):
        _fetch_bigg_model(getattr(comp, "provenance", None), rp)
    return rp


def _base_to_species(components):
    """Map a BiGG metabolite base id -> the microcosm shared-medium species, from the components' own
    port_map grounding (EX-reaction -> species). Lets the community transform speak microcosm's namespace."""
    b2s = {}
    for c in _fba_components(components):
        for rxn, sp in c.transfer.get("port_map", {}).items():
            b2s.setdefault(_base(rxn), sp)
    return b2s


def community_transform(components, available=None, vmax=10.0, feed=None, abundances=None, fair=False,
                        thermo=True, enzyme_budget=None):
    """Solve the joint community metabolism for the FBA members and return their net shared-medium exchange
    -- the community's TRANSFORM -- with competition and coexistence resolved inside one optimization.

    available: iterable of microcosm species offered as substrates (uptake allowed up to vmax); trace
               inorganics are always available. None -> the community's own default medium.
    Returns a dict: feasible, community_growth, members {id: growth}, transform {species: net_flux>0 produced
    / <0 consumed}, produced [species], consumed [species]. Returns feasible=False if it cannot grow."""
    import pandas as pd
    from micom import Community

    fba = _fba_components(components)
    if len(fba) < 1:
        return {"feasible": False, "reason": "no FBA members"}
    files = [_model_path(c) for c in fba]
    if enzyme_budget is not None:                      # GECKO-family: bounded proteome -> overflow/labor-division
        from . import enzyme
        files = [enzyme.write_constrained(f, enzyme_budget) for f in files]
    tax = pd.DataFrame({
        "id": [re.sub(r"[^a-zA-Z0-9]", "_", c.id) for c in fba],
        "file": files,
        "abundance": [1.0] * len(fba) if not abundances else [abundances.get(c.id, 1.0) for c in fba],
    })
    com = Community(tax, progress=False, solver=_SOLVER)
    # WEDGE GUARD: GLPK can hang FOREVER on some GEM-scale community pFBA LPs (observed as a freeze with idle
    # CPU). Cap every solve so no single LP can wedge a worker. HiGHS rarely needs it; GLPK
    # fallback does. Best-effort -- not every optlang interface exposes a timeout.
    try:
        com.solver.configuration.timeout = _SOLVE_TIMEOUT
    except Exception:
        pass
    if ABIOTIC:                                        # optional abiotic medium chemistry (precipitation, redox);
        from . import abiotic                          # off by default (MICROCOSM_ABIOTIC), so the demos are unaffected
        abiotic.inject(com)
    b2s = _base_to_species(components)
    s2base = {}
    for b, s in b2s.items():
        s2base.setdefault(s, b)

    if available is not None or feed is not None:
        feed = feed or {}                              # {species: per-substrate uptake rate} -- the FEED knob
        offered = set(available) if available is not None else set(feed)
        rate_by_base = {_FEED_ALIAS.get(s, s2base.get(s, s)): float(feed.get(s, vmax)) for s in offered}
        medium = {}
        # EXPLICITLY open each offered species' community exchange -- even if it is NOT in the default medium
        # (e.g. light/photon, which is closed by default), so an offered substrate is actually delivered.
        # Try BOTH namespace conventions (BiGG glc__D + AGORA2 glc_D) so AGORA2 communities are actually fed.
        for base, rate in rate_by_base.items():
            for rxn in _ex_ids(base):
                if com.reactions.has_id(rxn):
                    medium[rxn] = rate
                    break
        for r in list(com.medium):                     # trace nutrients stay generously open at default
            if _base(r) in TRACE:
                medium[r] = max(medium.get(r, 0.0), com.medium[r])
        if medium:
            com.medium = medium

    try:
        if fair:                                        # QP fairness refinement (needs a robust QP solver)
            sol = com.cooperative_tradeoff(fraction=0.5, fluxes=True)
        else:                                           # pFBA -> deterministic, reproducible representative
            sol = com.optimize(fluxes=True, pfba=True)
    except Exception as e:
        return {"feasible": False, "reason": f"{type(e).__name__}: {str(e)[:80]}"}
    # Same distinction the directed gate now draws: NO SOLUTION is "we could not answer" (the solver aborted
    # or returned nothing), which is NOT the biological claim "cannot grow on this medium". The except above
    # already names a raised failure; sol=None is the silent one, and it used to land in the same bucket as a
    # genuine zero-growth verdict.
    if sol is None:
        return {"feasible": False, "solver_failed": True, "community_growth": None,
                "reason": "solver returned no solution -- NOT a biological verdict"}
    if not (float(sol.growth_rate) > 1e-6):
        return {"feasible": False, "solver_failed": False,
                "reason": "community cannot grow on this medium",
                "community_growth": round(float(sol.growth_rate), 5)}

    fl = sol.fluxes
    transform = {}
    for col in fl.columns:
        if col.startswith("EX_"):
            v = float(fl.loc["medium", col])
            if abs(v) > 1e-4:
                sp = b2s.get(_base(col), _base(col))
                transform[sp] = transform.get(sp, 0.0) + v
    result = {
        "feasible": True,
        "community_growth": round(float(sol.growth_rate), 5),
        "members": {c.id: round(float(sol.members.growth_rate.get(
            re.sub(r"[^a-zA-Z0-9]", "_", c.id), 0.0)), 5) for c in fba},
        "transform": {k: round(v, 4) for k, v in transform.items()},
        "produced": sorted([k for k, v in transform.items() if v > 1e-3 and k not in ("H2O", "H_ext")]),
        "consumed": sorted([k for k, v in transform.items() if v < -1e-3]),
    }
    if thermo:                                # is the net transform thermodynamically REAL, not just balanced?
        from . import thermo as _t
        result["thermo"] = _t.transform_dg(transform, lambda sp: s2base.get(sp))
    return result


# Optional directed-production and screening solvers build on this substrate. They are an optional add-on;
# re-exported here so a caller has one import path whether or not they are installed.
try:  # noqa: SIM105
    from .directed import *  # noqa: F401,F403  (its __all__ lists the helpers re-exported here; absent in the open build)
except ImportError:
    pass
