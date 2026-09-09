"""
EMERGENCE DISCRIMINATOR -- what separates a genuine multi-organism capability from a trivial recombination
("organism B makes its known product once fed by A, and would make it anyway"). Without a COMPUTABLE
definition of novelty, every candidate is a curation opinion and results drown in rediscoveries. This makes
emergence measurable, so an emergent claim is defensible.

SUPER-ADDITIVITY (the labor-division / non-additivity signature): a community is super-additive for a
target product if it produces STRICTLY MORE of the target than the best of ALL its proper sub-communities
(every smaller subset, including singletons). If the whole exceeds every part, the interaction is
load-bearing -- a real emergent capability -- rather than a passenger. This is a NECESSARY condition for
calling a route a multi-organism discovery; obligate-intermediate (SMETANA-style) and mechanistic-distance
tests are complementary follow-ons.
"""
from __future__ import annotations

import itertools

from .community import community_transform
from .search import target_value


def _production(members, target, available, feed, vmax, enzyme_budget, consume, directed=False):
    """Return (target production, feasible). feasible=False means the community LP could not support all
    members (they cannot coexist on this medium) -- distinct from feasible-but-zero (they coexist but cannot
    make the target). The irreducibility metric must not conflate these.

    directed=False (default): SECRETION at the growth optimum, via community_transform. This is the path the
    ground-truth harness has always validated.
    directed=True: MAX EXPORT subject to a growth floor, via community_production. Both modes route target
    production through THIS one chokepoint so they stay comparable on the same ground truth, with the same
    significance, dynamic-survival, and redox-gate semantics rather than a divergent inline re-implementation.
    Nothing about the subset-isolation semantics below changes between the two.

    A proper subset is ALWAYS solved with the other members physically REMOVED from the LP -- so {A,B} here is
    A and B in ISOLATION (any C absent, none of C's fuel/detox/competition/redox/thermodynamic effect on A or
    B present). That is deliberately the correct baseline for "what can A,B do WITHOUT C". NOTE: do NOT
    memoize sub-community solves across candidates as a scaling shortcut. It is TECHNICALLY safe for this exact
    isolation semantics, but the platform's whole premise is that if C modifies A or B in ANY way the trio
    {A,B,C} is UNIQUE and not decomposable -- a cache that even looks like it assumes decomposability is the
    wrong instinct here, and exhaustive is cheap at current scale. Re-introduce only under real scale pressure,
    with an explicit feed-/modify-through-C test."""
    try:
        if directed and consume:
            # DESTRUCTION USE-CASES measure UPTAKE, not export. This branch used to not exist ("consume-mode is
            # a transform-only concept"), which forced every destruction problem to be scored as a production
            # problem -- and that is the wrong question. "Methane capture" is not CO2 output (every organism
            # makes CO2 from anything, so the metric barely moves); it is METHANE CONSUMED. Same for
            # thiocyanate destruction scored on sulfate, and phenol destruction scored on CO2.
            # community_screen reports per-species uptake off the SAME gate, so super-additivity on a
            # destruction target now asks the right question: does the consortium DESTROY more than any subset?
            # It also makes feedstock-dependence STRUCTURAL rather than a second solve -- uptake of a species
            # that is not offered is 0 by construction.
            from .community import community_screen
            r = community_screen(members, target, [target], available=available, feed=feed, vmax=vmax,
                                 enzyme_budget=enzyme_budget)
            if r.get("feasible"):
                return float((r.get("uptake") or {}).get(target, 0.0)), True
            if r.get("solver_failed"):
                return 0.0, False
            return 0.0, False
        if directed:
            from .community import community_production
            r = community_production(members, target, available=available, feed=feed, vmax=vmax,
                                     enzyme_budget=enzyme_budget)
            if r.get("feasible"):
                return float(r.get("production") or 0.0), True
            # A community that GREW but has no exchange for the target COEXISTS FINE -- it simply cannot make
            # the product. community_production only reports "no community exchange" AFTER the growth gate has
            # passed, so this is capability-gated: the STRONGEST basis, not the weakest. Collapsing it to
            # feasible=False made super_additivity report E.coli+methanogen as "all subsets cannot COEXIST
            # (medium-sensitive -- WEAK)" when in truth E. coli grows perfectly and structurally cannot make
            # methane. That is the exact distinction the basis logic below exists to draw, inverted.
            if "no community exchange" in str(r.get("reason", "")):
                return 0.0, True
            return 0.0, False
        r = community_transform(members, available=available, feed=feed, vmax=vmax, thermo=False,
                                enzyme_budget=enzyme_budget)
        return (target_value(r, target, consume=consume), True) if r.get("feasible") else (0.0, False)
    except Exception:
        return 0.0, False


def super_additivity(members, target, available=None, feed=None, vmax=10.0, max_enum=5, margin=0.02,
                     enzyme_budget=None, consume=False, directed=False):
    """Community target production minus the BEST proper sub-community's. synergy>0 (beyond margin) and a
    non-trivial community amount => emergent. For <=max_enum members every proper subset is enumerated;
    larger sets fall back to singletons + pairs (logged as partial). Returns the full comparison."""
    if directed:                                       # an optional directed-production path, used when present
        try:
            from . import directed as _directed_layer  # noqa: F401
        except ImportError:
            raise NotImplementedError(
                "directed=True is not available in this build; use the default super-additivity test (directed=False).")
    # Resolve the medium ONCE from the FULL member set, then hold it FIXED across every subset below. If it
    # were left None, each subset would re-derive its own medium and removing an organism would also remove
    # that organism's nutrients -- confounding "is this member load-bearing?" with "did the diet change?".
    # The medium is the ENVIRONMENT; the experiment only removes ORGANISMS from it.
    if available is None and directed:
        from .community import _derive_available
        available = _derive_available(members) or set()
    comm, _ = _production(members, target, available, feed, vmax, enzyme_budget, consume, directed)
    ids = list(range(len(members)))
    if len(members) <= max_enum:
        subsets = [c for k in range(1, len(members)) for c in itertools.combinations(ids, k)]
        exhaustive = True
    else:
        subsets = [(i,) for i in ids] + list(itertools.combinations(ids, 2))
        exhaustive = False
    best, best_members = 0.0, None                          # best production over ALL proper subsets
    best_feasible = 0.0                                     # best production over FEASIBLE proper subsets only
    n_sub = len(subsets); n_feasible = 0
    for idx in subsets:
        p, feas = _production([members[i] for i in idx], target, available, feed, vmax, enzyme_budget,
                              consume, directed)
        if feas:
            n_feasible += 1
            best_feasible = max(best_feasible, p)
        if p > best:
            best, best_members = p, [members[i].id for i in idx]
    synergy = comm - best
    single = len(members) < 2                              # a lone organism cannot be super-additive
    # HONESTY on WHY the subsets fail: is the irreducibility because smaller sets cannot
    # COEXIST (infeasible -- a coexistence-gated claim, medium-sensitive) or because they coexist but cannot
    # PRODUCE the target (capability-gated -- the stronger claim)? Distinguish them.
    n_infeasible = n_sub - n_feasible
    if single:
        basis = "single-member"
    elif best_feasible > 1e-3:
        basis = "reducible (a feasible subset already produces)"
    elif n_feasible == 0:
        basis = f"all {n_sub} subsets cannot COEXIST (coexistence-gated, medium-sensitive -- WEAK)"
    elif n_infeasible > n_feasible:
        basis = (f"mostly coexistence-gated (WEAK-MED): {n_infeasible}/{n_sub} subsets cannot coexist, "
                 f"{n_feasible} coexist but produce zero -- medium-sensitive")
    else:
        basis = (f"mostly capability-gated (STRONG): {n_feasible}/{n_sub} subsets coexist but cannot produce "
                 f"the target, {n_infeasible} cannot coexist")
    return {
        "target": target,
        "community": round(comm, 4),
        "best_subcommunity": round(best, 4),
        "best_subcommunity_members": best_members,
        "synergy": 0.0 if single else round(synergy, 4),
        "irreducibility": 0.0 if single else (round(synergy / comm, 4) if comm > 1e-9 else 0.0),
        "best_feasible_subcommunity": round(best_feasible, 4),
        "subsets_feasible": f"{n_feasible}/{n_sub}",
        "irreducibility_basis": basis,
        "emergent": (not single) and comm > best * (1.0 + margin) and comm > 1e-3,
        "produces": comm > 1e-3,                            # a direct (possibly single-organism) producer
        "subset_search": "exhaustive" if exhaustive else "singletons+pairs (partial)",
    }
