"""
ABIOTIC CHEMISTRY of the shared medium -- the reactions that happen WITHOUT an enzyme.

WHY THIS EXISTS: the inorganic and organic chemistry of the shared environment is easy to overlook, yet some
of it can react to produce a feedstock that was never supplied as an input. That gap was STRUCTURAL:
community._fba_components() keeps only `transfer.type == 'fba'` members, so every chemistry component
(carbonate_system, aeration_o2, ...) is silently dropped from every community LP. The platform's whole premise
is a netlist of heterogeneous boxes -- engine.py/Vivarium does stitch them -- but the community-FBA path has
always seen organisms only.

CHEMISTRY IS ENVIRONMENT, NOT A MEMBER. This is the same lesson the medium taught: the medium is resolved
once from the full member set and held FIXED while ORGANISMS are removed.
If the medium holds sulfide and Cu2+, precipitation HAPPENS -- it is not a member you select. So these
reactions are injected into the community's MEDIUM compartment when the abiotic layer is enabled
(MICROCOSM_ABIOTIC). Consequence that matters:
the tuple space stays C(organisms, n). Modelling chemistry as an optional member would multiply the search
space for nothing; modelling it as environment costs zero combinatorially.

It also enables a control that CANNOT EXIST without it: the ABIOTIC-ONLY BASELINE. If the chemistry alone
makes the target, the consortium is irrelevant and the hit is an artifact. super_additivity enumerates subsets
down to singletons but never to the empty set -- meaningless with organisms only, essential once the
environment can react on its own.

=== HONEST BOUNDARIES (each is a real limit, recorded not hidden) ===

1. FBA MODELS THIS CHEMISTRY AS *AVAILABLE*, NOT AS *SPONTANEOUS*. An LP carries flux through a reaction only
   if it serves the objective. Real precipitation happens whether or not it helps anything. So these reactions
   fire when the TARGET is downstream of them (the LP reaches for them) and sit idle otherwise. We may
   therefore UNDER-report chemistry, never over-report it -- the safe direction, but say so.
2. IRREVERSIBLE + STRONGLY FAVOURABLE ONLY. A reversible abiotic reaction in an LP is a free-energy loop
   waiting to happen (this is why loopless FBA exists). Every reaction here is thermodynamically one-way under
   the conditions of interest.
3. NO ACID/BASE EQUILIBRIA, NO pH. FBA sets FLUXES, not equilibrium ratios -- it would drive CO2/HCO3- to
   whatever the objective prefers rather than to the pKa 6.35 speciation, i.e. a free CO2 pump. A
   `carbonate_system` component is therefore deliberately NOT ported here. pH is out of reach of a
   stoichiometric method, exactly like regulation and kinetics. A "pH crash" toxicity concern remains
   unmodellable.
4. NO KINETICS. Sulfide-mediated azo reduction is real but slow and mediator-dependent; precipitation is fast.
   FBA cannot tell them apart. Rates are not claimed.
5. PRECIPITATES LEAVE THE SYSTEM. A solid is not a steady-state metabolite; each precipitate gets an explicit
   sink, or mass balance silently forces its flux to zero.
"""
from __future__ import annotations

# Each entry: (id, {medium_base: stoichiometric coeff}, note-with-citation)
# NEGATIVE coeff = consumed, POSITIVE = produced. Bases are BiGG medium bases (community compartment `_m`).
# Products suffixed `__precip` are SOLIDS: they get an automatic sink (see boundary 5).
ABIOTIC_REACTIONS = [
    # --- metal-sulfide precipitation: the actual product in acid-mine-drainage metal recovery. Biology
    # (sulfate reducers) makes the H2S; the recoverable solid (CuS/ZnS) is pure chemistry, which pure FBA
    # could not represent at all -- it could model the sulfide and never the recovery.
    # Ksp values are textbook (Stumm & Morgan, Aquatic Chemistry 3rd ed.): these are ~10^-25..10^-37, i.e.
    # quantitative and utterly one-way under any realistic condition.
    ("ABIO_CuS_precip", {"h2s": -1, "cu2": -1, "CuS__precip": 1, "h": 2},
     "H2S + Cu2+ -> CuS(s) + 2H+ ; Ksp(CuS) ~6e-37 (Stumm & Morgan). Quantitative, irreversible."),
    ("ABIO_ZnS_precip", {"h2s": -1, "zn2": -1, "ZnS__precip": 1, "h": 2},
     "H2S + Zn2+ -> ZnS(s) + 2H+ ; Ksp(ZnS) ~2e-25 (Stumm & Morgan). Quantitative, irreversible."),
    ("ABIO_FeS_precip", {"h2s": -1, "fe2": -1, "FeS__precip": 1, "h": 2},
     "H2S + Fe2+ -> FeS(s) + 2H+ ; Ksp(FeS) ~6e-19 (Stumm & Morgan). Irreversible; the classic AMD sink."),
    ("ABIO_CdS_precip", {"h2s": -1, "cd2": -1, "CdS__precip": 1, "h": 2},
     "H2S + Cd2+ -> CdS(s) + 2H+ ; Ksp(CdS) ~1e-27 (Stumm & Morgan)."),

    # --- abiotic redox by O2. Fe(II) oxidation is the reaction that REGENERATES Geobacter's electron acceptor
    # -- i.e. exactly the "chemistry produces a feedstock not present as an input" case.
    ("ABIO_Fe2_oxidation", {"fe2": -4, "o2": -1, "h": -4, "fe3": 4, "h2o": 2},
     "4Fe2+ + O2 + 4H+ -> 4Fe3+ + 2H2O ; abiotic ferrous oxidation (Stumm & Morgan). Fast at neutral pH."),
    ("ABIO_sulfide_oxidation", {"h2s": -2, "o2": -1, "s": 2, "h2o": 2},
     "2H2S + O2 -> 2S(0) + 2H2O ; abiotic sulfide oxidation (Stumm & Morgan). Irreversible."),

    # --- cyanide detox. CN- + S(0) -> SCN- is the documented abiotic route, and it MAKES THE FEEDSTOCK for
    # thiocyanate degradation out of a cyanide waste stream -- again chemistry producing a feedstock.
    ("ABIO_CN_to_SCN", {"cyan": -1, "s": -1, "tcynt": 1},
     "CN- + S(0) -> SCN- ; abiotic cyanide->thiocyanate (Luthy & Bruce 1979 ES&T 13:1481). One-way."),
]

PRECIPITATE_SUFFIX = "__precip"


def reachable_from(feed_bases, target, alias=None):
    """THE ABIOTIC-ONLY BASELINE, as a pure stoichiometric reachability question: can the CHEMISTRY ALONE make
    `target` from `feed_bases`, with no organism present? Returns (reachable, path) -- `path` is the list of
    abiotic reaction ids that get there, so the answer is auditable rather than a bare bool.

    THIS IS THE CONTROL THE MODULE DOCSTRING PROMISED: an abiotic-only baseline. Because ABIO_CN_to_SCN
    (CN- + S(0) -> SCN-) is injected into the medium and fires whether or not a cell is present, a screen that
    ignores this baseline can measure OUR OWN injected CHEMISTRY and report it as a biological capability -- and
    an identical score across many organisms is the tell-tale signature of that class of artifact.

    So a problem whose TARGET is abiotically reachable from its FEED is MIS-SPECIFIED -- it asks biology for
    something a bucket of reagents already does. Cyanide/thiocyanate is exactly that: the abiotic step makes
    SCN-, so the biology's real job is to DESTROY thiocyanate, not make it. Cheap enough (a fixed-point over ~7
    reactions, no LP) to run as a startup gate before any compute.

    NOTE the asymmetry, deliberately: this checks REACHABILITY, not yield. A target that chemistry can reach
    only partially is still suspect and still flagged; deciding it is fine is a human call, not a silent pass.
    """
    alias = alias or {}
    have = {alias.get(b, b) for b in feed_bases}
    tgt = alias.get(target, target)
    path = []
    for _ in range(len(ABIOTIC_REACTIONS) + 1):          # fixed point: each pass may unlock the next reaction
        grew = False
        for rid, stoich, _note in ABIOTIC_REACTIONS:
            if rid in path:
                continue
            need = [b for b, c in stoich.items() if c < 0]
            if all(b in have for b in need):
                for b, c in stoich.items():
                    if c > 0 and b not in have:
                        have.add(b)
                        grew = True
                path.append(rid)
        if not grew:
            break
    return (tgt in have), path


def _medium_met(com, base):
    """The community-compartment metabolite for a base, or None. micom names them `<base>_m`; AGORA2 drops the
    BiGG double underscore, so try both (the namespace split that silently starved every AGORA2 community --
    see community._ex_ids)."""
    for b in (base, base.replace("__", "_")):
        mid = f"{b}_m"
        if com.metabolites.has_id(mid):
            return com.metabolites.get_by_id(mid)
    return None


def applicable(com):
    """Which ABIOTIC_REACTIONS can actually be built in THIS community -- i.e. every non-precipitate species is
    present in the medium compartment. Returns [(id, stoich, note)]. Reported, never assumed: a reaction whose
    species are absent is silently impossible, and unreported silence is exactly the failure mode this guards against."""
    out = []
    for rid, stoich, note in ABIOTIC_REACTIONS:
        needed = [b for b in stoich if not b.endswith(PRECIPITATE_SUFFIX)]
        if all(_medium_met(com, b) is not None for b in needed):
            out.append((rid, stoich, note))
    return out


def inject(com):
    """Add the applicable abiotic reactions to the community's MEDIUM compartment, plus a sink per precipitate.
    Returns the list of reaction ids added (for logging/audit -- a silent modification of the model is exactly
    the kind of thing that turns into an unexplained result three days later).

    Injected into the community AFTER micom builds it, deliberately: these are not micom `members` (they have
    no biomass, and micom's growth machinery assumes one per member), they are the ENVIRONMENT the members
    share. This keeps the member machinery untouched.
    """
    from cobra import Metabolite, Reaction
    added = []
    for rid, stoich, note in applicable(com):
        if com.reactions.has_id(rid):
            continue
        rxn = Reaction(rid)
        rxn.name = note[:120]
        rxn.lower_bound = 0.0                 # IRREVERSIBLE -- see boundary 2 (no free-energy loops)
        rxn.upper_bound = 1000.0
        mets = {}
        for base, coeff in stoich.items():
            if base.endswith(PRECIPITATE_SUFFIX):
                m = Metabolite(f"{base}_m", name=base, compartment="m")   # the solid
                mets[m] = coeff
            else:
                mets[_medium_met(com, base)] = coeff
        rxn.add_metabolites(mets)
        com.add_reactions([rxn])
        added.append(rid)
        # every precipitate must LEAVE the system, or steady-state pins its flux to zero (boundary 5)
        for base, coeff in stoich.items():
            if base.endswith(PRECIPITATE_SUFFIX) and coeff > 0:
                # named EX_, not DM_: community_production resolves a target via _ex_ids() -> EX_<base>_m, so a
                # DM_-named sink would make the precipitate structurally unreachable AS A TARGET -- i.e. the
                # product the use-case exists to make would be invisible. Secretion-only (lb=0), so cobra's
                # medium setter cannot open it as a free input, and the defined-carbon closure ignores it
                # (it never appears in com.medium, which only lists exchanges with lb<0).
                sink_id = f"EX_{base}_m"
                if not com.reactions.has_id(sink_id):
                    sink = Reaction(sink_id)
                    sink.name = f"sink: {base} leaves as solid"
                    sink.lower_bound = 0.0
                    sink.upper_bound = 1000.0
                    sink.add_metabolites({com.metabolites.get_by_id(f"{base}_m"): -1})
                    com.add_reactions([sink])
                    added.append(sink_id)
    return added
