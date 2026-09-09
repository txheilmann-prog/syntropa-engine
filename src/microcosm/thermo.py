"""
THERMODYNAMIC FEASIBILITY GATE via eQuilibrator component-contribution (Noor et al. 2013; Beber et al.,
eQuilibrator 3.0, NAR 2022). The authoritative open dGf source -- no hand-rolled group contribution.

WHY: FBA and community FBA are stoichiometric and SECOND-LAW-BLIND -- they balance atoms and capacity, not
free energy -- so a transform that is "feasible" by yield can be a free-energy fantasy (a net-endergonic
transformation the community could not actually drive). This gate computes the NET transform's transformed
Gibbs energy dG'm and REJECTS what is CONFIDENTLY endergonic. It is a NECESSARY condition (net free energy
must be released), not a full pathway MDF -- internal driving-force bottlenecks are a deeper check (future).

Scope + honesty: eQuilibrator identifies compounds by BiGG metabolite id (verified: ChEBI accessions are
not indexed in the local cache, BiGG ids are), which is the native namespace of the community transform.
Species that do not map (non-chemical drivers like light/photons, or un-referenced metabolites) are
reported as `unmapped` and the verdict is "unknown" if a driver is missing -- we never fake a dG. The proton
is dropped (its contribution is absorbed into the pH-transformed potentials in the prime framework).

eQuilibrator stays BEHIND this seam (a sanctioned thermodynamics seam).
"""
from __future__ import annotations

_CC = None
_CC_FAILED = False
_COMPOUND_CACHE: dict = {}

# BiGG-id -> KEGG-id fallback for compounds whose BiGG entry in the local eQuilibrator cache is a STRUCTURELESS
# placeholder (inchi_key=None -> no reliable component-contribution dGf). This is not an approximation: each mapping points
# at the SAME chemical species in a namespace that DOES carry a structure/estimate, verified by inchi_key presence.
# Curated for the short-chain fatty acids + fermentation intermediates central to syntrophy, where BiGG silently
# fails (butyrate/propionate/isobutyrate/valerate...). Extend with verified mappings, never by guessing.
_KEGG_FALLBACK = {
    "but": "C00246", "ppa": "C00163", "isobut": "C02632", "ptrc": "C00134",
    "va": "C00803", "hxa": "C01585", "ibt": "C02632", "2mbut": "C18319", "3mbut": "C08262",
}


def _cc():
    """ComponentContribution singleton. The FIRST call may download the eQuilibrator compound cache (~1 GB);
    that fetch is BOUNDED -- a per-read socket-stall timeout plus a hard total cap, run in a daemon thread so
    the bound holds regardless of how the download library handles timeouts -- so a stalled or unreachable
    download RAISES instead of hanging forever. Callers (_compound, transform_dg, route_dg) already convert that
    into verdict 'unknown'/unmapped rather than blocking, and a failed init is remembered so later calls fail
    fast instead of re-hanging. Override the bounds with MICROCOSM_CC_TIMEOUT / MICROCOSM_CC_STALL (seconds)."""
    global _CC, _CC_FAILED
    if _CC is not None:
        return _CC
    if _CC_FAILED:
        raise RuntimeError("eQuilibrator ComponentContribution unavailable (a prior init failed or timed out)")
    import os
    import socket
    import threading
    from equilibrator_api import ComponentContribution, Q_
    total_s = float(os.environ.get("MICROCOSM_CC_TIMEOUT", "600"))   # hard cap on the first-run cache fetch
    stall_s = float(os.environ.get("MICROCOSM_CC_STALL", "90"))      # per-read socket-stall bound
    box: dict = {}

    def _build():
        prev = socket.getdefaulttimeout()
        socket.setdefaulttimeout(stall_s)              # a stalled SSL read raises instead of blocking forever
        try:
            box["cc"] = ComponentContribution()
        except Exception as e:                          # record ANY init failure so callers fail fast
            box["err"] = e
        finally:
            socket.setdefaulttimeout(prev)

    t = threading.Thread(target=_build, name="eq-cc-init", daemon=True)
    t.start()
    t.join(total_s)
    if "cc" not in box:                                 # still alive (timed out) or raised inside the thread
        _CC_FAILED = True
        if "err" in box:
            raise RuntimeError(f"eQuilibrator init failed: {type(box['err']).__name__}: {str(box['err'])[:80]}")
        raise TimeoutError(f"eQuilibrator cache init exceeded {total_s:.0f}s (stalled or unreachable download)")
    cc = box["cc"]
    cc.p_h = Q_(7.0)
    cc.ionic_strength = Q_("0.25M")
    cc.temperature = Q_("298.15K")
    _CC = cc
    return _CC


def _reliable(c):
    """A compound is usable for a dGf estimate only if it carries a structure (inchi_key). A BiGG placeholder
    (inchi_key=None) yields the eQuilibrator 1e5-uncertainty sentinel -> we treat it as UNMAPPED, never fake it."""
    return c is not None and getattr(c, "inchi_key", None) is not None


def _compound(base):
    """eQuilibrator compound for a BiGG metabolite base id, cached; None if unknown/non-chemical/structureless.
    Tries BiGG first; if that resolves to a STRUCTURELESS placeholder (no reliable dGf), falls back to the curated
    KEGG id for the SAME species. Returns None (reported as unmapped) rather than a compound we cannot estimate --
    the honesty rail: an unresolvable compound leaves its reaction UNCONSTRAINED, it is never assigned a fake dG."""
    if base in _COMPOUND_CACHE:
        return _COMPOUND_CACHE[base]
    c = None
    try:
        c = _cc().get_compound("bigg.metabolite:" + base)
    except Exception:
        c = None
    if not _reliable(c) and base in _KEGG_FALLBACK:
        try:
            alt = _cc().get_compound("kegg:" + _KEGG_FALLBACK[base])
            if _reliable(alt):
                c = alt
        except Exception:
            pass
    if not _reliable(c):
        c = None                               # structureless -> unmapped, not faked
    _COMPOUND_CACHE[base] = c
    return c


def transform_dg(transform, base_of):
    """Transformed Gibbs energy dG'm of a NET transform {species: net_flux} (produced>0, consumed<0).

    base_of: callable species -> BiGG metabolite base id (or None). Returns a dict:
      dg_prime_kj_mol, uncertainty_kj_mol, verdict (feasible|marginal|infeasible|unknown),
      thermo_ok (False only when CONFIDENTLY endergonic), unmapped [species].
    Only CONFIDENTLY endergonic transforms (dG'm - uncertainty > 0) are rejected -- we kill what is
    thermodynamically impossible, not what is merely uncertain."""
    from equilibrator_api import Reaction
    sparse, unmapped = {}, []
    for sp, flux in transform.items():
        if abs(flux) < 1e-6:
            continue
        base = base_of(sp)
        if base == "h":                       # proton: implicit in the prime (pH-transformed) framework
            continue
        c = _compound(base) if base else None
        if c is None:
            unmapped.append(sp)
            continue
        sparse[c] = sparse.get(c, 0.0) + flux
    sparse = {c: v for c, v in sparse.items() if abs(v) > 1e-9}
    if len(sparse) < 2:
        return {"dg_prime_kj_mol": None, "verdict": "unknown",
                "reason": "fewer than 2 mappable compounds", "unmapped": unmapped, "thermo_ok": True}
    if unmapped:                              # a missing driver (e.g. light) makes dG meaningless -- don't fake it
        return {"dg_prime_kj_mol": None, "verdict": "unknown",
                "reason": f"unmapped species ({','.join(unmapped[:4])})", "unmapped": unmapped, "thermo_ok": True}
    try:
        dg = _cc().physiological_dg_prime(Reaction(sparse))
        val = float(dg.value.to("kJ/mol").magnitude)
        err = float(dg.error.to("kJ/mol").magnitude)
    except Exception as e:
        return {"dg_prime_kj_mol": None, "verdict": "unknown",
                "reason": f"{type(e).__name__}: {str(e)[:60]}", "unmapped": unmapped, "thermo_ok": True}
    return _verdict(val, err, [])


# --- ELECTRON-CARRIER redox extension -------------------------------------------------------------------------------
# Component contribution cannot estimate PROTEIN-BOUND redox carriers (ferredoxin's Fe-S cluster has no small-molecule
# structure). eQuilibrator DOES already cover NAD(P), FAD, FMN, quinones -- verified -- so the ONLY missing electron
# carrier that matters for anaerobic syntrophy is FERREDOXIN, the low-potential carrier cells use to dispose reducing
# equivalents. We add it via its MEASURED standard redox potential E'0 (authoritative, cited), using eQuilibrator's
# OWN aqueous-H2 half-cell as the reference -- calibrated live from eQ's NAD couple so the two scales are consistent
# (this absorbs the H2 gas-vs-aqueous-1M standard-state offset exactly; validated to <1 kJ/mol against eQ-full NAD).
_F = 96.485                                    # Faraday constant, kJ/(mol*V)
_E_NAD = -0.320                                # E'0(NAD+/NADH), pH7 -- established constant (Thauer/Buckel; textbook)
_E_H2_AQ = None                                # eQ's effective aqueous-H2 half-cell potential (V), calibrated on demand
# each carrier: reduced/oxidized BiGG bases, electrons n, E'0 (V, pH7), citation. Fe2S2 ferredoxin is a 1-e- carrier
# (fdxrd charge +1 vs fdxo_2_2/fdxox charge +2, verified across genome-scale models).
CARRIER_POTENTIALS = [
    # E0_sd honestly spans the ferredoxin CLASS: plant 2Fe-2S ~ -0.42, bacterial/bifurcating [4Fe-4S] down to ~-0.50 V.
    # +/-0.05 V -> ~5 kJ/mol per electron of uncertainty, propagated into every ferredoxin reaction (never hidden).
    {"red": ("fdxrd",), "ox": ("fdxo_2_2", "fdxox"), "n": 1, "E0": -0.420, "E0_sd": 0.050,
     "cite": "2Fe-2S ferredoxin E'0 ~ -0.42 V, class spread -0.39..-0.50 V (Cammack 1992; Buckel & Thauer 2013, BBA 1827:94)"},
    # ModelSEED-namespace variant: a [8Fe-8S] (2x[4Fe-4S], clostridial-type) ferredoxin, distinct base ids
    # fdxr_42/fdxo_42 (formula Fe8S8X in the models). It is a 2-ELECTRON carrier: model charges are 0 (uninformative),
    # but the formula (two [4Fe-4S] clusters) AND the electron balance of its own hydrogenases (e.g. HYDFDN:
    # fdxr_42 + nadh -> fdxo_42 + 2 H2: 2 H2 = 4 e-, nadh = 2 e- => fdxr_42 must donate 2 e-) both require n=2. Same
    # ferredoxin-class E'0/uncertainty as above. Without this, the H2/formate-disposal reactions of clostridial-type
    # acetogens and syntrophs are left UNMAPPED.
    {"red": ("fdxr_42",), "ox": ("fdxo_42",), "n": 2, "E0": -0.420, "E0_sd": 0.050,
     "cite": "[8Fe-8S] 2x[4Fe-4S] 2-e- ferredoxin (ModelSEED fdx*_42, formula Fe8S8); E'0 ~ -0.42 V ferredoxin "
             "class (Cammack 1992; Buckel & Thauer 2013, BBA 1827:94); n=2 from formula + reaction electron balance"},
    # Coenzyme F420 (deazaflavin), the 2-electron hydride carrier central to methanogen H2/C1 metabolism. eQ covers
    # only the OXIDIZED form -> treat the whole couple as a custom carrier. BiGG (f420_2/f420h2_2) + ModelSEED (coF420/
    # coF420h) namespaces both listed. Well-characterized -> tighter uncertainty than ferredoxin.
    {"red": ("f420h2_2", "coF420h"), "ox": ("f420_2", "coF420"), "n": 2, "E0": -0.360, "E0_sd": 0.020,
     "cite": "coenzyme F420 2e- hydride carrier: E'0 standard midpoint -0.340 V (Jacobson & Walsh 1984, Biochemistry 23:979-988, "
             "doi 10.1021/bi00300a028 -- the primary source for F420's redox potential); -0.380 V physiological in "
             "hydrogenotrophic methanogens (10:1 ox:red). Central -0.360 +/- 0.020 V spans that standard-to-physiological "
             "range (within the reported E0_sd)."},
    # Methanophenazine (MP/MPH2), the membrane electron carrier of the Methanosarcinales (aceticlastic/methylotrophic
    # methanogens) -- a 2-hydroxyphenazine, 2-electron carrier. eQ cannot estimate it (polyisoprenoid phenazine). E'0
    # is well-determined and UNDISPUTED. (The CoB-S-S-CoM heterodisulfide is deliberately NOT added: its E'0 is disputed
    # in the literature -- -143 mV Tietze 2003 vs -281 mV Laird 2024 -- so hard-coding either would be a guess.)
    {"red": ("dhmphze",), "ox": ("mphze",), "n": 2, "E0": -0.165, "E0_sd": 0.006,
     "cite": "methanophenazine E'0 = -0.165 +/- 0.006 V, 2e- (Tietze et al. 2003, ChemBioChem 4:333)"},
]
_RED_TO_CARRIER = {b: c for c in CARRIER_POTENTIALS for b in c["red"]}
_OX_TO_CARRIER = {b: c for c in CARRIER_POTENTIALS for b in c["ox"]}


def _e_ref():
    """eQuilibrator's effective aqueous-H2 half-cell potential E(H2,aq) (V), so a carrier added via its textbook E'0
    lands on eQ's own scale. Calibrated from eQ's NAD couple: eQ(NAD+ + H2 -> NADH) = -2F(E_NAD - E_H2,aq)."""
    global _E_H2_AQ
    if _E_H2_AQ is None:
        d = reaction_dg0_prime({"nad": 1, "h2": 1}, {"nadh": 1}, _no_carriers=True)
        if not d["reliable"]:
            raise RuntimeError("cannot calibrate electron reference (eQ NAD couple unavailable)")
        _E_H2_AQ = _E_NAD - d["dg0_prime_kj_mol"] / (-2.0 * _F)
    return _E_H2_AQ


def _extract_carriers(reactants, products):
    """Split off electron-carrier species, returning (react2, prod2, correction_kJ, uncertainty_kJ, ok). The carrier's electrons are
    replaced with H2 (keeping the residual reaction eQ-mappable + balanced); correction accounts for the carrier's
    potential vs eQ's aqueous H2. ok=False if a carrier is not conserved (broken stoichiometry -> refuse, don't fake)."""
    react2 = {b: n for b, n in reactants.items()}
    prod2 = {b: n for b, n in products.items()}
    correction, var = 0.0, 0.0                 # var accumulates the carrier-potential uncertainty (variance, kJ^2)
    carriers = {id(c): c for b, c in {**_RED_TO_CARRIER, **_OX_TO_CARRIER}.items()
                if b in reactants or b in products}
    for c in carriers.values():
        s_red = sum(products.get(b, 0) for b in c["red"]) - sum(reactants.get(b, 0) for b in c["red"])
        s_ox = sum(products.get(b, 0) for b in c["ox"]) - sum(reactants.get(b, 0) for b in c["ox"])
        if abs(s_red + s_ox) > 1e-6:           # carrier must be conserved (red produced == ox consumed)
            return None, None, None, None, False
        a = s_red                              # net moles of carrier REDUCED going forward
        for b in c["red"] + c["ox"]:           # drop carrier species from the residual reaction
            react2.pop(b, None)
            prod2.pop(b, None)
        # electrons that reduced the carrier are supplied by H2 instead: add a*(n/2) H2 to the PRODUCT side
        h2_prod = a * c["n"] / 2.0
        if h2_prod >= 0:
            prod2["h2"] = prod2.get("h2", 0) + h2_prod
        else:
            react2["h2"] = react2.get("h2", 0) - h2_prod
        correction += a * (-c["n"] * _F * (c["E0"] - _e_ref()))
        var += (a * c["n"] * _F * c.get("E0_sd", 0.0)) ** 2   # E'0 uncertainty -> reaction dG uncertainty (honest)
    return react2, prod2, correction, var ** 0.5, True


def reaction_dg0_prime(reactants, products, _no_carriers=False):
    """STANDARD transformed Gibbs energy dG'^0 (1 M reference, pH 7, I=0.25 M, 298.15 K) of a reaction given as
    {bigg_base: stoich>0} reactant/product dicts. The clean building block for TFA (a per-reaction driving-force
    constraint needs dG'^0; the [c]-dependent term RT*sum(nu*ln c) is applied by the caller). Returns a dict:
      dg0_prime_kj_mol, uncertainty_kj_mol, unmapped [bases], reliable (bool).
    'reliable' is False (and dg0 None) if ANY species is unmappable/structureless or eQuilibrator returns its
    1e5 sentinel -- such a reaction is left thermodynamically UNCONSTRAINED and that is reported, never faked."""
    from equilibrator_api import Reaction
    correction, carrier_sd = 0.0, 0.0
    if not _no_carriers and (set(reactants) | set(products)) & (set(_RED_TO_CARRIER) | set(_OX_TO_CARRIER)):
        reactants, products, correction, carrier_sd, ok = _extract_carriers(reactants, products)
        if not ok:                             # a carrier was not conserved (broken stoichiometry) -> refuse
            return {"dg0_prime_kj_mol": None, "uncertainty_kj_mol": None, "unmapped": ["<unconserved-carrier>"],
                    "reliable": False}
    sparse, unmapped = {}, []
    for base, n in reactants.items():
        if base == "h" or abs(n) < 1e-9:
            continue
        c = _compound(base)
        (unmapped.append(base) if c is None else sparse.__setitem__(c, sparse.get(c, 0.0) - abs(n)))
    for base, n in products.items():
        if base == "h" or abs(n) < 1e-9:
            continue
        c = _compound(base)
        (unmapped.append(base) if c is None else sparse.__setitem__(c, sparse.get(c, 0.0) + abs(n)))
    sparse = {c: v for c, v in sparse.items() if abs(v) > 1e-9}
    if unmapped or len(sparse) < 2:
        return {"dg0_prime_kj_mol": None, "uncertainty_kj_mol": None, "unmapped": unmapped, "reliable": False}
    try:
        dg = _cc().standard_dg_prime(Reaction(sparse))
        val = float(dg.value.to("kJ/mol").magnitude)
        err = float(dg.error.to("kJ/mol").magnitude)
    except Exception as e:
        return {"dg0_prime_kj_mol": None, "uncertainty_kj_mol": None, "unmapped": unmapped,
                "reliable": False, "reason": f"{type(e).__name__}: {str(e)[:60]}"}
    if err > 1000:                             # eQuilibrator sentinel: a compound has no reliable dGf
        return {"dg0_prime_kj_mol": None, "uncertainty_kj_mol": None, "unmapped": [], "reliable": False}
    total_err = (err ** 2 + carrier_sd ** 2) ** 0.5   # combine eQ + carrier-potential uncertainty in quadrature
    return {"dg0_prime_kj_mol": round(val + correction, 2), "uncertainty_kj_mol": round(total_err, 2),
            "unmapped": [], "reliable": True}


def route_dg(reactants, products, base_of):
    """Thermodynamic feasibility of a DEFINED transformation A->B (the directed-search gate). reactants and
    products are {species: stoich>0} dicts. Same verdict semantics as transform_dg. This is the clean
    interface for gating a candidate route: reject a substrate->product claim that is confidently endergonic
    before it is ever ranked or taken to the bench."""
    transform = {}
    for sp, n in reactants.items():
        transform[sp] = transform.get(sp, 0.0) - abs(n)
    for sp, n in products.items():
        transform[sp] = transform.get(sp, 0.0) + abs(n)
    return transform_dg(transform, base_of)


def _verdict(val, err, unmapped):
    if err > 1000:                            # eQuilibrator sentinel (~1e5): a compound has NO reliable dGf
        return {"dg_prime_kj_mol": None, "verdict": "unknown",
                "reason": "contains a compound eQuilibrator cannot estimate",
                "unmapped": unmapped, "thermo_ok": True}
    if val - err > 0:
        verdict = "infeasible"                # confidently endergonic -> a free-energy fantasy, reject
    elif val + err < 0:
        verdict = "feasible"                  # confidently exergonic
    else:
        verdict = "marginal"                  # uncertainty genuinely straddles zero
    return {"dg_prime_kj_mol": round(val, 1), "uncertainty_kj_mol": round(err, 1),
            "verdict": verdict, "thermo_ok": verdict != "infeasible", "unmapped": []}
