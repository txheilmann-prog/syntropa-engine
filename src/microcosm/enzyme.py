"""
ENZYME-CONSTRAINED metabolism -- the sMOMENT lumped protein-pool constraint (Bekiaris & Klamt, BMC
Bioinformatics 2020; the lightweight form of the GECKO/MOMENT family). Plain FBA gives every reaction an
unlimited, free catalytic rate, so it has no reason to PREFER overflow metabolism / division-of-labor /
syntrophy -- the very emergent behaviors of interest in a multi-organism community. A finite
proteome budget forces the resource trade-off that makes those behaviors optimal -- the single biggest lever
for surfacing them. Verified: it reproduces AEROBIC ACETATE OVERFLOW (the Crabtree effect) that
plain FBA structurally cannot, and those overflow products are exactly the cross-feeding substrates.

Construction: split each reversible internal reaction so both directions draw enzyme, add a `prot_pool`
pseudo-metabolite every internal reaction consumes proportionally to (MW / kcat), and a bounded pool-supply
reaction -> sum((MW_j/kcat_j)|v_j|) <= budget, a hard catalytic ceiling. The constraint is INTRACELLULAR,
so enzyme-constrained models still compose in a community (each keeps its own proteome pool).

NO-FUDGE NOTE: kcat and MW here are a UNIFORM NOMINAL first cut (a documented catalytic efficiency,
order-of-magnitude per Bar-Even et al. 2011 "The moderately efficient enzyme"), NOT per-reaction measured
kcats. That is enough to introduce the bounded-proteome trade-off and reproduce overflow, and it is LABELED
nominal. The authoritative upgrade is per-reaction kcats from DLKcat/BRENDA -- a documented next step, not a
hidden approximation. cobra stays behind this seam (a sanctioned model-transformation seam).
"""
import hashlib
import os
import threading

import cobra

NOMINAL_KCAT_PER_S = 65.0    # nominal uniform turnover (1/s); central-carbon enzymes ~ tens-hundreds /s
NOMINAL_MW_KDA = 40.0        # nominal average enzyme molecular weight (kDa = g/mmol)

# THE BUDGET EVERY CALLER ACTUALLY USES, and it is not `enzyme_constrain`'s own default. The pilots, the use
# case search and every evidence harness pass 0.1 g protein / gDW; `enzyme_constrain(budget=0.04)` below is the
# signature default and NOTHING in the tree passes it. Keeping the live value as a bare literal at dozens of call
# sites is what let the two drift apart unnoticed, so it lives here and callers import it.
#
# THIS IS A CHOSEN SETTING, NOT A MEASURED PARAMETER, and it is labeled that way for the same reason the nominal
# kcat above is. Neither 0.1 nor 0.04 is derived from a measurement or a citation; the change that introduced 0.1
# named the mechanism and not the magnitude.
#
# WHAT IS MEASURED ABOUT IT, because a setting without a reason can still have a bound. Sweeping the budget over
# 0.04, 0.1, 0.2 and 0.4 on the methanogen syntrophy use case:
#   - EMERGENCE VERDICTS are invariant across the whole range, so a binary discovery claim does not rest on it.
#   - ABSOLUTE PRODUCTION tracks the budget below 0.1 and saturates above it, unevenly across candidates, so a
#     flux computed at one budget CANNOT be converted to another by arithmetic. Re-run instead.
#   - ABOVE ROUGHLY 0.1 THE POOL GOES SLACK AND THE CONSTRAINT STOPS BINDING. Six of eight abundance ratios give
#     the identical unblocked optimum at 0.2 and at 0.4, and a route dependence that measures about 15 percent at
#     0.1 measures exactly 0 at 0.4: block the reaction and the model pays for the alternative out of a budget it
#     was not spending. At that point the model is back to the unlimited free catalysis this module exists to
#     remove, so budgets much above 0.1 defeat the purpose of enzyme-constraining at all.
# So any claim whose magnitude depends on the budget must be quoted WITH the budget, and 0.1 is inside the regime
# where the constraint is active while 0.2 is already leaving it.
#
# THE AUTHORITATIVE UPGRADE is the same one the kcat note names: a budget derived from a measured proteome mass
# fraction and saturation for these organisms, which turns the setting into a parameter.
#
# WHICH CALL SITES IMPORT THIS, and why the rest deliberately do NOT. The ones that convert are the ones where a
# single source of truth is what you want: `outcome.py`'s two signature defaults (the only place under `src/`
# that baked a VALUE rather than taking one, every other module correctly defaulting to None) and the use case
# search's scorers. The ones that KEEP their literal are TESTS and RECORDED-EVIDENCE HARNESSES, on purpose:
# `tests/test_carbon_leak.py` and `tests/test_feedstock_dependence.py` PIN behavior at 0.1, so importing the
# constant would silently retarget what they assert the moment the constant moved, converting an assertion into a
# follower. Evidence harnesses are the same case, since their recorded numbers belong to the value that produced
# them. `tests/test_enzyme_budget_home.py` enforces exactly this split, both directions: no baked budget VALUE
# anywhere under `src/`, and the pinning tests must KEEP their literal, so a well-meant cleanup cannot quietly
# turn an assertion into a follower.
PILOT_BUDGET = 0.1           # g protein / gDW -- the budget the engine and every harness actually run

# On-disk cache of enzyme-constrained models. It persists across runs (avoids re-constraining the same member
# on every evaluation), but MUST be bounded with an eviction policy (below): keyed with no eviction, a cache
# like this grows without limit and can fill the system drive.
_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # microcosm project root
CACHE_DIR = os.environ.get("MICROCOSM_EC_CACHE_DIR", os.path.join(_ROOT, ".cache", "enzyme_constrained"))
# The cache is size-capped (default 10 GB) via the MICROCOSM_EC_CACHE_CAP_GB environment variable; LRU eviction
# (below) keeps it under the cap. 10 GB holds ~2,400 models warm (a ~1.52s cache hit vs a ~5.62s rebuild on a
# miss); a larger working set misses more often, but a full disk breaks everything. Raise the cap when more
# disk is available.
CACHE_CAP_BYTES = float(os.environ.get("MICROCOSM_EC_CACHE_CAP_GB", "10")) * (1024 ** 3)


def _evict_lru(cache_dir, cap_bytes):
    """Delete oldest-mtime cache entries until the directory is back under cap_bytes."""
    entries = [os.path.join(cache_dir, f) for f in os.listdir(cache_dir)]
    stats = [(p, os.path.getmtime(p), os.path.getsize(p)) for p in entries if os.path.isfile(p)]
    total = sum(size for _, _, size in stats)
    if total <= cap_bytes:
        return
    for path, _, size in sorted(stats, key=lambda s: s[1]):
        if total <= cap_bytes:
            break
        try:
            os.remove(path)
            total -= size
        except OSError:
            pass


def _coef(kcat_per_s, mw_kda):
    return mw_kda / (kcat_per_s * 3600.0)    # g_protein * h / mmol  == MW(g/mmol) / kcat(1/h)


def enzyme_constrain(model, budget=0.04, kcat_per_s=NOMINAL_KCAT_PER_S, mw_kda=NOMINAL_MW_KDA):
    """Return a copy of `model` with a sMOMENT lumped protein-pool constraint. `budget` = total enzyme
    budget (g protein / gDW); lower = tighter catalytic ceiling -> earlier overflow. Exchanges, biomass,
    and ATP maintenance draw no enzyme."""
    m = model.copy()
    coef = _coef(kcat_per_s, mw_kda)
    pool = cobra.Metabolite("prot_pool_c", name="protein pool (sMOMENT, nominal kcat)", compartment="c")
    supply = cobra.Reaction("PROT_POOL_SUPPLY")
    supply.bounds = (0.0, budget)
    supply.add_metabolites({pool: 1.0})
    m.add_reactions([supply])
    reverses = []
    for r in list(m.reactions):
        if r is supply or r.boundary or r.id.startswith("EX_") or "biomass" in r.id.lower() \
                or "ATPM" in r.id.upper():
            continue
        if r.lower_bound < 0.0:                          # split reversible: both directions cost enzyme
            rev = cobra.Reaction(r.id + "_ecrev")
            rev.add_metabolites({met: -c for met, c in r.metabolites.items()})
            rev.bounds = (0.0, -r.lower_bound)
            rev.add_metabolites({pool: -coef})
            reverses.append(rev)
            r.lower_bound = 0.0
        r.add_metabolites({pool: -coef})
    m.add_reactions(reverses)
    return m


def write_constrained(src_path, budget=0.04, kcat_per_s=NOMINAL_KCAT_PER_S, mw_kda=NOMINAL_MW_KDA):
    """Enzyme-constrain the SBML at src_path and write it to a cached temp SBML; return that path (so the
    community solver can compose enzyme-constrained members). Cached by (path, budget, kcat, mw) so the same
    member is not re-constrained on every evaluation."""
    # KEY ON THE MODEL'S CONTENT, NOT ITS PATH. Keying on abspath(src_path) makes this cache STRUCTURALLY
    # INCAPABLE OF EVER HITTING ACROSS PROCESSES on Windows when the loader materializes each model's payload to
    # a PER-PROCESS temp directory: the same model then has a different "path" -- and therefore a different key
    # -- in every worker, so each worker writes its own duplicate and the cache fills with per-PID copies of the
    # same few thousand models (which is what makes a path-keyed cache "thrash": not too small, just full of
    # garbage).
    # A content hash is stable across processes and across machines, and it is SAFER than a path key: if a
    # model's payload is ever re-generated or corrected, the key changes and the stale constrained model cannot
    # be silently reused. Cost is ~10-20ms of md5 on a ~4 MB SBML against a ~4 s rebuild -- noise.
    h = hashlib.md5()
    with open(src_path, "rb") as _f:
        for chunk in iter(lambda: _f.read(1 << 20), b""):
            h.update(chunk)
    key = hashlib.md5(f"{h.hexdigest()}|{budget}|{kcat_per_s}|{mw_kda}".encode()).hexdigest()[:12]
    os.makedirs(CACHE_DIR, exist_ok=True)
    out = os.path.join(CACHE_DIR, f"microcosm_ec_{key}.xml")
    if os.path.exists(out):
        os.utime(out, None)  # bump mtime so a hot entry survives LRU eviction
        return out
    cobra.io.sbml.LOGGER.disabled = True
    m = enzyme_constrain(cobra.io.read_sbml_model(src_path), budget, kcat_per_s, mw_kda)
    # ATOMIC WRITE (concurrency fix): parallel workers share this cache path (keyed by model+budget only, in
    # CACHE_DIR). Writing straight to `out` let a second worker see os.path.exists(out)==True while the file was
    # still half-written, then hand micom a TRUNCATED SBML that libSBML hangs parsing forever (a wedge: idle
    # CPU, solver-independent, only under concurrency).
    # Write to a per-writer-unique temp, then os.replace() -> `out` is published atomically and only ever
    # observed COMPLETE. Concurrent writers each rename their own temp; last wins, identical content.
    tmp = f"{out}.{os.getpid()}.{threading.get_ident()}.tmp"
    cobra.io.write_sbml_model(m, tmp)
    os.replace(tmp, out)
    _evict_lru(CACHE_DIR, CACHE_CAP_BYTES)
    return out
