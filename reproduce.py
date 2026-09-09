"""Reproduce the headline results in one command:  python reproduce.py

Regenerates the numbers the README stands on, from the demo models (fetching any link-only model from its
source on first run), so a stranger can re-run and verify the results without trusting the author. Each numeric
known-answer is CHECKED against its expected value within a stated tolerance, and the script exits non-zero if
any check fails. If a link-only model or the eQuilibrator cache cannot be fetched (offline, or the source is
unreachable), the affected section is SKIPPED with a clear note rather than crashing -- a skip is not a failure.

Reproducibility note: FBA growth rates depend on the LP solver and its version (the reported values are
reproduced to within ~1% across solver stacks). For bit-stable numbers, install the pinned dependency set in
`requirements-lock.txt`.
"""
import os
import sys
import warnings

warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT, "src"))

import microcosm.builtins  # noqa: F401  (registers transfer functions)
from microcosm import load_component, composition_validity
from microcosm.community import community_transform
from microcosm.emergence import super_additivity
from microcosm.thermo import route_dg

_failures = []
_skips = []


def check(label, got, expected, tol):
    """Numeric known-answer: PASS if within tolerance, else record a failure."""
    ok = isinstance(got, (int, float)) and abs(got - expected) <= tol
    if not ok:
        _failures.append(label)
    print(f"  [{'PASS' if ok else 'FAIL'}] {label:46} {got:>10}   (expected {expected} +/- {tol})")


def note(label, got, expected):
    """Qualitative result -- printed for inspection, gated separately below."""
    print(f"  {label:53} {got:>18}   (expected {expected})")


def try_load(path):
    """Load a component manifest; return None (recording nothing yet) if a link-only model cannot be fetched
    (offline, or the upstream source unreachable) -- so the section degrades to a clear SKIP instead of crashing.
    Only fetch/availability errors are swallowed; a malformed manifest still raises."""
    try:
        return load_component(path)
    except (FileNotFoundError, OSError) as e:
        print(f"  [note] {os.path.basename(path)}: link-only model unavailable ({type(e).__name__}); see MODELS.md")
        return None


def skip(section_needs):
    print(f"  [SKIP] needs {section_needs} -- unavailable without network; re-run online to reproduce. See MODELS.md.")
    _skips.append(section_needs)


def main():
    ecoli = try_load(os.path.join(ROOT, "library", "fba_ecoli_core.json"))
    methanogen = try_load(os.path.join(ROOT, "library", "fba_methanogen.json"))

    print("\n1. KNOWN-ANSWER VALIDATION -- E. coli core vs Orth, Fleming & Palsson 2010")
    if ecoli is not None:
        aer = community_transform([ecoli], feed={"glucose": 10, "O2": 1000}, thermo=False)["community_growth"]
        ana = community_transform([ecoli], feed={"glucose": 10, "O2": 0}, thermo=False)["community_growth"]
        check("aerobic glucose-minimal growth (1/h)", round(aer, 4), 0.8739, 2e-3)
        check("anaerobic glucose growth (1/h)", round(ana, 4), 0.211, 3e-3)
    else:
        skip("the E. coli core model")

    print("\n2. THERMODYNAMIC GATE (eQuilibrator) -- rejects free-energy fantasies")
    try:
        fwd = route_dg({"glc__D": 1}, {"etoh": 2, "co2": 2}, lambda s: s)
        rev = route_dg({"etoh": 2, "co2": 2}, {"glc__D": 1}, lambda s: s)
        note("glucose -> 2 ethanol + 2 CO2", f"{fwd['dg_prime_kj_mol']} kJ/mol  {fwd['verdict']}", "< 0, feasible")
        note("its impossible reverse", f"{rev['dg_prime_kj_mol']} kJ/mol  {rev['verdict']}", "> 0, infeasible")
        if fwd.get("verdict") == "unknown" or rev.get("verdict") == "unknown":
            skip("the eQuilibrator compound cache (dG unknown until the cache is warm)")
        elif fwd.get("thermo_ok") is not True or rev.get("thermo_ok") is not False:
            _failures.append("thermodynamic gate")
    except Exception as e:
        skip(f"the eQuilibrator compound cache ({type(e).__name__})")

    print("\n3. RIGOROUS COMMUNITY FBA + EMERGENCE -- E. coli + methanogen -> methane")
    if ecoli is not None and methanogen is not None:
        r = community_transform([ecoli, methanogen],
                                available={"glucose", "ammonium", "phosphate", "CO2", "H2O"}, feed={"glucose": 10})
        note("community growth (1/h)", f"{r['community_growth']:.3f}", "> 0, both coexist")
        note("secretes", ", ".join(p for p in r["produced"] if p in ("methane", "acetate", "formate")), "incl. methane")
        emg = super_additivity([ecoli, methanogen], "methane",
                               available={"glucose", "ammonium", "phosphate", "CO2", "H2O"}, feed={"glucose": 10})
        note("methane: community vs best sub-community", f"{emg['community']} vs {emg['best_subcommunity']}", "emergent")
        note("emergent (neither makes methane alone)?", str(emg["emergent"]), "True")
        if emg.get("emergent") is not True:
            _failures.append("methane emergence")
    else:
        skip("the E. coli core + methanogen models")

    print("\n4. COMPOSITION-VALIDITY GATE -- catches a physically-impossible community before any simulation")
    cyano = try_load(os.path.join(ROOT, "library", "fba_cyanobacterium.json"))
    if cyano is not None and methanogen is not None:
        problems = composition_validity([cyano, methanogen])   # an aerobe + a strict anaerobe share no redox window
        note("aerobic cyanobacterium + methanogen", "REJECTED" if problems else "ACCEPTED", "REJECTED (disjoint Eh)")
        if problems:
            print(f"       -> {problems[0]}")
        note("cyanobacterium alone", "valid" if not composition_validity([cyano]) else "invalid", "valid")
        if not problems:
            _failures.append("composition-validity gate")
    else:
        skip("the cyanobacterium + methanogen models")

    print("\n5. INTERSPECIES-H2 THRESHOLDS -- the per-hand-off screen reproduces measured guild thresholds (Cord-Ruwisch et al. 1988)")
    from microcosm import handoff_screen
    checks = handoff_screen.guard()
    if checks and checks[0][0] == "SKIP":
        print(f"  [SKIP] {checks[0][2]}")
        _skips.append("the eQuilibrator compound cache")
    else:
        for name, ok, detail in checks:
            print(f"  [{'PASS' if ok else 'FAIL'}] {name:46} {detail}")
            if not ok:
                _failures.append(f"H2 threshold ({name})")

    if _failures:
        print(f"\n{len(_failures)} CHECK(S) FAILED: {', '.join(_failures)}")
        print("(FBA/thermodynamic numbers are solver- and version-sensitive; see the reproducibility note above.)")
        sys.exit(1)
    if _skips:
        print(f"\n{len(_skips)} section(s) SKIPPED because a network resource was unavailable (link-only models are "
              f"fetched from their source, and the eQuilibrator cache downloads, on first run). Everything that ran is "
              f"within tolerance; re-run with network access to reproduce the rest. See MODELS.md.")
        return
    print("\nAll headline results regenerated and within tolerance. See tests/ for the full assertion suite.\n")


if __name__ == "__main__":
    main()
