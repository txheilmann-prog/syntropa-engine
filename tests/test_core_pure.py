"""Pure, deterministic tests of the substrate core -- formula/mass handling, whole-system conservation, and the
source-licensing policy. These run in seconds with NO solver stack (cobra/micom) and NO network, so a reviewer
can check the core directly; the FBA/thermodynamics path is covered by tests/test_public_api.py.
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

from microcosm import licensing, units                       # noqa: E402
from microcosm.conservation import ConservationMonitor        # noqa: E402


def test_molar_mass_and_formula_parsing():
    assert abs(units.molar_mass("H2O") - 18.015) < 1e-2
    assert 180.0 < units.molar_mass("C6H12O6") < 180.3        # glucose ~ 180.16 g/mol
    assert units.parse_formula("CO2") == {"C": 1, "O": 2}
    assert units.parse_formula("C6H12O6") == {"C": 6, "H": 12, "O": 6}


def test_conservation_passes_an_atom_balanced_reaction():
    mon = ConservationMonitor({"glucose": "C6H12O6", "o2": "O2", "co2": "CO2", "h2o": "H2O"})
    series = [{"glucose": 1.0, "o2": 6.0, "co2": 0.0, "h2o": 0.0},
              {"glucose": 0.0, "o2": 0.0, "co2": 6.0, "h2o": 6.0}]   # C6H12O6 + 6 O2 -> 6 CO2 + 6 H2O
    assert mon.check(series).conserved


def test_conservation_flags_an_atom_leak():
    mon = ConservationMonitor({"glucose": "C6H12O6", "o2": "O2", "co2": "CO2", "h2o": "H2O"})
    series = [{"glucose": 1.0, "o2": 6.0, "co2": 0.0, "h2o": 0.0},
              {"glucose": 0.0, "o2": 0.0, "co2": 7.0, "h2o": 6.0}]   # an extra CO2 -> carbon from nowhere
    report = mon.check(series)
    assert not report.conserved
    assert "C" in report.leaking


def test_heterogeneous_composition_conserves_mass():
    """Compose a first-principles function box with an empirical-lookup box over one shared medium, run the
    Vivarium trajectory, and confirm whole-system mass conservation -- the composition, validity, and
    conservation core exercised end to end, offline (no FBA solve, no network)."""
    import microcosm.builtins  # noqa: F401  (registers the function transfers, e.g. aerobic_respiration)
    from microcosm import load_component
    from microcosm.engine import run_netlist
    lib = os.path.join(ROOT, "library")
    resp = load_component(os.path.join(lib, "respiration.json"))       # first-principles function
    ferm = load_component(os.path.join(lib, "fermenter_lookup.json"))  # empirical lookup table
    res, report = run_netlist([resp, ferm],
                              {"initial_medium": {"glucose": 20.0, "O2": 10.0}, "duration": 20.0})
    assert res["composition_problems"] == []      # the two boxes can share one well-mixed medium
    assert res["diverged"] is None                # the trajectory did not blow up
    assert report.conserved                       # atoms balance over the whole trajectory
    assert not report.leaking


def test_licensing_is_fail_closed_and_correct():
    assert licensing.may_redistribute("BioModels-CC0 (BMDB)") is True           # CC0, freely re-servable
    assert licensing.may_redistribute("BiGG") is False                          # link-only; cite the source
    assert licensing.may_redistribute(None) is False                            # unknown -> deny (fail closed)
    assert licensing.may_redistribute("some brand new unvetted source") is False
    assert licensing.in_subset("BiGG", "commercial") is False
    assert licensing.policy_for("AGORA2 v2.01 (Heinken 2023)")["redistribute"] is False
