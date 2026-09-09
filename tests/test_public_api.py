"""Public-API tests + the known-answer benchmark. Runnable by a stranger:

    pip install -e ".[test]"
    pytest tests/test_public_api.py

Exercises the substrate's public surface (load_component, community_transform, super_additivity, route_dg) and
reproduces the published E. coli core growth rates. Link-only demo models are fetched from their source on first
run (see MODELS.md), so this runs from a clean checkout.
"""
import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

pytest.importorskip("cobra")
pytest.importorskip("micom")

import microcosm.builtins  # noqa: E402, F401  (registers transfer functions)
from microcosm import community_transform, load_component, route_dg, super_additivity  # noqa: E402


def _lib(name):
    return load_component(os.path.join(ROOT, "library", name))


def test_ecoli_core_known_answer_growth_rates():
    """Community FBA reproduces the canonical E. coli core rates (Orth, Fleming & Palsson 2010)."""
    ecoli = _lib("fba_ecoli_core.json")
    aerobic = community_transform([ecoli], feed={"glucose": 10, "O2": 1000}, thermo=False)["community_growth"]
    anaerobic = community_transform([ecoli], feed={"glucose": 10, "O2": 0}, thermo=False)["community_growth"]
    assert aerobic == pytest.approx(0.8739, abs=2e-3)
    assert anaerobic == pytest.approx(0.211, abs=3e-3)     # canonical ~0.211/h; FBA is solver-version-sensitive


def test_thermodynamic_gate_rejects_the_infeasible_reverse():
    """The thermodynamic gate accepts glucose -> 2 ethanol + 2 CO2 and rejects its confidently-endergonic reverse."""
    forward = route_dg({"glc__D": 1}, {"etoh": 2, "co2": 2}, lambda s: s)
    reverse = route_dg({"etoh": 2, "co2": 2}, {"glc__D": 1}, lambda s: s)
    assert forward["thermo_ok"] is True
    assert reverse["thermo_ok"] is False


def test_ecoli_methanogen_methane_is_emergent():
    """E. coli + methanogen produce methane that neither member makes alone: a super-additive (emergent) result."""
    ecoli = _lib("fba_ecoli_core.json")
    methanogen = _lib("fba_methanogen.json")
    emg = super_additivity([ecoli, methanogen], "methane",
                           available={"glucose", "ammonium", "phosphate", "CO2", "H2O"}, feed={"glucose": 10})
    assert emg["emergent"] is True
    assert emg["community"] > emg["best_subcommunity"]
