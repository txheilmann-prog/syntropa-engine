# Syntropa

A physics-gated substrate for composing **heterogeneous** biological models -- first-principles reactions, empirical
lookup tables, whole kinetic SBML models, and genome-scale flux-balance (FBA) models -- into micro-ecosystems that
exchange matter through a shared, chemically-grounded medium. Distributed as the `syntropa-engine` package; imported as `microcosm`.

Its distinguishing feature is a **layered validity stack**: every composed result is checked by independent gates
(composition-validity, rigorous community FBA, whole-system conservation, thermodynamic feasibility, per-box
validity domains), so a composition that is physically impossible, non-conserving, or thermodynamically infeasible
is caught before its result is trusted. Shared identities are grounded (ChEBI/BiGG), quantities are dimensioned,
and the pipeline is anchored to a published known-answer benchmark.

This is the **open substrate**. It composes *on top of* Vivarium, libRoadRunner, COBRApy, MICOM, and eQuilibrator;
the contribution is the verifiability discipline, not a new composition framework.

## Install

    pip install -e .

Python 3.11+; imported as `microcosm`. Core dependencies (see `pyproject.toml`): cobra, micom, equilibrator-api,
libroadrunner, vivarium-core, numpy, scipy, pandas.

This is a source-checkout project: **clone the repository and run from its root.** `pip install -e .` installs the
`microcosm` package for `import`, but the demo models (`library/`), `reproduce.py`, and the tests live in the
checkout, not inside the installed package -- so run them from the repository root, not from a bare index install.

## Reproduce the headline results

Run from the source checkout (`reproduce.py` and the `library/` paths resolve relative to it):

    python reproduce.py

Regenerates, from the bundled and link-only demo models: the known-answer *E. coli* core growth rates (aerobic 0.874/h, anaerobic
0.211/h), the thermodynamic gate rejecting a free-energy-infeasible reverse reaction, and an emergent *E. coli* +
methanogen cross-feeding that produces methane neither member makes alone, and the per-hand-off thermodynamic screen
reproducing the measured interspecies-H2 thresholds of the sulfate-reducer, methanogen, and acetogen guilds
(Cord-Ruwisch et al. 1988). Each numeric result is checked against its expected value within a stated tolerance. A few demo models are link-only (their upstream license does not
permit redistribution); `reproduce.py` fetches them from their source on first run, and, run without network, skips
just the sections whose model or eQuilibrator cache it cannot fetch (with a clear note) rather than failing. See `MODELS.md`.

FBA growth rates depend on the LP solver and its version (the reported values reproduce to within ~1% across
solver stacks). For bit-stable numbers, install the pinned environment in `requirements-lock.txt`. The first run
that uses the thermodynamic gate (including `reproduce.py`) downloads the eQuilibrator compound cache (a few
hundred MB) once; later runs reuse it.

## Public API

Run from the repository root (the `library/` paths are relative):

    import microcosm
    ecoli = microcosm.load_component("library/fba_ecoli_core.json")
    # thermo=False keeps this quick; the thermodynamic gate (below) downloads the eQuilibrator cache on first use
    result = microcosm.community_transform([ecoli], feed={"glucose": 10, "O2": 1000}, thermo=False)
    print(result["community_growth"])                                                     # ~ 0.874 /h

    # the other entry points (see reproduce.py and tests/ for complete, runnable calls):
    #   microcosm.composition_validity(components)                                -> physical-coexistence gate
    #   microcosm.super_additivity(components, target, available=..., feed=...)   -> emergence discriminator
    #   microcosm.route_dg(substrates, products, id_resolver)                     -> thermodynamic gate (downloads cache)

To compose your own models, author a manifest following `SCHEMA.md` (copy the closest example in `library/`).

## Security

`load_component` and the community-FBA path parse model files (SBML/JSON) named by a manifest, so load only
manifests and model files you trust. The runtime fetch of a link-only demo model is validated (a BiGG model id
of the expected form), restricted to its cited host (no cross-host redirects), size-capped, and written only
inside the `library/` tree.

## How to cite

Please cite the Syntropa preprint: https://doi.org/10.5281/zenodo.22129408 (a software paper is in preparation).
Each bundled demonstration model retains its own original source and citation; see `MODELS.md`.

## License

Apache-2.0 (`LICENSE` and `NOTICE`). Models under `library/` carry their own upstream licenses (`MODELS.md`).
Runtime dependencies carry their own licenses (including LGPL for libSBML and cobra); using them as libraries does
not affect the license of your own code.
