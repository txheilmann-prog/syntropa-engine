---
title: 'Syntropa: a physics-gated substrate for composing and validating heterogeneous cell and ecosystem models'
tags:
  - Python
  - systems biology
  - model composition
  - co-simulation
  - flux balance analysis
  - reproducibility
authors:
  - name: Craig Heilmann
    orcid: 0009-0007-6648-7076
    affiliation: 1
affiliations:
  - name: Independent researcher, USA
    index: 1
date: 9 September 2026
bibliography: paper.bib
---

# Summary

Syntropa is a Python substrate for composing *heterogeneous* biological models: first-principles reactions,
empirical lookup tables, kinetic SBML models, and genome-scale flux-balance (FBA) models, assembled into
micro-ecosystems that exchange matter through a shared, chemically grounded medium. It builds on Vivarium
[@agmon2022vivarium] for dynamic composition, libRoadRunner [@somogyi2015libroadrunner] for SBML kinetics,
COBRApy [@ebrahim2013cobrapy] for FBA, MICOM [@diener2020micom] for community FBA, and eQuilibrator
[@beber2022equilibrator] for thermodynamics. Its organizing idea is a **layered validity stack**: independent
gates, cheapest first, that reject a composition which cannot share one well-mixed medium, does not conserve
mass or charge, or whose net transform is thermodynamically infeasible, before its result is trusted. The stack
rests on a grounded-identity guard (ChEBI/BiGG identity, actual-atom formulas, charge) and a dimensioned-units
discipline. This paper describes the open composition-and-validation substrate; larger-scale search built on it
is reported elsewhere.

# Statement of need

Composing heterogeneous models is an active frontier: Vivarium [@agmon2022vivarium], whole-cell models
[@karr2012wholecell], and community-metabolism tools such as COMETS [@dukovski2021comets] all compose submodels
over shared state, but they generally validate only *computational* wiring (do ports connect, do types match)
and leave *physical* validity to the modeler. A heterogeneous composition can wire up cleanly, run to
completion, and still be physically inconsistent (incompatible phases or disjoint pH windows), silently
non-conserving (a boundary that leaks atoms or charge), or thermodynamically infeasible (a net transform that
runs uphill in free energy). Model-credibility tools such as MEMOTE [@lieven2020memote] check *single* models,
not the compatibility *between* composed ones. Syntropa targets that gap: a small, verifiable substrate in which
a heterogeneous composition is screened by physics before it is solved, identities and units are enforced, and
every model value is sourced.

# Software design and functionality

A model is a set of declarative **component manifests** (JSON, versioned in git). Each names the component's
*ports* (its interface to the medium, each a species with a grounded identity, role, and units), a *transfer*
(a first-principles function, an empirical lookup, a kinetic SBML model, or an FBA model), and a *validity
domain*. Components free-associate through the medium: two naming the same species must agree on its formula.
The public surface is small and versioned (`load_component`, `composition_validity`, `community_transform`,
`super_additivity`, `route_dg`, and the hand-off screen); each heavy engine sits behind one seam and is *used*,
not reinvented. Runs record the engine and schema versions and a declarative netlist, and the constraint-based
solve is deterministic; `requirements-lock.txt` pins solver versions.

The **validity stack**, cheapest first:

1. *Composition validity*: a static gate for whether the submodels can share one well-mixed medium. From each
   component's declared preconditions (phase; required constituents; ranges of pH, temperature, redox potential)
   it checks phase consistency, intersects the declared ranges, and closes implications (a charged species
   implies an aqueous phase). It screens the shared-medium assumption from modeler-declared metadata; it does
   **not** represent spatially or temporally structured communities (biofilms, sediments, diel redox cycling),
   where members that cannot share one well-mixed state still coexist. This applies contract reasoning
   [@benveniste2018contracts], distinct from conservation-by-construction composition [@shahidi2021bondgraph].
2. *Community FBA* [@diener2020micom]: members are solved in one steady-state optimization with the
   shared-medium mass balance as a hard constraint, so per-organism over-draw cannot occur within the joint
   solve (as in any joint formulation, e.g. SteadyCom [@chan2017steadycom], OptCom [@zomorrodi2012optcom]).
   Cross-feeding is resolved under the chosen community objective; the representative solution is regularized by
   parsimonious FBA [@lewis2010pfba] or MICOM's cooperative-tradeoff, with residual optima characterized by flux
   variability.
3. *Conservation*: whole-system atom and charge conservation over the trajectory, subtracting declared open
   boundaries; for a steady-state solve this mainly re-checks boundary arithmetic.
4. *Thermodynamic feasibility* [@noor2013componentcontribution; @beber2022equilibrator]: a net transform whose
   transformed Gibbs energy, at eQuilibrator's 1 mM physiological reference, is confidently endergonic is
   rejected. It is a necessary-condition filter on the *net* transform, so individual interspecies hand-offs
   cancel and are not themselves gated (per-hand-off feasibility is a separate screen). Evaluating at a fixed
   reference, it can miss a transform feasible only at other activities; and an unmapped species
   (terminal-methanogenesis C1 carriers are not yet covered) leaves its reaction unconstrained. The filter never
   fabricates a free energy: it can miss an infeasibility, not invent one. For per-reaction energies,
   low-potential carriers component contribution cannot estimate (ferredoxin, coenzyme F420, methanophenazine)
   come from literature redox potentials with propagated uncertainty (these cancel in a net transform).
5. *Per-box validity domains*: warn when a box is driven outside its characterized range.

An optional enzyme layer (sMOMENT [@bekiaris2020gecko]) adds a nominal protein-pool bound that can reproduce the
qualitative onset of overflow metabolism. A super-additivity ("emergence") test flags a candidate consortium
only when its target output exceeds the best sub-community: a *necessary* condition for synergy, not a sufficient
one (obligate-intermediate and mechanistic-distance tests are complementary, and beyond pairs the FBA walk tends
to collapse to single organisms). It is itself thermodynamics-blind; the bundled *E. coli* plus *M. barkeri*
methane example reproduces a known syntrophy [@hamilton2015], not a discovery.

**Known-answer validation.** `reproduce.py` regenerates each anchor in one command. The community-FBA path
reproduces the canonical *E. coli* core growth rates [@orth2010ecolicore] (aerobic 0.874/h, anaerobic 0.211/h)
within a per-number tolerance, validating the FBA-engine integration; a heterogeneous first-principles-plus-lookup
netlist runs offline and conserves atoms over its trajectory, exercising composition, the medium, and
conservation. The thermodynamic layer is *calibrated against, and consistent with*, the measured
interspecies-hydrogen thresholds of the three scavenging guilds (sulfate reducers below about 1 Pa, methanogens
1 to 10 Pa, acetogens 50 to 95 Pa) [@cordruwisch1988]: their ordering follows from the component-contribution
free energies alone, while the absolute thresholds also use the literature Schink quantum (about 20 kJ/mol
[@schink1997], applied to the energy-conserving scavenger) at nominal activities, so this is a consistency
check, not a parameter-free prediction. The FBA and thermodynamic demo models are fetched from BiGG and the
eQuilibrator cache on first run, so those sections require network access and skip with a note otherwise, and a
run that verifies nothing exits non-zero. An automated suite runs in continuous integration; MEMOTE
[@lieven2020memote] remains the standard single-model check.

# State of the field

Syntropa composes *on top of* Vivarium, COBRApy, MICOM, and eQuilibrator rather than replacing them, and is not a
substitute for production community tools (COMETS, MICOM, SteadyCom [@chan2017steadycom], OptCom
[@zomorrodi2012optcom]) or for thermodynamics-based flux analysis (TMFA [@henry2007tmfa]), of which its
net-transform check is a lightweight special case. Predicting coexistence and cross-feeding from genome-scale
models is the domain of SMETANA [@zelezniak2015smetana], MICOM, and COMETS, shown across many communities
[@machado2021polarization]; the methane coculture chemistry the example reproduces was analyzed by Hamilton and
colleagues [@hamilton2015]. Its contribution is a verifiability discipline for *heterogeneous* composition. The
least-precedented piece is the composition-validity gate, which, to the author's knowledge, is not integrated
into heterogeneous systems-biology composition tooling, though static precondition screening is established in
contract-based design [@benveniste2018contracts] and superstructure synthesis [@chen2021pyosyn]; identity
reconciliation follows the semantic-annotation tradition [@krause2010semanticsbml]. Anchoring composed results to
published known answers, rather than to internal consistency alone, is what makes them auditable by a stranger.

# AI usage disclosure

Large language models (Anthropic Claude, Opus and Sonnet models) were used in this work: to assist with drafting
and copy-editing the paper, and with code review, refactoring, and test scaffolding for the software. All
scientific choices, the software design, model-value sourcing, and validation were made and directed by the
author, who reviewed, edited, and verified every AI-assisted output; no result, number, or citation was accepted
without checking it against its source.

# Acknowledgements

Developed as an independent research project.

# References
