# Changelog

All notable changes to this project are documented here. The format follows Keep a Changelog, and the project
uses semantic versioning.

## [0.1.2] - 2026-10-20

Three of the four changes here are caveats and measured bounds rather than features. A validity stack is
only worth what its stated limits are worth, so the limits ship.

### Changed
- `handoff_screen`: the agreement with the published Cord-Ruwisch acetogen H2 threshold is now stated as
  CONDITIONAL ON THE SHIPPED SOLUTE CONCENTRATIONS, beside the claim rather than in a footnote. Swept
  across the cited acetate range this screen stops reproducing that threshold at BOTH ends, at 46.09 Pa
  low and 115.78 Pa high, against a published band of 50 to 95 Pa. Being inside a cited concentration
  window at a single point says nothing about the window. The agreement is real and it is not
  unconditional: a caller changing a concentration changes what the screen validates against.
- `handoff_screen`, `thermo`: documented, with measurements, that eQuilibrator does NOT
  temperature-correct dG'0. Setting `cc.temperature` is accepted and reported back while the value does
  not move: hydrogenotrophic methanogenesis reads -195.10 kJ/mol at 298.15, 313.15 and 328.15 K alike,
  and propionate oxidation shifts 0.01 kJ/mol over 30 K. Raising `RT` to match a warmer reactor is
  therefore a HALF correction, with the concentration term at the new temperature and dG'0 pinned at
  25 C, and that partial correction measurably worsened agreement with documented obligate and
  facultative lifestyles in our own testing. Treat every verdict as sign feasibility at 25 C.
- `handoff_screen`: every verdict is BULK PHASE and conditional on the members achieving local proximity.
  There is one shared concentration per metabolite and no spatial, diffusion or aggregation variable,
  while interspecies H2 and formate transfer is diffusion limited and the decisive concentration is at
  the cell to cell interface. The feasible windows on the obligate syntrophies span a factor of only
  about 1.5 to 4, which describes a controlled microenvironment and not a stirred vessel; the internal
  control is the not-obligate case, where lactate oxidation returns a window six orders of magnitude
  wide. Bulk phase feasibility is NECESSARY AND NOT SUFFICIENT for a hand-off to occur.

### Added
- `thermo.set_conditions()`, the spine entry point for eQuilibrator aqueous conditions, returning the
  conditions actually in effect so a caller can assert them. `ionic_strength_M` and `p_h` enter through
  the Legendre and Debye-Huckel terms, need no entropy data, and ARE honored; `temperature_K` is not,
  for the reason above, and the docstring says so where a caller will read it.
- `enzyme.PILOT_BUDGET`, naming the 0.1 g protein per gDW that the pilots, the use case search and every
  evidence harness actually pass. `enzyme_constrain(budget=0.04)` was the signature default and nothing
  in the tree passed it, so the live value lived as a bare literal at dozens of call sites. It is labeled
  a CHOSEN SETTING rather than a measured parameter, because neither number is derived from a measurement
  or a citation. What IS measured is its effect, swept over 0.04, 0.2 and 0.4: emergence verdicts are
  invariant across the whole range, so no binary discovery claim rests on it; absolute production tracks
  the budget below 0.1 and saturates above it, unevenly across candidates, so a flux computed at one
  budget cannot be converted to another by arithmetic; and above roughly 0.1 the pool goes slack and the
  constraint stops binding, with six of eight abundance ratios giving the identical unblocked optimum at
  0.2 and at 0.4. Any claim whose magnitude depends on the budget must be quoted with the budget.
- `emergence`: optional `abundances`, mapping component id to biomass fraction, passed through to the
  solver. The asymmetry is documented where it matters: micom renormalizes whatever fractions a subset
  is left with, so a singleton always solves at abundance 1.0 regardless of the community ratio. The
  singleton baselines are therefore abundance INDEPENDENT while the community term is not, which means
  synergy, as community minus best subset, absorbs the entire ratio effect.

### Fixed
- `handoff_screen`: the concentration scans evaluated dG'0 once per point when it is concentration
  independent, making n+1 identical eQuilibrator calls per scan. Hoisted out of the loop, and `_lnQ` is
  split out so the scan path and the single point path share one implementation and cannot drift.

## [0.1.1] - 2026-09-10

### Added
- The JOSS paper (paper.md, paper.bib).
- An offline heterogeneous-composition and whole-system mass-conservation check in reproduce.py, and a matching test.

### Changed
- reproduce.py exits non-zero when no numeric check ran (an all-skip run is no longer a silent pass).
- Clarified demo-manifest notes; corrected the E. coli core citation.

## [0.1.0] - 2026-09-09

First public release: the open simulator substrate for composing and physics-validating heterogeneous
biological models.

### Added
- Component-manifest composition over a shared, chemically-grounded medium.
- The layered validity stack: composition-validity, rigorous community FBA (MICOM/pFBA), whole-system
  conservation, thermodynamic feasibility (eQuilibrator + curated low-potential electron carriers), and per-box
  validity domains, plus an optional enzyme-constraint (sMOMENT) layer.
- Grounded identity (ChEBI/BiGG, actual-atom formulas, charge) and dimensioned units.
- A super-additivity emergence test and a per-hand-off thermodynamic feasibility screen.
- A published known-answer benchmark (E. coli core, Orth et al. 2010), reproduce.py, and a public-API test.
