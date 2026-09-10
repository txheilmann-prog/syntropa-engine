# Changelog

All notable changes to this project are documented here. The format follows Keep a Changelog, and the project
uses semantic versioning.

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
