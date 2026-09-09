# Component manifest format

A component is one JSON manifest in `library/` describing a box in the netlist: its shared-species ports, how it
transforms them, and the physical conditions under which it can coexist with other boxes. `load_component()` reads
and validates it (returning ALL problems at once if malformed). The shipped manifests are the canonical examples --
`fba_methanogen.json` is a full genome-scale FBA box, `respiration.json` a first-principles function box; copy the
closest one and adapt it.

## Fields

- `id` (str): unique id. `name` (str): human-readable name. `kind` (str): what the box is (e.g. `organism`,
  `chemistry`). `schema_version` (str): must equal the loader's `SCHEMA_VERSION` (`"1.0"`).
- `open_boundary` (bool): true if the box exchanges matter with an UNMODELLED compartment (feed, atmosphere,
  biomass); its net port atom-flux is then a declared boundary the conservation monitor subtracts.
- `provenance` (list of `{source, ref, retrieved, confidence?}`): where the model and its facts come from. The
  `source` string drives the licensing policy (see `MODELS.md` / `licensing.py`) -- e.g. `BioModels-CC0`, `BiGG`,
  or a first-principles tag (`first-principles`, `biomass`).
- `transfer` (obj): the box's behavior; `type` selects the engine:
  - `"fba"` -- flux-balance box: `path` (SBML file relative to `library/`, e.g. `data/imethanogen.xml`), `vmax`,
    and `port_map` mapping each SBML exchange-reaction id to a port species name.
  - `"function"` -- a registered transfer function: `name` is the registry key (see `builtins.py`).
- `params` (obj, optional): `{name: {default, range: [lo, hi]}}` tunable parameters.
- `ports` (obj): the shared-species interface -- `{species: {formula, role, chemical, entity_id, unit?}}`.
  `formula` is the actual-atom chemical formula (for mass/charge checks); `role` is `input` or `output`;
  `chemical` marks a real species vs an abstract signal; `entity_id` is the canonical id (ChEBI/BiGG) by which the
  same species is reconciled across boxes.
- `validity` (obj, optional): `{note}` free-text scope/caveats.
- `context` (obj, optional): the DECLARED environmental preconditions the composition-validity gate checks for
  mutual satisfiability BEFORE any simulation -- `phase` (e.g. `aqueous`), `requires` (species that must be
  present), `ranges` (`{intensive_var: [min, max]}`, e.g. `Eh_mV`, `pH`, `temperature`), and an optional
  `o2_relationship` descriptor. The gate intersects each variable's range across all boxes; a disjoint
  intersection means the community cannot physically coexist and is rejected (see `reproduce.py` section 4).

Species names in `feed` / `available` match the port species names (see `_FEED_ALIAS` in `community.py` for
accepted aliases).
