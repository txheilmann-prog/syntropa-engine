# Demo model licensing

Every demo manifest (our metadata + port map) ships in `library/`. A model's payload file ships directly when its upstream license permits us to re-serve the bytes -- curated CC0 BioModels and our own first-principles models -- and is **link-only** otherwise: the manifest ships, the model file does not, and `reproduce.py` fetches it from the source named in the manifest's `provenance` on first run (to place one by hand, obtain it from that source and save it at the path the manifest gives).

## Bundled (model file ships here)

| manifest | model file | source |
|---|---|---|
| `fba_methanotroph.json` | `library/data/imethanotroph_lumped.xml` | biomass composition, derivation, growth yield |
| `glycolysis_teusink.json` | `library/data/teusink_glycolysis.xml` | biomodels-cc0 |
| `pritchard_glycolysis.json` | `library/data/pritchard_glycolysis.xml` | biomodels-cc0 |
| `signaling_kholodenko.json` | `library/data/kholodenko_mapk.xml` | biomodels-cc0, kholodenko2000 |

## Link-only (manifest ships; model file fetched at runtime)

| manifest | model file | upstream source |
|---|---|---|
| `fba_cyano_auto.json` | `library/data/icyano_auto.xml` | BIGG (cite per the manifest) |
| `fba_cyanobacterium.json` | `library/data/icyano.xml` | BIGG (cite per the manifest) |
| `fba_ecoli_core.json` | `library/data/e_coli_core.xml` | BIGG (cite per the manifest) |
| `fba_methanogen.json` | `library/data/imethanogen.xml` | BIGG (cite per the manifest) |
