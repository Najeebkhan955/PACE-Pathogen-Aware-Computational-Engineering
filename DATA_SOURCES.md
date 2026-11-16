# PACE data sources and provenance

PACE integrates publicly available peptide sequences, measured antimicrobial and hemolysis activity, and computational benchmark outputs. See [the source-provenance table](DATA_SOURCE_PROVENANCE.csv) for upstream datasets, study-specific inputs, and associated citations.

## Research sources

| Dataset or method | Scientific role |
|---|---|
| AMPlify | AMP recognition model training and held-out evaluation |
| amPEPpy M-model corpus | Independent external classification evaluation after sequence-overlap filtering |
| GRAMPA | Experimentally measured MIC records for six ESKAPE-associated quantitative endpoints and activity controls |
| Hemolytik-derived dataset | Source erythrocyte hemolysis observations for HC₅₀ regression |
| AMPScanner v2, ESM-2, QMAP | Classification or quantitative prediction benchmarks using defined evaluation cohorts |
| HydrAMP, AMPGAN v2 | Exploratory peptide library comparisons under a common scoring surrogate |

The *E. faecium*-associated MIC endpoint pools *E. faecium* and *E. faecalis*; the *E. cloacae*-associated endpoint includes records historically named *E. aerogenes* (now *Klebsiella aerogenes*). A pooled endpoint should not be described as a strain-specific susceptibility estimate.

## Interpretation and data coverage

- **Experimentally measured source records** and **predictions for 20 designed PACE leads** are distinct evidence types. The designed peptide MIC, HC₅₀, and selectivity values are model-derived.
- Strongly active (≤10 µM) and target-specifically weak (≥64 µM) measured MIC labels do not define universal AMP or non-AMP behavior.
- Exact-sequence deduplication and overlap filtering do not necessarily eliminate homologous sequences, assay heterogeneity, strain differences, or changes in peptide modifications.
- The ESM-2, NSGA-II, and generative-library comparisons are computational evaluations; shared-surrogate scores are not experimentally verified antimicrobial hit rates.
- S01–S12 retain the numerical records, file provenance, and data-integrity information needed for inspection of reported analyses.

See the [supplementary index](SUPPLEMENTARY_INDEX.csv), [reproducibility guide](docs/03_REPRODUCIBILITY.md), and [source attribution guidance](docs/04_ATTRIBUTION_AND_LICENSES.md).
