# PACE reproducibility and computational records

## Data and archive integrity

From the repository root with Python 3:

```bash
python code/verify_checksums.py
python code/verify_repository_manifest.py
python code/reassemble_archives.py --verify-only
```

Checksum verification confirms data-file consistency, but does not by itself rerun model training or demonstrate experimental activity.

## Source scripts and analyses

- [PACE Python source scripts](../code/original_PACE/) and [declared run order](../code/original_PACE/RUN_ORDER.md) document the sequence curation, 90-descriptor encoding, AMP classifier, MIC/HC₅₀ regressors, GA search, and candidate selection stages.
- [S12](../S12/) contains additional computational scripts, numerical checks, and run records.
- [S01](../S01/) describes the input data and source provenance.

| Analysis | Supporting records |
|---|---|
| AMP recognition and external predictor comparison | [S02/DATA](../S02/DATA/) |
| Measured target-specific activity controls | [S03/DATA](../S03/DATA/) |
| Six MIC endpoints and predictive performance | [S04/DATA](../S04/DATA/) |
| Hemolysis prediction and uncertainty | [S05/DATA](../S05/DATA/) |
| ESM-2, hybrid, composition–length, and QMAP evaluation | [S06/DATA](../S06/DATA/) |
| Multitask learning, species pooling, and cross-pathogen MIC | [S07/DATA](../S07/DATA/) |
| PACE leads, diversity, and applicability domain | [S08/DATA](../S08/DATA/) |
| GA objective ablations | [S09/DATA](../S09/DATA/) |
| Controlled GA versus NSGA-II comparison | [S10/DATA](../S10/DATA/) |
| HydrAMP and reconstructed AMPGAN v2 libraries | [S11/DATA](../S11/DATA/) |
| Executable code and computational records | [S12/DATA](../S12/DATA/) |

## Scope of reproducibility

The deposited research code and retained output records enable inspection and reuse of the documented workflows. Some comparator models use specified refitted or reconstructed implementations. **Exact replay of the initial GA search may require fitted-model artifacts that are not included.** Complete provenance, model-specific settings, and the data files should be consulted when reproducing results.

The 20 designed peptide scores are model predictions, not experimentally determined MIC or HC₅₀ measurements. Comparing candidate libraries under a shared surrogate does not establish biological hit rates.

## Reassembling the large numerical archives

Complete ZIP contents are stored in ordered, byte-preserving parts in `S06/DATA/` and `S10/DATA/`. To reassemble and verify them:

```bash
python code/reassemble_archives.py
```

See [provenance](../DATA_SOURCES.md) and [licensing notes](04_ATTRIBUTION_AND_LICENSES.md).
