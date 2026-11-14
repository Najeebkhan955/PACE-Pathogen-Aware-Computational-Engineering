# Reproducibility guide

## 1. Verify scientific data checksums

From the repository root:

```bash
python code/verify_checksums.py
python code/verify_repository_manifest.py
```

These are **integrity checks**, not a claim that all model training environments run with zero configuration.

## 2. Source files and run order

- [Original PACE scripts](../code/original_PACE/): author-supplied Python code. The 14 original `.py` files have been separately extracted for code review.
- [Original run order](../code/original_PACE/RUN_ORDER.md): exact script sequence as supplied.
- [S12](../S12/): code and selected independent verification scripts; `.ipynb` notebooks may run in Google Colab.
- [S01](../S01/): curation source and measured MIC/HC50 provenance.

## 3. Numerical audit route

| Task | Evidence |
|---|---|
| Reproduce AMP classification comparisons | `S02/DATA/` |
| Inspect six MIC out-of-fold predictions | `S04/DATA/` |
| Inspect HC50 and uncertainty evaluation | `S05/DATA/` |
| Inspect ESM-2 and QMAP evaluations | `S06/DATA/` |
| Compare MTL and species-aware endpoints | `S07/DATA/` |
| Validate lead sequences and ranking criteria | `S08/DATA/` |
| Inspect GA ablations and controlled NSGA-II archives | `S09/DATA/`, `S10/DATA/` |
| Inspect external generative library comparisons | `S11/DATA/` |

## 4. Availability and computational caveats

Original model files and the original GA search trajectory were **not included** with the author's source inputs. Retained analysis outputs include new out-of-fold calculations, reconstructed external-model runs and score-based stress tests, with their provenance identified. These may not reproduce an unavailable historical run bit-for-bit. Consult model-specific notebooks and `DATA/` manifests for software versions before execution.

## 5. Reusing the data

- Do not conflate modeled MIC/HC50 with experimentally measured values.
- Treat sequence grouping, model reconstruction, and dataset overlap checks as mandatory analysis context.
- Consult [Attribution and licenses](04_ATTRIBUTION_AND_LICENSES.md) before redistribution.
