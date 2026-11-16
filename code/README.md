# PACE computational research code

The Python source files document the PACE sequence-curation, descriptor-generation, classification, quantitative modeling, evolutionary optimization, and lead-prioritization workflows. Complementary scientific analyses and numerical verification resources are organized in the supplementary collections.

| Resource | Purpose |
|---|---|
| [Source pipeline](original_PACE/) | Curation, 90-descriptor features, AMP recognition, MIC/HC₅₀ modeling, genetic search and lead selection |
| [Execution guide](original_PACE/RUN_ORDER.md) | Declared source-script sequence and software requirements |
| [Scientific analysis scripts](../S12/DATA/) | Model comparisons, uncertainty, objective evaluation and additional numerical records |
| [verify_checksums.py](verify_checksums.py) | Verify scientific files against the research SHA-256 manifest |
| [verify_repository_manifest.py](verify_repository_manifest.py) | Verify published documentation and supporting files |
| [reassemble_archives.py](reassemble_archives.py) | Restore and verify complete S06/S10 numerical archives |

## Verification commands

From the repository root:

```bash
python code/verify_checksums.py
python code/verify_repository_manifest.py
python code/reassemble_archives.py --verify-only
```

Reconstruct the large numerical archives as needed:

```bash
python code/reassemble_archives.py
```

The repository retains computational configurations and available result records. Reproducing particular model training runs depends on compatible source data, fitted artifacts, runtime environments and declared random seeds. Reported MIC/HC₅₀ values for the 20 designed PACE peptides remain model-derived predictions.
