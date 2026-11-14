# PACE — Pathogen-Aware Computational Engineering

**Computational antimicrobial peptide discovery and de novo design for ESKAPE-priority pathogens**

[![Data](https://img.shields.io/badge/Supplementary-S01%E2%80%93S12-235b80)](./supplementary)
[![Status](https://img.shields.io/badge/Status-data%20deposition%20in%20progress-e8a431)](./supplementary)
[![License](https://img.shields.io/badge/License-see%20data%20sources-lightgrey)](./DATA_SOURCES.md)

## Overview

PACE combines sequence-based machine-learning models, physicochemical descriptors, six ESKAPE-associated MIC prediction endpoints, HC50/hemolysis prediction, and genetic-algorithm candidate prioritization. Subsequent computational analyses evaluate pretrained ESM-2 representations, multi-task modeling, measured target-specific activity controls, species pooling, an NSGA-II comparator, and external AMP design approaches.

**Important distinction:** The 20 designed PACE peptides are *computational candidates*. Their predicted activity and selectivity are not experimental measurements.

## Repository status

| Collection | Current status | Link |
|---|---|---|
| S01–S12 Excel supplementary summaries | **Uploaded** | [Browse 12 workbooks](./supplementary) |
| Supplementary CSV index | **Uploaded** | [Index](./supplementary/SUPPLEMENTARY_INDEX.csv) |
| Full row-level CSV matrices | **Deposit in progress** | [Data landing page](./data/README.md) |
| Original and analysis Python scripts | **Deposit in progress** | [Code landing page](./code/README.md) |
| Peptide FASTA sequences | **Deposit in progress** | [Sequences landing page](./sequences/README.md) |

The summary workbooks **do not replace** full underlying row-level data. No complete-data or independent reproduction claim should be made until the remaining files are uploaded and verified.

## Supplementary data guide

| Group | Topic |
|---|---|
| **S01** | Dataset provenance and activity/MIC data curation |
| **S02** | AMP classification and independent AMPScanner v2 comparison |
| **S03** | Measured peptide–pathogen activity control benchmarks |
| **S04** | Six ESKAPE MIC prediction endpoints |
| **S05** | HC50/hemolysis modeling and uncertainty |
| **S06** | ESM-2 embeddings, hybrid models, and QMAP evaluation |
| **S07** | Multi-task learning, species pooling, measured specificity |
| **S08** | Original 20 leads, design filters and applicability domain |
| **S09** | Fitness-objective ablations |
| **S10** | Ten-seed GA versus NSGA-II analysis |
| **S11** | Exploratory HydrAMP and reconstructed AMPGAN v2 comparisons |
| **S12** | Scripts, inputs and source provenance |

## Analysis interpretation

- MIC and HC50 values for original PACE leads are **model predictions**, not laboratory-confirmed activity.
- The independent AMPScanner comparison concerns **AMP classification** on retained sequences; it is not an assay hit-rate comparison.
- The NSGA-II and external-generator comparisons involve model scores and differing methodological conditions; do not interpret them as experimentally demonstrated superiority.
- Species-group pooling and assay heterogeneity are evaluated in supporting sensitivity analyses.
- Some subsequent analyses use reconstructed/refitted models instead of original saved fitted artifacts. The distinction is documented.

## Getting started

1. Review the [supplementary index](./supplementary/SUPPLEMENTARY_INDEX.csv) and the S01–S12 workbooks.
2. Consult [data sources and attribution](./DATA_SOURCES.md).
3. When the full scripts and data have been uploaded, follow the [code guide](./code/README.md) to reproduce supported analyses.
4. Use SHA-256 manifests to verify the version-specific scientific files.

## Citation and availability

**Associated manuscript:** *PACE: Pathogen-Aware Computational Engineering Evolutionary Framework for ESKAPE-Directed Antimicrobial Peptide Discovery and De Novo Design* (under revision).

**Repository:** https://github.com/Najeebkhan955/PACE-Pathogen-Aware-Computational-Engineering

No permanent archival DOI is claimed. Do not cite this repository as a complete code-and-data deposit until all listed materials are available.

## Research integrity and privacy

Reviewer correspondence and unpublished manuscript drafts are intentionally **not** included. External datasets, pretrained weights and third-party code retain the licenses and attribution requirements of their original sources.
