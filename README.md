# PACE | Pathogen-Aware Computational Engineering

**A computational evolutionary framework for ESKAPE-directed antimicrobial peptide discovery and de novo design**

[![Supplementary](https://img.shields.io/badge/Research%20data-S01%E2%80%93S12-164b74)](supplementary/) [![Reproducibility](https://img.shields.io/badge/Code-Python-227d69)](code/) [![Integrity](https://img.shields.io/badge/Integrity-SHA--256-227d69)](MASTER_MANIFEST_SHA256.csv)

## Project overview

PACE combines curated antimicrobial peptide data, handcrafted sequence representations, AMP classification, quantitative MIC and HC50 prediction, and genetic-algorithm candidate prioritization. Supporting computational analyses evaluate pretrained ESM-2 representations, multitask models, experimentally measured target-specific activity controls, species-pooling sensitivity, NSGA-II versus GA, and external AMP generators.

**The 20 final PACE sequences are computationally designed candidates. Predicted antimicrobial activity and hemolysis are not laboratory-validated outcomes.**

## Start here

| Resource | Contents |
|---|---|
| [S01–S12 supplementary Excel workbooks](supplementary/) | Twelve readable scientific summary/index workbooks |
| [Full supporting datasets](S01/) | S01–S12 organized `DATA/` folders with row-level CSV, FASTA, and figures |
| [Original PACE Python code](code/original_PACE/) | Author-supplied source scripts and execution order |
| [Additional analysis scripts](S12/DATA/) | Independent calculation, validation, notebook and verification code |
| [Reproducibility guide](docs/03_REPRODUCIBILITY.md) | File integrity and archive recovery |
| [Data provenance](DATA_SOURCES.md) | Source citations, license caveats and data lineage |
| [Supplementary index](SUPPLEMENTARY_INDEX.csv) | Dataset-level summary |
| [Visually formatted overview](PACE_Repository_Overview.xlsx) | Workbook for navigating supplemental results |

## Supplementary research groups

| Group | Area | Location |
|---|---|---|
| S01 | Source curation and provenance | [S01](S01/) |
| S02 | AMP classification and independent AMPScanner benchmark | [S02](S02/) |
| S03 | Measured target-specific activity and negative-control analyses | [S03](S03/) |
| S04 | Six ESKAPE MIC prediction endpoints and held-out results | [S04](S04/) |
| S05 | HC50/hemolysis estimates and uncertainty | [S05](S05/) |
| S06 | ESM-2, hybrid prediction and QMAP benchmarks | [S06](S06/) |
| S07 | Multitask learning and species-level MIC sensitivity | [S07](S07/) |
| S08 | Original 20 leads, FASTA, diversity and applicability domain | [S08](S08/) |
| S09 | GA objective ablation experiments | [S09](S09/) |
| S10 | Controlled GA vs NSGA-II comparison, all candidate evaluations | [S10](S10/) |
| S11 | HydrAMP and reconstructed AMPGAN v2 comparisons | [S11](S11/) |
| S12 | Reproducibility notebooks, scripts and source-code archive | [S12](S12/) |

## Reproduce and validate

Requires Python 3 for integrity checks. From the root of a cloned repository:

```bash
python code/verify_checksums.py
python code/verify_repository_manifest.py
python code/reassemble_archives.py --verify-only
python code/reassemble_archives.py
```

**Large original result archives are available losslessly** as `.zip.part001` etc, each approximately 2 MiB, in `S06/DATA/` and `S10/DATA/`. The reassembly script joins these chunks, verifies per-part SHA-256 and verifies the reconstructed ZIP SHA-256. The ZIP archives are not shortened.

## Reproducibility scope

The repository preserves available source code, numerical model comparisons and full retained computational analysis outputs. Some analyses rely on refitted or reconstructed models; the original saved fitted-model artifacts and exact historical GA run were not provided. Benchmarks using a common scoring model must not be interpreted as experimentally demonstrated antimicrobial hit rates.

## Data availability and citation

Public repository: **https://github.com/Najeebkhan955/PACE-Pathogen-Aware-Computational-Engineering**. See [citation guidance](CITATION.md). A permanent archived release DOI has not been assigned. The associated manuscript is under revision, not claimed published.

## Responsible reuse

The source datasets and third-party benchmark software have distinct owners and reuse terms. See [data sources](DATA_SOURCES.md) and [license notes](docs/04_ATTRIBUTION_AND_LICENSES.md). No blanket license is asserted. Reviewer comments, editor correspondence and manuscript drafts are excluded from this public repository.
