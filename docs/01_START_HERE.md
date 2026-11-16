# PACE | Guide to the research repository

The PACE repository provides a navigable record of the computational framework, its data sources, model development, peptide-design results and supporting evaluations.

## Recommended entry points

1. Read the [research overview](../README.md) for the framework and key findings.
2. Browse the [supplementary workbooks](../supplementary/) for scientific summary tables.
3. Open the corresponding numbered [S01–S12 data folders](../SUPPLEMENTARY_INDEX.csv) for complete records.
4. Follow the [source Python scripts](../code/original_PACE/) and [run order](../code/original_PACE/RUN_ORDER.md) for the original computational pipeline.
5. Use the [reproducibility guide](03_REPRODUCIBILITY.md) to verify checksums and restore multipart archives.

## Locate the supporting evidence

| Research question | Collections |
|---|---|
| How were source peptides and labels curated? | [S01](../S01/), [S02](../S02/) |
| How was activity against individual target groupings modeled? | [S03](../S03/), [S04](../S04/), [S07](../S07/) |
| How were hemolysis and selectivity modeled? | [S05](../S05/) |
| How were pretrained representations evaluated? | [S06](../S06/) |
| What are the exact sequences and modeled characteristics of the 20 candidates? | [S08](../S08/) |
| How do evolutionary objectives and alternatives compare? | [S09](../S09/), [S10](../S10/), [S11](../S11/) |
| Where are computational scripts and integrity records? | [S12](../S12/), [code](../code/) |

Each collection contains scientific context and its underlying numerical files. For consistent interpretation, distinguish measured source activity from model-predicted characteristics of newly designed peptides.
