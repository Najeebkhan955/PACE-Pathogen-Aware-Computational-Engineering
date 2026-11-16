# PACE | Pathogen-Aware Computational Engineering Evolutionary Framework

### Computational discovery and *de novo* design of antimicrobial peptides against ESKAPE-associated pathogens

PACE combines machine-learning prediction with evolutionary sequence design to prioritize antimicrobial peptide (AMP) candidates across multiple bacterial targets. The study integrates AMP recognition, quantitative minimum inhibitory concentration (MIC) modeling, erythrocyte hemolysis prediction, and sequence optimization within a reproducible, data-driven research workflow.

**Research contribution:** an interpretable design objective combining overall and worst-target predicted antimicrobial potency, modeled hemolytic selectivity, reference-sequence similarity, and sequence diversity to generate a focused library of **20 candidate peptides**.

## Study highlights

| Research component | Main finding |
|---|---|
| Curated AMP classification dataset | **7,839** AMPlify sequences |
| Independent external AMP evaluation | **1,848** amPEPpy sequences after overlap filtering |
| Quantitative antimicrobial activity | **44,151** curated GRAMPA MIC observations |
| Hemolysis dataset | **1,271** unique peptide records |
| AMP recognition | **0.972 ROC–AUC** on the held-out AMPlify set |
| External classification | **0.915 ROC–AUC** for the soft-voting ensemble |
| Six pathogen-associated MIC models | Cross-validated **R² = 0.248–0.584** |
| Hemolysis modeling | Reported **HC₅₀ R² = 0.432** |
| Final computational design output | **20 peptides, 28–30 amino acids** |

These results describe specified predictive evaluations; activity and selectivity values for the 20 designed peptides are **computational predictions**, not measurements from prospective laboratory tests.

## PACE computational architecture

1. **AMP recognition.** Ninety sequence-level descriptors (28 physicochemical, 20 amino-acid composition, and 42 composition–transition descriptors) support a soft-voting AMP classifier.
2. **Quantitative antimicrobial activity.** Six separate MIC regressors model *Enterococcus faecium*-associated, *Staphylococcus aureus*, *Klebsiella pneumoniae*, *Acinetobacter baumannii*, *Pseudomonas aeruginosa*, and *Enterobacter cloacae*-associated endpoints. The first and last endpoint groupings incorporate specified additional species.
3. **Hemolytic selectivity.** An HC₅₀ model estimates erythrocyte hemolysis. The computational selectivity ratio is the predicted HC₅₀ divided by the geometric mean of six predicted MIC values.
4. **Evolutionary sequence selection.** A weighted genetic algorithm balances AMP probability, mean and worst-target predicted MIC, modeled selectivity, novelty/similarity penalties, and charge deviation, followed by eligibility screening and diversity-constrained selection.

### Computational evaluation

The research evaluates pretrained ESM-2 representations, AMPScanner v2, composition–length baselines, measured target-specific activity controls, multitask and pooled-species MIC models, GA objective ablations, matched-budget NSGA-II comparisons, generative peptide libraries (HydrAMP and AMPGAN v2), descriptor applicability domains, and prediction uncertainty.

The matched external representation benchmark reported ROC–AUC **0.912** (90D Extra-Trees), **0.932** (ESM-2), **0.932** (hybrid), and **0.803** (pretrained AMPScanner). These results use a different eligibility-filtered cohort from the full external soft-voting ensemble.

## Research materials and downloads

| Resource | Description |
|---|---|
| [Supplementary collections S01–S12](supplementary/) | Downloadable scientific workbooks with summary tables and results |
| [Full supplementary datasets](S01/) | Numbered research folders, each with its supporting `DATA/` content |
| [Dataset index](SUPPLEMENTARY_INDEX.csv) | Collection titles, counts, and workbook locations |
| [Source data and provenance](DATA_SOURCES.md) | Dataset origins, biological endpoint definitions, and interpretation |
| [PACE Python pipeline](code/original_PACE/) | Curation, model development, evolutionary search and figure scripts |
| [Analysis and verification code](S12/DATA/) | Computational benchmarks, numerical analysis and research records |
| [Reproducibility guide](docs/03_REPRODUCIBILITY.md) | Execution context, integrity verification and archive restoration |
| [Citation guidance](CITATION.md) | Study identity, version-specific referencing and attribution |

### Supplementary collections at a glance

| Collection | Research focus |
|---|---|
| [S01](S01/) | Data sources and peptide curation |
| [S02](S02/) | AMP recognition and external classification |
| [S03](S03/) | Measured peptide–pathogen activity controls |
| [S04](S04/) | Six MIC regression endpoints |
| [S05](S05/) | Hemolysis prediction and uncertainty |
| [S06](S06/) | ESM-2, hybrid and QMAP benchmarks |
| [S07](S07/) | Multitask and species-associated MIC analyses |
| [S08](S08/) | Designed peptides, diversity and applicability domains |
| [S09](S09/) | GA objective-component analyses |
| [S10](S10/) | GA versus NSGA-II comparison |
| [S11](S11/) | HydrAMP and AMPGAN v2 comparisons |
| [S12](S12/) | Computational scripts and data-integrity resources |

## Reproducibility

Use Python 3 from a local clone of this repository:

```bash
python code/verify_checksums.py
python code/verify_repository_manifest.py
python code/reassemble_archives.py --verify-only
```

Large numerical archives in S06 and S10 are preserved as verifiable multipart files. To reconstruct them, run `python code/reassemble_archives.py`. The repository documents model configurations, retained results, and scripts; an exact replay of the initial evolutionary search may additionally require fitted-model artifacts unavailable in the deposit.

## Citation and data use

**Research repository:** https://github.com/Najeebkhan955/PACE-Pathogen-Aware-Computational-Engineering

Cite the exact Git commit or tagged release alongside the study publication when available. A permanent archival DOI is not listed. The research uses third-party datasets and software with source-specific reuse terms; please consult the [citation guidance](CITATION.md), [provenance](DATA_SOURCES.md) and [licensing notes](docs/04_ATTRIBUTION_AND_LICENSES.md).

*The 20 PACE-designed peptides are computational research candidates. Their predicted MIC, HC₅₀, and selectivity have not been prospectively confirmed by laboratory testing.*
