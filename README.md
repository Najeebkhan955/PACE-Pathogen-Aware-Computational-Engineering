# PACE: Pathogen-Aware Computational Engineering Evolutionary Framework

**ESKAPE-directed antimicrobial peptide discovery and de novo design**

PACE is an integrated four-stage computational framework for antimicrobial peptide (AMP) recognition, pathogen-associated minimum inhibitory concentration (MIC) prediction, erythrocyte hemolysis (HC₅₀) modeling, and genetic-algorithm (GA) sequence optimization. It combines predicted potency, modeled selectivity, similarity constraints, and physicochemical properties to prioritize experimentally testable peptide candidates.

**The 20 designed PACE peptides are computational candidates. Their MIC, HC₅₀, and selectivity estimates are model predictions, not new laboratory measurements.**

[![Supplementary](https://img.shields.io/badge/Research%20data-S01%E2%80%93S12-164b74)](supplementary/) [![Code](https://img.shields.io/badge/Research%20code-Python-227d69)](code/) [![Integrity](https://img.shields.io/badge/Data%20integrity-SHA--256-227d69)](MASTER_MANIFEST_SHA256.csv)

## Computational framework

1. **AMP recognition:** a soft-voting classifier uses 90 sequence-level descriptors (28 physicochemical, 20 amino-acid composition, and 42 composition–transition features) to estimate AMP probability.
2. **Pathogen-associated MIC prediction:** six separate quantitative regressors cover *Enterococcus faecium*-associated, *Staphylococcus aureus*, *Klebsiella pneumoniae*, *Acinetobacter baumannii*, *Pseudomonas aeruginosa*, and *Enterobacter cloacae*-associated endpoints. The *Enterococcus* and *Enterobacter* labels include specified pooled species and should not be interpreted as strain-specific predictions.
3. **Hemolysis and selectivity:** a quantitative model predicts HC₅₀. Computational selectivity is **cSI = predicted HC₅₀ / GM(predicted MIC₁, …, predicted MIC₆)**, where GM is the geometric mean.
4. **Sequence optimization:** the fixed-weight GA objective includes AMP probability (+2.0), mean predicted log₁₀ MIC (−1.0), worst-target predicted log₁₀ MIC (−0.5), bounded log-selectivity (+0.6), sequence-similarity penalty (−2.0), and charge-deviation penalty (−0.35). Post-search screening and diversity selection yield 20 candidate peptides.

## Main data and results

| Item | Reported finding |
|---|---|
| AMP recognition corpus | 7,839 curated AMPlify sequences |
| External recognition corpus | 1,848 leakage-filtered amPEPpy sequences |
| Quantitative activity corpus | 44,151 curated GRAMPA MIC observations |
| Hemolysis corpus | 1,271 peptide-level HC₅₀ records |
| Held-out soft-voting AMP classifier | ROC–AUC **0.972** |
| External soft-voting AMP classifier | ROC–AUC **0.915** on the full independent external cohort |
| Six quantitative MIC models | Cross-validated R² **0.248–0.584** |
| HC₅₀ regression | Reported R² **0.432** |
| Designed PACE library | **20 peptides**, 28–30 residues per candidate |

On a distinct comparator-overlap-filtered external cohort, the 90-descriptor Extra-Trees, ESM-2, hybrid, and pretrained AMPScanner models achieved ROC–AUC values of **0.912, 0.932, 0.932, and 0.803**, respectively. These benchmark results have a different evaluation cohort from the full external soft-voting ensemble.

## Computational evaluation

The study includes independent AMP recognition and pretrained-model comparisons, ESM-2 representation benchmarks, experimentally measured target-specific strong/weak activity controls, multitask and species-aware sensitivity analyses, objective ablations, matched-budget GA versus NSGA-II comparisons, HydrAMP and AMPGAN v2 exploratory library scoring, and applicability-domain and prediction-uncertainty assessments. The ten-seed GA versus NSGA-II comparisons did not identify an overall BH-FDR-significant winner.

Measured activity-control observations originate from public source datasets; scores assigned to newly designed PACE peptides remain model-derived. Comparative common-surrogate scores are **not experimental antimicrobial hit rates**.

## Research materials

| Resource | Content |
|---|---|
| [S01–S12 supplementary workbooks](supplementary/) | Twelve navigable scientific collections |
| [Data collection index](SUPPLEMENTARY_INDEX.csv) | Group titles, counts, and workbook paths |
| [Data provenance](DATA_SOURCES.md) | Input origins, target definitions, and evidence interpretation |
| [Source Python scripts](code/original_PACE/) | Original algorithmic workflow and execution order |
| [Supporting analysis scripts](S12/DATA/) | Numerical evaluations and file integrity |
| [Reproducibility guide](docs/03_REPRODUCIBILITY.md) | Dataset access, SHA-256 checks, and large-archive assembly |
| [Citation guidance](CITATION.md) | Repository and manuscript identification |

### Supplementary collection overview

| Group | Content |
|---|---|
| [S01](S01/) | Peptide source curation and provenance |
| [S02](S02/) | AMP recognition and external classifier comparison |
| [S03](S03/) | Measured strong-versus-weak antimicrobial activity |
| [S04](S04/) | Six ESKAPE-associated MIC endpoints |
| [S05](S05/) | Hemolysis prediction and uncertainty |
| [S06](S06/) | ESM-2, composition–length, hybrid, and QMAP comparisons |
| [S07](S07/) | Multitask, pooled-species, and cross-pathogen analyses |
| [S08](S08/) | The 20 PACE leads, sequence diversity, and applicability domain |
| [S09](S09/) | GA fitness-component ablations |
| [S10](S10/) | Matched-budget GA versus NSGA-II |
| [S11](S11/) | HydrAMP and AMPGAN v2 exploratory scoring |
| [S12](S12/) | Research code and reproducibility records |

## File integrity and reproducibility

With Python 3 from the repository root:

```bash
python code/verify_checksums.py
python code/verify_repository_manifest.py
python code/reassemble_archives.py --verify-only
```

Complete large archives in S06 and S10 are distributed in verifiable parts; run `python code/reassemble_archives.py` to assemble them. Available scripts and preserved output records support the documented analyses, although exact replay of the initial GA trajectory may require fitted-model artifacts not included in the deposited materials.

## Data availability and citation

Repository: **https://github.com/Najeebkhan955/PACE-Pathogen-Aware-Computational-Engineering**. Cite the relevant commit or tagged release and the associated research paper once a publication citation is available. No permanent repository DOI is currently listed. Consult [source-specific licensing information](docs/04_ATTRIBUTION_AND_LICENSES.md) before redistributing third-party records or code.
