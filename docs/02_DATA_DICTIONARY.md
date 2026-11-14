# Dataset interpretation guide and reporting units

This guide identifies core *classes* of variables. The exact column-level definitions, units, and original selection criteria must be taken from the relevant supplementary workbook, group README, and source code.

| Field family | Interpretation | Caution |
|---|---|---|
| `sequence` / `peptide` | Amino-acid sequence identifier or canonical peptide sequence | Noncanonical residues and modifications may have separate curation treatment |
| `MIC` | Minimum inhibitory concentration, typically reported in µM for these analyses | Source study strain, conditions and transformation must be checked |
| `log10_MIC` | Log10-transformed MIC concentration, when explicitly labeled | Transform back using 10 to the power of the stored value |
| `HC50` | Hemolytic concentration, when stated, typically in µM | Censoring at source assay limits may affect comparisons |
| `log10_HC50` | Log10-transformed HC50 when explicitly labeled | A model prediction is not an assay result |
| `ROC_AUC` | Receiver operating characteristic area under the curve | Recomputed on a defined evaluation cohort; not comparable across different cohorts without care |
| `PR_AUC` | Precision–recall area | Depends strongly on case prevalence and evaluation definition |
| `OOF` | Out-of-fold predictions | Check whether splits are sequence-grouped / similarity-aware |
| `GA` / `NSGA-II` | Genetic algorithm / Pareto evolutionary baseline | Model-predicted fitness is not observed antimicrobial effect |
| `AD` | Applicability domain | Proximity to observed descriptor space is not a prediction guarantee |

**Record conventions:** Different CSVs may represent individual peptide records, peptide–pathogen pairs, experimental measurements, repeated cross-validation split rows, or evolutionary candidate evaluations. Therefore **rows should not be summed as independent biological observations** across files. Exact columns should be inspected rather than inferred from this guide.
