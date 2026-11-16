# PACE | Data dictionary and reporting conventions

This guide defines commonly used variables across the PACE research collections. Individual workbook codebooks and source scripts provide the exact column-level definitions, preprocessing decisions, and evaluation protocols.

| Variable or concept | Interpretation |
|---|---|
| Peptide sequence | Canonical amino-acid sequence used as a molecular or training-set identifier |
| Length | Number of amino-acid residues in a peptide |
| Estimated net charge | Sequence-derived charge under the specified pH and pKₐ assumptions |
| MIC (µM) | Minimum inhibitory concentration from a source experiment or model prediction, as indicated |
| log₁₀ MIC (µM) | Base-10 transformed concentration; the original µM scale is recovered by exponentiation |
| HC₅₀ (µM) | Concentration associated with 50% erythrocyte hemolysis, measured or predicted as labeled |
| log₁₀ HC₅₀ (µM) | Base-10 logarithm of HC₅₀ concentration |
| Computational selectivity ratio (cSI) | Predicted HC₅₀ divided by the geometric mean of six predicted MIC values |
| ROC–AUC | Area under the receiver operating characteristic curve for specified classification cohorts |
| PR–AUC | Area under the precision–recall curve, interpreted with the evaluation prevalence |
| R² | Coefficient of determination for quantitative prediction |
| RMSE and MAE | Root-mean-square error and mean absolute error, with endpoint-specific units |
| OOF | Out-of-fold prediction under a declared cross-validation protocol |
| GA | Weighted genetic-algorithm sequence optimization |
| NSGA-II | Nondominated sorting genetic algorithm II, used in comparative optimization |
| Applicability domain | Descriptor-space proximity to the source training data, evaluated using a defined distance threshold |
| 3-mer Jaccard similarity | Similarity between sets of contiguous amino-acid triplets |
| BH-FDR | Benjamini–Hochberg false-discovery-rate adjustment |

## Record structure

A row may represent a sequence, a peptide–pathogen observation, a measurement replicate, a model prediction, a resampled evaluation, or an evolutionary candidate. Consequently, row counts across different output files should not be treated as numbers of distinct biological samples without inspecting the dataset definition.

**Evidence labels:** experimental source MIC/hemolysis values and computational estimates for the designed PACE peptides are separate evidence types. Consult [S01–S12](../supplementary/) and [data provenance](../DATA_SOURCES.md) for their origins and appropriate interpretation.
