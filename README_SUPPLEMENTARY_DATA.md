# PACE — Supplementary Data S01–S12: Author/Editor Index

**This is an analysis-based supplementary data draft, not a new manuscript, a revised reviewer response, or experimental AMP validation.** All author source files are unchanged.

## Structure

Each of 12 groups contains its own styled XLSX navigation/results workbook, complete CSV records (no large-data truncation), original FASTA where available, and data source provenance. Any shorter Excel sample clearly states it is a view rather than a replacement for the full CSV. All supporting CSV records are included with the same values as their verified source files. Full 10-seed GA/NSGA-II candidate archives are stored losslessly as six Git-tracked binary ZIP parts under `S10/DATA/`; original bytes are reconstructed with `python code/reassemble_archives.py`.

**Groups:**

- **S01 — Source provenance and curation.** Data sources, original curated records, original input metadata and transparent curation checks.
- **S02 — AMP classification and independent external predictor.** AMP labeling, held-out predictions and overlap-filtered AMPScanner v2 comparator.
- **S03 — Measured potency and weak-activity negative controls.** Experimental, species-specific strong- versus weak-activity peptide observations and composition-matched controls.
- **S04 — Six ESKAPE MIC regression models.** Complete independently curated target MIC, original-model and extended OOF prediction evidence.
- **S05 — Hemolysis prediction and uncertainty.** Full curated hemolysis records, model residuals, similarity-group validation and uncertainty.
- **S06 — Protein language models and QMAP transfer.** Original 90D, 21D, pretrained ESM-2, hybrid comparisons and QMAP benchmark OOF predictions.
- **S07 — Multi-task species pooling and cross-pathogen evidence.** Masked multi-task OOF outputs, group-held-out species pooling, cross-pathogen paired measured MIC.
- **S08 — Original 20 leads and applicability domain.** Complete lead sequences, original source predictions, seven-target applicability domain and ranking diagnostics.
- **S09 — Evolutionary design and objective ablation.** Complete reoptimization diagnostics for full, no-hemolysis, no-novelty and reduced fitness terms.
- **S10 — Controlled fixed-weight GA versus NSGA-II.** Matched 9000-unique evaluations per seed, 10 seeds per algorithm, 40 selected lead sets and full candidate archives.
- **S11 — External generated libraries and comparative rescoring.** Published HydrAMP candidates, original-checkpoint reconstructed AMPGAN v2 peptides, PACE leads and length sensitivity.
- **S12 — Data availability and reproducible code.** Original PACE scripts, versioned source feature engine, independent metric calculators, run order and checksums.

## Interpretation and submission safeguards

1. Distinguish measured MIC/HC50 source data from **model-predicted** lead values throughout. The 20 original leads are not experimentally validated.
2. A high MIC for one pathogen is a target-specific weak-activity control, not evidence of universal biological inactivity.
3. AMPGAN v2 sequences came from original weights and a reconstructed TensorFlow forward pass; exact historical native inference equivalence remains unverified.
4. The fixed-weight GA and NSGA-II analyses are **new controlled/reconstructed** comparisons; do not describe their sequences as the original historical GA output.
5. QMAP is a different species/data benchmark and does not directly prove ESKAPE MIC generalization.
6. The original seven input files, including the original Reviewer 1 and Reviewer 2 comments, were not edited. The reviewer comments are not reproduced in this publication-focused package.
7. Original research figures, table numbering, references, and all supporting file citations must be reconciled during the subsequent manuscript revision. Verify the target journal's allowed supplement formats, third-party source reuse permissions and repository deposit requirements before submission.

## Provenance

- `MASTER_MANIFEST_SHA256.csv`: content hashes, source file paths, line counts and group affiliation.
- `ORIGINAL_AUTHOR_INPUTS_UNCHANGED_SHA256.csv`: immutable author-input checksums, **without** embedding private reviewer documents.
- `DATA_SOURCE_PROVENANCE.csv`: public source URLs and provenance limitations.
- `SUPPLEMENTARY_INDEX.csv`: navigation and total input data records.

## Reproducibility

- S12 contains original scripts, additional portable analysis scripts, and manifests. Original fitted PACE model pickles and complete historical GA trajectory were never provided; do not imply exact original-run reproduction.
- Data and figure provenance is preserved using unchanged copies of analyzed outputs.

## Suggested assembly workflow

1. Freeze and verify the sources and compare every included CSV header/row count with the manuscript tables.
2. Synchronize manuscript tables/figures and journal numbering to S01–S12 using the actual results, without changing original reviewer statements.
3. Public data/code are deposited at https://github.com/Najeebkhan955/PACE-Pathogen-Aware-Computational-Engineering; no permanent archived DOI is assigned.
4. Finalize manuscript and author responses only after scientific cross-file consistency checks.
