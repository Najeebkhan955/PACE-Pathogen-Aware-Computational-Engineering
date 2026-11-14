# PACE Phase 11E5 independent verification & Phase 11F exploratory comparative audit

**Scope:** ANALYSIS ONLY. **Status:** Reconstruction completed; numerical and sequence checks passed; historical native-TensorFlow forward equivalence **not verified**. Original PACE manuscript, supplementary data, reviewer comments and reproducibility code are **unchanged** (all six saved original-input hashes match prior Phase 10B snapshot).

## A. Exact evidence from user-uploaded Phase 11E5 ZIP

- ZIP CRC/structural integrity: **passed**. The uploaded archive contains a complete reconstructed inference run, source variable inventory, FASTA and CSV outputs, run log and status manifest.
- **2,024 generated outputs:** **1,000 AMP-positive-conditioned** and **1,024 non-AMP-conditioned** rows. The latter 1,024 must not enter the AMP-generation comparison.
- All **1,000/1,000 AMP-positive-conditioned sequences use canonical 20-amino-acid symbols**, and all **1,000 sequences are unique**. FASTA matched positive rows of CSV in exact order with no mismatch. None exactly matches the **12,979-sequence PACE Phase 8 reference**.
- Recalculated positive lengths: **mean 17.491, median 17.0, range 1–32 aa**. **197/1,000 (19.7%) are <12 aa** and therefore outside original PACE search range. **649/1,000** lie within 12–25 aa; **793/1,000** within 12–30 aa.
- Reconstructed generator uses the **original v2.0 architecture and 16 restored checkpoint weight tensors**. All sixteen shapes reported in the archived inventory match those expected by the generator, and generator source-section SHA-256 matches original `ampgan/model.py` recovered in Phase 11E3. Official `PairedAMPSequence` source explicitly concatenates positive conditions before negative conditions, validating the positive-half label assignment.
- The archive **does not contain original checkpoint weight bytes**; their previous SHA-256 validation is supported by earlier Colab logs, rather than independently repeated here.
- **Critical provenance:** Reconstruction implemented the published/original layer algebra with direct TensorFlow ops to bypass an obsolete `CoordinateChannel1D` SavedModel placeholder. Exact output equivalence with AMPGAN v2's native historical TF graph has **not** been tested. Label as *AMPGAN v2 reconstructed-from-official-weights candidates*; **not** confirmed native-generation samples, not ESKAPE-specific optimized candidates, and not experimentally active AMPs.

## B. Sequence diversity, descriptors, and common PACE quantitative-model rescoring

| Candidate library | N | Mean length (aa) | Mean net charge (pH 7.4) | Median within-library 3-mer Jaccard ↓ | Mean max 3-mer Jaccard vs same 12,979-reference ↓ | Mean **predicted** geometric MIC (µM) ↓ | Mean **predicted** worst-target MIC (µM) ↓ | Mean modeled SI ↑ | Predicted joint threshold pass (MIC max ≤32 µM, SI ≥10) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| PACE original top-20 | 20 | 29.45 | 8.90 | 0.068 | 0.101 | 2.37 | 3.80 | 68.11 | 20/20 |
| HydrAMP released top-100 | 100 | 20.62 | 3.63 | 0.000 | 0.148 | 17.07 | 31.62 | 9.46 | 32/100 |
| AMPGAN v2 reconstructed AMP-positive (all) | 1000 | 17.49 | 4.52 | 0.000 | 0.176 | 19.50 | 36.55 | 15.17 | 441/1000 |
| AMPGAN v2 reconstructed AMP-positive (12–25 aa) | 649 | 17.71 | 4.63 | 0.000 | 0.170 | 16.97 | 30.59 | 15.89 | 314/649 |
| PACE post-hoc 12–25 (10 seeds) | 200 | 24.51 | 7.82 | 0.022 | 0.124 | 3.58 | 5.73 | 48.19 | 200/200 |

**Definitions and limitations:** 3-mer Jaccard 0 means no shared 3-mers (not necessarily no homologous sequence); median Jaccard can be 0 in diverse libraries. Short sequences can have identical k-mer *sets* even when sequences differ: full AMPGAN set's maximum pairwise Jaccard is 1 despite all 1,000 primary sequences being distinct. No exact reference match does not prove evolutionary/mechanistic novelty. HydrAMP input is the **released starter-kit top-100** library, not necessarily the precise peptides experimentally tested in its original publication.

The seven quantitative predictors (six MIC endpoints and HC50) were refitted as 300-tree ExtraTrees models to **the same original PACE quantitative labels/90-descriptor pipeline**. For the original PACE-20 + HydrAMP-100 sequences, refitted predicted geometric MIC agrees with independently saved Phase 11 results to **max absolute deviation 4.6e-14 µM**. This confirms implementation consistency. The MIC and selectivity figures are **model predictions only**, not measured MIC/HC50 or true screening hit rates. The simple threshold does **not** include AMP-classifier acceptance.

**Major interpretive limitation — evaluator circularity:** PACE peptides were optimized by prediction models of this same family, whereas external generated libraries were not. The reward-like evaluator is therefore predictably favorable to PACE. Additionally, libraries differ in length, intended task, selection stringency and candidate counts; no unbiased model-superiority inference or p-value is warranted. Earlier independent QMAP assessment found HC50 transfer performance to be weak; modeled selectivity is especially uncertain.

## C. Common 22–25 aa length-band sensitivity (not a matched experimental design)

| Library | N within 22–25 aa | Mean actual length | Mean **predicted** MIC (µM) | Mean **predicted** worst-target MIC (µM) | Mean modeled SI |
|---|---:|---:|---:|---:|---:|
| PACE post-hoc 12–25 (10 seeds) | 200 | 24.51 | 3.58 | 5.73 | 48.19 |
| HydrAMP released top-100 | 40 | 23.70 | 15.40 | 30.64 | 10.42 |
| AMPGAN v2 reconstructed AMP-positive (all) | 121 | 23.78 | 12.34 | 20.69 | 17.29 |

This helps expose length effects but is **not** exact per-length matching, nor does it remove optimization-selection bias. PACE here is the **post-hoc selected** 200 sequences from 10 earlier evolutionary archives, while HydrAMP and reconstructed AMPGAN sequences are published/unaltered generation outputs, respectively. The original PACE 20 (length 28–30) have *no* sequences in this shared band.

## D. Verdict and next evidence needed

**What is now complete:** (i) the reconstructed AMPGAN v2 sampling output exists, (ii) FASTA, CSV, canonicality, uniqueness, source-architecture consistency and length checks pass, (iii) same-source PACE/HydrAMP baseline predictions reproduce, and (iv) an **exploratory three-library** dataset with transparent lengths, novelty/diversity and reward-model rescoring is available. Repeated full checkpoint downloads are unnecessary *for this exploratory descriptive analysis*.

**What remains unresolved:** (1) native historical generator equivalence to the reconstructed forward pass; (2) homology-aware **independent** ESKAPE MIC/HC50 evaluation (QMAP validation alone is not a directly comparable six-ESKAPE predictor for the generated sequences); (3) length/budget/selection/condition-matched generation under genuinely comparable objectives; (4) wet-lab MIC/HC50 if biological superiority is to be claimed.

**Defensible reviewer response at present:** We evaluated an external released HydrAMP library and outputs reconstructed from the released AMPGAN v2 checkpoint; comparisons reveal pronounced differences in length/charge/diversity and exploratory PACE-model scores. Because of evaluator circularity, different optimization objectives, and unconfirmed native equivalence for reconstructed AMPGAN outputs, this analysis does not establish superior antimicrobial potency or safety of PACE. Manuscript changes remain prohibited pending consolidated final reviewer audit and explicit user go-ahead.

## Files / validation
- `phase11f_library_summary.csv` and JSON: independent library metrics.
- `phase11f_all_sequence_model_rescoring.csv`: per-peptide seven predictor outputs.
- `phase11f_sequence_properties_and_novelty.csv`: descriptor and reference-overlap checks.
- `phase11f_common_length_22to25_sensitivity.csv`: length-band sensitivity.
- `AMPGANv2_reconstructed_positive1000.fasta`: verbatim from user ZIP; reconstructed-status, not verified native equivalence.
- `phase11f_original_source_integrity.csv`: six user-provided input-file SHA-256 rechecks; all unchanged.
- `phase11f_verification_manifest.json`: precise independent verification scope.
- `phase11f_analysis_run.py`: source analysis script, reconstructs outputs using user ZIPs.

**No original manuscript, supplement, reviewer response, or submitted source file was edited.**
