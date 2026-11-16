# S02 — AMP classification and independent external predictor

AMP classifier development, held-out evaluation, and external predictor benchmarks.

**Data organization:** The supplementary workbook provides scientific summary and index tables; deposited row-level numerical data, figure files, and scripts are in `DATA/` where applicable. File names, accession provenance, and checksums are preserved for traceability.

**Scientific interpretation:** Classifier performance is evaluated on the specified training, held-out, and independent external cohorts.

## Complete files

| File | Rows (where applicable) | Source archive | SHA256 |
|---|---:|---|---|
| composition_matched_auc_diagnostic.csv | 8 | pace_phase12_audit/composition_matched_auc_diagnostic.csv | `b96eedd7c45ace549376cfb10df72bac7d4f73debec728cf05249feee0e42b6d` |
| independently_recomputed_classification_metrics.csv | 30 | pace_phase14b_verified/independently_recomputed_classification_metrics.csv | `2497d1f091e976d2d0bdaf72f26eae1f5bf725a1fca6ee40f1d3cfb9220a0ea4` |
| independent_stratified_paired_bootstrap.csv | 24 | pace_phase14b_verified/independent_stratified_paired_bootstrap.csv | `e9af7d60af76a31af11578215d9a078e577579dfb6b1dc74d9c50aa70e91999b` |
| evaluated_and_excluded_counts.csv | 2 | pace_phase14b_verified/evaluated_and_excluded_counts.csv | `74766779840baba0e8096bc406850a81d396937a4e42b25efe58f64597d2ffcc` |
| phase14b_AMPScanner_evaluation_manifest.csv | 3416 | pace_phase14b_analysis/phase14b_AMPScanner_evaluation_manifest.csv | `ab6c47257f5bee73e56dd02dae29cd574a41602e460091a54ca6fe2ae8e72bb1` |
| classification_predictions.csv | 27328 | PACE_Phase8_ESM2_results.zip::classification_predictions.csv | `e08d43e40934ab1ba63d983029e6dfdc8f5e3fbd01f26860e2aa4ced8eb4279a` |
| classification_metrics.csv | 16 | PACE_Phase8_ESM2_results.zip::classification_metrics.csv | `0554836dd5278f96a53f38b416add929ab3ce6741a6e5d830296d866f77e6953` |
| AMPScanner_PACE_matched_predictions.csv | 3050 | PACE_Phase14B_AMPScanner_External_Results.zip::AMPScanner_PACE_matched_predictions.csv | `5fbd2a37b68b2cd5a2cb27c48053943f5ed638dbdcd1a62ced6433f7a450795b` |
| S02_Supplementary_Tables.xlsx |  | Generated with artifact_tool from original verified CSV files | `ab9b4ea8c33ac84b5ca8e0bb93609022c043b21651561d5f9141abd53480eafa` |
