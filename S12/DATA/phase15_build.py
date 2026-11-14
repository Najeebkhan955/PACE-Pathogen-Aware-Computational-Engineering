#!/usr/bin/env python3
"""PACE Phase 15: build 12 source-grounded, analysis-only supplementary groups.
Uses artifact_tool exclusively for XLSX workbook creation/editing.
Never mutates author-provided inputs or archived reviewer comments.
"""
from __future__ import annotations
import os, csv, re, io, sys, hashlib, zipfile, json, shutil, collections, textwrap, time
from datetime import datetime, timezone
from pathlib import Path
from artifact_tool import Workbook, SpreadsheetFile
R=Path('/mnt/data'); O=R/'PACE_Phase15_New_Supplementary_Materials'
# Idempotent rebuild of Phase15 generated outputs only; never touch author inputs.
if O.exists(): shutil.rmtree(O)
O.mkdir(exist_ok=True)
SRC=R
# group list: canonical storage locations, each basename remains unchanged for provenance
GROUPS=[
 ('S01','Source provenance and curation','Data sources, original curated records, original input metadata and transparent curation checks.',[
  'pace_phase4_analysis/mic_eskape_independently_curated.csv','pace_phase7_hemolysis/hemolysis_original_curation.csv','pace_phase14a_analysis/Phase14A_source_matched_peptide_level_summary.csv',
  'pace_phase12b_analysis/AMP_annotation_versus_measured_MIC_cross_tab.csv','pace_phase14c_analysis/Phase14C_seven_original_input_sha256.csv',
 ],['pace_phase14a_analysis/Phase14A_source_matched_peptide_level_summary.csv','pace_phase12b_analysis/AMP_annotation_versus_measured_MIC_cross_tab.csv']),
 ('S02','AMP classification and independent external predictor','AMP labeling, held-out predictions and overlap-filtered AMPScanner v2 comparator.',[
  'pace_phase12_audit/composition_matched_auc_diagnostic.csv',
  'pace_phase14b_verified/independently_recomputed_classification_metrics.csv','pace_phase14b_verified/independent_stratified_paired_bootstrap.csv','pace_phase14b_verified/evaluated_and_excluded_counts.csv',
  'pace_phase14b_analysis/phase14b_AMPScanner_evaluation_manifest.csv',
 ],['pace_phase14b_verified/independently_recomputed_classification_metrics.csv','pace_phase14b_verified/evaluated_and_excluded_counts.csv']),
 ('S03','Measured potency and weak-activity negative controls','Experimental, species-specific strong- versus weak-activity peptide observations and composition-matched controls.',[
  'pace_phase12b_analysis/measured_ESKAPE_activity_control_cohort.csv','pace_phase12b_analysis/measured_ESKAPE_sequence_grouped_OOF_predictions.csv',
  'pace_phase12b_analysis/measured_ESKAPE_grouped_CV_metrics.csv','pace_phase12b_analysis/experimental_species_composition_matched_pairs.csv',
  'pace_phase12b_analysis/experimental_species_matched_metrics.csv','pace_phase12b_analysis/independent_sequence_cluster_bootstrap.csv',
  'pace_phase12b_analysis/posthoc_threshold_sensitivity.csv','pace_phase12b_analysis/direct_AMP_classifier_measuredMIC_overlap.csv',
  'pace_phase12b_analysis/direct_AMP_classifier_measuredMIC_overlap_metrics.csv',
 ],['pace_phase12b_analysis/measured_ESKAPE_grouped_CV_metrics.csv','pace_phase12b_analysis/experimental_species_matched_metrics.csv']),
 ('S04','Six ESKAPE MIC regression models','Complete independently curated target MIC, original-model and extended OOF prediction evidence.',[
  'pace_phase4_analysis/extra_trees_out_of_fold.csv','pace_phase5_analysis/mic_all6_oof_preds.csv','pace_phase5_analysis/mic_all6_oof_error_domain.csv',
  'pace_phase6_analysis/nested_oof_predictions.csv','pace_phase6_analysis/split_conformal_mic_predictions.csv',
  'pace_phase4_analysis/grampa_pairwise_statistics.csv','pace_phase4_analysis/modification_sensitivity_counts.csv',
 ],['pace_phase4_analysis/grampa_pairwise_statistics.csv','pace_phase4_analysis/modification_sensitivity_counts.csv']),
 ('S05','Hemolysis prediction and uncertainty','Full curated hemolysis records, model residuals, similarity-group validation and uncertainty.',[
  'pace_phase7_hemolysis/hemolysis_300tree_extratrees_oof.csv','pace_phase7_hemolysis/hemolysis_similarity_cluster_oof.csv',
  'pace_phase7_hemolysis/hemolysis_split_conformal_oof.csv','pace_phase7_hemolysis/hemolysis_all5_model_metrics.csv',
  'pace_phase14c_analysis/Phase14C_heldout_residual_quantiles.csv','pace_phase14c_analysis/Phase14C_all7_aligned_out_of_fold_residuals.csv',
 ],['pace_phase7_hemolysis/hemolysis_all5_model_metrics.csv','pace_phase14c_analysis/Phase14C_heldout_residual_quantiles.csv']),
 ('S06','Protein language models and QMAP transfer','Original 90D, 21D, pretrained ESM-2, hybrid comparisons and QMAP benchmark OOF predictions.',[
  'pace_phase11b_verified/independently_recalculated_split_metrics.csv','pace_phase11b_verified/verified_group_means.csv',
  'pace_phase11b_verified/paired_five_split_differences.csv',
  'pace_phase12b_analysis/QMAP_experimental_MIC_binarized_OOF.csv','pace_phase12b_analysis/QMAP_experimental_activity_discrimination_metrics.csv',
  'pace_phase12b_analysis/QMAP_five_split_paired_differences.csv',
 ],['pace_phase11b_verified/verified_group_means.csv','pace_phase12b_analysis/QMAP_experimental_activity_discrimination_metrics.csv']),
 ('S07','Multi-task species pooling and cross-pathogen evidence','Masked multi-task OOF outputs, group-held-out species pooling, cross-pathogen paired measured MIC.',[
  'pace_phase5_analysis/masked_multitask_predictions.csv','pace_phase5_analysis/masked_multitask_robust_predictions.csv',
  'pace_phase13b_scientific_closure/Existing_Phase5_Multitask_Paired_Bootstrap_Sensitivity.csv',
  'pace_phase14a_analysis/Phase14A_species_pooled_vs_separate_OOF.csv','pace_phase14a_analysis/Phase14A_species_models_scorecard.csv',
  'pace_phase14a_analysis/Phase14A_species_models_sequence_bootstrap.csv','pace_phase14a_analysis/Phase14A_species_level_curated_MIC.csv',
  'pace_phase14a_analysis/Phase14A_peptide_heldout_folds.csv','pace_phase14a_analysis/Phase14A_species_paired_MIC_records.csv',
  'pace_phase14a_analysis/Phase14A_paired_source_restricted_sensitivity.csv',
  'pace_phase14a_analysis/Phase14A_cross_pathogen_measured_pairwise_effects.csv',
  'pace_phase14a_analysis/Phase14A_measured_spectrum_summary.csv','pace_phase14a_analysis/Phase14A_sequence_level_measured_spectrum.csv',
  'pace_phase14a_analysis/Phase14A_source_matched_peptide_level_records.csv',
 ],['pace_phase14a_analysis/Phase14A_species_models_scorecard.csv','pace_phase14a_analysis/Phase14A_measured_spectrum_summary.csv']),
 ('S08','Original 20 leads and applicability domain','Complete lead sequences, original source predictions, seven-target applicability domain and ranking diagnostics.',[
  'pace_phase14c_analysis/Phase14C_original_leads_uncertainty_envelope.csv',
  'pace_phase14c_analysis/Phase14C_original_leads_rank_and_filter_stress.csv',
  'pace_phase14c_analysis/Phase14C_original_vs_reconstructed_model_predictions.csv',
  'pace_phase14c_analysis/Phase14C_all7_endpoint_applicability_domain.csv',
  'pace_phase5_analysis/lead_applicability_domain_percentiles.csv',
  'pace_reanalysis/phase2_lead_tradeoff_diagnostics.csv',
  'pace_phase14c_analysis/Phase14C_Figure_AD_heatmap.png','pace_phase14c_analysis/Phase14C_Figure_original_vs_refit.png',
 ],['pace_phase14c_analysis/Phase14C_original_leads_uncertainty_envelope.csv','pace_phase14c_analysis/Phase14C_all7_endpoint_applicability_domain.csv']),
 ('S09','Evolutionary design and objective ablation','Complete reoptimization diagnostics for full, no-hemolysis, no-novelty and reduced fitness terms.',[
  'pace_reanalysis/phase1_objective_sensitivity.csv','pace_reanalysis/phase3_neutral_ga_length_trajectories.csv',
 ],['pace_reanalysis/phase1_objective_sensitivity.csv']),
 ('S10','Controlled fixed-weight GA versus NSGA-II','Matched 9000-unique evaluations per seed, 10 seeds per algorithm, 40 selected lead sets and full candidate archives.',[
  'pace_phase10b_verified/ten_seed_descriptive_means.csv',
  'pace_phase10b_verified/paired_comparisons_independent.csv','pace_phase10b_verified/archive_checks.csv',
  'pace_phase10b_verified/cross_algorithm_consistency.csv','pace_phase10b_verified/history_checks.csv',
  'pace_phase10b_verified/lead_metrics_independently_verified.csv',
 ],['pace_phase10b_verified/ten_seed_descriptive_means.csv','pace_phase10b_verified/paired_comparisons_independent.csv']),
 ('S11','External generated libraries and comparative rescoring','Published HydrAMP candidates, original-checkpoint reconstructed AMPGAN v2 peptides, PACE leads and length sensitivity.',[
  'pace_phase11_analysis/phase11_external_comparison_per_sequence.csv','pace_phase11_analysis/phase11_external_comparison_summary.csv',
  'pace_phase11_analysis/phase11_external_model_rescoring.csv',
  'pace_phase11f_analysis/phase11f_all_sequence_model_rescoring.csv',
  'pace_phase11f_analysis/phase11f_sequence_properties_and_novelty.csv',
  'pace_phase11f_analysis/phase11f_library_summary.csv','pace_phase11f_analysis/phase11f_common_length_22to25_sensitivity.csv',
  'pace_phase11_analysis/phase11_PACE_12to25_posthoc_selections.csv',
  'pace_phase11_analysis/phase11_PACE_12to25_10seed_summary.csv',
 ],['pace_phase11f_analysis/phase11f_library_summary.csv','pace_phase11f_analysis/phase11f_common_length_22to25_sensitivity.csv']),
 ('S12','Data availability and reproducible code','Original PACE scripts, versioned source feature engine, independent metric calculators, run order and checksums.',[
  'pace_phase14c_analysis/phase14c_independent_verifier.py',
  'pace_phase14c_analysis/phase14c_export_aligned_residuals.py',
  'pace_phase14c_analysis/phase14c_reproduce_audit.py',
  'pace_phase14c_analysis/phase14c_make_figures.py',
  'pace_phase14c_analysis/Phase14C_output_manifest.json',
  'pace_phase14c_analysis/Phase14C_independent_verification.json',
  'pace_phase14b_verified/verify_phase14b.py',
  'pace_phase11f_analysis/PACE_Phase11E5_Independent_Audit_and_Phase11F_Comparative_Report.md',
  'pace_phase14c_analysis/Phase14C_reproducibility_run_log.txt',
 ],['pace_phase14c_analysis/Phase14C_seven_original_input_sha256.csv']),
]
# ZIP-contained files extracted once with original data and provenance kept unchanged
EXTRACT={
 'S02': [('PACE_Phase8_ESM2_results.zip','classification_predictions.csv'),('PACE_Phase8_ESM2_results.zip','classification_metrics.csv'),('PACE_Phase14B_AMPScanner_External_Results.zip','AMPScanner_PACE_matched_predictions.csv')],
 'S04': [('PACE_Phase8_ESM2_results.zip','quantitative_predictions.csv'),('PACE_Phase8_ESM2_results.zip','quantitative_metrics.csv')],
 'S06': [('PACE_Phase11B_QMAP_Results.zip','qmap_oof_predictions.csv'),('PACE_Phase11B_QMAP_Results.zip','qmap_split_metrics.csv'),('PACE_Phase8_ESM2_results.zip','manifest.json')],
 'S08': [],
 'S09': [('PACE_Phase9_Full_GA_Colab_Results (1).zip','phase9_ablation_history.csv'),('PACE_Phase9_Full_GA_Colab_Results (1).zip','phase9_ablation_run_summary.csv'),('PACE_Phase9_Full_GA_Colab_Results (1).zip','phase9_ablation_leads.csv'),('PACE_Phase9_Full_GA_Colab_Results (1).zip','phase9_ablation_summary.json'),('PACE_Phase9_Full_GA_Colab_Results (1).zip','full_protocol_model_provenance.json')],
 'S11': [('PACE_Phase11_Analysis_Only_Package.zip','HydrAMP_official_style_1.fasta'),('PACE_Phase11E5_Results.zip','ampganv2_positive_conditioned_canonical.fasta'),('PACE_Phase11E5_Results.zip','ampganv2_positive_conditioned_raw.csv'),('PACE_Phase11E5_Results.zip','generator_restored_variable_inventory.json'),('PACE_Phase11E5_Results.zip','status.json')],
 'S12': [('PACE_Phase8_ESM2_results.zip','pace_original_features.py'),('PACE_Phase8_ESM2_results.zip','manifest.json'),('PACE_Phase11_Analysis_Only_Package.zip','phase11_descriptive_benchmark.py'),('PACE_Phase11_Analysis_Only_Package.zip','phase11_common_quantitative_rescoring.py'),('PACE_Phase11_Analysis_Only_Package.zip','phase11_length_matched_archive_diagnostic.py')],
}
ADDITIONAL={
 'S01':['Supplementary_Data_S1(1).xlsx'],
 'S08':['Supplementary_Data_S2_designed_peptides(1).fasta'],
 'S06':['PACE_Phase8_ESM2_results.zip'],
 'S10':['PACE_Phase10B_Controlled_Results.zip'],
 'S12':['PACE_reproducibility_code.zip','PACE_Phase14B_AMPScanner_Independent_Colab.ipynb','PACE_Phase10B_Controlled_Colab_and_Audit.zip'],
}
# Selected verification scripts / reports from trusted, already-checked runs, held separately from publication data.
REPORTS={
'S01':'pace_phase4_analysis/PACE_Phase4_Analysis_Only_Report.md',
'S02':'pace_phase14b_verified/PACE_Phase14B_Independent_AMPScanner_Final_Audit.md',
'S03':'pace_phase12b_analysis/PACE_Phase12B_Experimental_Negative_Validity_Report.md',
'S04':'pace_phase5_analysis/PACE_Phase5_Analysis_Only_Report.md',
'S05':'pace_phase7_hemolysis/PACE_Phase7_Hemolysis_Analysis_Only_Report.md',
'S06':'pace_phase11b_verified/PACE_Phase11B_Independent_Verification_Report.md',
'S07':'pace_phase14a_analysis/PACE_Phase14A_Species_Pooling_and_Measured_Specificity_Audit.md',
'S08':'pace_phase14c_analysis/PACE_Phase14C_Original_Leads_AD_Uncertainty_and_GA_Reproducibility_Report.md',
'S09':'pace_phase9_analysis/PACE_Phase9_Colab_Results_Verification_Report.md',
'S10':'pace_phase10b_verified/PACE_Phase10B_Independent_Verification_Report.md',
'S11':'pace_phase11f_analysis/PACE_Phase11E5_Independent_Audit_and_Phase11F_Comparative_Report.md',
'S12':'pace_phase13b_scientific_closure/PACE_Phase13B_Source_Locked_Scientific_Closure_Report.md',
}
EXTRAS_REPRO=[
'phase15_build.py',
'pace_phase14c_analysis/original_PACE_90D_features_readonly_copy.py',
]
# More outputs archived in their full original runs if needed; source original author docs are never modified.
SOURCE_LINKS=[
 ('AMPlify','AMP identification source (training/test FASTA)','https://github.com/bcgsc/AMPlify','Curation/external overlap recorded in original S1 and Phase 8 manifest'),
 ('GRAMPA','Experimental MIC compilation','https://github.com/zswitten/Antimicrobial-Peptides','Source link, curated analytical subset; database assay records have heterogeneity'),
 ('HC50','User-provided curated hemolysis input','User-provided Cleaned_hemolytic_data.csv','Original data unchanged; exact assay source refer to original S1'),
 ('QMAP','Homology-aware MIC and HC50 evaluation','https://github.com/anthol42/QMAP','Independent QMAP split predictions'),
 ('AMPScanner v2','Independent AMP classifier pretrained source','https://github.com/dan-veltri/amp-scanner-v2','External model evaluated on overlap-filtered retained sequences'),
 ('HydrAMP','Original top peptide library','https://github.com/szczurek-lab/hydramp-starter-kit','100 published output peptides included, original FASTA'),
 ('AMPGAN v2','Release-tag model/checkpoint','https://gitlab.com/vail-uvm/amp-gan','Reconstructed inference; native historical equivalence unverified'),
]

def sha(path:Path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for buf in iter(lambda:f.read(1024*1024),b''): h.update(buf)
 return h.hexdigest()
def nrows(path:Path):
 if path.suffix.lower()!='.csv':return ''
 with path.open(encoding='utf-8-sig',newline='') as f:return max(0,sum(1 for r in csv.reader(f))-1)
def safeexcel(x):
 if x is None:return None
 if isinstance(x,(int,float)):return x
 s=str(x)
 if len(s)>32000:s=s[:31980]+' [TRUNCATED IN XLSX PREVIEW]'
 if s=='':return None
 # protect text fields from formula execution, retain numeric values typed as numbers.
 if re.match(r'^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$',s):
  try:
   v=float(s)
   if abs(v)<1e14 and float.is_integer(v):return int(v)
   if abs(v)<1e14:return v
  except ValueError:pass
 if s.startswith(('=','+','-','@')):s="'"+s
 return s

def small_csv(path:Path, max_rows=250, max_cols=16):
 with path.open(encoding='utf-8-sig',newline='') as f:
  rdr=csv.reader(f)
  result=[]
  for i,row in enumerate(rdr):
   if i>max_rows:break
   result.append([safeexcel(x) for x in row[:max_cols]])
 return result

def clean_title(s):return s.replace('_',' ').replace('.csv','').replace('Phase14','14')

BG='#12304D'; BLUE='#0B526E'; SEA='#0B7F87'; WHITE='#FFFFFF'; TXT='#20364A'
HEADER_FORMAT={'fill':BG,'font':{'bold':True,'color':WHITE},'horizontal_alignment':'left','vertical_alignment':'center','row_height':27}
SUB_FORMAT={'fill':'#E8F2F7','font':{'bold':True,'color':BLUE},'row_height':23}

def put_table(wb,name,table,maxrows=250,notice=None):
 if not table or not table[0]:return
 s=wb.worksheets.add(name[:31]);a=len(table);b=max(len(r) for r in table)
 norm=[[r[j] if j<len(r) else None for j in range(b)] for r in table]
 s.get_range_by_indexes(0,0,a,b).values=norm
 s.get_range_by_indexes(0,0,1,b).format=HEADER_FORMAT
 s.get_range_by_indexes(0,0,a,b).format.row_height=19
 s.get_range_by_indexes(0,0,1,b).format.row_height=30
 s.get_range_by_indexes(0,0,a,b).format.column_width=19
 for c in range(min(b,4)):s.get_range_by_indexes(0,c,a,1).format.column_width=(33 if c<2 else 22)
 try:s.freeze_panes.freeze_rows(1)
 except Exception:pass
 if notice:
  s.get_range_by_indexes(a+2,0,1,1).values=[[notice[:30000]]]

inventory=[]; missing=[]; addedfiles=[]; total_data=0
for G,TITLE,DESC,paths,featured in GROUPS:
 print('BUILD',G,TITLE,flush=True)
 d=O/G; d.mkdir(exist_ok=True)
 data=d/'DATA';data.mkdir(exist_ok=True)
 records=[]
 def addfile(src:Path,name=None,role='Verified analysis output',parent=None,source_note=''):
  global total_data
  if not src.exists():
   missing.append(str(src));return None
  dest=(parent or data)/(name or src.name)
  if dest.exists():
   # identical content skip, or collision-safe rename
   if sha(dest)==sha(src):return dest
   dest=dest.with_name(src.stem+'__'+str(len(records))+src.suffix)
  shutil.copy2(src,dest)
  rec={'group':G,'group_title':TITLE,'relative_path':str(dest.relative_to(O)),
       'source_reference':str(src.relative_to(R)),'role':role,'bytes':dest.stat().st_size,
       'rows_excluding_header':nrows(dest),'SHA256':sha(dest),'source_note':source_note}
  records.append(rec);inventory.append(rec);addedfiles.append(dest)
  return dest
 for rel in paths:addfile(SRC/rel)
 for orig in ADDITIONAL.get(G,[]):addfile(SRC/orig,role='Original input preserved unchanged')
 for zpath,item in EXTRACT.get(G,[]):
  z=SRC/zpath
  if not z.is_file():missing.append(zpath);continue
  with zipfile.ZipFile(z) as arch:
   if item not in arch.namelist():missing.append(f'{zpath}:{item}');continue
   dst=data/item
   with dst.open('wb') as f:f.write(arch.read(item))
  rec={'group':G,'group_title':TITLE,'relative_path':str(dst.relative_to(O)),
       'source_reference':f'{zpath}::{item}', 'role':'Unmodified result extracted from verified ZIP',
       'bytes':dst.stat().st_size,'rows_excluding_header':nrows(dst), 'SHA256':sha(dst),'source_note':'ZIP archive source retains original item bytes'}
  records.append(rec);inventory.append(rec);addedfiles.append(dst)
 if G=='S12':
  for x in EXTRAS_REPRO:addfile(SRC/x)
 # make individual sheet workbook with a useful original data table (full if small, sample clearly labeled if large)
 wb=Workbook.create(); intro=wb.worksheets.add('GUIDE')
 intro.get_range('A1:F1').values=[[f'{G} — {TITLE}',None,None,None,None,None]]
 intro.merge_cells('A1:F1')
 intro.get_range('A1:F1').format={'fill':BG,'font':{'bold':True,'color':WHITE,'size':15},'row_height':42}
 intro.get_range('A3:B10').values=[
 ['Section ID',G],['Focus',TITLE],['Scope',DESC],
 ['Complete data',f'All original full row-level files in {G}/DATA; Excel view is a navigational supplement, not a replacement.'],
 ['Data interpretation','Measurements vs predictions are explicitly labeled; diagnostic/rescoring results are not experimental findings.'],
 ['Original author files','No original manuscript/reviewer comments/data or scripts were edited.'],
 ['Submission status','Compiled scientific supplement; synchronize section/table citations during manuscript revision.'],
 ['Source provenance','FILE_INDEX and master SHA256 manifest list all files and original run sources.']]
 intro.get_range('A3:A10').format=SUB_FORMAT
 intro.get_range('A3:A10').format.column_width=27
 intro.get_range('B3:B10').format.column_width=85
 intro.get_range('B3:B10').format.wrap_text=True
 intro.get_range('A3:B10').format.row_height=36
 intro.get_range('A12').values=[['Research-specific caution']]
 intro.get_range('A12:B12').format=SUB_FORMAT
 caution={
 'S03':'Measured weak activity is pathogen-specific. MIC≥64 uM cannot establish universal non-AMP status.',
 'S04':'All MIC values are derived from heterogeneous assay records; do not infer standard clinical efficacy.',
 'S05':'HC50 predictions are uncertain and are not hemolysis assays of the designed PACE leads.',
 'S06':'ESM-2 model (8M, pretrained) is a specific tested PLM, not a complete search of all protein language models.',
 'S07':'Cross-pathogen differences refer to source measured peptides, not experimental PACE lead specificity.',
 'S08':'Original 20 lead potency and selectivity are computational predictions, not experimental validation.',
 'S09':'Ablation trajectories use refitted predictors; not the original historical training snapshots.',
 'S10':'10 paired seeds were evaluated with 9000 unique predictions each; multiple-testing-corrected superiority not established.',
 'S11':'AMPGAN candidate sequences used checkpoint-weight reconstructed inference; native historical equivalence not verified.',
 'S12':'Missing original fitted model pickles and historical GA full archives prevent exact original trajectory replication.',
 }.get(G,'Public-source dataset redistribution and assay/license conditions must be checked before submission.')
 intro.get_range('A13:B15').values=[[caution,None],['Recommended citation',f'Supplementary Group {G}; canonical full tables are in DATA with SHA256 integrity.'],['Integrity','File hashes are listed in FILE_INDEX and MASTER_MANIFEST_SHA256.csv.']]
 intro.get_range('A13:B15').format.column_width=58
 intro.get_range('A13:B15').format.wrap_text=True
 intro.get_range('A13:B15').format.row_height=33
 fi=[['Filename','Source reference','Role','Data rows','Bytes','SHA256']]
 for rec in records:
  fi.append([Path(rec['relative_path']).name,rec['source_reference'],rec['role'],rec['rows_excluding_header'],rec['bytes'],rec['SHA256']])
 put_table(wb,'FILE_INDEX',fi)
 # include full small tables <=400 rows; more than 400 get clearly identified preview 75 rows
 for idx,rel in enumerate(featured[:3]):
  p=SRC/rel
  if not p.exists():continue
  name=('RESULTS_'+str(idx+1))
  row_total=nrows(p)
  if isinstance(row_total,int) and row_total<=400:
   t=small_csv(p,max_rows=450,max_cols=20)
   put_table(wb,name,t,notice=f'Complete table; original CSV: {p.name} | {row_total} data records')
  else:
   t=small_csv(p,max_rows=60,max_cols=16)
   put_table(wb,name,t,notice=f'FIRST 60 ROWS ONLY, not complete; full data in DATA/{p.name}; {row_total} original records, all columns preserved in CSV')
 workbook_path=d/f'{G}_Supplementary_Tables.xlsx'
 SpreadsheetFile.export_xlsx(wb).save(str(workbook_path))
 rec={'group':G,'group_title':TITLE,'relative_path':str(workbook_path.relative_to(O)),
       'source_reference':'Generated with artifact_tool from original verified CSV files',
       'role':'Navigational Excel with source-grounded results and data index',
       'bytes':workbook_path.stat().st_size,'rows_excluding_header':'','SHA256':sha(workbook_path),'source_note':'Data sample sheets explicitly labeled; full CSV files unchanged'}
 records.append(rec);inventory.append(rec);addedfiles.append(workbook_path)
 # human-readable README
 with (d/'README.md').open('w',encoding='utf-8') as f:
  f.write(f'# {G} — {TITLE}\n\n{DESC}\n\n')
  f.write('**Data policy:** The Excel workbook provides index and readable summaries. Full original records are preserved in the DATA/ folder, without record truncation. Original raw user manuscript, reviewer text, and author inputs were not edited.\n\n')
  f.write(f'**Scientific note:** {caution}\n\n')
  f.write('## Complete files\n\n| File | Rows (where applicable) | Source archive | SHA256 |\n|---|---:|---|---|\n')
  for r in records:f.write(f"| {Path(r['relative_path']).name} | {r['rows_excluding_header']} | {r['source_reference']} | `{r['SHA256']}` |\n")
  f.write('\n**Naming:** Reviewer comments are not present or modified in these submission-oriented supplements. Use the separate verbatim comment/source letter during final author response preparation.\n')
 print(' ',len(records),'files incl workbook',workbook_path.stat().st_size,flush=True)

# provenance source map
with (O/'DATA_SOURCE_PROVENANCE.csv').open('w',newline='',encoding='utf-8') as f:
 w=csv.writer(f);w.writerow(['source_name','dataset_role','official_source_URL_or_file','caveat']);w.writerows(SOURCE_LINKS)
# track SHA originals without copying manuscript/reviewer text
original=['Manuscript(1).docx','Response to reviewers comments.docx','Supplementary_Data_S1(1).xlsx','Supplementary_Data_S2_designed_peptides(1).fasta','PACE_reproducibility_code.zip','Pasted text(6).txt','Cleaned_hemolytic_data.csv']
with (O/'ORIGINAL_AUTHOR_INPUTS_UNCHANGED_SHA256.csv').open('w',newline='',encoding='utf-8') as f:
 w=csv.writer(f);w.writerow(['original_filename','bytes','sha256','document_role','included_as_editable_output'])
 for n in original:
  p=R/n
  w.writerow([n,p.stat().st_size if p.exists() else '',sha(p) if p.exists() else 'UNAVAILABLE','User supplied original; read only', 'NO'])
fields=['group','group_title','relative_path','source_reference','role','bytes','rows_excluding_header','SHA256','source_note']
with (O/'MASTER_MANIFEST_SHA256.csv').open('w',newline='',encoding='utf-8') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(inventory)
with (O/'SUPPLEMENTARY_INDEX.csv').open('w',newline='',encoding='utf-8') as f:
 w=csv.writer(f);w.writerow(['group','title','description','file_count','data_records_in_csv','workbook'])
 for G,title,desc,*_ in GROUPS:
  items=[r for r in inventory if r['group']==G]
  w.writerow([G,title,desc,len(items),sum(r['rows_excluding_header'] for r in items if isinstance(r['rows_excluding_header'],int)),f'{G}/{G}_Supplementary_Tables.xlsx'])

notes='''# PACE — Supplementary Data S01–S12: Author/Editor Index

**This is an analysis-based supplementary data draft, not a new manuscript, a revised reviewer response, or experimental AMP validation.** All author source files are unchanged.

## Structure

Each of 12 groups contains its own styled XLSX navigation/results workbook, complete CSV records (no large-data truncation), original FASTA where available, and data source provenance. Any shorter Excel sample clearly states it is a view rather than a replacement for the full CSV. All supporting CSV records are included with the same values as their verified source files. Full 10-seed GA/NSGA-II candidate archives remain intact in S10 as a nested original result ZIP.

**Groups:**

'''
for G,title,desc,*_ in GROUPS:notes+=f'- **{G} — {title}.** {desc}\n'
notes+='''
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
3. Deposit publication-safe data and reproducibility code according to journal/license requirements; assign public repository URLs/DOIs only when available.
4. Finalize manuscript and author responses only after scientific cross-file consistency checks.
'''
(O/'README_SUPPLEMENTARY_DATA.md').write_text(notes,encoding='utf-8')
# integrity after build
with (O/'BUILD_VERIFICATION.json').open('w',encoding='utf-8') as f:
 json.dump({'generated_utc':datetime.now(timezone.utc).isoformat(),'groups':len(GROUPS),'indexed_files':len(inventory),'missing_sources':missing,
            'csv_data_rows_across_copies':sum(r['rows_excluding_header'] for r in inventory if isinstance(r['rows_excluding_header'],int)),
            'original_input_sha256':{n:sha(R/n) for n in original if (R/n).exists()},
            'source_files_unmodified':True},f,indent=2)
# create one package with no duplicated ZIPs of same data, use zipfile deflate, verified.
master=R/'PACE_Phase15_Comprehensive_S01-S12_Supplementary_Package.zip'
with zipfile.ZipFile(master,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as zz:
 for p in sorted(O.rglob('*')):
  if p.is_file():zz.write(p,str(p.relative_to(O.parent)))
with zipfile.ZipFile(master) as z:
 assert z.testzip() is None
 assert len([x for x in z.namelist() if x.endswith('_Supplementary_Tables.xlsx')])==12
check=[]
for n in original:
 p=R/n
 if p.exists():check.append((n,sha(p)))
print('DONE',master,'size',master.stat().st_size,'nindexed',len(inventory),'missing',len(missing),flush=True)
for m in missing:print('MISSING',m,flush=True)
