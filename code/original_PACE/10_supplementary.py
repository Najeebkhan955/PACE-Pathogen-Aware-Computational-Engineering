import numpy as np, pandas as pd, json, sys, os
sys.path.insert(0, os.path.dirname(__file__))
from features import FEATURE_NAMES, GLOBAL_NAMES, AA, CTD_ATTR
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

W = '/home/claude/work/'
OUT = '/mnt/user-data/outputs/'
os.makedirs(OUT, exist_ok=True)
TARGETS = ['E. faecium', 'S. aureus', 'K. pneumoniae', 'A. baumannii', 'P. aeruginosa', 'E. cloacae']

HDR_FILL = PatternFill('solid', fgColor='1A3350')
HDR_FONT = Font(name='Arial', bold=True, color='FFFFFF', size=10)
BODY = Font(name='Arial', size=10)
TITLE = Font(name='Arial', bold=True, size=12, color='1A3350')
NOTE = Font(name='Arial', italic=True, size=9, color='555555')
THIN = Side(style='thin', color='D0D5DC')
BORD = Border(left=THIN, right=THIN, top=THIN, bottom=THIN)

wb = openpyxl.Workbook()
wb.remove(wb.active)


def sheet(name, title, note, df, widths=None, numfmt=None):
    ws = wb.create_sheet(name[:31])
    ws['A1'] = title; ws['A1'].font = TITLE
    ws['A2'] = note; ws['A2'].font = NOTE
    r0 = 4
    for j, c in enumerate(df.columns, 1):
        cell = ws.cell(r0, j, str(c))
        cell.fill = HDR_FILL; cell.font = HDR_FONT
        cell.alignment = Alignment(horizontal='center', vertical='center', wrap_text=True)
        cell.border = BORD
    for i, (_, row) in enumerate(df.iterrows(), r0 + 1):
        for j, v in enumerate(row, 1):
            if isinstance(v, (np.floating, float)) and not isinstance(v, bool):
                v = None if pd.isna(v) else float(v)
            elif isinstance(v, (np.integer,)):
                v = int(v)
            cell = ws.cell(i, j, v)
            cell.font = BODY; cell.border = BORD
            if isinstance(v, float):
                cell.number_format = (numfmt or {}).get(df.columns[j - 1], '0.000')
    for j, c in enumerate(df.columns, 1):
        w = (widths or {}).get(c)
        if w is None:
            w = min(34, max(11, int(df[c].astype(str).str.len().max() if len(df) else 10) + 3, len(str(c)) + 3))
        ws.column_dimensions[get_column_letter(j)].width = w
    ws.freeze_panes = ws.cell(r0 + 1, 1)
    ws.row_dimensions[r0].height = 30
    return ws


# ---------------------------------------------------------------- S1 sources
log = json.load(open(W + 'curation_log.json'))
src = pd.DataFrame([
    dict(Dataset='AMPlify benchmark', Role='Classifier training / held-out test',
         Repository='github.com/bcgsc/AMPlify', File='data/AMPlify_{AMP,non_AMP}_{train,test}_*.fa',
         Raw_records=sum(log['amplify_raw'].values()), After_curation=log['cls_train']['n'] + log['cls_test']['n'],
         Reference='Li et al., BMC Genomics 2022'),
    dict(Dataset='amPEPpy M-model set', Role='External validation',
         Repository='github.com/tlawrence3/amPEPpy', File='training_data/M_model_*.fasta',
         Raw_records=log['external']['before_leakage_filter'], After_curation=log['external']['n'],
         Reference='Lawrence et al., Bioinformatics 2021'),
    dict(Dataset='GRAMPA', Role='Pathogen-resolved MIC regression',
         Repository='github.com/zswitten/Antimicrobial-Peptides', File='data/grampa.csv',
         Raw_records=log['grampa_raw'], After_curation=log['grampa_clean'],
         Reference='Witten & Witten, bioRxiv 2019 (aggregates APD, DBAASP, DRAMP, YADAMP, DADP)'),
    dict(Dataset='Hemolytik-derived HC50', Role='Haemolysis regression',
         Repository='github.com/zswitten/Antimicrobial-Peptides', File='Hemolysis/Cleaned_hemolytic_data.csv',
         Raw_records=1377, After_curation=log['hemolysis']['n'],
         Reference='Gautam et al., Database 2014 (Hemolytik)'),
])
sheet('S1_Data_sources', 'Table S1. Data provenance and curation',
      'Curation: uppercase, 20 standard residues only, length 5-60, exact de-duplication, removal of '
      'sequences with unusual chemical modifications, and removal of any external sequence also present in training.',
      src, widths={'File': 42, 'Reference': 52, 'Repository': 40, 'Role': 32, 'Dataset': 24})

curation = pd.DataFrame([
    dict(Step='Classifier training set', n=log['cls_train']['n'], Positives=log['cls_train']['pos']),
    dict(Step='Classifier held-out test set', n=log['cls_test']['n'], Positives=log['cls_test']['pos']),
    dict(Step='External validation set', n=log['external']['n'], Positives=log['external']['pos']),
    dict(Step='Train/test exact overlaps removed', n=log['train_test_exact_overlap_removed'], Positives=None),
    dict(Step='External sequences removed for overlap with training',
         n=log['external']['before_leakage_filter'] - log['external']['n'], Positives=None),
    dict(Step='Peptides with MIC on >=2 ESKAPE targets', n=log['mic_matrix']['with_ge2_targets'], Positives=None),
    dict(Step='Peptides with MIC on all six ESKAPE targets', n=log['mic_matrix']['with_all6'], Positives=None),
])
sheet('S2_Curation_log', 'Table S2. Curation audit trail',
      'Every count is produced by scripts/01_curate.py and is reproducible from the raw archives.', curation)

# ---------------------------------------------------------------- S3 feature dictionary
desc = {}
for n in GLOBAL_NAMES:
    desc[n] = 'Global physicochemical descriptor'
desc.update({
    'length': 'Number of residues', 'net_charge_pH7': 'Net charge at pH 7.4 (Henderson-Hasselbalch, EMBOSS pKa)',
    'charge_density': 'Net charge divided by length', 'pI': 'Isoelectric point (Bio.SeqUtils.IsoelectricPoint)',
    'gravy_kd': 'Mean Kyte-Doolittle hydropathy', 'hydrophobicity_eis': 'Mean Eisenberg consensus hydrophobicity',
    'hydrophobic_moment': 'Maximum Eisenberg mean hydrophobic moment, 100 deg, 11-residue window',
    'amphipathicity_ratio': 'Hydrophobic moment divided by |mean hydrophobicity|',
    'aliphatic_index': 'Ikai aliphatic index', 'boman_index': 'Mean Boman protein-binding potential (kcal/mol)',
    'instability_like': 'Mean absolute Boman difference between adjacent residues (relative index)',
    'flexibility': 'Mean Vihinen flexibility index',
    'max_hydrophobic_run': 'Longest consecutive run of hydrophobic residues',
    'max_cationic_run': 'Longest consecutive run of K/R', 'kr_ratio': '(count K + 1) / (count R + 1)',
    'hydromoment_x_charge': 'Hydrophobic moment multiplied by net charge',
})
rows = []
for f in FEATURE_NAMES:
    if f.startswith('AAC_'):
        d, fam = f'Fractional composition of {f[4:]}', 'Amino-acid composition'
    elif '_C' in f or '_T' in f:
        attr = f.rsplit('_', 1)[0]; code = f.rsplit('_', 1)[1]
        d = (f'CTD composition, group {code[1]} of {attr}' if code.startswith('C')
             else f'CTD transition frequency between groups {code[1]} and {code[2]} of {attr}')
        fam = 'CTD'
    else:
        d, fam = desc.get(f, 'Global physicochemical descriptor'), 'Physicochemical'
    rows.append(dict(Feature=f, Family=fam, Description=d))
sheet('S3_Feature_dictionary', 'Table S3. The 90-descriptor feature set',
      'CTD attribute groupings follow Dubchak et al. (1995). All scales are hard-coded in scripts/features.py.',
      pd.DataFrame(rows), widths={'Description': 62, 'Feature': 24, 'Family': 24})

# ---------------------------------------------------------------- S4-S5 classifier
cm = pd.read_csv(W + 'cls_metrics.csv')
sheet('S4_Classifier_metrics', 'Table S4. AMP recognition performance',
      'Threshold 0.5 for all threshold-dependent metrics. CV = 5-fold stratified cross-validation on the training set.',
      cm, widths={'Model': 24, 'Set': 22})
sheet('S5_CV_folds', 'Table S5. Per-fold cross-validation results',
      'Fold-level metrics underlying the stability analysis in Figure 3d.',
      pd.read_csv(W + 'cls_folds.csv'), widths={'Model': 24})

# ---------------------------------------------------------------- S6 regression
rm = pd.read_csv(W + 'reg_metrics.csv')
sheet('S6_Regression_metrics', 'Table S6. Potency and haemolysis regression performance',
      'All values are 5-fold cross-validated on log10-transformed MIC (uM) or HC50 (uM).',
      rm, widths={'Target': 22, 'Model': 20})

cc = pd.read_csv(W + 'mic_cross_corr.csv', index_col=0).reset_index().rename(columns={'index': 'Target'})
sheet('S7_Cross_pathogen', 'Table S7. Cross-pathogen Spearman correlation of measured log10 MIC',
      'Computed on peptides with measurements against both members of each pair (minimum 25 shared peptides).',
      cc, widths={'Target': 20})
npc = pd.read_csv(W + 'mic_pair_counts.csv', index_col=0).reset_index().rename(columns={'index': 'Target'})
sheet('S8_Pair_counts', 'Table S8. Shared-peptide counts per pathogen pair',
      'Denominators for Table S7.', npc, widths={'Target': 20},
      numfmt={c: '0' for c in npc.columns})

# ---------------------------------------------------------------- S9 importance
pi = pd.read_csv(W + 'perm_importance.csv')
sh = pd.read_csv(W + 'shap_importance.csv')
imp = pi.merge(sh, on='feature', how='outer').sort_values('importance', ascending=False)
imp.columns = ['Feature', 'Permutation_importance_AUC_drop', 'SD', 'Mean_abs_SHAP']
sheet('S9_Feature_importance', 'Table S9. Feature importance',
      'Permutation importance from the gradient-boosting model on the held-out test set (10 repeats); '
      'SHAP values from the extra-trees model on 400 randomly selected test peptides.',
      imp, widths={'Feature': 28}, numfmt={'Permutation_importance_AUC_drop': '0.00000',
                                           'SD': '0.00000', 'Mean_abs_SHAP': '0.00000'})

sheet('S10_Similarity_bins', 'Table S10. Performance stratified by sequence similarity',
      'Test peptides binned by their maximum 3-mer Jaccard similarity to any training peptide.',
      pd.read_csv(W + 'similarity_bins.csv'))

# ---------------------------------------------------------------- S11 reference peptides
ref = pd.read_csv(W + 'reference_peptides.csv')
ws = sheet('S11_Reference_peptides', 'Table S11. Model behaviour on characterised reference peptides',
           'Predicted values for the MIC models are recall checks, not independent predictions: these peptides are '
           'present in the GRAMPA training data for the regressors. The P(AMP) column is the classifier ensemble.',
           ref, widths={'peptide': 34, 'sequence': 52})

# ---------------------------------------------------------------- S12 leads
lead = pd.read_csv(W + 'designed_leads.csv')
cols = ['peptide_id', 'sequence', 'length', 'charge', 'p_amp', 'mic_mean', 'mic_worst',
        'hc50', 'ti', 'sim'] + ['mic_' + t for t in TARGETS]
lead = lead[cols]
ws = sheet('S12_Designed_leads', 'Table S12. Prioritised de novo designed peptides',
           'All MIC and HC50 values are model predictions on the log10 (uM) scale and have NOT been experimentally '
           'validated. Columns U-V are computed in-sheet from the predictions.',
           lead, widths={'sequence': 36, 'peptide_id': 13})
r0 = 5
ws.cell(4, len(cols) + 1, 'Mean MIC (uM)').fill = HDR_FILL
ws.cell(4, len(cols) + 1).font = HDR_FONT
ws.cell(4, len(cols) + 2, 'Selectivity index').fill = HDR_FILL
ws.cell(4, len(cols) + 2).font = HDR_FONT
mic_col = get_column_letter(cols.index('mic_mean') + 1)
ti_col = get_column_letter(cols.index('ti') + 1)
for i in range(len(lead)):
    r = r0 + i
    c1 = ws.cell(r, len(cols) + 1, f'=POWER(10,{mic_col}{r})'); c1.font = BODY; c1.number_format = '0.00'; c1.border = BORD
    c2 = ws.cell(r, len(cols) + 2, f'=POWER(10,{ti_col}{r})'); c2.font = BODY; c2.number_format = '0.0'; c2.border = BORD
ws.column_dimensions[get_column_letter(len(cols) + 1)].width = 15
ws.column_dimensions[get_column_letter(len(cols) + 2)].width = 16

# ---------------------------------------------------------------- S13-S15
sheet('S13_Design_funnel', 'Table S13. Candidate prioritisation funnel',
      'Filters applied in the order shown.', pd.read_csv(W + 'design_funnel.csv'),
      widths={'step': 56}, numfmt={'n': '#,##0'})
sheet('S14_GA_history', 'Table S14. Genetic-algorithm optimisation history',
      'Population statistics per generation (population 240, elite 24, 45 generations).',
      pd.read_csv(W + 'ga_history.csv'))
sheet('S15_Property_comparison', 'Table S15. Physicochemical comparison',
      'Means and standard deviations across natural AMPs, non-AMPs and the 20 designed leads.',
      pd.read_csv(W + 'property_comparison.csv'), widths={'property': 26})

wb.save(OUT + 'Supplementary_Data_S1.xlsx')
print('workbook saved with', len(wb.sheetnames), 'sheets:', ', '.join(wb.sheetnames))
