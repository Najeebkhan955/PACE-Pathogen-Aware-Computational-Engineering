"""Exports exact phase5/7 random/group OOF residuals as a self-contained audit table."""
from pathlib import Path
import csv
D=Path('/mnt/data/pace_phase14c_analysis');R=D.parent
T=['E. faecium','S. aureus','K. pneumoniae','A. baumannii','P. aeruginosa','E. cloacae']
K=['E._faecium','S._aureus','K._pneumoniae','A._baumannii','P._aeruginosa','E._cloacae']
def records(path):
 with open(path,encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
all_random=records(R/'pace_phase5_analysis/mic_all6_oof_error_domain.csv');out=[]
for t,k in zip(T,K):
 group=records(R/f'pace_phase5_analysis/cluster_cv_{k}.csv')
 lookup={x['sequence']:x for x in group}
 rows=[x for x in all_random if x['target']==t]
 assert len(rows)==len(lookup)==len(set(x['sequence'] for x in rows))
 for row in rows:
  v=lookup[row['sequence']]
  assert abs(float(row['observed'])-float(v['observed']))<1e-8
  out.append(dict(target=t,sequence=row['sequence'],observed_log10=float(row['observed']),ordinary_cv_pred_log10=float(row['predicted']),ordinary_residual_observed_minus_pred=float(row['residual']),similarity_group_cv_pred_log10=float(v['prediction_group_cv']),similarity_group_residual_observed_minus_pred=float(v['observed'])-float(v['prediction_group_cv']),similarity_cluster_id=v['cluster_id'],ordinary_cv_fold=row['fold'],source='PACE Phase5 300tree OOF + similarity-group 150tree OOF'))
hr=records(R/'pace_phase7_hemolysis/hemolysis_300tree_extratrees_oof.csv');hg={x['sequence']:x for x in records(R/'pace_phase7_hemolysis/hemolysis_similarity_cluster_oof.csv')};assert len(hr)==len(hg)==1271
for row in hr:
 v=hg[row['sequence']]
 assert abs(float(row['observed_log10_hc50'])-float(v['observed']))<1e-8
 out.append(dict(target='HC50',sequence=row['sequence'],observed_log10=float(row['observed_log10_hc50']),ordinary_cv_pred_log10=float(row['oof_pred_log10_hc50']),ordinary_residual_observed_minus_pred=float(row['residual']),similarity_group_cv_pred_log10=float(v['group_prediction']),similarity_group_residual_observed_minus_pred=float(v['observed'])-float(v['group_prediction']),similarity_cluster_id=v['cluster'],ordinary_cv_fold=row['fold'],source='PACE Phase7 300tree OOF + similarity-group 150tree OOF'))
p=D/'Phase14C_all7_aligned_out_of_fold_residuals.csv'
with p.open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
print('Aligned residual rows:',len(out),'per endpoint:',{t:sum(x['target']==t for x in out) for t in T+['HC50']})
