#!/usr/bin/env python3
"""Checks key Phase 14C claims independently from exported CSVs."""
from pathlib import Path
import csv,json,numpy as np,hashlib,zipfile
D=Path('/mnt/data/pace_phase14c_analysis');ROOT=D.parent

def get(s):
 with (D/s).open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f))
L=get('Phase14C_original_leads_uncertainty_envelope.csv');R=get('Phase14C_original_leads_rank_and_filter_stress.csv');A=get('Phase14C_all7_endpoint_applicability_domain.csv');C=get('Phase14C_original_vs_reconstructed_model_predictions.csv')
assert len(L)==len(R)==20;assert len(A)==140;assert len(C)==140
assert set(x['peptide_id'] for x in L)==set(x['peptide_id'] for x in R)
T=['E._faecium','S._aureus','K._pneumoniae','A._baumannii','P._aeruginosa','E._cloacae']
maxmean=maxworst=0
for row in L:
 mic=np.array([float(row['original_'+t+'_MIC_uM']) for t in T]);pmean=mic.prod()**(1/6)
 maxmean=max(maxmean,abs(pmean-float(row['original_mean_MIC_uM'])))
 maxworst=max(maxworst,abs(mic.max()-float(row['original_worst_MIC_uM'])))
assert maxmean<1e-10 and maxworst<1e-10,(maxmean,maxworst)
def km(s):return {s[i:i+3] for i in range(len(s)-2)}
S=[x['sequence'] for x in L]
J=[len(km(x)&km(y))/len(km(x)|km(y)) for i,x in enumerate(S) for y in S[i+1:]]
assert max(J)<=0.45
assert not any(x['above_reference_p95']=='True' for x in A)
assert not any(x['worst_MIC_stress_below_32']=='True' or x['SI_stress_above_10']=='True' for x in L)
assert all(0<=float(x['winner_frequency_independent'])<=1 for x in R)
assert abs(sum(float(x['winner_frequency_independent']) for x in R)-1)<1e-10
assert abs(sum(float(x['winner_frequency_shared_systematic']) for x in R)-1)<1e-10
assert min(float(x['original_pAMP']) for x in L)>=0.9
assert max(float(x['original_worst_MIC_uM']) for x in L)<=32
assert min(float(x['original_SI']) for x in L)>=10
f={x['file']:x for x in get('Phase14C_seven_original_input_sha256.csv')}
for name,d in f.items():
 assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==d['sha256_before']==d['sha256_after']
report={'verified':True,'original_leads':20,'endpoints':7,'prediction_differences':140,'max_mean_MIC_uM_recalculation_error':maxmean,'max_worst_MIC_uM_recalculation_error':maxworst,'max_pairwise_3mer_jaccard':max(J),'leads_beyond_AD_p95':sum(x['above_reference_p95']=='True' for x in A),'original_source_files_sha256_matched':len(f),'note':'Stress metrics are not prospective probabilities.'}
(D/'Phase14C_independent_verification.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
