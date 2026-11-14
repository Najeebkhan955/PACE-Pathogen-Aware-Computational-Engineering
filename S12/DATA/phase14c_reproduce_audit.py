#!/usr/bin/env python3
"""Analysis-only reproduction of PACE Phase 14C.
Never rewrites original PACE inputs, original reviewers, or original supplementary materials.
Uses original S12 point predictions, phase5/7 heldout residuals and 90D feature engine.
Important: residual stress tests are NOT calibrated prediction probabilities for new leads.
"""
import os, io, csv, json, hashlib, zipfile, importlib.util
from pathlib import Path
from collections import Counter
import numpy as np
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler
from artifact_tool import Blob, SpreadsheetFile
ROOT=Path('/mnt/data'); OUT=ROOT/'pace_phase14c_analysis'; OUT.mkdir(exist_ok=True)
TARGETS=['E. faecium','S. aureus','K. pneumoniae','A. baumannii','P. aeruginosa','E. cloacae']
KEYS=['E._faecium','S._aureus','K._pneumoniae','A._baumannii','P._aeruginosa','E._cloacae']
N_DRAW=12000
SEED=20261114
SOURCE=['Manuscript(1).docx','Supplementary_Data_S1(1).xlsx','Supplementary_Data_S2_designed_peptides(1).fasta','PACE_reproducibility_code.zip','Pasted text(6).txt','Cleaned_hemolytic_data.csv','Response to reviewers comments.docx']
def hashes():return {s:hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in SOURCE}
before=hashes()
def readcsv(f):
 with open(f,encoding='utf-8-sig',newline='') as z:return list(csv.DictReader(z))
def csvout(name,rows,fields=None):
 p=OUT/name
 if fields is None: fields=list(rows[0]) if rows else []
 with open(p,'w',encoding='utf-8',newline='') as z:
  wr=csv.DictWriter(z,fieldnames=fields,extrasaction='ignore');wr.writeheader();wr.writerows(rows)
 return p
# Read original user workbook through approved artifact_tool (read-only)
w=SpreadsheetFile.import_xlsx(Blob.load(str(ROOT/'Supplementary_Data_S1(1).xlsx')))
v=w.worksheets.get_item('S12_Designed_leads').get_range('A4:R24').values
header=[str(x) for x in v[0]]; leads=[dict(zip(header,r)) for r in v[1:] if r[0] and str(r[0]).startswith('PACE-')]
assert len(leads)==20,len(leads)
seqs=[r['sequence'] for r in leads]
fasta=(ROOT/'Supplementary_Data_S2_designed_peptides(1).fasta').read_text()
fa=[]
for line in fasta.splitlines():
 if line.startswith('>'):fa.append(['', '']) ;fa[-1][0]=line.split()[0][1:]
 elif line.strip():fa[-1][1]+=line.strip()
assert [(r['peptide_id'],r['sequence']) for r in leads]==[(i,s) for i,s in fa]
assert len(set(seqs))==20
log_32=np.log10(32.)
MIC=np.array([[float(x['mic_'+t]) for t in TARGETS] for x in leads]); H=np.array([float(x['hc50']) for x in leads]);
P=np.array([float(x['p_amp']) for x in leads]); CH=np.array([float(x['charge']) for x in leads]); SIM=np.array([float(x['sim']) for x in leads]);
mean=MIC.mean(1);worst=MIC.max(1);ti=H-mean
fitness=lambda mean,worst,ti:(2*P -mean -0.5*worst +.6*np.clip(ti,-2,4)-2*np.maximum(0,SIM-.55)-.35*np.maximum(0,np.abs(CH-6)-3))
fit=fitness(mean,worst,ti)
summary_recon={
 'n_leads':len(leads),
 'max_abs_mean_logMIC_delta':float(np.max(np.abs(mean-np.array([float(x['mic_mean']) for x in leads])))),
 'max_abs_worst_logMIC_delta':float(np.max(np.abs(worst-np.array([float(x['mic_worst']) for x in leads])))),
 'max_abs_selectivity_log10_delta':float(np.max(np.abs(ti-np.array([float(x['ti']) for x in leads])))),
 'all_final_hard_filter_pass':bool(np.all((P>=.9)&(worst<=log_32)&(ti>=1)&(SIM<.55))),
 'all_no_novelty_penalty':bool(np.all(SIM<=.55))}
def km(s):return {s[i:i+3] for i in range(len(s)-2)}
jacs=[]
for i in range(20):
 for j in range(i+1,20):
  a,b=km(seqs[i]),km(seqs[j]);jacs.append(len(a&b)/len(a|b));
summary_recon.update(max_pairwise_jaccard=float(max(jacs)),median_pairwise_jaccard=float(np.median(jacs)),count_length_ge29=int(sum(len(s)>=29 for s in seqs)),original_fitness_sort_matches=bool(np.all(np.diff(fit)<=1e-8)))
# Independent source CV error estimates
source_oof={}
source_group={}
for t,key in zip(TARGETS,KEYS):
 oo=readcsv(ROOT/'pace_phase5_analysis'/'mic_all6_oof_error_domain.csv');rr=[x for x in oo if x['target']==t]
 source_oof[t]=np.array([float(x['residual']) for x in rr]) # observed minus predicted
 gg=readcsv(ROOT/'pace_phase5_analysis'/('cluster_cv_'+key+'.csv'))
 source_group[t]=np.array([float(x['observed'])-float(x['prediction_group_cv']) for x in gg])
source_oof['HC50']=np.array([float(x['residual']) for x in readcsv(ROOT/'pace_phase7_hemolysis'/'hemolysis_300tree_extratrees_oof.csv')])
source_group['HC50']=np.array([float(x['observed'])-float(x['group_prediction']) for x in readcsv(ROOT/'pace_phase7_hemolysis'/'hemolysis_similarity_cluster_oof.csv')])
q90={k:float(np.quantile(np.abs(v),.9,method='higher')) for k,v in source_group.items()}
q90_oof={k:float(np.quantile(np.abs(v),.9,method='higher')) for k,v in source_oof.items()}
calrows=[]
for t in TARGETS+['HC50']:
 for kind,source in [('random_oof',source_oof),('similarity_group_oof',source_group)]:
  e=source[t];q=float(np.quantile(np.abs(e),.9,method='higher'))
  calrows.append(dict(target=t,source=kind,n=len(e),rmse_log10=float(np.sqrt(np.mean(e**2))),mae_log10=float(np.mean(np.abs(e))),abs_error_q90_log10=q,concentration_upper_lower_fold=10**(2*q),comment='Diagnostic empirical residual quantile; not a calibrated coverage guarantee for PACE leads'))
csvout('Phase14C_heldout_residual_quantiles.csv',calrows)
# Position in measured MIC 90D AD (existing) and new HC50 90D feature-space AD
lead_ad=readcsv(ROOT/'pace_phase5_analysis'/'lead_applicability_domain_percentiles.csv')
ad_by={(r['peptide'],r['target']):r for r in lead_ad}
with zipfile.ZipFile(ROOT/'PACE_Phase8_ESM2_results.zip') as z:
 ndata=np.load(io.BytesIO(z.read('all_embeddings_and_features.npz')))
 seq_ref=np.asarray(ndata['sequences']); Xref=np.asarray(ndata['X90'],dtype=np.float64)
 fsource=OUT/'original_PACE_90D_features_readonly_copy.py';fsource.write_bytes(z.read('pace_original_features.py'))
spec=importlib.util.spec_from_file_location('pace14c_original_feats',str(fsource));mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
Xlead=mod.featurize_many(seqs)
refmap={str(s):i for i,s in enumerate(seq_ref)}
hemo_rows=readcsv(ROOT/'pace_phase7_hemolysis'/'hemolysis_original_curation.csv')
hemo_seqs=[r['sequence'] for r in hemo_rows if r['sequence'] in refmap];Xh=Xref[[refmap[s] for s in hemo_seqs]]
assert len(hemo_seqs)==len(hemo_rows)==1271,(len(hemo_seqs),len(hemo_rows))
scale=StandardScaler().fit(Xh);Hx=scale.transform(Xh);Hl=scale.transform(Xlead)
kn=NearestNeighbors(n_neighbors=6,metric='euclidean',n_jobs=1).fit(Hx)
selfd=kn.kneighbors(Hx,n_neighbors=6,return_distance=True)[0][:,1:].mean(axis=1)
leadd=kn.kneighbors(Hl,n_neighbors=5,return_distance=True)[0].mean(axis=1)
hemo_p95=float(np.quantile(selfd,.95)); hemo_pct=[float(100*np.mean(selfd<=d)) for d in leadd]
adrows=[]
for i,l in enumerate(leads):
 for t in TARGETS:
  r=ad_by[(l['peptide_id'],t)]
  adrows.append(dict(peptide_id=l['peptide_id'],target=t,lead_mean5NN_distance=float(r['lead_mean5NN_distance']),reference_percentile=float(r['training_self_excluded_mean5NN_percentile']),reference_p95=float(r['reference_p95']),above_reference_p95=bool(str(r['above_reference_p95'])=='True'),distance_protocol='Existing Phase5 endpoint-specific 90D scaler and 5NN'))
 adrows.append(dict(peptide_id=l['peptide_id'],target='HC50',lead_mean5NN_distance=float(leadd[i]),reference_percentile=hemo_pct[i],reference_p95=hemo_p95,above_reference_p95=bool(leadd[i]>hemo_p95),distance_protocol='New Phase14C hemolysis full-data standardized 90D and 5NN'))
csvout('Phase14C_all7_endpoint_applicability_domain.csv',adrows)
# Score reconstruction versus previous refit, for transparency
refit={r['sequence']:r for r in readcsv(ROOT/'pace_phase11f_analysis'/'phase11f_all_sequence_model_rescoring.csv') if r['library']=='PACE_original20'}
shift=[]
for i,l in enumerate(leads):
 r=refit[l['sequence']]
 for j,t in enumerate(TARGETS):
  shift.append(dict(peptide_id=l['peptide_id'],endpoint=t,original_log10=float(MIC[i,j]),independent_refit_log10=float(r['pred_'+KEYS[j]+'_log10']),refit_minus_original_log10=float(r['pred_'+KEYS[j]+'_log10'])-MIC[i,j]))
 shift.append(dict(peptide_id=l['peptide_id'],endpoint='HC50',original_log10=H[i],independent_refit_log10=float(r['pred_HC50_log10']),refit_minus_original_log10=float(r['pred_HC50_log10'])-H[i]))
csvout('Phase14C_original_vs_reconstructed_model_predictions.csv',shift)
# 90th-percentile absolute-error envelope, NO joint or lead-level coverage claim
qmic=np.array([q90[t] for t in TARGETS]);qhc=q90['HC50']
envrows=[]
for i,l in enumerate(leads):
 worstU=float(np.max(MIC[i]+qmic));worstL=float(np.max(MIC[i]-qmic));
 lowSI=float(H[i]-qhc-np.mean(MIC[i]+qmic))
 hiSI=float(H[i]+qhc-np.mean(MIC[i]-qmic))
 obj=dict(peptide_id=l['peptide_id'],sequence=l['sequence'],length=len(l['sequence']),original_pAMP=P[i],original_charge=float(CH[i]),original_known_3mer_similarity=float(SIM[i]),original_pred_mean_log10MIC=float(mean[i]),original_pred_worst_log10MIC=float(worst[i]),original_pred_HC50_log10=float(H[i]),original_fitness=float(fit[i]),original_mean_MIC_uM=float(10**mean[i]),original_worst_MIC_uM=float(10**worst[i]),original_pred_HC50_uM=float(10**H[i]),original_SI=float(10**ti[i]),worst_MIC_upper_stress_uM=float(10**worstU),worst_MIC_lower_stress_uM=float(10**worstL),SI_lower_stress=float(10**lowSI),SI_upper_stress=float(10**hiSI),worst_MIC_stress_below_32=bool(worstU<=log_32),SI_stress_above_10=bool(lowSI>=1),max_MIC_domain_reference_percentile=max(float(ad_by[(l['peptide_id'],t)]['training_self_excluded_mean5NN_percentile']) for t in TARGETS),HC50_domain_reference_percentile=hemo_pct[i]);
 for j,t in enumerate(TARGETS):
  tag=KEYS[j];obj['original_'+tag+'_MIC_uM']=float(10**MIC[i,j]);obj['upper90_residual_stress_'+tag+'_MIC_uM']=float(10**(MIC[i,j]+qmic[j]))
 envrows.append(obj)
csvout('Phase14C_original_leads_uncertainty_envelope.csv',envrows)
# Two deliberately extreme residual correlation assumptions for lead ranking stability
rng=np.random.default_rng(SEED);stats={}
for mode in ['independent_peptide_endpoint','shared_endpoint_systematic']:
 if mode=='independent_peptide_endpoint':
  mm=np.stack([rng.choice(source_group[t],size=(N_DRAW,20),replace=True) for t in TARGETS],axis=2)+MIC[None,:,:]
  hh=H[None,:]+rng.choice(source_group['HC50'],size=(N_DRAW,20),replace=True)
 else:
  mm=MIC[None,:,:]+np.stack([rng.choice(source_group[t],size=(N_DRAW,1),replace=True) for t in TARGETS],axis=2)
  hh=H[None,:]+rng.choice(source_group['HC50'],size=(N_DRAW,1),replace=True)
 avg=mm.mean(2);wm=mm.max(2);tiv=hh-avg
 fs=(2*P[None,:]-avg-.5*wm+.6*np.clip(tiv,-2,4)-2*np.maximum(0,SIM[None,:]-.55)-.35*np.maximum(0,np.abs(CH[None,:]-6)-3))
 args=np.argsort(-fs,axis=1)
 winner=np.bincount(args[:,0],minlength=20)/N_DRAW
 top5=np.array([(args[:,:5]==i).any(axis=1).mean() for i in range(20)])
 all_filters=((P[None,:]>=.9)&(wm<=log_32)&(tiv>=1)&(SIM[None,:]<.55))
 passfreq=all_filters.mean(0)
 base_order=np.argsort(-fit)
 base_top5=set(base_order[:5]);top5_overlap=np.mean([len(base_top5.intersection(set(a[:5]))) for a in args])
 stats[mode]={'winner_freq':winner,'top5_freq':top5,'hard_filter_stress_freq':passfreq,'mean_top5_overlap':float(top5_overlap),'fraction_all_20_pass':float(all_filters.all(axis=1).mean()),'mean_pass_count':float(all_filters.sum(axis=1).mean()),'reference_top1_stress_win_freq':float(winner[base_order[0]])}
rankrows=[]
for i,l in enumerate(leads):
 rankrows.append(dict(peptide_id=l['peptide_id'],original_rank=int(np.where(np.argsort(-fit)==i)[0][0])+1,original_fitness=float(fit[i]),winner_frequency_independent=float(stats['independent_peptide_endpoint']['winner_freq'][i]),top5_frequency_independent=float(stats['independent_peptide_endpoint']['top5_freq'][i]),hard_filter_pass_frequency_independent=float(stats['independent_peptide_endpoint']['hard_filter_stress_freq'][i]),winner_frequency_shared_systematic=float(stats['shared_endpoint_systematic']['winner_freq'][i]),top5_frequency_shared_systematic=float(stats['shared_endpoint_systematic']['top5_freq'][i]),hard_filter_pass_frequency_shared_systematic=float(stats['shared_endpoint_systematic']['hard_filter_stress_freq'][i]),note='Uncalibrated residual stress-test frequencies; not real-world probabilities'))
csvout('Phase14C_original_leads_rank_and_filter_stress.csv',rankrows)
# Original GA source code checks and provenance
with zipfile.ZipFile(ROOT/'PACE_reproducibility_code.zip') as z:
 ga=z.read('code/05_design_ga.py').decode(); sel=z.read('code/06_leads.py').decode()
 source_check={'stage5_seed7':('np.random.default_rng(7)' in ga),'stage5_population240_generations45_elite24':('POP, GENS, ELITE = 240, 45, 24' in ga),'stage5_original_ranked_top25':('f4.head(25)' in ga),'stage6_original_top20_greedy_0.45':("if len(sel) == 20:" in sel and '0.45' in sel),'stage5_five_opt_terms':all(s in ga for s in ['2.0 * pa[i]','0.5 * mic_worst','0.6 * float(np.clip(ti','max(0.0, sim - 0.55)','abs(q - 6.0)']),'stage5_requires_external_model_pickle':("model_cls_final.pkl" in ga and 'reg_Haemolysis_HC50.pkl' in ga),'original_model_pickle_in_zip':any(zf.endswith('.pkl') for zf in z.namelist()),'original_ga_candidate_archive_in_zip':any(zf.endswith('ga_all_candidates.csv') for zf in z.namelist())}
phase9z=ROOT/'PACE_Phase9_Full_GA_Colab_Results (1).zip';recon_intersection=None
if phase9z.exists():
 with zipfile.ZipFile(phase9z) as z:
  leads9=list(csv.DictReader(io.StringIO(z.read('phase9_ablation_leads.csv').decode())))
  full_set={x['sequence'] for x in leads9 if x.get('arm')=='full' or x.get('configuration')=='full'}
  recon_intersection={'rows':len(leads9),'full_arm_n':len(full_set),'full_arm_overlap_original20':len(set(seqs)&full_set),'all_arms_overlap_original20':len(set(seqs)&{x['sequence'] for x in leads9})}
summary_recon.update(source_check)
summary_recon.update({'phase9_refit_comparison':recon_intersection,'median_hemolysis_domain_percentile':float(np.median(hemo_pct)),'hemolysis_leads_above_p95':int(sum(np.asarray(leadd)>hemo_p95)), 'all7_AD_above_p95':int(sum(r['above_reference_p95'] for r in adrows)),'n_group_error_cv_sources':{k:len(v) for k,v in source_group.items()},'group_q90_log10':q90,'random_q90_log10':q90_oof,'group_uncertainty_stress_worst_MIC_all20_below32':int(sum(r['worst_MIC_stress_below_32'] for r in envrows)),'group_uncertainty_stress_SI_all20_above10':int(sum(r['SI_stress_above_10'] for r in envrows)),'fitness_original_top_id':leads[int(np.argmax(fit))]['peptide_id'],'stress_draws_per_scenario':N_DRAW,'stress_seed':SEED,'ranking_stress_summary':{k:{kk:v for kk,v in val.items() if not isinstance(v,np.ndarray)} for k,val in stats.items()},'source_input_sha256':before})
after=hashes();assert before==after,'SOURCES MUTATED'
summary_recon['source_hashes_unchanged']=True
(OUT/'Phase14C_integrity_and_summary.json').write_text(json.dumps(summary_recon,indent=2,ensure_ascii=False))
csvout('Phase14C_seven_original_input_sha256.csv',[{'file':k,'sha256_before':v,'sha256_after':after[k],'same':True} for k,v in before.items()])
print(json.dumps({k:v for k,v in summary_recon.items() if k in ['n_leads','max_pairwise_jaccard','median_pairwise_jaccard','all7_AD_above_p95','median_hemolysis_domain_percentile','hemolysis_leads_above_p95','group_uncertainty_stress_worst_MIC_all20_below32','group_uncertainty_stress_SI_all20_above10','ranking_stress_summary','phase9_refit_comparison','original_fitness_sort_matches','source_hashes_unchanged']},indent=2))
