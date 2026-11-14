"""Independent verification of user-returned PACE Phase14B external pretrained AMPScanner run.
Analysis-only. Original user inputs are opened read-only.
"""
import csv,io,json,zipfile,hashlib,collections, pathlib
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score,balanced_accuracy_score,matthews_corrcoef,brier_score_loss,confusion_matrix

ROOT=pathlib.Path('/mnt/data')
OUT=ROOT/'pace_phase14b_verified';OUT.mkdir(exist_ok=True)
Z=ROOT/'PACE_Phase14B_AMPScanner_External_Results.zip'
MODS=['PACE90_ET','ESM2_ET','Hybrid_ET','Composition21_ET','AMPScanner_original2018']
def readz(z,p): return list(csv.DictReader(io.StringIO(z.read(p).decode('utf-8'))))
def savecsv(name,rows,cols=None):
 p=OUT/name
 if not rows: return
 with p.open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=cols or list(rows[0]));w.writeheader();w.writerows(rows)
 return str(p)
def metr(y,p):
 y=np.asarray(y,dtype=int);p=np.asarray(p,dtype=float);pred=(p>=0.5).astype(int)
 tn,fp,fn,tp=confusion_matrix(y,pred,labels=[0,1]).ravel()
 return dict(n=len(y),n_amp=int(sum(y)),n_nonamp=int(len(y)-sum(y)),ROC_AUC=roc_auc_score(y,p),PR_AUC=average_precision_score(y,p),Balanced_Accuracy=balanced_accuracy_score(y,pred),MCC=matthews_corrcoef(y,pred),Brier=brier_score_loss(y,p),TP=int(tp),TN=int(tn),FP=int(fp),FN=int(fn),positive_prediction_rate=float(np.mean(pred)))

with zipfile.ZipFile(Z) as z:
 assert z.testzip() is None
 rows=readz(z,'AMPScanner_PACE_matched_predictions.csv'); removed=readz(z,'source_overlap_exclusion_log.csv'); rep=readz(z,'AMPScanner_PACE_verified_metric_comparison.csv'); bs=readz(z,'AMPScanner_PACE_paired_bootstrap.csv'); meta=json.loads(z.read('provenance_and_readiness.json'));run=json.loads(z.read('run_status.json'))

assert run['stage']=='complete' and meta['stage']=='complete'
assert len(rows)==3050 and len(removed)==366 and meta['n_predictions']==3050
assert len({(r['split'],r['sequence']) for r in rows})==len(rows)
assert set(r['split'] for r in rows)=={'heldout','external'}
assert set(r['label'] for r in rows)=={'0','1'}
assert all(int(r['length'])==len(r['sequence']) and int(r['length'])>=10 for r in rows)
assert all(r['exclusion_reason']=='' and r['known_AMPScanner_trainvaltest_overlap']=='0' for r in rows)
assert all(r['AMPScanner_original2018']!='' for r in rows)
assert all(0<=float(r[m])<=1 and np.isfinite(float(r[m])) for r in rows for m in MODS)
assert all(r['AMPScanner_original2018'] is not None for r in rows)
assert not(set((r['split'],r['sequence']) for r in rows) & set((r['split'],r['sequence']) for r in removed))
assert len(set(r['sequence'] for r in rows if r['split']=='external') & set(r['sequence'] for r in rows if r['split']=='heldout'))==0
assert len(rows)+len(removed)==3416

# Verify exact exported scores and labels against original Phase8 payload.
phase8_map={}
with zipfile.ZipFile(ROOT/'PACE_Phase8_ESM2_results.zip') as z:
 for r in readz(z,'classification_predictions.csv'):
  if r['model']=='ExtraTrees200':
   k=(r['split'],r['sequence'],r['label'])
   phase8_map.setdefault(k,{})[r['representation']]=float(r['score'])
maprepr={'PACE90_ET':'90D_PACE','ESM2_ET':'ESM2_320D','Hybrid_ET':'Hybrid_90D_ESM2','Composition21_ET':'21D_length_AAC'}
# Validate exact actual representation names from source
names=collections.Counter(k for v in phase8_map.values() for k in v)
print('Phase8 rep names',names)
for key in maprepr:
 if maprepr[key] not in names:
  matches=[s for s in names if key.split('_')[0].lower() in s.lower()];print('No exact name',key,'candidates',matches)

match_count=0;error_max=0.;notfound=[]
for r in rows:
 k=(r['split'],r['sequence'],r['label'])
 if k not in phase8_map: notfound.append(k);continue
 match_count+=1
 for m,src in maprepr.items():
  if src in phase8_map[k]:
   error_max=max(error_max,abs(float(r[m])-phase8_map[k][src]))
print('Phase8 matched rows',match_count,'missing',len(notfound),'maxscore absdiff',error_max)
assert not notfound

# Verify previously published 30 metric records exact.
metrics=[];errs=[]
for split in ['heldout','external']:
 for subset in ['all_eligible','length_10to30','length_31to60']:
  sub=[r for r in rows if r['split']==split and (subset=='all_eligible' or (10<=int(r['length'])<=30 if subset=='length_10to30' else 31<=int(r['length'])<=60))]
  y=np.array([int(r['label']) for r in sub]);assert len(set(y))==2
  for model in MODS:
   p=np.array([float(r[model]) for r in sub]); d={'split':split,'subset':subset,'model':model,**metr(y,p)};metrics.append(d)
   matches=[r for r in rep if (r['split'],r['subset'],r['model'])==(split,subset,model)]
   assert len(matches)==1
   for field in ['ROC_AUC','PR_AUC','Balanced_Accuracy','MCC','Brier']:
    e=abs(d[field]-float(matches[0][field]));errs.append((e,split,subset,model,field))
   assert d['n']==int(matches[0]['n'])
assert len(metrics)==30
print('Metric worst deviation',max(errs))

# Independently recompute stratified paired AUROC intervals using vectorized average ranks.
# This is the same Mann-Whitney AUROC statistic as roc_auc_score.
from scipy.stats import rankdata
N=600;rng=np.random.default_rng(1414);out=[]
for split in ['heldout','external']:
 for subset in ['all_eligible','length_10to30','length_31to60']:
  sub=[r for r in rows if r['split']==split and (subset=='all_eligible' or (10<=int(r['length'])<=30 if subset=='length_10to30' else 31<=int(r['length'])<=60))]
  y=np.array([int(r['label']) for r in sub]);s=np.array([float(r['AMPScanner_original2018']) for r in sub]);i0=np.where(y==0)[0];i1=np.where(y==1)[0];n0=len(i0);n1=len(i1)
  boot_indices=np.concatenate([rng.choice(i0,(N,n0),replace=True),rng.choice(i1,(N,n1),replace=True)],axis=1)
  def boot_auc(scores):
   ranks=rankdata(scores[boot_indices],axis=1,method='average')
   return (ranks[:,n0:].sum(axis=1)-n1*(n1+1)/2)/(n0*n1)
  auc_base=roc_auc_score(y,s); auc_base_boot=boot_auc(s)
  for model in MODS[:-1]:
   p=np.array([float(r[model]) for r in sub]);delta=roc_auc_score(y,p)-auc_base
   arr=boot_auc(p)-auc_base_boot
   lo,hi=np.quantile(arr,[.025,.975]);
   out.append(dict(split=split,subset=subset,model=model,comparator='AMPScanner_original2018',n=len(y),AUC_delta=delta,bootstrap_CI95_low=lo,bootstrap_CI95_high=hi,bootstrap_reps=N,two_sided_sign_p=2*min((np.sum(arr<=0)+1)/(N+1),(np.sum(arr>=0)+1)/(N+1)),significance_at_95='yes' if lo>0 or hi<0 else 'no'))

# compare length-conditional difference of joint splits; proxy for heterogeneity.
exclusions=collections.Counter((r['split'],r['exclusion_reason']) for r in removed)
counts=[]
for split in ['heldout','external']:
 counts.append(dict(split=split,n_original=1568 if split=='heldout' else 1848,n_evaluated=sum(r['split']==split for r in rows),n_removed=sum(r['split']==split for r in removed),n_short=exclusions[(split,'unsupported_length_lt10')],n_dev_overlap=exclusions[(split,'known_AMPScanner_2018_development_overlap')]))

# Check full original source file integrity snapshot, do not modify originals.
snapshot=list(csv.DictReader((ROOT/'pace_phase14b_analysis/phase14b_original_source_sha256_snapshot.csv').open()))
hashes=[]
for r in snapshot:
 p=ROOT/r['file'];actual=hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None
 hashes.append(dict(filename=r['file'],expected_sha256=r['sha256'],actual_sha256=actual,match='yes' if actual==r['sha256'] else 'NO'))
assert all(r['match']=='yes' for r in hashes)

savecsv('independently_recomputed_classification_metrics.csv',metrics)
savecsv('independent_stratified_paired_bootstrap.csv',out)
savecsv('evaluated_and_excluded_counts.csv',counts)
savecsv('original_source_integrity_checked.csv',hashes)
sumry={'zip_readable':True,'stage':run['stage'],'n_matched_peptides':len(rows),'n_excluded':len(removed),'per_split':counts,'max_abs_recalculated_metric_discrepancy':max(x[0] for x in errs),'phase8_reconciled_rows':match_count,'phase8_max_abs_prediction_difference':error_max,'checkpoint_sha256_reported':meta['source_urls']['OriginalPaper_081917_FULL_MODEL.h5']['sha256'],'original_inputs_integrity':all(x['match']=='yes' for x in hashes),'statistical_scope':'Exploratory stratified sequence bootstrap only; no homology-aware train/test overlap resolution beyond exact matching; different model vintages.','n_models':len(MODS),'n_stratification_rows':len(metrics)}
(OUT/'verification_manifest.json').write_text(json.dumps(sumry,indent=2))
print('SUMMARY',json.dumps(sumry,indent=2));print('BOOTSTRAP')
for x in out:
 if x['model'] in ['PACE90_ET','Hybrid_ET']:print(x)
