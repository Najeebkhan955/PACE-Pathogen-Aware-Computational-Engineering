"""ANALYSIS ONLY. Full 300-tree PACE 90D quantitative predictors refit on published Phase 8 labels.
Common PACE scoring of external HydrAMP top-100 and original PACE-20. Not an independent
unbiased efficacy test because PACE leads were optimized with related reward models.
"""
import zipfile, io, csv, importlib.util, time, json
from pathlib import Path
from collections import defaultdict
import numpy as np
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.metrics import r2_score
ROOT=Path('/mnt/data');OUT=ROOT/'pace_phase11_analysis'
z=zipfile.ZipFile(ROOT/'PACE_Phase8_ESM2_results.zip')
db=np.load(io.BytesIO(z.read('all_embeddings_and_features.npz')))
seqs=db['sequences'].tolist();XX=db['X90'].astype(np.float32)
lookup={s:i for i,s in enumerate(seqs)}
ft_spec=importlib.util.spec_from_file_location('fengine',OUT/'features_original_source.py')
f=importlib.util.module_from_spec(ft_spec);ft_spec.loader.exec_module(f)
external=[s.strip() for s in (OUT/'HydrAMP_top100_sequences.txt').read_text().splitlines() if s.strip()]
pace=[s.strip() for s in (ROOT/'Supplementary_Data_S2_designed_peptides(1).fasta').read_text().splitlines() if s and not s.startswith('>')]
matched_rows=list(csv.DictReader((OUT/'phase11_PACE_12to25_posthoc_selections.csv').open()))
matched=[r['sequence'] for r in matched_rows]
checks=pace+external+matched
Xcheck=f.featurize_many(checks).astype(np.float32)
assert Xcheck.shape==(320,90)
q=list(csv.DictReader(io.StringIO(z.read('quantitative_predictions.csv').decode())))
metrics=[];preds={}
for endpoint in ['E. faecium','S. aureus','K. pneumoniae','A. baumannii','P. aeruginosa','E. cloacae','HC50']:
  selected=[r for r in q if r['endpoint']==endpoint and r['representation']=='90D_PACE']
  ind=[lookup[r['sequence']] for r in selected]
  labels=np.asarray([float(r['observed_log10']) for r in selected],dtype=float)
  t=time.time()
  model=ExtraTreesRegressor(n_estimators=300,random_state=2024,n_jobs=4,
                            min_samples_split=2,min_samples_leaf=1,max_features=1.0)
  model.fit(XX[ind],labels)
  out=model.predict(Xcheck)
  preds[endpoint]=out.tolist()
  print('Trained',endpoint,'n',len(selected),'seconds',round(time.time()-t,2),'range',round(out.min(),3),round(out.max(),3),flush=True)
  metrics.append({'endpoint':endpoint,'n_train':len(selected),'model':'ExtraTreesRegressor(300, random_state=2024, n_jobs=4)','n_features':90})
rows=[]
for i,s in enumerate(checks):
  mean=np.mean([preds[k][i] for k in ['E. faecium','S. aureus','K. pneumoniae','A. baumannii','P. aeruginosa','E. cloacae']])
  worst=np.max([preds[k][i] for k in ['E. faecium','S. aureus','K. pneumoniae','A. baumannii','P. aeruginosa','E. cloacae']])
  ti=preds['HC50'][i]-mean
  r={'library':('PACE_original20' if i<20 else ('HydrAMP_top100' if i<120 else 'PACE_posthoc_12to25')),'candidate_id':1+i if i<20 else (i-19 if i<120 else i-119),'sequence':s,'length':len(s),
     'pred_mean_log10_MIC':mean,'pred_geometric_mean_MIC_uM':10**mean,'pred_worst_log10_MIC':worst,
     'pred_worst_MIC_uM':10**worst,'pred_log10_HC50':preds['HC50'][i],
     'pred_HC50_uM':10**preds['HC50'][i],'pred_log10_SI':ti,'pred_fold_SI':10**ti}
  for k in preds:r['pred_'+k.replace(' ','_')+'_log10']=preds[k][i]
  rows.append(r)
with (OUT/'phase11_external_model_rescoring.csv').open('w',newline='') as h:
 w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
summary={}
for lib in ['PACE_original20','HydrAMP_top100','PACE_posthoc_12to25']:
 v=[r for r in rows if r['library']==lib]
 summary[lib]={
 'n':len(v),'mean_pred_geomMIC_uM':float(np.mean([r['pred_geometric_mean_MIC_uM'] for r in v])),
 'median_pred_geomMIC_uM':float(np.median([r['pred_geometric_mean_MIC_uM'] for r in v])),
 'mean_pred_worstMIC_uM':float(np.mean([r['pred_worst_MIC_uM'] for r in v])),
 'median_pred_worstMIC_uM':float(np.median([r['pred_worst_MIC_uM'] for r in v])),
 'mean_pred_SI':float(np.mean([r['pred_fold_SI'] for r in v])),
 'median_pred_SI':float(np.median([r['pred_fold_SI'] for r in v])),
 'pass_pred_MIC_32':int(sum(r['pred_worst_MIC_uM']<=32 for r in v)),
 'pass_pred_SI_10':int(sum(r['pred_fold_SI']>=10 for r in v)),
 'pass_both':int(sum(r['pred_worst_MIC_uM']<=32 and r['pred_fold_SI']>=10 for r in v)),
 }
(OUT/'phase11_common_scorer_summary.json').write_text(json.dumps({'training':metrics,'summary':summary},indent=2))
print(json.dumps(summary,indent=2))
