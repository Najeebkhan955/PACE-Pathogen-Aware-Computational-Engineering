"""Analysis-only post-hoc PACE length<=25 archive selection, not a new length-restricted GA."""
import zipfile,csv,io,json
from pathlib import Path
from itertools import combinations
from collections import defaultdict
import numpy as np
OUT=Path('/mnt/data/pace_phase11_analysis')
z=zipfile.ZipFile('/mnt/data/PACE_Phase10B_Controlled_Results.zip')

def grams(s):return frozenset(s[i:i+3] for i in range(len(s)-2))
def sim(a,b):
 a,b=grams(a),grams(b);return len(a&b)/len(a|b)

def fitness(x):
 a=2*x['p_amp']-x['mic_mean']
 return a-.5*x['mic_worst']+.6*np.clip(x['ti'],-2,4)-2*max(0,x['sim']-.55)-.35*max(0,abs(x['charge']-6)-3)

saved=[];stats=[]
for seed in (7,29,41,53,67,79,97,113,131,149):
 rows=list(csv.DictReader(io.StringIO(z.read(f'Fixed-weight_GA_{seed}_full_archive.csv').decode())))
 rec=[]
 for r in rows:
  v={k:float(r[k]) for k in ('p_amp','mic_mean','mic_worst','hc50','ti','sim','charge','length')};v['sequence']=r['sequence']
  if v['length']<=25 and v['p_amp']>=.9 and v['mic_worst']<=np.log10(32) and v['ti']>=1 and v['sim']<.55:rec.append(v)
 rec.sort(key=fitness,reverse=True)
 selected=[]
 for r in rec:
  if all(sim(r['sequence'],v['sequence'])<=.45 for v in selected):selected.append(r)
  if len(selected)==20:break
 if len(selected)!=20:raise RuntimeError(f'seed {seed}: only {len(selected)}')
 for i,r in enumerate(selected,1):saved.append({'seed':seed,'rank':i,**r})
 stats.append({'seed':seed,'eligible_12to25':len(rec),'selected_n':len(selected),
 'mean_length':np.mean([r['length'] for r in selected]),
 'mean_pred_gmMIC_uM':np.mean([10**r['mic_mean'] for r in selected]),
 'mean_pred_worstMIC_uM':np.mean([10**r['mic_worst'] for r in selected]),
 'mean_pred_SI':np.mean([10**r['ti'] for r in selected]),
 'mean_pAMP':np.mean([r['p_amp'] for r in selected]),
 'median_pairwise_jaccard':np.median([sim(a['sequence'],b['sequence']) for a,b in combinations(selected,2)])})
for name,records in [('phase11_PACE_12to25_posthoc_selections.csv',saved),('phase11_PACE_12to25_10seed_summary.csv',stats)]:
 with (OUT/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
print(json.dumps({'10seed_mean':{k:float(np.mean([r[k] for r in stats])) for k in stats[0] if k!='seed'},'example_seed7':stats[0]},indent=2))
