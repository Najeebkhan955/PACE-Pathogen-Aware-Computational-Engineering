"""Analysis only: compare two published candidate sequence libraries under one descriptor/novelty protocol."""
from pathlib import Path
from collections import defaultdict
from itertools import combinations
import zipfile, io, hashlib, json, csv, statistics
import numpy as np
import importlib.util

ROOT=Path('/mnt/data'); OUT=ROOT/'pace_phase11_analysis'; OUT.mkdir(exist_ok=True)
SOURCES={
    'HydrAMP_top100': OUT/'HydrAMP_official_style_1.fasta',
    'PACE_original20': ROOT/'Supplementary_Data_S2_designed_peptides(1).fasta',
}

def fasta(path):
    ss=[]
    for ln in Path(path).read_text().splitlines():
        ln=ln.strip()
        if ln and not ln.startswith('>'): ss.append(ln)
    assert all(5 <= len(s) <= 60 and set(s) <= set('ACDEFGHIKLMNPQRSTVWY') for s in ss)
    assert len(ss)==len(set(ss)), 'not unique'
    return ss

with zipfile.ZipFile(ROOT/'PACE_Phase8_ESM2_results.zip') as z:
    feat_path=OUT/'features_original_source.py'
    feat_path.write_bytes(z.read('pace_original_features.py'))
    source_seqs=np.load(io.BytesIO(z.read('all_embeddings_and_features.npz')))['sequences'].tolist()

spec=importlib.util.spec_from_file_location('features_original_source',str(feat_path));ft=importlib.util.module_from_spec(spec);spec.loader.exec_module(ft)
print('module functions:',[s for s in dir(ft) if not s.startswith('_') and callable(getattr(ft,s))][-15:])

def ngrams(seq,n=3): return frozenset(seq[i:i+n] for i in range(len(seq)-n+1))

def pairj(sa,sb):
    a,b=ngrams(sa),ngrams(sb)
    return len(a&b)/len(a|b) if a|b else 0

ref_set=set(source_seqs)
ref_tri=[ngrams(s) for s in source_seqs]
idx=defaultdict(list)
for i,gr in enumerate(ref_tri):
    for x in gr:idx[x].append(i)


def maxnov(seq):
    g=ngrams(seq)
    candidates=set()
    for x in g:candidates.update(idx.get(x,()))
    best=0
    for i in candidates:
        t=ref_tri[i]
        val=len(g&t)/(len(g)+len(t)-len(g&t))
        if val>best:best=val
    return best

def features(seqs):
    for name in ('featurize_many','features','extract_features','sequence_features','compute_features','feature_vector'):
        if hasattr(ft,name):
            fn=getattr(ft,name)
            try:
                xx=fn(seqs)
                if np.asarray(xx).shape==(len(seqs),90):return np.asarray(xx,float)
            except Exception:pass
    raise RuntimeError('Feature extraction API not recognized')

records=[];metrics={}
for name,path in SOURCES.items():
    seqs=fasta(path)
    X=features(seqs)
    dists=np.array([pairj(a,b) for a,b in combinations(seqs,2)])
    nn=np.array([maxnov(s) for s in seqs])
    idx_name=ft.FEATURE_NAMES if hasattr(ft,'FEATURE_NAMES') else None
    print('names',idx_name[:12] if idx_name else 'none')
    lens=np.array([len(s) for s in seqs])
    # Original engine: descriptor first columns length, charge, charge density, pI, GRAVY, ...
    vals={
        'N':len(seqs),'length_mean':float(lens.mean()),'length_min':int(lens.min()),'length_max':int(lens.max()),
        'length_sd':float(lens.std(ddof=1)),
        'net_charge_mean':float(X[:,1].mean()),'net_charge_sd':float(X[:,1].std(ddof=1)),
        'hydrophobic_moment_mean':float(X[:,6].mean()),'gravy_mean':float(X[:,4].mean()),
        'fraction_C_present':float(np.mean(['C' in s for s in seqs])),
        'fraction_M_present':float(np.mean(['M' in s for s in seqs])),
        'fraction_exact_match_reference':float(np.mean([s in ref_set for s in seqs])),
        'pairwise_3mer_Jaccard_median':float(np.median(dists)),
        'pairwise_3mer_Jaccard_mean':float(np.mean(dists)),
        'pairwise_3mer_Jaccard_max':float(np.max(dists)),
        'fraction_pairs_above_045':float(np.mean(dists>0.45)),
        'known_corpus_max3mer_Jaccard_mean':float(nn.mean()),
        'known_corpus_max3mer_Jaccard_median':float(np.median(nn)),
        'known_corpus_max3mer_Jaccard_max':float(nn.max()),
        'fraction_sequences_above_055_reference':float(np.mean(nn>=.55)),
        'fraction_cationic_gt8':float(np.mean(X[:,1]>8)),
    }
    metrics[name]=vals
    for i,(s,row,nv) in enumerate(zip(seqs,X,nn),1):
        records.append({'library':name,'candidate_id':i,'sequence':s,'length':len(s),
          'net_charge_pH7_4':row[1],'pI_bounded':row[3],
          'hydrophobic_moment':row[6],'gravy':row[4],
          'max_jaccard3_to_PACE_reference':nv,'exact_PACE_reference_match':int(s in ref_set),
          'contains_C':int('C' in s),'contains_M':int('M' in s)})

with (OUT/'phase11_external_comparison_per_sequence.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
with (OUT/'phase11_external_comparison_summary.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['library']+list(next(iter(metrics.values()))));w.writeheader();w.writerows([{'library':k,**v} for k,v in metrics.items()])
(OUT/'phase11_external_comparison_summary.json').write_text(json.dumps(metrics,indent=2))
print(json.dumps(metrics,indent=2))
