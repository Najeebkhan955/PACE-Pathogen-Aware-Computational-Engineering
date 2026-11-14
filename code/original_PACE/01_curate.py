"""Dataset curation for the ESKAPE-AMP study.

Sources (all public, version-controlled):
  AMPlify      : github.com/bcgsc/AMPlify  -- data/*.fa            (Li et al. 2022)
  amPEPpy      : github.com/tlawrence3/amPEPpy -- training_data/   (Lawrence et al. 2021)
  GRAMPA       : github.com/zswitten/Antimicrobial-Peptides -- data/grampa.csv
                 (Witten & Witten 2019; aggregates APD, DBAASP, DRAMP, YADAMP, DADP)
  Hemolytik    : same repo, Hemolysis/Cleaned_hemolytic_data.csv
"""
import pandas as pd, numpy as np, json, os, sys
sys.path.insert(0, os.path.dirname(__file__))
from features import clean

D = '/home/claude/data/'
OUT = '/home/claude/work/'
os.makedirs(OUT, exist_ok=True)
MINLEN, MAXLEN = 5, 60
log = {}


def readfa(p):
    seqs, name, buf = [], None, []
    for line in open(p):
        line = line.strip()
        if line.startswith('>'):
            if name is not None:
                seqs.append(''.join(buf))
            name, buf = line[1:], []
        else:
            buf.append(line)
    if name is not None:
        seqs.append(''.join(buf))
    return seqs


def prep(seqs):
    out = []
    for s in seqs:
        c = clean(s)
        if c and MINLEN <= len(c) <= MAXLEN:
            out.append(c)
    return out


# ---------------------------------------------------------------- classification
A = D + 'AMPlify-master/data/'
raw = {k: readfa(A + v) for k, v in {
    'amp_tr': 'AMPlify_AMP_train_common.fa', 'amp_te': 'AMPlify_AMP_test_common.fa',
    'non_tr': 'AMPlify_non_AMP_train_balanced.fa', 'non_te': 'AMPlify_non_AMP_test_balanced.fa'}.items()}
log['amplify_raw'] = {k: len(v) for k, v in raw.items()}
cl = {k: prep(v) for k, v in raw.items()}
log['amplify_after_filter'] = {k: len(v) for k, v in cl.items()}

# de-duplicate within class, then drop any train sequence also present in test
tr = pd.DataFrame({'sequence': cl['amp_tr'] + cl['non_tr'],
                   'label': [1] * len(cl['amp_tr']) + [0] * len(cl['non_tr'])}).drop_duplicates('sequence')
te = pd.DataFrame({'sequence': cl['amp_te'] + cl['non_te'],
                   'label': [1] * len(cl['amp_te']) + [0] * len(cl['non_te'])}).drop_duplicates('sequence')
overlap = set(tr.sequence) & set(te.sequence)
te = te[~te.sequence.isin(overlap)]
log['train_test_exact_overlap_removed'] = len(overlap)
tr.to_csv(OUT + 'cls_train.csv', index=False)
te.to_csv(OUT + 'cls_test.csv', index=False)
log['cls_train'] = dict(n=len(tr), pos=int(tr.label.sum()))
log['cls_test'] = dict(n=len(te), pos=int(te.label.sum()))

# ---------------------------------------------------------------- external set
P = D + 'amPEPpy-master/training_data/'
ext = pd.DataFrame({
    'sequence': prep(readfa(P + 'M_model_train_AMP_sequence.numbered.fasta')) +
                prep(readfa(P + 'M_model_train_nonAMP_sequence.numbered.proplen.subsample.fasta')),
    'label': None})
n_amp_ext = len(prep(readfa(P + 'M_model_train_AMP_sequence.numbered.fasta')))
ext['label'] = [1] * n_amp_ext + [0] * (len(ext) - n_amp_ext)
ext = ext.drop_duplicates('sequence')
seen = set(tr.sequence) | set(te.sequence)
before = len(ext)
ext = ext[~ext.sequence.isin(seen)]
log['external'] = dict(before_leakage_filter=before, n=len(ext), pos=int(ext.label.sum()))
ext.to_csv(OUT + 'cls_external.csv', index=False)

# ---------------------------------------------------------------- ESKAPE MIC
g = pd.read_csv(D + 'Antimicrobial-Peptides-master/data/grampa.csv')
log['grampa_raw'] = len(g)
g = g[(~g.has_unusual_modification.astype(bool))]
g['seq'] = g.sequence.map(clean)
g = g.dropna(subset=['seq'])
g = g[g.seq.str.len().between(MINLEN, MAXLEN)]
log['grampa_clean'] = len(g)

ESKAPE = {
    'E. faecium':    ['E. faecium', 'E. faecalis'],
    'S. aureus':     ['S. aureus'],
    'K. pneumoniae': ['K. pneumoniae'],
    'A. baumannii':  ['A. baumannii'],
    'P. aeruginosa': ['P. aeruginosa'],
    'E. cloacae':    ['E. cloacae', 'E. aerogenes'],
}
mic = {}
rows = []
for target, names in ESKAPE.items():
    sub = g[g.bacterium.isin(names)]
    agg = sub.groupby('seq').value.agg(['median', 'count', 'std']).reset_index()
    agg.columns = ['sequence', 'log10_MIC_uM', 'n_records', 'sd_records']
    agg['target'] = target
    mic[target] = agg
    rows.append(agg)
    log.setdefault('mic', {})[target] = dict(records=len(sub), unique_peptides=len(agg))
micall = pd.concat(rows, ignore_index=True)
micall.to_csv(OUT + 'mic_eskape.csv', index=False)

# multi-target matrix (peptides assayed against >=2 ESKAPE species)
piv = micall.pivot_table(index='sequence', columns='target', values='log10_MIC_uM')
piv.to_csv(OUT + 'mic_matrix.csv')
log['mic_matrix'] = dict(peptides=len(piv),
                         with_ge2_targets=int((piv.notna().sum(1) >= 2).sum()),
                         with_all6=int((piv.notna().sum(1) == 6).sum()))

# ---------------------------------------------------------------- hemolysis
h = pd.read_csv(D + 'Antimicrobial-Peptides-master/Hemolysis/Cleaned_hemolytic_data.csv')
h['sequence'] = h.Sequence.map(clean)
h = h.dropna(subset=['sequence'])
h = h[h.sequence.str.len().between(MINLEN, MAXLEN)]
h = h.groupby('sequence').log10_HC50.median().reset_index()
h.to_csv(OUT + 'hemolysis.csv', index=False)
log['hemolysis'] = dict(n=len(h), median_log10_HC50=round(float(h.log10_HC50.median()), 3))

json.dump(log, open(OUT + 'curation_log.json', 'w'), indent=2)
print(json.dumps(log, indent=2))
