"""Stage 6: final lead selection with intra-library diversity control,
benchmark predictions for reference AMPs, and property comparison."""
import numpy as np, pandas as pd, pickle, os, sys, collections, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(__file__))
from features import featurize_many, FEATURE_NAMES

W = '/home/claude/work/'
TARGETS = ['E. faecium', 'S. aureus', 'K. pneumoniae', 'A. baumannii', 'P. aeruginosa', 'E. cloacae']
K = 3


def kmers(s):
    return {s[i:i + K] for i in range(len(s) - K + 1)}


def jac(a, b):
    return len(a & b) / len(a | b) if (a | b) else 0.0


df = pd.read_csv(W + 'ga_all_candidates.csv')
cand = df[(df.p_amp >= 0.90) & (df.mic_worst <= np.log10(32)) &
          (df.ti >= 1.0) & (df.sim < 0.55)].sort_values('fitness', ascending=False)

# greedy diversity selection: no two leads may share >0.45 3-mer Jaccard
sel, sel_k = [], []
for _, r in cand.iterrows():
    ks = kmers(r.sequence)
    if all(jac(ks, k2) <= 0.45 for k2 in sel_k):
        sel.append(r); sel_k.append(ks)
    if len(sel) == 20:
        break
lead = pd.DataFrame(sel).reset_index(drop=True)
lead.insert(0, 'peptide_id', ['PACE-%02d' % i for i in range(1, len(lead) + 1)])
lead.to_csv(W + 'designed_leads.csv', index=False)
print('diverse leads selected:', len(lead))
print(lead[['peptide_id', 'sequence', 'length', 'charge', 'p_amp', 'mic_mean', 'ti', 'sim']].round(3).to_string(index=False))

fun = pd.read_csv(W + 'design_funnel.csv')
fun.loc[fun.step == 'Prioritised leads', 'n'] = len(lead)
fun = pd.concat([fun.iloc[:-1], pd.DataFrame([
    dict(step='Diversity filter (pairwise 3-mer Jaccard <= 0.45)', n=len(lead))])], ignore_index=True)
fun.to_csv(W + 'design_funnel.csv', index=False)

# ------------------------------------------------------------ reference peptides
REF = {
    'Magainin-2':   'GIGKFLHSAKKFGKAFVGEIMNS',
    'LL-37':        'LLGDFFRKSKEKIGKEFKRIVQRIKDFLRNLVPRTES',
    'Melittin':     'GIGAVLKVLTTGLPALISWIKRKRQQ',
    'Pexiganan':    'GIGKFLKKAKKFGKAFVKILKK',
    'Indolicidin':  'ILPWKWPWWPWRR',
    'Buforin-II':   'TRSSRAGLQFPVGRVHRLLRK',
    'Cecropin-B':   'KWKVFKKIEKMGRNIRNGIVKAGPAIAVLGEAKAL',
    'Aurein-1.2':   'GLFDIIKKIAESF',
    'Insulin B-chain (non-AMP control)': 'FVNQHLCGSHLVEALYLVCGERGFFYTPKA',
    'Ubiquitin 1-30 (non-AMP control)':  'MQIFVKTLTGKTITLEVEPSDTIENVKAKI',
}
cls = pickle.load(open(W + 'model_cls_final.pkl', 'rb'))['members']
REG = {t: pickle.load(open(W + 'reg_' + t.replace('. ', '_').replace(' ', '_') + '.pkl', 'rb'))['model'] for t in TARGETS}
HEM = pickle.load(open(W + 'reg_Haemolysis_HC50.pkl', 'rb'))['model']

names, seqs = list(REF), list(REF.values())
X = featurize_many(seqs)
p = np.mean([m.predict_proba(X)[:, 1] for m in cls.values()], axis=0)
mic = np.column_stack([REG[t].predict(X) for t in TARGETS])
ref = pd.DataFrame({'peptide': names, 'sequence': seqs, 'P_AMP': p,
                    'pred_HC50_log10uM': HEM.predict(X)})
for i, t in enumerate(TARGETS):
    ref['pred_' + t] = mic[:, i]

# observed MIC where the peptide exists in GRAMPA (these peptides are in training
# data for the regressors, so this is a recall check, not an independent test)
obs = pd.read_csv(W + 'mic_eskape.csv')
om = obs.pivot_table(index='sequence', columns='target', values='log10_MIC_uM')
for t in TARGETS:
    ref['obs_' + t] = [om[t].get(s, np.nan) if t in om.columns else np.nan for s in seqs]
ref.to_csv(W + 'reference_peptides.csv', index=False)
print('\nReference-peptide check (P(AMP) and predicted S. aureus / P. aeruginosa log10 MIC):')
print(ref[['peptide', 'P_AMP', 'pred_S. aureus', 'obs_S. aureus', 'pred_P. aeruginosa', 'obs_P. aeruginosa']].round(3).to_string(index=False))

# ------------------------------------------------------------ property comparison
tr = pd.read_csv(W + 'cls_train.csv')
nat = tr[tr.label == 1].sequence.tolist()
non = tr[tr.label == 0].sequence.tolist()
Xn = featurize_many(nat); Xo = featurize_many(non); Xl = featurize_many(lead.sequence.tolist())
keys = ['length', 'net_charge_pH7', 'hydrophobic_moment', 'gravy_kd', 'pI',
        'frac_basic', 'frac_aromatic', 'aliphatic_index', 'boman_index', 'charge_density']
rows = []
for k in keys:
    i = FEATURE_NAMES.index(k)
    rows.append(dict(property=k,
                     natural_AMP_mean=Xn[:, i].mean(), natural_AMP_sd=Xn[:, i].std(),
                     non_AMP_mean=Xo[:, i].mean(), non_AMP_sd=Xo[:, i].std(),
                     designed_mean=Xl[:, i].mean(), designed_sd=Xl[:, i].std()))
comp = pd.DataFrame(rows)
comp.to_csv(W + 'property_comparison.csv', index=False)
np.save(W + 'X_natural.npy', Xn); np.save(W + 'X_nonamp.npy', Xo); np.save(W + 'X_leads.npy', Xl)
print('\nProperty comparison:')
print(comp.round(2).to_string(index=False))

with open(W + 'designed_leads.fasta', 'w') as fh:
    for _, r in lead.iterrows():
        fh.write('>%s | P_AMP=%.3f | mean_log10MIC=%.3f | SI=%.1f | novelty=%.3f\n%s\n'
                 % (r.peptide_id, r.p_amp, r.mic_mean, 10 ** r.ti, r.sim, r.sequence))
print('\nFASTA written')
