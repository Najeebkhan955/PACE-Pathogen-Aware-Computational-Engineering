"""Stage 4: model interpretation and robustness analysis."""
import numpy as np, pandas as pd, json, pickle, os, sys, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(__file__))
from features import FEATURE_NAMES
from sklearn.inspection import permutation_importance
from sklearn.metrics import roc_auc_score, matthews_corrcoef

W = '/home/claude/work/'
rng = np.random.default_rng(2024)

Xtr, ytr = np.load(W + 'Xtr.npy'), np.load(W + 'ytr.npy')
Xte, yte = np.load(W + 'Xte.npy'), np.load(W + 'yte.npy')
tr = pd.read_csv(W + 'cls_train.csv'); te = pd.read_csv(W + 'cls_test.csv')
pte = np.load(W + 'pte_ENS.npy')

# ---------------------------------------------------------- permutation importance
gb = pickle.load(open(W + 'clsmodel_GB.pkl', 'rb'))
pi = permutation_importance(gb, Xte, yte, n_repeats=10, random_state=2024, scoring='roc_auc')
imp = pd.DataFrame({'feature': FEATURE_NAMES, 'importance': pi.importances_mean,
                    'sd': pi.importances_std}).sort_values('importance', ascending=False)
imp.to_csv(W + 'perm_importance.csv', index=False)
print('Top 12 features by permutation importance (AUC drop):')
print(imp.head(12).round(4).to_string(index=False), flush=True)

# ---------------------------------------------------------- SHAP
try:
    import shap
    et = pickle.load(open(W + 'clsmodel_ET.pkl', 'rb'))
    idx = rng.choice(len(Xte), 400, replace=False)
    ex = shap.TreeExplainer(et)
    sv = ex.shap_values(Xte[idx], check_additivity=False)
    sv = sv[..., 1] if np.ndim(sv) == 3 else sv
    np.save(W + 'shap_values.npy', sv)
    np.save(W + 'shap_X.npy', Xte[idx])
    ms = pd.DataFrame({'feature': FEATURE_NAMES, 'mean_abs_shap': np.abs(sv).mean(0)}
                      ).sort_values('mean_abs_shap', ascending=False)
    ms.to_csv(W + 'shap_importance.csv', index=False)
    print('\nTop 10 by mean |SHAP|:', ', '.join(ms.feature.head(10)), flush=True)
except Exception as e:
    print('SHAP failed:', e)

# ---------------------------------------------------------- redundancy-stratified test
def kmers(s, k=3):
    return {s[i:i + k] for i in range(len(s) - k + 1)}

train_sets = [kmers(s) for s in tr.sequence]
train_pos = [kmers(s) for s, l in zip(tr.sequence, tr.label) if l == 1]


def max_jaccard(s, pool):
    a = kmers(s)
    if not a:
        return 0.0
    best = 0.0
    for b in pool:
        if not b:
            continue
        inter = len(a & b)
        if inter == 0:
            continue
        j = inter / len(a | b)
        if j > best:
            best = j
            if best > 0.999:
                break
    return best


sim = np.array([max_jaccard(s, train_sets) for s in te.sequence])
np.save(W + 'test_similarity.npy', sim)
bins = [(0, .2), (.2, .4), (.4, .6), (.6, .8), (.8, 1.01)]
rows = []
for lo, hi in bins:
    m = (sim >= lo) & (sim < hi)
    if m.sum() < 20 or len(np.unique(yte[m])) < 2:
        continue
    rows.append(dict(bin=f'{lo:.1f}-{hi:.1f}', n=int(m.sum()),
                     ROC_AUC=roc_auc_score(yte[m], pte[m]),
                     MCC=matthews_corrcoef(yte[m], (pte[m] >= .5).astype(int))))
sb = pd.DataFrame(rows)
sb.to_csv(W + 'similarity_bins.csv', index=False)
print('\nPerformance vs. max 3-mer similarity to training set:')
print(sb.round(3).to_string(index=False), flush=True)

# ---------------------------------------------------------- cross-pathogen correlation
piv = pd.read_csv(W + 'mic_matrix.csv', index_col=0)
corr = piv.corr(min_periods=25, method='spearman')
corr.to_csv(W + 'mic_cross_corr.csv')
n_pairs = pd.DataFrame(index=piv.columns, columns=piv.columns, dtype=float)
for a in piv.columns:
    for b in piv.columns:
        n_pairs.loc[a, b] = int(piv[[a, b]].dropna().shape[0])
n_pairs.to_csv(W + 'mic_pair_counts.csv')
print('\nCross-pathogen Spearman of measured log10 MIC:')
print(corr.round(2).to_string())
