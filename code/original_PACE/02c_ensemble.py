import numpy as np, pandas as pd, json, pickle, os, sys, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(__file__))
from sklearn.metrics import (roc_auc_score, average_precision_score, matthews_corrcoef,
                             accuracy_score, f1_score, confusion_matrix, roc_curve,
                             precision_recall_curve)
W = '/home/claude/work/'
MEMBERS = ['GB', 'ET', 'SVM', 'RF']
LABEL = 'Soft-voting ensemble'


def metrics(y, p, thr=0.5):
    yh = (p >= thr).astype(int)
    tn, fp, fn, tp = confusion_matrix(y, yh).ravel()
    return dict(Accuracy=accuracy_score(y, yh), Sensitivity=tp / (tp + fn),
                Specificity=tn / (tn + fp), Precision=tp / (tp + fp),
                F1=f1_score(y, yh), MCC=matthews_corrcoef(y, yh),
                ROC_AUC=roc_auc_score(y, p), PR_AUC=average_precision_score(y, p))


ytr, yte, yex = np.load(W + 'ytr.npy'), np.load(W + 'yte.npy'), np.load(W + 'yex.npy')
Xte, Xex = np.load(W + 'Xte.npy'), np.load(W + 'Xex.npy')
oof = np.mean([np.load(W + f'oof_{k}.npy') for k in MEMBERS], axis=0)
pte = np.mean([np.load(W + f'pte_{k}.npy') for k in MEMBERS], axis=0)
mdls = {k: pickle.load(open(W + f'clsmodel_{k}.pkl', 'rb')) for k in MEMBERS}
pex = np.mean([mdls[k].predict_proba(Xex)[:, 1] for k in MEMBERS], axis=0)

rows = [dict(Model=LABEL, Set='5-fold CV (train)', **metrics(ytr, oof)),
        dict(Model=LABEL, Set='Held-out test', **metrics(yte, pte)),
        dict(Model=LABEL, Set='External validation', **metrics(yex, pex))]
pd.DataFrame(rows).to_csv(W + 'cls_metrics.csv', mode='a', header=False, index=False)
np.save(W + 'oof_ENS.npy', oof); np.save(W + 'pte_ENS.npy', pte); np.save(W + 'pex_ENS.npy', pex)

fpr, tpr, _ = roc_curve(yte, pte); pr, rc, _ = precision_recall_curve(yte, pte)
c = json.load(open(W + 'cls_curves.json'))
c[LABEL] = dict(fpr=fpr[::max(1, len(fpr)//400)].tolist(), tpr=tpr[::max(1, len(tpr)//400)].tolist(),
                prec=pr[::max(1, len(pr)//400)].tolist(), rec=rc[::max(1, len(rc)//400)].tolist())
json.dump(c, open(W + 'cls_curves.json', 'w'))

res = pd.read_csv(W + 'cls_metrics.csv')
print(res[res.Set == 'Held-out test'].round(3).to_string(index=False))
print()
print(res[res.Set == 'External validation'].round(3).to_string(index=False))

# ---- final production model: refit ensemble members on train+test
Xtr = np.load(W + 'Xtr.npy')
Xfull, yfull = np.vstack([Xtr, Xte]), np.concatenate([ytr, yte])
import importlib
mod = importlib.import_module('02b_classify_one') if False else None
sys.path.insert(0, '/home/claude/scripts')
from importlib.machinery import SourceFileLoader
b = SourceFileLoader('b', '/home/claude/scripts/02b_classify_one.py').load_module()
final = {}
for k in MEMBERS:
    m = b.build(k); m.fit(Xfull, yfull); final[k] = m
    print('refit', k, flush=True)
pickle.dump({'members': final, 'names': MEMBERS}, open(W + 'model_cls_final.pkl', 'wb'))
print('saved production ensemble')
