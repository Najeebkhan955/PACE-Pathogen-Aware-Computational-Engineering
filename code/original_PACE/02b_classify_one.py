"""Stage 1 runner: evaluate ONE classifier and checkpoint its results."""
import numpy as np, pandas as pd, json, os, sys, pickle, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(__file__))

from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier, HistGradientBoostingClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import (roc_auc_score, average_precision_score, matthews_corrcoef,
                             accuracy_score, f1_score, confusion_matrix, roc_curve,
                             precision_recall_curve)

W = '/home/claude/work/'
SEED = 2024
NFOLD = 5


def build(name):
    return {
        'LR':   Pipeline([('s', StandardScaler()), ('m', LogisticRegression(max_iter=3000, random_state=SEED))]),
        'SVM':  Pipeline([('s', StandardScaler()), ('m', SVC(C=4, gamma='scale', probability=True, random_state=SEED))]),
        'RF':   RandomForestClassifier(n_estimators=400, n_jobs=-1, random_state=SEED),
        'ET':   ExtraTreesClassifier(n_estimators=400, n_jobs=-1, random_state=SEED),
        'GB':   HistGradientBoostingClassifier(max_iter=400, learning_rate=0.07, random_state=SEED),
        'MLP':  Pipeline([('s', StandardScaler()), ('m', MLPClassifier(hidden_layer_sizes=(128, 64), alpha=1e-3,
                                                                       max_iter=400, early_stopping=True, random_state=SEED))]),
    }[name]


LABEL = {'LR': 'Logistic regression', 'SVM': 'SVM (RBF)', 'RF': 'Random forest',
         'ET': 'Extra trees', 'GB': 'Gradient boosting', 'MLP': 'Neural network (MLP)'}


def metrics(y, p, thr=0.5):
    yh = (p >= thr).astype(int)
    tn, fp, fn, tp = confusion_matrix(y, yh).ravel()
    return dict(Accuracy=accuracy_score(y, yh), Sensitivity=tp / (tp + fn),
                Specificity=tn / (tn + fp), Precision=tp / (tp + fp) if (tp + fp) else 0.0,
                F1=f1_score(y, yh), MCC=matthews_corrcoef(y, yh),
                ROC_AUC=roc_auc_score(y, p), PR_AUC=average_precision_score(y, p))


def main(key):
    Xtr, ytr = np.load(W + 'Xtr.npy'), np.load(W + 'ytr.npy')
    Xte, yte = np.load(W + 'Xte.npy'), np.load(W + 'yte.npy')
    Xex, yex = np.load(W + 'Xex.npy'), np.load(W + 'yex.npy')
    name = LABEL[key]
    cv = StratifiedKFold(n_splits=NFOLD, shuffle=True, random_state=SEED)

    p = cross_val_predict(build(key), Xtr, ytr, cv=cv, method='predict_proba')[:, 1]
    rows = [dict(Model=name, Set=f'{NFOLD}-fold CV (train)', **metrics(ytr, p))]

    # per-fold spread
    folds = []
    for fi, (a, b) in enumerate(cv.split(Xtr, ytr), 1):
        m = build(key); m.fit(Xtr[a], ytr[a])
        folds.append(dict(Model=name, fold=fi, **metrics(ytr[b], m.predict_proba(Xtr[b])[:, 1])))

    mdl = build(key); mdl.fit(Xtr, ytr)
    pte = mdl.predict_proba(Xte)[:, 1]
    pex = mdl.predict_proba(Xex)[:, 1]
    rows.append(dict(Model=name, Set='Held-out test', **metrics(yte, pte)))
    rows.append(dict(Model=name, Set='External validation', **metrics(yex, pex)))

    fpr, tpr, _ = roc_curve(yte, pte)
    pr, rc, _ = precision_recall_curve(yte, pte)
    curves = json.load(open(W + 'cls_curves.json')) if os.path.exists(W + 'cls_curves.json') else {}
    curves[name] = dict(fpr=fpr[::max(1, len(fpr)//400)].tolist(), tpr=tpr[::max(1, len(tpr)//400)].tolist(),
                        prec=pr[::max(1, len(pr)//400)].tolist(), rec=rc[::max(1, len(rc)//400)].tolist())
    json.dump(curves, open(W + 'cls_curves.json', 'w'))

    np.save(W + f'oof_{key}.npy', p)
    np.save(W + f'pte_{key}.npy', pte)
    for fn, df in [('cls_metrics.csv', pd.DataFrame(rows)), ('cls_folds.csv', pd.DataFrame(folds))]:
        df.to_csv(W + fn, mode='a', header=not os.path.exists(W + fn), index=False)
    pickle.dump(mdl, open(W + f'clsmodel_{key}.pkl', 'wb'))
    print(key, 'CV MCC %.3f AUC %.3f | test MCC %.3f | ext MCC %.3f' %
          (rows[0]['MCC'], rows[0]['ROC_AUC'], rows[1]['MCC'], rows[2]['MCC']), flush=True)


if __name__ == '__main__':
    main(sys.argv[1])
