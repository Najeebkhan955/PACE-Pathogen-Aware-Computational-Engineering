"""Stage 2/3 runner: fit potency (log10 MIC) or haemolysis (log10 HC50) models for ONE endpoint."""
import numpy as np, pandas as pd, json, os, sys, pickle, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(__file__))
from features import featurize_many

from sklearn.ensemble import RandomForestRegressor, HistGradientBoostingRegressor, ExtraTreesRegressor
from sklearn.svm import SVR
from sklearn.linear_model import RidgeCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
from scipy.stats import pearsonr, spearmanr

W = '/home/claude/work/'
SEED, NFOLD = 2024, 5
TARGETS = ['E. faecium', 'S. aureus', 'K. pneumoniae', 'A. baumannii', 'P. aeruginosa', 'E. cloacae']


def build():
    return {
        'Ridge':             Pipeline([('s', StandardScaler()), ('m', RidgeCV(alphas=np.logspace(-3, 3, 25)))]),
        'SVR (RBF)':         Pipeline([('s', StandardScaler()), ('m', SVR(C=8, gamma='scale', epsilon=0.1))]),
        'Random forest':     RandomForestRegressor(n_estimators=300, n_jobs=-1, random_state=SEED),
        'Extra trees':       ExtraTreesRegressor(n_estimators=300, n_jobs=-1, random_state=SEED),
        'Gradient boosting': HistGradientBoostingRegressor(max_iter=400, learning_rate=0.06, random_state=SEED),
    }


def run(tag, seqs, y):
    X = featurize_many(seqs)
    cv = KFold(n_splits=NFOLD, shuffle=True, random_state=SEED)
    rows, best, best_r2, best_pred = [], None, -9e9, None
    for name, mdl in build().items():
        p = cross_val_predict(mdl, X, y, cv=cv)
        r2 = r2_score(y, p)
        rows.append(dict(Target=tag, Model=name, n=len(y), R2=r2,
                         RMSE=float(np.sqrt(mean_squared_error(y, p))), MAE=mean_absolute_error(y, p),
                         Pearson_r=pearsonr(y, p)[0], Spearman_rho=spearmanr(y, p)[0]))
        if r2 > best_r2:
            best, best_r2, best_pred = name, r2, p
        print(f'  {name:18s} R2={r2:.3f}', flush=True)
    pd.DataFrame(rows).to_csv(W + 'reg_metrics.csv', mode='a',
                              header=not os.path.exists(W + 'reg_metrics.csv'), index=False)
    preds = json.load(open(W + 'reg_preds.json')) if os.path.exists(W + 'reg_preds.json') else {}
    preds[tag] = dict(y=[round(v, 4) for v in y.tolist()],
                      pred=[round(v, 4) for v in best_pred.tolist()], model=best)
    json.dump(preds, open(W + 'reg_preds.json', 'w'))
    m = build()[best]; m.fit(X, y)
    key = tag.replace('. ', '_').replace(' ', '_').replace('(', '').replace(')', '')
    pickle.dump({'model': m, 'name': best, 'r2': best_r2}, open(W + f'reg_{key}.pkl', 'wb'))
    print(f'{tag}: best={best} R2={best_r2:.3f} n={len(y)}', flush=True)


if __name__ == '__main__':
    which = sys.argv[1]
    if which == 'hem':
        h = pd.read_csv(W + 'hemolysis.csv')
        run('Haemolysis (HC50)', h.sequence, h.log10_HC50.values)
    else:
        t = TARGETS[int(which)]
        mic = pd.read_csv(W + 'mic_eskape.csv')
        sub = mic[mic.target == t].reset_index(drop=True)
        run(t, sub.sequence, sub.log10_MIC_uM.values)
