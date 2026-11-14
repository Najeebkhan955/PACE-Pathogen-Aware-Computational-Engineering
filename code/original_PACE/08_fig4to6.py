import numpy as np, pandas as pd, json, sys, os
sys.path.insert(0, os.path.dirname(__file__))
from style import *
from features import FEATURE_NAMES
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from sklearn.metrics import confusion_matrix

W = '/home/claude/work/'
TARGETS = ['E. faecium', 'S. aureus', 'K. pneumoniae', 'A. baumannii', 'P. aeruginosa', 'E. cloacae']
GRAM = {'E. faecium': 'pos', 'S. aureus': 'pos', 'K. pneumoniae': 'neg',
        'A. baumannii': 'neg', 'P. aeruginosa': 'neg', 'E. cloacae': 'neg'}
res = pd.read_csv(W + 'cls_metrics.csv')

# ============================================================ FIGURE 4 robustness
sb = pd.read_csv(W + 'similarity_bins.csv')
yte, pte = np.load(W + 'yte.npy'), np.load(W + 'pte_ENS.npy')
yex, pex = np.load(W + 'yex.npy'), np.load(W + 'pex_ENS.npy')
sim = np.load(W + 'test_similarity.npy')

fig, axes = plt.subplots(2, 2, figsize=(11, 8.4))
a = axes[0, 0]
x = np.arange(len(sb))
a.bar(x, sb.n, color='#dfe6ee', width=.66, label='peptides in bin')
a.set_ylabel('Peptides', color=MUTED); a.set_xticks(x); a.set_xticklabels(sb.bin)
a.set_xlabel('Maximum 3-mer Jaccard similarity to any training peptide')
a2 = a.twinx(); a2.grid(False)
a2.plot(x, sb.ROC_AUC, 'o-', color=PALETTE[0], label='ROC-AUC', ms=6)
a2.plot(x, sb.MCC, 's-', color=PALETTE[1], label='MCC', ms=5.5)
for xi, v in zip(x, sb.ROC_AUC):
    a2.annotate(f'{v:.3f}', (xi, v), textcoords='offset points', xytext=(0, 8),
                ha='center', fontsize=7.6, color=PALETTE[0])
a2.set_ylim(.55, 1.06); a2.set_ylabel('Performance')
h1, l1 = a.get_legend_handles_labels(); h2, l2 = a2.get_legend_handles_labels()
a.legend(h1 + h2, l1 + l2, loc='lower right')
panel(a, 'a'); a.set_title('Performance is retained on low-similarity peptides', loc='left')

a = axes[0, 1]
order = ['Logistic regression', 'Neural network (MLP)', 'Random forest', 'SVM (RBF)',
         'Extra trees', 'Gradient boosting', 'Soft-voting ensemble']
t_auc = [res[(res.Model == m) & (res.Set == 'Held-out test')].ROC_AUC.iloc[0] for m in order]
e_auc = [res[(res.Model == m) & (res.Set == 'External validation')].ROC_AUC.iloc[0] for m in order]
y = np.arange(len(order))
a.barh(y - .2, t_auc, height=.38, color=PALETTE[0], label='Held-out test')
a.barh(y + .2, e_auc, height=.38, color=PALETTE[3], label='External validation')
for yi, (t, e) in enumerate(zip(t_auc, e_auc)):
    a.text(t + .004, yi - .2, f'{t:.3f}', va='center', fontsize=7.4, color=MUTED)
    a.text(e + .004, yi + .2, f'{e:.3f}', va='center', fontsize=7.4, color=MUTED)
a.set_yticks(y); a.set_yticklabels(order, fontsize=8.2); a.invert_yaxis()
a.set_xlim(.75, 1.02); a.set_xlabel('ROC-AUC'); a.legend(loc='lower right')
panel(a, 'b'); a.set_title('Generalisation to an independent corpus', loc='left')

a = axes[1, 0]
cm = confusion_matrix(yte, (pte >= .5).astype(int))
im = a.imshow(cm, cmap=SEQ_CMAP)
for i in range(2):
    for j in range(2):
        a.text(j, i, f'{cm[i, j]:,}', ha='center', va='center', fontsize=15,
               fontweight='bold', color='white' if cm[i, j] > cm.max() * .55 else INK)
a.set_xticks([0, 1]); a.set_xticklabels(['predicted\nnon-AMP', 'predicted\nAMP'])
a.set_yticks([0, 1]); a.set_yticklabels(['true\nnon-AMP', 'true\nAMP'])
a.grid(False); panel(a, 'c')
a.set_title('Ensemble confusion matrix (held-out test)', loc='left')

a = axes[1, 1]
bins = np.linspace(0, 1, 46)
a.hist(pex[yex == 0], bins=bins, color=PALETTE[1], alpha=.75, label='non-AMP')
a.hist(pex[yex == 1], bins=bins, color=PALETTE[0], alpha=.75, label='AMP')
a.axvline(.5, color=INK, ls='--', lw=1.2)
a.set_xlabel('Ensemble P(AMP)'); a.set_ylabel('Count'); a.legend()
panel(a, 'd'); a.set_title('Score separation on the external set', loc='left')
fig.suptitle('Robustness and independent validation', fontsize=13, fontweight='bold', x=.06, ha='left', y=1.0)
fig.tight_layout(); save(fig, 'Figure4_robustness')

# ============================================================ FIGURE 5 interpretation
imp = pd.read_csv(W + 'perm_importance.csv').head(18).iloc[::-1]
sv = np.load(W + 'shap_values.npy'); Xs = np.load(W + 'shap_X.npy')
ms = np.abs(sv).mean(0)
top = np.argsort(ms)[::-1][:14][::-1]

fig, axes = plt.subplots(2, 2, figsize=(11.4, 9.2))
a = axes[0, 0]
a.barh(range(len(imp)), imp.importance, xerr=imp.sd, color=PALETTE[0],
       error_kw=dict(ecolor=MUTED, lw=.9, capsize=2), height=.68)
a.set_yticks(range(len(imp))); a.set_yticklabels(imp.feature, fontsize=7.8)
a.set_xlabel('Decrease in ROC-AUC when permuted'); panel(a, 'a')
a.set_title('Permutation importance (gradient boosting)', loc='left')

a = axes[0, 1]
rng = np.random.default_rng(3)
for row, fi in enumerate(top):
    v = sv[:, fi]; f = Xs[:, fi]
    lo, hi = np.percentile(f, [2, 98])
    col = np.clip((f - lo) / (hi - lo + 1e-9), 0, 1)
    a.scatter(v, row + rng.normal(0, .12, len(v)), c=col, cmap=DIV_CMAP, s=5.5, lw=0, alpha=.8)
a.axvline(0, color=MUTED, lw=.9)
a.set_yticks(range(len(top))); a.set_yticklabels([FEATURE_NAMES[i] for i in top], fontsize=7.8)
a.set_xlabel('SHAP value (impact on P(AMP))'); panel(a, 'b')
sm = plt.cm.ScalarMappable(cmap=DIV_CMAP); sm.set_array([])
cb = fig.colorbar(sm, ax=a, fraction=.035, pad=.02, ticks=[0, 1])
cb.ax.set_yticklabels(['low', 'high']); cb.set_label('feature value', fontsize=8)
a.set_title('SHAP attribution (extra trees, 400 test peptides)', loc='left')

a = axes[1, 0]
qi = FEATURE_NAMES.index('net_charge_pH7'); mi = FEATURE_NAMES.index('hydrophobic_moment')
s = a.scatter(Xs[:, qi], sv[:, qi], c=Xs[:, mi], cmap=DIV_CMAP, s=18, lw=.2, edgecolor='white')
a.axhline(0, color=MUTED, lw=.9)
a.set_xlabel('Net charge at pH 7.4'); a.set_ylabel('SHAP value for net charge')
cb = fig.colorbar(s, ax=a, fraction=.035, pad=.02); cb.set_label('$\\mu_H$', fontsize=8)
panel(a, 'c'); a.set_title('Charge contribution is modulated by amphipathicity', loc='left')

a = axes[1, 1]
block = ['Physicochemical' if f in FEATURE_NAMES[:28] else ('Amino-acid composition'
         if f.startswith('AAC_') else 'CTD') for f in FEATURE_NAMES]
bd = pd.DataFrame({'block': block, 'ms': ms}).groupby('block').ms.agg(['sum', 'mean', 'count'])
bd = bd.reindex(['Physicochemical', 'Amino-acid composition', 'CTD'])
xb = np.arange(3)
a.bar(xb, bd['sum'], color=[PALETTE[0], PALETTE[2], PALETTE[3]], width=.6)
for xi, (v, c) in enumerate(zip(bd['sum'], bd['count'])):
    a.text(xi, v, f'{v:.3f}\n({c} features)', ha='center', va='bottom', fontsize=8, color=MUTED)
a.set_xticks(xb); a.set_xticklabels(bd.index, fontsize=8.5)
a.set_ylabel('Total mean |SHAP|'); a.set_ylim(0, bd['sum'].max() * 1.25)
panel(a, 'd'); a.set_title('Contribution by descriptor family', loc='left')
fig.suptitle('What the model has learned', fontsize=13, fontweight='bold', x=.06, ha='left', y=1.0)
fig.tight_layout(); save(fig, 'Figure5_interpretation')

# ============================================================ FIGURE 6 potency
preds = json.load(open(W + 'reg_preds.json'))
rm = pd.read_csv(W + 'reg_metrics.csv')
corr = pd.read_csv(W + 'mic_cross_corr.csv', index_col=0)

fig = plt.figure(figsize=(12, 9.6))
gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 1.05], hspace=.55, wspace=.52)
for k, t in enumerate(TARGETS):
    a = fig.add_subplot(gs[k // 3, k % 3])
    d = preds[t]; y, p = np.array(d['y']), np.array(d['pred'])
    col = GRAM_NEG if GRAM[t] == 'neg' else GRAM_POS
    a.scatter(y, p, s=11, color=col, alpha=.4, lw=0)
    lim = [min(y.min(), p.min()) - .2, max(y.max(), p.max()) + .2]
    a.plot(lim, lim, ls='--', color=INK, lw=1)
    r2 = rm[(rm.Target == t) & (rm.Model == d['model'])].R2.iloc[0]
    rho = rm[(rm.Target == t) & (rm.Model == d['model'])].Spearman_rho.iloc[0]
    a.text(.04, .96, f'$R^2$ = {r2:.3f}\n$\\rho$ = {rho:.3f}\nn = {len(y):,}', transform=a.transAxes,
           va='top', fontsize=8.2, bbox=dict(fc='white', ec='#d5dae1', lw=.8, boxstyle='round,pad=0.4'))
    a.set_title(t, loc='left', style='italic', fontsize=10)
    a.set_xlim(lim); a.set_ylim(lim)
    if k // 3 == 1:
        a.set_xlabel('Measured $\\log_{10}$ MIC')
    if k % 3 == 0:
        a.set_ylabel('Predicted (5-fold CV)')
    if k == 0:
        panel(a, 'a', dx=-.22, dy=1.22)

a = fig.add_subplot(gs[2, 0])
best = rm.loc[rm.groupby('Target').R2.idxmax()].set_index('Target').reindex(TARGETS)
cols = [GRAM_NEG if GRAM[t] == 'neg' else GRAM_POS for t in TARGETS]
b = a.bar(range(6), best.R2, color=cols, width=.66)
for r, (v, n) in zip(b, zip(best.R2, best.n)):
    a.text(r.get_x() + r.get_width() / 2, v + .012, f'{v:.3f}\nn={int(n):,}', ha='center',
           fontsize=7.4, color=MUTED)
a.set_xticks(range(6)); a.set_xticklabels([t.replace('. ', '.') for t in TARGETS],
                                          style='italic', fontsize=7.4, rotation=38, ha='right')
a.set_ylabel('Cross-validated $R^2$'); a.set_ylim(0, best.R2.max() * 1.32)
panel(a, 'b'); a.set_title('Predictability by target', loc='left')

a = fig.add_subplot(gs[2, 1])
o = ['S. aureus', 'E. faecium', 'K. pneumoniae', 'A. baumannii', 'E. cloacae', 'P. aeruginosa']
C = corr.reindex(index=o, columns=o)
im = a.imshow(C.values, cmap=SEQ_CMAP, vmin=.3, vmax=1)
for i in range(6):
    for j in range(6):
        a.text(j, i, f'{C.values[i, j]:.2f}', ha='center', va='center', fontsize=7.2,
               color='white' if C.values[i, j] > .82 else INK)
a.set_xticks(range(6)); a.set_xticklabels([t.replace('. ', '.') for t in o], rotation=45,
                                          ha='right', style='italic', fontsize=7.4)
a.set_yticks(range(6)); a.set_yticklabels([t.replace('. ', '.') for t in o], style='italic', fontsize=7.4)
a.add_patch(plt.Rectangle((-.5, -.5), 2, 2, fill=False, ec='#e05252', lw=2.4))
a.add_patch(plt.Rectangle((1.5, 1.5), 4, 4, fill=False, ec='#f2b705', lw=2.4))
a.text(.5, -1.15, 'Gram +', ha='center', fontsize=7.4, color='#e05252', fontweight='bold')
a.text(3.5, -1.15, 'Gram \u2212', ha='center', fontsize=7.4, color='#c99a04', fontweight='bold')
a.grid(False); fig.colorbar(im, ax=a, fraction=.04, pad=.02)
panel(a, 'c'); a.set_title('Gram-class structure', loc='left')

a = fig.add_subplot(gs[2, 2])
piv = rm.pivot_table(index='Model', columns='Target', values='R2')
piv = piv.reindex(columns=TARGETS + ['Haemolysis (HC50)'])
piv = piv.reindex(['Ridge', 'SVR (RBF)', 'Gradient boosting', 'Random forest', 'Extra trees'])
V = np.clip(piv.values, 0, None)
im = a.imshow(V, cmap=SEQ_CMAP, vmin=0, vmax=.62, aspect='auto')
for i in range(piv.shape[0]):
    for j in range(piv.shape[1]):
        raw = piv.values[i, j]
        lab = '< 0' if raw < 0 else f'{raw:.2f}'
        a.text(j, i, lab, ha='center', va='center', fontsize=7,
               color='white' if V[i, j] > .42 else INK)
a.set_xticks(range(piv.shape[1]))
a.set_xticklabels([c.replace('. ', '.').replace(' (HC50)', '') for c in piv.columns],
                  rotation=45, ha='right', fontsize=7.2)
a.set_yticks(range(piv.shape[0])); a.set_yticklabels(piv.index, fontsize=7.8)
a.grid(False); fig.colorbar(im, ax=a, fraction=.04, pad=.02)
panel(a, 'd'); a.set_title('Algorithm $\\times$ endpoint', loc='left')
fig.suptitle('Pathogen-resolved potency modelling', fontsize=13, fontweight='bold', x=.06, ha='left', y=.98)
save(fig, 'Figure6_potency')
