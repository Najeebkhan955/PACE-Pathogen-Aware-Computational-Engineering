import numpy as np, pandas as pd, json, sys, os
sys.path.insert(0, os.path.dirname(__file__))
from style import *
from features import FEATURE_NAMES, featurize_many
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

W = '/home/claude/work/'
TARGETS = ['E. faecium', 'S. aureus', 'K. pneumoniae', 'A. baumannii', 'P. aeruginosa', 'E. cloacae']
GRAM = {'E. faecium': 'pos', 'S. aureus': 'pos', 'K. pneumoniae': 'neg',
        'A. baumannii': 'neg', 'P. aeruginosa': 'neg', 'E. cloacae': 'neg'}

# ============================================================ FIGURE 1 workflow
fig, ax = plt.subplots(figsize=(11.5, 6.2))
ax.set_xlim(0, 100); ax.set_ylim(0, 60); ax.axis('off')

def box(x, y, w, h, title, lines, fc, ec, tc='white'):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.6,rounding_size=1.4',
                                fc=fc, ec=ec, lw=1.4))
    ax.text(x + w / 2, y + h - 2.6, title, ha='center', va='top', fontsize=9.6,
            fontweight='bold', color=tc)
    for i, ln in enumerate(lines):
        ax.text(x + w / 2, y + h - 6.6 - i * 3.1, ln, ha='center', va='top',
                fontsize=7.9, color=tc, alpha=0.95)

def arrow(x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle='-|>', mutation_scale=13,
                                 lw=1.5, color='#8a94a2', shrinkA=2, shrinkB=2))

box(1, 33, 21, 24, 'PUBLIC DATA SOURCES',
    ['AMPlify benchmark', 'n = 7,839 peptides', '', 'amPEPpy M-model',
     'n = 1,848 (external)', '', 'GRAMPA MIC records', 'n = 44,151 / 6 ESKAPE',
     '', 'Hemolytik HC50', 'n = 1,271'], '#1a3350', '#1a3350')
box(26, 37, 18, 16, 'DESCRIPTORS',
    ['90 features', '', '28 physicochemical', '20 AAC', '42 CTD (7 attributes)'],
    '#2b6cb0', '#2b6cb0')
box(48, 41, 22, 15, 'STAGE 1 — RECOGNITION',
    ['6 algorithms + ensemble', 'AMP vs non-AMP', 'external validation'], '#2f9e6f', '#2f9e6f')
box(48, 23, 22, 15, 'STAGE 2 — POTENCY',
    ['6 pathogen-specific', 'log10 MIC regressors', 'Gram-resolved'], '#0e7c86', '#0e7c86')
box(48, 5, 22, 15, 'STAGE 3 — SAFETY',
    ['log10 HC50 model', 'selectivity index', 'from erythrocyte data'], '#d98324', '#d98324')
box(75, 23, 23, 26, 'STAGE 4 — DE NOVO DESIGN',
    ['multi-objective', 'genetic algorithm', '', '9,315 candidates', '45 generations',
     '', 'novelty + diversity', 'filtering', '', '20 prioritised leads'], '#7b5ea7', '#7b5ea7')
arrow(22, 45, 26, 45); arrow(44, 45, 48, 48); arrow(44, 45, 48, 30); arrow(44, 45, 48, 12)
arrow(70, 48, 75, 40); arrow(70, 30, 75, 34); arrow(70, 12, 75, 28)
ax.text(50, 59.5, 'Integrated in-silico pipeline for ESKAPE-directed antimicrobial peptides',
        ha='center', va='top', fontsize=12, fontweight='bold', color=INK)
ax.text(50, 1.5, 'All datasets are public and version-controlled; all models are trained, validated and '
                 'applied within a single reproducible workflow',
        ha='center', fontsize=8, color=MUTED, style='italic')
save(fig, 'Figure1_workflow')

# ============================================================ FIGURE 2 landscape
tr = pd.read_csv(W + 'cls_train.csv')
Xn = np.load(W + 'X_natural.npy'); Xo = np.load(W + 'X_nonamp.npy')
mic = pd.read_csv(W + 'mic_eskape.csv')
iL, iQ, iM = (FEATURE_NAMES.index(k) for k in ['length', 'net_charge_pH7', 'hydrophobic_moment'])

fig, axes = plt.subplots(2, 2, figsize=(11, 8.4))
a = axes[0, 0]
bins = np.arange(5, 62, 2)
a.hist(Xn[:, iL], bins=bins, color=PALETTE[0], alpha=.85, label='AMP (n=3,124)')
a.hist(Xo[:, iL], bins=bins, color=PALETTE[1], alpha=.55, label='non-AMP (n=3,147)')
a.set_xlabel('Peptide length (residues)'); a.set_ylabel('Count'); a.legend(); panel(a, 'a')
a.set_title('Length distribution of the training corpus', loc='left')

a = axes[0, 1]
a.scatter(Xo[:, iQ], Xo[:, iM], s=7, c=PALETTE[1], alpha=.28, lw=0, label='non-AMP')
a.scatter(Xn[:, iQ], Xn[:, iM], s=7, c=PALETTE[0], alpha=.35, lw=0, label='AMP')
a.set_xlabel('Net charge at pH 7.4'); a.set_ylabel('Hydrophobic moment  $\\mu_H$')
a.set_xlim(-15, 22); a.legend(markerscale=2.2); panel(a, 'b')
a.set_title('Cationic–amphipathic separation', loc='left')

a = axes[1, 0]
cnt = mic.groupby('target').size().reindex(TARGETS)
cols = [GRAM_NEG if GRAM[t] == 'neg' else GRAM_POS for t in TARGETS]
b = a.barh(range(len(TARGETS)), cnt.values, color=cols, height=.66)
a.set_yticks(range(len(TARGETS))); a.set_yticklabels(TARGETS, style='italic')
a.invert_yaxis(); a.set_xlabel('Unique peptides with measured MIC')
for r, v in zip(b, cnt.values):
    a.text(v + 40, r.get_y() + r.get_height() / 2, f'{v:,}', va='center', fontsize=8, color=MUTED)
a.set_xlim(0, cnt.max() * 1.16)
from matplotlib.lines import Line2D
a.legend(handles=[Line2D([], [], color=GRAM_NEG, lw=6, label='Gram-negative'),
                  Line2D([], [], color=GRAM_POS, lw=6, label='Gram-positive')], loc='lower right')
panel(a, 'c'); a.set_title('ESKAPE coverage after curation', loc='left')

a = axes[1, 1]
data = [mic[mic.target == t].log10_MIC_uM.values for t in TARGETS]
vp = a.violinplot(data, showextrema=False, widths=.85)
for body, t in zip(vp['bodies'], TARGETS):
    body.set_facecolor(GRAM_NEG if GRAM[t] == 'neg' else GRAM_POS)
    body.set_alpha(.55); body.set_edgecolor('white')
bp = a.boxplot(data, widths=.13, showfliers=False, patch_artist=True)
for p in bp['boxes']:
    p.set_facecolor('white'); p.set_edgecolor(INK); p.set_linewidth(.9)
for k in ['medians', 'whiskers', 'caps']:
    for p in bp[k]:
        p.set_color(INK); p.set_linewidth(.9)
a.set_xticks(range(1, 7)); a.set_xticklabels([t.replace('. ', '.\n') for t in TARGETS],
                                             style='italic', fontsize=7.6)
a.set_ylabel('Measured $\\log_{10}$ MIC ($\\mu$M)'); panel(a, 'd')
a.set_title('Potency distribution by pathogen', loc='left')
fig.suptitle('Curated data landscape', fontsize=13, fontweight='bold', x=.06, ha='left', y=1.0)
fig.tight_layout()
save(fig, 'Figure2_data_landscape')

# ============================================================ FIGURE 3 classifier
res = pd.read_csv(W + 'cls_metrics.csv')
curves = json.load(open(W + 'cls_curves.json'))
order = ['Logistic regression', 'Neural network (MLP)', 'Random forest', 'SVM (RBF)',
         'Extra trees', 'Gradient boosting', 'Soft-voting ensemble']
fig, axes = plt.subplots(2, 2, figsize=(11, 8.6))

a = axes[0, 0]
for i, m in enumerate(order):
    c = curves[m]; auc = res[(res.Model == m) & (res.Set == 'Held-out test')].ROC_AUC.iloc[0]
    lw = 2.6 if 'ensemble' in m else 1.5
    a.plot(c['fpr'], c['tpr'], color=PALETTE[i % len(PALETTE)], lw=lw,
           label=f'{m} ({auc:.3f})', zorder=5 if 'ensemble' in m else 2)
a.plot([0, 1], [0, 1], ls=':', color=MUTED, lw=1)
a.set_xlabel('False positive rate'); a.set_ylabel('True positive rate')
a.legend(title='ROC-AUC', title_fontsize=8.5, loc='lower right'); panel(a, 'a')
a.set_title('Receiver operating characteristic (held-out test)', loc='left')

a = axes[0, 1]
for i, m in enumerate(order):
    c = curves[m]; ap = res[(res.Model == m) & (res.Set == 'Held-out test')].PR_AUC.iloc[0]
    lw = 2.6 if 'ensemble' in m else 1.5
    a.plot(c['rec'], c['prec'], color=PALETTE[i % len(PALETTE)], lw=lw, label=f'{m} ({ap:.3f})',
           zorder=5 if 'ensemble' in m else 2)
a.set_xlabel('Recall'); a.set_ylabel('Precision'); a.set_ylim(.5, 1.02)
a.legend(title='PR-AUC', title_fontsize=8.5, loc='lower left'); panel(a, 'b')
a.set_title('Precision–recall', loc='left')

a = axes[1, 0]
met = ['Accuracy', 'Sensitivity', 'Specificity', 'F1', 'MCC', 'ROC_AUC']
M = res[res.Set == 'Held-out test'].set_index('Model').reindex(order)[met]
im = a.imshow(M.values, cmap=SEQ_CMAP, vmin=.55, vmax=1.0, aspect='auto')
a.set_xticks(range(len(met))); a.set_xticklabels([m.replace('_', '-') for m in met], rotation=30, ha='right')
a.set_yticks(range(len(order))); a.set_yticklabels(order, fontsize=8)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        v = M.values[i, j]
        a.text(j, i, f'{v:.3f}', ha='center', va='center', fontsize=7.6,
               color='white' if v > .87 else INK)
a.grid(False); panel(a, 'c'); a.set_title('Held-out performance matrix', loc='left')
fig.colorbar(im, ax=a, fraction=.035, pad=.02)

a = axes[1, 1]
fo = pd.read_csv(W + 'cls_folds.csv')
grp = [fo[fo.Model == m].MCC.values for m in order[:-1]]
bp = a.boxplot(grp, patch_artist=True, widths=.6, showfliers=False)
for p, c in zip(bp['boxes'], PALETTE):
    p.set_facecolor(c); p.set_alpha(.6); p.set_edgecolor(c)
for k in ['medians', 'whiskers', 'caps']:
    for p in bp[k]:
        p.set_color(INK); p.set_linewidth(1.0)
for i, g in enumerate(grp, 1):
    a.scatter(np.random.normal(i, .06, len(g)), g, s=13, color=INK, alpha=.65, zorder=4, lw=0)
ens = res[(res.Model == 'Soft-voting ensemble') & (res.Set == '5-fold CV (train)')].MCC.iloc[0]
a.axhline(ens, color=PALETTE[6], ls='--', lw=1.6, label=f'Ensemble CV MCC = {ens:.3f}')
a.set_xticks(range(1, len(order))); a.set_xticklabels([m.split(' (')[0] for m in order[:-1]],
                                                      rotation=28, ha='right', fontsize=8)
a.set_ylabel('Matthews correlation coefficient'); a.legend(loc='lower right'); panel(a, 'd')
a.set_title('Cross-validation stability (5 folds)', loc='left')
fig.suptitle('AMP recognition benchmark', fontsize=13, fontweight='bold', x=.06, ha='left', y=1.0)
fig.tight_layout()
save(fig, 'Figure3_classifier_benchmark')
