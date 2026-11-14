import numpy as np, pandas as pd, json, sys, os
sys.path.insert(0, os.path.dirname(__file__))
from style import *
from features import FEATURE_NAMES, featurize_many, EIS
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors

W = '/home/claude/work/'
TARGETS = ['E. faecium', 'S. aureus', 'K. pneumoniae', 'A. baumannii', 'P. aeruginosa', 'E. cloacae']
GRAM = {'E. faecium': 'pos', 'S. aureus': 'pos', 'K. pneumoniae': 'neg',
        'A. baumannii': 'neg', 'P. aeruginosa': 'neg', 'E. cloacae': 'neg'}

hist = pd.read_csv(W + 'ga_history.csv')
cand = pd.read_csv(W + 'ga_all_candidates.csv')
lead = pd.read_csv(W + 'designed_leads.csv')
fun = pd.read_csv(W + 'design_funnel.csv')
Xn = np.load(W + 'X_natural.npy'); Xo = np.load(W + 'X_nonamp.npy'); Xl = np.load(W + 'X_leads.npy')
Xc = featurize_many(cand.sequence.tolist())

# ============================================================ FIGURE 7
fig, axes = plt.subplots(2, 3, figsize=(13.4, 8.2))
a = axes[0, 0]
a.plot(hist.generation, hist.best_fitness, color=PALETTE[4], label='best individual')
a.plot(hist.generation, hist.mean_fitness, color=PALETTE[0], label='population mean')
a.fill_between(hist.generation, hist.mean_fitness, hist.best_fitness, color=PALETTE[4], alpha=.12)
a.set_xlabel('Generation'); a.set_ylabel('Multi-objective fitness'); a.legend(loc='lower right')
panel(a, 'a'); a.set_title('Optimisation trajectory', loc='left')

a = axes[0, 1]
a.plot(hist.generation, hist.mean_p_amp, color=PALETTE[2], label='mean P(AMP)')
a.set_ylabel('P(AMP)', color=PALETTE[2]); a.set_xlabel('Generation'); a.set_ylim(0, 1.02)
a2 = a.twinx(); a2.grid(False)
a2.plot(hist.generation, hist.mean_mic, color=PALETTE[1], label='mean predicted $\\log_{10}$ MIC')
a2.plot(hist.generation, hist.mean_ti, color=PALETTE[3], label='mean selectivity ($\\log_{10}$)')
a2.set_ylabel('$\\log_{10}$ units')
h1, l1 = a.get_legend_handles_labels(); h2, l2 = a2.get_legend_handles_labels()
a.legend(h1 + h2, l1 + l2, loc='center right', fontsize=7.8)
panel(a, 'b'); a.set_title('Objectives improve jointly', loc='left')

a = axes[0, 2]
y = np.arange(len(fun))
a.barh(y, fun.n, color=[PALETTE[0]] * (len(fun) - 1) + [PALETTE[4]], height=.62)
for yi, v in zip(y, fun.n):
    a.text(v * 1.08, yi, f'{v:,}', va='center', fontsize=7.8, color=MUTED)
a.set_yticks(y); a.set_yticklabels([s.replace(' against all six ESKAPE targets', '\nvs all six targets')
                                    .replace(' (pairwise 3-mer Jaccard <= 0.45)', '\n(3-mer Jaccard \u2264 0.45)')
                                    for s in fun.step], fontsize=7.4)
a.invert_yaxis(); a.set_xscale('log'); a.set_xlim(8, 3e4)
a.set_xlabel('Candidates remaining (log scale)')
panel(a, 'c'); a.set_title('Prioritisation funnel', loc='left')

# --- chemical space
sc = StandardScaler().fit(np.vstack([Xn, Xo]))
pca = PCA(n_components=2).fit(sc.transform(np.vstack([Xn, Xo])))
Pn, Po, Pc, Pl = (pca.transform(sc.transform(x)) for x in (Xn, Xo, Xc, Xl))
a = axes[1, 0]
a.scatter(Po[:, 0], Po[:, 1], s=6, c='#c9d2dc', lw=0, alpha=.6, label='non-AMP')
a.scatter(Pn[:, 0], Pn[:, 1], s=6, c=PALETTE[0], lw=0, alpha=.35, label='natural AMP')
a.scatter(Pc[:, 0], Pc[:, 1], s=5, c=PALETTE[2], lw=0, alpha=.25, label='GA candidates')
a.scatter(Pl[:, 0], Pl[:, 1], s=52, c=PALETTE[4], marker='*', edgecolor='white', lw=.6,
          label='prioritised leads', zorder=6)
a.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
a.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
a.set_xlim(np.percentile(Po[:, 0], .3), np.percentile(Pn[:, 0], 99.9) + 3)
a.legend(markerscale=1.8, loc='upper left', fontsize=7.6)
panel(a, 'd'); a.set_title('Descriptor space occupancy', loc='left')

# --- applicability domain
nn = NearestNeighbors(n_neighbors=5).fit(sc.transform(Xn))
dn = nn.kneighbors(sc.transform(Xn))[0][:, 1:].mean(1)
dl = nn.kneighbors(sc.transform(Xl))[0].mean(1)
dc = nn.kneighbors(sc.transform(Xc[:3000]))[0].mean(1)
a = axes[1, 1]
bins = np.linspace(0, max(dl.max(), np.percentile(dn, 99.5)) * 1.05, 45)
a.hist(dn, bins=bins, color=PALETTE[0], alpha=.7, density=True, label='natural AMPs')
a.hist(dc, bins=bins, color=PALETTE[2], alpha=.55, density=True, label='GA candidates')
for v in dl:
    a.axvline(v, color=PALETTE[4], lw=1.1, alpha=.8)
a.axvline(np.percentile(dn, 95), color=INK, ls='--', lw=1.3, label='95th pct. of natural AMPs')
a.plot([], [], color=PALETTE[4], lw=1.1, label='prioritised leads')
a.set_xlabel('Mean distance to 5 nearest natural AMPs\n(standardised descriptor space)')
a.set_ylabel('Density'); a.set_ylim(0, 1.15); a.legend(fontsize=7.2, loc='upper right')
panel(a, 'e'); a.set_title('Applicability domain', loc='left')

a = axes[1, 2]
bins = np.linspace(0, .6, 40)
a.hist(cand.sim, bins=bins, color=PALETTE[2], alpha=.7, label='all candidates')
a.hist(lead.sim, bins=bins, color=PALETTE[4], alpha=.95, label='prioritised leads')
a.axvline(.55, color=INK, ls='--', lw=1.2, label='novelty threshold')
a.set_xlabel('Max 3-mer Jaccard similarity to known peptides')
a.set_ylabel('Count'); a.legend(fontsize=7.6, loc='upper right')
panel(a, 'f'); a.set_title('Sequence novelty', loc='left')
fig.suptitle('De novo design under multi-objective optimisation', fontsize=13,
             fontweight='bold', x=.055, ha='left', y=1.0)
fig.tight_layout(); save(fig, 'Figure7_design')

# ============================================================ FIGURE 8 leads
fig = plt.figure(figsize=(13.4, 8.8))
gs = fig.add_gridspec(2, 3, height_ratios=[1.12, 1], hspace=.42, wspace=.34)

a = fig.add_subplot(gs[0, :2])
M = lead[['mic_' + t for t in TARGETS]].values
im = a.imshow(M.T, cmap=SEQ_CMAP.reversed(), aspect='auto', vmin=0, vmax=1.2)
a.set_yticks(range(6)); a.set_yticklabels(TARGETS, style='italic', fontsize=8.4)
a.set_xticks(range(len(lead))); a.set_xticklabels(lead.peptide_id, rotation=90, fontsize=7)
for i in range(6):
    for j in range(len(lead)):
        a.text(j, i, f'{10**M[j, i]:.0f}', ha='center', va='center', fontsize=6.2,
               color='white' if M[j, i] < .45 else INK)
a.grid(False)
cb = fig.colorbar(im, ax=a, fraction=.028, pad=.015)
cb.set_label('predicted $\\log_{10}$ MIC ($\\mu$M)\ndarker = more potent', fontsize=7.6)
panel(a, 'a', dx=-.07)
a.set_title('Predicted potency of prioritised leads (cell labels: MIC in $\\mu$M)', loc='left')

# helical wheel
a = fig.add_subplot(gs[0, 2])
seq = lead.sequence.iloc[0]
CLASSC = {'basic': '#2b6cb0', 'acidic': '#e05252', 'polar': '#9aa6b4', 'hydrophobic': '#d98324'}
def cls_of(c):
    if c in 'KRH': return 'basic'
    if c in 'DE': return 'acidic'
    if c in 'AILVFWMYC': return 'hydrophobic'
    return 'polar'
n = min(len(seq), 18)
ang = np.deg2rad(100)
hx = hy = 0.0
for i in range(n):
    th = np.pi / 2 - i * ang
    r = 1 - .028 * i
    x, y = r * np.cos(th), r * np.sin(th)
    a.add_patch(plt.Circle((x, y), .115, fc=CLASSC[cls_of(seq[i])], ec='white', lw=1.3, zorder=3))
    a.text(x, y, seq[i], ha='center', va='center', fontsize=8.4, color='white',
           fontweight='bold', zorder=4)
    hx += EIS[seq[i]] * np.cos(th); hy += EIS[seq[i]] * np.sin(th)
norm = np.hypot(hx, hy)
a.arrow(0, 0, hx / norm * .62, hy / norm * .62, width=.022, color=INK, zorder=5,
        length_includes_head=True, head_width=.09)
a.set_xlim(-1.3, 1.3); a.set_ylim(-1.3, 1.3); a.set_aspect('equal'); a.axis('off')
a.legend(handles=[Line2D([], [], marker='o', ls='', color=v, label=k, markersize=7)
                  for k, v in CLASSC.items()], loc='lower center', ncol=2,
         bbox_to_anchor=(.5, -.13), fontsize=7.6)
panel(a, 'b', dx=-.10, dy=1.18)
a.set_title(f'{lead.peptide_id.iloc[0]} helical wheel  (arrow: hydrophobic moment)',
            loc='center', fontsize=9.0, pad=18)

# radar
a = fig.add_subplot(gs[1, 0], projection='polar')
keys = ['net_charge_pH7', 'hydrophobic_moment', 'frac_basic', 'aliphatic_index',
        'frac_aromatic', 'charge_density']
labels = ['net\ncharge', '$\\mu_H$', 'basic\nfraction', 'aliphatic\nindex', 'aromatic\nfraction', 'charge\ndensity']
idx = [FEATURE_NAMES.index(k) for k in keys]
ref = np.vstack([Xn, Xo])
lo = np.percentile(ref[:, idx], 2, axis=0); hi = np.percentile(ref[:, idx], 98, axis=0)
def norm_row(v):
    return np.clip((v - lo) / (hi - lo + 1e-9), 0, 1.35)
th = np.linspace(0, 2 * np.pi, len(keys), endpoint=False).tolist(); th += th[:1]
for vals, c, lab in [(Xo[:, idx].mean(0), '#9aa6b4', 'non-AMP'),
                     (Xn[:, idx].mean(0), PALETTE[0], 'natural AMP'),
                     (Xl[:, idx].mean(0), PALETTE[4], 'designed leads')]:
    v = norm_row(vals).tolist(); v += v[:1]
    a.plot(th, v, color=c, lw=2); a.fill(th, v, color=c, alpha=.16)
a.set_xticks(th[:-1]); a.set_xticklabels(labels, fontsize=7.4)
a.set_yticks([0, .5, 1]); a.set_yticklabels([], fontsize=6.6)
a.set_ylim(0, 1.38); a.grid(color=GRID)
a.legend(handles=[Line2D([], [], color='#9aa6b4', lw=2.4, label='non-AMP'),
                  Line2D([], [], color=PALETTE[0], lw=2.4, label='natural AMP'),
                  Line2D([], [], color=PALETTE[4], lw=2.4, label='designed leads')],
         loc='lower center', bbox_to_anchor=(.5, -.30), fontsize=7.4, ncol=1)
a.set_title('Property profile', loc='left', pad=18)
a.text(-.24, 1.12, 'c', transform=a.transAxes, fontsize=12.5, fontweight='bold', color=INK)

# position-specific residue frequency
a = fig.add_subplot(gs[1, 1:])
L = 28
counts = np.zeros((L, 20))
AAs = list("ACDEFGHIKLMNPQRSTVWY")
for s in lead.sequence:
    for i in range(min(L, len(s))):
        counts[i, AAs.index(s[i])] += 1
freq = counts / counts.sum(1, keepdims=True)
grp = {'basic': 'KRH', 'acidic': 'DE', 'hydrophobic': 'AILVFWMYC', 'polar': 'GNQSTP'}
bottom = np.zeros(L)
for g, mem in grp.items():
    v = freq[:, [AAs.index(c) for c in mem if c in AAs]].sum(1)
    a.bar(range(1, L + 1), v, bottom=bottom, width=.82,
          color=CLASSC[g], label=g)
    bottom += v
a.set_xlabel('Position'); a.set_ylabel('Residue frequency across the 20 leads')
a.set_xlim(.4, L + .6); a.set_ylim(0, 1.18)
a.legend(ncol=4, fontsize=7.8, loc='upper center', framealpha=.95, frameon=True)
panel(a, 'd', dx=-.05)
a.set_title('Positional composition of the designed library', loc='left')
fig.suptitle('Characterisation of the prioritised peptide library', fontsize=13,
             fontweight='bold', x=.055, ha='left', y=.99)
save(fig, 'Figure8_leads')
