import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from style import *
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

fig, ax = plt.subplots(figsize=(12.6, 7.0))
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis('off')
LH = 4.4   # line height


def box(x, y, w, title, lines, fc):
    h = 9.5 + LH * len(lines)
    ax.add_patch(FancyBboxPatch((x, y - h), w, h, boxstyle='round,pad=0.7,rounding_size=1.6',
                                fc=fc, ec=fc, lw=1.2))
    ax.text(x + w / 2, y - 4.4, title, ha='center', va='center', fontsize=9.8,
            fontweight='bold', color='white')
    for i, ln in enumerate(lines):
        ax.text(x + w / 2, y - 9.8 - i * LH, ln, ha='center', va='center',
                fontsize=8.1, color='white', alpha=.95)
    return (x, y, w, h)


def arrow(p1, p2):
    ax.add_patch(FancyArrowPatch(p1, p2, arrowstyle='-|>', mutation_scale=14,
                                 lw=1.6, color='#96a0ad', shrinkA=3, shrinkB=3))


ax.text(0, 98, 'An integrated in-silico pipeline for ESKAPE-directed antimicrobial peptides',
        ha='left', va='top', fontsize=13.5, fontweight='bold', color=INK)

d = box(0, 88, 22, 'PUBLIC DATA SOURCES',
        ['AMPlify benchmark', '7,839 peptides', '', 'amPEPpy M-model',
         '1,848 external', '', 'GRAMPA MIC records', '44,151 across ESKAPE',
         '', 'Hemolytik HC50', '1,271 peptides'], '#1a3350')
f = box(27, 74, 18, 'DESCRIPTOR SPACE',
        ['90 features', '', '28 physicochemical', '20 amino-acid composition',
         '42 CTD (7 attributes)'], '#2b6cb0')
s1 = box(50, 91, 22, 'STAGE 1 · RECOGNITION',
         ['6 algorithms + ensemble', 'AMP vs non-AMP', 'external validation'], '#2f9e6f')
s2 = box(50, 68, 22, 'STAGE 2 · POTENCY',
         ['6 pathogen-specific', '$\\log_{10}$ MIC regressors', 'Gram-resolved'], '#0e7c86')
s3 = box(50, 45, 22, 'STAGE 3 · SAFETY',
         ['$\\log_{10}$ HC$_{50}$ model', 'selectivity index', 'erythrocyte data'], '#d98324')
s4 = box(77, 82, 23, 'STAGE 4 · DE NOVO DESIGN',
         ['multi-objective', 'genetic algorithm', '', '9,315 candidates',
          'over 45 generations'], '#7b5ea7')
s5 = box(77, 44, 23, 'PRIORITISED OUTPUT',
         ['novelty + diversity filter', '', '20 leads: PACE-01 … 20',
          '', 'no wet-lab validation'], '#b03060')

arrow((22, 76), (27, 70))
arrow((45, 70), (50, 84))
arrow((45, 68), (50, 61))
arrow((45, 66), (50, 38))
arrow((72, 84), (77, 74))
arrow((72, 61), (77, 66))
arrow((72, 38), (77, 58))
arrow((88.5, 50.5), (88.5, 44))

ax.text(0, 4, 'Every dataset is public and version-controlled; every model is trained, validated and applied '
              'within one reproducible workflow.',
        ha='left', fontsize=8.2, color=MUTED, style='italic')
save(fig, 'Figure1_workflow')
