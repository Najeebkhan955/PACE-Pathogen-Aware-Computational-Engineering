"""Shared figure style: clean, high-contrast, publication-grade."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.colors as mcolors

FIG = '/home/claude/figures/'
import os
os.makedirs(FIG, exist_ok=True)

INK = '#12161c'
MUTED = '#6b7684'
GRID = '#e3e7ec'
PALETTE = ['#2b6cb0', '#e05252', '#2f9e6f', '#d98324', '#7b5ea7', '#0e7c86', '#b03060', '#5a6b7c']
GRAM_NEG = '#2b6cb0'
GRAM_POS = '#e05252'
SEQ_CMAP = mcolors.LinearSegmentedColormap.from_list(
    'seq', ['#f7fafc', '#c5dced', '#7fb3d8', '#3f86bd', '#1a4e80'])
DIV_CMAP = mcolors.LinearSegmentedColormap.from_list(
    'div', ['#1a4e80', '#7fb3d8', '#f4f6f8', '#eb9b8f', '#b5342a'])

rcParams.update({
    'figure.dpi': 120, 'savefig.dpi': 600,
    'font.family': 'DejaVu Sans', 'font.size': 9,
    'axes.titlesize': 10.5, 'axes.titleweight': 'bold', 'axes.labelsize': 9.5,
    'axes.edgecolor': '#c6ccd4', 'axes.linewidth': 0.9, 'axes.labelcolor': INK,
    'axes.facecolor': 'white', 'figure.facecolor': 'white',
    'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.7, 'grid.alpha': 1.0,
    'axes.axisbelow': True, 'axes.spines.top': False, 'axes.spines.right': False,
    'xtick.color': MUTED, 'ytick.color': MUTED, 'text.color': INK,
    'xtick.labelsize': 8.5, 'ytick.labelsize': 8.5,
    'legend.frameon': False, 'legend.fontsize': 8.5,
    'lines.linewidth': 1.9, 'lines.solid_capstyle': 'round',
    'savefig.bbox': 'tight', 'savefig.pad_inches': 0.25,
})


def panel(ax, letter, dx=-0.085, dy=1.06):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=12.5,
            fontweight='bold', va='top', ha='left', color=INK)


def save(fig, name):
    fig.savefig(FIG + name + '.png', facecolor='white')
    fig.savefig(FIG + name + '.pdf', facecolor='white')
    plt.close(fig)
    print('saved', name, flush=True)
