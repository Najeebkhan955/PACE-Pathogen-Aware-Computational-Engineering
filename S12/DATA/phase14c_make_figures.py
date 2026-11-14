from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
D=Path('/mnt/data/pace_phase14c_analysis')
def read(name):
 with open(D/name,encoding='utf-8',newline='') as f:return list(csv.DictReader(f))
A=read('Phase14C_all7_endpoint_applicability_domain.csv');T=['E. faecium','S. aureus','K. pneumoniae','A. baumannii','P. aeruginosa','E. cloacae','HC50'];leadids=[f'PACE-{i:02}' for i in range(1,21)]
m=np.array([[float(next(x for x in A if x['peptide_id']==p and x['target']==t)['reference_percentile']) for t in T] for p in leadids]);assert m.shape==(20,7)
fig,ax=plt.subplots(figsize=(11,7),facecolor='white');im=ax.imshow(m,interpolation='nearest',aspect='auto',vmin=0,vmax=100,cmap='viridis')
ax.set_xticks(np.arange(len(T)),labels=[x.replace('. ','\n') for x in T]);ax.set_yticks(np.arange(20),labels=leadids);ax.set_title('Original PACE-20: descriptor-space 5NN distance percentiles',fontsize=13,pad=15)
for i in range(20):
 for j in range(7):
  ax.text(j,i,f'{m[i,j]:.0f}',ha='center',va='center',fontsize=7,color='white' if m[i,j]<55 else 'black')
cbar=fig.colorbar(im,ax=ax,shrink=.83);cbar.set_label('Percentile versus self-excluded measured-peptide 5NN distances')
fig.text(.51,.02,'All 140 lead–endpoint checks below the 95th-percentile reference distance cutoff. Not a predictive-accuracy guarantee.',ha='center',fontsize=9)
fig.tight_layout(rect=(0,.05,1,1));fig.savefig(D/'Phase14C_Figure_AD_heatmap.png',dpi=250,bbox_inches='tight');plt.close(fig)
P=read('Phase14C_original_vs_reconstructed_model_predictions.csv')
fig,ax=plt.subplots(figsize=(7.5,7),facecolor='white')
for t in T:
 x=[float(r['original_log10']) for r in P if r['endpoint']==t];y=[float(r['independent_refit_log10']) for r in P if r['endpoint']==t];ax.scatter(x,y,s=24,alpha=.7,label=t)
xx=[min(float(r['original_log10']) for r in P)-.15,max(float(r['original_log10']) for r in P)+.15];ax.plot(xx,xx,ls='--',lw=1.3,color='black',label='Identity')
ax.set_xlabel('Original S12 point predictions (log₁₀ µM)');ax.set_ylabel('Independent refitted model predictions (log₁₀ µM)');ax.set_title('Reconstructed quantitative predictors versus original PACE lead scores');ax.legend(fontsize=8,ncol=2,frameon=False);ax.grid(alpha=.2);fig.tight_layout();fig.savefig(D/'Phase14C_Figure_original_vs_refit.png',dpi=250,bbox_inches='tight');plt.close(fig)
print('Figures created: heatmap and model agreement')
