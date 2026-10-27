from smcbench.paths import RESULTS, FIGS
import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from scipy import stats
plt.rcParams.update({'font.family': 'serif', 'font.size': 9, 'pdf.fonttype': 42})
COL = {'CSMC': '#2a78d6', 'ISMC': '#eb6834', 'TSMC': '#1baf7a', 'ASMC': '#eda100', 'FOSMC': '#e87ba4', 'STSMC': '#008300'}
ORDER = ['CSMC','ISMC','TSMC','ASMC','FOSMC','CSMC-BL','ISMC-BL','TSMC-BL','ASMC-BL','FOSMC-BL','STSMC']
ev = json.load(open(f'{RESULTS}/eval_results_lyap_final.json'))
J = np.array([[r['J'] for r in ev[n]['mc']] for n in ORDER]); rk = np.mean(stats.rankdata(J, axis=0), axis=1)
fig = plt.figure(figsize=(10, 4))
ax = fig.add_axes([0.01, 0, 0.58, 1]); ax.set_xlim(0, 10); ax.set_ylim(0, 4); ax.axis('off')
boxes = [(0.2, 'Plant & drive', 'Siemens 1FK7063 SPMSM\n(data sheet), PI current\nloops, $V_{dc}/\\sqrt{3}$ and\n$I_{\\max}$ limits, delay,\n20-bit encoder, $T_s$=100 µs'),
         (3.5, 'Six SMC laws', 'CSMC · ISMC · TSMC\nASMC · FOSMC · STSMC\nsgn and boundary layer\nno load in control law\nLyapunov gain bounds'),
         (6.8, 'Fair evaluation', 'identical DE tuning\n(1920 evals, 3 seeds)\n7 unseen scenarios\n200-run Monte Carlo\nFriedman test')]
for x, title, body in boxes:
    ax.add_patch(FancyBboxPatch((x, 0.9), 3.0, 2.3, boxstyle='round,pad=0.05,rounding_size=0.15', fc='#eef3f8', ec='#1F4E79', lw=1))
    ax.text(x + 1.5, 2.8, title, ha='center', va='center', weight='bold', color='#1F4E79', fontsize=10)
    ax.text(x + 1.5, 1.95, body, ha='center', va='center', fontsize=9.2, linespacing=1.3)
for x in [3.22, 6.52]:
    ax.annotate('', (x + 0.26, 1.95), (x, 1.95), arrowprops=dict(arrowstyle='-|>', color='#1F4E79', lw=1.2))
b = fig.add_axes([0.70, 0.17, 0.28, 0.71])
o = np.argsort(rk)
names = [ORDER[i] for i in o][::-1]; vals = [rk[i] for i in o][::-1]
b.barh(names, vals, color=[COL[n.replace('-BL', '')] for n in names], height=0.65,
       alpha=1.0)
for i, n in enumerate(names):
    if not (n.endswith('-BL') or n == 'STSMC'):
        b.patches[i].set_alpha(0.4)
b.set_xlabel('Monte Carlo mean rank (1 = best)\nlight bars: discontinuous versions', fontsize=8); b.tick_params(labelsize=7.5)
b.set_title('Integral SMC + boundary layer\nranks first', fontsize=9, weight='bold', color='#1F4E79')
for s in ['top', 'right']: b.spines[s].set_visible(False)
b.set_xlim(0, 11)
with plt.rc_context({'savefig.bbox': 'standard'}):
    fig.savefig(f'{FIGS}/graphical_abstract.pdf'); fig.savefig(f'{FIGS}/graphical_abstract.png', dpi=300)
