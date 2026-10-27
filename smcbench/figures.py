"""Publication figures (vector PDF) from the final constrained-protocol results.
Colour = controller family (fixed order, validated palette); line style = sgn (dashed) vs BL (solid)."""
from smcbench.paths import RESULTS, FIGS
import json, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from smcbench import tune as T, evaluate as E
from smcbench.sim import gl_weights

OUT = FIGS
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.family': 'serif', 'font.size': 8, 'axes.labelsize': 8, 'legend.fontsize': 7,
                     'xtick.labelsize': 7, 'ytick.labelsize': 7, 'axes.linewidth': 0.6, 'lines.linewidth': 1.0,
                     'axes.grid': True, 'grid.color': '#dddddd', 'grid.linewidth': 0.4, 'savefig.bbox': 'tight',
                     'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42})
COL = {'CSMC': '#2a78d6', 'ISMC': '#eb6834', 'TSMC': '#1baf7a', 'ASMC': '#eda100', 'FOSMC': '#e87ba4', 'STSMC': '#008300'}
MK = {'CSMC': 'o', 'ISMC': 's', 'TSMC': '^', 'ASMC': 'D', 'FOSMC': 'v', 'STSMC': 'P'}
ORDER = ['CSMC', 'ISMC', 'TSMC', 'ASMC', 'FOSMC', 'CSMC-BL', 'ISMC-BL', 'TSMC-BL', 'ASMC-BL', 'FOSMC-BL', 'STSMC']
base = lambda n: n.replace('-BL', '')
ls = lambda n: '-' if (n.endswith('-BL') or n == 'STSMC') else '--'
W1, W2 = 3.5, 7.16

best = E.best_params()
R = json.load(open(f'{RESULTS}/eval_results_{T.PROTOCOL}{T.TAG}.json'))

# ---- Fig: time responses around the load step (nominal scenario) ----
fig, ax = plt.subplots(3, 2, figsize=(W2, 5.2), sharex=True)
for j, grp in enumerate([[n for n in ORDER if not n.endswith('-BL') and n != 'STSMC'],
                         [n for n in ORDER if n.endswith('-BL')] + ['STSMC']]):
    for n in sorted(grp, key=lambda k: 0 if k == 'STSMC' else 1):
        L = T.run(n, best[n]['x'], S=E.SCEN['E0 nominal'])
        t = L[:, 0]; m = (t >= 0.28) & (t <= 0.40)
        kw = dict(color=COL[base(n)], ls=ls(n), lw=0.6 if n == 'STSMC' else 0.9, label=n,
                  alpha=0.55 if n == 'STSMC' else 1.0, zorder=1 if n == 'STSMC' else 2)
        ax[0, j].plot(t[m] * 1e3, L[m, 1] - L[m, 2], **kw)
        ax[1, j].plot(t[m] * 1e3, L[m, 4], **kw)
        ax[2, j].plot(t[m] * 1e3, L[m, 7], **kw)
    ax[0, j].set_title(['(a) switching $\\mathrm{sgn}(s)$', '(b) boundary layer and STSMC'][j], fontsize=8)
    h_, l_ = ax[0, j].get_legend_handles_labels()
    ax[2, j].legend(h_, l_, ncol=3, frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.32))
    ax[2, j].set_xlabel('time (ms)')
ax[0, 0].set_ylabel('speed error $e$ (rad/s)'); ax[1, 0].set_ylabel('$i_q^\\star$ (A)'); ax[2, 0].set_ylabel('$T_e$ (N·m)')
for a in ax[0]: a.set_ylim(-3, 3.5)
for a in ax.flat: a.set_xlim(280, 400)
fig.subplots_adjust(hspace=0.15)
fig.savefig(f'{OUT}/fig_load_step.pdf'); plt.close(fig)

# ---- Fig: Monte Carlo cost (box) + accuracy vs smoothness ----
NS = {'CSMC': 6.22, 'ISMC': 7.50, 'TSMC': 32.53, 'ASMC': 16.18, 'FOSMC': 3005, 'STSMC': 7.18}   # medians of 3 runs
J = [np.clip([r['J'] for r in R[n]['mc']], 1e-2, 1e3) for n in ORDER]
fig, ax = plt.subplots(1, 2, figsize=(W2, 2.7), gridspec_kw={'width_ratios': [1.35, 1]})
bp = ax[0].boxplot(J, tick_labels=ORDER, showfliers=False, patch_artist=True, widths=0.6,
                   medianprops=dict(color='black', lw=1))
for patch, n in zip(bp['boxes'], ORDER):
    patch.set_facecolor(COL[base(n)]); patch.set_alpha(0.35 if not n.endswith('-BL') and n != 'STSMC' else 0.8)
    patch.set_edgecolor('#444444'); patch.set_linewidth(0.6)
ax[0].set_yscale('log'); ax[0].set_ylabel('Monte Carlo objective $\\mathcal{J}$'); ax[0].tick_params(axis='x', rotation=55)
ax[0].set_title('(a) objective over 200 random plants (lower is better)', fontsize=8)
for n in ORDER:
    mc = R[n]['mc']
    x = np.nanmedian([r.get('IAE_load', np.nan) for r in mc]); y = np.nanmedian([r.get('Te_p2p', np.nan) for r in mc])
    ax[1].scatter(x, y, s=22, marker=MK[base(n)], facecolor=COL[base(n)] if ls(n) == '-' else 'white',
                  edgecolor=COL[base(n)], linewidth=1.1, zorder=3, label=n)
ax[1].set_xscale('log'); ax[1].set_yscale('log')
ax[1].set_xlabel('median load IAE (rad)'); ax[1].set_ylabel('median torque ripple (N·m p–p)')
ax[1].set_title('(b) accuracy vs smoothness', fontsize=8)
ax[1].legend(ncol=2, frameon=False, fontsize=6, loc='center left', bbox_to_anchor=(1.0, 0.5), handletextpad=0.2)
fig.tight_layout(); fig.savefig(f'{OUT}/fig_montecarlo.pdf'); plt.close(fig)

# ---- Fig: fractional memory accuracy and execution time ----
Ts, Tend = 1e-4, 0.6
n = int(Tend / Ts) + 1; t = np.arange(n) * Ts
from math import gamma
fig, ax = plt.subplots(1, 2, figsize=(W2, 2.3))
Lms = [100, 300, 1000, 2000, 4000, 6001]
for a_, mk in [(0.5, 'o'), (0.89, 's')]:
    errs = []
    for Lm in Lms:
        w = gl_weights(-a_, Lm); j = np.arange(min(n, Lm))
        val = Ts ** a_ * np.sum(w[j] * 1.0)
        errs.append(abs(val - Tend ** a_ / gamma(1 + a_)) / (Tend ** a_ / gamma(1 + a_)) * 100)
    ax[0].plot([L * Ts * 1e3 for L in Lms], errs, marker=mk, ms=4, color='#e87ba4' if a_ == 0.5 else '#4a3aa7',
               label=f'$\\alpha={a_}$')
ax[0].set_xscale('log'); ax[0].set_yscale('log')
ax[0].set_xlabel('memory window $L_m T_s$ (ms)'); ax[0].set_ylabel('error of ${}_0D_t^{-\\alpha}1$ at 0.6 s (%)')
ax[0].legend(frameon=False); ax[0].set_title('(a) short-memory truncation error', fontsize=8)
ax[0].annotate('full history', (600, 1e-3), xytext=(150, 1e-2), fontsize=7, arrowprops=dict(arrowstyle='->', lw=0.5))
names = ['CSMC', 'ISMC', 'STSMC', 'ASMC', 'TSMC', 'FOSMC']
ax[1].barh(names, [NS[k] for k in names], color=[COL[k] for k in names], height=0.6)
for i, k in enumerate(names):
    ax[1].text(NS[k] * 1.15, i, f'{NS[k]:.0f} ns' if NS[k] < 100 else f'{NS[k]/1000:.1f} µs', va='center', fontsize=7)
ax[1].set_xscale('log'); ax[1].set_xlim(1, 3e4); ax[1].set_xlabel('execution time per control step (ns)')
ax[1].set_title('(b) measured cost (C, -O2, 2.8 GHz)', fontsize=8); ax[1].grid(axis='y', visible=False)
fig.tight_layout(); fig.savefig(f'{OUT}/fig_cost.pdf'); plt.close(fig)
print('figures written to', OUT)
