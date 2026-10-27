"""Additional publication figures from the final (constrained, 2 us) results."""
from smcbench.paths import RESULTS, FIGS
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, SymLogNorm
from smcbench import tune as T, evaluate as E

W = RESULTS
OUT = FIGS
plt.rcParams.update({'font.family': 'serif', 'font.size': 8, 'axes.labelsize': 8, 'legend.fontsize': 7,
                     'xtick.labelsize': 7, 'ytick.labelsize': 7, 'axes.linewidth': 0.6, 'axes.grid': True,
                     'grid.color': '#dddddd', 'grid.linewidth': 0.4, 'savefig.bbox': 'tight',
                     'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42})
COL = {'CSMC': '#2a78d6', 'ISMC': '#eb6834', 'TSMC': '#1baf7a', 'ASMC': '#eda100', 'FOSMC': '#e87ba4', 'STSMC': '#008300'}
MK = {'CSMC': 'o', 'ISMC': 's', 'TSMC': '^', 'ASMC': 'D', 'FOSMC': 'v', 'STSMC': 'P'}
ORDER = ['CSMC', 'ISMC', 'TSMC', 'ASMC', 'FOSMC', 'CSMC-BL', 'ISMC-BL', 'TSMC-BL', 'ASMC-BL', 'FOSMC-BL', 'STSMC']
CONT = ['CSMC-BL', 'ISMC-BL', 'TSMC-BL', 'ASMC-BL', 'FOSMC-BL', 'STSMC']
base = lambda n: n.replace('-BL', '')
ls = lambda n: '-' if (n.endswith('-BL') or n == 'STSMC') else '--'
W2 = 7.16
best = E.best_params()
ev = json.load(open(f'{W}/eval_results_lyap_final.json'))

# ---- test scenarios: low speed, reversal, saturation (continuous variants) ----
sc = [('E1 low speed', '(a) E1: 10% speed, 50% load', (0.25, 0.45)),
      ('E4 reversal', '(b) E4: speed reversal', (0.25, 0.45)),
      ('E5 saturation', '(c) E5: near-step reference', (0.0, 0.12))]
fig, ax = plt.subplots(2, 3, figsize=(W2, 3.7))
DRAW = ['STSMC', 'CSMC-BL', 'ASMC-BL', 'FOSMC-BL', 'TSMC-BL', 'ISMC-BL']   # best-ranked drawn last (on top)
for j, (k, title, (a, b)) in enumerate(sc):
    for n in DRAW:
        L = T.run(n, best[n]['x'], S=E.SCEN[k]); t = L[:, 0]; m = (t >= a) & (t <= b)
        kw = dict(color=COL[base(n)], lw=0.55 if n == 'STSMC' else 0.9, alpha=0.55 if n == 'STSMC' else 1, label=n)
        if k == 'E5 saturation':
            ax[0, j].plot(t[m] * 1e3, L[m, 2], **kw)
        else:
            ax[0, j].plot(t[m] * 1e3, L[m, 1] - L[m, 2], **kw)
        ax[1, j].plot(t[m] * 1e3, L[m, 4], **kw)
    for r in range(2):
        ax[r, j].set_xlim(a * 1e3, b * 1e3)
    ax[0, j].set_title(title, fontsize=8); ax[1, j].set_xlabel('time (ms)')
ax[0, 0].set_ylabel('speed error (rad/s)'); ax[0, 1].set_ylabel('speed error (rad/s)'); ax[0, 2].set_ylabel('speed (rad/s)')
ax[1, 0].set_ylabel('$i_q^\\star$ (A)')
ax[0, 0].set_ylim(-1.5, 1.5); ax[0, 1].set_ylim(-3, 3)
h_, l_ = ax[0, 0].get_legend_handles_labels()
order = [l_.index(n) for n in CONT]
fig.legend([h_[i] for i in order], [l_[i] for i in order], ncol=6, frameon=False, loc='lower center', bbox_to_anchor=(0.5, 0.0))
fig.subplots_adjust(left=0.07, right=0.99, top=0.93, bottom=0.2, wspace=0.3, hspace=0.32)
fig.savefig(f'{OUT}/fig_scenarios.pdf'); plt.close(fig)

# ---- parameter-variation heatmap ----
PV = ['Rs x1.39 (hot)', 'Rs x0.7', 'L x0.8', 'L x1.2', 'psi x0.947 (hot)', 'psi x0.9', 'J x0.5', 'J x2', 'J x3', 'B=0.003']
HD = ['$R_s$×1.39', '$R_s$×0.7', '$L_s$×0.8', '$L_s$×1.2', '$\\psi_f$×0.947', '$\\psi_f$×0.9', '$J$×0.5', '$J$×2', '$J$×3', '$B$=0.003']
M = np.array([[100 * (ev[n]['pvar'][k]['IAE_load'] / ev[n]['scen']['E0 nominal']['IAE_load'] - 1) for k in PV] for n in ORDER])
cmap = LinearSegmentedColormap.from_list('div', ['#2a78d6', '#f0efec', '#e34948'])
fig, axh = plt.subplots(figsize=(W2 * 0.62, 3.0))
im = axh.imshow(M, cmap=cmap, norm=SymLogNorm(linthresh=10, vmin=-3000, vmax=3000), aspect='auto')
axh.set_xticks(range(len(PV)), HD, rotation=45, ha='right'); axh.set_yticks(range(len(ORDER)), ORDER)
axh.grid(False)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        axh.text(j, i, f'{M[i, j]:+.0f}', ha='center', va='center', fontsize=5.5, color='#111111')
cb = fig.colorbar(im, ax=axh, fraction=0.04, pad=0.02); cb.set_label('change of load IAE (%)', fontsize=7)
fig.savefig(f'{OUT}/fig_pvar.pdf'); plt.close(fig)

# ---- sensitivity: sampling time and load-rise time ----
sens = json.load(open(f'{W}/sensitivity_lyap.json'))
fig, ax = plt.subplots(1, 2, figsize=(W2, 2.4))
for n in ORDER:
    Ts = [r[0] * 1e6 for r in sens[n]['Ts']]; J = [min(r[1], 400) for r in sens[n]['Ts']]
    ax[0].plot(Ts, J, ls=ls(n), marker=MK[base(n)], ms=3.5, color=COL[base(n)], lw=0.9, label=n,
               mfc=COL[base(n)] if ls(n) == '-' else 'white')
ax[0].set_yscale('log'); ax[0].set_xticks([50, 100, 200]); ax[0].set_xlabel('sampling time $T_s$ (µs)')
ax[0].set_ylabel('objective $\\mathcal{J}$ (gains tuned at 100 µs)'); ax[0].set_title('(a) sampling-time sensitivity', fontsize=8)
ax[0].legend(ncol=2, frameon=False, fontsize=6, loc='center left', bbox_to_anchor=(1.0, 0.5))
tr = {5: json.load(open(f'{W}/tuning_results_lyap_final.json')), 20: json.load(open(f'{W}/tuning_results_lyap_tr20.json')),
      50: json.load(open(f'{W}/tuning_results_lyap_tr50.json'))}
for n in ['ISMC-BL', 'TSMC-BL', 'STSMC']:
    y = [min(r['J'] for r in tr[k] if r['name'] == n) for k in [5, 20, 50]]
    ax[1].plot([5, 20, 50], y, marker=MK[base(n)], color=COL[base(n)], lw=1, label=n)
ax[1].set_xlabel('load-rise time $t_r$ (ms)'); ax[1].set_ylabel('tuned objective $\\mathcal{J}$')
ax[1].set_yscale('log'); ax[1].set_xticks([5, 20, 50]); ax[1].legend(frameon=False)
ax[1].set_title('(b) load-rise sensitivity (re-tuned)', fontsize=8)
fig.tight_layout(); fig.savefig(f'{OUT}/fig_sensitivity.pdf'); plt.close(fig)

# ---- current spectra in steady state (nominal, loaded) ----
fig, axs = plt.subplots(figsize=(W2 * 0.5, 2.6))
for n in ['CSMC', 'CSMC-BL', 'ISMC-BL', 'STSMC']:
    L = T.run(n, best[n]['x'], S=E.SCEN['E0 nominal']); t = L[:, 0]; m = t >= 0.5
    x = L[m, 5] - L[m, 5].mean(); dt = t[1] - t[0]
    X = np.abs(np.fft.rfft(x * np.hanning(len(x)))) * 2 / len(x); f = np.fft.rfftfreq(len(x), dt)
    axs.semilogy(f[1:], X[1:] + 1e-6, color=COL[base(n)], ls=ls(n), lw=0.8, label=n)
axs.set_xlim(0, 5000); axs.set_ylim(1e-5, 30); axs.set_xlabel('frequency (Hz)'); axs.set_ylabel('$|I_q(f)|$ (A)')
axs.legend(frameon=False, ncol=2, fontsize=6, loc='upper right'); axs.set_title('(b) steady-state $i_q$ spectrum (rated load)', fontsize=8)
fig.subplots_adjust(left=0.17, right=0.97, bottom=0.17, top=0.9)
with plt.rc_context({'savefig.bbox': 'standard'}): fig.savefig(f'{OUT}/fig_spectrum.pdf'); plt.close(fig)   # same canvas as fig_sliding

# ---- sliding-variable trajectories around the load step ----
fig, axp = plt.subplots(figsize=(W2 * 0.5, 2.6))
for n in CONT:
    L = T.run(n, best[n]['x'], S=E.SCEN['E0 nominal']); t = L[:, 0]; m = (t >= 0.29) & (t <= 0.36)
    axp.plot(t[m] * 1e3, L[m, 9], color=COL[base(n)], lw=0.55 if n == 'STSMC' else 0.9,
             alpha=0.55 if n == 'STSMC' else 1, label=n)
axp.set_xlim(290, 360); axp.set_ylim(-1.2, 7.4); axp.set_xlabel('time (ms)'); axp.set_ylabel('sliding variable $s$ (rad/s)')
axp.set_title('(a) $s(t)$ around the rated-load step', fontsize=8)
axp.legend(frameon=False, ncol=3, fontsize=6, loc='upper center')
fig.subplots_adjust(left=0.17, right=0.97, bottom=0.17, top=0.9)
with plt.rc_context({'savefig.bbox': 'standard'}): fig.savefig(f'{OUT}/fig_sliding.pdf'); plt.close(fig)   # same canvas as fig_spectrum

# ---- non-constant disturbances ----
d = json.load(open(f'{W}/disturb_lyap.json'))
names = ['CSMC-BL', 'ISMC-BL', 'TSMC-BL', 'FOSMC-BL', 'ASMC-BL', 'STSMC']
fig, ax = plt.subplots(1, 3, figsize=(W2, 2.3))
for j, w in enumerate([0.02, 0.1]):
    for n in names:
        y = [d[f'ripple w={w} a={a}'][n][0] for a in [0.01, 0.03, 0.05]]
        ax[j].plot([1, 3, 5], y, marker=MK[base(n)], ms=3.5, color=COL[base(n)], lw=0.9, label=n)
    ax[j].set_yscale('log'); ax[j].set_xlabel('ripple amplitude (% of $T_{\\mathrm{rated}}$)')
    ax[j].set_title(f'({"ab"[j]}) {int(w*100)}% speed, {int(6*4*w*T.M["w_rated"]/2/np.pi)} Hz ripple', fontsize=8)
ax[0].set_ylabel('speed ripple (rad/s p–p)')
for n in names:
    y = [d[f'coulomb {tf}'][n] for tf in [0.0, 0.02, 0.05]]
    ax[2].plot([0, 2, 5], y, marker=MK[base(n)], ms=3.5, color=COL[base(n)], lw=0.9, label=n)
ax[2].set_yscale('log'); ax[2].set_xlabel('Coulomb friction (% of $T_{\\mathrm{rated}}$)'); ax[2].set_ylabel('reversal IAE (rad)')
ax[2].set_title('(c) Coulomb friction, ±5% reversal', fontsize=8)
ax[2].legend(frameon=False, fontsize=6, loc='center left', bbox_to_anchor=(1.0, 0.5))
fig.tight_layout(); fig.savefig(f'{OUT}/fig_disturb.pdf'); plt.close(fig)
print('done')
