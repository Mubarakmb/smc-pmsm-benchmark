"""Constrained ('lyap') versus unconstrained ('free') tuning protocol: tuned gains, paired
Monte Carlo comparison (same 200 plants/operating points, same noise seeds), figure and summary."""
from smcbench.paths import RESULTS, FIGS, TABLES
import json
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker
from smcbench import tune as T

W = RESULTS
OUT = FIGS
plt.rcParams.update({'font.family': 'serif', 'font.size': 8, 'axes.labelsize': 8, 'legend.fontsize': 7,
                     'xtick.labelsize': 7, 'ytick.labelsize': 7, 'axes.linewidth': 0.6, 'axes.grid': True,
                     'grid.color': '#dddddd', 'grid.linewidth': 0.4, 'savefig.bbox': 'tight',
                     'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42})
ORDER = ['CSMC', 'ISMC', 'TSMC', 'ASMC', 'FOSMC', 'CSMC-BL', 'ISMC-BL', 'TSMC-BL', 'ASMC-BL', 'FOSMC-BL', 'STSMC']
D = T.D          # rated-load disturbance bound (4834.4 rad/s^2)
C = 1518783.203391092   # load-rise bound for t_r = 5 ms (T.C_DOT under the default LOAD_RISE)


def best(prot):
    res = json.load(open(f'{W}/tuning_results_{prot}_final.json'))
    b = {}
    for r in res:
        if r['name'] not in b or r['J'] < b[r['name']]['J']:
            b[r['name']] = r
    return b, res


def gain_ratio(n, p):
    """Switching gain normalised by the bound the reaching condition requires."""
    if n == 'STSMC':
        return p['k2'] / C
    if n.startswith('ASMC'):
        return p['Kp'] / D
    return p['k'] / D


B = {pr: best(pr) for pr in ('lyap', 'free')}
EV = {pr: json.load(open(f'{W}/eval_results_{pr}_final.json')) for pr in ('lyap', 'free')}
Jm = {pr: np.array([[r['J'] for r in EV[pr][n]['mc']] for n in ORDER]) for pr in EV}
rank = {pr: np.mean(stats.rankdata(Jm[pr], axis=0), axis=1) for pr in EV}
joint = stats.rankdata(np.vstack([Jm['lyap'], Jm['free']]), axis=0).mean(axis=1)
med = lambda pr, k: np.array([np.nanmedian([r.get(k, np.nan) for r in EV[pr][n]['mc']]) for n in ORDER])

lines = ['variant  | J_train lyap  free | gain/bound lyap  free | MC rank lyap free | joint rank lyap free | '
         'med Te_p2p lyap free | med IAE_load lyap free | med dw lyap free | fail% lyap free | med J ratio free/lyap, Wilcoxon p']
summary = {}
for i, n in enumerate(ORDER):
    bl, bf = B['lyap'][0][n], B['free'][0][n]
    jl, jf = Jm['lyap'][i], Jm['free'][i]
    ok = (jl < 100) & (jf < 100)
    ratio = np.median(jf[ok] / jl[ok])
    try:
        pw = stats.wilcoxon(jf[ok], jl[ok]).pvalue
    except ValueError:
        pw = np.nan
    seedsJ = {pr: sorted(r['J'] for r in B[pr][1] if r['name'] == n) for pr in B}
    s = dict(Jtr=(bl['J'], bf['J']), seedsJ=seedsJ, g=(gain_ratio(n, bl['params']), gain_ratio(n, bf['params'])),
             params=(bl['params'], bf['params']), rank=(rank['lyap'][i], rank['free'][i]),
             joint=(joint[i], joint[i + len(ORDER)]),
             Te=(med('lyap', 'Te_p2p')[i], med('free', 'Te_p2p')[i]),
             IAE=(med('lyap', 'IAE_load')[i], med('free', 'IAE_load')[i]),
             dw=(med('lyap', 'dw_max')[i], med('free', 'dw_max')[i]),
             TV=(med('lyap', 'TV_iqref')[i], med('free', 'TV_iqref')[i]),
             fail=(100 * np.mean(jl >= 100), 100 * np.mean(jf >= 100)), ratio=ratio, p=pw,
             frac_better=float(np.mean(jf[ok] < jl[ok])))
    summary[n] = s
    lines.append(f"{n:9s}| {s['Jtr'][0]:.4f} {s['Jtr'][1]:.4f} | {s['g'][0]:7.3f} {s['g'][1]:7.3f} | "
                 f"{s['rank'][0]:5.2f} {s['rank'][1]:5.2f} | {s['joint'][0]:6.2f} {s['joint'][1]:6.2f} | "
                 f"{s['Te'][0]:7.3f} {s['Te'][1]:7.3f} | {s['IAE'][0]:.4f} {s['IAE'][1]:.4f} | "
                 f"{s['dw'][0]:.3f} {s['dw'][1]:.3f} | {s['fail'][0]:4.1f} {s['fail'][1]:4.1f} | "
                 f"{ratio:.3f} p={pw:.1e} better={s['frac_better']:.2f}")
fr = stats.friedmanchisquare(*Jm['free'])
lines.append(f'Friedman (free): chi2={fr.statistic:.1f} p={fr.pvalue:.2e}')
for n in ORDER:
    lines.append(f"{n}: lyap {B['lyap'][0][n]['params']}\n      free {B['free'][0][n]['params']}  seedsJ free {summary[n]['seedsJ']['free']}")
txt = '\n'.join(lines)
open(f'{W}/protocol_compare.txt', 'w').write(txt)
json.dump(summary, open(f'{W}/protocol_compare.json', 'w'), indent=1, default=float)
print(txt)

# ---------------- figure ----------------
x = np.arange(len(ORDER)); w = 0.38
cL, cF = '#2a78d6', '#eb6834'
fig, ax = plt.subplots(1, 3, figsize=(7.16, 2.5))
g = np.array([summary[n]['g'] for n in ORDER])
ax[0].bar(x - w / 2, g[:, 0], w, color=cL, label='constrained')
ax[0].bar(x + w / 2, g[:, 1], w, color=cF, label='unconstrained')
ax[0].axhline(1.0, color='k', lw=0.7, ls='--'); ax[0].axhline(0.1, color='0.4', lw=0.6, ls=':')
ax[0].set_yscale('log'); ax[0].set_ylabel(r'switching gain / required bound')
ax[0].set_title('(a) tuned switching gain', fontsize=8)
rt = np.array([summary[n]['ratio'] for n in ORDER])
ax[1].bar(x, rt, 0.6, color=['#1baf7a' if v < 1 else '#c0392b' for v in rt])
ax[1].axhline(1.0, color='k', lw=0.7, ls='--'); ax[1].set_yscale('log'); ax[1].set_ylim(0.1, 2.5); ax[1].set_yticks([0.1, 0.2, 0.5, 1, 2]); ax[1].set_yticklabels(['0.1', '0.2', '0.5', '1', '2']); ax[1].yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax[1].set_ylabel(r'median of $\mathcal{J}_{\mathrm{unc}}/\mathcal{J}_{\mathrm{con}}$')
ax[1].set_title('(b) paired Monte Carlo change', fontsize=8)
te = np.array([summary[n]['Te'] for n in ORDER])
ax[2].bar(x - w / 2, te[:, 0], w, color=cL); ax[2].bar(x + w / 2, te[:, 1], w, color=cF)
ax[2].set_yscale('log'); ax[2].set_ylim(0.03, 40); ax[2].set_ylabel(r'median $T_e$ ripple (N$\cdot$m)')
ax[2].set_title('(c) Monte Carlo torque ripple', fontsize=8)
for a in ax:
    a.set_xticks(x); a.set_xticklabels(ORDER, rotation=60, ha='right')
ax[0].set_ylim(0.03, 30); ax[0].legend(loc='upper left', frameon=False, ncol=2, fontsize=6.5)
fig.tight_layout(); fig.savefig(f'{OUT}/fig_protocol.pdf'); plt.close(fig)

# ---------------- table ----------------
def f2(v): return f'{v:.2f}'
def fr(v): return f'{v:.3g}' if v < 10 else f'{v:.1f}'
rows = []
for n in ORDER:
    s = summary[n]
    gl, gu = s['g']
    cond = 'yes' if gu >= 1.0 else 'no'
    if n == 'STSMC':
        k1, k2 = s['params'][1]['k1'], s['params'][1]['k2']
        cond = 'no' if (k2 < C or k1 < np.sqrt(k2 + C)) else 'yes'
    rows.append(f"{n} & {gl:.2f} & {gu:.2f} & {cond} & {s['Jtr'][0]:.3f} & {s['Jtr'][1]:.3f} & "
                f"{f2(s['rank'][0])} & {f2(s['rank'][1])} & {fr(s['Te'][0])} & {fr(s['Te'][1])} & "
                f"{s['ratio']:.2f} & {100*s['frac_better']:.0f} & {s['fail'][1]:.1f} \\\\")
tex = r"""\begin{table*}[!t]
\centering
\caption{Stability-constrained (con) versus unconstrained (unc) tuning. Switching gain relative to the bound required by the reaching condition ($k/D$, $K^+/D$, or $k_2/C$ for STSMC); whether the unconstrained gains still satisfy the condition of Table~\ref{tab:conditions}; training objective; Monte Carlo mean rank within each protocol; median torque ripple; median of the paired ratio $\mathcal J_{\text{unc}}/\mathcal J_{\text{con}}$ over the same 200 Monte Carlo runs (runs failed under either protocol excluded); share of runs in which the unconstrained gains are better; failure rate with unconstrained gains.}
\label{tab:protocol}
\small
\setlength{\tabcolsep}{4pt}
\begin{tabular}{@{}lcccccccccccc@{}}
\toprule
 & \multicolumn{3}{c}{gain / bound} & \multicolumn{2}{c}{training $\mathcal J$} & \multicolumn{2}{c}{MC rank} & \multicolumn{2}{c}{$T_e$ ripple (N·m)} & median & unc better & unc fail \\
\cmidrule(lr){2-4}\cmidrule(lr){5-6}\cmidrule(lr){7-8}\cmidrule(lr){9-10}
\textbf{Variant} & con & unc & cond.\ met & con & unc & con & unc & con & unc & $\mathcal J_{\text{unc}}/\mathcal J_{\text{con}}$ & (\%) & (\%) \\
\midrule
""" + '\n'.join(rows) + r"""
\bottomrule
\end{tabular}
\end{table*}
"""
open(f'{TABLES}/tab_protocol.tex', 'w').write(tex)
