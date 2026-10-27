"""Generate LaTeX tables for the manuscript directly from result files (no manual transcription)."""
from smcbench.paths import RESULTS, TABLES
import json, os
import numpy as np
from scipy import stats

W = RESULTS
OUT = TABLES
os.makedirs(OUT, exist_ok=True)
ORDER = ['CSMC', 'ISMC', 'TSMC', 'ASMC', 'FOSMC', 'CSMC-BL', 'ISMC-BL', 'TSMC-BL', 'ASMC-BL', 'FOSMC-BL', 'STSMC']
tun = json.load(open(f'{W}/tuning_results_lyap_final.json'))
ev = json.load(open(f'{W}/eval_results_lyap_final.json'))


def g(x, d=3):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return 'NR'
    ax = abs(x)
    if ax == 0: return '0'
    if ax >= 1000: return f'{x:.3g}'.replace('e+0', 'e').replace('e+', 'e') if ax >= 1e5 else f'{x:.0f}'
    if ax >= 100: return f'{x:.0f}'
    if ax >= 10: return f'{x:.1f}'
    if ax >= 0.01: return f'{x:.{d}f}'.rstrip('0').rstrip('.') if d > 2 else f'{x:.2f}'
    return f'{x:.1e}'


PNAME = {'k': 'k', 'q': 'q', 'lam': '\\lambda', 'nu': '\\nu', 'al': '\\alpha', 'Kp': 'K^+', 'g0': '\\gamma_0',
         'tau': '\\tau', 'alf': '\\bar\\alpha', 'ki': 'k_i', 'phi': '\\phi', 'rho': '\\rho', 'k1': 'k_1', 'k2': 'k_2'}

# ---- Table: tuned parameters and seed spread ----
rows = []
for n in ORDER:
    rs = sorted([r for r in tun if r['name'] == n], key=lambda r: r['J'])
    p = rs[0]['params']
    ps = ', '.join(f'${PNAME.get(k, k)}={g(v)}$' for k, v in p.items() if k != 'rho')
    Js = [r['J'] for r in rs]
    rows.append(f"{n} & {ps} & {Js[0]:.4f} & {np.median(Js):.4f} & {Js[-1]:.4f} \\\\")
open(f'{OUT}/tab_gains.tex', 'w').write(r"""\begin{table*}[!t]
\centering
\caption{Tuned parameters (best of three seeds, constrained protocol, 1920 evaluations per run) and training objective $\mathcal J$ over the three seeds.}
\label{tab:gains}
\small
\begin{tabular}{@{}lp{10.2cm}ccc@{}}
\toprule
\textbf{Variant} & \textbf{Tuned parameters} & \textbf{$\mathcal J$ best} & \textbf{median} & \textbf{worst} \\
\midrule
""" + '\n'.join(rows) + r"""
\bottomrule
\end{tabular}
\end{table*}
""")

# ---- Table: nominal scenario incl. chattering indices ----
rows = []
for n in ORDER:
    m = ev[n]['scen']['E0 nominal']
    rows.append(f"{n} & {g(m['IAE_trk'])} & {g(m['IAE_load'])} & {g(m['dw_max'],2)} & {g(m['t_rec'] if m['t_rec']==m['t_rec'] else None)} & "
                f"{g(m['ess_end'])} & {g(m['TV_iqref'])} & {g(m['HF_iq'])} & {g(m['sw_rate'])} & {g(m['w_p2p'])} & {g(m['Te_p2p'])} \\\\")
open(f'{OUT}/tab_nominal.tex', 'w').write(r"""\begin{table*}[!t]
\centering
\caption{Nominal scenario (rated speed, rated-load step at 0.3~s): tracking, disturbance rejection and chattering indices (definitions in Table~\ref{tab:metrics}).}
\label{tab:nominal}
\small
\begin{tabular}{@{}lcccccccccc@{}}
\toprule
& $\mathrm{IAE}_{\text{trk}}$ & $\mathrm{IAE}_{\text{load}}$ & $\Delta\omega$ & $t_{\text{rec}}$ & $e_{ss}$ & TV($i_q^\star$) & HF $i_q$ & switch rate & $\omega$ ripple & $T_e$ ripple \\
\textbf{Variant} & (rad) & (rad) & (rad/s) & (s) & (rad/s) & (A/s) & (A rms) & (1/s) & (rad/s p--p) & (N·m p--p) \\
\midrule
""" + '\n'.join(rows) + r"""
\bottomrule
\end{tabular}
\end{table*}
""")

# ---- Table: test scenarios (objective) ----
SC = ['E1 low speed', 'E2 mid speed', 'E3 load 25%', 'E3 load 150%', 'E4 reversal', 'E5 saturation']
rows = []
for n in ORDER:
    vals = [ev[n]['scen'][s]['J'] for s in SC]
    rows.append(n + ' & ' + ' & '.join(('F' if v >= 100 else f'{v:.3f}') for v in vals) + ' \\\\')
open(f'{OUT}/tab_scen.tex', 'w').write(r"""\begin{table}[!t]
\centering
\caption{Objective $\mathcal J$ on test scenarios not used for tuning (F: steady-state error above 5\,\%).}
\label{tab:scen}
\resizebox{\columnwidth}{!}{%
\begin{tabular}{@{}lcccccc@{}}
\toprule
\textbf{Variant} & E1 & E2 & E3 (25\%) & E3 (150\%) & E4 & E5 \\
\midrule
""" + '\n'.join(rows) + r"""
\bottomrule
\end{tabular}}
\end{table}
""")

# ---- Table: parameter variation (% change in load IAE) ----
PV = ['Rs x1.39 (hot)', 'Rs x0.7', 'L x0.8', 'L x1.2', 'psi x0.947 (hot)', 'psi x0.9', 'J x0.5', 'J x2', 'J x3', 'B=0.003']
HD = ['$R_s{\\times}1.39$', '$R_s{\\times}0.7$', '$L_s{\\times}0.8$', '$L_s{\\times}1.2$', '$\\psi_f{\\times}0.947$', '$\\psi_f{\\times}0.9$',
      '$J{\\times}0.5$', '$J{\\times}2$', '$J{\\times}3$', '$B{=}0.003$']
rows = []
for n in ORDER:
    b = ev[n]['scen']['E0 nominal']['IAE_load']
    pc = lambda v: '0' if round(v) == 0 else f"{v:+.0f}"
    rows.append(n + ' & ' + ' & '.join(pc(100*(ev[n]['pvar'][k]['IAE_load']/b-1)) for k in PV) + ' \\\\')
open(f'{OUT}/tab_pvar.tex', 'w').write(r"""\begin{table}[!t]
\centering
\caption{Plant parameter variations (controllers keep the nominal model): change of $\mathrm{IAE}_{\text{load}}$ relative to the nominal plant (\%).}
\label{tab:pvar}
\resizebox{\columnwidth}{!}{%
\begin{tabular}{@{}lrrrrrrrrrr@{}}
\toprule
\textbf{Variant} & """ + ' & '.join(HD) + r""" \\
\midrule
""" + '\n'.join(rows) + r"""
\bottomrule
\end{tabular}}
\end{table}
""")

# ---- Table: Monte Carlo summary ----
def mranks(E):
    J = np.array([[r['J'] for r in E[n]['mc']] for n in ORDER])
    return dict(zip(ORDER, np.mean(stats.rankdata(J, axis=0), axis=1))), stats.friedmanchisquare(*J)
rk, fr = mranks(ev)
rows = []
for n in sorted(ORDER, key=lambda k: rk[k]):
    mc = ev[n]['mc']
    med = lambda k: np.nanmedian([r.get(k, np.nan) for r in mc])
    q = lambda k: np.nanpercentile([r.get(k, np.nan) for r in mc], [25, 75])
    J = [r['J'] for r in mc]
    rows.append(f"{n} & {rk[n]:.2f} & {np.median(J):.3f} [{np.percentile(J,25):.3f}, {np.percentile(J,75):.3f}] & "
                f"{g(med('IAE_load'))} & {g(med('dw_max'),2)} & {g(med('Te_p2p'))} & {100*np.mean(np.array(J)>=100):.1f} \\\\")
open(f'{OUT}/tab_mc.tex', 'w').write(r"""\begin{table*}[!t]
\centering
\caption{Monte Carlo study ($N=200$ random plants and operating points, identical for all variants): Friedman mean rank of $\mathcal J$ (1 = best), median [interquartile range] of $\mathcal J$, medians of key metrics, and failure rate. Results for the unconstrained tuning protocol are given in Table~\ref{tab:protocol}. Friedman test (constrained): $\chi^2=""" + f"{fr.statistic:.0f}" + r"""$, $p<10^{-300}$.}
\label{tab:mc}
\small
\begin{tabular}{@{}lcccccc@{}}
\toprule
\textbf{Variant} & \textbf{mean rank} & $\mathcal J$ median [IQR] & $\mathrm{IAE}_{\text{load}}$ (rad) & $\Delta\omega$ (rad/s) & $T_e$ ripple (N·m) & fail (\%) \\
\midrule
""" + '\n'.join(rows) + r"""
\bottomrule
\end{tabular}
\end{table*}
""")
print('tables written; Friedman chi2 =', fr.statistic, 'p =', fr.pvalue, )
