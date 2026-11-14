"""Monte Carlo ranking within subgroups of the operating conditions (second revision).

Every Monte Carlo run varies all plant parameters, the speed and the load together; this script groups
the 200 constrained-protocol runs by drawn speed, load and inertia, repeats the Friedman ranking in each
group and writes tables/tab_subgroups.tex. It uses only results/eval_results_lyap_final.json.

    python -m smcbench.subgroups
"""
import json, os
import numpy as np
from scipy import stats
from smcbench.paths import RESULTS, TABLES

ORDER = ['CSMC', 'ISMC', 'TSMC', 'ASMC', 'FOSMC', 'CSMC-BL', 'ISMC-BL', 'TSMC-BL', 'ASMC-BL', 'FOSMC-BL', 'STSMC']
W_RATED = 314.159


def main():
    ev = json.load(open(os.path.join(RESULTS, 'eval_results_lyap_final.json')))
    J = np.array([[r['J'] for r in ev[n]['mc']] for n in ORDER])
    R = stats.rankdata(J, axis=0)
    allr = R.mean(1)
    mc = ev['ISMC-BL']['mc']              # the draws are identical for all variants
    w = np.array([r['w_ref'] for r in mc]) / W_RATED
    L = np.array([r['load'] for r in mc])
    Jx = np.array([r['f_J'] for r in mc])
    groups = [('All runs', np.ones(len(mc), bool)),
              ('Speed $<40\\,\\%$', w < .4), ('Speed $40$--$70\\,\\%$', (w >= .4) & (w < .7)), ('Speed $\\ge70\\,\\%$', w >= .7),
              ('Load $<50\\,\\%$', L < .5), ('Load $50$--$100\\,\\%$', (L >= .5) & (L < 1)), ('Load $\\ge100\\,\\%$', L >= 1),
              ('Inertia $<1\\times$', Jx < 1), ('Inertia $1$--$2\\times$', (Jx >= 1) & (Jx < 2)), ('Inertia $\\ge2\\times$', Jx >= 2)]
    rows = []
    for g, m in groups:
        mr = R[:, m].mean(1)
        o = np.argsort(mr)
        tau = stats.kendalltau(allr, mr)[0]
        rows.append(f"{g} & {m.sum()} & " + ' & '.join(f'{ORDER[i]} ({mr[i]:.2f})' for i in o[:3]) + f" & {tau:.2f} \\\\")
        print(f"{g:24s} n={m.sum():3d}  " + ', '.join(f'{ORDER[i]} {mr[i]:.2f}' for i in o[:3]) + f"  tau={tau:.2f}")
    hi = Jx >= 2
    a, b = J[ORDER.index('FOSMC-BL'), hi], J[ORDER.index('ISMC-BL'), hi]
    print(f"inertia >= 2x: FOSMC-BL better than ISMC-BL in {100*np.mean(a < b):.0f} % of runs, "
          f"Wilcoxon p = {stats.wilcoxon(a, b).pvalue:.2g}")
    tex = r"""\begin{table*}[!t]
\centering
\caption{Monte Carlo ranking within subgroups of the 200 runs, defined by the drawn speed, load and inertia (all other parameters vary simultaneously in every run). Mean rank in parentheses; $\tau$: Kendall rank correlation of the full eleven-variant ranking with the ranking over all runs.}
\label{tab:subgroups}
\small
\begin{tabular}{@{}lccccc@{}}
\toprule
\textbf{Subgroup} & \textbf{Runs} & \textbf{First} & \textbf{Second} & \textbf{Third} & $\tau$ \\
\midrule
""" + '\n'.join(rows) + r"""
\bottomrule
\end{tabular}
\end{table*}
"""
    open(os.path.join(TABLES, 'tab_subgroups.tex'), 'w').write(tex)


if __name__ == '__main__':
    main()
