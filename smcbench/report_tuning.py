"""Summarise tuning results: best/median/worst cost per variant over seeds, tuned
parameters, nominal metrics, and response figures for the best seed."""
from smcbench.paths import RESULTS, FIGS
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from smcbench import tune as T

res = json.load(open(f'{RESULTS}/tuning_results_{T.PROTOCOL}{T.TAG}.json'))
by = {}
for r in res:
    by.setdefault(r['name'], []).append(r)
order = ['CSMC', 'ISMC', 'TSMC', 'ASMC', 'FOSMC', 'CSMC-BL', 'ISMC-BL', 'TSMC-BL', 'ASMC-BL', 'FOSMC-BL', 'STSMC']
order = [n for n in order if n in by]
print(f"{'variant':9s} {'J_best':>8s} {'J_med':>8s} {'J_worst':>8s} {'nfev':>5s}  params(best)")
best = {}
for n in order:
    rs = sorted(by[n], key=lambda r: r['J'])
    best[n] = rs[0]
    Js = [r['J'] for r in rs]
    ps = ', '.join(f"{k}={v:.4g}" for k, v in rs[0]['params'].items())
    print(f"{n:9s} {Js[0]:8.4f} {np.median(Js):8.4f} {Js[-1]:8.4f} {rs[0]['nfev']:5d}  {ps}")
keys = ['IAE_trk', 'IAE_load', 'dw_max', 't_rec', 'ess_end', 'TV_iqref', 'HF_iq', 'w_p2p', 'Te_p2p', 'sat_%']
print('\n' + 'variant'.ljust(10) + ''.join(k.rjust(10) for k in keys))
for n in order:
    m = best[n]['metrics']
    print(n.ljust(10) + ''.join(f"{m[k]:10.4g}" for k in keys))

fig, ax = plt.subplots(3, 2, figsize=(13, 10))
groups = [[n for n in order if not n.endswith('-BL') and n != 'STSMC'] + ['STSMC'],
          [n for n in order if n.endswith('-BL')] + ['STSMC']]
for j, g in enumerate(groups):
    for n in g:
        L = T.run(n, best[n]['x'])
        t = L[:, 0]
        ax[0, j].plot(t, L[:, 1] - L[:, 2], lw=0.7, label=n)
        ax[1, j].plot(t, L[:, 4], lw=0.6)
        ax[2, j].plot(t, L[:, 7], lw=0.6)
    ax[0, j].set_title(['switching sgn(s) + STSMC', 'boundary layer + STSMC'][j])
    ax[0, j].set_ylim(-3, 3); ax[0, j].set_ylabel('e (rad/s)'); ax[0, j].legend(fontsize=7)
    ax[1, j].set_ylabel('iq* (A)'); ax[2, j].set_ylabel('Te (N m)'); ax[2, j].set_xlabel('t (s)')
    for a in ax[:, j]:
        a.axvline(0.3, color='k', lw=0.4, ls=':')
plt.tight_layout()
plt.savefig(f'{FIGS}/diag_tuned_nominal_{T.PROTOCOL}{T.TAG}.png', dpi=80)
