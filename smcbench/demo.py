"""Quick start: simulate one or more tuned controllers in a chosen test scenario.

    python -m smcbench.demo                      # all 11 variants, nominal scenario E0
    python -m smcbench.demo ISMC-BL,STSMC "E4 reversal"
    PROTOCOL=free python -m smcbench.demo        # gains from the unconstrained protocol
    python -m smcbench.demo ISMC-BL,STSMC "E0 nominal" --show   # also open an interactive window

Scenarios: "E0 nominal", "E1 low speed", "E2 mid speed", "E3 load 25%", "E3 load 50%",
"E3 load 150%", "E4 reversal", "E5 saturation".
Uses the best seed of results/tuning_results_{PROTOCOL}_final.json, prints the main metrics
and writes figures/demo_<scenario>.png (speed error, current command, torque).
"""
import sys
import numpy as np
import matplotlib
SHOW = '--show' in sys.argv
if SHOW:
    sys.argv.remove('--show')
else:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
from smcbench.paths import FIGS
from smcbench import tune as T, evaluate as E

ORDER = ['CSMC', 'ISMC', 'TSMC', 'ASMC', 'FOSMC', 'CSMC-BL', 'ISMC-BL', 'TSMC-BL', 'ASMC-BL', 'FOSMC-BL', 'STSMC']

names = sys.argv[1].split(',') if len(sys.argv) > 1 else ORDER
scen = sys.argv[2] if len(sys.argv) > 2 else 'E0 nominal'
if scen not in E.SCEN:
    sys.exit(f'unknown scenario {scen!r}; choose one of {list(E.SCEN)}')
S = E.SCEN[scen]
best = E.best_params()
keys = ['J', 'IAE_trk', 'IAE_load', 'dw_max', 't_rec', 'ess_end', 'TV_iqref', 'Te_p2p']
print(f'protocol={T.PROTOCOL}  scenario={scen}')
print('variant'.ljust(10) + ''.join(k.rjust(11) for k in keys))
fig, ax = plt.subplots(4, 1, figsize=(8, 9), sharex=True)
ref_drawn = False
for n in names:
    r, L = E.evaluate(n, best[n]['x'], S)
    print(n.ljust(10) + ''.join(f"{r.get(k, np.nan):11.4g}" for k in keys))
    t = L[:, 0]
    if not ref_drawn:
        ax[0].plot(t, L[:, 1], 'k--', lw=0.8, label='reference'); ref_drawn = True
    ax[0].plot(t, L[:, 2], lw=0.7, label=n)
    ax[1].plot(t, L[:, 1] - L[:, 2], lw=0.7)
    ax[2].plot(t, L[:, 4], lw=0.6)
    ax[3].plot(t, L[:, 7], lw=0.6)
ax[0].set_ylabel('speed (rad/s)'); ax[0].legend(fontsize=7, ncol=3)
ax[1].set_ylabel('speed error (rad/s)'); ax[1].set_ylim(-5, 5)
ax[2].set_ylabel('i_q* (A)'); ax[3].set_ylabel('T_e (N m)'); ax[3].set_xlabel('t (s)')
ax[0].set_title(f'{scen} ({T.PROTOCOL} protocol)')
fig.tight_layout()
out = f"{FIGS}/demo_{scen.split()[0]}.png"
fig.savefig(out, dpi=100)
print('figure written to', out)
if SHOW:
    plt.show()
