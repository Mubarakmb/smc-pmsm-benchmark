"""Summarise out-of-sample evaluation: scenario table, parameter-variation degradation,
Monte-Carlo statistics with Friedman test and mean ranks, and two figures."""
from smcbench.paths import RESULTS, FIGS
import json, os
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PROT = os.environ.get('PROTOCOL', 'lyap')
R = json.load(open(f'{RESULTS}/eval_results_{PROT}' + os.environ.get('TAG', '_final') + '.json'))
order = [n for n in ['CSMC', 'ISMC', 'TSMC', 'ASMC', 'FOSMC', 'CSMC-BL', 'ISMC-BL', 'TSMC-BL',
                     'ASMC-BL', 'FOSMC-BL', 'STSMC'] if n in R]
NS = {'CSMC': 6.1, 'ISMC': 8.1, 'TSMC': 33.6, 'ASMC': 16.9, 'FOSMC': 3046, 'STSMC': 6.9}   # mean of 3 timing runs

print(f'=== protocol: {PROT} ===\nScenario cost J (lower is better; >=100 means failed steady-state test)')
sc = list(R[order[0]]['scen'])
print('variant'.ljust(10) + ''.join(s[:13].rjust(14) for s in sc))
for n in order:
    print(n.ljust(10) + ''.join(f"{R[n]['scen'][s]['J']:14.3f}" for s in sc))

print('\nTorque ripple p2p at end of run (N m) per scenario')
for n in order:
    print(n.ljust(10) + ''.join(f"{R[n]['scen'][s].get('Te_p2p', np.nan):14.3f}" for s in sc))

print('\nParameter variation: % change of load-IAE relative to nominal (E0)')
pv = list(R[order[0]]['pvar'])
print('variant'.ljust(10) + ''.join(s[:12].rjust(13) for s in pv))
for n in order:
    b = R[n]['scen']['E0 nominal']['IAE_load']
    print(n.ljust(10) + ''.join(f"{100*(R[n]['pvar'][s].get('IAE_load', np.nan)/b-1):13.1f}" for s in pv))

print('\nMonte Carlo (N=200): median [IQR]')
keys = ['J', 'IAE_load', 'dw_max', 'Te_p2p', 'TV_iqref']
print('variant'.ljust(10) + ''.join(k.rjust(24) for k in keys) + '   fail%')
Jmat = []
for n in order:
    mc = R[n]['mc']
    row = n.ljust(10)
    for k in keys:
        v = np.array([r.get(k, np.nan) for r in mc], float)
        q1, q2, q3 = np.nanpercentile(v, [25, 50, 75])
        row += f"{q2:10.3g} [{q1:5.3g},{q3:5.3g}]".rjust(24)
    J = np.array([r['J'] for r in mc])
    row += f"{100*np.mean(J >= 100):8.1f}"
    print(row)
    Jmat.append(J)
Jmat = np.array(Jmat)
chi2, p = stats.friedmanchisquare(*Jmat)
ranks = np.mean(stats.rankdata(Jmat, axis=0), axis=1)
print(f'\nFriedman test over MC cost: chi2={chi2:.1f}, p={p:.2e}')
print('Mean rank (1 = best): ' + ', '.join(f"{n} {r:.2f}" for n, r in sorted(zip(order, ranks), key=lambda z: z[1])))

fig, ax = plt.subplots(1, 2, figsize=(14, 5))
ax[0].boxplot([np.clip(j, 1e-3, 1e3) for j in Jmat], tick_labels=order, showfliers=False)
ax[0].set_yscale('log'); ax[0].set_ylabel('MC cost J (log)'); ax[0].tick_params(axis='x', rotation=45)
ax[0].set_title(f'Monte Carlo cost, N=200 ({PROT})')
for n in order:
    mc = R[n]['mc']
    x = np.nanmedian([r.get('IAE_load', np.nan) for r in mc])
    y = np.nanmedian([r.get('Te_p2p', np.nan) for r in mc])
    base = n.replace('-BL', '')
    ax[1].scatter(x, y, s=20 + 30 * np.log10(NS[base]), marker='o' if n.endswith('-BL') else 's')
    ax[1].annotate(n, (x, y), fontsize=8)
ax[1].set_xscale('log'); ax[1].set_yscale('log')
ax[1].set_xlabel('median load IAE (rad)'); ax[1].set_ylabel('median torque ripple p2p (N m)')
ax[1].set_title('Accuracy vs smoothness (marker size ~ log exec time)')
plt.tight_layout(); plt.savefig(f'{FIGS}/diag_mc_{PROT}' + os.environ.get('TAG', '_final') + '.png', dpi=80)

# ---------- weight sensitivity of the ranking (R3-8) ----------
if 'TV_both' in R[order[0]]['mc'][0]:
    from smcbench.tune import M
    N3 = 2 * M['I_rated_pk'] / 1e-3
    grid = [(a, b, 1 - a - b) for a in np.linspace(0.05, 0.9, 18) for b in np.linspace(0.05, 0.9, 18) if 1 - a - b >= 0.05]
    top = {n: 0 for n in order}
    rk = {n: [] for n in order}
    for a, b, c in grid:
        Jw = []
        for n in order:
            v = []
            for r in R[n]['mc']:
                if r['J'] >= 100 or 'IAE_trk' not in r:
                    v.append(1e6); continue
                band = 0.005 * abs(r['w_ref'])
                v.append(3 * (a * r['IAE_trk'] / (band * 0.3) + b * r['IAE_load'] / (band * 0.3) + c * r['TV_both'] / N3))
            Jw.append(v)
        mr = np.mean(stats.rankdata(np.array(Jw), axis=0), axis=1)
        for n, x in zip(order, mr):
            rk[n].append(x)
        top[order[int(np.argmin(mr))]] += 1
    print(f'\nWeight sensitivity over {len(grid)} weight triples (w_trk, w_load, w_chatter), each >= 0.05:')
    for n in sorted(order, key=lambda n: np.mean(rk[n])):
        print(f"  {n:9s} mean rank {np.mean(rk[n]):5.2f}  best {min(rk[n]):5.2f}  worst {max(rk[n]):5.2f}  "
              f"ranked #1 in {100*top[n]/len(grid):5.1f} % of weightings")
