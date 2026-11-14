"""Dedicated anti-windup for the FOSMC fractional memory (second revision, Reviewer comment on
saturation and anti-windup).

With FO_AW=1, the error history of FOSMC pauses (no new samples are stored) while i_q* is saturated
and the error drives further into saturation; this is the counterpart of the conditional integration
used for the other integral states. The script re-evaluates FOSMC and FOSMC-BL with the same tuned
gains on all test scenarios, parameter variations and the 200 Monte Carlo runs, writes
results/eval_results_lyap_foaw.json and prints the comparison with the common rule, including the
Monte Carlo ranking with the anti-windup variants.

    python -m smcbench.foaw_study          # ~2 min
"""
import os
os.environ['FO_AW'] = '1'
os.environ['PROTOCOL'] = 'lyap'      # the study uses the constrained-protocol gains
import json
import numpy as np
from scipy import stats
from smcbench.paths import RESULTS
from smcbench import tune as T, evaluate as E

ORDER = ['CSMC', 'ISMC', 'TSMC', 'ASMC', 'FOSMC', 'CSMC-BL', 'ISMC-BL', 'TSMC-BL', 'ASMC-BL', 'FOSMC-BL', 'STSMC']
NAMES = ['FOSMC', 'FOSMC-BL']


def run():
    assert T.FO_AW == 1
    best = E.best_params()
    res = {}
    for n in NAMES:
        x = best[n]['x']
        res[n] = {'scen': {k: E.evaluate(n, x, S)[0] for k, S in E.SCEN.items()},
                  'pvar': {k: E.evaluate(n, x, E.SCEN['E0 nominal'], Pp=E.plant(**f))[0] for k, f in E.PVAR.items()},
                  'mc': E.monte_carlo(n, x)}
    json.dump(res, open(os.path.join(RESULTS, 'eval_results_lyap_foaw.json'), 'w'), indent=1, default=float)
    return res


def report(res):
    base = json.load(open(os.path.join(RESULTS, 'eval_results_lyap_final.json')))
    def ranks(ev):
        J = np.array([[r['J'] for r in ev[n]['mc']] for n in ORDER])
        return dict(zip(ORDER, stats.rankdata(J, axis=0).mean(1)))
    aw = dict(base); aw.update(res)
    r0, r1 = ranks(base), ranks(aw)
    for n in NAMES:
        a, b = base[n]['scen']['E5 saturation'], res[n]['scen']['E5 saturation']
        nsat = sum(1 for m in base[n]['mc'] if m.get('sat_%', 0) > 0)
        print(f"{n}: E5 J {a['J']:.2f} -> {b['J']:.2f}; MC mean rank {r0[n]:.2f} -> {r1[n]:.2f}; "
              f"MC runs with saturation: {nsat} of {len(base[n]['mc'])}")
    print('Monte Carlo ranking with the anti-windup variants:')
    for n in sorted(ORDER, key=lambda k: r1[k]):
        print(f'  {n:9s} {r1[n]:.2f}')


if __name__ == '__main__':
    report(run())
