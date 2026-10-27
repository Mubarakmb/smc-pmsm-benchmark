"""Evaluation on scenarios NOT used for tuning. Uses the best seed per variant.

Test set
  E0  nominal      : rated speed, rated load step (training conditions, reported for reference)
  E1  low speed    : 10 % rated speed, 50 % load
  E2  mid speed    : 50 % rated speed, 100 % load
  E3  load sweep   : rated speed, load 25 / 50 / 100 / 150 % of rated
  E4  reversal     : +rated -> -rated at 0.3 s, no load step
  E5  saturation   : quasi-step reference (1e6 rad/s^2) to rated speed + rated load -> windup test
  P   one-at-a-time parameter variations of the PLANT (controller keeps nominal model)
  MC  Monte Carlo  : N random plants/operating points
"""
from smcbench.paths import RESULTS
import json, sys
import numpy as np
from smcbench import tune as T, config as cf, metrics as mt

M = T.M
WR, TR = M['w_rated'], M['T_rated']
LR = T.LOAD_RISE


def best_params():
    res = json.load(open(f'{RESULTS}/tuning_results_{T.PROTOCOL}{T.TAG}.json'))
    best = {}
    for r in res:
        if r['name'] not in best or r['J'] < best[r['name']]['J']:
            best[r['name']] = r
    return best


def scen(w=WR, load=1.0, acc=5000.0, T=0.6, w2=0.0, t2=-1.0):
    return cf.scenario_vec(T=T, wref=w, acc=acc, TL1=load * TR, TLt=0.3, TLr=LR, wref2=w2, tref2=t2)


def plant(**f):
    m = dict(M)
    for k, v in f.items():
        m[k] = m[k] * v if k != 'B' else v
    return cf.motor_vec(m)


def evaluate(name, x, S, Pp=None, quant=20, noise=0.0, seed=1, protocol=None):
    Dv = cf.drive_vec(dict(cf.DRIVE, Vdc=M['Vdc'], Imax=M['I_max_pk'], quant=quant, noise=noise, seed=seed), M)
    log = T.run(name, x, S=S, Pp=T.P if Pp is None else Pp, Dv=Dv, protocol=protocol)
    J, r = T.cost_of(log, S)
    if r is None:
        r = {}
    r['J'] = J
    return r, log


SCEN = {
    'E0 nominal': scen(),
    'E1 low speed': scen(w=0.1 * WR, load=0.5),
    'E2 mid speed': scen(w=0.5 * WR, load=1.0),
    'E3 load 25%': scen(load=0.25), 'E3 load 50%': scen(load=0.5),
    'E3 load 150%': scen(load=1.5),
    'E4 reversal': scen(load=0.0, w2=-WR, t2=0.3, T=0.6),
    'E5 saturation': scen(acc=1e6),
}
PVAR = {'Rs x1.39 (hot)': dict(Rs=1.39), 'Rs x0.7': dict(Rs=0.7), 'L x0.8': dict(Ld=0.8, Lq=0.8),
        'L x1.2': dict(Ld=1.2, Lq=1.2), 'psi x0.947 (hot)': dict(psi=0.947), 'psi x0.9': dict(psi=0.9),
        'J x0.5': dict(J=0.5), 'J x2': dict(J=2.0), 'J x3': dict(J=3.0), 'B=0.003': dict(B=0.003)}


def monte_carlo(name, x, N=200, seed=0, protocol=None):
    rng = np.random.default_rng(seed)
    out = []
    for i in range(N):
        f = dict(Rs=rng.uniform(0.8, 1.5), Ld=(l := rng.uniform(0.8, 1.2)), Lq=l,
                 psi=rng.uniform(0.9, 1.0), J=rng.uniform(0.5, 3.0), B=rng.uniform(0, 0.003))
        S = scen(w=rng.uniform(0.1, 1.0) * WR, load=rng.uniform(0.0, 1.5))
        r, log = evaluate(name, x, S, Pp=plant(**f), noise=0.02, seed=i + 1, protocol=protocol)
        r.update({'w_ref': float(S[1]), 'load': float(S[4] / TR), **{'f_' + k: float(v) for k, v in f.items()}})
        if r.get('IAE_trk') is not None:
            t = log[:, 0]; iqr = log[:, 4]
            w1 = (t >= S[5] - mt.W_SS) & (t < S[5]); w2 = t >= S[0] - mt.W_SS
            r['TV_both'] = (np.abs(np.diff(iqr[w1])).sum() + np.abs(np.diff(iqr[w2])).sum()) / (2 * mt.W_SS)
        out.append(r)
    return out


if __name__ == '__main__':
    best = best_params()
    names = sys.argv[1].split(',') if len(sys.argv) > 1 else list(best)
    res = {}
    for n in names:
        x = best[n]['x']
        res[n] = {'scen': {k: evaluate(n, x, S)[0] for k, S in SCEN.items()},
                  'pvar': {k: evaluate(n, x, SCEN['E0 nominal'], Pp=plant(**f))[0] for k, f in PVAR.items()}}
        res[n]['mc'] = monte_carlo(n, x)
        print(n, 'done', flush=True)
        json.dump(res, open(f'{RESULTS}/eval_results_{T.PROTOCOL}{T.TAG}.json', 'w'), indent=1, default=float)
