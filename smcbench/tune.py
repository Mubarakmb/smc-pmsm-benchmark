"""Common tuning protocol (identical for every controller variant).

Optimiser   : scipy differential evolution (best1bin, dithered mutation 0.5-1, crossover 0.7,
              deferred updating), polish off, Sobol init.
Budget      : exactly N_EVAL = 1920 cost evaluations per run (Sobol population of POP = 32,
              60 generations), identical for every variant; NSEED = 3 seeds per variant.
              (A run can stop earlier only if its population collapses to a single point.)
Protocols   : PROTOCOL=lyap (default) enforces the reaching/stability conditions at rated
              load as lower bounds; PROTOCOL=free lowers the switching-gain bounds to 0.1 D.
Usage       : python -m smcbench.tune [CSMC,ISMC-BL,...]   (env: PROTOCOL, TAG, LOAD_RISE, N_WORKERS)
              Finished (variant, seed) jobs found in the result file are skipped (resumable).
Search space: log-uniform for gains, linear for exponents; ranges scaled by the rated-load
              disturbance bound D = T_rated / J (same rule for every controller).
Training scenario (not used for reporting): ramp 0 -> rated speed at 5000 rad/s^2,
              rated-load raised-cosine step at 0.3 s (5 ms rise), 20-bit encoder speed,
              one-sample delay, current/voltage limits; T = 0.6 s.
Cost        : J = IAE_trk/N1 + IAE_load/N2 + TV/N3 + 10*sat_frac + 100*[unstable]
              N1 = N2 = (0.5 % of rated speed) x window length  (error sitting on the band edge)
              N3 = 2 * I_rated_pk / 1 ms                         (rated-current swing every ms)
              TV = total variation of iq* per second over both steady windows.
"""
from smcbench.paths import RESULTS
import sys, json, time, math, os
import numpy as np
from scipy.optimize import differential_evolution
from smcbench.sim import simulate, NG, gl_weights
from smcbench import config as cf, metrics as mt

M = cf.SIEMENS_1FK7063
P = cf.motor_vec(M)
DRV = cf.drive_vec(dict(cf.DRIVE, Vdc=M['Vdc'], Imax=M['I_max_pk'], quant=20), M)
TAG = os.environ.get('TAG', '_final')   # result-file suffix; '_final' = the results reported in the paper
LOAD_RISE = float(os.environ.get('LOAD_RISE', '5e-3'))   # raised-cosine load rise time (s)
S_TRAIN = cf.scenario_vec(T=0.6, wref=M['w_rated'], acc=5000.0, TL1=M['T_rated'], TLt=0.3, TLr=LOAD_RISE)
D = M['T_rated'] / M['J']
LM = 2000
N_EVAL, POP, NSEED = 1920, 32, 3
PROTOCOL = os.environ.get('PROTOCOL', 'lyap')   # 'lyap' (stability conditions enforced) or 'free'
C_DOT = D * np.pi / 2 / LOAD_RISE                       # |d_dot| bound for the 5 ms raised-cosine load step
KLO = 1.05 * D if PROTOCOL == 'lyap' else 0.1 * D

# name: (ctrl id, boundary layer, [(param, lo, hi, log)])
SPECS = {
    'CSMC':     (1, 0, [('k', KLO, 10*D, 1), ('q', 10, 3000, 1)]),
    'ISMC':     (2, 0, [('lam', 5, 1000, 1), ('k', KLO, 10*D, 1), ('q', 10, 3000, 1)]),
    'TSMC':     (3, 0, [('nu', 5, 1000, 1), ('al', .5, .95, 0), ('k', KLO, 10*D, 1), ('q', 10, 3000, 1)]),
    'ASMC':     (4, 0, [('q', 10, 3000, 1), ('Kp', KLO, 10*D, 1), ('g0', 1, 1000, 1), ('tau', 1e-4, 1e-2, 1), ('alf', .5, .95, 0)]),
    'FOSMC':    (5, 0, [('ki', 5, 1000, 1), ('al', .1, .9, 0), ('k', KLO, 10*D, 1), ('q', 10, 3000, 1)]),
}
# STSMC: under 'lyap' the conditions k2 > C and k1 > sqrt(k2 + C) (Seeber & Horn 2017) are
# enforced by searching rho = k1 / sqrt(k2 + C) >= 1.01 and k2 >= 1.05 C; under 'free' k1, k2 directly.
# The parameter vector x therefore has a different meaning for STSMC in the two protocols.
STSMC_SPEC = {'lyap': (6, 0, [('rho', 1.01, 10, 1), ('k2', 1.05 * C_DOT, 100 * C_DOT, 1)]),
              'free': (6, 0, [('k1', 1, 1e4, 1), ('k2', 1e2, 1e8, 1)])}
SPECS['STSMC'] = STSMC_SPEC[PROTOCOL]
for n in ['CSMC', 'ISMC', 'TSMC', 'ASMC', 'FOSMC']:   # boundary-layer variants: + phi
    c, _, sp = SPECS[n]
    SPECS[n + '-BL'] = (c, 1, sp + [('phi', 1e-2, 10, 1)])


def build(name, x, protocol=None):
    """Map a parameter vector x (physical values, order of SPECS) to the simulator gain vector."""
    c, bl, sp = STSMC_SPEC[protocol or PROTOCOL] if name == 'STSMC' else SPECS[name]
    p = {k: v for (k, *_), v in zip(sp, x)}
    G = np.zeros(NG)
    if c == 1:   G[:2] = [p['k'], p['q']]
    elif c == 2: G[[0, 1, 3]] = [p['lam'], p['k'], p['q']]
    elif c == 3: G[[0, 1, 3, 4]] = [p['nu'], p['al'], p['k'], p['q']]
    elif c == 4: G[:10] = [0, p['q'], 0.01 * p['Kp'], p['tau'], p['alf'], p['g0'], 0.0, 1e3, 0.01 * p['Kp'], p['Kp']]
    elif c == 5: G[:4] = [p['ki'], p['al'], p['k'], p['q']]
    else:
        k1 = p['k1'] if 'k1' in p else p['rho'] * math.sqrt(p['k2'] + C_DOT)
        G[:2] = [k1, p['k2']]
        p['k1'] = k1
    phi = p.get('phi', 1.0)
    return c, bl, G, phi, p


def to_x(name, u):
    _, _, sp = SPECS[name]
    return [math.exp(math.log(lo) + ui * (math.log(hi) - math.log(lo))) if lg else lo + ui * (hi - lo)
            for (k, lo, hi, lg), ui in zip(sp, u)]


def run(name, x, S=S_TRAIN, Pp=P, Dv=DRV, protocol=None):
    """Simulate variant `name` with parameters x (as stored in the tuning results) in scenario S.
    Pp: true plant vector (the controller always uses the nominal P); Dv: drive vector.
    protocol: 'lyap' or 'free' (only matters for STSMC); default = environment PROTOCOL."""
    c, bl, G, phi, _ = build(name, x, protocol)
    al = G[1] if c == 5 else 0.5
    wa, wb = (gl_weights(-al, LM), gl_weights(1 - al, LM)) if c == 5 else (np.ones(1), np.ones(1))
    return simulate(c, G, bl, phi, Pp, P.copy(), Dv, S, wa, wb, len(wa))


def cost_of(log, S):
    if not np.all(np.isfinite(log)):
        return 1e3, None
    r = mt.compute(log, S)
    band = 0.005 * abs(S[1])
    N1 = band * (S[5] - 0.0)
    N2 = band * (S[0] - S[5])
    N3 = 2 * M['I_rated_pk'] / 1e-3
    t = log[:, 0]; iqr = log[:, 4]
    w1 = (t >= S[5] - mt.W_SS) & (t < S[5]); w2 = t >= S[0] - mt.W_SS
    TV = (np.abs(np.diff(iqr[w1])).sum() + np.abs(np.diff(iqr[w2])).sum()) / (2 * mt.W_SS)
    J = r['IAE_trk'] / N1 + r['IAE_load'] / N2 + TV / N3 + 10 * r['sat_%'] / 100
    if r['ess_end'] > 0.05 * abs(S[1]) or r['ess_pre'] > 0.05 * abs(S[1]):
        J += 100
    return J, r


def objective(u, name):
    J, _ = cost_of(run(name, to_x(name, u)), S_TRAIN)
    return J


def tune(name, seed):
    d = len(SPECS[name][2])
    popsize = POP // d                    # popsize*d in (16, 32] -> Sobol init gives exactly POP members
    maxiter = N_EVAL // POP - 1           # identical budget: POP * (maxiter + 1) = N_EVAL
    t0 = time.time()
    hist = []
    res = differential_evolution(objective, [(0, 1)] * d, args=(name,), popsize=popsize, maxiter=maxiter,
                                 seed=seed, polish=False, init='sobol', tol=0, mutation=(0.5, 1.0),
                                 recombination=0.7, updating='deferred', workers=1,
                                 atol=0)
    x = to_x(name, res.x)
    J, r = cost_of(run(name, x), S_TRAIN)
    return dict(name=name, seed=seed, J=J, nfev=int(res.nfev), x=x,
                params=build(name, x)[4], metrics=r, secs=time.time() - t0, hist=hist)


def _tune_star(a):
    return tune(*a)


if __name__ == '__main__':
    from multiprocessing import Pool
    names = sys.argv[1].split(',') if len(sys.argv) > 1 else list(SPECS)
    jobs = [(n, s) for n in names for s in range(NSEED)]
    out = []
    fn = f'{RESULTS}/tuning_results_{PROTOCOL}{TAG}.json'
    if os.path.exists(fn):                      # resume: skip finished (name, seed) jobs
        out = json.load(open(fn))
        done = {(r['name'], r['seed']) for r in out}
        jobs = [j for j in jobs if j not in done]
    with Pool(int(os.environ.get('N_WORKERS', os.cpu_count() or 1))) as pool:
        for r in pool.imap_unordered(_tune_star, jobs):
            print(f"{r['name']:9s} seed={r['seed']} J={r['J']:.4f} nfev={r['nfev']} {r['secs']:.0f}s", flush=True)
            out.append(r)
            json.dump(out, open(f'{RESULTS}/tuning_results_{PROTOCOL}{TAG}.json', 'w'), indent=1, default=float)
