"""(a) Solver step-size convergence (R3-7): nominal scenario, h = 10, 5, 2, 1 us.
(b) Sampling-time sensitivity (R1-2): Ts = 50, 100, 200 us with the gains tuned at 100 us.
Uses the constrained-protocol tuned gains (best seed)."""
from smcbench.paths import RESULTS
import json
import numpy as np
from smcbench import tune as T, config as cf, evaluate as E

best = E.best_params()
S = E.SCEN['E0 nominal']
keys = ['IAE_trk', 'IAE_load', 'dw_max', 'Te_p2p', 'TV_iqref']
out = {}
print('(a) solver step: max relative change of metrics vs h = 1 us')
for n, b in best.items():
    vals = {}
    for h in [10e-6, 5e-6, 2e-6, 1e-6]:
        Dv = cf.drive_vec(dict(cf.DRIVE, Vdc=T.M['Vdc'], Imax=T.M['I_max_pk'], quant=20, h=h), T.M)
        _, r = T.cost_of(T.run(n, b['x'], S=S, Dv=Dv), S)
        vals[h] = r
    ref = vals[1e-6]
    dev = {h: max(abs(vals[h][k] - ref[k]) / max(abs(ref[k]), 1e-12) for k in keys) for h in vals}
    out.setdefault(n, {})['h'] = dev
    print(f"  {n:9s} " + '  '.join(f"h={h*1e6:3.0f}us: {100*d:7.4f}%" for h, d in dev.items()))
print('\n(b) sampling time: cost J and torque ripple (N m p-p), gains tuned at 100 us')
for n, b in best.items():
    row = []
    for Ts in [50e-6, 100e-6, 200e-6]:
        Dv = cf.drive_vec(dict(cf.DRIVE, Vdc=T.M['Vdc'], Imax=T.M['I_max_pk'], quant=20, Tsc=Ts), T.M)
        J, r = T.cost_of(T.run(n, b['x'], S=S, Dv=Dv), S)
        row.append((Ts, J, r['Te_p2p'] if r else np.nan))
    out[n]['Ts'] = row
    print(f"  {n:9s} " + '  '.join(f"Ts={Ts*1e6:3.0f}us: J={J:7.3f} ripple={rp:6.3f}" for Ts, J, rp in row))
json.dump(out, open(f'{RESULTS}/sensitivity_lyap.json', 'w'), indent=1, default=float)
