"""Non-constant disturbances that integral action cannot remove exactly:
  (i)  periodic torque ripple at 6x electrical frequency, amplitude 1/3/5 % of rated torque,
       evaluated at low speeds where its frequency falls inside the speed-loop bandwidth;
  (ii) Coulomb friction 2/5 % of rated torque during a low-speed reversal (+/-5 % rated speed).
Amplitudes are ASSUMPTIONS (not in the datasheet) and are swept rather than fixed.
Metric: speed ripple p-p and mean |e| over the last 0.1 s (ripple) / IAE around the reversal."""
from smcbench.paths import RESULTS
import json
import numpy as np
from smcbench import tune as T, config as cf, evaluate as E

best = E.best_params()
names = ['ISMC-BL', 'TSMC-BL', 'FOSMC-BL', 'ASMC-BL', 'STSMC', 'ISMC', 'CSMC-BL']
TR, WR = T.M['T_rated'], T.M['w_rated']
Dv = cf.drive_vec(dict(cf.DRIVE, Vdc=T.M['Vdc'], Imax=T.M['I_max_pk'], quant=20), T.M)
out = {}
print('(i) periodic ripple: speed ripple p-p (rad/s) over last 0.1 s, 50 % load')
for wfrac in [0.02, 0.1]:
    for amp in [0.01, 0.03, 0.05]:
        Pp = cf.motor_vec(dict(T.M, Trip=amp * TR))
        S = cf.scenario_vec(T=0.8, wref=wfrac * WR, acc=5000.0, TL1=0.5 * TR, TLt=0.3, TLr=5e-3)
        row = []
        for n in names:
            L = T.run(n, best[n]['x'], S=S, Pp=Pp, Dv=Dv)
            m = L[:, 0] >= S[0] - 0.1
            e = L[m, 1] - L[m, 2]
            row.append((n, np.ptp(L[m, 2]), np.mean(np.abs(e))))
            out.setdefault(f'ripple w={wfrac} a={amp}', {})[n] = (float(np.ptp(L[m, 2])), float(np.mean(np.abs(e))))
        f_hz = 6 * T.M['pn'] * wfrac * WR / (2 * np.pi)
        print(f"  speed {wfrac*100:3.0f}% (ripple {f_hz:5.0f} Hz) amp {amp*100:.0f}%: " +
              '  '.join(f"{n}={p:.4f}" for n, p, _ in row))
print('\n(ii) Coulomb friction, reversal +/-5 % rated speed: IAE over [0.25, 0.45] s (rad)')
for tf in [0.0, 0.02, 0.05]:
    Pp = cf.motor_vec(dict(T.M, Tf=tf * TR))
    S = cf.scenario_vec(T=0.6, wref=0.05 * WR, acc=500.0, TL1=0.0, TLt=0.59, TLr=5e-3,
                        wref2=-0.05 * WR, tref2=0.3)
    row = []
    for n in names:
        L = T.run(n, best[n]['x'], S=S, Pp=Pp, Dv=Dv)
        m = (L[:, 0] >= 0.25) & (L[:, 0] < 0.45)
        iae = np.sum(np.abs(L[m, 1] - L[m, 2])) * (L[1, 0] - L[0, 0])
        row.append((n, iae))
        out.setdefault(f'coulomb {tf}', {})[n] = float(iae)
    print(f"  Tf = {tf*100:.0f}% rated: " + '  '.join(f"{n}={v:.4f}" for n, v in row))
json.dump(out, open(f'{RESULTS}/disturb_lyap.json', 'w'), indent=1)
