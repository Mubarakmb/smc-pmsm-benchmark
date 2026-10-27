"""Accuracy of the short-memory Grunwald-Letnikov fractional integral vs memory length.
Test signal: unit step e(t)=1 (t>=0), exact D^{-a} 1 = t^a / Gamma(1+a); also e(t)=t:
exact t^{1+a}/Gamma(2+a). Sampling Ts = 100 us, horizon 0.6 s."""
import numpy as np
from math import gamma
from smcbench.sim import gl_weights
Ts, T = 1e-4, 0.6
n = int(T / Ts) + 1
t = np.arange(n) * Ts
for a in [0.5, 0.89]:
    for name, f, exact in [('step', np.ones(n), t ** a / gamma(1 + a)),
                           ('ramp', t.copy(), t ** (1 + a) / gamma(2 + a))]:
        print(f"alpha={a}  signal={name}")
        for Lm in [100, 300, 1000, 2000, 6001]:
            w = gl_weights(-a, Lm)
            k = n - 1
            j = np.arange(min(k + 1, Lm))
            val = Ts ** a * np.sum(w[j] * f[k - j])
            err = abs(val - exact[-1]) / abs(exact[-1]) * 100
            print(f"   Lm={Lm:5d} ({Lm*Ts*1e3:6.1f} ms window)  rel. error at t=0.6 s = {err:7.3f} %  "
                  f"cost = {2*Lm} MAC/step (2 operators), memory = {8*Lm/1024:.1f} kB/buffer (double)")
