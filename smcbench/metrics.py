"""Benchmark metrics with explicit, event-tied evaluation windows.

Log columns: t, w_ref, w, w_meas, iq_ref, iq, id, Te, TL, s, vd, vq, sat, aux
All windows are defined relative to scenario events:
  t_ramp = end of reference ramp, t_L = load-step start, T = end of run.
"""
import numpy as np

C = dict(t=0, wr=1, w=2, wm=3, iqr=4, iq=5, id=6, Te=7, TL=8, s=9, vd=10, vq=11, sat=12, aux=13)
BAND = 0.005          # +/-0.5 % of the final reference
W_SS = 0.05           # steady-state window length (s)
F_HP = 1000.0         # high-pass cutoff for HF current content (Hz)


def _settle(t, e, band, t0):
    """First time >= t0 after which |e| stays within band; NaN if never (NR)."""
    m = t >= t0
    tt, ee = t[m], np.abs(e[m])
    out = ee > band
    if not out.any():
        return 0.0
    last = np.nonzero(out)[0][-1]
    if last == len(ee) - 1:
        return np.nan
    return tt[last + 1] - t0


def _hp_rms(x, dt, fc):
    X = np.fft.rfft(x - x.mean())
    f = np.fft.rfftfreq(len(x), dt)
    X[f < fc] = 0
    return np.sqrt(np.mean(np.fft.irfft(X, len(x)) ** 2))


def compute(log, S):
    t = log[:, C['t']]
    dt = t[1] - t[0]
    wr, w = log[:, C['wr']], log[:, C['w']]
    e = wr - w
    wfin = abs(S[1])
    band = BAND * wfin
    t_ramp = wfin / S[2]
    t_L, T = S[5], S[0]

    def win(a, b):
        return (t >= a) & (t < b)

    trk = win(0, t_L)
    post_ramp = win(t_ramp, t_L)
    ld = win(t_L, T)
    ss1 = win(t_L - W_SS, t_L)          # steady window before load
    ss2 = win(T - W_SS, T)              # steady window at end
    iqr, iq, s, Te = log[:, C['iqr']], log[:, C['iq']], log[:, C['s']], log[:, C['Te']]

    r = {}
    r['IAE_trk'] = np.sum(np.abs(e[trk])) * dt
    r['ITAE_trk'] = np.sum(t[trk] * np.abs(e[trk])) * dt
    r['OS_%'] = max(0.0, (w[post_ramp].max() - wfin) / wfin * 100)
    r['t_set'] = _settle(t[trk], e[trk], band, t_ramp)
    r['IAE_load'] = np.sum(np.abs(e[ld])) * dt
    r['dw_max'] = np.max(np.abs(e[ld]))
    r['t_rec'] = _settle(t, e, band, t_L)
    r['ess_pre'] = np.mean(np.abs(e[ss1]))
    r['ess_end'] = np.mean(np.abs(e[ss2]))
    r['E_iq'] = np.sum(iq ** 2) * dt                       # A^2 s
    r['sat_%'] = 100 * np.mean(log[:, C['sat']])
    # chattering indices, evaluated on the final steady window (loaded)
    r['TV_iqref'] = np.sum(np.abs(np.diff(iqr[ss2]))) / W_SS   # A/s
    r['HF_iq'] = _hp_rms(iq[ss2], dt, F_HP)                    # A rms
    sg = np.sign(s[ss2])
    r['sw_rate'] = np.sum(sg[1:] != sg[:-1]) / W_SS            # 1/s
    r['w_p2p'] = np.ptp(w[ss2])                                # rad/s
    r['Te_p2p'] = np.ptp(Te[ss2])                              # N m
    return r
