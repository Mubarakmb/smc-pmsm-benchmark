"""
Unified PMSM speed-control benchmark: numba-accelerated simulation core.

Architecture (identical for all six controllers):
  speed SMC (Ts_w) -> i_q* (saturated, common anti-windup)
  -> PI current loops (Ts_c), i_d* = 0, decoupling + back-EMF feedforward,
     back-calculation anti-windup
  -> voltage limit |v_dq| <= Vdc/sqrt(3) (d-axis priority)
  -> one-sample computational delay -> average-value inverter (ZOH)
  -> continuous dq PMSM model integrated with fixed-step RK4 (step h).

The disturbance d = T_L/J is NOT available to any controller; all controllers
use only the nominal model (chi_hat, a_hat) and the known reference.

Controller output structure (common): chi_hat * i_q* = dw_ref + a_hat*w + N(.)
  1 CSMC : s = e                               N = k*sw(s) + q*s
  2 ISMC : s = e - e(0) + lam*int(e)           N = lam*e + k*sw(s) + q*s
  3 TSMC : s = e + nu*int(|e|^al sgn e)        N = nu|e|^al sgn e + k*sw(s) + q*s
           (integral terminal surface: no |e|^(al-1) term -> nonsingular)
  4 ASMC : s = e                               N = K(t)*sw(s) + q*s,
           K adapted by equivalent-control (LPF of sw) law with bounds
  5 FOSMC: s = e + ki*D^{-al} e  (GL, short memory Lm)
                                               N = ki*D^{1-al} e + k*sw(s) + q*s
  6 STSMC: s = e                               N = k1*|s|^0.5 sgn s + v,  dv/dt = k2*sgn(s)
           (gain vector G[0] = k1, G[1] = k2)
  sw(s) = sgn(s)  (bl = 0)   or  sat(s/phi)  (bl = 1, boundary layer phi)
"""
import numpy as np
from numba import njit

# ---- parameter vector layout ------------------------------------------------
# plant P (true) and nominal Pn (used by controllers) share layout:
P_RS, P_LD, P_LQ, P_PSI, P_PN, P_J, P_B, P_TRIP, P_TF = range(9)
# P_TRIP: amplitude of periodic torque ripple at 6x electrical frequency (N m)
# P_TF  : Coulomb friction torque (N m), smoothed as Tf*tanh(w/0.5)
# drive D
(D_VDC, D_IMAX, D_TSC, D_NW, D_H, D_WC, D_KAW, D_DELAY,
 D_NOISE, D_QUANT, D_SEED) = range(11)
# scenario S
(S_T, S_WREF, S_ACC, S_TL0, S_TL1, S_TLT, S_TLR, S_WREF2, S_TREF2) = range(9)
# gains G (meaning depends on controller; see controller_gains())
NG = 12


@njit(cache=True)
def _sgn(x):
    return 1.0 if x > 0 else (-1.0 if x < 0 else 0.0)


@njit(cache=True)
def _sw(s, bl, phi):
    if bl == 0:
        return _sgn(s)
    y = s / phi
    return 1.0 if y > 1.0 else (-1.0 if y < -1.0 else y)


@njit(cache=True)
def _ref(t, S):
    """Speed reference: ramp to S_WREF at S_ACC, optional second ramp to
    S_WREF2 starting at S_TREF2 (same |acc|). Returns (w*, dw*)."""
    w1, acc = S[S_WREF], S[S_ACC]
    t1 = abs(w1) / acc
    if t < t1:
        w, dw = acc * t * _sgn(w1), acc * _sgn(w1)
    else:
        w, dw = w1, 0.0
    if S[S_TREF2] > 0 and t >= S[S_TREF2]:
        w2 = S[S_WREF2]
        dt = t - S[S_TREF2]
        span = abs(w2 - w1) / acc
        if dt < span:
            w, dw = w1 + _sgn(w2 - w1) * acc * dt, _sgn(w2 - w1) * acc
        else:
            w, dw = w2, 0.0
    return w, dw


@njit(cache=True)
def _load(t, S):
    """Load torque: base TL0, raised-cosine step of TL1 at TLT over rise TLR
    (bounded dT/dt and d2T/dt2)."""
    TL = S[S_TL0]
    if t >= S[S_TLT]:
        x = (t - S[S_TLT]) / S[S_TLR] if S[S_TLR] > 0 else 1.0
        if x >= 1.0:
            TL += S[S_TL1]
        else:
            TL += S[S_TL1] * 0.5 * (1.0 - np.cos(np.pi * x))
    return TL


@njit(cache=True)
def _deriv(x, vd, vq, TL, P):
    id_, iq, w = x[0], x[1], x[2]
    we = P[P_PN] * w
    did = (vd - P[P_RS] * id_ + we * P[P_LQ] * iq) / P[P_LD]
    diq = (vq - P[P_RS] * iq - we * (P[P_LD] * id_ + P[P_PSI])) / P[P_LQ]
    Te = 1.5 * P[P_PN] * (P[P_PSI] * iq + (P[P_LD] - P[P_LQ]) * id_ * iq)
    Td = 0.0
    if P.shape[0] > P_TF:
        Td = P[P_TRIP] * np.sin(6.0 * P[P_PN] * x[3]) + P[P_TF] * np.tanh(w / 0.5)
    dw = (Te - P[P_B] * w - TL - Td) / P[P_J]
    out = np.empty(4)
    out[0], out[1], out[2], out[3] = did, diq, dw, w
    return out


@njit(cache=True)
def gl_weights(r, Lm):
    w = np.empty(Lm)
    w[0] = 1.0
    for j in range(1, Lm):
        w[j] = w[j - 1] * (1.0 - (r + 1.0) / j)
    return w


@njit(cache=True)
def simulate(ctrl, G, bl, phi, P, Pn, D, S, wa, wb, Lm):
    """Run one simulation. wa/wb: GL weights for orders -al and 1-al (FOSMC).
    Returns log array sampled at the speed-loop rate with columns:
    t, w_ref, w, w_meas, iq_ref, iq, id, Te, TL, s, vd, vq, sat, aux"""
    Tsc = D[D_TSC]
    nw = int(D[D_NW])
    Tsw = Tsc * nw
    h = D[D_H]
    nsub = int(round(Tsc / h))
    Ttot = S[S_T]
    nc = int(round(Ttot / Tsc))
    nlog = nc // nw + 1
    log = np.zeros((nlog, 14))

    np.random.seed(int(D[D_SEED]))
    # nominal model constants for controllers
    chi = 1.5 * Pn[P_PN] * Pn[P_PSI] / Pn[P_J]
    ah = Pn[P_B] / Pn[P_J]
    Imax = D[D_IMAX]
    Vmax = D[D_VDC] / np.sqrt(3.0)
    # current PI (pole-zero cancellation, bandwidth wc)
    wc = D[D_WC]
    Kpd, Kid = Pn[P_LD] * wc, Pn[P_RS] * wc
    Kpq, Kiq = Pn[P_LQ] * wc, Pn[P_RS] * wc
    Kaw = D[D_KAW]

    x = np.zeros(4)           # id, iq, w, theta (true, mech)
    thq_prev = 0.0
    xd = 0.0
    xq = 0.0                  # PI integrators
    vd_app = 0.0
    vq_app = 0.0              # applied (delayed) voltages
    vd_next = 0.0
    vq_next = 0.0
    iq_ref = 0.0
    # controller states
    I1 = 0.0                  # integral state (ISMC/TSMC)
    e0 = 0.0
    first = True
    K = G[2]                  # ASMC initial gain
    z = 0.0                   # ASMC LPF
    v = 0.0                   # STSMC integrator
    ebuf = np.zeros(Lm)       # FOSMC error history (circular)
    bidx = 0
    nb = 0
    s = 0.0
    sat = 0.0
    aux = 0.0
    il = 0

    for kc in range(nc):
        t = kc * Tsc
        # ---------------- measurement ----------------
        if D[D_QUANT] > 0:
            # speed from an N-bit encoder by backward difference at Ts_c
            res = 2.0 * np.pi / 2.0 ** D[D_QUANT]
            thq = res * np.floor(x[3] / res)
            w_meas = (thq - thq_prev) / Tsc if kc > 0 else 0.0
            thq_prev = thq
        else:
            w_meas = x[2]
        if D[D_NOISE] > 0:
            w_meas += D[D_NOISE] * np.random.randn()
        # ---------------- speed loop -----------------
        if kc % nw == 0:
            wr, dwr = _ref(t, S)
            e = wr - w_meas
            if first:
                e0 = e
                first = False
            if ctrl == 1:  # CSMC
                k, q = G[0], G[1]
                s = e
                N = k * _sw(s, bl, phi) + q * s
            elif ctrl == 2:  # ISMC
                lam, k, q = G[0], G[1], G[3]
                s = e - e0 + lam * I1
                N = lam * e + k * _sw(s, bl, phi) + q * s
            elif ctrl == 3:  # TSMC (integral terminal)
                nu, al, k, q = G[0], G[1], G[3], G[4]
                ft = nu * abs(e) ** al * _sgn(e)
                s = e + I1
                N = ft + k * _sw(s, bl, phi) + q * s
            elif ctrl == 4:  # ASMC
                q = G[1]
                s = e
                N = K * _sw(s, bl, phi) + q * s
            elif ctrl == 5:  # FOSMC
                ki, k, q = G[0], G[2], G[3]
                ebuf[bidx] = e
                nb = min(nb + 1, Lm)
                accA = 0.0
                accB = 0.0
                for j in range(nb):
                    ej = ebuf[(bidx - j) % Lm]
                    accA += wa[j] * ej
                    accB += wb[j] * ej
                al = G[1]
                Dm = accA * Tsw ** al            # D^{-al} e
                D1 = accB * Tsw ** (al - 1.0)    # D^{1-al} e
                bidx = (bidx + 1) % Lm
                s = e + ki * Dm
                N = ki * D1 + k * _sw(s, bl, phi) + q * s
                aux = Dm
            else:  # STSMC
                k1 = G[0]
                s = e
                N = k1 * np.sqrt(abs(s)) * _sgn(s) + v
            u = (dwr + ah * w_meas + N) / chi
            uc = min(max(u, -Imax), Imax)
            sat = 1.0 if uc != u else 0.0
            # common conditional-integration anti-windup:
            # freeze integral states if saturated and e drives further into it
            freeze = sat > 0 and _sgn(u) == _sgn(e)
            if ctrl == 2 and not freeze:
                I1 += e * Tsw
            elif ctrl == 3 and not freeze:
                I1 += G[0] * abs(e) ** G[1] * _sgn(e) * Tsw
            elif ctrl == 4:
                tau, alf, g0, g1, M, mu, Kp_ = G[3], G[4], G[5], G[6], G[7], G[8], G[9]
                z += Tsw / tau * (_sw(s, bl, phi) - z)
                dl = abs(z) - alf
                dK = (g0 + g1 * abs(e)) * K * _sgn(dl) \
                    - M * max(K - Kp_, 0.0) + M * max(mu - K, 0.0)
                if not (freeze and dK > 0):
                    K += dK * Tsw
                aux = K
            elif ctrl == 6 and not freeze:
                v += G[1] * _sgn(s) * Tsw
                aux = v
            iq_ref = uc
            if il < nlog:
                log[il, 0] = t
                log[il, 1] = wr
                log[il, 2] = x[2]
                log[il, 3] = w_meas
                log[il, 4] = iq_ref
                log[il, 5] = x[1]
                log[il, 6] = x[0]
                log[il, 7] = 1.5 * P[P_PN] * (P[P_PSI] * x[1] + (P[P_LD] - P[P_LQ]) * x[0] * x[1])
                log[il, 8] = _load(t, S)
                log[il, 9] = s
                log[il, 10] = vd_app
                log[il, 11] = vq_app
                log[il, 12] = sat
                log[il, 13] = aux
                il += 1
        # ---------------- current loop ---------------
        we = Pn[P_PN] * w_meas
        ed = 0.0 - x[0]
        eq = iq_ref - x[1]
        vd_u = Kpd * ed + xd - we * Pn[P_LQ] * x[1]
        vq_u = Kpq * eq + xq + we * (Pn[P_LD] * x[0] + Pn[P_PSI])
        vd_c = min(max(vd_u, -Vmax), Vmax)
        vqm = np.sqrt(max(Vmax * Vmax - vd_c * vd_c, 0.0))
        vq_c = min(max(vq_u, -vqm), vqm)
        xd += Tsc * (Kid * ed + Kaw * (vd_c - vd_u))
        xq += Tsc * (Kiq * eq + Kaw * (vq_c - vq_u))
        if D[D_DELAY] > 0:
            vd_app, vq_app = vd_next, vq_next
            vd_next, vq_next = vd_c, vq_c
        else:
            vd_app, vq_app = vd_c, vq_c
        # ---------------- plant (RK4) ----------------
        for ks in range(nsub):
            tt = t + ks * h
            TL1 = _load(tt, S)
            TL2 = _load(tt + 0.5 * h, S)
            TL3 = _load(tt + h, S)
            k1 = _deriv(x, vd_app, vq_app, TL1, P)
            k2 = _deriv(x + 0.5 * h * k1, vd_app, vq_app, TL2, P)
            k3 = _deriv(x + 0.5 * h * k2, vd_app, vq_app, TL2, P)
            k4 = _deriv(x + h * k3, vd_app, vq_app, TL3, P)
            x = x + h / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
    return log[:il]
