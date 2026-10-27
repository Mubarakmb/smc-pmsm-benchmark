"""Benchmark motor, drive settings and scenario definition.

Motor: Siemens SIMOTICS S-1FK7063-5AF71 surface-mounted PMSM; every value below is taken
from (or derived from) the manufacturer data sheet, see the notes next to SIEMENS_1FK7063.
"""
import numpy as np

# Siemens SIMOTICS S-1FK7063-5AF71 (manufacturer data sheet, Article 1FK7063-5AF71-1EA2):
#   8 poles, 3000 rpm rated, 7.3 Nm / 5.6 A rated (100 K), 35 Nm / 28 A max, kT = 1.37 Nm/A,
#   kE(20C) = 87.5 V/1000 rpm, R(20C) = 0.65 ohm, rotating-field L = 7.7 mH, J = 15.1 kg cm^2,
#   tau_el = 11.8 ms, tau_mech = 1.56 ms (both reproduced from R, L, J, kT -> R, L are per phase).
#   psi_f from kE (20 C): 0.1706 Wb;  from kT (100 K): 0.1615 Wb (-5.3 %, hot magnets).
#   Viscous friction is NOT given in the data sheet -> B = 0 nominal (varied in robustness study).
#   Rated data of the 1FK7 family are specified for a 600 V DC link.
SIEMENS_1FK7063 = dict(Rs=0.65, Ld=7.7e-3, Lq=7.7e-3, psi=0.1706, pn=4, J=15.1e-4, B=0.0,
                       I_rated_pk=5.6 * 2 ** 0.5, I_max_pk=28.0 * 2 ** 0.5, T_rated=7.3,
                       w_rated=3000 * 3.141592653589793 / 30, Vdc=600.0)


def motor_vec(m):
    return np.array([m['Rs'], m['Ld'], m['Lq'], m['psi'], m['pn'], m['J'], m['B'],
                     m.get('Trip', 0.0), m.get('Tf', 0.0)], float)


# Drive defaults. All benchmark scripts override Vdc, Imax and quant explicitly with the
# data-sheet values (Vdc = 600 V, Imax = 28 A rms * sqrt(2), 20-bit encoder).
#   Tsc  : sampling period of the speed and current loops (s)
#   nw   : speed-loop period in current-loop periods
#   h    : RK4 integration step of the continuous plant (s)
#   wc   : current-loop bandwidth (rad/s), PI by pole-zero cancellation
#   Kaw  : back-calculation anti-windup gain of the current PI (None -> wc)
#   delay: one-sample computational delay (1) or none (0)
#   noise: std of additive speed-measurement noise (rad/s)
#   quant: encoder resolution in bits (0 -> ideal speed measurement)
#   seed : noise seed
DRIVE = dict(Vdc=600.0, Imax=28.0 * 2 ** 0.5, Tsc=100e-6, nw=1, h=2e-6,
             wc=2 * np.pi * 400, Kaw=None, delay=1, noise=0.0, quant=0.0, seed=1)


def drive_vec(d, m):
    Kaw = d['Kaw'] if d['Kaw'] is not None else d['wc']   # back-calculation gain, 1/s
    return np.array([d['Vdc'], d['Imax'], d['Tsc'], d['nw'], d['h'], d['wc'], Kaw,
                     d['delay'], d['noise'], d['quant'], d['seed']], float)


def scenario_vec(T=0.8, wref=300.0, acc=3000.0, TL0=0.0, TL1=3.0, TLt=0.4, TLr=5e-3,
                 wref2=0.0, tref2=-1.0):
    """T: run length (s); wref: speed reference (rad/s, mech.) reached by a ramp of slope acc
    (rad/s^2); TL0: initial load (N m); TL1: load step (N m) starting at TLt (s) with a
    raised-cosine rise of TLr (s); optional second ramp to wref2 starting at tref2 (s, <0 = none)."""
    return np.array([T, wref, acc, TL0, TL1, TLt, TLr, wref2, tref2], float)


CTRL = {1: 'CSMC', 2: 'ISMC', 3: 'TSMC', 4: 'ASMC', 5: 'FOSMC', 6: 'STSMC'}
