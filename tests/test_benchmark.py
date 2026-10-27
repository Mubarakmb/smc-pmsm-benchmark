"""Regression and physics checks. Run from the repository root:  python -m pytest -q tests

The regression tests re-simulate the tuned controllers and require the metrics stored in
results/ to be reproduced; the physics tests check closed-form predictions of the paper.
"""
import json
import math
import os
import numpy as np
import pytest
from smcbench.paths import RESULTS
from smcbench import tune as T, evaluate as E
from smcbench.sim import gl_weights

ORDER = ['CSMC', 'ISMC', 'TSMC', 'ASMC', 'FOSMC', 'CSMC-BL', 'ISMC-BL', 'TSMC-BL', 'ASMC-BL', 'FOSMC-BL', 'STSMC']
RTOL = 1e-6


def _best(prot):
    best = {}
    for r in json.load(open(os.path.join(RESULTS, f'tuning_results_{prot}_final.json'))):
        if r['name'] not in best or r['J'] < best[r['name']]['J']:
            best[r['name']] = r
    return best


BEST = {p: _best(p) for p in ('lyap', 'free')}
EVAL = {p: json.load(open(os.path.join(RESULTS, f'eval_results_{p}_final.json'))) for p in ('lyap', 'free')}


@pytest.mark.parametrize('name', ORDER)
def test_training_cost_reproduced(name):
    """Training objective of the best constrained seed (tuned gains -> same J)."""
    b = BEST['lyap'][name]
    J, _ = T.cost_of(T.run(name, b['x'], protocol='lyap'), T.S_TRAIN)
    assert J == pytest.approx(b['J'], rel=RTOL)


@pytest.mark.parametrize('name', ORDER)
def test_training_cost_reproduced_free(name):
    b = BEST['free'][name]
    J, _ = T.cost_of(T.run(name, b['x'], protocol='free'), T.S_TRAIN)
    assert J == pytest.approx(b['J'], rel=RTOL)


@pytest.mark.parametrize('prot', ['lyap', 'free'])
@pytest.mark.parametrize('scen', ['E1 low speed', 'E4 reversal'])
@pytest.mark.parametrize('name', ['ISMC-BL', 'STSMC', 'FOSMC'])
def test_test_scenarios_reproduced(prot, scen, name):
    r, _ = E.evaluate(name, BEST[prot][name]['x'], E.SCEN[scen], protocol=prot)
    ref = EVAL[prot][name]['scen'][scen]
    for k in ('J', 'IAE_load', 'Te_p2p', 'TV_iqref'):
        assert r[k] == pytest.approx(ref[k], rel=RTOL, abs=1e-12), k


def test_monte_carlo_first_runs_reproduced():
    """First three Monte Carlo draws (random plant, operating point and noise) of ISMC-BL."""
    mc = E.monte_carlo('ISMC-BL', BEST['lyap']['ISMC-BL']['x'], N=3)
    for r, ref in zip(mc, EVAL['lyap']['ISMC-BL']['mc'][:3]):
        assert r['J'] == pytest.approx(ref['J'], rel=RTOL)
        assert r['f_J'] == ref['f_J'] and r['load'] == ref['load']


@pytest.mark.parametrize('name', ['CSMC-BL', 'ISMC-BL', 'TSMC-BL', 'FOSMC-BL'])
def test_boundary_layer_equilibrium(name):
    """Inside the boundary layer the sliding variable settles at s* = D / (k/phi + q)."""
    p = BEST['lyap'][name]['params']
    L = T.run(name, BEST['lyap'][name]['x'], S=E.SCEN['E0 nominal'])
    s_end = L[L[:, 0] > 0.55, 9].mean()
    assert s_end == pytest.approx(T.D / (p['k'] / p['phi'] + p['q']), rel=2e-3)


@pytest.mark.parametrize('name', ['ISMC', 'TSMC', 'FOSMC'])
def test_unconstrained_gain_no_sliding(name):
    """With k < D (unconstrained protocol) s stays at (D - k)/q and never changes sign
    (paper, Section 5.5: 6.14, 6.56 and 5.33 for ISMC, TSMC and FOSMC)."""
    p = BEST['free'][name]['params']
    assert p['k'] < T.D
    L = T.run(name, BEST['free'][name]['x'], S=E.SCEN['E0 nominal'], protocol='free')
    s = L[L[:, 0] > 0.55, 9]
    assert s.min() > 0
    assert s.mean() == pytest.approx((T.D - p['k']) / p['q'], rel=5e-3)   # FOSMC: 0.2 % memory effect


def test_steady_state_torque_balance():
    """Rated load, no friction: mean electromagnetic torque equals the load torque."""
    L = T.run('ISMC-BL', BEST['lyap']['ISMC-BL']['x'], S=E.SCEN['E0 nominal'])
    m = L[:, 0] > 0.55
    assert L[m, 7].mean() == pytest.approx(T.M['T_rated'], rel=1e-3)


def test_limits_respected():
    """|i_q*| never exceeds the data-sheet maximum current, also in the saturation scenario."""
    L = T.run('ISMC-BL', BEST['lyap']['ISMC-BL']['x'], S=E.SCEN['E5 saturation'])
    assert np.abs(L[:, 4]).max() <= T.M['I_max_pk'] + 1e-9
    v = np.hypot(L[:, 10], L[:, 11])
    assert v.max() <= T.M['Vdc'] / math.sqrt(3) + 1e-9


def test_gl_weights():
    """Grunwald-Letnikov weights: order 0 is the identity, order 1 the backward difference."""
    assert np.allclose(gl_weights(0.0, 5), [1, 0, 0, 0, 0])
    assert np.allclose(gl_weights(1.0, 5), [1, -1, 0, 0, 0])
