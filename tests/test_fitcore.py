"""The SIM155 optimiser (common/fitcore.py) must find the minimum of noisy, badly conditioned problems (ROADMAP 3.4).

Synthetic −log(post): Gaussian with correlations up to 0.995 and widths spanning five decades, plus the evaluation
noise measured with hi_class + Planck/DESI/Pantheon+ (0.035 RMS in −log(post); SIM155/diagnostics.py) and random
hi_class failures. Pass = within 0.15–0.3 of the true minimum (≤ 0.6 in χ²), far inside the pre-registered
thresholds (Δχ² ≥ 6; "competitive" = within 3 of w0wa). Observed: 0.01–0.08.
"""
import numpy as np
from scipy.optimize import minimize as sp_minimize

import fitcore


def problem(n, seed, noise=0.035, fail_rate=0.0, edge=False):
    rng = np.random.default_rng(seed)
    sig = 10.0 ** rng.uniform(-4, 1, n)
    Q, _ = np.linalg.qr(rng.normal(size=(n, n)))
    C = Q @ np.diag(np.geomspace(1.0, 0.005, n)) @ Q.T
    d = np.sqrt(np.diag(C))
    Hinv = np.linalg.inv(C / np.outer(d, d) * np.outer(sig, sig))
    xs = rng.normal(size=n) * sig
    lo, hi = xs - 30 * sig, xs + 30 * sig
    if edge:                                   # true minimum just outside the box: the fit must end on the edge
        lo[-1] = xs[-1] + 0.5 * sig[-1]
    prop = sig * rng.uniform(0.5, 1.5, n)
    exact = lambda x: 0.5 * (x - xs) @ Hinv @ (x - xs)

    def g(x):
        if fail_rate and rng.uniform() < fail_rate:
            return np.nan
        return exact(x) + noise * rng.normal() + 1200.0

    x0 = np.clip(xs + 3.0 * sig * rng.normal(size=n), lo, hi)
    res = sp_minimize(exact, np.clip(xs, lo, hi), jac=lambda x: Hinv @ (x - xs), bounds=list(zip(lo, hi)),
                      method='L-BFGS-B', options=dict(ftol=1e-15, gtol=1e-12, maxiter=10000))
    return g, x0, lo, hi, prop, exact, res.fun


def excess(n, seed, **kw):
    g, x0, lo, hi, prop, exact, fmin = problem(n, seed, **kw)
    x, fx, info = fitcore.minimize(g, x0, lo, hi, prop, seed=seed)
    assert np.all(x >= lo) and np.all(x <= hi)
    return exact(x) - fmin, info


def test_interior_minimum():
    for seed in range(2):
        e, info = excess(7, seed)
        assert e < 0.15, (seed, e)
        assert info['evaluations'] < 700


def test_minimum_on_prior_edge():
    e, _ = excess(9, 0, edge=True)
    assert e < 0.3


def test_noise_and_failures():
    """5% random hi_class failures (rescued by the jittered retry) on top of the measured noise."""
    for seed in range(2):
        e, info = excess(9, seed, fail_rate=0.05)
        assert e < 0.3, (seed, e)
        assert info['rescued_by_retry'] > 0
