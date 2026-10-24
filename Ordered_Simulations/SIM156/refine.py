#!/usr/bin/env python3
"""SIM156: unconstrained PT best fit in the strong-coupling region found by the grid (log10 Λ0 ≈ 0, v1 ≈ 0.45).
The free fits (run.py fits) started near log10 Λ0 = −1 and converged to the weak-coupling (TQ-like) optimum; the grid
shows a separate region with larger Δχ² that violates LLR. Three extra starts there report the true unconstrained
optimum, as SPEC.md requires ("the unconstrained PT best fit and its α_M(0)"). Writes refine.json."""
import json, os, sys
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run
import fitcore

def one(seed):
    names = run.BASE['names'] + ['log10_L0', 'v1']
    lo = np.array(run.BASE['lo'] + [-3.0, 0.0]); hi = np.array(run.BASE['hi'] + [1.5, 1.0])
    prop = np.array(run.BASE['prop'] + [0.15, 0.03])
    g = json.load(open(os.path.join(HERE, 'grid.json')))
    pt = [r for r in g if r['log10_L0'] is not None and r['chi2'] is not None]
    starts = sorted(pt, key=lambda r: r['chi2']['total'])[:3]           # the three best grid points
    s = starts[seed - 1]
    x0 = np.array([s['x']['omega_b'], s['x']['omega_cdm'], s['x']['h'], s['log10_L0'], s['v1']])
    def f(x):
        c = run.evaluate('pt', x)
        return np.nan if c is None else 0.5 * c['total']
    x, fx, info = fitcore.minimize(f, x0, lo, hi, prop, seed=seed, maxfun1=300, maxfun3=250)
    c = run.evaluate('pt', x)
    return dict(seed=seed, start=dict(zip(names, map(float, x0))), x=dict(zip(names, map(float, x))), chi2=c,
                evaluations=info['evaluations'])

if __name__ == '__main__':
    with Pool(3) as p:
        res = p.map(one, [1, 2, 3])
    json.dump(res, open(os.path.join(HERE, 'refine.json'), 'w'), indent=1, default=float)
    ref = json.load(open(os.path.join(HERE, 'fits.json')))['best']['lcdm']['chi2']['total']
    for r in res:
        c = r['chi2']
        print(f"start {r['start']['log10_L0']:+.1f}/{r['start']['v1']:.2f} -> Δχ² {ref - c['total']:+.3f}  log10Λ0 {r['x']['log10_L0']:+.3f} "
              f"v1 {r['x']['v1']:.3f} h {r['x']['h']:.4f}  α_M0 {c['alpha_M0']:+.4f} F0 {c['F0']:.4f} S {c['S']:.4f}  "
              f"cmb {c['cmb']:.2f} bao {c['bao']:.2f} sn {c['sn']:.2f}")
