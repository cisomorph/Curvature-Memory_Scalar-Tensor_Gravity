#!/usr/bin/env python3
"""
SIM157 runner (SPEC.md). Joint likelihood = background (SIM156: compressed Planck + DESI DR2 + Pantheon+) + eBOSS DR16
RSD-only fσ8 (6 points). The best fit of each model is the joint best fit; background and RSD parts are reported.
    python run.py fits      -> fits.json   (3 seeded starts per model; DPT also from the best grid points)
    python run.py grid      -> grid.json   (DPT (log10 Λ0, v1) grid, profiled over ω_b, ω_c, h)
"""
import json, os, sys, time
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'common'))
import cmstg_obs as obs
import fitcore
import dark as D

RSD_Z = np.array([0.15, 0.38, 0.51, 0.70, 0.85, 1.48])               # Alam et al. 2021, Table III, RSD-only
RSD_FS8 = np.array([0.53, 0.500, 0.455, 0.448, 0.315, 0.462])
RSD_ERR = np.array([0.16, 0.047, 0.039, 0.043, 0.095, 0.045])
LIK = None


def lik():
    global LIK
    if LIK is None:
        LIK = D.L.BgLikelihood()
    return LIK


BASE = dict(names=['omega_b', 'omega_cdm', 'h'], lo=[0.019, 0.10, 0.55], hi=[0.025, 0.14, 0.80],
            prop=[1.5e-4, 1e-3, 0.005], start=[0.02240, 0.1190, 0.680])
EXTRA = {'lcdm': ([], [], [], [], []),
         'w0wa': (['w0', 'wa'], [-2.0, -3.0], [0.0, 2.0], [0.06, 0.25], [-0.85, -0.5]),
         'tq': (['v1'], [0.0], [1.0], [0.05], [0.15]),
         'dpt': (['log10_L0', 'v1'], [-3.0, 0.0], [1.5, 1.0], [0.3, 0.05], [-0.5, 0.25])}


class _Timeout(Exception):
    pass


def _alarm(signum, frame):
    raise _Timeout()


def evaluate(model, x, limit=10):
    """χ² parts at x, or None (failed or unphysical point, or > `limit` s: a safety guard)."""
    import signal
    old = signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(limit)
    try:
        return _evaluate(model, x)
    except _Timeout:
        return None
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)


def _evaluate(model, x):
    ob, oc, h = x[:3]
    info = {}
    if model == 'w0wa':
        H = obs.H_lcdm(h, ob + oc + D.OMEGA_NU, D.OMEGA_R, w0=x[3], wa=x[4])
        d, dp = D.growth_lcdm_like(H, ob, oc, np.r_[RSD_Z, 0.0])
    else:
        L0 = 10.0 ** x[3] if model == 'dpt' else 0.0
        v1 = {'lcdm': 0.0, 'tq': x[3] if model == 'tq' else 0.0, 'dpt': x[4] if model == 'dpt' else 0.0}[model]
        m = D.DarkBG(ob, oc, L0, v1)
        bg = m.solve(h)
        if bg is None:
            return None
        H = m.H_func(bg)
        d, dp = m.growth_table(bg, np.r_[RSD_Z, 0.0])
        info = dict(A0=bg['A0'], phi0=bg['phi0'], y0=bg['y0'], beta_hi=bg['beta_hi'])
        if model == 'dpt':
            zw, w = D.w_eff(bg, ob, oc)
            ww = w[zw < 2.5]                              # ±0.005 tolerance: ΛCDM's w_eff = −1 up to numerical noise
            info.update(w_eff={str(z): float(np.interp(z, zw, w)) for z in (0.0, 0.3, 0.5, 1.0, 1.5, 2.0)},
                        w_eff_min=float(ww.min()), w_eff_max=float(ww.max()),
                        crosses=bool(np.any(ww < -1.005) and np.any(ww > -0.995)))
    c = lik().chi2(H, ob, oc, h)
    fs8 = D.SIGMA8_REF * dp[:-1] / D.ref_delta0()
    s8 = D.SIGMA8_REF * d[-1] / D.ref_delta0()
    om_today = (ob + oc * info.get('A0', 1.0) + D.OMEGA_NU) / h**2
    c.update(rsd=float(np.sum(((fs8 - RSD_FS8) / RSD_ERR) ** 2)), sigma8=float(s8), S8=float(s8 * np.sqrt(om_today / 0.3)),
             fsigma8=dict(zip(map(str, RSD_Z), map(float, fs8))), **info)
    c['background'] = c['total']
    c['total'] = c['background'] + c['rsd']
    return c


def fit(model, seed, x0=None):
    names, lo, hi, prop, start = EXTRA[model]
    names = BASE['names'] + names
    lo = np.array(BASE['lo'] + lo); hi = np.array(BASE['hi'] + hi)
    prop = np.array(BASE['prop'] + prop)
    if x0 is None:
        rng = np.random.default_rng(seed)
        x0 = np.clip(np.array(BASE['start'] + start) + 2 * prop * rng.normal(size=len(prop)), lo, hi)

    def f(x):
        c = evaluate(model, x)
        return np.nan if c is None else 0.5 * c['total']
    try:
        x, fx, info = fitcore.minimize(f, np.asarray(x0, float), lo, hi, prop, seed=seed, maxfun1=300, maxfun3=250)
    except RuntimeError:
        return dict(model=model, seed=seed, start=dict(zip(names, map(float, x0))), x=None, chi2=None, evaluations=None)
    return dict(model=model, seed=seed, start=dict(zip(names, map(float, x0))), x=dict(zip(names, map(float, x))),
                chi2=evaluate(model, x), evaluations=info['evaluations'])


def grid_point(args):
    lg, v1 = args
    lo, hi, prop = np.array(BASE['lo']), np.array(BASE['hi']), np.array(BASE['prop'])

    def f(x):
        c = evaluate('dpt', np.r_[x, [lg, v1]])
        return np.nan if c is None else 0.5 * c['total']
    try:
        x, fx, info = fitcore.minimize(f, np.array(BASE['start']), lo, hi, prop, seed=0, maxfun1=200, maxfun3=150)
    except RuntimeError:                     # every evaluation failed: the field always reaches the singular point of A
        return dict(log10_L0=lg, v1=v1, x=None, chi2=None, note='unphysical for all (ω_b, ω_c, h): A diverges')
    return dict(log10_L0=lg, v1=v1, x=dict(zip(BASE['names'], map(float, x))), chi2=evaluate('dpt', np.r_[x, [lg, v1]]))


def main_grid():
    v1s = [0.0, 0.1, 0.2, 0.3, 0.45, 0.6, 0.8]
    lgs = [-3.0, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5]
    t0 = time.time()
    with Pool(12) as pool:
        res = pool.map(grid_point, [(lg, v) for lg in lgs for v in v1s])
    json.dump(res, open(os.path.join(HERE, 'grid.json'), 'w'), indent=1, default=float)
    print(f'{len(res)} grid points in {(time.time() - t0) / 60:.1f} min')


def main_fits():
    jobs = [(m, s, None) for m in ('lcdm', 'w0wa', 'tq', 'dpt') for s in (1, 2, 3)]
    g = os.path.join(HERE, 'grid.json')
    if os.path.exists(g):                                   # extra DPT starts at the three best grid points
        grid = [r for r in json.load(open(g)) if r['chi2'] is not None]
        for k, r in enumerate(sorted(grid, key=lambda r: r['chi2']['total'])[:3]):
            jobs.append(('dpt', 10 + k, [r['x']['omega_b'], r['x']['omega_cdm'], r['x']['h'], r['log10_L0'], r['v1']]))
    with Pool(min(12, len(jobs))) as pool:
        res = pool.starmap(fit, jobs)
    best = {}
    for r in res:
        if r['chi2'] is not None and (r['model'] not in best or r['chi2']['total'] < best[r['model']]['chi2']['total']):
            best[r['model']] = r
    json.dump(dict(all=res, best=best), open(os.path.join(HERE, 'fits.json'), 'w'), indent=1, default=float)
    ref = best['lcdm']['chi2']
    print(f"{'model':5} {'joint χ²':>9} {'Δjoint':>7} {'Δbg':>7} {'ΔRSD':>6}  {'S8':>6}  params")
    for k, r in best.items():
        c = r['chi2']
        print(f"{k:5} {c['total']:9.3f} {ref['total'] - c['total']:+7.2f} {ref['background'] - c['background']:+7.2f} "
              f"{ref['rsd'] - c['rsd']:+6.2f}  {c['S8']:.4f}  " + ' '.join(f'{n}={v:.4g}' for n, v in r['x'].items())
              + (f"  A0={c['A0']:.4f} β_hi={c['beta_hi']:.1e} w_eff(0,0.5,1)=" + ','.join(f"{c['w_eff'][z]:.3f}" for z in ('0.0', '0.5', '1.0'))
                 + f" crosses={c['crosses']}" if k == 'dpt' else ''))


if __name__ == '__main__':
    os.environ.setdefault('OMP_NUM_THREADS', '1')
    {'fits': main_fits, 'grid': main_grid}[sys.argv[1]]()
