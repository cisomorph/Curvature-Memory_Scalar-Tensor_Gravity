#!/usr/bin/env python3
"""
SIM156 runner (SPEC.md): validation V3, then best fits of ΛCDM, w0wa, TQ and PT, plus the PT (Λ0, v1) grid.
    python run.py fits        # V3 and the best fits (3 seeds each)  -> fits.json
    python run.py grid        # PT grid with profiled (ω_b, ω_cdm, h) -> grid.json
"""
import json, os, sys, time
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'common'))
import cmstg_obs as obs
import fitcore
from st_solver import ScalarTensorBG
import bglik as L

LIK = None


def lik():
    global LIK
    if LIK is None:
        LIK = L.BgLikelihood()
    return LIK


BASE = dict(names=['omega_b', 'omega_cdm', 'h'], lo=[0.019, 0.10, 0.55], hi=[0.025, 0.14, 0.80],
            prop=[1.5e-4, 1e-3, 0.005], start=[0.02240, 0.1190, 0.680])
EXTRA = {'lcdm': ([], [], [], [], []),
         'w0wa': (['w0', 'wa'], [-2.0, -3.0], [0.0, 2.0], [0.06, 0.25], [-0.85, -0.5]),
         'tq': (['v1'], [0.0], [1.0], [0.05], [0.15]),
         'pt': (['log10_L0', 'v1'], [-3.0, 0.0], [1.5, 1.0], [0.3, 0.05], [-1.0, 0.15])}


def scalar_bg(om, h, L0, v1):
    """TQ (L0 = 0) or PT background: F = ½ − L0 δΨ², V = 3 v1 δΨ, δΨ = δΨ' = 0 at z = 3000."""
    m = ScalarTensorBG(F=lambda p, N: 0.5 - L0 * p * p, F_psi=lambda p, N: -2.0 * L0 * p,
                       V=lambda p: 3.0 * v1 * p, V_psi=lambda p: 3.0 * v1,
                       omega_m=om, omega_r=obs.OMEGA_R, z_ini=3000.0, n_points=1500, max_step=0.05)
    bg = m.calibrated(0.0, 0.0, h)
    if bg is None:
        return None, None
    return m.H_func(bg), bg


def evaluate(model, x):
    """Returns (chi2 parts dict incl. alpha_M0, F0, H) or None."""
    ob, oc, h = x[:3]
    om = ob + oc + L.OMEGA_NU
    info = {}
    if model == 'lcdm':
        H = obs.H_lcdm(h, om, obs.OMEGA_R)
    elif model == 'w0wa':
        H = obs.H_lcdm(h, om, obs.OMEGA_R, w0=x[3], wa=x[4])
    else:
        L0 = 0.0 if model == 'tq' else 10.0 ** x[3]
        v1 = x[3] if model == 'tq' else x[4]
        H, bg = scalar_bg(om, h, L0, v1)
        if H is None:
            return None
        F0, p0, y0 = bg['F0'], bg['psi0'], bg['y0']
        info = dict(F0=F0, dpsi0=p0, y0=y0, alpha_M0=float(-2.0 * L0 * p0 * y0 / F0),
                    S=float(np.sqrt(F0 / 0.5)))
    c = lik().chi2(H, ob, oc, h)
    c.update(info)
    return c


def fit(model, seed, llr=False):
    names, lo, hi, prop, start = EXTRA[model]
    names = BASE['names'] + names
    lo = np.array(BASE['lo'] + lo); hi = np.array(BASE['hi'] + hi)
    prop = np.array(BASE['prop'] + prop); x0 = np.array(BASE['start'] + start)
    rng = np.random.default_rng(seed)
    x0 = np.clip(x0 + 2 * prop * rng.normal(size=len(x0)), lo, hi)

    def f(x):
        c = evaluate(model, x)
        if c is None:
            return np.nan
        v = 0.5 * c['total']
        if llr and 'alpha_M0' in c and abs(c['alpha_M0']) > 0.002:      # LLR as a stiff, continuous wall
            v += 0.5 * ((abs(c['alpha_M0']) - 0.002) / 2e-4) ** 2
        return v

    x, fx, info = fitcore.minimize(f, x0, lo, hi, prop, seed=seed, maxfun1=300, maxfun3=200)
    c = evaluate(model, x)
    return dict(model=model, seed=seed, llr_constrained=llr, x=dict(zip(names, map(float, x))), chi2=c,
                evaluations=info['evaluations'])


def main_fits():
    jobs = [(m, s, False) for m in ('lcdm', 'w0wa', 'tq', 'pt') for s in (1, 2, 3)] + [('pt', s, True) for s in (1, 2, 3)]
    with Pool(min(len(jobs), 12)) as pool:
        res = pool.starmap(fit, jobs)
    best = {}
    for r in res:
        k = r['model'] + ('_llr' if r['llr_constrained'] else '')
        if k not in best or r['chi2']['total'] < best[k]['chi2']['total']:
            best[k] = r
    json.dump(dict(all=res, best=best), open(os.path.join(HERE, 'fits.json'), 'w'), indent=1, default=float)
    ref = best['lcdm']['chi2']['total']
    print(f"{'model':8} {'χ²':>9} {'Δχ² vs ΛCDM':>11}  cmb    bao     sn       params")
    for k, r in best.items():
        c = r['chi2']
        extra = f" α_M0 {c['alpha_M0']:+.4f} F0 {c['F0']:.5f}" if 'alpha_M0' in c else ''
        print(f"{k:8} {c['total']:9.3f} {ref - c['total']:+11.3f}  {c['cmb']:.2f}  {c['bao']:.2f}  {c['sn']:.2f}  "
              + ' '.join(f'{n}={v:.4g}' for n, v in r['x'].items()) + extra)
    d = ref - best['w0wa']['chi2']['total']
    print(f"\nV3: Δχ²(w0wa vs ΛCDM) = {d:+.2f}; required in [6.2, 12.2] (SIM155 full likelihood: +9.19) -> "
          f"{'PASS' if 6.2 <= d <= 12.2 else 'FAIL'}")


def grid_point(args):
    lg, v1 = args
    names_fixed = (lg, v1)

    def f(x):
        c = evaluate('pt' if lg is not None else 'tq', np.r_[x, [lg, v1]] if lg is not None else np.r_[x, [v1]])
        return np.nan if c is None else 0.5 * c['total']
    lo, hi, prop = np.array(BASE['lo']), np.array(BASE['hi']), np.array(BASE['prop'])
    x, fx, info = fitcore.minimize(f, np.array(BASE['start']), lo, hi, prop, seed=0, maxfun1=200, maxfun3=150)
    model = 'pt' if lg is not None else 'tq'
    c = evaluate(model, np.r_[x, [lg, v1]] if lg is not None else np.r_[x, [v1]])
    return dict(log10_L0=lg, v1=v1, x=dict(zip(BASE['names'], map(float, x))), chi2=c)


def main_grid():
    v1s = [0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.45, 0.6]
    lgs = [None, -3.0, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0]          # None = TQ (Λ0 = 0)
    jobs = [(lg, v) for lg in lgs for v in v1s]
    t0 = time.time()
    with Pool(12) as pool:
        res = pool.map(grid_point, jobs)
    json.dump(res, open(os.path.join(HERE, 'grid.json'), 'w'), indent=1, default=float)
    print(f'{len(res)} grid points in {(time.time() - t0) / 60:.1f} min')


if __name__ == '__main__':
    os.environ.setdefault('OMP_NUM_THREADS', '1')
    {'fits': main_fits, 'grid': main_grid}[sys.argv[1]]()
