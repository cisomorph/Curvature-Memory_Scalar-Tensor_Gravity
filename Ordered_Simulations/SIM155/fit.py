#!/usr/bin/env python3
"""
SIM155 fit driver (v3): python fit.py MODEL DATA SEED   (MODEL: lcdm | w0wa | cmstgA; DATA: D1 | D2 | D3)

Maximum-posterior fit with common/fitcore.py (tested in tests/test_fitcore.py). Every evaluation goes through
Cobaya's model (theory + likelihoods + priors). From a seeded random start drawn from `ref`:
  1. coarse : Py-BOBYQA in the prior box (≤ 400 evals)
  2. Hessian: central differences at the current point (steps ½ proposal width), regularised in proposal units
  3. fine   : Py-BOBYQA in whitened coordinates (≤ 400 evals); 2–3 repeat until a round gains < 0.05 in −log(post)
Failed hi_class evaluations are retried at a 1e-3-proposal-width jitter; the prior box is a continuous wall.
Progress is written to fits_v3/<run>.json after every stage; `sigma` holds Laplace (inverse-Hessian) errors.

History:
  - v1 (`run_fit.py`, Cobaya `minimize` defaults) stopped 2–5 in χ² above the minimum and crashed on −inf.
  - v2 (BOBYQA rhoend 2e-5 then repeated MIGRAD) was correct but too slow (> 50 min per ΛCDM fit).
  See RESEARCH_LOG.md. Pass criteria are unchanged (SPEC.md).
"""
import json, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'common'))
from cobaya.model import get_model
from cobaya_hiclass import HiClassCMSTG
import fitcore

OUTDIR = os.path.join(HERE, os.environ.get('FIT_OUTDIR', 'fits_v3'))


def build_info(model, data):
    sn = {'D1': 'sn.pantheonplus', 'D2': 'sn.union3', 'D3': 'sn.desy5'}[data]
    like = {'planck_2018_lowl.TT': None, 'planck_2018_lowl.EE': None,
            'planck_2018_highl_plik.TTTEEE_lite_native': None, 'planck_2018_lensing.native': None,
            'bao.desi_dr2.desi_bao_all': None, sn: None}
    P = lambda lo, hi, loc, sc, prop, **kw: dict(prior={'min': lo, 'max': hi},
                                                 ref={'dist': 'norm', 'loc': loc, 'scale': sc}, proposal=prop, **kw)
    params = {
        'omega_b': P(0.019, 0.025, 0.02237, 0.00015, 1.3e-4),
        'omega_cdm': P(0.10, 0.14, 0.1200, 0.0012, 9e-4),
        'H0': P(50., 85., 67.5, 0.6, 0.4),
        'logA': P(2.9, 3.2, 3.044, 0.014, 0.013, drop=True),
        'A_s': {'value': 'lambda logA: 1e-10*np.exp(logA)'},
        'n_s': P(0.9, 1.05, 0.9649, 0.004, 0.0035),
        'tau_reio': P(0.02, 0.10, 0.054, 0.007, 0.007),
    }
    extra = {'non_linear': 'halofit', 'lensing': 'yes'}
    cm = False
    if model == 'w0wa':
        extra.update({'Omega_Lambda': 0.0, 'use_ppf': 'yes'})
        params['w0_fld'] = P(-2.0, 0.0, -0.8, 0.1, 0.06)
        params['wa_fld'] = P(-3.0, 2.0, -0.6, 0.3, 0.25)
    elif model == 'cmstgA':
        cm = True
        params['cmstg_L0'] = dict(prior={'min': 0.02, 'max': 1.0}, ref={'dist': 'uniform', 'loc': 0.10, 'scale': 0.20}, proposal=0.03)
        params['cmstg_xi'] = dict(prior={'min': 0.0, 'max': 0.3}, ref={'dist': 'uniform', 'loc': 0.05, 'scale': 0.07}, proposal=0.01)
    return {'packages_path': '/home/aion/cobaya_packages',
            'theory': {'hiclass': {'external': HiClassCMSTG, 'cmstg_model': cm, 'cmstg_family': 0, 'extra_args': extra}},
            'likelihood': like, 'params': params, 'debug': False}


def main(model, data, seed):
    t0 = time.time()
    run = f'{model}_{data}_s{seed}'
    os.makedirs(OUTDIR, exist_ok=True)
    out = os.path.join(OUTDIR, run + '.json')
    m = get_model(build_info(model, data))
    names = list(m.parameterization.sampled_params())
    lo, hi = m.prior.bounds(confidence_for_unbounded=0.9999995).T
    pinfo = m.parameterization.sampled_params_info()
    prop = np.array([pinfo[n].get('proposal', 0.01 * (b - a)) for n, a, b in zip(names, lo, hi)])
    x0 = m.prior.reference(random_state=np.random.default_rng(seed))
    rec = dict(model=model, data=data, seed=seed, names=names, start=dict(zip(names, map(float, x0))),
               proposal=dict(zip(names, map(float, prop))))

    def save(stage, x, fx, info):
        lls = m.loglikes(x, as_dict=True)[0]
        chi2 = {k: -2.0 * float(v) for k, v in lls.items()}
        rec.update(stage=stage, seconds=time.time() - t0, minimum=dict(zip(names, map(float, x))),
                   minus_log_post=float(fx), chi2_total=float(sum(chi2.values())), chi2=chi2, optimiser=info)
        json.dump(rec, open(out, 'w'), indent=2)
        print(f"[{run}] {stage}: -logpost {fx:.4f}, chi2 {rec['chi2_total']:.3f}, {rec['seconds'] / 60:.1f} min",
              flush=True)

    x, fx, info = fitcore.minimize(lambda v: -m.logpost(v), x0, lo, hi, prop, save=save, seed=seed,
                                   maxfun1=int(os.environ.get('MAXFUN1', 400)),
                                   maxfun3=int(os.environ.get('MAXFUN3', 400)))
    if info.get('sigma') is not None:
        info['sigma'] = dict(zip(names, info['sigma']))
    rec['converged'] = bool(info['rounds'] and info['rounds'][-1]['gain'] < 0.05)
    save('final', x, fx, info)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
