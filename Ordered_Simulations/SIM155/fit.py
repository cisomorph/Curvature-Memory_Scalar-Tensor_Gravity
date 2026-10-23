#!/usr/bin/env python3
"""
SIM155 fit driver (v2): python fit.py MODEL DATA SEED   (MODEL: lcdm | w0wa | cmstgA; DATA: D1 | D2 | D3)

Two-stage maximum-posterior fit, using Cobaya's model (theory + likelihoods + priors) for every evaluation:
  1. BOBYQA (Py-BOBYQA), scaled to the prior box, from a seeded random start drawn from `ref`.
     Failed evaluations (hi_class errors) return a large finite penalty, not +inf, so the solver keeps working.
  2. MIGRAD (iminuit) from the BOBYQA point, restarted until two passes agree to 0.01 in −log(post).
v1 (`run_fit.py`, Cobaya's `minimize`) stopped up to 2–5 χ² above the minimum and crashed on −inf;
see ROADMAP log 2026-10-23. Criteria are unchanged (SPEC.md).
"""
import json, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'common'))
from cobaya.model import get_model
from cobaya_hiclass import HiClassCMSTG
import pybobyqa
from iminuit import Minuit

PENALTY = 1e6


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
    m = get_model(build_info(model, data))
    names = list(m.parameterization.sampled_params())
    lo, hi = m.prior.bounds(confidence_for_unbounded=0.9999995).T
    info_p = m.parameterization.sampled_params_info()
    steps = np.array([info_p[n].get('proposal', 0.01 * (h - l)) for n, l, h in zip(names, lo, hi)])
    nev = {'n': 0, 'fail': 0}

    def f(x):
        nev['n'] += 1
        x = np.clip(x, lo, hi)
        lp = m.logpost(x)
        if not np.isfinite(lp):
            nev['fail'] += 1
            return PENALTY
        return -lp

    rng = np.random.default_rng(seed)
    x0 = m.prior.reference(random_state=rng)
    f0 = f(x0)
    # stage 1: BOBYQA in the scaled box
    sol = pybobyqa.solve(f, x0, bounds=(lo, hi), scaling_within_bounds=True, rhobeg=0.05, rhoend=2e-5,
                         maxfun=int(os.environ.get('MAXFUN', 2500)), seek_global_minimum=False, do_logging=False,
                         user_params={'slow.max_slow_iters': 60})
    x1, f1 = sol.x, sol.f
    # stage 2: MIGRAD polish, repeated until two passes agree
    passes = []
    x, fx = x1.copy(), f1
    for k in range(4):
        mi = Minuit(lambda v: f(np.asarray(v)), x, name=names)
        mi.errordef = Minuit.LIKELIHOOD
        mi.limits = list(zip(lo, hi))
        mi.errors = steps
        mi.strategy = 1
        mi.migrad(ncall=int(os.environ.get('NCALL', 1500)))
        xn, fn = np.array(mi.values), mi.fval
        passes.append(dict(fval=float(fn), valid=bool(mi.valid), edm=float(mi.fmin.edm)))
        improved = fx - fn
        if fn < fx:
            x, fx = xn, fn
        if improved < 0.01:
            break
    lls = m.loglikes(x, as_dict=True)[0]
    chi2 = {k: -2.0 * float(v) for k, v in lls.items()}
    res = dict(model=model, data=data, seed=seed, seconds=time.time() - t0, evaluations=nev['n'], failed_evaluations=nev['fail'],
               start=dict(zip(names, map(float, x0))), f_start=float(f0), f_bobyqa=float(f1), bobyqa_flag=int(sol.flag),
               migrad_passes=passes, converged=bool(len(passes) and passes[-1]['valid'] and (len(passes) < 4)),
               minimum=dict(zip(names, map(float, x))), minus_log_post=float(fx),
               chi2_total=float(sum(chi2.values())), chi2=chi2)
    out = os.path.join(HERE, 'fits_v2', f'{model}_{data}_s{seed}.json')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    json.dump(res, open(out, 'w'), indent=2)
    print(json.dumps({k: res[k] for k in ('model', 'data', 'seed', 'seconds', 'evaluations', 'failed_evaluations',
                                          'f_bobyqa', 'minus_log_post', 'chi2_total', 'converged')}, default=float))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
