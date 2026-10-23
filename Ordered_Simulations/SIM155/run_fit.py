#!/usr/bin/env python3
"""SIM155 — one best-fit run: python run_fit.py MODEL DATA SEED  (MODEL: lcdm | w0wa | cmstgA; DATA: D1 | D2 | D3)."""
import json, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'common'))
from cobaya.run import run
from cobaya_hiclass import HiClassCMSTG

model, data, seed = sys.argv[1], sys.argv[2], int(sys.argv[3])
SN = {'D1': 'sn.pantheonplus', 'D2': 'sn.union3', 'D3': 'sn.desy5'}[data]
like = {'planck_2018_lowl.TT': None, 'planck_2018_lowl.EE': None,
        'planck_2018_highl_plik.TTTEEE_lite_native': None, 'planck_2018_lensing.native': None,
        'bao.desi_dr2.desi_bao_all': None, SN: None}
params = {
    'omega_b':   {'prior': {'min': 0.019, 'max': 0.025}, 'ref': {'dist': 'norm', 'loc': 0.02237, 'scale': 0.00015}, 'proposal': 1e-4},
    'omega_cdm': {'prior': {'min': 0.10, 'max': 0.14}, 'ref': {'dist': 'norm', 'loc': 0.1200, 'scale': 0.0012}, 'proposal': 1e-3},
    'H0':        {'prior': {'min': 50., 'max': 85.}, 'ref': {'dist': 'norm', 'loc': 67.5, 'scale': 0.6}, 'proposal': 0.5},
    'logA':      {'prior': {'min': 2.9, 'max': 3.2}, 'ref': {'dist': 'norm', 'loc': 3.044, 'scale': 0.014}, 'proposal': 0.01, 'drop': True},
    'A_s':       {'value': 'lambda logA: 1e-10*np.exp(logA)'},
    'n_s':       {'prior': {'min': 0.9, 'max': 1.05}, 'ref': {'dist': 'norm', 'loc': 0.9649, 'scale': 0.004}, 'proposal': 0.003},
    'tau_reio':  {'prior': {'min': 0.02, 'max': 0.10}, 'ref': {'dist': 'norm', 'loc': 0.054, 'scale': 0.007}, 'proposal': 0.005},
}
extra = {'non_linear': 'halofit', 'lensing': 'yes'}
cm = False
if model == 'w0wa':
    extra.update({'Omega_Lambda': 0.0, 'use_ppf': 'yes'})
    params['w0_fld'] = {'prior': {'min': -2.0, 'max': 0.0}, 'ref': {'dist': 'norm', 'loc': -0.8, 'scale': 0.1}, 'proposal': 0.05}
    params['wa_fld'] = {'prior': {'min': -3.0, 'max': 2.0}, 'ref': {'dist': 'norm', 'loc': -0.6, 'scale': 0.3}, 'proposal': 0.1}
elif model == 'cmstgA':
    cm = True
    params['cmstg_L0'] = {'prior': {'min': 0.02, 'max': 1.0}, 'ref': {'dist': 'uniform', 'loc': 0.10, 'scale': 0.20}, 'proposal': 0.02}
    params['cmstg_xi'] = {'prior': {'min': 0.0, 'max': 0.3}, 'ref': {'dist': 'uniform', 'loc': 0.05, 'scale': 0.07}, 'proposal': 0.01}

out = os.path.join(HERE, 'fits', f'{model}_{data}_s{seed}')
info = {'packages_path': '/home/aion/cobaya_packages',
        'theory': {'hiclass': {'external': HiClassCMSTG, 'cmstg_model': cm, 'cmstg_family': 0, 'extra_args': extra}},
        'likelihood': like, 'params': params,
        'sampler': {'minimize': {'method': 'bobyqa', 'seed': seed, 'max_evals': int(os.environ.get('MAX_EVALS', 3000))}},
        'output': out, 'force': True}
t0 = time.time()
try:
    run(info)
    status = 'converged'
except Exception as e:          # e.g. max_evals reached; the last point is still informative
    status = f'not converged: {type(e).__name__}'


def read_minimum(path):
    """Parse Cobaya's <out>.minimum.txt (GetDist format: '#' header with names, one data row)."""
    if not os.path.exists(path):
        return None
    lines = open(path).read().splitlines()
    names = lines[0].lstrip('#').split()
    vals = [float(x) for x in lines[1].split()]
    return dict(zip(names, vals))


mini = read_minimum(out + '.minimum.txt')
res = dict(model=model, data=data, seed=seed, seconds=time.time() - t0, status=status,
           minimum=mini,
           chi2_total=None if mini is None else mini.get('chi2'),
           chi2={k: v for k, v in (mini or {}).items() if k.startswith('chi2__')})
json.dump(res, open(out + '_result.json', 'w'), indent=2, default=float)
print(json.dumps(res, indent=1, default=float))
