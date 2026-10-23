#!/usr/bin/env python3
"""
SIM155 optimiser diagnostics on the real D1 likelihood (Planck 2018 + DESI DR2 + Pantheon+). Writes diagnostics.json.

  noise : −log(post) along tiny steps (±0.02 proposal widths, 11 points) in H0 and omega_cdm around a ΛCDM
          point; the RMS residual of a quadratic fit is the evaluation noise that the optimiser must tolerate.
  fails : fraction of CMSTG-A points drawn from the fit's `ref` distribution that hi_class rejects, and whether
          a 1e-3-proposal-width jitter (fitcore's retry rule) rescues them.
Usage: python diagnostics.py noise | fails WORKER NWORKERS NPOINTS
"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'common'))
from cobaya.model import get_model
from fit import build_info


def prop_of(m, names):
    pinfo = m.parameterization.sampled_params_info()
    return np.array([pinfo[n]['proposal'] for n in names])


def noise():
    m = get_model(build_info('lcdm', 'D1'))
    names = list(m.parameterization.sampled_params())
    x0 = np.array([m.parameterization.sampled_params_info()[n]['ref']['loc'] if isinstance(
        m.parameterization.sampled_params_info()[n]['ref'], dict) else 1.0 for n in names])
    prop = prop_of(m, names)
    out = {}
    for p in ('H0', 'omega_cdm'):
        i = names.index(p)
        t = np.linspace(-0.02, 0.02, 11)
        f, parts = [], []
        for ti in t:
            x = x0.copy(); x[i] += ti * prop[i]
            f.append(-m.logpost(x))
            parts.append({k: -2 * float(v) for k, v in m.loglikes(x, as_dict=True)[0].items()})
        f = np.array(f)
        c = np.polyfit(t, f, 2)
        res = f - np.polyval(c, t)
        per = {}
        for k in parts[0]:                     # noise of each likelihood's χ² (quadratic-fit residual RMS)
            y = np.array([q[k] for q in parts])
            per[k] = float((y - np.polyval(np.polyfit(t, y, 2), t)).std(ddof=3))
        out[p] = dict(step_in_proposals=t.tolist(), f=f.tolist(), rms_residual=float(res.std(ddof=3)),
                      slope_per_proposal=float(c[1]), chi2_noise_by_likelihood=per)
        print(p, 'noise RMS in −log(post):', out[p]['rms_residual'], ' slope:', c[1], flush=True)
        print('   χ² noise by likelihood:', {k: round(v, 4) for k, v in per.items()}, flush=True)
    return out


def fails(worker, nworkers, npoints):
    m = get_model(build_info('cmstgA', 'D1'))
    names = list(m.parameterization.sampled_params())
    prop = prop_of(m, names)
    rng = np.random.default_rng(500 + worker)
    rows = []
    for k in range(npoints):
        x = m.prior.reference(random_state=rng)
        lp = m.logpost(x)
        row = dict(point=dict(zip(names, map(float, x))), ok=bool(np.isfinite(lp)), f=float(-lp) if np.isfinite(lp) else None)
        if not row['ok']:
            res = []
            for j in range(3):
                xj = x + 1e-3 * prop * rng.normal(size=len(x))
                lj = m.logpost(xj)
                res.append(bool(np.isfinite(lj)))
                if res[-1]:
                    row['rescued_f'] = float(-lj)
                    break
            row['retries_ok'] = res
        rows.append(row)
        print(worker, k, row['ok'], row.get('retries_ok'), flush=True)
    return rows


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'noise':
        json.dump(noise(), open(os.path.join(HERE, 'diag_noise.json'), 'w'), indent=1)
    else:
        w, nw, npnt = map(int, sys.argv[2:5])
        json.dump(fails(w, nw, npnt), open(os.path.join(HERE, f'diag_fails_{w}.json'), 'w'), indent=1)
