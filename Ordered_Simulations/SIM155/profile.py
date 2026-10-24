#!/usr/bin/env python3
"""
SIM155 profile likelihood of CMSTG-A in the attractor strength ΔF = F* − ½ = ξ²/(4Λ0):
    python profile.py DATA DF        (e.g. python profile.py D1 0.012)

Why: the free CMSTG fits (fit.py) end in a long, curved, nearly flat valley in (Λ0, ξ, H0). Along it the attractor
field φ* = ξ/(2Λ0) stays ≈ 0.28–0.35 while ΔF and H0 change, and χ² changes by < 1. Local optimisers stop at different
points of such a valley (seed spreads > 0.5 in χ²). Fixing ΔF and minimising over everything else is the standard
remedy: min over the ΔF grid of the profile is the best fit, and the profile shows directly how much improvement
survives where BBN allows it (C3 needs ΔF ≲ 0.012). This is a method for finding the best fit; data, model, priors
and the pre-registered criteria are unchanged (the pre-registered box ξ ≤ 0.3 is kept: Λ0 ≤ 0.0225/ΔF).

At fixed ΔF the sampled CMSTG parameter is Λ0 (flat prior); ξ = 2√(Λ0 ΔF) is computed. Starts: the six ΛCDM
parameters from the D1 ΛCDM best fit, Λ0 on the valley guess ΔF/φ*(ΔF)².
Output: profile/cmstgA_<DATA>_dF<DF>.json (fit.py format, plus 'dF').
"""
import json, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'common'))
from cobaya.model import get_model
import fitcore
from fit import build_info

OUTDIR = os.path.join(HERE, 'profile')
def phi_guess(dF):
    """Valley guess for φ* = ξ/(2Λ0) from the free D1 fits: (ΔF, φ*) = (0.009, 0.29), (0.011, 0.28), (0.024, 0.30),
    (0.035, 0.35), i.e. ≈ 0.28 up to ΔF ≈ 0.018 and rising beyond. Only the start depends on it."""
    return 0.28 + 4.0 * max(0.0, dF - 0.018)


def setup(data, dF):
    """Cobaya model at fixed ΔF, its sampled names, box, proposals and start point."""
    info = build_info('cmstgA', data)
    L_hi = min(1.0, 0.0225 / dF)                       # keeps ξ = 2√(Λ0 ΔF) ≤ 0.3 (pre-registered box)
    L_ref = float(np.clip(dF / phi_guess(dF)**2, 0.03, 0.9 * L_hi))
    info['params']['cmstg_L0'] = dict(prior={'min': 0.02, 'max': L_hi},
                                      ref={'dist': 'norm', 'loc': L_ref, 'scale': 0.1 * L_ref}, proposal=0.15 * L_ref)
    info['params']['cmstg_xi'] = {'value': f'lambda cmstg_L0: 2*np.sqrt(cmstg_L0*{dF!r})'}
    m = get_model(info)
    names = list(m.parameterization.sampled_params())
    lo, hi = m.prior.bounds(confidence_for_unbounded=0.9999995).T
    pinfo = m.parameterization.sampled_params_info()
    prop = np.array([pinfo[n].get('proposal', 0.01 * (b - a)) for n, a, b in zip(names, lo, hi)])
    # start: ΛCDM D1 best fit for the shared parameters, the valley guess for Λ0
    ref = min((json.load(open(os.path.join(HERE, 'fits_v3', f'lcdm_D1_s{s}.json'))) for s in (1, 2, 3)),
              key=lambda r: r['minus_log_post'])['minimum']
    x0 = np.array([ref.get(n, L_ref) for n in names])
    return m, names, lo, hi, prop, x0


def main(data, dF):
    t0 = time.time()
    run = f'cmstgA_{data}_dF{dF:.4f}'
    os.makedirs(OUTDIR, exist_ok=True)
    out = os.path.join(OUTDIR, run + '.json')
    m, names, lo, hi, prop, x0 = setup(data, dF)
    rec = dict(model='cmstgA', data=data, dF=dF, phi_star_guess=phi_guess(dF), names=names,
               start=dict(zip(names, map(float, x0))), proposal=dict(zip(names, map(float, prop))))

    def save(stage, x, fx, oinfo):
        lls = m.loglikes(x, as_dict=True)[0]
        chi2 = {k: -2.0 * float(v) for k, v in lls.items()}
        mini = dict(zip(names, map(float, x)))
        mini['cmstg_xi'] = float(2 * np.sqrt(mini['cmstg_L0'] * dF))
        rec.update(stage=stage, seconds=time.time() - t0, minimum=mini, minus_log_post=float(fx),
                   chi2_total=float(sum(chi2.values())), chi2=chi2, optimiser=oinfo)
        json.dump(rec, open(out, 'w'), indent=2)
        print(f"[{run}] {stage}: -logpost {fx:.4f}, chi2 {rec['chi2_total']:.3f}, L0 {mini['cmstg_L0']:.4f}, "
              f"H0 {mini['H0']:.2f}, {rec['seconds'] / 60:.1f} min", flush=True)

    x, fx, oinfo = fitcore.minimize(lambda v: -m.logpost(v), x0, lo, hi, prop, save=save, seed=0,
                                    maxfun1=int(os.environ.get('MAXFUN1', 300)))
    if oinfo.get('sigma') is not None:
        oinfo['sigma'] = dict(zip(names, oinfo['sigma']))
    rec['converged'] = bool(oinfo['rounds'] and oinfo['rounds'][-1]['gain'] < 0.1)
    save('final', x, fx, oinfo)


if __name__ == '__main__':
    main(sys.argv[1], float(sys.argv[2]))
