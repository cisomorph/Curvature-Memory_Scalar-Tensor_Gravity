#!/usr/bin/env python3
"""
SIM155 systematic check: χ² offset of the hi_class CMSTG model in its exact GR limit (ξ = 0, so φ ≡ 0 and F ≡ ½)
against plain ΛCDM at the same cosmological parameters, on the full D1 likelihood. hi_class_model/README.md
records a 2e-4 … 9e-4 (of peak) C_ℓ residual in this limit from hi_class's scalar-sector numerics; this turns it into
a Δχ², the size of any artificial CMSTG-vs-ΛCDM offset. ξ = 1e-4 checks continuity. Writes gr_limit_chi2.json.
"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'common'))
from cobaya.model import get_model
from fit import build_info

POINT = dict(omega_b=0.02237, omega_cdm=0.1200, H0=67.5, logA=3.044, n_s=0.9649, tau_reio=0.054, A_planck=1.0)


def chi2_parts(m, extra=None):
    names = list(m.parameterization.sampled_params())
    p = dict(POINT, **(extra or {}))
    x = np.array([p[n] for n in names])
    lls = m.loglikes(x, as_dict=True)[0]
    return {k: -2.0 * float(v) for k, v in lls.items()}


out = {'point': POINT}
ml = get_model(build_info('lcdm', 'D1'))
out['lcdm'] = ref = chi2_parts(ml)
print('ΛCDM', round(sum(ref.values()), 3), flush=True)
mc = get_model(build_info('cmstgA', 'D1'))
for L0, xi in ((0.15, 0.0), (0.5, 0.0), (0.15, 1e-4)):
    c = chi2_parts(mc, dict(cmstg_L0=L0, cmstg_xi=xi))
    d = {k: c[k] - ref[k] for k in ref}
    out[f'cmstg_L0={L0}_xi={xi}'] = dict(chi2=c, delta_vs_lcdm=d, delta_total=float(sum(d.values())))
    print(f'CMSTG L0={L0} xi={xi}: Δχ² vs ΛCDM = {sum(d.values()):+.4f}', {k: round(v, 4) for k, v in d.items()}, flush=True)
# the same comparison at three more cosmologies (±1σ moves), to separate a systematic offset from evaluation noise
moves = [dict(H0=68.1, omega_cdm=0.1188), dict(H0=66.9, n_s=0.9689), dict(omega_b=0.02252, tau_reio=0.060, logA=3.056)]
dl = [out['cmstg_L0=0.15_xi=0.0']['delta_total']]
for mv in moves:
    base = dict(POINT, **mv)
    names_l = list(ml.parameterization.sampled_params()); names_c = list(mc.parameterization.sampled_params())
    cl = {k: -2.0 * float(v) for k, v in ml.loglikes(np.array([base[n] for n in names_l]), as_dict=True)[0].items()}
    bc = dict(base, cmstg_L0=0.15, cmstg_xi=0.0)
    cc = {k: -2.0 * float(v) for k, v in mc.loglikes(np.array([bc[n] for n in names_c]), as_dict=True)[0].items()}
    dl.append(float(sum(cc.values()) - sum(cl.values())))
    print('move', mv, f'Δχ² (GR-limit CMSTG − ΛCDM) = {dl[-1]:+.4f}', flush=True)
out['multi_point'] = dict(moves=moves, delta_total=dl, mean=float(np.mean(dl)), std=float(np.std(dl, ddof=1)))
print(f'mean offset {np.mean(dl):+.3f} ± {np.std(dl, ddof=1) / np.sqrt(len(dl)):.3f} (sd {np.std(dl, ddof=1):.3f}, n = {len(dl)})')
json.dump(out, open(os.path.join(HERE, 'gr_limit_chi2.json'), 'w'), indent=1)
