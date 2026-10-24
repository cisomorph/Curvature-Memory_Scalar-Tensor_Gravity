#!/usr/bin/env python3
"""
SIM155 follow-up (ROADMAP 5.2 and 5.3 at the full-CMB best fits): G(z) history and growth of structure.

For the D1 best fits of ΛCDM and w0waCDM, the CMSTG-A best fit, and the best CMSTG-A profile point allowed by BBN
(C3), hi_class computes:
  - σ8, S8 = σ8 √(Ωm/0.3), and fσ8(z) = dσ8/d ln a (central differences of σ8(z) from the linear P(k, z), which
    includes the scalar field's effect on growth);
  - for CMSTG: the G history G(z)/G0 = F0/F(z) from hi_class's M*² = 2F (the local G today is G0; the attractor
    keeps F' ≈ 0 today, so the Brans–Dicke correction to G_N is < 1e-7).
Reads output.json (summarize.py); writes growth.json and prints a table.
"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'common'))
from hiclassy import HiClass
from cobaya_hiclass import SMG_SETTINGS, NEUTRINOS, yhe, F_attractor

out = json.load(open(os.path.join(HERE, 'output.json')))
ZS = [0.0, 0.3, 0.5, 0.8, 1.0, 1.5]


def run(kind, mn):
    p = dict(NEUTRINOS, H0=mn['H0'], omega_b=mn['omega_b'], omega_cdm=mn['omega_cdm'], A_s=1e-10 * np.exp(mn['logA']),
             n_s=mn['n_s'], tau_reio=mn['tau_reio'], output='mPk', **{'P_k_max_1/Mpc': 3.0, 'z_max_pk': 2.0})
    S = 1.0
    if kind == 'cmstgA':
        p.update(SMG_SETTINGS)
        p['parameters_smg'] = f"0.7, {mn['cmstg_L0']:.12g}, {mn['cmstg_xi']:.12g}, 0, 0., 0."
        p['output_background_smg'] = 10
        S = float(np.sqrt(2 * F_attractor(mn['cmstg_L0'], mn['cmstg_xi'], 0)))
    elif kind == 'w0wa':
        p.update(Omega_Lambda=0.0, use_ppf='yes', w0_fld=mn['w0_fld'], wa_fld=mn['wa_fld'])
    p['YHe'] = yhe(mn['omega_b'], S)
    c = HiClass(); c.set(p); c.compute()
    h = c.h()
    s8 = lambda z: c.sigma(8.0 / h, z)
    res = dict(sigma8=float(s8(0.0)), Omega_m=float(c.Omega_m()))
    res['S8'] = res['sigma8'] * np.sqrt(res['Omega_m'] / 0.3)
    dz = 0.01
    res['fsigma8'] = {}
    for z in ZS:
        za, zb = max(z - dz, 0.0), z + dz
        res['fsigma8'][str(z)] = float(-(1 + z) * (s8(zb) - s8(za)) / (zb - za))
    if kind == 'cmstgA':
        bg = c.get_background()
        z, F = bg['z'], bg['M*^2_smg'] / 2
        F0 = float(F[np.argmin(np.abs(z))])
        res['G_over_G0'] = {str(zz): float(F0 / np.interp(zz, z[::-1], F[::-1])) for zz in (0.0, 0.5, 2.0, 10.0, 1090.0, 1e6)}
        res['G_BBN_over_G0'] = float(F0 / 0.5)                    # φ = 0, F = ½ before the field moves
        res['F0'] = F0
    c.struct_cleanup()
    return res


cases = {'ΛCDM best fit': ('lcdm', out['table']['lcdm_D1']['minimum']),
         'w0waCDM best fit': ('w0wa', out['table']['w0wa_D1']['minimum']),
         'CMSTG-A best fit': ('cmstgA', out['table']['cmstgA_D1']['minimum'])}
bb = out.get('bbn_allowed_best')
if bb:
    cases[f"CMSTG-A best BBN-allowed (ΔF = {bb['dF']:.4f})"] = ('cmstgA', bb['minimum'])

res = {name: dict(kind=k, minimum=mn, **run(k, mn)) for name, (k, mn) in cases.items()}
ref = res['ΛCDM best fit']
print(f"{'case':38} {'σ8':>7} {'S8':>7} {'Ωm':>7}  " + ' '.join(f'fσ8({z})' for z in ZS))
for name, r in res.items():
    print(f"{name:38} {r['sigma8']:7.4f} {r['S8']:7.4f} {r['Omega_m']:7.4f}  " +
          ' '.join(f"{r['fsigma8'][str(z)]:8.4f}" for z in ZS))
for name, r in res.items():
    if name != 'ΛCDM best fit':
        print(f"{name:38} vs ΛCDM: σ8 {100 * (r['sigma8'] / ref['sigma8'] - 1):+.2f}%, S8 {100 * (r['S8'] / ref['S8'] - 1):+.2f}%, "
              + ', '.join(f"fσ8({z}) {100 * (r['fsigma8'][str(z)] / ref['fsigma8'][str(z)] - 1):+.2f}%" for z in (0.0, 0.5, 1.0)))
    if 'G_over_G0' in r:
        print(f"   G/G0: BBN {r['G_BBN_over_G0']:.4f}; " + ', '.join(f"z={k}: {v:.4f}" for k, v in r['G_over_G0'].items()))
json.dump(res, open(os.path.join(HERE, 'growth.json'), 'w'), indent=1, default=float)
