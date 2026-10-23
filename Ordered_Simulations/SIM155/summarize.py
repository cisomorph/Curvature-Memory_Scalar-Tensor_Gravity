#!/usr/bin/env python3
"""SIM155 — aggregate the best fits, apply the pre-registered criteria (SPEC.md), write RESULT/output."""
import glob, json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'common'))
import audit_driver as drv

K = {'lcdm': 0, 'w0wa': 2, 'cmstgA': 2}            # extra parameters relative to ΛCDM
rows = [json.load(open(f)) for f in glob.glob(os.path.join(HERE, 'fits', '*_result.json'))]
rows = [r for r in rows if r.get('chi2_total') is not None]
best = {}
for r in rows:
    key = (r['model'], r['data'])
    if key not in best or r['chi2_total'] < best[key]['chi2_total']:
        best[key] = r


def cmstg_checks(mini):
    """Local gravity and BBN at a CMSTG best fit, from hi_class's own background."""
    from hiclassy import HiClass
    from cobaya_hiclass import SMG_SETTINGS, NEUTRINOS, yhe, F_attractor
    L0, xi = mini['cmstg_L0'], mini['cmstg_xi']
    pars = dict(SMG_SETTINGS, **NEUTRINOS)
    pars.update(H0=mini['H0'], omega_b=mini['omega_b'], omega_cdm=mini['omega_cdm'],
                parameters_smg=f"0.7, {L0}, {xi}, 0, 0., 0.", output='',
                YHe=yhe(mini['omega_b'], np.sqrt(2 * F_attractor(L0, xi, 0))), output_background_smg=10)
    c = HiClass(); c.set(pars); c.compute()
    bg = c.get_background()
    i0 = int(np.argmin(np.abs(bg['z'])))
    phi0, dphi0, M2 = bg['phi_smg'][i0], bg["phi'"][i0], bg['M*^2_smg'][i0]
    H0 = bg['H [1/Mpc]'][i0]
    F0 = M2 / 2
    Fp0 = -2 * L0 * phi0 + xi
    y0 = dphi0 / H0                                      # dφ/dN at a = 1 (conformal φ' / (aH))
    lg = drv.local_gravity(F0, [Fp0], Fdot_over_HF=Fp0 * y0 / F0, H0=mini['H0'])
    S = np.sqrt(2 * F0)
    eta10 = 273.9 * mini['omega_b']
    yp = 0.2485 + 0.0016 * ((eta10 - 6) + 100 * (S - 1))
    dh = 2.64 * (6 / (eta10 - 6 * (S - 1)))**1.6
    pyp = (yp - 0.245) / np.hypot(0.003, 0.0006)
    pdh = (dh - 2.547) / np.hypot(0.025, 0.03 * dh)
    return dict(F0=F0, F_attractor=F_attractor(L0, xi, 0), phi0=phi0, Fpsi0=Fp0, local_gravity=lg,
                bbn=dict(S=S, Yp=yp, DH=dh, pull_Yp=pyp, pull_DH=pdh),
                C2=bool(lg['cassini_ok'] and lg['llr_ok']), C3=bool(abs(pyp) <= 2 and abs(pdh) <= 2))


table = {}
print(f"{'data':4} {'model':7} {'χ²':>9} {'Δχ² vs ΛCDM':>12} {'ΔAIC vs ΛCDM':>13} {'Δχ² vs w0wa':>12}  status")
for d in ('D1', 'D2', 'D3'):
    ref = best.get(('lcdm', d))
    w = best.get(('w0wa', d))
    for m in ('lcdm', 'w0wa', 'cmstgA'):
        b = best.get((m, d))
        if b is None:
            continue
        dchi = None if ref is None else ref['chi2_total'] - b['chi2_total']          # > 0: better than ΛCDM
        daic = None if dchi is None else -dchi + 2 * K[m]                              # < 0: preferred
        dw = None if w is None else w['chi2_total'] - b['chi2_total']
        table[(m, d)] = dict(chi2=b['chi2_total'], dchi2_vs_lcdm=dchi, dAIC_vs_lcdm=daic, dchi2_vs_w0wa=dw,
                             status=b['status'], minimum=b['minimum'], chi2_parts=b['chi2'])
        print(f"{d:4} {m:7} {b['chi2_total']:9.2f} {dchi if dchi is not None else float('nan'):+12.2f} "
              f"{daic if daic is not None else float('nan'):+13.2f} {dw if dw is not None else float('nan'):+12.2f}  {b['status']}")

verdict, checks = 'INCOMPLETE', {}
if ('cmstgA', 'D1') in table and ('lcdm', 'D1') in table:
    t1 = table[('cmstgA', 'D1')]
    checks = cmstg_checks(t1['minimum'])
    a = t1['dchi2_vs_lcdm'] >= 6
    c = any(table.get(('cmstgA', d), {}).get('dchi2_vs_lcdm', -1) >= 6 for d in ('D2', 'D3'))
    verdict = 'PASS' if (a and checks['C2'] and checks['C3'] and c) else 'FAIL'
    print(f"\nCMSTG-A best fit (D1): Λ0 = {t1['minimum']['cmstg_L0']:.4f}, ξ = {t1['minimum']['cmstg_xi']:.4f}, "
          f"H0 = {t1['minimum']['H0']:.2f}; F0 = {checks['F0']:.5f} (attractor {checks['F_attractor']:.5f})")
    print(f"  (a) Δχ²(D1) ≥ 6: {a}   (b) C2 {checks['C2']} (|γ−1| = {checks['local_gravity']['gamma_minus_1']:.1e}), "
          f"C3 {checks['C3']} (Y_p pull {checks['bbn']['pull_Yp']:+.2f})   (c) D2/D3: {c}")
    print(f"VERDICT (pre-registered criteria): {verdict}")

json.dump(dict(verdict=verdict, table={f'{m}_{d}': v for (m, d), v in table.items()}, cmstg_checks=checks,
               n_runs=len(rows)), open(os.path.join(HERE, 'output.json'), 'w'), indent=2, default=float)
