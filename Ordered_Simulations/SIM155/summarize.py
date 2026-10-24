#!/usr/bin/env python3
"""SIM155 — aggregate the v3 best fits (fits_v3/*.json), apply the pre-registered criteria (SPEC.md), write
output.json and results.md. Runs on partial output too (unfinished fits are listed as 'running')."""
import glob, json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'common'))
import audit_driver as drv

K = {'lcdm': 0, 'w0wa': 2, 'cmstgA': 2}            # extra parameters relative to ΛCDM
FITDIR = os.path.join(HERE, os.environ.get('FIT_OUTDIR', 'fits_v3'))
runs = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(FITDIR, '*.json')))]
done = [r for r in runs if r.get('stage') == 'final']
SIG_APL = 0.0025                                   # plik-lite's Gaussian prior on A_planck (1 ± 0.0025)
for r in runs:
    # χ² = sum of likelihood χ² (as pre-registered and as reported by Cobaya). The fit minimises −log(post), which
    # also contains A_planck's Gaussian prior; χ²_eff = χ² + ((A_planck − 1)/0.0025)² is the quantity actually
    # minimised (up to model-dependent constants from flat-prior widths) and is reported alongside.
    r['chi2_eff'] = r['chi2_total'] + ((r['minimum'].get('A_planck', 1.0) - 1.0) / SIG_APL) ** 2
    # converged = the optimiser stopped on its own criterion (last round gained < fitcore's tol = 0.1), not because
    # the rounds ran out. fit.py's own flag used 0.05, inconsistent with tol; recomputed here from the stored gains.
    rr = (r.get('optimiser') or {}).get('rounds') or []
    r['converged'] = bool(rr and rr[-1]['gain'] < 0.1)
best, seeds = {}, {}
for r in done:
    key = (r['model'], r['data'])
    seeds.setdefault(key, []).append(r)
    if key not in best or r['chi2_eff'] < best[key]['chi2_eff']:     # best = lowest minimised objective
        best[key] = r
for r in best.values():
    r['source'] = f"free fit, seed {r['seed']}"

# Profile likelihood in ΔF = ξ²/(4Λ0) (RESEARCH_LOG decision 6). Compared on χ²_eff, because −log(post) carries
# model-dependent flat-prior constants (the profile samples Λ0 only, over a ΔF-dependent range).
prof = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(HERE, 'profile', 'cmstgA_*.json')))]
for r in prof:
    r['chi2_eff'] = r['chi2_total'] + ((r['minimum'].get('A_planck', 1.0) - 1.0) / SIG_APL) ** 2
    rr = (r.get('optimiser') or {}).get('rounds') or []
    r['converged'] = bool(rr and rr[-1]['gain'] < 0.1)
    r['source'] = f"profile, ΔF = {r['dF']:.3f}"
prof_done = sorted((r for r in prof if r.get('stage') == 'final'), key=lambda r: (r['data'], r['dF']))
for r in prof_done:
    key = ('cmstgA', r['data'])
    if key in best and r['chi2_eff'] < best[key]['chi2_eff']:
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


def fmt(v, f='+.2f'):
    return '—' if v is None else format(v, f)


table, md = {}, []
hdr = (f"{'data':4} {'model':7} {'χ²':>9} {'Δχ² vs ΛCDM':>12} {'ΔAIC vs ΛCDM':>13} {'Δχ² vs w0wa':>12} "
       f"{'Δχ²_eff vs ΛCDM':>16} {'seed spread':>11}  conv  evals")
print(hdr)
md += ['| data | model | best χ² | Δχ² vs ΛCDM | ΔAIC vs ΛCDM | Δχ² vs w0wa | Δχ²_eff vs ΛCDM | seed spread (χ²_eff) | converged | fits |',
       '|---|---|---|---|---|---|---|---|---|---|']
for d in ('D1', 'D2', 'D3'):
    ref, w = best.get(('lcdm', d)), best.get(('w0wa', d))
    for m in ('lcdm', 'w0wa', 'cmstgA'):
        b = best.get((m, d))
        if b is None:
            continue
        dchi = None if ref is None else ref['chi2_total'] - b['chi2_total']          # > 0: better than ΛCDM
        daic = None if dchi is None else -dchi + 2 * K[m]                              # < 0: preferred
        dw = None if w is None else w['chi2_total'] - b['chi2_total']                  # > 0: better than w0wa
        c2 = [r['chi2_total'] for r in seeds[(m, d)]]
        ce = [r['chi2_eff'] for r in seeds[(m, d)]]
        spread = max(ce) - min(ce)                                                     # in the minimised quantity
        nconv = sum(bool(r.get('converged')) for r in seeds[(m, d)])
        deff = None if ref is None else ref['chi2_eff'] - b['chi2_eff']
        table[(m, d)] = dict(chi2=b['chi2_total'], dchi2_vs_lcdm=dchi, dAIC_vs_lcdm=daic, dchi2_vs_w0wa=dw,
                             chi2_eff=b['chi2_eff'], dchi2_eff_vs_lcdm=deff,
                             seed_chi2=sorted(c2), seed_chi2_eff=sorted(ce), seed_spread=spread, n_converged=nconv, n_fits=len(c2),
                             best_source=b['source'], minimum=b['minimum'], chi2_parts=b['chi2'],
                             sigma=b['optimiser'].get('sigma'), evaluations=b['optimiser']['evaluations'],
                             failed_evaluations=b['optimiser']['failed_evaluations'])
        print(f"{d:4} {m:7} {b['chi2_total']:9.2f} {fmt(dchi):>12} {fmt(daic):>13} {fmt(dw):>12} {fmt(deff):>16} "
              f"{spread:11.2f}  {nconv}/{len(c2)}  {b['optimiser']['evaluations']}")
        md.append(f"| {d} | {m} | {b['chi2_total']:.2f} | {fmt(dchi)} | {fmt(daic)} | {fmt(dw)} | {fmt(deff)} | {spread:.2f} | "
                  f"{nconv}/{len(c2)} | {len(c2)} |")

running = sorted(f"{r['model']}_{r['data']}_s{r['seed']} ({r.get('stage')}: χ² {r['chi2_total']:.2f})"
                 for r in runs if r.get('stage') != 'final')
if running:
    print('running:', ', '.join(running))

profile_rows = []
if prof:
    ref_eff = {d: best[('lcdm', d)]['chi2_eff'] for d in ('D1', 'D2', 'D3') if ('lcdm', d) in best}
    ref_chi = {d: best[('lcdm', d)]['chi2_total'] for d in ('D1', 'D2', 'D3') if ('lcdm', d) in best}
    print('\nProfile likelihood in ΔF = ξ²/(4Λ0) (CMSTG-A; Δ vs the ΛCDM best fit of the same data):')
    print(f"{'data':4} {'ΔF':>6} {'stage':7} {'χ²':>9} {'Δχ²':>7} {'Δχ²_eff':>8} {'Λ0':>7} {'ξ':>7} {'φ*':>6} {'H0':>6} "
          f"{'ω_b':>8} {'Y_p pull':>8} {'D/H pull':>8}  C2 C3")
    md += ['', '## Profile likelihood in ΔF = ξ²/(4Λ0) (CMSTG-A)', '',
           '| data | ΔF | stage | χ² | Δχ² vs ΛCDM | Δχ²_eff vs ΛCDM | Λ0 | ξ | φ* | H0 | ω_b | Y_p pull | D/H pull | C2 | C3 |',
           '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|']
    for r in sorted(prof, key=lambda r: (r['data'], r['dF'])):
        mn, ck = r['minimum'], cmstg_checks(r['minimum'])
        row = dict(data=r['data'], dF=r['dF'], stage=r['stage'], converged=r['converged'], chi2=r['chi2_total'],
                   chi2_eff=r['chi2_eff'], dchi2_vs_lcdm=ref_chi[r['data']] - r['chi2_total'],
                   dchi2_eff_vs_lcdm=ref_eff[r['data']] - r['chi2_eff'], minimum=mn, checks=ck)
        profile_rows.append(row)
        phi = mn['cmstg_xi'] / (2 * mn['cmstg_L0'])
        print(f"{r['data']:4} {r['dF']:6.3f} {r['stage']:7} {r['chi2_total']:9.2f} {row['dchi2_vs_lcdm']:+7.2f} "
              f"{row['dchi2_eff_vs_lcdm']:+8.2f} {mn['cmstg_L0']:7.4f} {mn['cmstg_xi']:7.4f} {phi:6.3f} {mn['H0']:6.2f} "
              f"{mn['omega_b']:8.5f} {ck['bbn']['pull_Yp']:+8.2f} {ck['bbn']['pull_DH']:+8.2f}  {str(ck['C2'])[0]}  {str(ck['C3'])[0]}")
        md.append(f"| {r['data']} | {r['dF']:.3f} | {r['stage']} | {r['chi2_total']:.2f} | {row['dchi2_vs_lcdm']:+.2f} | "
                  f"{row['dchi2_eff_vs_lcdm']:+.2f} | {mn['cmstg_L0']:.4f} | {mn['cmstg_xi']:.4f} | {phi:.3f} | {mn['H0']:.2f} | "
                  f"{mn['omega_b']:.5f} | {ck['bbn']['pull_Yp']:+.2f} | {ck['bbn']['pull_DH']:+.2f} | {ck['C2']} | {ck['C3']} |")

# Exploratory (not a pre-registered criterion): the best D1 CMSTG point that passes C2 and C3, among the finished
# profile points and the free fits.
bbn_best = None
if ('lcdm', 'D1') in best:
    cands = [dict(source=f"profile, ΔF = {x['dF']:.3f}", dF=x['dF'], chi2=x['chi2'], chi2_eff=x['chi2_eff'],
                  minimum=x['minimum'], checks=x['checks'])
             for x in profile_rows if x['data'] == 'D1' and x['stage'] == 'final']
    for r in seeds.get(('cmstgA', 'D1'), []):
        mn = r['minimum']
        cands.append(dict(source=f"free fit, seed {r['seed']}", dF=mn['cmstg_xi'] ** 2 / (4 * mn['cmstg_L0']),
                          chi2=r['chi2_total'], chi2_eff=r['chi2_eff'], minimum=mn, checks=cmstg_checks(mn)))
    ok = [x for x in cands if x['checks']['C2'] and x['checks']['C3']]
    if ok:
        bbn_best = min(ok, key=lambda x: x['chi2_eff'])
        bbn_best['dchi2_vs_lcdm'] = best[('lcdm', 'D1')]['chi2_total'] - bbn_best['chi2']
        bbn_best['dchi2_eff_vs_lcdm'] = best[('lcdm', 'D1')]['chi2_eff'] - bbn_best['chi2_eff']
        line = (f"Exploratory (not a pre-registered criterion): the best D1 CMSTG point that passes C2 and C3 "
                f"({bbn_best['source']}, ΔF = {bbn_best['dF']:.4f}, H0 = {bbn_best['minimum']['H0']:.2f}, Y_p pull "
                f"{bbn_best['checks']['bbn']['pull_Yp']:+.2f}) has Δχ² = {bbn_best['dchi2_vs_lcdm']:+.2f} "
                f"(Δχ²_eff {bbn_best['dchi2_eff_vs_lcdm']:+.2f}) vs ΛCDM.")
        print(line)
        md += ['', f'- {line}']

verdict, checks, crit = 'INCOMPLETE', {}, {}
if ('cmstgA', 'D1') in table and ('lcdm', 'D1') in table:
    t1 = table[('cmstgA', 'D1')]
    checks = cmstg_checks(t1['minimum'])
    crit['a'] = bool(t1['dchi2_vs_lcdm'] >= 6)
    crit['b'] = bool(checks['C2'] and checks['C3'])
    crit['c'] = bool(any((table.get(('cmstgA', d)) or {}).get('dchi2_vs_lcdm') is not None
                         and table[('cmstgA', d)]['dchi2_vs_lcdm'] >= 6 for d in ('D2', 'D3')))
    complete = all(k in table for k in [(m, d) for m in K for d in ('D1', 'D2', 'D3')])
    verdict = ('PASS' if all(crit.values()) else 'FAIL') if complete else 'INCOMPLETE'
    mn = t1['minimum']
    lines = [f"CMSTG-A best fit (D1, {t1['best_source']}): Λ0 = {mn['cmstg_L0']:.4f}, ξ = {mn['cmstg_xi']:.4f}, ξ²/(4Λ0) = "
             f"{mn['cmstg_xi']**2 / (4 * mn['cmstg_L0']):.5f}, H0 = {mn['H0']:.2f}; F0 = {checks['F0']:.5f} "
             f"(attractor {checks['F_attractor']:.5f})",
             f"(a) Δχ²(D1) = {t1['dchi2_vs_lcdm']:+.2f} ≥ 6: {crit['a']}",
             f"(b) C2 {checks['C2']} (|γ−1| = {abs(checks['local_gravity']['gamma_minus_1']):.1e}), C3 {checks['C3']} "
             f"(Y_p pull {checks['bbn']['pull_Yp']:+.2f}, D/H pull {checks['bbn']['pull_DH']:+.2f})",
             f"(c) Δχ² ≥ 6 on D2 or D3: {crit['c']}",
             f"VERDICT (pre-registered criteria): {verdict}"]
    print('\n' + '\n'.join(lines))
    md += [''] + [f'- {x}' for x in lines]

open(os.path.join(HERE, 'results.md'), 'w').write('# SIM155 best fits (generated by summarize.py)\n\n' + '\n'.join(md) + '\n')
json.dump(dict(verdict=verdict, criteria=crit, table={f'{m}_{d}': v for (m, d), v in table.items()},
               cmstg_checks=checks, profile=profile_rows, bbn_allowed_best=bbn_best, n_runs=len(runs),
               n_finished=len(done), running=running),
          open(os.path.join(HERE, 'output.json'), 'w'), indent=2, default=float)
