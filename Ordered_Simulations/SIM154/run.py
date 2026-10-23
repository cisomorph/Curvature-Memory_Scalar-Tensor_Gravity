#!/usr/bin/env python3
"""
SIM154 — Attractor scalar-tensor models: wide scan + local gravity + BBN (ROADMAP Phase 2).
Implements SPEC.md (pre-registered 2026-10-22). Families:
  A (subtractive, SIM132): F = ½(1 − 2Λ0Ψ²) + ξΨ
  B (exponential, SIM134): F = ½e^{−2Λ0Ψ²} + ξΨ
Ψ = Ψ' = 0 at z = 1e5; H0 solved from θ*; DESI DR1 + DR2; criteria C1–C4 from SPEC.md.
"""
import json, os, sys, warnings
from datetime import datetime
from multiprocessing import Pool
import numpy as np
warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'common'))
import cmstg_obs as obs
from st_solver import ScalarTensorBG
import audit_driver as drv

LAMBDA0 = [0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5]
XI = [1 / 24, 1 / 12, 1 / 8, 1 / 6, 1 / 4]

# BBN (SPEC): Kneller & Steigman 2004 fits; observations with fit errors added in quadrature
ETA10 = 273.9 * obs.OMEGA_B
YP_OBS, YP_ERR, YP_FIT = 0.245, 0.003, 0.0006
DH_OBS, DH_ERR, DH_FIT_FRAC = 2.547, 0.025, 0.03


def family(name, L0, xi):
    if name == 'A':
        F = lambda p: 0.5 * (1 - 2 * L0 * p * p) + xi * p
        Fp = lambda p: -2 * L0 * p + xi
        attractor = xi / (2 * L0)
    else:
        F = lambda p: 0.5 * np.exp(-2 * L0 * p * p) + xi * p
        Fp = lambda p: -2 * L0 * p * np.exp(-2 * L0 * p * p) + xi
        attractor = None
        if xi < np.sqrt(L0) * np.exp(-0.5):
            from scipy.optimize import brentq
            attractor = brentq(Fp, 1e-9, 0.5 / np.sqrt(L0))
    return F, Fp, attractor


def bbn(S):
    yp = 0.2485 + 0.0016 * ((ETA10 - 6) + 100 * (S - 1))
    dh = 2.64 * (6 / (ETA10 - 6 * (S - 1)))**1.6
    pull_yp = (yp - YP_OBS) / np.hypot(YP_ERR, YP_FIT)
    pull_dh = (dh - DH_OBS) / np.hypot(DH_ERR, DH_FIT_FRAC * dh)
    return dict(S=S, Yp=yp, DH=dh, pull_Yp=pull_yp, pull_DH=pull_dh,
                dNeff_equiv=(S * S - 1) * 43 / 7, ok=bool(abs(pull_yp) <= 2 and abs(pull_dh) <= 2))


def run_point(args):
    fam, L0, xi = args
    F, Fp, att = family(fam, L0, xi)
    m = ScalarTensorBG(F=lambda p, N: F(p), F_psi=lambda p, N: Fp(p))
    rec = dict(family=fam, Lambda0=L0, xi=xi, attractor=att)
    bg = m.theta_matched(0.0, 0.0, h_guess=0.674)
    if bg is None:
        probe = m.calibrated(0.0, 0.0, 0.674)
        rec['status'] = 'no theta* solution' if probe is not None else 'integration breakdown (F -> 0 or constraint)'
        rec.update(C1=False, C2=False, C3=False, C4=probe is not None, survives=False)
        return rec
    o = m.observe(bg, n_params=2)
    p0, y0, F0 = bg['psi0'], bg['y0'], bg['F0']
    lg = drv.local_gravity(F0, [Fp(p0)], Fdot_over_HF=Fp(p0) * y0 / F0, H0=bg['H0'])
    b = bbn(np.sqrt(F0 / F(0.0)))
    c1 = o['DR2']['gof_sigma'] < 2.0 and o['DR2']['delta_chi2_tot_vs_LCDM'] >= 4.0
    c2 = lg['cassini_ok'] and lg['llr_ok']
    c4 = bool(np.min(bg['F']) > 0)
    rec.update(status='ok', H0=bg['H0'], psi0=p0, y0=y0, F0=F0, Fpsi0=float(Fp(p0)),
               F_zstar=float(np.interp(obs.Z_STAR, bg['z'], bg['F'])),
               theta_star_100=o['theta_star_100'], r_drag=o['r_drag'],
               DR1=o['DR1'], DR2=o['DR2'], local_gravity=lg, bbn=b,
               C1=bool(c1), C2=bool(c2), C3=b['ok'], C4=c4)
    rec['survives'] = bool(c1 and c2 and b['ok'] and c4)
    return rec


if __name__ == '__main__':
    t0 = datetime.utcnow()
    jobs = [(f, L0, xi) for f in ('A', 'B') for L0 in LAMBDA0 for xi in XI]
    with Pool(10) as pool:
        rows = pool.map(run_point, jobs, chunksize=1)

    print("=" * 100)
    print("SIM154 — attractor models: θ*-matched DESI + local gravity + BBN   (SPEC.md, pre-registered)")
    print("=" * 100)
    hdr = (f"{'fam':>3} {'Λ0':>5} {'ξ':>6} {'H0':>6} {'F0':>7} {'Fψ(Ψ0)':>9} {'Δχ²DR1':>7} {'Δχ²DR2':>7} "
           f"{'GOF2':>5} {'|γ−1|':>8} {'Ġ/G':>9} {'Yp':>6} {'pYp':>5} {'pDH':>5} {'ΔNeff':>6}  C1 C2 C3 C4")
    print(hdr); print('-' * len(hdr))
    for r in rows:
        if r['status'] != 'ok':
            print(f"{r['family']:>3} {r['Lambda0']:>5} {r['xi']:6.4f}  {r['status']}")
            continue
        lg, b = r['local_gravity'], r['bbn']
        flag = lambda x: ' ✓' if x else ' ✗'
        print(f"{r['family']:>3} {r['Lambda0']:>5} {r['xi']:6.4f} {r['H0']:6.2f} {r['F0']:7.4f} "
              f"{r['Fpsi0']:9.2e} {r['DR1']['delta_chi2_tot_vs_LCDM']:+7.2f} {r['DR2']['delta_chi2_tot_vs_LCDM']:+7.2f} "
              f"{r['DR2']['gof_sigma']:5.2f} {lg['gamma_minus_1']:8.1e} {lg['Gdot_over_G_per_yr']:9.1e} "
              f"{b['Yp']:6.4f} {b['pull_Yp']:+5.1f} {b['pull_DH']:+5.1f} {b['dNeff_equiv']:6.2f} "
              f"{flag(r['C1'])}{flag(r['C2'])}{flag(r['C3'])}{flag(r['C4'])}")
    surv = [r for r in rows if r.get('survives')]
    print("=" * 100)
    for f in ('A', 'B'):
        sub = [r for r in rows if r['family'] == f and r['status'] == 'ok']
        n = len(sub)
        print(f"Family {f}: {n} solved; C1 {sum(r['C1'] for r in sub)}, C2 {sum(r['C2'] for r in sub)}, "
              f"C3 {sum(r['C3'] for r in sub)}, C1&C3 {sum(r['C1'] and r['C3'] for r in sub)}, "
              f"C1&C2&C3&C4 {sum(r['survives'] for r in sub)}")
    print(f"SURVIVORS (C1–C4): {len(surv)}")
    for r in surv:
        print(f"  family {r['family']} Λ0={r['Lambda0']} ξ={r['xi']:.4f}: Δχ² DR1 {r['DR1']['delta_chi2_tot_vs_LCDM']:+.2f}, "
              f"DR2 {r['DR2']['delta_chi2_tot_vs_LCDM']:+.2f}, H0 {r['H0']:.2f}, Yp {r['bbn']['Yp']:.4f}")
    verdict = 'SURVIVORS FOUND' if surv else 'NO SURVIVORS'
    print(f"VERDICT: {verdict}")

    # figures: Δχ²(DR2) maps with C2/C3 overlays
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for ax, f in zip(axes, ('A', 'B')):
        M = np.full((len(XI), len(LAMBDA0)), np.nan)
        for r in rows:
            if r['family'] == f and r['status'] == 'ok':
                M[XI.index(r['xi']), LAMBDA0.index(r['Lambda0'])] = r['DR2']['delta_chi2_tot_vs_LCDM']
        im = ax.imshow(M, origin='lower', cmap='RdBu', vmin=-25, vmax=25, aspect='auto')
        for r in rows:
            if r['family'] == f and r['status'] == 'ok':
                i, j = XI.index(r['xi']), LAMBDA0.index(r['Lambda0'])
                ax.text(j, i, ('C' if r['C2'] else '·') + ('B' if r['C3'] else '·'), ha='center', va='center', fontsize=9)
        ax.set_xticks(range(len(LAMBDA0))); ax.set_xticklabels(LAMBDA0)
        ax.set_yticks(range(len(XI))); ax.set_yticklabels([f'{x:.3f}' for x in XI])
        ax.set_xlabel('Λ0'); ax.set_ylabel('ξ')
        ax.set_title(f'Family {f}: Δχ²(DESI DR2 + θ*) vs ΛCDM\n"C" = passes Solar System, "B" = passes BBN')
        fig.colorbar(im, ax=ax)
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(HERE, 'figures', f'sim154_scan.{ext}'), dpi=150)
    plt.close(fig)

    json.dump(dict(sim='SIM154', started=t0.isoformat() + 'Z', finished=datetime.utcnow().isoformat() + 'Z',
                   spec='SPEC.md (pre-registered)', verdict=verdict, n_points=len(rows), survivors=surv, rows=rows,
                   lcdm_reference={rel: obs.lcdm_reference(rel) for rel in ('DR1', 'DR2')}),
              open(os.path.join(HERE, 'output.json'), 'w'), indent=2, default=float)
    print("SIM154 complete.")
