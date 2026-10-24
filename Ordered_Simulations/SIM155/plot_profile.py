#!/usr/bin/env python3
"""SIM155 figure: CMSTG-A profile likelihood in ΔF = ξ²/(4Λ0) on D1, with the BBN pulls. Reads output.json
(written by summarize.py) and fits_v3/; writes figures/sim155_profile.png."""
import glob, json, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
out = json.load(open(os.path.join(HERE, 'output.json')))
rows = sorted((r for r in out['profile'] if r['data'] == 'D1'), key=lambda r: r['dF'])
ref = out['table']['lcdm_D1']['chi2_eff']
w0wa = out['table']['w0wa_D1']['dchi2_eff_vs_lcdm']
free = []
for f in sorted(glob.glob(os.path.join(HERE, 'fits_v3', 'cmstgA_D1_s*.json'))):
    r = json.load(open(f))
    m = r['minimum']
    ce = r['chi2_total'] + ((m.get('A_planck', 1.0) - 1.0) / 0.0025) ** 2
    free.append((m['cmstg_xi'] ** 2 / (4 * m['cmstg_L0']), ref - ce, r['seed']))

dF = np.array([r['dF'] for r in rows])
d_eff = np.array([r['dchi2_eff_vs_lcdm'] for r in rows])
fin = np.array([r['stage'] == 'final' for r in rows])
yp = np.array([r['checks']['bbn']['pull_Yp'] for r in rows])
dh = np.array([r['checks']['bbn']['pull_DH'] for r in rows])
h0 = np.array([r['minimum']['H0'] for r in rows])

fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.2, 7.6), sharex=True, gridspec_kw=dict(height_ratios=[3, 2]))
# BBN-excluded region: where the Y_p or D/H pull exceeds 2 (interpolated between profile points)
worst = np.maximum(np.abs(yp), np.abs(dh))
edge = None
for i in range(1, len(dF)):
    if worst[i - 1] <= 2 < worst[i]:
        edge = dF[i - 1] + (2 - worst[i - 1]) * (dF[i] - dF[i - 1]) / (worst[i] - worst[i - 1])
        break
for ax in (a1, a2):
    if edge is not None:
        ax.axvspan(edge, max(dF.max(), 0.06) + 0.005, color='0.88', zorder=0)
a1.axhline(6, ls='--', color='C3', lw=1.2, label='pre-registered pass line, Δχ² = 6')
a1.axhline(w0wa, ls=':', color='C2', lw=1.5, label=f'w0waCDM best fit (Δχ²_eff = {w0wa:.1f})')
a1.axhline(0, color='k', lw=0.8)
a1.plot(dF[fin], d_eff[fin], 'o-', color='C0', label='CMSTG-A profile (fixed ΔF, all else minimised)')
if (~fin).any():
    a1.plot(dF[~fin], d_eff[~fin], 'o', mfc='none', color='C0', label='profile point, not finished')
for x, y, s in free:
    a1.plot(x, y, 'x', color='C1', ms=8, mew=2, label='free fits (3 seeds)' if s == 1 else None)
a1.set_ylabel('Δχ²_eff vs ΛCDM  (D1: Planck 2018 + DESI DR2 + Pantheon+)')
a1.set_title('SIM155: CMSTG attractor (family A) — profile in the attractor strength')
if edge is not None:
    a1.text(edge + 0.001, a1.get_ylim()[0] + 0.3, 'BBN excluded\n(Y_p or D/H > 2σ)', fontsize=9, color='0.35')
a1.legend(fontsize=8, loc='upper left')
a2.axhline(2, ls='--', color='0.4', lw=1)
a2.plot(dF, yp, 's-', color='C4', label='Y_p pull (σ)')
a2.plot(dF, dh, 'd-', color='C5', label='D/H pull (σ)')
a2.set_ylabel('BBN pull at the profile point')
a2.set_xlabel('ΔF = F* − ½ = ξ²/(4Λ0)   (fractional weakening of G since BBN ≈ ΔF/F*)')
t = a2.twinx()
t.plot(dF, h0, '^:', color='C7', ms=4, label='H0')
t.set_ylabel('H0 at the profile point [km/s/Mpc]', color='C7')
a2.legend(fontsize=8, loc='upper left')
fig.tight_layout()
os.makedirs(os.path.join(HERE, 'figures'), exist_ok=True)
fig.savefig(os.path.join(HERE, 'figures', 'sim155_profile.png'), dpi=150)
print('wrote figures/sim155_profile.png; BBN edge at ΔF ≈', None if edge is None else round(float(edge), 4))
