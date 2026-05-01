"""
Regenerate fig6_polarization.pdf for the master paper.
SIM95: CMSTG vs LCDM CMB EE and TE power spectra.
Top panel: fractional deviations for CMSTG joint and BAO-only vs LCDM.
Bottom panel: bar chart of chi2 values (EE and TE) for all three cosmologies.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json, os

matplotlib.rcParams.update({
    'font.size': 11, 'axes.labelsize': 12, 'axes.titlesize': 12,
    'legend.fontsize': 10, 'xtick.labelsize': 10, 'ytick.labelsize': 10,
    'text.usetex': False, 'figure.dpi': 150,
    'axes.linewidth': 1.2, 'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.top': True, 'ytick.right': True,
})

SIMDIR = os.path.expanduser('~/Ordered_Simulations/SIM95/Outputs')
HERE   = os.path.dirname(os.path.abspath(__file__))

# ── Load CLASS C_ell files ──────────────────────────────────────────────────
def load_cl(subdir):
    path = os.path.join(SIMDIR, subdir, '00_00_cl_lensed.dat')
    d = np.loadtxt(path)
    ell = d[:, 0]
    EE  = d[:, 2]   # col 2 = EE
    TE  = d[:, 3]   # col 3 = TE
    return ell, EE, TE

ell_lcdm, EE_lcdm, TE_lcdm = load_cl('class_lcdm_planck')
ell_joint, EE_joint, TE_joint = load_cl('class_cmstg_joint')
ell_bao,   EE_bao,   TE_bao   = load_cl('class_cmstg_bao')

# Trim to common ell range 2-1996 (matches chi2 computation in SIM95)
mask = (ell_lcdm >= 2) & (ell_lcdm <= 1996)
ell  = ell_lcdm[mask]
EE_l = EE_lcdm[mask];  TE_l = TE_lcdm[mask]
EE_j = EE_joint[mask]; TE_j = TE_joint[mask]
EE_b = EE_bao[mask];   TE_b = TE_bao[mask]

# Fractional deviations vs LCDM (avoid division by zero for TE)
eps = 1e-30
dEE_joint = (EE_j - EE_l) / (np.abs(EE_l) + eps) * 100   # percent
dTE_joint = (TE_j - TE_l) / (np.abs(TE_l) + eps) * 100
dEE_bao   = (EE_b - EE_l) / (np.abs(EE_l) + eps) * 100
dTE_bao   = (TE_b - TE_l) / (np.abs(TE_l) + eps) * 100

# Running mean smoother
def smooth(x, w=30):
    kernel = np.ones(w) / w
    return np.convolve(x, kernel, mode='same')

# ── Load chi2 summary ────────────────────────────────────────────────────────
with open(os.path.join(os.path.expanduser('~/Ordered_Simulations/SIM95/Outputs'),
                       'sim95_diagnostics.json')) as f:
    diag = json.load(f)

pl = diag['planck_likelihood']
chi2_ee = [pl['EE']['chi2_lcdm_planck'],
           pl['EE']['chi2_cmstg_joint'],
           pl['EE']['chi2_cmstg_bao_only']]
chi2_te = [pl['TE']['chi2_lcdm_planck'],
           pl['TE']['chi2_cmstg_joint'],
           pl['TE']['chi2_cmstg_bao_only']]
labels  = [r'$\Lambda$CDM', 'CMSTG\njoint', 'CMSTG\nBAO-only']
colors  = ['C2', 'C0', 'C3']

# ── Plot ─────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(12, 9),
                         gridspec_kw={'height_ratios': [2, 1]})
fig.suptitle(r'SIM95 — CMSTG vs $\Lambda$CDM CMB EE and TE polarization', fontsize=13)

ax_ee, ax_te = axes[0, 0], axes[0, 1]
ax_cb_ee, ax_cb_te = axes[1, 0], axes[1, 1]

# Top-left: EE fractional deviation
ax_ee.plot(ell, smooth(dEE_joint), color='C0', lw=1.5,
           label=r'CMSTG joint best-fit $(\Delta\chi^2_{\rm EE}=+31.5)$')
ax_ee.plot(ell, smooth(dEE_bao),   color='C3', lw=1.5, ls='--',
           label=r'CMSTG BAO-only $(\chi^2_{\rm EE}=4353)$')
ax_ee.axhline(0, color='C2', lw=1.2, ls='-', label=r'$\Lambda$CDM (reference)')
ax_ee.set_xlim(2, 1996)
ax_ee.set_ylim(-5, 5)
ax_ee.set_xlabel(r'Multipole $\ell$')
ax_ee.set_ylabel(r'$\Delta C_\ell^{EE}/C_\ell^{EE}$ [%]')
ax_ee.set_title(r'EE polarization (RMS$_{\rm joint}=0.19\%$)')
ax_ee.legend(fontsize=9, loc='upper right')

# Top-right: TE fractional deviation (mask near TE=0 crossings)
te_denom_ok = np.abs(TE_l) > 0.01 * np.max(np.abs(TE_l))
ell_te  = ell[te_denom_ok]
dTE_j_m = dTE_joint[te_denom_ok]
dTE_b_m = dTE_bao[te_denom_ok]

ax_te.plot(ell_te, smooth(dTE_j_m, 30), color='C0', lw=1.5,
           label=r'CMSTG joint $(\Delta\chi^2_{\rm TE}=+22.8)$')
ax_te.plot(ell_te, smooth(dTE_b_m, 30), color='C3', lw=1.5, ls='--',
           label=r'CMSTG BAO-only $(\chi^2_{\rm TE}=3110)$')
ax_te.axhline(0, color='C2', lw=1.2, ls='-', label=r'$\Lambda$CDM (reference)')
ax_te.set_xlim(2, 1996)
ax_te.set_ylim(-10, 10)
ax_te.set_xlabel(r'Multipole $\ell$')
ax_te.set_ylabel(r'$\Delta C_\ell^{TE}/C_\ell^{TE}$ [%]')
ax_te.set_title(r'TE cross-spectrum (RMS$_{\rm joint}=0.43\%$)')
ax_te.legend(fontsize=9, loc='upper right')

# Bottom-left: EE chi2 bar chart
x = np.arange(3)
bars_ee = ax_cb_ee.bar(x, chi2_ee, color=colors, edgecolor='k', linewidth=0.8,
                       width=0.55)
ax_cb_ee.set_yscale('log')
ax_cb_ee.set_xticks(x)
ax_cb_ee.set_xticklabels(labels, fontsize=10)
ax_cb_ee.set_ylabel(r'$\chi^2_{\rm EE}$ (approx.\ Gaussian, 1995 modes)')
ax_cb_ee.set_title('EE polarization likelihood')
for bar, val in zip(bars_ee, chi2_ee):
    ax_cb_ee.text(bar.get_x() + bar.get_width()/2, val * 1.5,
                  f'{val:.1f}', ha='center', va='bottom', fontsize=9.5)
# Mark dof reference
ax_cb_ee.axhline(pl['EE']['n_modes'], color='gray', lw=1, ls=':', alpha=0.7)
ax_cb_ee.text(2.4, pl['EE']['n_modes'] * 1.1, 'dof=1995',
              color='gray', fontsize=8.5, ha='right')

# Bottom-right: TE chi2 bar chart
bars_te = ax_cb_te.bar(x, chi2_te, color=colors, edgecolor='k', linewidth=0.8,
                       width=0.55)
ax_cb_te.set_yscale('log')
ax_cb_te.set_xticks(x)
ax_cb_te.set_xticklabels(labels, fontsize=10)
ax_cb_te.set_ylabel(r'$\chi^2_{\rm TE}$ (approx.\ Gaussian, 1995 modes)')
ax_cb_te.set_title('TE cross-spectrum likelihood')
for bar, val in zip(bars_te, chi2_te):
    ax_cb_te.text(bar.get_x() + bar.get_width()/2, val * 1.5,
                  f'{val:.1f}', ha='center', va='bottom', fontsize=9.5)
ax_cb_te.axhline(pl['TE']['n_modes'], color='gray', lw=1, ls=':', alpha=0.7)
ax_cb_te.text(2.4, pl['TE']['n_modes'] * 1.1, 'dof=1995',
              color='gray', fontsize=8.5, ha='right')

plt.tight_layout(rect=[0, 0, 1, 0.96])
out = os.path.join(HERE, 'fig6_polarization.pdf')
plt.savefig(out, bbox_inches='tight')
plt.close()
print(f"Saved {out}")
