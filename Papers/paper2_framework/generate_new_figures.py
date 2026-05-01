"""
Generate new figures for Paper 2 (CMSTG Framework):
  fig8_w_of_z.pdf  -- w0-wa plane (SIM109)
  fig9_sound_horizon.pdf -- r_d vs Lambda0 with N_eff correction (SIM101)
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse
import os

OUT = os.path.join(os.path.dirname(__file__), 'figures')
os.makedirs(OUT, exist_ok=True)

# ── Matplotlib style ──────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'DejaVu Serif',
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.dpi': 150,
    'lines.linewidth': 1.5,
    'axes.grid': True,
    'grid.alpha': 0.3,
    'text.usetex': False,
})

# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 8: w0-wa plane (SIM109)
# ─────────────────────────────────────────────────────────────────────────────

def draw_ellipse(ax, w0_c, wa_c, sw0, swa, rho, n_sigma=1, **kwargs):
    """Draw a confidence ellipse for (w0, wa) with correlation rho."""
    # Eigenvalue decomposition of the covariance matrix
    cov = np.array([[sw0**2, rho*sw0*swa],
                    [rho*sw0*swa, swa**2]])
    vals, vecs = np.linalg.eigh(cov)
    # Width and height in units of sigma
    scale = np.sqrt(vals * (2 * n_sigma))  # chi2=2 for 1-sigma in 2D
    angle = np.degrees(np.arctan2(vecs[1, -1], vecs[0, -1]))
    e = Ellipse(xy=(w0_c, wa_c), width=2*scale[0], height=2*scale[1],
                angle=angle, **kwargs)
    ax.add_patch(e)

fig, ax = plt.subplots(figsize=(6.5, 5.5))

# CMSTG (SIM109): w0 = -0.992, wa = -0.082
w0_cmstg = -0.992
wa_cmstg = -0.082

# Planck 2018 wCDM: w0 = -1.03 +/- 0.03, wa = 0 (fixed)
# Shown as a band: w0 in [-1.06, -1.00], wa = 0
ax.axhspan(-0.3, 0.3, color='grey', alpha=0.15, zorder=0,
           label='$\Lambda$CDM ($w_0=-1$, $w_a=0$)')
ax.axvspan(-1.06, -1.00, color='grey', alpha=0.25, zorder=1)

# DESI Y1 1-sigma and 2-sigma ellipses
# Central values: w0 = -0.827, wa = -0.75
# Uncertainties: sigma_w0 = 0.060, sigma_wa = 0.29, rho ~ -0.7 (typical)
w0_desi = -0.827
wa_desi = -0.75
sw0_desi = 0.060
swa_desi = 0.29
rho_desi = -0.70

draw_ellipse(ax, w0_desi, wa_desi, sw0_desi, swa_desi, rho_desi,
             n_sigma=2, facecolor='#FFD580', edgecolor='orange',
             linewidth=1.2, zorder=2, label='DESI Y1 $2\sigma$')
draw_ellipse(ax, w0_desi, wa_desi, sw0_desi, swa_desi, rho_desi,
             n_sigma=1, facecolor='#FFA500', edgecolor='darkorange',
             linewidth=1.5, zorder=3, label='DESI Y1 $1\sigma$')
ax.plot(w0_desi, wa_desi, 'o', color='darkorange', ms=5, zorder=5)

# Planck 2018 wCDM point
ax.errorbar(-1.03, 0, xerr=0.03, fmt='s', color='dimgrey', ms=7,
            capsize=4, zorder=6, label='Planck 2018 wCDM')

# CMSTG point
ax.plot(w0_cmstg, wa_cmstg, '*', color='royalblue', ms=16, zorder=7,
        label='CMSTG Phase 1 locked action', markeredgecolor='navy', markeredgewidth=0.5)

# Lambda CDM point
ax.plot(-1.0, 0.0, 'D', color='black', ms=8, zorder=7,
        label='$\Lambda$CDM', markeredgecolor='black')

# Annotation arrows showing sigma distances
ax.annotate('', xy=(-1.03, 0), xytext=(w0_cmstg, wa_cmstg),
            arrowprops=dict(arrowstyle='->', color='dimgrey',
                           lw=1.2, connectionstyle='arc3,rad=0.15'))
ax.text(-1.012, -0.041, '$1.3\sigma$', fontsize=9, color='dimgrey', ha='center')

ax.set_xlabel('$w_0$')
ax.set_ylabel('$w_a$')
ax.set_title('Dark Energy Equation of State (SIM109)', pad=8)
ax.set_xlim(-1.22, -0.65)
ax.set_ylim(-1.45, 0.45)
ax.axhline(0, color='black', lw=0.5, ls='--', alpha=0.4)
ax.axvline(-1.0, color='black', lw=0.5, ls='--', alpha=0.4)
ax.legend(loc='lower left', framealpha=0.9, fontsize=9)
ax.text(0.98, 0.97, 'CMSTG vs Planck: $1.3\sigma$\nCMSTG vs DESI Y1: $3.6\sigma$',
        transform=ax.transAxes, ha='right', va='top', fontsize=9,
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

fig.tight_layout()
fig.savefig(os.path.join(OUT, 'fig8_w_of_z.pdf'), bbox_inches='tight')
fig.savefig(os.path.join(OUT, 'fig8_w_of_z.png'), bbox_inches='tight')
plt.close(fig)
print("fig8 done")

# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 9: r_d vs Lambda0 with N_eff correction (SIM101)
# ─────────────────────────────────────────────────────────────────────────────

# SIM101 data: shift is dr_d / r_d = Delta_r / r0
# r0 (LCDM with proper N_eff) = 147.1 Mpc (Planck value)
# Without N_eff: r0_approx = 153.05 Mpc
# The fractional shifts are the same; just anchor to 147.1 Mpc
r0_planck = 147.1  # Mpc, N_eff-corrected

Lambda0_vals = np.array([0.000, 0.003, 0.008, 0.050, 0.095])
# Fractional shifts from SIM101 (in ppm): 0, 0.68, 1.8, 11, 21
frac_shifts_ppm = np.array([0, 0.68, 1.8, 11.0, 21.0])
r_d_cmstg = r0_planck * (1 - frac_shifts_ppm * 1e-6)

# Dense curve for plotting
lam_dense = np.linspace(0, 0.10, 300)
# Linear interpolation of fractional shift vs Lambda0^2 (approximately linear in Lambda0)
# From the table: shift grows roughly as Lambda0 * constant
# Fit: frac_ppm ~ 220 * Lambda0 (approximately)
frac_interp = np.interp(lam_dense, Lambda0_vals, frac_shifts_ppm)
r_d_dense = r0_planck * (1 - frac_interp * 1e-6)

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

# Left panel: r_d vs Lambda0
ax = axes[0]
ax.plot(lam_dense, r_d_dense, 'b-', lw=2, label='CMSTG $r_d$')
ax.plot(Lambda0_vals, r_d_cmstg, 'bo', ms=7, zorder=5)
ax.axhline(r0_planck, color='red', ls='--', lw=1.5,
           label='Planck $r_d = 147.1$ Mpc')
ax.axhline(r0_planck * 0.997, color='red', ls=':', lw=1.0, alpha=0.6)
ax.axhline(r0_planck * 1.003, color='red', ls=':', lw=1.0, alpha=0.6,
           label='BAO precision $(\pm 0.3\%)$')
ax.fill_between(lam_dense, r0_planck * 0.997, r0_planck * 1.003,
                color='red', alpha=0.07, label='_nolegend_')
ax.axvline(0.003, color='grey', ls='--', lw=1.0, alpha=0.7,
           label='$\Lambda_0 = 0.003$ (locked)')
ax.axvline(0.095, color='orange', ls='--', lw=1.0, alpha=0.7,
           label='$\Lambda_0 = 0.095$ (95% bound)')
ax.set_xlabel('$\Lambda_0$ [$M_{\mathrm{Pl}}^{-2}$]')
ax.set_ylabel('$r_d$ [Mpc]')
ax.set_title('Sound Horizon vs Coupling Strength')
ax.set_xlim(0, 0.10)
ax.set_ylim(147.05, 147.20)
ax.legend(fontsize=9, loc='upper right')

# Right panel: fractional shift in ppm
ax = axes[1]
ax.plot(lam_dense, frac_interp, 'b-', lw=2, label='CMSTG $|\Delta r_d / r_d|$')
ax.plot(Lambda0_vals, frac_shifts_ppm, 'bo', ms=7, zorder=5)
ax.axhline(3000, color='red', ls='--', lw=1.5,
           label='BAO precision ($\sim 3000$ ppm)')
ax.fill_between(lam_dense, 0, 3000, color='green', alpha=0.06)
ax.axvline(0.003, color='grey', ls='--', lw=1.0, alpha=0.7,
           label='$\Lambda_0 = 0.003$')
ax.axvline(0.095, color='orange', ls='--', lw=1.0, alpha=0.7,
           label='95% upper bound')
ax.text(0.050, 600, 'CMSTG $r_d$ degenerate\nwith $\Lambda$CDM\n(all couplings)',
        fontsize=9, ha='center', va='bottom', color='darkgreen',
        bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8))
ax.set_xlabel('$\Lambda_0$ [$M_{\mathrm{Pl}}^{-2}$]')
ax.set_ylabel('$|\Delta r_d / r_d|$ [ppm]')
ax.set_title('Fractional Sound Horizon Shift')
ax.set_xlim(0, 0.10)
ax.set_ylim(0, 3200)
ax.legend(fontsize=9, loc='upper left')

fig.suptitle('Sound Horizon from First Principles (SIM101) — $N_{\mathrm{eff}}$-corrected baseline',
             fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(OUT, 'fig9_sound_horizon.pdf'), bbox_inches='tight')
fig.savefig(os.path.join(OUT, 'fig9_sound_horizon.png'), bbox_inches='tight')
plt.close(fig)
print("fig9 done")

print("All Paper 2 new figures generated.")
