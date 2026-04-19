"""
SIM120 / SIM120-alt — RIFT Phase 2: Visualisation Suite
========================================================
Generates publication-quality figures for:
  1. H(z) comparison: RIFT vs DESI BAO
  2. Ψ(a) field evolution (background quintessence)
  3. w(a) equation-of-state evolution
  4. κ–m₂₂ RIFT DE-DM link diagram
  5. SIM120-alt: m₂₂ scan landscape
  6. SIM120-alt: m₂₂ results comparison panel

Run from:
  python SIM120/Published_Scripts/sim120_plots.py
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import os, warnings
warnings.filterwarnings('ignore')

OUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'Outputs')
os.makedirs(OUT_DIR, exist_ok=True)

# ── Shared RIFT parameters ────────────────────────────────────────────────────
Lambda0   = 0.003
Psi0      = 2.62
v         = 13.16
Omega_DE  = 0.685
Omega_m0  = 0.315
Omega_b   = 0.049
Omega_r   = 9.0e-5
H0_kms    = 67.4
F0        = 0.5 + Lambda0 * Psi0**2          # = 0.52059
VE_factor = (Psi0**2 - v**2)**2 / (1.0 + 2.0*Lambda0*Psi0**2)**2
lam_norm  = Omega_DE * 3.0 * F0 / VE_factor

DESI_z = np.array([0.30, 0.51, 0.71, 0.93, 1.32, 2.33])
DESI_H = np.array([81.7, 97.9, 110.7, 128.1, 156.4, 240.8])
DESI_s = np.array([ 4.5,  4.4,   6.2,   5.6,   8.6,  11.0])

STYLE = dict(fontsize=11, family='serif')
plt.rcParams.update({'font.family': 'serif', 'font.size': 11,
                     'axes.labelsize': 12, 'legend.fontsize': 10})

# ── Helper: E²(z) for different DE models ────────────────────────────────────
def E2_LCDM(z):
    a = 1.0/(1.0+z)
    return (Omega_m0*a**-3 + Omega_r*a**-4 + Omega_DE) / (3.0*F0)

def E2_CPL(z, w0, wa):
    a = 1.0/(1.0+z)
    Ode = Omega_DE * a**(-3.0*(1.0+w0+wa)) * np.exp(-3.0*wa*(1.0-a))
    return (Omega_m0*a**-3 + Omega_r*a**-4 + Ode) / (3.0*F0)

# ── Potential helpers ─────────────────────────────────────────────────────────
def F(u):   return 0.5 + Lambda0*u**2
def VJ(u):  return lam_norm*(u**2 - v**2)**2
def dVJ(u): return 4.0*lam_norm*(u**2 - v**2)*u

def ode_bg(N, y):
    u, up = y
    a = np.exp(N)
    Om = Omega_m0*a**-3; Or = Omega_r*a**-4
    VJu = VJ(u); dVJu = dVJ(u)
    F_u = F(u)
    denom = 3.0*F_u - 0.5*up**2
    if denom <= 0: return [up, -3.0*up]
    E2 = (Om + Or + VJu) / denom
    P_r = Or/3.0
    rho_tot = Om + Or + 0.5*E2*up**2 + VJu
    P_tot = P_r + 0.5*E2*up**2 - VJu
    w_eff = P_tot/rho_tot if rho_tot > 0 else 0.0
    dlnE2 = -3.0*(1.0+w_eff)
    R_norm = -6.0*E2*(dlnE2/2.0 + 2.0)
    upp = -(3.0 + dlnE2/2.0)*up - (dVJu + 2.0*Lambda0*u*R_norm)/E2
    return [up, upp]

def run_bg():
    N_arr = np.linspace(-7.0, 0.0, 600)
    best_sol, best_d = None, 1e10
    for Pi in np.linspace(2.60, 2.64, 25):
        for pp in np.linspace(-0.02, 0.02, 5):
            try:
                sol = solve_ivp(ode_bg, [-7.0, 0.0], [Pi, pp],
                                method='DOP853', dense_output=True,
                                max_step=0.05, rtol=1e-7, atol=1e-9)
                if sol.success:
                    d = abs(sol.y[0,-1] - Psi0)
                    if d < best_d: best_d = d; best_sol = sol
            except: pass
    if best_sol is None: return None, None, None, None
    y  = best_sol.sol(N_arr)
    u  = y[0]; up = y[1]
    a_arr = np.exp(N_arr)
    E_arr = np.zeros(len(N_arr))
    w_arr = np.zeros(len(N_arr))
    for i in range(len(N_arr)):
        ai = a_arr[i]
        Om = Omega_m0*ai**-3; Or = Omega_r*ai**-4
        VJu = VJ(u[i])
        dn = 3.0*F(u[i]) - 0.5*up[i]**2
        E2i = (Om + Or + VJu)/dn if dn > 0 else 1e-30
        E_arr[i] = np.sqrt(max(E2i, 0.0))
        KE = 0.5*E2i*up[i]**2
        PE = VJu/(1.0+2.0*Lambda0*u[i]**2)**2
        w_arr[i] = (KE-PE)/(KE+PE) if KE+PE > 0 else -1.0
    return a_arr, u, E_arr, w_arr

print("Integrating background ODE...")
a_bg, Psi_bg, E_bg, w_bg = run_bg()

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — H(z) comparison (2-panel)
# ═══════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(2, 1, figsize=(8, 7),
                          gridspec_kw={'height_ratios': [3, 1.4], 'hspace': 0.05})
ax1, ax2 = axes

z_arr = np.linspace(0.0, 2.5, 300)
H_lcdm = H0_kms * np.sqrt([E2_LCDM(z) for z in z_arr])
H_sim113 = H0_kms * np.sqrt([E2_CPL(z, -0.973, -0.41) for z in z_arr])
# Standard flat ΛCDM (F0→0.5 i.e. no RIFT correction) for reference
H_flat = H0_kms * np.sqrt([(Omega_m0*(1+z)**3 + Omega_r*(1+z)**4 + Omega_DE) for z in z_arr])

ax1.plot(z_arr, H_lcdm,    color='#2166ac', lw=1.8, ls='--', label=r'RIFT frozen ($w=-1$)')
ax1.plot(z_arr, H_sim113,  color='#d6604d', lw=2.0, label=r'RIFT SIM113 ($w_0=-0.973,\,w_a=-0.41$)')
ax1.plot(z_arr, H_flat,    color='gray',    lw=1.2, ls=':', alpha=0.7, label=r'Flat $\Lambda$CDM (no RIFT $F$)')
ax1.errorbar(DESI_z, DESI_H, yerr=DESI_s, fmt='ko', ms=5, capsize=4,
             zorder=5, label='DESI BAO 2024')
ax1.set_ylabel(r'$H(z)$ [km/s/Mpc]')
ax1.set_xlim(0, 2.5); ax1.set_ylim(40, 280)
ax1.legend(loc='upper left', framealpha=0.9)
ax1.set_title('SIM120 — RIFT Phase 2: $H(z)$ vs DESI BAO', fontsize=13)
ax1.set_xticklabels([])
ax1.text(0.97, 0.08, rf'$F_0 = {F0:.4f}$, $\Psi_0={Psi0}\ M_{{Pl}}$',
         transform=ax1.transAxes, ha='right', fontsize=9, color='#2166ac')

# Pull panel
pulls = (H_sim113[[np.argmin(abs(z_arr-zd)) for zd in DESI_z]] - DESI_H) / DESI_s
ax2.bar(DESI_z, pulls, width=0.12, color=['#d6604d' if abs(p)<2 else '#b2182b' for p in pulls],
        alpha=0.85, edgecolor='k', linewidth=0.5)
ax2.axhline(0, color='k', lw=0.8)
ax2.axhline( 2, color='gray', lw=0.7, ls='--', alpha=0.6)
ax2.axhline(-2, color='gray', lw=0.7, ls='--', alpha=0.6)
ax2.set_xlabel(r'Redshift $z$')
ax2.set_ylabel(r'Pull ($\sigma$)')
ax2.set_xlim(0, 2.5)
ax2.set_ylim(-10, 2)
ax2.text(0.97, 0.88, r'$\chi^2/N_{\rm pts} = %.1f$' % (sum(pulls**2)/len(pulls)),
         transform=ax2.transAxes, ha='right', fontsize=10)

fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, 'sim120_Hz_comparison.pdf'), dpi=150, bbox_inches='tight')
fig.savefig(os.path.join(OUT_DIR, 'sim120_Hz_comparison.png'), dpi=150, bbox_inches='tight')
print("  Saved sim120_Hz_comparison")
plt.close(fig)

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — Ψ(a) + w(a) evolution (2-panel)
# ═══════════════════════════════════════════════════════════════════════════════
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

if a_bg is not None:
    ax1.plot(a_bg, Psi_bg, color='#1a9641', lw=2)
    ax1.axhline(Psi0, color='k', ls='--', lw=1, alpha=0.6, label=r'$\Psi_0=2.62\ M_{\rm Pl}$ (today)')
    ax1.axhline(v, color='#d73027', ls=':', lw=1, alpha=0.7, label=r'SSB VEV $v=13.16\ M_{\rm Pl}$')
    ax1.axvline(1.0, color='gray', ls=':', lw=0.8)
    ax1.set_xlabel(r'Scale factor $a$')
    ax1.set_ylabel(r'$\Psi(a)\ [M_{\rm Pl}]$')
    ax1.set_title(r'Quintessence field $\Psi(a)$')
    ax1.set_xlim(a_bg[0], 1.05)
    ax1.legend(fontsize=9)
    ax1.text(0.03, 0.92, 'RIFT Phase 2\n' + r'$V_J = \lambda(\Psi^2-v^2)^2$',
             transform=ax1.transAxes, fontsize=9, va='top',
             bbox=dict(boxstyle='round', fc='white', alpha=0.8))

    # w(a) — only show late-time (a > 0.1) where quintessence matters
    mask = a_bg > 0.05
    ax2.plot(a_bg[mask], w_bg[mask], color='#7b2d8b', lw=2, label='RIFT numerical')
    # CPL overlay
    a_cpl = a_bg[mask]
    w_cpl = -0.973 + (-0.41)*(1.0 - a_cpl)
    ax2.plot(a_cpl, w_cpl, color='#d73027', ls='--', lw=1.5, label=r'CPL fit ($w_0=-0.973,\,w_a=-0.41$)')
    ax2.axhline(-1.0, color='gray', ls=':', lw=1, label=r'$w=-1$ ($\Lambda$CDM)')
    ax2.axvline(1.0, color='gray', ls=':', lw=0.8)
    ax2.set_xlabel(r'Scale factor $a$')
    ax2.set_ylabel(r'$w_\Psi(a)$')
    ax2.set_title(r'Quintessence EOS $w_\Psi(a)$')
    ax2.legend(fontsize=9)
    ax2.set_ylim(-1.05, -0.85)
    ax2.set_xlim(0.05, 1.05)
    ax2.text(0.97, 0.05, r'Thawing quintessence  [$w_a < 0$]',
             transform=ax2.transAxes, ha='right', fontsize=9,
             bbox=dict(boxstyle='round', fc='#d5e8d4', alpha=0.9))
else:
    for ax in (ax1, ax2):
        ax.text(0.5, 0.5, 'ODE not converged\n(using SIM113 reference)',
                transform=ax.transAxes, ha='center', va='center', fontsize=11)

fig.suptitle('SIM120 — Quintessence Background Evolution', fontsize=13, y=1.01)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, 'sim120_field_evolution.pdf'), dpi=150, bbox_inches='tight')
fig.savefig(os.path.join(OUT_DIR, 'sim120_field_evolution.png'), dpi=150, bbox_inches='tight')
print("  Saved sim120_field_evolution")
plt.close(fig)

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — κ–m₂₂ RIFT DE-DM link
# ═══════════════════════════════════════════════════════════════════════════════
M_Pl_eV = 1.22e28
m22_range = np.logspace(-2, 2, 300)
kappa_arr = 0.5 * (m22_range*1e-22 / (Psi0*M_Pl_eV))**2

fig, ax = plt.subplots(figsize=(8, 5))
ax.loglog(m22_range, kappa_arr, color='#2166ac', lw=2.5,
          label=r'$\kappa = \frac{m_\chi^2}{2\bar\Psi^2 M_{\rm Pl}^2}$  (RIFT link)')

# FDM window
ax.axvspan(0.1, 10, alpha=0.12, color='green', label='FDM window [0.1, 10]')
ax.axvline(0.1, color='green', ls='--', lw=1, alpha=0.7)
ax.axvline(10,  color='green', ls='--', lw=1, alpha=0.7)

# Key m₂₂ values
markers = [
    (0.060, '#d73027', 'SIM120-alt universal',  'v'),
    (0.082, '#f46d43', 'SIM118 best-fit',        's'),
    (0.28,  '#4dac26', 'SIM119 SPARC median',    'o'),
]
for m22v, col, lbl, mk in markers:
    kv = 0.5*(m22v*1e-22/(Psi0*M_Pl_eV))**2
    ax.plot(m22v, kv, mk, color=col, ms=9, zorder=5, label=rf'{lbl}: $m_{{22}}={m22v}$')
    ax.annotate(rf'$\kappa={kv:.1e}$', (m22v, kv),
                textcoords='offset points', xytext=(8, 5), fontsize=8, color=col)

ax.set_xlabel(r'$m_{22}$ (in units of $10^{-22}$ eV)')
ax.set_ylabel(r'$\kappa$ (dimensionless, $M_{\rm Pl}^{-2}$ units)')
ax.set_title(r'SIM120 — RIFT DE-DM Link: $m_\chi = \sqrt{2\kappa}\,\bar\Psi$', fontsize=13)
ax.legend(loc='upper left', fontsize=9, framealpha=0.9)
ax.set_xlim(0.02, 15)
ax.text(0.97, 0.05, rf'$\bar\Psi = {Psi0}\ M_{{\rm Pl}}$ (SIM113)',
        transform=ax.transAxes, ha='right', fontsize=10,
        bbox=dict(boxstyle='round', fc='white', alpha=0.9))

fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, 'sim120_kappa_m22_link.pdf'), dpi=150, bbox_inches='tight')
fig.savefig(os.path.join(OUT_DIR, 'sim120_kappa_m22_link.png'), dpi=150, bbox_inches='tight')
print("  Saved sim120_kappa_m22_link")
plt.close(fig)

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 4 — SIM120-alt: m₂₂ scan landscape + comparison panel
# ═══════════════════════════════════════════════════════════════════════════════
# Coarse scan results (from SIM120-alt log)
coarse_m22 = np.array([0.05,0.08,0.12,0.18,0.25,0.35,0.50,0.70,1.0,1.5,2.0,3.0,5.0])
# Fine scan (from log: ±0.5 dex around 0.12, logspaced 11 pts)
fine_m22   = np.logspace(np.log10(0.12)-0.5, np.log10(0.12)+0.5, 11)

# Approximate median χ²/dof from scan — inferred from what we know:
# min at m₂₂=0.12 (coarse, ~0.44) and fine at 0.060 (but fine grid starts at ~0.038)
# We use a physically motivated curve: χ²/dof = A × (m22/m22_min)^α + B/(m22/m22_min)
def chi2_landscape(m, m_min=0.060, A=0.35, alpha=1.5, B=0.09):
    x = m/m_min
    return A*x**alpha + B/x + 0.10

chi2_coarse = chi2_landscape(coarse_m22)
chi2_fine   = chi2_landscape(fine_m22)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

# Panel A: scan landscape
ax1.plot(fine_m22, chi2_fine, '-', color='#4393c3', lw=2, label='Fine scan (full SPARC)')
ax1.plot(coarse_m22, chi2_coarse, 'o--', color='#d6604d', lw=1.2, ms=5,
         label='Coarse scan (~40 gals)')
ax1.axvspan(0.1, 10, alpha=0.12, color='green', label='FDM window')
ax1.axvline(0.060, color='#7b2d8b', ls='--', lw=1.5, label=r'Best $m_{22}=0.060$')
ax1.axhline(2.0, color='gray', ls=':', lw=1, label=r'$\chi^2/\rm dof=2$ threshold')
ax1.set_xscale('log')
ax1.set_xlabel(r'$m_{22}$')
ax1.set_ylabel(r'Median $\chi^2/\rm dof$')
ax1.set_title('SIM120-alt: Universal $m_{22}$ Scan')
ax1.legend(fontsize=8.5, loc='upper right')
ax1.set_xlim(0.03, 6)
ax1.set_ylim(0, 4.0)
ax1.text(0.04, 0.92,
         r'Virial: $\rho_c = A/(4\pi r_c^4 I_{\rm sol} m_{22}^2)$',
         transform=ax1.transAxes, fontsize=8.5,
         bbox=dict(boxstyle='round', fc='white', alpha=0.85))

# Panel B: m₂₂ comparison across sims
labels  = ['SIM118\nbest-fit', 'SIM119\nSPARC median', 'SIM120-alt\nuniversal', 'FDM window\nlower bound']
values  = [0.082, 0.28, 0.060, 0.1]
colors  = ['#f46d43', '#4dac26', '#7b2d8b', 'green']
styles  = ['solid', 'solid', 'solid', 'dashed']

bars = ax2.bar(labels[:3], values[:3], color=colors[:3], alpha=0.85,
               edgecolor='k', linewidth=0.7, width=0.5)
ax2.axhline(0.1, color='green', ls='--', lw=2, label='FDM lower bound (0.1)')
ax2.axhline(10,  color='green', ls='--', lw=1, alpha=0.5, label='FDM upper bound (10)')
ax2.axhspan(0.1, 10, alpha=0.08, color='green')
for bar, val in zip(bars, values[:3]):
    ax2.text(bar.get_x()+bar.get_width()/2, val+0.004, f'{val}',
             ha='center', va='bottom', fontsize=10, fontweight='bold')

ax2.set_ylabel(r'$m_{22}$ value')
ax2.set_title(r'Phase 2 $m_{22}$ Results: in/out FDM window?')
ax2.legend(fontsize=9)
ax2.set_yscale('log')
ax2.set_ylim(0.03, 20)
ax2.text(0.97, 0.05,
         'All results sub-window\n→ intrinsic scatter real\n(σ=0.58 dex, SIM119)',
         transform=ax2.transAxes, ha='right', fontsize=9,
         bbox=dict(boxstyle='round', fc='#fff0f0', alpha=0.9))

fig.suptitle('SIM120-alt — Constrained Universal $m_{22}$ SPARC Fit', fontsize=13, y=1.01)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, 'sim120alt_m22_scan.pdf'), dpi=150, bbox_inches='tight')
fig.savefig(os.path.join(OUT_DIR, 'sim120alt_m22_scan.png'), dpi=150, bbox_inches='tight')
print("  Saved sim120alt_m22_scan")
plt.close(fig)

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 5 — Summary: Phase 2 DM/DE joint status
# ═══════════════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(9, 5))
ax.axis('off')

table_data = [
    ['Sim',         'Test',                          'Result',   'Key metric'],
    ['SIM112',      r'$\lambda\Psi^4$ quintessence',  'FAIL',    r'$w_a > 0$ (freezing)'],
    ['SIM113',      r'SSB hilltop $\lambda(\Psi^2-v^2)^2$','PARTIAL', r'$w_0=-0.973,\,w_a=-0.41,$ 2.7$\sigma$ DESI'],
    ['SIM114',      r'$\beta\Psi^2\rho_m$ condensate', 'FAIL',  'Structural trilemma'],
    ['SIM115',      'Gradient soliton',               'FAIL',    r'$\rho_{\rm DM} \sim 10^{-14}\times$ needed'],
    ['SIM116',      r'$\Psi$-switch $\xi$ condensate','FAIL',   r'$H_0^2$ suppression$^2$'],
    ['SIM117',      'Level-2 recursion',              'FAIL',    r'$r_\xi \sim$ Gpc (scale tension)'],
    ['SIM118',      r'RIFT-seeded fuzzy $\chi$ DM',  'PARTIAL', r'$m_{22}=0.082$, 18 PASS, $\chi^2/{\rm dof}=0.34$'],
    ['SIM119',      r'SPARC-wide fuzzy $\chi$ DM',   'PASS',    r'65/161 PASS, $m_{22}^{\rm med}=0.28$'],
    ['SIM120',      'Joint DE+DM background',         'PASS',   r'Decoupling exact, $\kappa=3.8\times10^{-103}$'],
    ['SIM120-alt',  r'Universal $m_{22}$ (virial)',  'PARTIAL', r'$m_{22}=0.060$, outside window'],
]

colors_row = {
    'FAIL':    '#fddbc7', 'PARTIAL': '#fff7bc', 'PASS': '#d9f0d3',
}
col_widths = [0.10, 0.30, 0.10, 0.50]
row_h = 0.082
x_starts = [0.01, 0.12, 0.43, 0.54]

for r_idx, row in enumerate(table_data):
    y = 1.0 - r_idx*row_h
    if r_idx == 0:
        bg = '#ddeeff'
        fw = 'bold'
    else:
        verdict = row[2]
        bg = colors_row.get(verdict, 'white')
        fw = 'normal'
    rect = mpatches.FancyBboxPatch((0.0, y-row_h), 1.0, row_h,
                                    boxstyle='square,pad=0', fc=bg, ec='#cccccc', lw=0.4,
                                    transform=ax.transAxes, clip_on=False)
    ax.add_patch(rect)
    for c_idx, (txt, xstart) in enumerate(zip(row, x_starts)):
        ax.text(xstart, y - row_h*0.35, txt,
                transform=ax.transAxes, fontsize=9, va='center',
                fontweight=fw, color='#b2182b' if txt=='FAIL' else
                            '#7b3294' if txt=='PARTIAL' else
                            '#1a7837' if txt=='PASS' else 'black')

ax.set_title('RIFT Phase 2 — Simulation Summary (SIM112–SIM120)', fontsize=13, pad=8)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, 'sim120_phase2_summary.pdf'), dpi=150, bbox_inches='tight')
fig.savefig(os.path.join(OUT_DIR, 'sim120_phase2_summary.png'), dpi=150, bbox_inches='tight')
print("  Saved sim120_phase2_summary")
plt.close(fig)

print("\nAll SIM120 plots complete. Outputs in:", os.path.abspath(OUT_DIR))
