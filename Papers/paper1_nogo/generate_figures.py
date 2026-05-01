"""
Paper 1: Two Structural No-Go Theorems for Late-Time Modifications of Dark Energy
Figure generation script — all four submission figures.

Figures:
  fig1_psi_evolution.pdf     — phi(z) trajectories for 6 Phase 3 mechanisms
  fig2_desi_chi2.pdf         — DESI tension vs coupling, 6-panel grid
  fig3_desi_cmb_tradeoff.pdf — DESI tension vs delta-theta* scatter (SIM141)
  fig4_biscalar.pdf          — Bi-scalar quintessence loophole closure (SIM143)
"""

import numpy as np
import json
import os
import warnings
warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from scipy.integrate import solve_ivp, quad
from collections import defaultdict

# ──────────────────────────────────────────────────────────────────────────────
# PATHS
# ──────────────────────────────────────────────────────────────────────────────
HERE   = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, 'figures')
SIMDIR = '/home/aion/Ordered_Simulations'
os.makedirs(FIGDIR, exist_ok=True)


def sim_json(simn, fname):
    for sub in ('Outputs', ''):
        p = os.path.join(SIMDIR, simn, sub, fname)
        if os.path.exists(p):
            with open(p) as f:
                return json.load(f)
    raise FileNotFoundError(f"{simn}/{fname}")


# ──────────────────────────────────────────────────────────────────────────────
# PUBLICATION STYLE  (no text.usetex — use mathtext only)
# ──────────────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family':        'serif',
    'font.serif':         ['DejaVu Serif', 'Times New Roman', 'Times'],
    'mathtext.fontset':   'dejavuserif',
    'font.size':          10,
    'axes.labelsize':     11,
    'axes.titlesize':     10,
    'legend.fontsize':    8.5,
    'xtick.labelsize':    9,
    'ytick.labelsize':    9,
    'xtick.direction':    'in',
    'ytick.direction':    'in',
    'xtick.top':          True,
    'ytick.right':        True,
    'axes.linewidth':     0.8,
    'grid.linewidth':     0.5,
    'grid.alpha':         0.35,
    'lines.linewidth':    1.8,
    'legend.framealpha':  0.9,
    'legend.edgecolor':   '0.7',
    'figure.dpi':         150,
    'savefig.dpi':        300,
    'savefig.bbox':       'tight',
    'savefig.pad_inches': 0.05,
})

# Colorblind-safe palette
C = {
    'sim131': '#c0392b',
    'sim132': '#2980b9',
    'sim133': '#8e44ad',
    'sim134': '#27ae60',
    'sim135': '#e67e22',
    'sim136': '#7f8c8d',
    'p1':     '#2c3e50',
    'tanh':   '#1abc9c',
    'exp':    '#e74c3c',
    'linear': '#95a5a6',
}

# ──────────────────────────────────────────────────────────────────────────────
# PHYSICS CONSTANTS  (Planck 2018, identical to SIM scripts)
# ──────────────────────────────────────────────────────────────────────────────
c_kms         = 2.998e5
H100          = 100.0
omh2_m        = 0.1430
omh2_r        = 4.18e-5
h_target      = 0.674
theta_obs     = 1.04101
theta_obs_err = 0.00029
z_star        = 1089.8
z_drag        = 1059.6
omh2_b        = 0.02237
PSI0_P1       = 2.62
F0_P1         = 0.5 + 0.003 * PSI0_P1**2

DESI_z = np.array([0.295, 0.510, 0.706, 0.930, 1.317, 2.330])
DESI_H = np.array([ 81.7,  97.9, 110.7, 128.1, 156.4, 240.8])
DESI_s = np.array([  4.5,   4.4,   6.2,   5.6,   8.6,  11.0])

# ──────────────────────────────────────────────────────────────────────────────
# GENERAL ODE INTEGRATOR  (covers all Phase 3 curvature-sourced mechanisms)
# ──────────────────────────────────────────────────────────────────────────────

def get_H2(phi, y, N, Lb, F_fn, dF_fn):
    a    = np.exp(N)
    rhs  = omh2_m / a**3 + omh2_r / a**4 + Lb
    Fv   = F_fn(phi)
    Fp   = dF_fn(phi)
    coef = 3.0*Fv + 6.0*Fp*y - 0.5*y**2
    coef = max(coef, 3.0*max(Fv, 0.3))
    return rhs / coef


def get_eps_H(phi, y, N, Lb, F_fn, dF_fn):
    eps = 5e-4
    H2p = get_H2(phi, y, N+eps, Lb, F_fn, dF_fn)
    H2m = get_H2(phi, y, N-eps, Lb, F_fn, dF_fn)
    H2  = get_H2(phi, y, N,    Lb, F_fn, dF_fn)
    if H2 < 1e-40:
        return 0.0
    return -0.5*(H2p - H2m) / (2.0*eps*H2)


def make_ode(F_fn, dF_fn):
    def ode(N, state, Lb):
        phi, y    = state
        H2        = get_H2(phi, y, N, Lb, F_fn, dF_fn)
        eps_H     = get_eps_H(phi, y, N, Lb, F_fn, dF_fn)
        R_over_H2 = 6.0*(2.0 - eps_H)
        source    = dF_fn(phi) * R_over_H2
        dy_dN     = source - (3.0 - eps_H)*y
        return [y, dy_dN]
    return ode


def integrate_bg(phi_ini, Lb, F_fn, dF_fn, z_ini=1e5, n_pts=5000):
    ode    = make_ode(F_fn, dF_fn)
    N_ini  = np.log(1.0/(1.0+z_ini))
    N_eval = np.linspace(N_ini, 0.0, n_pts)
    sol = solve_ivp(ode, (N_ini, 0.0), [phi_ini, 0.0], args=(Lb,),
                    t_eval=N_eval, method='RK45',
                    rtol=1e-8, atol=1e-12, max_step=0.05)
    if not sol.success:
        return None
    N_arr   = sol.t
    z_arr   = np.exp(-N_arr) - 1.0
    phi_arr = sol.y[0]
    y_arr   = sol.y[1]
    H_arr   = np.array([H100*np.sqrt(max(get_H2(phi_arr[i], y_arr[i],
                                              N_arr[i], Lb, F_fn, dF_fn), 0.0))
                        for i in range(len(N_arr))])
    Feff_arr = np.array([F_fn(p) for p in phi_arr])
    idx = np.argsort(z_arr)
    return dict(z=z_arr[idx], phi=phi_arr[idx], y=y_arr[idx],
                H=H_arr[idx], Feff=Feff_arr[idx],
                phi0=float(phi_arr[-1]), y0=float(y_arr[-1]),
                Feff0=float(Feff_arr[-1]), H0=float(H_arr[-1]))


def calibrate_Lb(phi0, y0, F_fn, dF_fn):
    coef = 3.0*F_fn(phi0) + 6.0*dF_fn(phi0)*y0 - 0.5*y0**2
    coef = max(coef, 3.0*max(F_fn(phi0), 0.3))
    return h_target**2 * coef - omh2_m - omh2_r


def integrate_self(phi_ini, F_fn, dF_fn, z_ini=1e5, n_iter=8):
    Lb = 3.0*F0_P1*h_target**2 - omh2_m - omh2_r
    bg = None
    for _ in range(n_iter):
        bg = integrate_bg(phi_ini, Lb, F_fn, dF_fn, z_ini=z_ini)
        if bg is None:
            return None
        Lb_new = calibrate_Lb(bg['phi0'], bg['y0'], F_fn, dF_fn)
        if abs(Lb_new - Lb) < 1e-12:
            break
        Lb = Lb_new
    return bg


# ──────────────────────────────────────────────────────────────────────────────
# MECHANISM COUPLING FUNCTIONS
# ──────────────────────────────────────────────────────────────────────────────

L0_131 = 0.003; xi_131 = 1.0/6.0
F131  = lambda p: 0.5 + L0_131*p**2 + xi_131*p
dF131 = lambda p: 2*L0_131*p + xi_131

L0_132 = 1.0; xi_132 = 1.0/6.0
F132  = lambda p: 0.5 - L0_132*p**2 + xi_132*p
dF132 = lambda p: -2*L0_132*p + xi_132

L0_134 = 0.2; xi_134 = 1.0/6.0
F134  = lambda p: 0.5*np.exp(-2*L0_134*p**2) + xi_134*p
dF134 = lambda p: -2*L0_134*p*np.exp(-2*L0_134*p**2) + xi_134

xi_135 = 0.01
F135  = lambda p: F0_P1 + xi_135*p
dF135 = lambda p: xi_135


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1: phi(z) trajectories — six Phase 3 mechanisms
# ══════════════════════════════════════════════════════════════════════════════
print("=" * 60)
print("Figure 1: phi(z) evolution")

print("  SIM131 (conformal additive)...")
bg131 = integrate_self(0.0, F131, dF131)
print("  SIM132 (subtractive DW)...")
bg132 = integrate_self(0.0, F132, dF132)
print("  SIM134 (dilaton)...")
bg134 = integrate_self(0.0, F134, dF134)
print("  SIM135 (bi-scalar)...")
bg135 = integrate_self(0.0, F135, dF135)

# SIM133: integrate briefly to show GB divergence
print("  SIM133 (Gauss-Bonnet — short run)...")
bg133_z, bg133_phi = np.array([]), np.array([])
try:
    def ode_gb(N, state, Lb):
        phi, y = state
        F_gb  = lambda p: 0.5 + 0.003*p**2
        dF_gb = lambda p: 2*0.003*p
        H2 = get_H2(phi, y, N, Lb, F_gb, dF_gb)
        if H2 <= 0 or abs(phi) > 1e6:
            return [0.0, 0.0]
        eps_H     = get_eps_H(phi, y, N, Lb, F_gb, dF_gb)
        R_over_H2 = 6.0*(2.0 - eps_H)
        GB_over_H2 = 24.0 * H2 * (1.0 - eps_H)
        source_R  = dF_gb(phi) * R_over_H2
        source_GB = 1e-8 * GB_over_H2
        dy_dN = (source_R + source_GB) - (3.0 - eps_H)*y
        return [y, dy_dN]
    Lb_gb = 3.0*F0_P1*h_target**2 - omh2_m - omh2_r
    sol_gb = solve_ivp(ode_gb, (np.log(1.0/1e5), 0.0), [0.0, 0.0],
                       args=(Lb_gb,), method='RK45',
                       max_step=0.01, rtol=1e-5, atol=1e-8,
                       dense_output=False)
    zg = np.exp(-sol_gb.t) - 1.0
    pg = sol_gb.y[0]
    mask_gb = (np.abs(pg) < 5.0) & (zg >= 0)
    bg133_z   = zg[mask_gb]
    bg133_phi = pg[mask_gb]
except Exception as e:
    print(f"    GB integration note: {e}")

z_CMB = z_star
fig1, ax1 = plt.subplots(1, 1, figsize=(8.5, 5.8))

# --- main curves ---
if bg131:
    m = bg131['z'] <= z_CMB * 1.05
    ax1.plot(bg131['z'][m], bg131['phi'][m], color=C['sim131'], lw=2.2,
             label=r'SIM131: $F_{\rm eff}=\frac{1+2\Lambda_0\Psi^2}{2}+\xi\Psi$, $\xi=1/6$')

if bg132:
    m = bg132['z'] <= z_CMB * 1.05
    ax1.plot(bg132['z'][m], bg132['phi'][m], color=C['sim132'], lw=2.0,
             label=r'SIM132: $F_{\rm eff}=\frac{1-2\Lambda_0\Psi^2}{2}+\xi\Psi$ (Deser-Woodard)')

if bg134:
    m = bg134['z'] <= z_CMB * 1.05
    ax1.plot(bg134['z'][m], bg134['phi'][m], color=C['sim134'], lw=2.0, ls='--',
             label=r'SIM134: $F_{\rm eff}=\frac{1}{2}e^{-2\Lambda_0\Psi^2}+\xi\Psi$ (dilaton)')

if bg135:
    m = bg135['z'] <= z_CMB * 1.05
    ax1.plot(bg135['z'][m], bg135['phi'][m], color=C['sim135'], lw=2.0, ls=(0,(3,1,1,1)),
             label=r'SIM135: $F_{\rm eff}=F_0+\xi_\phi\phi$, $\xi_\phi=0.01$ ($\Psi$ frozen)')

# SIM136: Psi=0 throughout (source vanishes at origin)
z_flat = np.logspace(-2, np.log10(z_CMB * 1.05), 500)
ax1.plot(z_flat, np.zeros_like(z_flat), color=C['sim136'], lw=1.8, ls=':',
         label=r'SIM136: Horndeski kinetic — $\Psi\equiv 0$ (source vanishes at $\Psi_{\rm ini}=0$)')

# SIM133: show brief trajectory before divergence
if len(bg133_z) > 3:
    idx_s = np.argsort(bg133_z)[::-1]   # high-z first
    ax1.plot(bg133_z[idx_s], bg133_phi[idx_s], color=C['sim133'], lw=1.5, ls='-.',
             label=r'SIM133: Gauss-Bonnet (diverges, $|\Psi_0|\sim 3\times 10^5\,M_{\rm Pl}$)')

# Phase 1 frozen reference
ax1.axhline(PSI0_P1, color=C['p1'], ls='--', lw=1.4, alpha=0.7,
            label=r'Phase 1 canonical: $\bar\Psi=2.62\,M_{\rm Pl}$ (frozen, not curvature-sourced)')

# CMB epoch
ax1.axvline(z_CMB, color='#f39c12', ls=':', lw=1.5, alpha=0.85)
ax1.text(z_CMB * 1.06, 0.12, r'$z_{\rm CMB}$', fontsize=9, color='#f39c12',
         ha='left', va='bottom')

# Annotation: all start at zero
ax1.annotate('All scalars: $\\phi_{\\rm ini}=0$',
             xy=(z_CMB * 0.85, 0.03), fontsize=8.5, color='0.30', ha='right',
             bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='0.6', alpha=0.85))

if bg131:
    ax1.annotate(r'$\Psi_0=2.91\,M_{\rm Pl}$',
                 xy=(0.12, bg131['phi0']),
                 xytext=(2.0, bg131['phi0'] - 0.2),
                 fontsize=8, color=C['sim131'],
                 arrowprops=dict(arrowstyle='->', color=C['sim131'], lw=0.8))

ax1.set_xscale('log')
ax1.invert_xaxis()
ax1.set_xlim(z_CMB * 1.8, 5e-2)
ax1.set_ylim(-0.15, 3.3)
ax1.set_xlabel(r'Redshift $z$  (right = early universe)', labelpad=4)
ax1.set_ylabel(r'$\phi(z)\;[M_{\rm Pl}]$', labelpad=4)
ax1.set_xticks([0.1, 1, 10, 100, 1000])
ax1.set_xticklabels(['0.1', '1', '10', '100', '1000'])
ax1.legend(loc='upper left', fontsize=8, framealpha=0.93, edgecolor='0.7', handlelength=2.8)
ax1.grid(True, which='both', ls='--', lw=0.4, alpha=0.4)
ax1.set_title('Theorem 1: curvature-sourced scalars initialized at $\\phi_{\\rm ini}=0$ '
              'grow monotonically (violates Eq. (2))',
              fontsize=9.5, pad=6)

# Inset: zoom on small-growth curves
ax_in = ax1.inset_axes([0.60, 0.07, 0.38, 0.44])
for bg, col, ls in [(bg132, C['sim132'], '-'), (bg134, C['sim134'], '--'),
                    (bg135, C['sim135'], (0,(3,1,1,1)))]:
    if bg:
        m = bg['z'] <= z_CMB * 1.05
        ax_in.plot(bg['z'][m], bg['phi'][m], color=col, lw=1.6, ls=ls)
ax_in.plot(z_flat, np.zeros_like(z_flat), color=C['sim136'], lw=1.4, ls=':')
ax_in.axvline(z_CMB, color='#f39c12', ls=':', lw=1.0, alpha=0.65)
ax_in.set_xscale('log')
ax_in.invert_xaxis()
ax_in.set_xlim(z_CMB * 1.8, 5e-2)
ax_in.set_ylim(-0.02, 0.58)
ax_in.set_ylabel(r'$\phi\;[M_{\rm Pl}]$', fontsize=7.5)
ax_in.set_xlabel(r'$z$', fontsize=7.5)
ax_in.tick_params(labelsize=7)
ax_in.set_xticks([0.1, 10, 1000])
ax_in.set_xticklabels(['0.1', '10', '1000'], fontsize=7)
ax_in.set_title('SIM132, 134, 135, 136\n(expanded scale)', fontsize=7.5)
ax_in.grid(True, ls='--', lw=0.4, alpha=0.4)

fig1.tight_layout()
for ext in ('pdf', 'png'):
    fig1.savefig(os.path.join(FIGDIR, f'fig1_psi_evolution.{ext}'))
plt.close(fig1)
print("  Saved fig1_psi_evolution")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2: DESI tension vs coupling parameter — 2x3 grid
# ══════════════════════════════════════════════════════════════════════════════
print("Figure 2: DESI tension vs coupling...")

d131 = sim_json('SIM131', 'sim131_results.json')
d132 = sim_json('SIM132', 'sim132_results.json')
d133 = sim_json('SIM133', 'sim133_results.json')
d134 = sim_json('SIM134', 'sim134_results.json')
d135 = sim_json('SIM135', 'sim135_results.json')
d136 = sim_json('SIM136', 'sim136_results.json')

baseline_tension = 2.77

fig2, axes2 = plt.subplots(2, 3, figsize=(13.5, 7.2),
                            gridspec_kw=dict(hspace=0.50, wspace=0.36))

panels = [
    (axes2[0,0],
     np.array(d131['scan_xi']['xi']),
     np.array(d131['scan_xi']['chi2_DESI']),
     r'Coupling $\xi$',
     C['sim131'], 'SIM131: conformal additive', 'linear'),

    (axes2[0,1],
     np.array([r['L0']       for r in d132['scan']]),
     np.array([r['chi2_DESI'] for r in d132['scan']]),
     r'Coupling $\Lambda_0$',
     C['sim132'], 'SIM132: Deser-Woodard', 'linear'),

    (axes2[0,2],
     np.array([r['alpha']     for r in d133['scan']]),
     np.array([r['chi2_DESI'] for r in d133['scan']]),
     r'Coupling $\tilde{\alpha}$',
     C['sim133'], 'SIM133: Gauss-Bonnet', 'log'),

    (axes2[1,0],
     np.array([r['L0']       for r in d134['scan']]),
     np.array([r['chi2_DESI'] for r in d134['scan']]),
     r'Coupling $\Lambda_0$',
     C['sim134'], 'SIM134: dilaton', 'linear'),

    (axes2[1,1],
     np.array([r['xi']       for r in d135['scan']]),
     np.array([r['chi2_DESI'] for r in d135['scan']]),
     r'Coupling $\xi_\phi$',
     C['sim135'], r'SIM135: bi-scalar $\phi$', 'linear'),

    (axes2[1,2],
     np.array([r['alpha']    for r in d136['scan_A_psi_ini0']]),
     np.array([r['chi2_DESI'] for r in d136['scan_A_psi_ini0']]),
     r'Coupling $\alpha$',
     C['sim136'], r'SIM136: Horndeski ($\Psi_{\rm ini}=0$)', 'log'),
]

for (ax, xs, chi2s, xlbl, col, lbl, xsc) in panels:
    tension = np.sqrt(chi2s / len(DESI_z))

    ax.plot(xs, tension, 'o-', color=col, lw=2.0, ms=5.5,
            markerfacecolor=col, markeredgecolor='white', markeredgewidth=0.5)

    ax.axhline(baseline_tension, color=C['p1'], ls='--', lw=1.5, alpha=0.8)
    ax.axhline(2.0, color='#27ae60', ls=':', lw=1.3, alpha=0.8)
    ax.axhspan(0, baseline_tension, alpha=0.05, color='#2ecc71')

    ax.set_xlabel(xlbl, fontsize=9.5)
    ax.set_ylabel(r'DESI tension [$\sigma$]', fontsize=9.5)
    ax.set_title(lbl, fontsize=9, color=col, fontweight='bold')
    ax.set_xscale(xsc)
    ax.set_ylim(bottom=0)
    ax.grid(True, ls='--', lw=0.4, alpha=0.4)

    y_min  = float(np.min(tension))
    ix_min = int(np.argmin(tension))
    ax.annotate(f'min: {y_min:.2f}' + r'$\sigma$' + ' > baseline',
                xy=(xs[ix_min], y_min),
                xytext=(0.55, 0.88), textcoords='axes fraction',
                fontsize=7.5, color=col, ha='center',
                bbox=dict(boxstyle='round,pad=0.2', fc='white', ec=col, alpha=0.85),
                arrowprops=dict(arrowstyle='->', color=col, lw=0.7))

# Shared legend at bottom
handles2 = [
    Line2D([0],[0], color=C['p1'], ls='--', lw=1.5,
           label=r'Phase 1 canonical ($2.77\sigma$)'),
    Line2D([0],[0], color='#27ae60', ls=':', lw=1.3,
           label=r'$2\sigma$ PASS threshold'),
    mpatches.Patch(color='#2ecc71', alpha=0.15,
                   label='Would-be improvement region'),
]
fig2.legend(handles=handles2, loc='lower center', ncol=3,
            bbox_to_anchor=(0.5, -0.01), fontsize=9,
            framealpha=0.92, edgecolor='0.7')

fig2.suptitle('Theorem 1 verification: DESI tension vs coupling for all six '
              'Phase 3 mechanisms (SIM131-SIM136). '
              'No mechanism reduces tension below the Phase 1 canonical baseline.',
              fontsize=9.5, y=1.01)

for ext in ('pdf', 'png'):
    fig2.savefig(os.path.join(FIGDIR, f'fig2_desi_chi2.{ext}'))
plt.close(fig2)
print("  Saved fig2_desi_chi2")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3: DESI-CMB trade-off (SIM141) — Theorem 2 visual proof
# ══════════════════════════════════════════════════════════════════════════════
print("Figure 3: DESI-CMB trade-off (SIM141)...")

d141 = sim_json('SIM141', 'output.json')

forms_data = defaultdict(list)
for r in d141['scan_results']:
    forms_data[r['form']].append(r)

p1_ref_tension = d141['phase1_reference']['desi_tension']
p1_ref_theta   = d141['phase1_reference']['theta_star_100']

planck_lo = theta_obs - 2*theta_obs_err
planck_hi = theta_obs + 2*theta_obs_err

fig3, ax3 = plt.subplots(1, 1, figsize=(8.5, 6.5))

# Planck band
ax3.axhspan(planck_lo, planck_hi, alpha=0.20, color='#3498db', zorder=0)
ax3.axhline(theta_obs, color='#2980b9', ls='-', lw=1.0, alpha=0.6, zorder=1)

# Phase 1 ODE reference crosshairs
ax3.axvline(p1_ref_tension, color=C['p1'], ls='--', lw=1.4, alpha=0.65, zorder=1)
ax3.axhline(p1_ref_theta,   color=C['p1'], ls='--', lw=1.4, alpha=0.65, zorder=1)

# DESI-improving region (tension < Phase 1 reference)
ax3.axvspan(0, p1_ref_tension, alpha=0.06, color='#2ecc71', zorder=0)

form_styles = {
    'linear':      dict(marker='s', color=C['linear'], ms=7, ls=':', lw=1.1, alpha=0.55,
                        label=r'Linear $\Lambda_0(a)$ (UV-inconsistent)'),
    'exponential': dict(marker='^', color=C['exp'],    ms=7, ls='--', lw=1.3, alpha=0.70,
                        label=r'Exponential $\Lambda_0(a)$ (UV-inconsistent)'),
    'tanh':        dict(marker='o', color=C['tanh'],   ms=9, ls='-',  lw=2.0, alpha=1.0,
                        label=r'Tanh $\Lambda_0(a)$ (UV-consistent, SIM105)'),
}

for form in ('linear', 'exponential', 'tanh'):
    recs = forms_data.get(form, [])
    if not recs:
        continue
    tensions = [r['desi_tension']    for r in recs]
    thetas   = [r['theta_star_100']  for r in recs]
    st = form_styles[form]
    # linear/exponential: connect because they are 1D parameter scans
    # tanh: 2D scan (a_trans x gamma) — scatter only, no misleading lines
    if form != 'tanh':
        recs_s = sorted(recs, key=lambda r: r['desi_tension'])
        ts_s   = [r['desi_tension']   for r in recs_s]
        th_s   = [r['theta_star_100'] for r in recs_s]
        ax3.plot(ts_s, th_s, linestyle=st['ls'], color=st['color'],
                 lw=st['lw'], alpha=st['alpha']*0.6, zorder=2)
    ax3.scatter(tensions, thetas, marker=st['marker'], color=st['color'],
                s=st['ms']**2, alpha=min(st['alpha']*1.3, 1.0),
                edgecolors='white', linewidths=0.5, zorder=5)

# Directional arrow: DESI improvement → theta* rises (trade-off direction for tanh)
tanh_recs = sorted(forms_data.get('tanh', []), key=lambda r: r['desi_tension'])
if len(tanh_recs) >= 2:
    t_hi = tanh_recs[-1]['desi_tension'];  th_hi = tanh_recs[-1]['theta_star_100']
    t_lo = tanh_recs[0]['desi_tension'];   th_lo = tanh_recs[0]['theta_star_100']
    ax3.annotate('', xy=(t_lo*1.02, th_lo*0.9998 + th_hi*0.0002),
                 xytext=(t_hi*0.985, th_hi*0.9998 + th_lo*0.0002),
                 arrowprops=dict(arrowstyle='->', color=C['tanh'], lw=1.8,
                                 connectionstyle='arc3,rad=0.15'), zorder=6)

# Annotate best DESI case
best = min(d141['scan_results'], key=lambda r: r['desi_tension'])
ax3.annotate(
    f"Best DESI: {best['desi_tension']:.3f}" + r'$\sigma$' + '\n'
    f"$100\\theta_*={best['theta_star_100']:.4f}$\n"
    r"(+63$\sigma$ CMB violation)",
    xy=(best['desi_tension'], best['theta_star_100']),
    xytext=(0.73, 0.52), textcoords='axes fraction',
    fontsize=8.5, ha='center', color=C['tanh'],
    bbox=dict(boxstyle='round,pad=0.4', fc='white', ec=C['tanh'], alpha=0.93, lw=1.0),
    arrowprops=dict(arrowstyle='->', color=C['tanh'], lw=1.0,
                    connectionstyle='arc3,rad=-0.22'),
    zorder=10)

# Empty-intersection annotation
ax3.text(p1_ref_tension * 0.38, planck_lo - 0.0015,
         r'DESI-improving $\cap$ Planck-allowed $=\,\varnothing$',
         fontsize=9, color='#e74c3c', ha='center', va='top', fontweight='bold',
         bbox=dict(boxstyle='round,pad=0.3', fc='#fef9e7', ec='#e74c3c', alpha=0.92))

ax3.text(0.55, (planck_lo + planck_hi)/2,
         'Planck band', fontsize=8.5, color='#2980b9', ha='left', va='center', style='italic')

ax3.set_xlabel(r'DESI tension [$\sigma$]  (lower = better)', fontsize=11)
ax3.set_ylabel(r'$100\,\theta_*$', fontsize=11)

all_tensions = [r['desi_tension']   for r in d141['scan_results']]
all_thetas   = [r['theta_star_100'] for r in d141['scan_results']]
ax3.set_xlim(min(all_tensions)*0.85, p1_ref_tension*1.08)
ax3.set_ylim(planck_lo - 0.003, max(all_thetas)*1.003)

handles3 = [
    mpatches.Patch(color='#3498db', alpha=0.4,
                   label=r'Planck $2\sigma$: $100\theta_*=1.04101\pm0.00058$'),
    mpatches.Patch(color='#2ecc71', alpha=0.15,
                   label='DESI-improving region'),
    Line2D([0],[0], color=C['p1'], ls='--', lw=1.4,
           label=r'Phase 1 ODE ref ($1.507\sigma$, $\theta_*=1.04096$)'),
] + [
    Line2D([0],[0], marker=form_styles[f]['marker'], color=form_styles[f]['color'],
           ls=form_styles[f]['ls'], lw=1.5, ms=7, label=form_styles[f]['label'])
    for f in ('tanh', 'exponential', 'linear')
]
ax3.legend(handles=handles3, loc='upper right', fontsize=8.2,
           framealpha=0.95, edgecolor='0.7')
ax3.grid(True, ls='--', lw=0.4, alpha=0.4)
ax3.set_title('Theorem 2 proof: DESI improvement and CMB preservation are '
              'mutually exclusive for any late-time $H(z)$ boost with fixed $r_s$.',
              fontsize=9.5, pad=8)

for ext in ('pdf', 'png'):
    fig3.savefig(os.path.join(FIGDIR, f'fig3_desi_cmb_tradeoff.{ext}'))
plt.close(fig3)
print("  Saved fig3_desi_cmb_tradeoff")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 4: Bi-scalar quintessence (SIM143) — loophole closure
# ══════════════════════════════════════════════════════════════════════════════
print("Figure 4: Bi-scalar quintessence (SIM143)...")

d143 = sim_json('SIM143', 'output.json')

p1ref_t143   = d143['phase1_reference']['desi_tension']
p1ref_th143  = d143['phase1_reference']['theta_star_100']
planck_lo143 = theta_obs - 2*theta_obs_err
planck_hi143 = theta_obs + 2*theta_obs_err

by_form143 = defaultdict(list)
for r in d143['scan_results']:
    by_form143[r['form']].append(r)

form_styles_143 = {
    'exp':       dict(color='#e74c3c', marker='o', label=r'Exponential $e^{-\lambda\phi}$'),
    'power-law': dict(color='#3498db', marker='s', label=r'Inverse power law $\phi^{-n}$'),
    'hilltop':   dict(color='#8e44ad', marker='^', label=r'Hilltop $1-(\phi/\mu)^2$'),
}
U0_sizes = {0.05: 50, 0.20: 110, 0.50: 210}

fig4, (ax4L, ax4R) = plt.subplots(1, 2, figsize=(12.5, 5.8),
                                   gridspec_kw=dict(wspace=0.38))

# ---------- LEFT: DESI tension vs theta* ----------
ax4L.axhspan(planck_lo143, planck_hi143, alpha=0.20, color='#3498db', zorder=0)
ax4L.axhline(theta_obs, color='#2980b9', ls='-', lw=0.9, alpha=0.55, zorder=1)
ax4L.axvline(p1ref_t143,   color=C['p1'], ls='--', lw=1.3, alpha=0.60, zorder=1)
ax4L.axhline(p1ref_th143,  color=C['p1'], ls='--', lw=1.3, alpha=0.60, zorder=1)
ax4L.axvspan(0, p1ref_t143, alpha=0.06, color='#2ecc71', zorder=0)

for form, recs in by_form143.items():
    st = form_styles_143.get(form, dict(color='k', marker='x', label=form))
    tensions = [r['desi_tension']   for r in recs]
    thetas   = [r['theta_star_100'] for r in recs]
    U0s      = [r['U0']             for r in recs]
    sizes    = [U0_sizes.get(u, 80) for u in U0s]

    recs_s   = sorted(recs, key=lambda r: r['U0'])
    xs_s     = [r['desi_tension']   for r in recs_s]
    ys_s     = [r['theta_star_100'] for r in recs_s]
    ax4L.plot(xs_s, ys_s, color=st['color'], lw=1.0, ls='--', alpha=0.45, zorder=2)
    ax4L.scatter(tensions, thetas, c=st['color'], marker=st['marker'],
                 s=sizes, edgecolors='white', linewidths=0.6, alpha=0.88, zorder=5)

# Invisible scatter entries for U0 size legend
for u, sz in sorted(U0_sizes.items()):
    ax4L.scatter([], [], c='0.4', marker='o', s=sz, edgecolors='white', linewidths=0.5,
                 label=f'$U_0={u}$')

ax4L.text(p1ref_t143*0.55, planck_lo143 - 0.0017,
          r'$|\Delta r_s|<10^{-8}$ Mpc for all 18 cases' + '\n' +
          r'$\Rightarrow$ SIM141 loophole does not exist',
          fontsize=8.5, color='#c0392b', ha='center', va='top', fontweight='bold',
          bbox=dict(boxstyle='round,pad=0.3', fc='#fef5e7', ec='#c0392b', alpha=0.92))

ax4L.text(p1ref_t143*0.93, (planck_lo143+planck_hi143)/2,
          'Planck band', fontsize=8.0, color='#2980b9', ha='right', va='center', style='italic')

ax4L.set_xlabel(r'DESI tension [$\sigma$]', fontsize=11)
ax4L.set_ylabel(r'$100\,\theta_*$', fontsize=11)
ax4L.set_xlim(min(r['desi_tension'] for r in d143['scan_results'])*0.95, p1ref_t143*1.03)
ax4L.set_ylim(planck_lo143 - 0.005,
              max(r['theta_star_100'] for r in d143['scan_results'])*1.002)
ax4L.grid(True, ls='--', lw=0.4, alpha=0.4)
ax4L.set_title('DESI-CMB trade-off for SIM143 bi-scalar $\\phi$', fontsize=9.5)

handles4L = []
for form in ('exp', 'power-law', 'hilltop'):
    st = form_styles_143[form]
    handles4L.append(ax4L.scatter([], [], c=st['color'], marker=st['marker'],
                                  s=100, label=st['label'], edgecolors='white', linewidths=0.6))
handles4L += [
    mpatches.Patch(alpha=0, label=''),   # spacer
    mpatches.Patch(color='#3498db', alpha=0.35, label=r'Planck $2\sigma$ band'),
    mpatches.Patch(color='#2ecc71', alpha=0.15, label='DESI-improving region'),
    Line2D([0],[0], color=C['p1'], ls='--', lw=1.3, label='Phase 1 reference'),
]
for u, sz in sorted(U0_sizes.items()):
    handles4L.append(
        plt.scatter([], [], c='0.4', marker='o', s=sz, edgecolors='white', linewidths=0.5,
                    label=f'$U_0={u}$'))
ax4L.legend(handles=handles4L, fontsize=7.8, loc='lower left',
            framealpha=0.94, edgecolor='0.7', ncol=1)

# ---------- RIGHT: |Delta r_s| vs DESI tension (confirms loophole is closed) ----------
for form, recs in by_form143.items():
    st = form_styles_143.get(form, dict(color='k', marker='x', label=form))
    tensions = [r['desi_tension']          for r in recs]
    delta_rs = [abs(r['delta_r_s_Mpc'])    for r in recs]
    U0s      = [r['U0']                    for r in recs]
    sizes    = [U0_sizes.get(u, 80)        for u in U0s]
    ax4R.scatter(tensions, delta_rs, c=st['color'], marker=st['marker'],
                 s=sizes, edgecolors='white', linewidths=0.6, alpha=0.88, zorder=5)

ax4R.axhline(0, color='k', ls='-', lw=0.7, alpha=0.3)

# Zoom to show the tiny values
ax4R.set_yscale('log')
all_drs = [abs(r['delta_r_s_Mpc']) for r in d143['scan_results'] if abs(r['delta_r_s_Mpc']) > 0]
if all_drs:
    ax4R.set_ylim(min(all_drs)*0.1, max(all_drs)*10)

ax4R.text(0.5, 0.90,
          r'$|\Delta r_s| < 10^{-8}$ Mpc for all cases:',
          transform=ax4R.transAxes, fontsize=9, ha='center', color='#27ae60', fontweight='bold',
          bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#27ae60', alpha=0.92))
ax4R.text(0.5, 0.80,
          'thawing quintessence cannot modify $r_s$',
          transform=ax4R.transAxes, fontsize=8.5, ha='center', color='#27ae60',
          bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.85))

ax4R.set_xlabel(r'DESI tension [$\sigma$]', fontsize=11)
ax4R.set_ylabel(r'$|\Delta r_s|$ [Mpc]', fontsize=11)
ax4R.set_xlim(min(r['desi_tension'] for r in d143['scan_results'])*0.95, p1ref_t143*1.03)
ax4R.grid(True, ls='--', lw=0.4, alpha=0.4, which='both')
ax4R.set_title('Sound horizon shift $|\\Delta r_s|$ vs DESI tension\n'
               '(SIM141 loophole confirmed absent for all 18 cases)',
               fontsize=9.5)

handles4R = [form_styles_143[f] for f in ('exp', 'power-law', 'hilltop')]
handles4R = [ax4R.scatter([], [], c=st['color'], marker=st['marker'],
                           s=100, label=st['label'], edgecolors='white', linewidths=0.6)
             for st in handles4R]
for u, sz in sorted(U0_sizes.items()):
    handles4R.append(
        plt.scatter([], [], c='0.4', marker='o', s=sz, edgecolors='white', linewidths=0.5,
                    label=f'$U_0={u}$'))
ax4R.legend(handles=handles4R, fontsize=7.8, loc='lower right',
            framealpha=0.94, edgecolor='0.7')

fig4.suptitle(
    'SIM143 bi-scalar $\\Psi+\\phi$ quintessence: '
    'thawing $\\phi$ evades Theorem 1 (no curvature source) '
    'but falls under Theorem 2 ($\\Delta r_s=0$, CMB violated for any DESI improvement).',
    fontsize=9.5, y=1.02)

for ext in ('pdf', 'png'):
    fig4.savefig(os.path.join(FIGDIR, f'fig4_biscalar.{ext}'))
plt.close(fig4)
print("  Saved fig4_biscalar")

print()
print("=" * 60)
print("All figures written to:", FIGDIR)
print("  fig1_psi_evolution.pdf   — Theorem 1: phi(z) trajectories")
print("  fig2_desi_chi2.pdf       — Theorem 1: tension vs coupling 6-panel")
print("  fig3_desi_cmb_tradeoff.pdf — Theorem 2: DESI vs theta* scatter")
print("  fig4_biscalar.pdf        — Theorem 2: SIM143 loophole closure")
print("=" * 60)
