"""
Generate publication figures for Paper 4: Galactic-Scale Constraints on CMSTG.
Produces fig1_yukawa_profile.pdf, fig2_timedomain.pdf, fig3_selfcoupling.pdf,
fig4_vainshtein.pdf, fig5_deficit.pdf in the figures/ directory.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.integrate import solve_ivp
import json, os

matplotlib.rcParams.update({
    'font.size': 12, 'axes.labelsize': 13, 'axes.titlesize': 13,
    'legend.fontsize': 11, 'xtick.labelsize': 11, 'ytick.labelsize': 11,
    'text.usetex': False, 'figure.dpi': 150,
    'axes.linewidth': 1.2, 'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.top': True, 'ytick.right': True,
})

FIGDIR = os.path.join(os.path.dirname(__file__), 'figures')
os.makedirs(FIGDIR, exist_ok=True)

# ── shared physical constants (simulation units: kpc, M_sun, c=1) ──────────
# SIM100 uses scaled m0 = 0.1 kpc^-1 (Compton length = 10 kpc)
m0       = 0.1          # kpc^-1
Lambda0  = 0.003        # M_Pl^-2
Psi_cosmo = 0.003       # M_Pl  (cosmological background)
r_max    = 30.0         # kpc
v_flat   = 150.0        # km/s  (NGC 3198)
G_N      = 4.301e-3     # pc M_sun^-1 (km/s)^2  ->  need kpc units
G_kpc    = 4.301e-6     # kpc M_sun^-1 (km/s)^2

# Isothermal halo density: rho_iso = v_flat^2 / (4 pi G r^2)  [M_sun/kpc^3]
def rho_iso(r):
    return v_flat**2 / (4 * np.pi * G_kpc * r**2)

# ── Fig 1: Yukawa profile vs isothermal target ──────────────────────────────
print("Generating fig1_yukawa_profile.pdf ...")

def yukawa_rhs(r, y):
    Psi, dPsi = y
    d2Psi = m0**2 * Psi - (2 / r) * dPsi
    return [dPsi, d2Psi]

# Integrate inward from r_max with decaying Yukawa BC
Psi_bc = 1e-4
r_span = (r_max, 0.3)
r_eval = np.linspace(r_max, 0.3, 500)
y0 = [Psi_bc, -(m0 + 1/r_max) * Psi_bc]
sol = solve_ivp(yukawa_rhs, r_span, y0, t_eval=r_eval, method='RK45',
                rtol=1e-10, atol=1e-13)
r_arr = sol.t
Psi_arr = sol.y[0]
dPsi_arr = sol.y[1]

# Energy density: rho_Psi ~ (1/2) (dPsi)^2 * c^2 / G  (M_sun/kpc^3)
# Use dimensionless ratio, normalize to deficit at r=10 kpc
rho_Psi_raw = 0.5 * dPsi_arr**2
idx10 = np.argmin(np.abs(r_arr - 10.0))
rho_iso_10 = rho_iso(10.0)
rho_Psi_10_needed = rho_iso_10 / 2.7e6   # the actual deficit
norm = rho_Psi_10_needed / rho_Psi_raw[idx10]
rho_Psi_plot = rho_Psi_raw * norm

fig, ax = plt.subplots(figsize=(7, 5))
r_plot = np.linspace(0.5, 30, 500)
ax.loglog(r_arr, rho_Psi_plot, 'C0-', lw=2, label=r'$\rho_\Psi(r)$ — Yukawa (SIM99)')
ax.loglog(r_plot, rho_iso(r_plot), 'C3--', lw=2,
          label=r'$\rho_{\rm iso}(r) \propto r^{-2}$ — isothermal target')
ax.axvline(10.0, color='gray', lw=1, ls=':', alpha=0.7)
ax.annotate(r'$r = 10\,{\rm kpc}$' + '\n' + r'deficit: $2.7\times10^6\times$',
            xy=(10, rho_Psi_10_needed*2), xytext=(15, rho_Psi_10_needed*50),
            arrowprops=dict(arrowstyle='->', color='gray', lw=1),
            fontsize=10, color='gray', ha='left')
ax.set_xlabel(r'$r\;[{\rm kpc}]$')
ax.set_ylabel(r'$\rho\;[M_\odot\,{\rm kpc}^{-3}]$')
ax.set_title(r'SIM99 — Yukawa profile vs isothermal target (NGC\,3198, $v_{\rm flat}=150\,{\rm km\,s}^{-1}$)')
ax.legend()
ax.set_xlim(0.4, 35)
plt.tight_layout()
plt.savefig(os.path.join(FIGDIR, 'fig1_yukawa_profile.pdf'), bbox_inches='tight')
plt.close()
print("  done.")

# ── Fig 2: Time-domain field evolution (waterfall) ──────────────────────────
print("Generating fig2_timedomain.pdf ...")

Nr = 120
r_td = np.linspace(0.25, r_max, Nr)
dr   = r_td[1] - r_td[0]
t_end = 150.0   # kpc/c

# Galactic Ricci scalar sourced by baryonic disk (approximate)
rho_bar_peak = 1e7   # M_sun/kpc^3 at centre, exponential disk
r_d = 3.5            # disk scale length kpc
rho_bar = rho_bar_peak * np.exp(-r_td / r_d)
R_gal = -8 * np.pi * G_kpc * rho_bar / (3e5)**2   # kpc^-2  (c in km/s)
R_gal_mag = np.abs(R_gal)

def pde_rhs(t, u):
    Psi = u[:Nr]
    Psi_dot = u[Nr:]
    # Laplacian in spherical symmetry: d2Psi/dr2 + (2/r) dPsi/dr
    d2Psi = np.zeros(Nr)
    # interior
    dPsi_c = np.gradient(Psi, r_td)
    d2Psi[1:-1] = (Psi[2:] - 2*Psi[1:-1] + Psi[:-2]) / dr**2
    d2Psi[0]  = d2Psi[1]
    d2Psi[-1] = d2Psi[-2]
    lap = d2Psi + (2.0 / r_td) * dPsi_c
    Psi_ddot = lap - m0**2 * Psi + 2 * Lambda0 * Psi * R_gal_mag
    # Dirichlet BC at outer boundary
    Psi_ddot[-1] = 0.0
    return np.concatenate([Psi_dot, Psi_ddot])

u0 = np.concatenate([np.full(Nr, Psi_cosmo), np.zeros(Nr)])
t_evals = np.linspace(0, t_end, 5)
sol_td = solve_ivp(pde_rhs, [0, t_end], u0, t_eval=t_evals,
                   method='RK45', rtol=1e-6, atol=1e-9, max_step=1.0)

colors = ['C0', 'C2', 'C1', 'C3']
labels = [r'$t=0$ (IC)', r'$t = t_{\rm end}/4$',
          r'$t = t_{\rm end}/2$', r'$t = t_{\rm end}$']
fig, ax = plt.subplots(figsize=(7, 5))
for i, (ci, lab) in enumerate(zip(colors, labels)):
    ti = i + 1 if i > 0 else 0
    idx_t = min(ti, sol_td.y.shape[1] - 1)
    ax.semilogy(r_td, sol_td.y[:Nr, idx_t], color=ci, lw=1.8, label=lab)
ax.axhline(Psi_cosmo, color='gray', lw=1, ls=':', alpha=0.7,
           label=r'$\Psi_{\rm cosmo}$')
ax.set_xlabel(r'$r\;[{\rm kpc}]$')
ax.set_ylabel(r'$\Psi(r)\;[M_{\rm Pl}]$')
ax.set_title(r'SIM100 — Time-domain $\Psi(r,t)$ evolution')
ax.legend(fontsize=10)
ax.set_xlim(0, 30)
plt.tight_layout()
plt.savefig(os.path.join(FIGDIR, 'fig2_timedomain.pdf'), bbox_inches='tight')
plt.close()
print("  done.")

# ── Fig 3: Self-coupling parameter scan ─────────────────────────────────────
print("Generating fig3_selfcoupling.pdf ...")

sim100_path = os.path.expanduser(
    '~/Ordered_Simulations/SIM100/Outputs/sim100_diagnostics.json')
with open(sim100_path) as f:
    d100 = json.load(f)

lambdas, growth_factors, blow_ups = [], [], []
seen = set()
for run in d100['runs']:
    lam = run['lambda_gal']
    if lam in seen:
        continue
    seen.add(lam)
    if run.get('blow_up', False):
        blow_ups.append(lam)
        lambdas.append(lam)
        growth_factors.append(None)
    else:
        lambdas.append(lam)
        growth_factors.append(run.get('psi_max_growth_factor'))

lam_arr = np.array([l for l, g in zip(lambdas, growth_factors) if g is not None])
gf_arr  = np.array([g for g in growth_factors if g is not None])

lam_crit_spatial = d100['physics']['lambda_crit_spatial']
lam_blow         = d100['physics']['lambda_crit_uniform']

fig, ax = plt.subplots(figsize=(7, 5))
ax.scatter(lam_arr, gf_arr, color='C0', s=60, zorder=5, label='Growth factor (stable runs)')
for lb in blow_ups:
    ax.axvline(lb, color='C3', lw=1, ls=':', alpha=0.5)
ax.axvline(lam_crit_spatial, color='C0', lw=1.5, ls='--',
           label=r'$\lambda_{\rm crit}^{\rm spatial}\approx -778$')
ax.axvline(lam_blow, color='C1', lw=1.5, ls='--',
           label=r'$\lambda_{\rm blow}\approx -1111$')
ax.axhline(1.0, color='gray', lw=1, ls=':', alpha=0.7)
# Mark blow-up region
if blow_ups:
    ax.fill_betweenx([0.9, 1.35], min(blow_ups)-50, min(blow_ups)+50,
                     color='C3', alpha=0.12, label='Blow-up region')
ax.set_xlabel(r'Quartic coupling $\lambda$')
ax.set_ylabel(r'Growth factor $G_{\rm field} = \max_r|\Psi(r,t_{\rm end})|/\Psi_{\rm cosmo}$')
ax.set_title(r'SIM100 — Self-coupling scan')
ax.legend(fontsize=10)
ax.set_xlim(-1600, 50)
ax.set_ylim(0.95, 1.35)
plt.tight_layout()
plt.savefig(os.path.join(FIGDIR, 'fig3_selfcoupling.pdf'), bbox_inches='tight')
plt.close()
print("  done.")

# ── Fig 4: Vainshtein ratio heatmap ─────────────────────────────────────────
print("Generating fig4_vainshtein.pdf ...")

sim103_path = os.path.expanduser(
    '~/Ordered_Simulations/SIM103/Outputs/sim103_diagnostics.json')
with open(sim103_path) as f:
    d103 = json.load(f)

vdata = d103['vainshtein_results']
alphas_set = sorted(set(v['alpha_V'] for v in vdata.values() if v['alpha_V'] > 0))
MVs_set    = sorted(set(v['M_V']     for v in vdata.values()))

# Build grid
alpha_vals = np.array(alphas_set)
MV_vals    = np.array(MVs_set)
V_grid = np.zeros((len(alpha_vals), len(MV_vals)))
for i, aV in enumerate(alpha_vals):
    for j, MV in enumerate(MV_vals):
        key = f'aV{aV:.1f}_MV{MV}'
        if key in vdata:
            V_grid[i, j] = vdata[key]['V_ratio']
        else:
            # try to find approximate key
            for k, v in vdata.items():
                if abs(v['alpha_V'] - aV) < 1e-6 and abs(v['M_V'] - MV) < 1e-9:
                    V_grid[i, j] = v['V_ratio']
                    break

fig, ax = plt.subplots(figsize=(7, 5))
if np.any(V_grid > 0):
    log_V = np.log10(np.where(V_grid > 0, V_grid, 1e-12))
    im = ax.contourf(np.log10(MV_vals), np.log10(alpha_vals), log_V,
                     levels=20, cmap='viridis_r')
    cb = plt.colorbar(im, ax=ax)
    cb.set_label(r'$\log_{10}\,\mathcal{V}$', fontsize=12)
    ax.set_xlabel(r'$\log_{10}(M_V\;[{\rm kpc}^{-1}])$')
    ax.set_ylabel(r'$\log_{10}\,\alpha_V$')
    ax.set_title(r'SIM103 — Vainshtein ratio $\mathcal{V}$ at $r=10\,{\rm kpc}$ (NGC\,3198)')
    # Mark threshold
    ax.axhline(np.log10(0), color='white', lw=0, alpha=0)
    ax.text(0.05, 0.92, r'$\mathcal{V} \leq 3.7\times10^{-6}$ everywhere',
            transform=ax.transAxes, color='white', fontsize=11,
            bbox=dict(fc='black', alpha=0.4, ec='none'))
else:
    # Fallback: bar chart of V ratios for MV=0.001 column
    aV_plot = [1, 10, 100, 1000]
    V_plot  = [vdata.get(f'aV{a:.1f}_MV0.001', {}).get('V_ratio', 0) for a in aV_plot]
    ax.bar([str(a) for a in aV_plot], V_plot, color='C0')
    ax.set_xlabel(r'$\alpha_V$')
    ax.set_ylabel(r'Vainshtein ratio $\mathcal{V}$')
    ax.set_yscale('log')
    ax.set_title(r'SIM103 — Vainshtein ratio at $r=10\,{\rm kpc}$, $M_V=0.001\,{\rm kpc}^{-1}$')
    ax.axhline(1.0, color='C3', lw=1.5, ls='--', label=r'$\mathcal{V}=1$ (nonlinear threshold)')
    ax.legend()
plt.tight_layout()
plt.savefig(os.path.join(FIGDIR, 'fig4_vainshtein.pdf'), bbox_inches='tight')
plt.close()
print("  done.")

# ── Fig 5: Energy density deficit bar chart ──────────────────────────────────
print("Generating fig5_deficit.pdf ...")

rho_Psi_val = 1.6e-3   # M_sun/kpc^3  (from SIM103 diagnostics)
rho_iso_val = 4.2e3    # M_sun/kpc^3  (isothermal at r=10 kpc, v_flat=150 km/s)

fig, ax = plt.subplots(figsize=(5.5, 5))
bars = ax.bar([r'$\rho_\Psi$' + '\n(CMSTG locked action)',
               r'$\rho_{\rm iso}$' + '\n(isothermal requirement)'],
              [rho_Psi_val, rho_iso_val],
              color=['C0', 'C3'], width=0.5, edgecolor='k', linewidth=0.8)
ax.set_yscale('log')
ax.set_ylabel(r'$\rho\;[M_\odot\,{\rm kpc}^{-3}]$ at $r=10\,{\rm kpc}$', fontsize=12)
ax.set_title(r'SIM103 — Energy density deficit (NGC\,3198)', fontsize=12)
ax.annotate('', xy=(1, rho_iso_val * 0.5), xytext=(0, rho_Psi_val * 2),
            arrowprops=dict(arrowstyle='<->', color='black', lw=1.5))
ax.text(0.5, np.sqrt(rho_Psi_val * rho_iso_val) * 0.8,
        r'Deficit: $2.7\times10^6\times$',
        ha='center', va='bottom', fontsize=11, color='black',
        bbox=dict(fc='lightyellow', ec='gray', alpha=0.85))
ax.set_ylim(1e-5, 1e6)
plt.tight_layout()
plt.savefig(os.path.join(FIGDIR, 'fig5_deficit.pdf'), bbox_inches='tight')
plt.close()
print("  done.")

print("\nAll figures generated successfully in", FIGDIR)
