"""
Generate publication figures for Paper 3: UV Finiteness and the Lz Fixed Point of CMSTG.
Produces fig1_sigma_km.pdf, fig2_diagrams.pdf, fig3_rg_flow.pdf, fig4_hierarchy.pdf
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import json, os

matplotlib.rcParams.update({
    'font.size': 12, 'axes.labelsize': 13, 'axes.titlesize': 13,
    'legend.fontsize': 11, 'xtick.labelsize': 11, 'ytick.labelsize': 11,
    'text.usetex': False, 'figure.dpi': 150,
    'axes.linewidth': 1.2, 'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.top': True, 'ytick.right': True,
})

FIGDIR = os.path.join(os.path.dirname(__file__), 'figures')
SIMDIR = os.path.expanduser('~/Ordered_Simulations')

# ── Fig 1: One-loop Σ(0) vs memory cutoff k_m ──────────────────────────────
print("Generating fig1_sigma_km.pdf ...")

with open(os.path.join(SIMDIR, 'SIM102/Outputs/sim102_diagnostics.json')) as f:
    d102 = json.load(f)

km_arr    = np.array(d102['k_m_scan'])
sigma_num = np.array(d102['Sigma_num_arr'])
sigma_ana = np.array(d102['Sigma_ana_arr'])
knat      = d102['k_m_naturalness_Mpc_inv']

# Bare quartic: extrapolate from analytic formula with no memory (sigma ~ km^4)
Lambda0 = d102['Lambda0']
km_plot = np.logspace(-3, 1.5, 200)
sigma_analytic = Lambda0**2 * km_plot**4 / (64 * np.pi**2)
# Bare theory line (same formula but labelled as bare)
sigma_bare = sigma_analytic  # same formula, just different label

fig, ax = plt.subplots(figsize=(7, 5.5))
ax.loglog(km_plot, sigma_analytic, color='C3', lw=1.5, ls='--',
          label=r'Bare theory: $\propto k_{\rm UV}^4$ (slope 4)')
ax.loglog(km_plot, sigma_analytic, color='C1', lw=1.5, ls=':',
          label=r'Analytic: $\Lambda_0^2 k_m^4/(64\pi^2)$')
ax.loglog(km_arr, sigma_num, 'C0o', ms=7, zorder=5,
          label=r'Numerical $\Sigma(0)$ (SIM102)')
ax.loglog(km_arr, sigma_ana, 'C1-', lw=1.5, alpha=0.7)
ax.axvline(knat, color='gray', lw=1.2, ls='--', alpha=0.8)
ax.text(knat*1.15, sigma_analytic[np.argmin(np.abs(km_plot-knat))]*0.3,
        r'$k_m^{\rm nat} = 9.15\,{\rm Mpc}^{-1}$',
        color='gray', fontsize=10, ha='left', va='top', rotation=90)
ax.set_xlabel(r'Memory cutoff $k_m\;[{\rm Mpc}^{-1}]$')
ax.set_ylabel(r'$\Sigma(0)$ (one-loop $\Psi$ self-energy)')
ax.set_title(r'SIM102 — Memory regulation of one-loop $\Psi$ self-energy')
ax.legend(fontsize=10)
ax.set_xlim(5e-4, 20)
plt.tight_layout()
plt.savefig(os.path.join(FIGDIR, 'fig1_sigma_km.pdf'), bbox_inches='tight')
plt.close()
print("  done.")

# ── Fig 2: Schematic diagrammatic summary ────────────────────────────────────
print("Generating fig2_diagrams.pdf ...")

fig, axes = plt.subplots(2, 2, figsize=(9, 7))
fig.suptitle(r'One-loop 1PI diagrams in CMSTG (SIM104)', fontsize=13, y=0.98)

diagram_info = [
    {
        'ax': axes[0, 0], 'label': 'A',
        'title': r'$\Psi$ self-energy',
        'desc': r'$\Sigma_A = \Lambda_0^2 k_m^4/(64\pi^2)$',
        'div': r'Bare: $k_{\rm UV}^4$ (quartic)',
        'mech': 'Memory damping',
        'color': 'C0', 'result': r'$\delta m^2/m_0^2 \approx 1$ at $k_m^{\rm nat}$',
    },
    {
        'ax': axes[0, 1], 'label': 'B',
        'title': r'Graviton wfn renorm',
        'desc': r'$\Pi_{hh}(0) = 0$ (Ward identity)',
        'div': r'Protected, $dI/dp^2|_0$ finite',
        'mech': 'Diffeomorphism invariance',
        'color': 'C2', 'result': r'$\Delta Z_h/Z_h \sim 10^{-6}$',
    },
    {
        'ax': axes[1, 0], 'label': 'C',
        'title': r'Vertex correction',
        'desc': r'$\delta\Lambda_0/\Lambda_0 = \Lambda_0^2 k_m^2/(16\pi^2 m_0^2)$',
        'div': r'Bare: $k_{\rm UV}^2$ (quadratic)',
        'mech': 'Memory damping',
        'color': 'C1', 'result': r'$\delta\Lambda_0/\Lambda_0 \approx 2.4\times10^{-2}$',
    },
    {
        'ax': axes[1, 1], 'label': 'D',
        'title': r'Two-loop $\Psi$ self-energy',
        'desc': r'$O(\Lambda_0^4)$',
        'div': r'$O(\Lambda_0^4)$ suppressed',
        'mech': 'Memory damping',
        'color': 'C3', 'result': r'$\Sigma_D/\Sigma_A \approx 3.4\times10^{-6}$',
    },
]

for info in diagram_info:
    ax = info['ax']
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.axis('off')
    # Draw loop circle
    theta = np.linspace(0, 2*np.pi, 100)
    cx, cy, r = 0.5, 0.62, 0.18
    ax.plot(cx + r*np.cos(theta), cy + r*np.sin(theta),
            color=info['color'], lw=2.5)
    # Incoming/outgoing lines
    if info['label'] in ('A', 'D'):
        ax.annotate('', xy=(cx-r-0.12, cy), xytext=(cx-r, cy),
                    arrowprops=dict(arrowstyle='->', color='k', lw=1.5))
        ax.annotate('', xy=(cx+r+0.12, cy), xytext=(cx+r, cy),
                    arrowprops=dict(arrowstyle='->', color='k', lw=1.5))
        ax.text(cx-r-0.15, cy, r'$\Psi$', ha='right', va='center', fontsize=11)
        ax.text(cx+r+0.14, cy, r'$\Psi$', ha='left', va='center', fontsize=11)
    else:
        ax.annotate('', xy=(cx-r-0.12, cy), xytext=(cx-r, cy),
                    arrowprops=dict(arrowstyle='->', color='C4', lw=2))
        ax.annotate('', xy=(cx+r+0.12, cy), xytext=(cx+r, cy),
                    arrowprops=dict(arrowstyle='->', color='C4', lw=2))
        ax.text(cx-r-0.15, cy, r'$h$', ha='right', va='center', fontsize=11,
                color='C4')
        ax.text(cx+r+0.14, cy, r'$h$', ha='left', va='center', fontsize=11,
                color='C4')
    # Label in loop
    ax.text(cx, cy, info['label'], ha='center', va='center', fontsize=12,
            fontweight='bold', color=info['color'])
    # Text below
    ax.text(0.5, 0.38, info['title'], ha='center', va='top', fontsize=11,
            fontweight='bold')
    ax.text(0.5, 0.28, info['desc'], ha='center', va='top', fontsize=9.5)
    ax.text(0.5, 0.18, info['div'], ha='center', va='top', fontsize=9,
            color='gray')
    ax.text(0.5, 0.09, r'$\rightarrow$ ' + info['result'],
            ha='center', va='top', fontsize=9, color=info['color'])

plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.savefig(os.path.join(FIGDIR, 'fig2_diagrams.pdf'), bbox_inches='tight')
plt.close()
print("  done.")

# ── Fig 3: RG flow 1/Λ₀² vs k_m ─────────────────────────────────────────────
print("Generating fig3_rg_flow.pdf ...")

with open(os.path.join(SIMDIR, 'SIM105/Outputs/sim105_diagnostics.json')) as f:
    d105 = json.load(f)

Lambda0_obs = d105['parameters']['Lambda0_obs']
m0_val      = d105['parameters']['m0']
km_nat_105  = d105['parameters']['k_m_nat']
coeff       = d105['running_coupling']['coefficient_Mpc2']

# Build data from numerical_vs_analytic
km_num, L0_num, L0_ana = [], [], []
for km_str, vals in d105['numerical_vs_analytic'].items():
    km_num.append(float(km_str))
    L0_num.append(vals['num'])
    L0_ana.append(vals['analytic'])
idx = np.argsort(km_num)
km_num = np.array(km_num)[idx]
L0_num = np.array(L0_num)[idx]
L0_ana_arr = np.array(L0_ana)[idx]

km_plot2 = np.logspace(-1, 4, 500)
inv_L0_sq_obs = 1.0 / Lambda0_obs**2
inv_L0_sq = inv_L0_sq_obs + (km_plot2**2 - km_nat_105**2) / (16 * np.pi**2 * m0_val**2)
L0_plot2 = 1.0 / np.sqrt(np.maximum(inv_L0_sq, 1e-6))

fig, ax = plt.subplots(figsize=(7, 5.5))
ax.loglog(km_num, 1.0/L0_num**2, 'C0o', ms=7, zorder=5, label='Numerical (SIM105)')
ax.loglog(km_plot2, 1.0/L0_plot2**2, 'C1-', lw=2,
          label=r'Analytic: $1/\Lambda_0^2 = 1/\Lambda_{0,{\rm obs}}^2 + (k_m^2 - k_{\rm nat}^2)/(16\pi^2 m_0^2)$')
ax.axvline(km_nat_105, color='gray', lw=1.2, ls='--', alpha=0.7)
ax.axhline(inv_L0_sq_obs, color='C3', lw=1, ls=':', alpha=0.8)
ax.text(km_nat_105*1.2, inv_L0_sq_obs*0.6,
        r'$k_m^{\rm nat} = 9.15\,{\rm Mpc}^{-1}$',
        color='gray', fontsize=9.5, rotation=90, va='top', ha='left')
ax.text(0.22, inv_L0_sq_obs*0.7,
        r'$\Lambda_{0,{\rm obs}} = 0.003$', color='C3', fontsize=10)
ax.set_xlabel(r'Memory scale $k_m\;[{\rm Mpc}^{-1}]$')
ax.set_ylabel(r'$1/\Lambda_0^2\;[M_{\rm Pl}^4]$')
ax.set_title(r'SIM105 — RG flow of $\Lambda_0$ ($\Lambda_0$ fixed point)')
ax.legend(fontsize=9.5)
# Annotate UV fixed point
ax.text(2e3, 1.0/Lambda0_obs**2 * 1.5,
        r'$\Lambda_0 \to 0$ (GR limit) as $k_m\to\infty$',
        fontsize=10, color='C1', ha='center')
ax.set_xlim(0.08, 2e4)
plt.tight_layout()
plt.savefig(os.path.join(FIGDIR, 'fig3_rg_flow.pdf'), bbox_inches='tight')
plt.close()
print("  done.")

# ── Fig 4: Perturbative hierarchy bar chart ──────────────────────────────────
print("Generating fig4_hierarchy.pdf ...")

with open(os.path.join(SIMDIR, 'SIM104/Outputs/sim104_diagnostics.json')) as f:
    d104 = json.load(f)

km_nat_key = '9.15'
sigma_A_over_m02 = d104['diagram_A']['results'][km_nat_key]['dm2_over_m02']
delta_L0_over_L0 = d104['diagram_C']['results'][km_nat_key]['delta_Lambda0_over_Lambda0']
sigma_D_over_A   = d104['diagram_D']['results'][km_nat_key]['ratio']

labels  = [r'One-loop $\Sigma_A/m_0^2$',
           r'Vertex $\delta\Lambda_0/\Lambda_0$',
           r'Two-loop $\Sigma_D/\Sigma_A$']
values  = [sigma_A_over_m02, delta_L0_over_L0, sigma_D_over_A]
colors  = ['C0', 'C1', 'C2']
hatches = ['', '//', 'xx']

fig, ax = plt.subplots(figsize=(6.5, 5.5))
bars = ax.bar(range(3), values, color=colors, edgecolor='k', linewidth=0.8,
              width=0.55)
for bar, h in zip(bars, hatches):
    bar.set_hatch(h)
ax.set_yscale('log')
ax.set_xticks(range(3))
ax.set_xticklabels(labels, fontsize=11)
ax.set_ylabel('Magnitude (dimensionless)', fontsize=12)
ax.set_title(r'Perturbative hierarchy at $k_m^{\rm nat} = 9.15\,{\rm Mpc}^{-1}$ (SIM104/SIM106)',
             fontsize=11)
# Annotate values
for i, (val, bar) in enumerate(zip(values, bars)):
    ax.text(bar.get_x() + bar.get_width()/2, val * 3,
            f'{val:.2e}', ha='center', va='bottom', fontsize=10.5)
# Power hierarchy annotation
ax.text(0.97, 0.92, r'Hierarchy: $1\ :\ 10^{-2}\ :\ 10^{-6}$',
        transform=ax.transAxes, ha='right', va='top', fontsize=11,
        bbox=dict(fc='lightyellow', ec='gray', alpha=0.85))
ax.set_ylim(1e-7, 20)
plt.tight_layout()
plt.savefig(os.path.join(FIGDIR, 'fig4_hierarchy.pdf'), bbox_inches='tight')
plt.close()
print("  done.")

print("\nAll figures generated in", FIGDIR)
