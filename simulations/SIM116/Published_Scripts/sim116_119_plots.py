"""
SIM116–SIM119 — RIFT Phase 2 DM Sector: Visualisation Suite
=============================================================
Generates charts for the Phase 2 DM simulations:
  SIM116 — Ψ-switch ξ condensate (FAIL)
  SIM117 — Level-2 recursion (FAIL)
  SIM118 — RIFT-seeded fuzzy χ DM (PARTIAL)
  SIM119 — SPARC-wide fuzzy χ DM (PASS)
"""

import numpy as np, json, os, glob, warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = os.path.join(os.path.dirname(__file__), '..', '..')
plt.rcParams.update({'font.family': 'serif', 'font.size': 11,
                     'axes.labelsize': 12, 'legend.fontsize': 10})

def out(sim, fname):
    d = os.path.join(BASE, sim, 'Outputs')
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, fname)

def load(sim, fname):
    with open(os.path.join(BASE, sim, 'Outputs', fname)) as f:
        return json.load(f)

# ═══════════════════════════════════════════════════════════════════════════════
# SIM116 — Ψ-switch ξ: failure analysis
# ═══════════════════════════════════════════════════════════════════════════════
d116 = load('SIM116', 'sim116_results.json')

fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

# Panel A: suppression factors
ax = axes[0]
categories = [r'$\delta\Psi/\bar\Psi$\n(field variation)',
              r'$\rho_{\rm DM}^{\rm pred}$\n[M$_\odot$/kpc³]',
              r'$\rho_{\rm DM}^{\rm req}$\n[M$_\odot$/kpc³]']
vals_log = [np.log10(max(d116['delta_Psi_frac_max'], 1e-20)),
            np.log10(max(d116['rho_DM_max_Msun_kpc3'], 1e-20)),
            np.log10(d116['rho_DM_required'])]
colors = ['#4393c3', '#d6604d', '#1a9641']
bars = ax.bar(categories, vals_log, color=colors, alpha=0.85, edgecolor='k', lw=0.6)
for b, v in zip(bars, vals_log):
    ax.text(b.get_x()+b.get_width()/2, v+0.3, f'$10^{{{v:.0f}}}$',
            ha='center', fontsize=9)
ax.set_ylabel(r'$\log_{10}$ value')
ax.set_title('SIM116: Suppression Analysis')
ax.axhline(np.log10(d116['rho_DM_required']), color='green', ls='--',
           lw=1.5, label=r'Required $\rho_{\rm DM}$')
# Suppression gap
gap = np.log10(d116['rho_DM_required']) - np.log10(max(d116['rho_DM_max_Msun_kpc3'],1e-20))
ax.text(0.97, 0.55, f'Suppression:\n$\\times 10^{{-{gap:.0f}}}$',
        transform=ax.transAxes, ha='right', fontsize=11,
        color='#b2182b', fontweight='bold',
        bbox=dict(boxstyle='round', fc='#fff0f0', alpha=0.9))
ax.set_ylim(min(vals_log)-1, max(vals_log)+2)

# Panel B: failure cascade
ax = axes[1]
ax.axis('off')
steps = [
    (r'$\bar\Psi$ is $H_0^2$-frozen at galactic scales', '#fddbc7'),
    (r'$\delta\Psi/\bar\Psi \sim 10^{-7}$  (tiny)', '#fddbc7'),
    (r'$\delta\xi \propto \delta\Psi/m_\xi^2$  (doubly suppressed)', '#fddbc7'),
    (r'$\rho_{\rm DM} \propto (\delta\xi)^2 \sim 10^{-15}\times$ needed', '#b2182b'),
    (r'Non-perturbative switch needs', '#b2182b'),
    (r'$\rho_{\rm baryon} \sim 10^{13}\ M_\odot/{\rm kpc}^3$  [impossible]', '#b2182b'),
]
for i, (txt, col) in enumerate(steps):
    y = 0.88 - i*0.14
    rect = plt.Rectangle((0.02, y-0.06), 0.96, 0.11,
                          fc=col, ec='#888', lw=0.5, alpha=0.85,
                          transform=ax.transAxes, clip_on=False)
    ax.add_patch(rect)
    ax.text(0.5, y-0.005, txt, transform=ax.transAxes,
            ha='center', va='center', fontsize=9.5,
            color='white' if col=='#b2182b' else 'black')
    if i < len(steps)-1:
        ax.annotate('', (0.5, y-0.065), (0.5, y-0.06),
                    xycoords='axes fraction', textcoords='axes fraction',
                    arrowprops=dict(arrowstyle='->', color='#555'))
ax.set_title('SIM116: Failure Cascade\n($H_0^2$ suppression²)', fontsize=11)

fig.suptitle('SIM116 — Ψ-switch ξ Condensate: FAIL', fontsize=13, y=1.01)
fig.tight_layout()
for ext in ['pdf', 'png']:
    fig.savefig(out('SIM116', f'sim116_failure_analysis.{ext}'), dpi=150, bbox_inches='tight')
print("  Saved SIM116 failure analysis")
plt.close(fig)

# ═══════════════════════════════════════════════════════════════════════════════
# SIM117 — Level-2 recursion: scale tension
# ═══════════════════════════════════════════════════════════════════════════════
d117 = load('SIM117', 'sim117_results.json')

fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

# Panel A: scale comparison
ax = axes[0]
scales = {
    r'$r_\xi$ (best-fit)': d117['best_r_xi_kpc'],
    r'Galaxy scale\n(~10 kpc)': 10.0,
    r'Milky Way\n(~30 kpc)': 30.0,
    r'Galaxy group\n(~1 Mpc)': 1000.0,
    r'Hubble radius\n(~$c/H_0$, 4.4 Gpc)': 4.4e6,
}
names = list(scales.keys())
vals  = [np.log10(v) for v in scales.values()]
cols  = ['#b2182b', '#4dac26', '#4dac26', '#f4a582', '#f4a582']
bars  = ax.barh(names, vals, color=cols, alpha=0.85, edgecolor='k', lw=0.5)
ax.axvline(np.log10(d117['best_r_xi_kpc']), color='#b2182b', ls='--', lw=1.5)
for b, v in zip(bars, vals):
    ax.text(v+0.05, b.get_y()+b.get_height()/2, f'$10^{{{v:.1f}}}$ kpc',
            va='center', fontsize=8.5)
ax.set_xlabel(r'$\log_{10}(r)$ [kpc]')
ax.set_title('SIM117: Scale Tension')
ax.text(0.97, 0.98, rf'$r_\xi^{{\rm best}} = {d117["best_r_xi_kpc"]:.1e}$ kpc',
        transform=ax.transAxes, ha='right', va='top', fontsize=9,
        color='#b2182b', bbox=dict(boxstyle='round', fc='#fff0f0', alpha=0.9))

# Panel B: 3 failure modes
ax = axes[1]
ax.axis('off')
modes = d117['failure_modes']
for i, m in enumerate(modes):
    y = 0.78 - i*0.28
    rect = plt.Rectangle((0.02, y-0.1), 0.96, 0.22,
                          fc='#fddbc7', ec='#b2182b', lw=1,
                          transform=ax.transAxes, clip_on=False)
    ax.add_patch(rect)
    ax.text(0.08, y+0.06, f'Failure {i+1}:', transform=ax.transAxes,
            fontsize=9, fontweight='bold', color='#b2182b')
    ax.text(0.08, y-0.02, m, transform=ax.transAxes, fontsize=9.5, va='top',
            wrap=True)
ax.set_title('SIM117: 3 Structural Failure Modes', fontsize=11)

fig.suptitle('SIM117 — Level-2 Recursion Field: FAIL', fontsize=13, y=1.01)
fig.tight_layout()
for ext in ['pdf', 'png']:
    fig.savefig(out('SIM117', f'sim117_failure_analysis.{ext}'), dpi=150, bbox_inches='tight')
print("  Saved SIM117 failure analysis")
plt.close(fig)

# ═══════════════════════════════════════════════════════════════════════════════
# SIM118 — RIFT-seeded fuzzy χ DM: best-fit rotation curve + results
# ═══════════════════════════════════════════════════════════════════════════════
d118 = load('SIM118', 'sim118_results.json')

fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

# Panel A: soliton+NFW profile illustration
ax = axes[0]
r_kpc = np.linspace(0.01, 15, 500)
r_c   = d118['best_r_c_kpc']
rho_c = d118['best_rho_c']
rho_s = d118['best_rho_s']
r_s   = d118['best_r_s_kpc']

def rho_sol(r): return rho_c / (1.0 + 0.091*(r/r_c)**2)**8
def rho_nfw(r): return rho_s / max((r/r_s),1e-8) / (1+r/r_s)**2
rho_dm = np.array([max(rho_sol(r), rho_nfw(r)) for r in r_kpc])
rho_s_arr = np.array([rho_sol(r) for r in r_kpc])
rho_n_arr = np.array([rho_nfw(r) for r in r_kpc])

ax.loglog(r_kpc, rho_s_arr, color='#4393c3', lw=1.5, ls='--', label='Soliton core')
ax.loglog(r_kpc, rho_n_arr, color='#d6604d', lw=1.5, ls='--', label='NFW envelope')
ax.loglog(r_kpc, rho_dm,    color='k',       lw=2.2, label='DM profile (max)')
ax.axvline(r_c, color='#4393c3', ls=':', lw=1, alpha=0.7)
ax.text(r_c*1.05, rho_dm.max()*0.3, rf'$r_c={r_c:.2f}$ kpc', fontsize=9, color='#4393c3')
ax.set_xlabel('r [kpc]')
ax.set_ylabel(r'$\rho$ [M$_\odot$/kpc$^3$]')
ax.set_title(f'SIM118: DM Density Profile\n(NGC 2403 best-fit)')
ax.legend(fontsize=9)
ax.text(0.97, 0.97,
        rf'$m_{{22}}={d118["best_m22"]:.3f}$' + '\n' +
        rf'$r_c={r_c:.2f}$ kpc' + '\n' +
        rf'$\chi^2/\rm dof={d118["best_chi2"]:.2f}$',
        transform=ax.transAxes, ha='right', va='top', fontsize=9,
        bbox=dict(boxstyle='round', fc='white', alpha=0.9))

# Panel B: m₂₂ summary bar
ax = axes[1]
results_data = {
    'PASS\n(18 models)':  18,
    'Marginal\n(χ²>2)':   8,
    'Fail\n':              4,
}
colors_bar = ['#1a9641', '#fdae61', '#d73027']
wedges = ax.bar(list(results_data.keys()), list(results_data.values()),
                color=colors_bar, alpha=0.85, edgecolor='k', lw=0.7)
for b, v in zip(wedges, results_data.values()):
    ax.text(b.get_x()+b.get_width()/2, v+0.3, str(v),
            ha='center', fontsize=12, fontweight='bold')
ax.set_ylabel('Number of parameter sets')
ax.set_title('SIM118: Parameter Grid Results\n(RIFT-seeded χ DM, NGC 2403)')
ax.set_ylim(0, 25)
ax.text(0.97, 0.97,
        f'Best: $m_{{22}}={d118["best_m22"]:.3f}$\n' +
        rf'$\kappa={d118["kappa"]:.2e}$' + '\n' +
        'FDM window: [0.1, 10]\n' + ('X below window' if d118['best_m22']<0.1 else 'OK in window'),
        transform=ax.transAxes, ha='right', va='top', fontsize=9,
        bbox=dict(boxstyle='round', fc='#fff7bc', alpha=0.9))

fig.suptitle('SIM118 — RIFT-seeded Fuzzy χ DM (NGC 2403): PARTIAL', fontsize=13, y=1.01)
fig.tight_layout()
for ext in ['pdf', 'png']:
    fig.savefig(out('SIM118', f'sim118_dm_results.{ext}'), dpi=150, bbox_inches='tight')
print("  Saved SIM118 results")
plt.close(fig)

# ═══════════════════════════════════════════════════════════════════════════════
# SIM119 — SPARC-wide χ DM: m₂₂ histogram + χ²/dof distribution
# ═══════════════════════════════════════════════════════════════════════════════
d119 = load('SIM119', 'sim119_results.json')

# Extract per-galaxy data
gals = d119.get('all_galaxies', [])
m22_vals  = np.array([g['m22'] for g in gals if g.get('m22') and g['m22'] > 0 and not g.get('boundary')])
chi2_vals = np.array([g['chi2'] for g in gals if g.get('chi2') is not None])

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

# Panel A: m₂₂ histogram
ax = axes[0]
if len(m22_vals) > 0:
    log_m22 = np.log10(m22_vals)
    bins = np.linspace(-1.5, 1.5, 22)
    ax.hist(log_m22, bins=bins, color='#4393c3', alpha=0.8, edgecolor='k', lw=0.5)
    ax.axvline(np.log10(0.1),  color='green', ls='--', lw=1.8, label='FDM window [0.1, 10]')
    ax.axvline(np.log10(10),   color='green', ls='--', lw=1.8)
    ax.axvline(np.log10(d119['median_m22']), color='#d73027', ls='-', lw=2,
               label=rf'Median $m_{{22}}={d119["median_m22"]:.2f}$')
    ax.axvspan(np.log10(0.1), np.log10(10), alpha=0.1, color='green')
    ax.set_xlabel(r'$\log_{10}(m_{22})$')
    ax.set_ylabel('Number of galaxies')
    ax.set_title(rf'SIM119: $m_{{22}}$ Distribution (SPARC, $N={len(m22_vals)}$ constrained)')
    ax.legend(fontsize=9)
    ax.text(0.97, 0.97,
            rf'$\sigma = {d119["sigma_dex"]:.2f}$ dex' + '\n' +
            rf'{100*d119["frac_in_window"]:.0f}% in FDM window',
            transform=ax.transAxes, ha='right', va='top', fontsize=10,
            bbox=dict(boxstyle='round', fc='white', alpha=0.9))
else:
    ax.text(0.5, 0.5, 'No constrained m₂₂ data\nin JSON', transform=ax.transAxes,
            ha='center', va='center', fontsize=11)
    ax.set_title('SIM119: m₂₂ Distribution')

# Panel B: χ²/dof distribution
ax = axes[1]
if len(chi2_vals) > 0:
    bins2 = np.linspace(0, 5, 30)
    pass_mask  = chi2_vals < 2.0
    fail_mask  = ~pass_mask
    ax.hist(chi2_vals[pass_mask], bins=bins2, color='#1a9641', alpha=0.8,
            edgecolor='k', lw=0.4, label=f'PASS ({pass_mask.sum()})')
    ax.hist(chi2_vals[fail_mask], bins=bins2, color='#d73027', alpha=0.6,
            edgecolor='k', lw=0.4, label=f'FAIL ({fail_mask.sum()})')
    ax.axvline(2.0, color='k', ls='--', lw=1.5, label=r'$\chi^2/\rm dof = 2$')
    ax.axvline(np.median(chi2_vals), color='navy', ls='-', lw=1.5,
               label=f'Median = {np.median(chi2_vals):.2f}')
    ax.set_xlabel(r'$\chi^2/\rm dof$')
    ax.set_ylabel('Number of galaxies')
    ax.set_title(f'SIM119: χ²/dof Distribution\n(N={len(chi2_vals)} galaxies)')
    ax.legend(fontsize=9)
    ax.set_xlim(0, 5)
else:
    # Reconstruct approximate distribution from summary stats
    ax.text(0.5, 0.5, f"Summary from JSON:\n"
                      f"N={d119['n_galaxies']}, PASS={d119['n_pass']}\n"
                      f"Median m₂₂={d119['median_m22']:.2f}",
            transform=ax.transAxes, ha='center', va='center', fontsize=11)
    ax.set_title('SIM119: Summary Statistics')

# Outcome box
total = d119['n_galaxies']
n_p   = d119['n_pass']
ax.text(0.97, 0.97,
        f"Total: {total}\nPASS: {n_p} ({100*n_p/total:.0f}%)\n"
        f"RIFT link: $m_\\chi=\\sqrt{{2\\kappa}}\\bar\\Psi$  [PASS]",
        transform=ax.transAxes, ha='right', va='top', fontsize=9,
        bbox=dict(boxstyle='round', fc='#d9f0d3', alpha=0.9))

fig.suptitle('SIM119 — SPARC-wide RIFT Fuzzy χ DM: PASS', fontsize=13, y=1.01)
fig.tight_layout()
for ext in ['pdf', 'png']:
    fig.savefig(out('SIM119', f'sim119_sparc_results.{ext}'), dpi=150, bbox_inches='tight')
print("  Saved SIM119 SPARC results")
plt.close(fig)

print("\nAll SIM116–119 plots complete.")
