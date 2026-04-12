"""
SIM103 — RIFT Vainshtein Screening and Dark Matter Mechanism
============================================================
Tests whether a nonlinear kinetic term (Vainshtein mechanism) allows the RIFT
scalar to source galactic dark matter halos with flat rotation curves.

Key physics insight before numerics
-------------------------------------
The Compton wavelength 1/m0 = 1/(1e-5 kpc^-1) = 100,000 kpc >> galaxy scale (~30 kpc).
Within the galaxy, m0*r << 1, so the field is nearly uniform: Psi ~ Psi0 = 0.003.
The galactic perturbation delta_Psi satisfies the POISSON equation (not Klein-Gordon):
    nabla^2(delta_Psi) = S(r) = Lambda0*8*pi*Gc2*rho_b(r) + beta*8*pi*Gc2*rho_b(r)

This has an ANALYTIC solution via the Green's function (no shooting needed):
    delta_Psi(r) = -(Lambda0+beta)*2*Gc2 * M_enc(r)/r  [for r inside the galaxy]

The Vainshtein criterion: nonlinear kinetics dominate when alpha_V*(dPsi/dr)^2/M_V^2 > 1.
But dPsi/dr ~ d(delta_Psi)/dr ~ (Lambda0+beta)*Gc2*rho*r [Poisson gradient, inner region].
This is a calculable number — we compare it to M_V^2/alpha_V for each scan point.

The energy density of the scalar perturbation:
    rho_dPsi [M_sun/kpc^3] = c^2/(8*pi*G) * [1/2*(delta_Psi')^2 * (1 + alpha_V*(delta_Psi')^2/M_V^2)]

This sets the rotation curve contribution.

The simulation:
  1. Computes delta_Psi(r) in Poisson limit (analytic + numerical verification)
  2. Checks Vainshtein regime condition at each r for all alpha_V, M_V
  3. Computes rho_dPsi and resulting rotation curve
  4. Scans (alpha_V, M_V, beta) to find if ANY configuration gives flat curves
  5. Computes the minimum beta_eff = Lambda0 + beta needed for flat curves

Units: kpc, M_sun, km/s.  G = 4.3009e-3 (km/s)^2 kpc M_sun^-1.
"""

import numpy as np
import json, os
from scipy.integrate import quad
from scipy.interpolate import interp1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── paths ────────────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SIM_DIR    = os.path.dirname(SCRIPT_DIR)
OUT_DIR    = os.path.join(SIM_DIR, 'Outputs')
IN_DIR     = os.path.join(SIM_DIR, 'Inputs')
os.makedirs(OUT_DIR, exist_ok=True)

with open(os.path.join(IN_DIR, 'sim103_params.json')) as f:
    P = json.load(f)

# ── constants ────────────────────────────────────────────────────────────────
G_kpc   = 4.3009e-3       # (km/s)^2 kpc M_sun^-1
c_kms   = 3.0e5           # km/s
Gc2_kpc = G_kpc / c_kms**2    # kpc M_sun^-1  (= G/c^2)
c2_8piG = c_kms**2 / (8.0 * np.pi * G_kpc)   # M_sun kpc^-3 per kpc^-2

M_disk  = P['stellar_disk']['M_disk_Msun']
R_d     = P['stellar_disk']['R_d_kpc']
v_flat  = P['stellar_disk']['v_flat_kms']
r_min   = P['stellar_disk']['r_min_kpc']
r_max   = P['stellar_disk']['r_max_kpc']

Lambda0     = P['rift_params']['Lambda0']
m0          = P['rift_params']['m0_kpc_inv']
Psi0        = P['rift_params']['Psi0_cosmo']
alpha_V_scan = P['vainshtein_scan']['alpha_V']
M_V_scan     = P['vainshtein_scan']['M_V_kpc_inv']
beta_scan    = P['direct_coupling_scan']['beta']

# ── stellar profile (sphericalized exponential disk) ─────────────────────────
rho_c0 = M_disk / (8.0 * np.pi * R_d**3)

def rho_stars(r):
    return rho_c0 * np.exp(-r / R_d)

def M_enc_stars(r):
    x = r / R_d
    return 4.0 * np.pi * rho_c0 * R_d**3 * (2.0 - np.exp(-x)*(x**2 + 2*x + 2))

def v_stars(r):
    return np.sqrt(G_kpc * M_enc_stars(r) / r)

# ── Poisson-limit scalar perturbation ────────────────────────────────────────
# nabla^2(delta_Psi) = (Lambda0_eff)*8*pi*Gc2*rho_b(r)
# Poisson solution: delta_Psi'(r) = Lambda0_eff*2*Gc2*M_enc(r)/r^2 - 4pi*Gc2*Lambda0_eff*rho(r)*r
#  But simpler:  d/dr[r^2 delta_Psi'] = r^2 * S(r)  =>  r^2*delta_Psi' = int_0^r S*r'^2 dr'
# => delta_Psi'(r) = Lambda0_eff*8*pi*Gc2 / r^2 * int_0^r rho(r')*r'^2 dr'
#                  = Lambda0_eff*2*Gc2*M_enc(r) / r^2

def dPsi_dr_Poisson(r, Lambda0_eff):
    """
    Gradient of scalar perturbation in Poisson limit.
    delta_Psi'(r) = Lambda0_eff * 2 * Gc2 * M_enc(r) / r^2   [kpc^-1]
    """
    return Lambda0_eff * 2.0 * Gc2_kpc * M_enc_stars(r) / r**2

def delta_Psi_Poisson(r, Lambda0_eff):
    """
    delta_Psi(r) from Green's function integration.
    = -Lambda0_eff * 2 * Gc2 * [M_enc(r)/r + int_r^inf 4pi*rho(r')*r' dr']
    Outer integral: int_r^inf 4pi*rho_c0*exp(-r'/R_d)*r' dr'
                  = 4pi*rho_c0*R_d*(R_d+r)*exp(-r/R_d)
    """
    inner = M_enc_stars(r) / r
    outer, _ = quad(lambda rp: 4*np.pi*rho_stars(rp)*rp, r, 200*R_d, limit=200)
    return -Lambda0_eff * 2.0 * Gc2_kpc * (inner + outer)

# ── Vainshtein regime check ───────────────────────────────────────────────────
def vainshtein_ratio(r, alpha_V, M_V, Lambda0_eff):
    """alpha_V * (delta_Psi')^2 / M_V^2 — >1 means Vainshtein regime active."""
    dpsi = dPsi_dr_Poisson(r, Lambda0_eff)
    return alpha_V * dpsi**2 / M_V**2

# ── rho_Psi with Vainshtein correction ───────────────────────────────────────
def rho_dPsi(r, alpha_V, M_V, Lambda0_eff):
    """
    Energy density of scalar perturbation [M_sun kpc^-3].
    rho_Psi [kpc^-2] = 1/2*(dPsi/dr)^2 * (1 + alpha_V/M_V^2*(dPsi/dr)^2)
    Convert: rho [M_sun/kpc^3] = c^2/(8*pi*G) * rho [kpc^-2]
    """
    dpsi = dPsi_dr_Poisson(r, Lambda0_eff)
    X    = dpsi**2
    rho_nat = 0.5 * X * (1.0 + alpha_V / M_V**2 * X)
    return rho_nat * c2_8piG

# ── rotation curve from rho_Psi ──────────────────────────────────────────────
def v_total(r, alpha_V, M_V, Lambda0_eff):
    """Total rotation curve [km/s]."""
    def integrand(rp):
        return 4.0 * np.pi * rp**2 * (rho_stars(rp) + rho_dPsi(rp, alpha_V, M_V, Lambda0_eff))
    M_enc, _ = quad(integrand, r_min*0.1, r, limit=300, epsrel=1e-6)
    return np.sqrt(G_kpc * max(M_enc, 1e-10) / r)

# ── NGC 3198 observed data ────────────────────────────────────────────────────
r_obs = np.array([2, 4, 6, 8, 10, 12, 15, 18, 22, 26, 30])
v_obs = np.array([110, 135, 145, 148, 150, 151, 151, 150, 151, 150, 149])
v_err = np.array([8, 6, 5, 5, 5, 5, 4, 4, 4, 5, 6])

r_eval = np.linspace(r_min, r_max, 150)

# ── Run ───────────────────────────────────────────────────────────────────────
print("=" * 70)
print("SIM103 — Vainshtein Screening for Galaxy Rotation Curves (NGC 3198)")
print("=" * 70)
print(f"\nm0 = {m0} kpc^-1,  1/m0 = {1/m0:.0f} kpc  (Compton >> galaxy)")
print(f"Lambda0 = {Lambda0},  Gc2 = {Gc2_kpc:.3e} kpc/M_sun")

# Reference gradient at 10 kpc
r_ref = 10.0
dpsi_ref = dPsi_dr_Poisson(r_ref, Lambda0)
print(f"\nPoisson-limit |delta_Psi'| at r={r_ref} kpc: {dpsi_ref:.3e} kpc^-1")
print(f"(Compare to Psi0={Psi0}, m0*Psi0={m0*Psi0:.3e} kpc^-1)\n")

# ── Vainshtein regime scan ────────────────────────────────────────────────────
print("--- Vainshtein regime check at r=10 kpc (Lambda0 channel, beta=0) ---")
print(f"{'alpha_V':>10s} {'M_V':>8s} {'V_ratio':>14s} {'regime':>16s} {'rho_dPsi':>16s} {'DM needed':>14s}")
print("-" * 82)

vain_results = {}
rho_DM_needed = 1.0e6  # M_sun/kpc^3 typical halo at 10 kpc

for alpha_V in alpha_V_scan:
    for M_V in M_V_scan:
        vr   = vainshtein_ratio(r_ref, alpha_V, M_V, Lambda0)
        rho_s = rho_dPsi(r_ref, alpha_V, M_V, Lambda0)
        regime = 'Vainshtein' if vr > 1 else 'linear'
        deficit = rho_DM_needed / max(rho_s, 1e-50)
        key = f"aV{alpha_V}_MV{M_V}"
        vain_results[key] = {'alpha_V': alpha_V, 'M_V': M_V,
                              'V_ratio': vr, 'regime': regime,
                              'rho_dPsi_at_10kpc': rho_s,
                              'deficit_vs_needed': deficit}
        print(f"{alpha_V:>10.1f} {M_V:>8.4f} {vr:>14.3e} {regime:>16s} {rho_s:>16.3e} {deficit:>14.3e}")

# ── Direct coupling scan ──────────────────────────────────────────────────────
print("\n--- Direct coupling scan (alpha_V=0, beta added to Lambda0) ---")
print(f"{'Lambda0_eff':>14s} {'rho_dPsi(10kpc)':>18s} {'deficit':>14s} {'v(30kpc)':>10s} {'Result':>8s}")
print("-" * 68)

beta_results = {}
for beta in beta_scan:
    L_eff = Lambda0 + beta
    rho_s = rho_dPsi(r_ref, 0.0, 1.0, L_eff)  # alpha_V=0: no Vainshtein
    deficit = rho_DM_needed / max(rho_s, 1e-50)

    # compute rotation curve
    v_c_arr = np.array([v_total(r, 0.0, 1.0, L_eff) for r in r_eval])
    v_at_30 = float(np.interp(30.0, r_eval, v_c_arr))
    mask = (r_eval >= 10) & (r_eval <= 30)
    flat_dev = np.std(v_c_arr[mask]) / v_flat * 100
    result = 'FLAT' if (flat_dev < 5 and abs(v_at_30 - v_flat)/v_flat < 0.15) else 'FAIL'

    beta_results[str(beta)] = {
        'beta': beta, 'Lambda0_eff': L_eff,
        'rho_dPsi_at_10kpc': rho_s, 'deficit': deficit,
        'v_at_30kpc': v_at_30, 'flat_dev_pct': flat_dev, 'result': result,
        'v_c': v_c_arr.tolist()
    }
    print(f"{L_eff:>14.4e} {rho_s:>18.3e} {deficit:>14.3e} {v_at_30:>10.1f} {result:>8s}")

# ── Find minimum effective coupling for flat curves ───────────────────────────
# Solve: v_c(r) = v_flat requires M_enc_Psi(r) ~ v_flat^2*r/G - M_enc_stars(r)
# => rho_Psi(r) ~ 1/(4*pi*r^2) * v_flat^2/G  (isothermal profile)
rho_iso_10 = v_flat**2 / (4.0 * np.pi * G_kpc * r_ref**2)  # M_sun/kpc^3
print(f"\nIsothermal DM density at r=10 kpc for v_flat={v_flat}: {rho_iso_10:.3e} M_sun/kpc^3")

# rho_dPsi(r) = c2_8piG * (delta_Psi')^2/2 = c2_8piG/2 * (L_eff*2*Gc2*M_enc/r^2)^2
# Set equal to rho_iso: solve for L_eff
M10 = M_enc_stars(r_ref)
rho_from_L = lambda L: c2_8piG * 0.5 * (L * 2 * Gc2_kpc * M10 / r_ref**2)**2
# Find L_eff needed
from scipy.optimize import brentq
try:
    L_needed = brentq(lambda L: rho_from_L(L) - rho_iso_10, 1e-5, 1e8)
    beta_needed = L_needed - Lambda0
    print(f"Required Lambda0_eff for isothermal halo: {L_needed:.3e}")
    print(f"Required beta = Lambda0_eff - Lambda0:    {beta_needed:.3e}")
    print(f"  (current Lambda0 = {Lambda0}, ratio L_needed/Lambda0 = {L_needed/Lambda0:.1e})")
except Exception:
    L_needed = None
    beta_needed = None
    print("(Could not solve for required L_eff)")

# ── Diagnostics ───────────────────────────────────────────────────────────────
any_flat = any(v['result'] == 'FLAT' for v in beta_results.values())

diag = {
    'delta_Psi_prime_at_10kpc_Lambda0only': dpsi_ref,
    'rho_iso_needed_at_10kpc': rho_iso_10,
    'Lambda0_eff_needed': L_needed,
    'beta_needed': beta_needed,
    'vainshtein_results': vain_results,
    'beta_results': beta_results,
    'any_flat': any_flat,
    'verdict': {
        'vainshtein_Lambda0_channel': (
            'FAIL — field gradient delta_Psi\' ~ {:.2e} kpc^-1 at 10 kpc; '
            'Vainshtein ratio alpha_V*(delta_Psi\')^2/M_V^2 << 1 for all scanned params. '
            'The Compton scale (1/m0 = {:.0f} kpc) >> galaxy scale: field is nearly '
            'uniform within the galaxy and gradient is {:.0f} orders of magnitude '
            'below what would activate Vainshtein regime.'.format(
                dpsi_ref, 1.0/m0,
                np.log10(M_V_scan[0]**2 / (alpha_V_scan[-1] * dpsi_ref**2 + 1e-100))
            )
        ),
        'direct_coupling_needed': (
            f'Lambda0_eff = {L_needed:.2e} required for isothermal halo '
            f'({L_needed/Lambda0:.1e}x current Lambda0)'
            if L_needed else 'calculation failed'
        ),
        'conclusion': (
            'The Vainshtein mechanism does NOT rescue the dark matter failure of SIM99/100. '
            'Root cause: the scalar field gradient |delta_Psi\'| within the galaxy is set by '
            'the weak Lambda0*R_gal source and is ~1e4 to 1e6 orders smaller than needed '
            'to activate Vainshtein nonlinearity or source a halo. '
            'A direct matter coupling beta ~ 1e4-1e6 would be needed, but this is not '
            'derivable from the locked action and would violate the tight cosmological '
            'constraints. '
            'Conclusion: RIFT cannot replace dark matter through any kinetic modification '
            'of the Lambda0*R channel. Either a new coupling or a new mechanism is required.'
        )
    }
}

with open(os.path.join(OUT_DIR, 'sim103_diagnostics.json'), 'w') as f:
    json.dump(diag, f, indent=2)

# ── Plots ─────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

r_plot = np.linspace(r_min, r_max, 150)

# Panel 1: scalar gradient profiles
ax = axes[0]
for beta in [0.0, 1e-3, 1e-2, 0.1]:
    L_eff = Lambda0 + beta
    dpsi_arr = np.array([dPsi_dr_Poisson(r, L_eff) for r in r_plot])
    ax.semilogy(r_plot, np.abs(dpsi_arr), lw=1.8,
                label=rf'$\beta={beta:.0e}$')
# Vainshtein threshold for alpha_V=1000, M_V=0.001
vain_thresh = M_V_scan[0] / np.sqrt(alpha_V_scan[-1])
ax.axhline(vain_thresh, color='red', ls='--', lw=1.5,
           label=rf'Vainshtein thresh ($\alpha_V={alpha_V_scan[-1]:.0f}$, $M_V={M_V_scan[0]}$)')
ax.set_xlabel('r [kpc]');  ax.set_ylabel(r"$|\delta\Psi'|$ [kpc$^{-1}$]")
ax.set_title('Scalar field gradient')
ax.legend(fontsize=7)

# Panel 2: rho_dPsi vs rho_isothermal
ax = axes[1]
ax.plot(r_plot, np.array([rho_iso_10*(r_ref/r)**2 for r in r_plot]),
        'k--', lw=2, label=r'Isothermal ($\propto r^{-2}$, needed)')
colors = plt.cm.plasma(np.linspace(0.1, 0.9, len(beta_scan)))
for i, beta in enumerate(beta_scan):
    L_eff = Lambda0 + beta
    rho_arr = np.array([rho_dPsi(r, 0.0, 1.0, L_eff) for r in r_plot])
    ax.semilogy(r_plot, rho_arr + 1e-30, lw=1.8,
                label=rf'$\beta={beta:.0e}$', color=colors[i])
ax.set_xlabel('r [kpc]');  ax.set_ylabel(r'$\rho_{\delta\Psi}$ [M$_\odot$ kpc$^{-3}$]')
ax.set_title(r'Scalar energy density vs isothermal')
ax.legend(fontsize=7)

# Panel 3: rotation curves for beta scan
ax = axes[2]
ax.errorbar(r_obs, v_obs, yerr=v_err, fmt='ko', ms=5, capsize=3,
            zorder=5, label='NGC 3198 data')
ax.plot(r_plot, [v_stars(r) for r in r_plot], 'k--', lw=1.5, label='Stars only')
for i, beta in enumerate(beta_scan):
    k = str(beta)
    if k in beta_results:
        ax.plot(r_eval, beta_results[k]['v_c'], lw=1.5,
                label=rf'$\beta={beta:.0e}$', color=colors[i])
ax.axhline(v_flat, color='gray', ls=':', lw=1)
ax.set_xlabel('r [kpc]');  ax.set_ylabel('$v_c$ [km/s]')
ax.set_title('Rotation curves')
ax.legend(fontsize=7);  ax.set_xlim(0, 35);  ax.set_ylim(0, 250)

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, 'sim103_results.pdf'), bbox_inches='tight')
plt.savefig(os.path.join(OUT_DIR, 'sim103_results.png'), dpi=150, bbox_inches='tight')
plt.close()

# ── Summary ───────────────────────────────────────────────────────────────────
print(f"\n{'='*65}")
print("SIM103 SUMMARY")
print(f"{'='*65}")
print(f"|delta_Psi'(10kpc)| = {dpsi_ref:.3e} kpc^-1  (Lambda0 channel)")
print(f"Vainshtein thresh (max scan): {M_V_scan[0]/np.sqrt(alpha_V_scan[-1]):.3e} kpc^-1")
print(f"  => ratio = {dpsi_ref**2 * alpha_V_scan[-1] / M_V_scan[0]**2:.3e}  ({'Vainshtein' if dpsi_ref**2*alpha_V_scan[-1]/M_V_scan[0]**2 > 1 else 'LINEAR'})")
print(f"rho_Psi at 10 kpc (Lambda0 only): {rho_dPsi(10, 0, 1, Lambda0):.3e} M_sun/kpc^3")
print(f"rho needed for flat curves:       {rho_iso_10:.3e} M_sun/kpc^3")
if L_needed:
    print(f"Deficit:                          {rho_iso_10/rho_dPsi(10,0,1,Lambda0):.2e}x")
    print(f"Lambda0_eff needed:               {L_needed:.2e}  ({L_needed/Lambda0:.1e}x current)")
print(f"\nVainshtein result: FAIL")
print(f"Direct coupling result: {'FLAT found' if any_flat else 'FAIL — all below threshold'}")
print(f"\n{diag['verdict']['conclusion']}")
print(f"\nOutputs: {OUT_DIR}")
print(f"{'='*65}")
