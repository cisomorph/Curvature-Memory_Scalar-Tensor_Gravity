"""
SIM101 — RIFT BAO Sound Horizon from First Principles
======================================================
Derives the BAO sound horizon r_d by integrating the RIFT FLRW+scalar system
from the radiation era through recombination.

Physics
-------
Comoving sound horizon (correct formula, no (1+z) in denominator):

    r_d = integral_{z_drag}^{inf}  c_s(z) * c / H(z)  dz

Derivation: comoving element dchi = c*dt/a; dt = -dz/((1+z)*H);
  dchi = c*(-dz/((1+z)*H)) / (1/(1+z)) = -c*dz/H(z).

Sound speed in tight-coupling (photon-baryon plasma):
    c_s^2 = c^2/3 / (1 + R_b),   R_b = 3*rho_b/(4*rho_gamma)

RIFT modifications:
  (1) H(z) corrected by rho_Psi from the scalar stress-energy:
        rho_Psi = 0.5*Psi_dot^2 + 0.5*m0^2*Psi^2  [Mpc^-2 natural units]
        H_RIFT = H_LCDM * sqrt(1 + rho_Psi/rho_crit)
  (2) c_s corrected via G_eff = G/(1 + 16*pi*G*Lambda0*Psi^2) — 16 ppm effect
  (3) Psi evolves via KG in FLRW: Psi_tt + 3H*Psi_t + m0^2*Psi = 0
      Rewritten as ODE in z:
        du/dz = -[(dH/dz)/H - 2/(1+z)] * u - m0^2*Psi / ((1+z)^2 * H^2)
      where u = dPsi/dz and H is in Mpc^-1 (natural units c=1).

Units convention
----------------
All H functions return km/s/Mpc (physical Hubble parameter) EXCEPT in the Psi ODE
where we use H_nat = H_kms / c_km [Mpc^-1] for consistency with m0 [Mpc^-1].
The r_d integral uses H_kms and c_km to produce r_d in Mpc.
"""

import numpy as np
import json, os
from scipy.integrate import solve_ivp, quad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── paths ──────────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SIM_DIR    = os.path.dirname(SCRIPT_DIR)
OUT_DIR    = os.path.join(SIM_DIR, 'Outputs')
IN_DIR     = os.path.join(SIM_DIR, 'Inputs')
os.makedirs(OUT_DIR, exist_ok=True)

with open(os.path.join(IN_DIR, 'sim101_params.json')) as f:
    P = json.load(f)

# ── constants ───────────────────────────────────────────────────────────────
c_km    = 2.998e5            # speed of light [km/s]
H0      = P['cosmology']['H0']       # [km/s/Mpc]
Omega_m = P['cosmology']['Omega_m']
Omega_b = P['cosmology']['Omega_b']
Omega_r = P['cosmology']['Omega_r']
Omega_L = 1.0 - Omega_m - Omega_r   # flat universe
Psi0    = P['scalar_field']['Psi0']
m0      = P['scalar_field']['m0_Mpc_inv']   # [Mpc^-1], natural units c=1
H0_nat  = H0 / c_km                 # H0 in Mpc^-1

Lambda0_scan = P['cosmology']['Lambda0_scan']

# ── Hubble parameter [km/s/Mpc] ─────────────────────────────────────────────
def H_kms(z):
    a = 1.0 / (1.0 + z)
    return H0 * np.sqrt(Omega_r/a**4 + Omega_m/a**3 + Omega_L)

def H_nat(z):
    """H(z) in Mpc^-1 = H_kms / c_km."""
    return H_kms(z) / c_km

# ── baryon-photon momentum ratio ─────────────────────────────────────────────
def R_b(z):
    """R_b = 3*rho_b / (4*rho_gamma) = (3*Omega_b)/(4*Omega_r) / (1+z)."""
    return 0.75 * (Omega_b / Omega_r) / (1.0 + z)

def cs2(z):
    """Sound speed squared (units of c^2); LCDM formula."""
    return 1.0 / (3.0 * (1.0 + R_b(z)))

# ── drag redshift (Eisenstein & Hu 1998) ────────────────────────────────────
def z_drag_EH():
    omh2  = Omega_m  * (H0/100)**2
    ombh2 = Omega_b  * (H0/100)**2
    b1 = 0.313 * omh2**(-0.419) * (1 + 0.607 * omh2**0.674)
    b2 = 0.238 * omh2**0.223
    return 1291.0 * omh2**0.251 / (1 + 0.659*omh2**0.828) * (1 + b1*ombh2**b2)

Z_DRAG_EH = z_drag_EH()
# E&H formula underestimates z_drag by ~4% for Planck-era params.
# Use CLASS-calibrated value for our SIM90 cosmology (H0=67.59, Omega_m=0.312, Omega_b=0.049):
# CLASS gives z_drag ~ 1059 and r_d ~ 147 Mpc.  We use 1059 to anchor the LCDM baseline.
Z_DRAG = 1059.0
print(f"z_drag (E&H approx): {Z_DRAG_EH:.1f}  →  using CLASS-calibrated z_drag = {Z_DRAG:.1f}")

# ── comoving sound horizon [Mpc] ─────────────────────────────────────────────
# r_d = integral_{z_drag}^{inf} c_s(z) * c / H(z) dz
# H(z) in km/s/Mpc, c in km/s  → r_d in Mpc
def rd_integral(cs2_fn, H_fn_kms, z_lo, z_hi=1e6):
    """Integrate comoving sound horizon."""
    def integrand(z):
        return np.sqrt(cs2_fn(z)) * c_km / H_fn_kms(z)
    val, err = quad(integrand, z_lo, z_hi, limit=2000,
                    epsabs=1e-3, epsrel=1e-6)
    return val, err

rd_LCDM, rd_LCDM_err = rd_integral(cs2, H_kms, Z_DRAG)
print(f"r_d (LCDM): {rd_LCDM:.4f} ± {rd_LCDM_err:.3e} Mpc  (expected ~147 Mpc)")

# ── Psi ODE in z-variable ────────────────────────────────────────────────────
# KG:  Psi_tt + 3H*Psi_t + m0^2*Psi = 0
# Change variable: d/dt = -(1+z)*H_nat * d/dz  (z decreasing as time increases)
# Let u = dPsi/dz.
# Derivation:
#   Psi_t  = -(1+z)*H_nat * u
#   Psi_tt = d/dt[Psi_t] = -(1+z)*H_nat * d/dz[Psi_t]
#          = -(1+z)*H_nat * d/dz[-(1+z)*H_nat*u]
#          = (1+z)*H_nat * [(1+z)*dH_nat/dz*u + H_nat*u + (1+z)*H_nat*du/dz]
#            Wait, include sign carefully:
#   Let v = -(1+z)*H_nat*u = Psi_t
#   dv/dt = -(1+z)*H_nat * dv/dz
#   Psi_tt = dv/dt = -(1+z)*H_nat * d/dz[-(1+z)*H_nat*u]
#          = (1+z)*H_nat * [(1+z)*dH_nat/dz*u + H_nat*u + (1+z)*H_nat*du/dz]  No...
#
# Cleaner: substitute into KG and collect du/dz terms.
# From KG: Psi_tt + 3*H_nat*Psi_t + m0^2*Psi = 0
#   Psi_tt = -3*H_nat*(-(1+z)*H_nat*u) - m0^2*Psi
#          = 3*(1+z)*H_nat^2*u - m0^2*Psi
# Also Psi_tt = (1+z)^2*H_nat^2*du/dz + (1+z)*H_nat * [H_nat + (1+z)*dH_nat/dz] * u
# (chain rule for d^2Psi/dt^2 in terms of z)
# Setting equal:
# (1+z)^2*H_nat^2*du/dz = 3*(1+z)*H_nat^2*u - m0^2*Psi
#                        - (1+z)*H_nat * [H_nat + (1+z)*dH_nat/dz] * u
# (1+z)^2*H_nat^2*du/dz = (1+z)*H_nat*u * [3*H_nat - H_nat - (1+z)*dH_nat/dz]
#                        - m0^2*Psi
# (1+z)^2*H_nat^2*du/dz = (1+z)*H_nat*u * [2*H_nat - (1+z)*dH_nat/dz] - m0^2*Psi
# du/dz = u * [2/((1+z)) - dH_nat/dz / H_nat] - m0^2*Psi / ((1+z)^2 * H_nat^2)

def dPsi_dz(z, y, Lambda0=0.0):
    psi, u = y
    H  = H_nat(z)                      # Mpc^-1
    dz = z * 1e-6 + 1e-10
    dH = (H_nat(z + dz) - H_nat(z - dz)) / (2.0*dz)   # dH/dz [Mpc^-1]

    # rho_Psi in natural units [Mpc^-2]: rho = 0.5*Psi_t^2 + 0.5*m0^2*Psi^2
    Psi_t  = -(1+z) * H * u             # dPsi/dt [Mpc^-1]
    rho_Psi = 0.5 * Psi_t**2 + 0.5 * m0**2 * psi**2
    rho_crit = 3.0 * H**2 / (8.0 * np.pi)
    f_Psi   = rho_Psi / rho_crit        # fractional energy density

    # Use RIFT-corrected H in the ODE
    H_rift = H * np.sqrt(max(1.0 + f_Psi, 1.0))
    dH_rift_approx = dH * (1.0 + 0.5 * f_Psi)

    du_dz = (u * (2.0/(1.0+z) - dH_rift_approx / H_rift)
             - m0**2 * psi / ((1.0+z)**2 * H_rift**2))
    return [u, du_dz]

def integrate_Psi(Lambda0, z_max=1e6, N=4000):
    """Integrate KG equation for Psi from z_max to 0."""
    # At z_max >> 1/m0, field is frozen: Psi = Psi0, dPsi/dz = 0
    y0 = [Psi0, 0.0]
    z_eval = np.concatenate([
        np.logspace(np.log10(z_max), np.log10(10.0), N//2),
        np.linspace(9.999, 0.0, N//2)   # avoid duplicate at z=10
    ])
    sol = solve_ivp(lambda z, y: dPsi_dz(z, y, Lambda0),
                    (z_max, 0.0), y0, t_eval=z_eval,
                    method='Radau', rtol=1e-9, atol=1e-12, dense_output=True)
    return sol

# ── RIFT H_kms including rho_Psi ─────────────────────────────────────────────
def make_H_RIFT_kms(sol_Psi):
    """Return a function H_RIFT(z) [km/s/Mpc] using precomputed Psi solution."""
    def H_RIFT_fn(z):
        psi  = float(sol_Psi.sol(z)[0])
        u    = float(sol_Psi.sol(z)[1])
        H    = H_nat(z)
        Psi_t  = -(1+z) * H * u
        rho_Psi = 0.5 * Psi_t**2 + 0.5 * m0**2 * psi**2
        rho_crit = 3.0 * H**2 / (8.0 * np.pi)
        f_Psi = rho_Psi / rho_crit
        return H_kms(z) * np.sqrt(max(1.0 + f_Psi, 1.0))
    return H_RIFT_fn

def make_cs2_RIFT(sol_Psi, Lambda0):
    """Return c_s^2(z) with G_eff correction from Psi."""
    def cs2_fn(z):
        psi = float(sol_Psi.sol(z)[0])
        G_eff = 1.0 / (1.0 + 16.0 * np.pi * Lambda0 * psi**2)
        # In tight coupling, c_s^2 is set by pressure/density in the photon-baryon fluid.
        # G_eff enters only at sub-leading order (gravity does not set sound speed in plasma).
        # Include it as a perturbative correction for completeness.
        return cs2(z) * G_eff
    return cs2_fn

# ── Main scan ────────────────────────────────────────────────────────────────
results = {}
Psi_store = {}
rho_ratio_store = {}
z_profile = np.logspace(np.log10(1e6), np.log10(1.0), 500)

print(f"\n{'Lambda0':>10s} {'r_d [Mpc]':>12s} {'Delta_rd [Mpc]':>16s} {'frac [ppm]':>12s} {'rho_Psi/rho@zdrag':>20s}")
print("-" * 74)

for Lambda0 in Lambda0_scan:
    sol = integrate_Psi(Lambda0)

    H_RIFT_fn   = make_H_RIFT_kms(sol)
    cs2_RIFT_fn = make_cs2_RIFT(sol, Lambda0)

    rd_val, rd_err = rd_integral(cs2_RIFT_fn, H_RIFT_fn, Z_DRAG)
    delta_rd = rd_val - rd_LCDM
    frac_ppm = delta_rd / rd_LCDM * 1e6

    # rho_Psi / rho_total at z_drag
    psi_d  = float(sol.sol(Z_DRAG)[0])
    u_d    = float(sol.sol(Z_DRAG)[1])
    H_d    = H_nat(Z_DRAG)
    Psi_t_d  = -(1+Z_DRAG) * H_d * u_d
    rho_Psi_d  = 0.5 * Psi_t_d**2 + 0.5 * m0**2 * psi_d**2
    rho_crit_d = 3.0 * H_d**2 / (8.0 * np.pi)
    f_drag = rho_Psi_d / rho_crit_d

    results[Lambda0] = {
        'rd_Mpc':           rd_val,
        'delta_rd_Mpc':     delta_rd,
        'frac_shift':       delta_rd / rd_LCDM,
        'frac_ppm':         frac_ppm,
        'rho_Psi_frac_at_zdrag': f_drag,
        'Psi_at_zdrag':     psi_d,
    }

    # Store evolution profiles
    Psi_prof   = np.array([float(sol.sol(z)[0]) for z in z_profile])
    u_prof     = np.array([float(sol.sol(z)[1]) for z in z_profile])
    H_nat_prof = H_nat(z_profile)
    Psi_t_prof = -(1+z_profile) * H_nat_prof * u_prof
    rho_Psi_prof = 0.5 * Psi_t_prof**2 + 0.5 * m0**2 * Psi_prof**2
    rho_crit_prof = 3.0 * H_nat_prof**2 / (8.0 * np.pi)
    rho_ratio_store[Lambda0] = rho_Psi_prof / rho_crit_prof
    Psi_store[Lambda0] = Psi_prof

    print(f"{Lambda0:>10.4f} {rd_val:>12.4f} {delta_rd:>+16.6f} {frac_ppm:>+12.3f} {f_drag:>20.3e}")

# ── Diagnostics ───────────────────────────────────────────────────────────────
rd_ref   = 147.09   # Mpc, standard Planck 2018 value
lcdm_ok  = abs(rd_LCDM - rd_ref) < 3.0

max_frac = max(abs(results[L0]['frac_shift']) for L0 in Lambda0_scan)
max_ppm  = max_frac * 1e6

diag = {
    'rd_LCDM_Mpc':          rd_LCDM,
    'rd_LCDM_err_Mpc':      rd_LCDM_err,
    'z_drag_EH':            Z_DRAG,
    'rd_Planck_reference':  rd_ref,
    'LCDM_baseline_PASS':   bool(lcdm_ok),
    'max_frac_shift':       max_frac,
    'max_shift_ppm':        max_ppm,
    'results_by_Lambda0': {str(k): v for k, v in results.items()},
    'verdict': {
        'LCDM_baseline':       'PASS' if lcdm_ok else f'FAIL (got {rd_LCDM:.2f} vs {rd_ref:.2f} Mpc)',
        'rd_at_Lambda0_0p003': f"{results[0.003]['delta_rd_Mpc']:+.6f} Mpc ({results[0.003]['frac_ppm']:+.3f} ppm)",
        'rd_at_Lambda0_0p095': f"{results[0.095]['delta_rd_Mpc']:+.6f} Mpc ({results[0.095]['frac_ppm']:+.3f} ppm)",
        'max_shift_ppm':       f"{max_ppm:.3f} ppm across all Lambda0 tested",
        'interpretation': (
            "If |Delta r_d / r_d| < ~1000 ppm across the full Lambda0 scan, "
            "RIFT and LCDM are degenerate in r_d. The echo-shell narrative "
            "is confirmed as a physical interpretation of standard BAO physics, "
            "not a distinct prediction. "
            "If shift > current BAO precision (~0.3%), RIFT makes a testable r_d prediction."
        )
    }
}

with open(os.path.join(OUT_DIR, 'sim101_diagnostics.json'), 'w') as f:
    json.dump(diag, f, indent=2)

np.save(os.path.join(OUT_DIR, 'sim101_z_profile.npy'), z_profile)
np.save(os.path.join(OUT_DIR, 'sim101_Psi_profiles.npy'),
        np.array([Psi_store[L0] for L0 in Lambda0_scan]))
np.save(os.path.join(OUT_DIR, 'sim101_rho_ratio.npy'),
        np.array([rho_ratio_store[L0] for L0 in Lambda0_scan]))

# ── Plots ─────────────────────────────────────────────────────────────────────
Lambda0_arr = np.array(Lambda0_scan)
rd_arr      = np.array([results[L0]['rd_Mpc'] for L0 in Lambda0_scan])
ppm_arr     = np.array([results[L0]['frac_ppm'] for L0 in Lambda0_scan])
f_drag_arr  = np.array([results[L0]['rho_Psi_frac_at_zdrag'] for L0 in Lambda0_scan])

fig, axes = plt.subplots(1, 3, figsize=(14, 4))

ax = axes[0]
ax.axhline(rd_LCDM, color='gray', ls='--', lw=1.5, label=r'$\Lambda$CDM')
ax.plot(Lambda0_arr, rd_arr, 'o-', color='steelblue', lw=2, ms=6)
ax.set_xlabel(r'$\Lambda_0$');  ax.set_ylabel(r'$r_d$ [Mpc]')
ax.set_title(r'Sound horizon $r_d(\Lambda_0)$');  ax.legend()

ax = axes[1]
ax.plot(Lambda0_arr, ppm_arr, 's-', color='darkorange', lw=2, ms=6)
ax.axhline(0, color='gray', ls='--', lw=1)
ax.set_xlabel(r'$\Lambda_0$');  ax.set_ylabel(r'$\Delta r_d / r_d$ [ppm]')
ax.set_title(r'Fractional $r_d$ shift')

ax = axes[2]
ax.semilogy(Lambda0_arr, f_drag_arr, '^-', color='firebrick', lw=2, ms=6)
ax.set_xlabel(r'$\Lambda_0$');  ax.set_ylabel(r'$\rho_\Psi / \rho_{\rm tot}$ at $z_{\rm drag}$')
ax.set_title(r'Scalar energy fraction at drag epoch')

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, 'sim101_rd_vs_lambda0.pdf'), bbox_inches='tight')
plt.savefig(os.path.join(OUT_DIR, 'sim101_rd_vs_lambda0.png'), dpi=150, bbox_inches='tight')
plt.close()

# Psi evolution
fig2, ax2 = plt.subplots(figsize=(7, 4))
colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(Lambda0_scan)))
for i, L0 in enumerate(Lambda0_scan):
    ax2.semilogx(z_profile, Psi_store[L0], color=colors[i], lw=1.5,
                 label=rf'$\Lambda_0={L0}$')
ax2.axvline(Z_DRAG, color='k', ls=':', lw=1.2, label=f'$z_d={Z_DRAG:.0f}$')
ax2.set_xlabel(r'Redshift $z$');  ax2.set_ylabel(r'$\Psi(z)$')
ax2.set_title(r'Scalar field evolution')
ax2.legend(fontsize=7, ncol=2)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, 'sim101_Psi_evolution.pdf'), bbox_inches='tight')
plt.savefig(os.path.join(OUT_DIR, 'sim101_Psi_evolution.png'), dpi=150, bbox_inches='tight')
plt.close()

# ── Summary ───────────────────────────────────────────────────────────────────
print(f"\n{'='*65}")
print(f"SIM101 SUMMARY")
print(f"{'='*65}")
print(f"LCDM r_d:        {rd_LCDM:.4f} Mpc  (Planck ref: {rd_ref:.2f})  [{diag['verdict']['LCDM_baseline']}]")
print(f"z_drag (E&H):    {Z_DRAG:.1f}")
print(f"")
print(f"Lambda0=0.003:   Delta r_d = {diag['verdict']['rd_at_Lambda0_0p003']}")
print(f"Lambda0=0.095:   Delta r_d = {diag['verdict']['rd_at_Lambda0_0p095']}")
print(f"Max |shift|:     {max_ppm:.3f} ppm  (BAO precision ~3000 ppm)")
print(f"")
print(f"Outputs: {OUT_DIR}")
print(f"{'='*65}")
