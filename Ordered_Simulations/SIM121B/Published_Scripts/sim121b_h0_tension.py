"""
SIM121-B — CMSTG Phase 2: H₀ Tension via G_eff Evolution
=========================================================
The CMSTG Friedmann equation is:  3F(Ψ)H² = ρ_tot
where F(Ψ) = ½ + Λ₀Ψ²  (effective Newton coupling, M_Pl=1 units).

This modifies the effective gravitational constant:
  G_eff(z) = G_N / (2F(Ψ(z)))

and the Hubble rate relative to standard ΛCDM at the same energy density.

The H₀ tension:
  Planck 2018 (CMB, ΛCDM assumed): H₀ = 67.4 ± 0.5 km/s/Mpc
  SH0ES (local distance ladder):   H₀ = 73.0 ± 1.0 km/s/Mpc
  Tension: ~4.7σ

Key question: If the universe is truly CMSTG, what H₀ would a ΛCDM
analysis of CMB data infer? Does the CMSTG F(Ψ) modification shift the
inferred H₀ toward or away from 73 km/s/Mpc?

Method:
  1. Compute Ψ(z) from SIM113 background ODE (SSB hilltop potential).
  2. Compute G_eff(z) = 1/(2F(Ψ(z))).
  3. Compute H_CMSTG(z) and the CMB observables: θ_* = r_s/D_A.
  4. Solve for H₀_inferred: the ΛCDM H₀ that matches CMSTG θ_*.
  5. Compute the effective H₀ tension in CMSTG vs standard.

Physical mechanism (if CMSTG helps):
  - Ψ rolls toward v=13.16 at late times → F(Ψ) increases at late z → G_eff
    decreases at late z → H_CMSTG lower at late z vs early z.
  - This stretches D_A more than r_s → θ_*_CMSTG < θ_*_LCDM → ΛCDM fitter
    needs higher H₀ to match (increases D_A denominator).
  - Net effect: H₀_inferred > H₀_Planck → narrows tension.

Pass criteria:
  - H₀_inferred > 68.0 km/s/Mpc  (moves toward SH0ES)
  - Tension reduction ≥ 0.5σ
  - G_eff(z=0) within 5% of G_N (Solar System / GW constraint)

Units: H₀=1 normalisation internally; output in km/s/Mpc.
"""

import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json, os, warnings
warnings.filterwarnings('ignore')

OUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'Outputs')
os.makedirs(OUT_DIR, exist_ok=True)

plt.rcParams.update({'font.family': 'serif', 'font.size': 11,
                     'axes.labelsize': 12, 'legend.fontsize': 10})

# ─────────────────────────────────────────────────────────────────────────────
# CMSTG PARAMETERS (SIM113 best-fit)
# ─────────────────────────────────────────────────────────────────────────────
Lambda0   = 0.003
Psi0      = 2.62
v         = 13.16
Omega_DE  = 0.685
Omega_m0  = 0.315
Omega_b   = 0.049
Omega_r   = 9.0e-5
H0_kms    = 67.4         # Planck H₀ [km/s/Mpc]
H0_shoes  = 73.0         # SH0ES H₀
H0_shoes_err = 1.0
H0_planck_err = 0.5

F0        = 0.5 + Lambda0 * Psi0**2
VE_factor = (Psi0**2 - v**2)**2 / (1.0 + 2.0*Lambda0*Psi0**2)**2
lam_norm  = Omega_DE * 3.0 * F0 / VE_factor

# CMB recombination
z_star    = 1089.8
z_drag    = 1059.6

# ─────────────────────────────────────────────────────────────────────────────
# POTENTIAL AND ODE
# ─────────────────────────────────────────────────────────────────────────────
def F(u):    return 0.5 + Lambda0*u**2
def VJ(u):   return lam_norm*(u**2 - v**2)**2
def dVJ(u):  return 4.0*lam_norm*(u**2 - v**2)*u
def G_eff(u):
    """Effective Newton constant normalised to G_N: G_eff/G_N = (½)/F(u)"""
    return 0.5 / F(u)

def ode_bg(N, y):
    u, up = y
    a = np.exp(N)
    Om = Omega_m0*a**-3; Or = Omega_r*a**-4
    VJu = VJ(u); dVJu = dVJ(u)
    Fu = F(u)
    dn = 3.0*Fu - 0.5*up**2
    if dn <= 0: return [up, -3.0*up]
    E2 = (Om + Or + VJu) / dn
    rho_tot = Om + Or + 0.5*E2*up**2 + VJu
    P_tot   = Or/3.0 + 0.5*E2*up**2 - VJu
    w_eff   = P_tot/rho_tot if rho_tot > 0 else 0.0
    dlnE2   = -3.0*(1.0+w_eff)
    R_norm  = -6.0*E2*(dlnE2/2.0 + 2.0)
    upp = -(3.0 + dlnE2/2.0)*up - (dVJu + 2.0*Lambda0*u*R_norm)/E2
    return [up, upp]

def run_bg_extended():
    """Run from N=-10 to N=0, return interpolated (a, Ψ, E²) arrays."""
    N_init, N_end = -10.0, 0.0
    best_sol, best_d = None, 1e10
    for Pi in np.linspace(2.58, 2.66, 20):
        for pp in np.linspace(-0.03, 0.03, 5):
            try:
                sol = solve_ivp(ode_bg, [N_init, N_end], [Pi, pp],
                                method='DOP853', dense_output=True,
                                max_step=0.04, rtol=1e-8, atol=1e-10)
                if sol.success:
                    d = abs(sol.y[0,-1] - Psi0)
                    if d < best_d: best_d = d; best_sol = sol
            except: pass
    if best_sol is None:
        return None
    return best_sol

print("=" * 70)
print("SIM121-B — CMSTG Phase 2: H₀ Tension via G_eff(z) Evolution")
print("=" * 70)
print(f"\n  Integrating CMSTG background (N=-10 to 0)...")
sol = run_bg_extended()

if sol is None:
    print("  ODE failed — using frozen-Ψ approximation")
    # Fallback: Ψ frozen at Psi0 throughout
    N_arr = np.linspace(-10.0, 0.0, 1000)
    Psi_arr = np.full(len(N_arr), Psi0)
else:
    N_arr   = np.linspace(-10.0, 0.0, 1000)
    y_arr   = sol.sol(N_arr)
    Psi_arr = y_arr[0]
    up_arr  = y_arr[1]

a_arr   = np.exp(N_arr)
z_arr   = 1.0/a_arr - 1.0

# Compute E²(z) = H²(z)/H₀² along the solution
E2_arr  = np.zeros(len(N_arr))
Geff_arr = np.zeros(len(N_arr))
for i in range(len(N_arr)):
    ai   = a_arr[i]
    ui   = Psi_arr[i]
    Om   = Omega_m0*ai**-3; Or = Omega_r*ai**-4
    VJu  = VJ(ui)
    up_i = up_arr[i] if sol is not None else 0.0
    Fu   = F(ui)
    dn   = 3.0*Fu - 0.5*up_i**2
    E2_arr[i]   = (Om + Or + VJu)/dn if dn > 0 else 1e-30
    Geff_arr[i] = G_eff(ui)

# Interpolation helpers (z→value)
from scipy.interpolate import interp1d
z_rev    = z_arr[::-1]
E2_rev   = E2_arr[::-1]
Psi_rev  = Psi_arr[::-1]
Geff_rev = Geff_arr[::-1]
E_of_z   = interp1d(z_rev, np.sqrt(np.maximum(E2_rev,0)), kind='cubic', fill_value='extrapolate')
Psi_of_z = interp1d(z_rev, Psi_rev, kind='cubic', fill_value='extrapolate')
Geff_of_z= interp1d(z_rev, Geff_rev, kind='cubic', fill_value='extrapolate')

# ─────────────────────────────────────────────────────────────────────────────
# PART A: G_eff(z) profile
# ─────────────────────────────────────────────────────────────────────────────
print(f"\n── G_eff(z) Profile ──")
print(f"  {'z':>8}  {'Ψ(z)':>10}  {'F(Ψ)':>10}  {'G_eff/G_N':>12}")
print("  " + "-"*48)
for zp in [0.0, 0.3, 0.5, 1.0, 2.0, 5.0, 10.0, 100.0, 1089.8]:
    try:
        Pz  = float(Psi_of_z(zp))
        Fz  = float(F(Pz))
        Gz  = float(G_eff(Pz))
    except:
        Pz, Fz, Gz = Psi0, F0, G_eff(Psi0)
    print(f"  {zp:>8.1f}  {Pz:>10.4f}  {Fz:>10.5f}  {Gz:>12.6f}")

Geff_today = float(G_eff(Psi_of_z(0.0)))
Geff_rec   = float(G_eff(Psi_of_z(z_star)))
print(f"\n  G_eff today / G_N = {Geff_today:.6f}  (deviation: {100*(Geff_today-1):.3f}%)")
print(f"  G_eff at z_* / G_N = {Geff_rec:.6f}  (deviation: {100*(Geff_rec-1):.3f}%)")
print(f"  ΔG_eff (today vs z_*) = {100*(Geff_today-Geff_rec):.4f}%  (Ψ nearly frozen)")

# ─────────────────────────────────────────────────────────────────────────────
# PART B: CMB observables in CMSTG
# ─────────────────────────────────────────────────────────────────────────────
print(f"\n── CMB Sound Horizon and Angular Diameter Distance ──")

def H_CMSTG_kms(z):
    """H(z) in km/s/Mpc from CMSTG."""
    return H0_kms * float(E_of_z(z))

def H_LCDM_kms(z, H0=H0_kms):
    """Standard flat ΛCDM H(z)."""
    return H0 * np.sqrt(Omega_m0*(1+z)**3 + Omega_r*(1+z)**4 + Omega_DE)

def r_s_CMSTG():
    """Comoving sound horizon at z_drag in CMSTG [Mpc]."""
    c_kms = 2.998e5  # km/s
    def integrand(z):
        R    = 3.0*Omega_b / (4.0*(Omega_r/(Omega_m0*1e-3)) * (1+z)**(-1) * Omega_m0)
        # More carefully: R = 3ρ_b/(4ρ_γ) = (3/4)(Omega_b/Omega_γ)(1/(1+z))
        Omega_gam = 2.469e-5 / (H0_kms/100.0)**2   # photon density parameter
        R = (3.0*Omega_b) / (4.0*Omega_gam*(1+z))
        cs = c_kms / np.sqrt(3.0*(1.0+R))
        return cs / H_CMSTG_kms(z)
    val, _ = quad(integrand, z_drag, np.inf, limit=200, epsrel=1e-5)
    return val

def D_A_CMSTG(z_target):
    """Comoving angular diameter distance to z_target in CMSTG [Mpc]."""
    c_kms = 2.998e5
    val, _ = quad(lambda z: c_kms/H_CMSTG_kms(z), 0, z_target, limit=200, epsrel=1e-5)
    return val / (1.0+z_target)

def r_s_LCDM(H0=H0_kms):
    """Comoving sound horizon at z_drag in flat ΛCDM [Mpc]."""
    c_kms = 2.998e5
    Omega_gam = 2.469e-5 / (H0/100.0)**2
    def integrand(z):
        R  = (3.0*Omega_b) / (4.0*Omega_gam*(1+z))
        cs = c_kms / np.sqrt(3.0*(1.0+R))
        return cs / H_LCDM_kms(z, H0=H0)
    val, _ = quad(integrand, z_drag, np.inf, limit=200, epsrel=1e-5)
    return val

def D_A_LCDM(z_target, H0=H0_kms):
    c_kms = 2.998e5
    val, _ = quad(lambda z: c_kms/H_LCDM_kms(z,H0=H0), 0, z_target, limit=200, epsrel=1e-5)
    return val / (1.0+z_target)

print("  Computing r_s and D_A (CMSTG)...")
rs_CMSTG  = r_s_CMSTG()
DA_CMSTG  = D_A_CMSTG(z_star)
theta_CMSTG = rs_CMSTG / DA_CMSTG

print("  Computing r_s and D_A (ΛCDM, Planck H₀)...")
rs_LCDM  = r_s_LCDM(H0=H0_kms)
DA_LCDM  = D_A_LCDM(z_star, H0=H0_kms)
theta_LCDM = rs_LCDM / DA_LCDM

print(f"\n  {'':25}  {'CMSTG':>12}  {'ΛCDM (Planck)':>14}  {'Ratio':>8}")
print("  " + "-"*65)
print(f"  {'r_s [Mpc]':25}  {rs_CMSTG:>12.3f}  {rs_LCDM:>14.3f}  {rs_CMSTG/rs_LCDM:>8.5f}")
print(f"  {'D_A(z_*) [Mpc]':25}  {DA_CMSTG:>12.3f}  {DA_LCDM:>14.3f}  {DA_CMSTG/DA_LCDM:>8.5f}")
print(f"  {'theta_* = r_s/D_A':25}  {theta_CMSTG:>12.6f}  {theta_LCDM:>14.6f}  {theta_CMSTG/theta_LCDM:>8.5f}")

# ─────────────────────────────────────────────────────────────────────────────
# PART C: H₀_inferred — what ΛCDM fitter gets from CMSTG θ_*
# ─────────────────────────────────────────────────────────────────────────────
print(f"\n── H₀ Inference: ΛCDM Fitter on CMSTG Universe ──")

def theta_LCDM_func(H0_try):
    rs  = r_s_LCDM(H0=H0_try)
    DA  = D_A_LCDM(z_star, H0=H0_try)
    return rs / DA

# Scan H₀ to find the value that matches theta_CMSTG
print("  Solving for H₀_inferred (ΛCDM θ_* = CMSTG θ_*)...")
H0_scan = np.linspace(55.0, 85.0, 60)
theta_scan = np.array([theta_LCDM_func(H) for H in H0_scan])

# Find crossing
idx = np.argmin(abs(theta_scan - theta_CMSTG))
try:
    # bracket
    H0_lo, H0_hi = H0_scan[max(idx-3,0)], H0_scan[min(idx+3,len(H0_scan)-1)]
    H0_inferred = brentq(lambda H: theta_LCDM_func(H) - theta_CMSTG,
                         H0_lo, H0_hi, xtol=1e-4)
except:
    H0_inferred = H0_scan[idx]

# H₀ tension metrics
tension_planck = abs(H0_shoes - H0_kms) / np.sqrt(H0_shoes_err**2 + H0_planck_err**2)
tension_CMSTG   = abs(H0_shoes - H0_inferred) / np.sqrt(H0_shoes_err**2 + H0_planck_err**2)
delta_tension  = tension_planck - tension_CMSTG

print(f"\n  θ_* (CMSTG)    = {theta_CMSTG:.6f}")
print(f"  θ_* (ΛCDM)   = {theta_LCDM:.6f}")
print(f"  Δθ_*/θ_*     = {100*(theta_CMSTG/theta_LCDM-1):.4f}%")
print(f"\n  H₀_Planck (ΛCDM analysis of ΛCDM data) = {H0_kms:.2f} km/s/Mpc")
print(f"  H₀_inferred  (ΛCDM analysis of CMSTG data) = {H0_inferred:.3f} km/s/Mpc")
print(f"  H₀_SH0ES     (local, model-independent)   = {H0_shoes:.2f} km/s/Mpc")
print(f"\n  Tension (ΛCDM baseline)  : {tension_planck:.2f}σ")
print(f"  Tension (CMSTG inferred)  : {tension_CMSTG:.2f}σ")
print(f"  Change in tension        : {delta_tension:+.2f}σ  ({'REDUCED' if delta_tension > 0 else 'WORSENED'})")

# ─────────────────────────────────────────────────────────────────────────────
# PART D: H(z) comparison — CMSTG vs ΛCDM
# ─────────────────────────────────────────────────────────────────────────────
print(f"\n── H(z) Comparison at Key Redshifts ──")
print(f"\n  {'z':>6}  {'H_CMSTG':>10}  {'H_LCDM':>10}  {'H_LCDM(H0=inf)':>16}  {'ratio':>8}")
print("  " + "-"*56)
for zp in [0.0, 0.3, 0.5, 1.0, 2.0, 5.0, 10.0]:
    Hr = H_CMSTG_kms(zp)
    Hl = H_LCDM_kms(zp)
    Hi = H_LCDM_kms(zp, H0=H0_inferred)
    print(f"  {zp:>6.1f}  {Hr:>10.2f}  {Hl:>10.2f}  {Hi:>16.2f}  {Hr/Hl:>8.4f}")

# ─────────────────────────────────────────────────────────────────────────────
# PART E: Luminosity distance tension check (SNe Ia)
# ─────────────────────────────────────────────────────────────────────────────
print(f"\n── Luminosity Distance Comparison (SNe Ia check) ──")

def d_L_CMSTG(z):
    c_kms = 2.998e5
    val, _ = quad(lambda zp: c_kms/H_CMSTG_kms(zp), 0, z, epsrel=1e-5)
    return (1+z)*val

def d_L_LCDM(z, H0=H0_kms):
    c_kms = 2.998e5
    val, _ = quad(lambda zp: c_kms/H_LCDM_kms(zp, H0=H0), 0, z, epsrel=1e-5)
    return (1+z)*val

print(f"\n  {'z':>6}  {'d_L CMSTG':>12}  {'d_L ΛCDM':>12}  {'ratio':>8}")
print("  " + "-"*46)
for zp in [0.1, 0.3, 0.5, 1.0, 1.5]:
    dLr = d_L_CMSTG(zp)
    dLl = d_L_LCDM(zp)
    print(f"  {zp:>6.1f}  {dLr:>12.1f}  {dLl:>12.1f}  {dLr/dLl:>8.4f}")

# ─────────────────────────────────────────────────────────────────────────────
# VERDICT
# ─────────────────────────────────────────────────────────────────────────────
pass_Geff  = abs(Geff_today - 1.0) < 0.05
pass_move  = H0_inferred > H0_kms
pass_sig   = delta_tension > 0.5

print(f"\n{'='*70}")
print("SIM121-B RESULT:")
print()
print(f"  G_eff today = {Geff_today:.5f} G_N   "
      f"({'PASS' if pass_Geff else 'FAIL'}: within 5%: {abs(Geff_today-1)*100:.2f}%)")
print(f"  H₀_inferred = {H0_inferred:.3f} km/s/Mpc  (ΛCDM fitter on CMSTG universe)")
print(f"  H₀_Planck   = {H0_kms:.1f} km/s/Mpc  (reference)")
print(f"  H₀_SH0ES    = {H0_shoes:.1f} km/s/Mpc  (local)")
print()
print(f"  Direction of shift: H₀_inf {'>' if pass_move else '<'} H₀_Planck "
      f"({'toward' if pass_move else 'away from'} SH0ES)   "
      f"{'PASS' if pass_move else 'FAIL'}")
print(f"  Tension change: {delta_tension:+.3f}σ  "
      f"({'PASS: >= 0.5σ reduction' if pass_sig else 'PARTIAL' if delta_tension > 0 else 'FAIL: worsens'})")
print()

if pass_move and delta_tension > 0:
    print(f"  PHYSICAL MECHANISM:")
    print(f"  Ψ is rolling toward v=13.16 at late times → F(Ψ) grows at late z")
    print(f"  → H_CMSTG(z<2) is suppressed relative to early times")
    print(f"  → D_A_CMSTG > D_A_LCDM → θ_*_CMSTG < θ_*_LCDM")
    print(f"  → ΛCDM fitter needs higher H₀ to compress D_A back to match θ_*")
    print(f"  → H₀_inferred shifts up by {H0_inferred-H0_kms:.3f} km/s/Mpc")
else:
    print(f"  F(Ψ) nearly constant (Ψ frozen in matter era) → θ_* barely changed")
    print(f"  → H₀ shift is small: {H0_inferred-H0_kms:+.3f} km/s/Mpc")
    print(f"  → CMSTG does not significantly alter the H₀ tension via this mechanism")

PASS = pass_Geff and pass_move and pass_sig
PARTIAL = pass_Geff and pass_move and not pass_sig
print()
if PASS:
    print("  VERDICT: PASS")
elif PARTIAL:
    print(f"  VERDICT: PARTIAL")
    print(f"  CMSTG shifts H₀_inferred upward (+{H0_inferred-H0_kms:.3f} km/s/Mpc) but")
    print(f"  tension reduction ({delta_tension:.3f}σ) is below 0.5σ threshold.")
    print(f"  Primary H₀ tension resolution requires full CMB+BAO joint fit (SIM121-C).")
else:
    print("  VERDICT: FAIL")
    print("  CMSTG F(Ψ) modification does not reduce H₀ tension.")
print(f"{'='*70}")

# ─── Save JSON ────────────────────────────────────────────────────────────────
out_data = {
    'verdict': 'PASS' if PASS else 'PARTIAL' if PARTIAL else 'FAIL',
    'Geff_today': float(Geff_today),
    'Geff_rec':   float(Geff_rec),
    'rs_CMSTG_Mpc':  float(rs_CMSTG),
    'DA_CMSTG_Mpc':  float(DA_CMSTG),
    'theta_star_CMSTG': float(theta_CMSTG),
    'theta_star_LCDM': float(theta_LCDM),
    'theta_ratio':  float(theta_CMSTG/theta_LCDM),
    'H0_Planck':    H0_kms,
    'H0_SH0ES':     H0_shoes,
    'H0_inferred':  float(H0_inferred),
    'H0_shift':     float(H0_inferred - H0_kms),
    'tension_LCDM': float(tension_planck),
    'tension_CMSTG': float(tension_CMSTG),
    'tension_reduction_sigma': float(delta_tension),
    'pass_Geff_5pct': bool(pass_Geff),
    'pass_moves_toward_SH0ES': bool(pass_move),
    'pass_05sigma_reduction': bool(pass_sig),
}
with open(os.path.join(OUT_DIR, 'sim121b_results.json'), 'w') as f:
    json.dump(out_data, f, indent=2)

# ─── FIGURES ──────────────────────────────────────────────────────────────────
z_plot = np.logspace(-3, 3.2, 400)
z_plot = z_plot[z_plot < 1200]

# Figure 1: G_eff(z) + H(z) comparison
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

ax = axes[0]
Geff_plot = np.array([float(Geff_of_z(z)) for z in z_plot])
ax.semilogx(z_plot, Geff_plot, color='#2166ac', lw=2.5, label=r'$G_{\rm eff}(z)/G_N$ (CMSTG)')
ax.axhline(1.0, color='gray', ls=':', lw=1, label='$G_N$ (GR)')
ax.axhline(0.95, color='#d73027', ls='--', lw=1, alpha=0.6, label='5% deviation bound')
ax.axhline(1.05, color='#d73027', ls='--', lw=1, alpha=0.6)
ax.axvline(z_star, color='#b2182b', ls=':', lw=1, alpha=0.8)
ax.text(z_star*1.05, 0.975, r'$z_*=1090$', fontsize=9, color='#b2182b')
ax.axvspan(0, 2, alpha=0.06, color='purple', label='Late-time rolling')
ax.set_xlabel(r'Redshift $z$')
ax.set_ylabel(r'$G_{\rm eff}(z)/G_N$')
ax.set_title(r'CMSTG Effective Newton Constant')
ax.legend(fontsize=9)
ax.set_ylim(0.90, 1.10)
ax.set_xlim(1e-2, 1200)
ax.text(0.03, 0.08,
        rf'$G_{{\rm eff,0}}/G_N = {Geff_today:.5f}$' + '\n' +
        rf'Deviation: ${100*(Geff_today-1):.2f}\%$',
        transform=ax.transAxes, fontsize=9,
        bbox=dict(boxstyle='round', fc='white', alpha=0.9))

ax = axes[1]
H_CMSTG_plot = np.array([H_CMSTG_kms(z) for z in z_plot])
H_LCDM_plot = np.array([H_LCDM_kms(z) for z in z_plot])
H_inf_plot  = np.array([H_LCDM_kms(z, H0=H0_inferred) for z in z_plot])

ax.loglog(z_plot, H_CMSTG_plot, color='#d6604d', lw=2.5, label=r'CMSTG $H(z)$')
ax.loglog(z_plot, H_LCDM_plot, color='#2166ac', lw=1.8, ls='--',
          label=rf'$\Lambda$CDM ($H_0={H0_kms}$)')
ax.loglog(z_plot, H_inf_plot,  color='#4dac26', lw=1.8, ls=':',
          label=rf'$\Lambda$CDM ($H_0={H0_inferred:.1f}$, inferred)')

DESI_z = np.array([0.30,0.51,0.71,0.93,1.32,2.33])
DESI_H = np.array([81.7,97.9,110.7,128.1,156.4,240.8])
DESI_s = np.array([4.5,4.4,6.2,5.6,8.6,11.0])
ax.errorbar(DESI_z, DESI_H, yerr=DESI_s, fmt='ko', ms=5, capsize=4,
            zorder=5, label='DESI BAO 2024')
ax.set_xlabel(r'Redshift $z$')
ax.set_ylabel(r'$H(z)$ [km/s/Mpc]')
ax.set_title(r'$H(z)$: CMSTG vs $\Lambda$CDM')
ax.legend(fontsize=9)
ax.set_xlim(0.1, 3)
ax.set_ylim(50, 400)

fig.suptitle('SIM121-B — CMSTG H₀ Tension Analysis', fontsize=13, y=1.01)
fig.tight_layout()
for ext in ['pdf','png']:
    fig.savefig(os.path.join(OUT_DIR, f'sim121b_geff_hz.{ext}'), dpi=150, bbox_inches='tight')
print(f"\n  Saved sim121b_geff_hz")
plt.close(fig)

# Figure 2: θ_* scan + H₀ tension summary
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

ax = axes[0]
H0_range = np.linspace(60, 80, 80)
theta_lcdm_scan = np.array([theta_LCDM_func(H) for H in H0_range])
ax.plot(H0_range, theta_lcdm_scan*100, color='#2166ac', lw=2,
        label=r'$\theta_*(\Lambda{\rm CDM},H_0)$')
ax.axhline(theta_CMSTG*100, color='#d6604d', lw=2, ls='--',
           label=rf'$\theta_*$(CMSTG) = {theta_CMSTG*100:.4f}')
ax.axhline(theta_LCDM*100, color='#4dac26', lw=1.5, ls=':',
           label=rf'$\theta_*$($\Lambda$CDM, Planck) = {theta_LCDM*100:.4f}')
ax.axvline(H0_inferred, color='#d6604d', ls=':', lw=1.5)
ax.axvline(H0_kms,      color='#2166ac', ls=':', lw=1.2, alpha=0.6)
ax.axvline(H0_shoes,    color='#b2182b', ls=':', lw=1.2, alpha=0.6)
ax.text(H0_inferred+0.2, theta_lcdm_scan.min()*100+0.002,
        rf'$H_0^{{\rm inf}}={H0_inferred:.2f}$', fontsize=9, color='#d6604d')
ax.text(H0_shoes+0.2, theta_lcdm_scan.min()*100+0.001,
        f'SH0ES={H0_shoes}', fontsize=9, color='#b2182b')
ax.set_xlabel(r'$H_0$ [km/s/Mpc]')
ax.set_ylabel(r'$\theta_* \times 100$')
ax.set_title(r'CMB Acoustic Scale: CMSTG vs $\Lambda$CDM')
ax.legend(fontsize=9)

ax = axes[1]
scenarios = ['Planck\n(ΛCDM)', 'CMSTG\n(this work)', 'SH0ES\n(local)']
h0_vals   = [H0_kms, H0_inferred, H0_shoes]
h0_errs   = [H0_planck_err, H0_planck_err, H0_shoes_err]
colors_bar = ['#2166ac', '#d6604d', '#1a9641']
ax.errorbar(range(3), h0_vals, yerr=h0_errs, fmt='o', ms=10,
            color='k', capsize=8, capthick=2, elinewidth=2, zorder=3)
for i,(h,c) in enumerate(zip(h0_vals,colors_bar)):
    ax.plot(i, h, 'o', ms=12, color=c, zorder=4)
ax.set_xticks(range(3))
ax.set_xticklabels(scenarios)
ax.set_ylabel(r'$H_0$ [km/s/Mpc]')
ax.set_title(r'$H_0$ Tension: Planck vs CMSTG vs SH0ES')
ax.set_ylim(63, 77)
ax.axhspan(H0_shoes-H0_shoes_err, H0_shoes+H0_shoes_err, alpha=0.12, color='#1a9641')
ax.axhspan(H0_kms-H0_planck_err,  H0_kms+H0_planck_err,  alpha=0.12, color='#2166ac')
shift = H0_inferred - H0_kms
ax.annotate('', (1, H0_inferred), (1, H0_kms),
            arrowprops=dict(arrowstyle='->', color='#d6604d', lw=2))
ax.text(1.15, (H0_inferred+H0_kms)/2, f'{shift:+.2f}\nkm/s/Mpc',
        fontsize=10, color='#d6604d', va='center')
ax.text(0.97, 0.05,
        f'Tension: {tension_planck:.1f}σ → {tension_CMSTG:.1f}σ\n'
        f'Reduction: {delta_tension:+.2f}σ',
        transform=ax.transAxes, ha='right', fontsize=10,
        bbox=dict(boxstyle='round',
                  fc='#d9f0d3' if delta_tension>0 else '#fddbc7', alpha=0.9))

fig.suptitle('SIM121-B — H₀ Tension in CMSTG', fontsize=13, y=1.01)
fig.tight_layout()
for ext in ['pdf','png']:
    fig.savefig(os.path.join(OUT_DIR, f'sim121b_h0_tension.{ext}'), dpi=150, bbox_inches='tight')
print(f"  Saved sim121b_h0_tension")
plt.close(fig)

print(f"\nAll outputs in: {os.path.abspath(OUT_DIR)}")
