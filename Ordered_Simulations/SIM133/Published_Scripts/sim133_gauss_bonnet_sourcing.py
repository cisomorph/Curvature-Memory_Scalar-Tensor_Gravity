"""
SIM133 — CMSTG Phase 3: Gauss-Bonnet Sourcing (Direction B)
===========================================================
SIM131/132 show that sourcing Ψ from R always drives F_eff the wrong way.
Alternative: source Ψ from the Gauss-Bonnet term 𝒢 instead.

Key physics: in FRW, 𝒢 = 24H²(Ḣ + H²) = 24H⁴(1−ε_H)
  - Matter domination: ε_H = 3/2 → (1−ε_H) = −½ → 𝒢 < 0 → Ψ driven negative
  - Λ domination:      ε_H → 0  → (1−ε_H) = +1  → 𝒢 > 0 → Ψ turns around

𝒢 CHANGES SIGN at matter-Λ equality (z ≈ 0.3) — exactly where DESI measures.

With Ψ negative at intermediate z: F_eff = ½(1+2Λ₀Ψ²)+ξΨ dips below Phase 1
via the ξΨ < 0 term → G_eff spikes → H(z) enhanced right where DESI wants it.

Action:
  S = ∫d⁴x√(-g) [ ½(1+2Λ₀Ψ²)R − ½(∂Ψ)² + α_GB·Ψ·𝒢 ] + S_SM
  F_eff = ½(1+2Λ₀Ψ²)  (Phase 1 backbone, unchanged)

Field EOM (FRW, N = ln a):
  Ψ'' + (3−ε_H)Ψ' = 2Λ₀Ψ·R/H² + α_GB·𝒢/H²

where:
  R/H²  = 6(2−ε_H)
  𝒢/H²  = 24H²(1−ε_H)  [dimensional: H² in H100 units]

Modified Friedmann: same Phase 1 form
  H²[3F_eff + 6·2Λ₀Ψ·Ψ' − ½Ψ'²] = ω_m a⁻³ + ω_r a⁻⁴ + Λ_bare

Note on α_GB units:
  □Ψ has units M_Pl·H². R·Ψ term: [Λ₀]·[Ψ]·[H²] ✓ (Λ₀ dimensionless).
  𝒢 has units H⁴. So [α_GB·𝒢] must have units M_Pl·H² → [α_GB] = M_Pl/H² = 1/H₀²
  In ODE (N-deriv): α_GB·𝒢/H² = α_GB·24H²(1−ε_H) → define α̃ = α_GB·H₀² (dimensionless)
  Then: α_GB·𝒢/H² = α̃·24·(H/H₀)²·(1−ε_H) = α̃·24·E²·(1−ε_H)

CMB risk: GB source ∝ E²·(1−ε_H). At z=1089: E ~ 1090^{1.5}/h² >> 1.
  Must use very small α̃ to keep |Ψ(z_CMB)| small.
  Scan: α̃ ∈ [1e-8, 1e-7, 1e-6, 1e-5, 1e-4]

Baseline: SIM121C DESI tension = 2.77σ, Λ₀ = 0.003 (Phase 1 locked)
"""

import numpy as np
from scipy.integrate import quad, solve_ivp
import json, os, warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

OUT = os.path.join(os.path.dirname(__file__), '..', 'Outputs')
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({'font.family': 'serif', 'font.size': 11,
                     'axes.labelsize': 12, 'legend.fontsize': 10})

# ─── Constants ───────────────────────────────────────────────────────────────
c_kms    = 2.998e5
H100     = 100.0
omh2_m   = 0.1430
omh2_r   = 4.18e-5
h_target = 0.674
omh2_b   = 0.02237

theta_obs     = 1.04101
theta_obs_err = 0.00029
z_drag        = 1059.6
z_star        = 1089.8

DESI_z = np.array([0.295, 0.510, 0.706, 0.930, 1.317, 2.330])
DESI_H = np.array([ 81.7,  97.9, 110.7, 128.1, 156.4, 240.8])
DESI_s = np.array([  4.5,   4.4,   6.2,   5.6,   8.6,  11.0])

SIM121C_chi2    = 18.26
SIM121C_tension = 2.77

Lambda0  = 0.003   # Phase 1 locked value
XI       = 1.0/6.0

F0_phase1    = 0.5 + Lambda0 * 2.62**2
Lbare_phase1 = 3.0*F0_phase1*h_target**2 - omh2_m - omh2_r

# ─── Physics functions ───────────────────────────────────────────────────────

def F_eff(Psi):
    """Phase 1 backbone — unchanged"""
    return 0.5 + Lambda0*Psi**2

def F_eff_prime(Psi):
    return 2.0*Lambda0*Psi

def get_H2(Psi, y, N, Lbare):
    a    = np.exp(N)
    rhs  = omh2_m/a**3 + omh2_r/a**4 + Lbare
    Fv   = F_eff(Psi)
    Fp   = F_eff_prime(Psi)
    coef = 3.0*Fv + 6.0*Fp*y - 0.5*y**2
    if coef < 0.3:
        coef = 3.0*Fv
    return rhs/coef

def get_eps_H(Psi, y, N, Lbare):
    eps = 5e-4
    H2p = get_H2(Psi, y, N+eps, Lbare)
    H2m = get_H2(Psi, y, N-eps, Lbare)
    H2  = get_H2(Psi, y, N, Lbare)
    if H2 < 1e-40:
        return 0.0
    return -0.5*(H2p-H2m)/(2.0*eps*H2)

def ode(N, state, Lbare, alpha_tilde):
    """
    α̃ = α_GB·H₀² (dimensionless).
    GB source in N-deriv: α̃·24·(H²/H₀²)·(1−ε_H) = α̃·24·E²·(1−ε_H)
    where H₀ = h_target·H100 = 67.4 km/s/Mpc.
    In H100 units: H₀² = h_target².
    """
    Psi, y  = state
    H2      = get_H2(Psi, y, N, Lbare)
    eps_H   = get_eps_H(Psi, y, N, Lbare)

    # R/H² source (Phase 1 backbone)
    R_H2    = 6.0*(2.0-eps_H)
    src_R   = F_eff_prime(Psi) * R_H2   # = 2Λ₀Ψ·R/H²

    # GB source: α̃·(𝒢/H²) = α̃·24H²(1−ε_H)
    # H² is in H100² units; α̃ is dimensionless (α_GB·H₀² = α_GB·h_target²·H100²)
    # But in our ODE everything is in H100 units, so use α_GB·H100² = α̃/h_target²
    # Simpler: define α̃ as the coefficient of 24H²(1−ε_H)/H_today²
    #   → src_GB = α̃ * 24 * (H2/h_target**2) * (1−ε_H)
    E2      = H2 / h_target**2          # (H/H₀)²
    src_GB  = alpha_tilde * 24.0 * E2 * (1.0 - eps_H)

    dy_dN   = src_R + src_GB - (3.0-eps_H)*y
    return [y, dy_dN]

def integrate(psi_ini, Lbare, alpha_tilde, z_ini=1e5, n_pts=4000):
    N_ini  = np.log(1.0/(1.0+z_ini))
    N_eval = np.linspace(N_ini, 0.0, n_pts)
    sol = solve_ivp(ode, (N_ini, 0.0), [psi_ini, 0.0],
                    args=(Lbare, alpha_tilde), t_eval=N_eval,
                    method='RK45', rtol=1e-8, atol=1e-12, max_step=0.03)
    if not sol.success:
        return None
    N_arr  = sol.t; z_arr = np.exp(-N_arr)-1.0
    ps_arr = sol.y[0]; y_arr = sol.y[1]
    H_arr, Fv_arr = [], []
    for i in range(len(N_arr)):
        H2 = get_H2(ps_arr[i], y_arr[i], N_arr[i], Lbare)
        H_arr.append(H100*np.sqrt(max(H2,0.0)))
        Fv_arr.append(F_eff(ps_arr[i]))
    H_arr = np.array(H_arr); Fv_arr = np.array(Fv_arr)
    idx = np.argsort(z_arr)
    return dict(z=z_arr[idx], psi=ps_arr[idx], y=y_arr[idx],
                H=H_arr[idx], Feff=Fv_arr[idx],
                psi0=float(ps_arr[-1]), y0=float(y_arr[-1]),
                Feff0=float(Fv_arr[-1]), H0=float(H_arr[-1]),
                Lbare=Lbare, alpha=alpha_tilde)

def calibrate(psi0, y0):
    Fv   = F_eff(psi0)
    Fp   = F_eff_prime(psi0)
    coef = 3.0*Fv + 6.0*Fp*y0 - 0.5*y0**2
    return h_target**2*coef - omh2_m - omh2_r

def run(alpha_tilde, n_iter=6):
    Lbare = Lbare_phase1
    bg = None
    for _ in range(n_iter):
        bg = integrate(0.0, Lbare, alpha_tilde)
        if bg is None:
            return None
        Lbare_new = calibrate(bg['psi0'], bg['y0'])
        if abs(Lbare_new-Lbare) < 1e-10:
            break
        Lbare = Lbare_new
    return integrate(0.0, Lbare, alpha_tilde)

def H_interp(z, bg):
    return float(np.interp(z, bg['z'], bg['H']))

def theta_star(bg):
    H0 = bg['H0']; h = H0/100.0
    Ogam = 2.469e-5/h**2
    def rs_int(z):
        R  = (3.0*omh2_b/h**2)/(4.0*Ogam*(1+z))
        cs = c_kms/np.sqrt(3.0*(1.0+R))
        Hz = H_interp(z, bg)
        return cs/Hz if Hz>0 else 0.0
    rs, _ = quad(rs_int, z_drag, 1e4, limit=200, epsrel=1e-5)
    def DC_int(z):
        Hz = H_interp(z, bg)
        return c_kms/Hz if Hz>0 else 0.0
    DC, _ = quad(DC_int, 0, z_star, limit=200, epsrel=1e-5)
    return 100.0*rs/DC if DC>0 else np.nan

def chi2_DESI(bg):
    return float(np.sum(((np.array([H_interp(z,bg) for z in DESI_z])-DESI_H)/DESI_s)**2))

# ─── Main ────────────────────────────────────────────────────────────────────
print("="*72)
print("SIM133 — Phase 3: Gauss-Bonnet Sourcing (Direction B)")
print("="*72)
print(f"  Action: ½(1+2Λ₀Ψ²)R - ½(∂Ψ)² + α̃·Ψ·𝒢  (Λ₀=0.003, Phase 1 backbone)")
print(f"  𝒢 changes sign at matter-Λ transition → non-monotonic Ψ(z)")
print(f"  Baseline: SIM121C {SIM121C_tension:.2f}σ")
print()

alpha_grid = [1e-8, 1e-7, 1e-6, 1e-5, 1e-4]

results_scan = []

print(f"  {'α̃':>10}  {'Ψ(z*)':>9}  {'Ψ₀':>8}  {'F₀':>8}  {'100θ*':>8}  {'tension':>8}")

for at in alpha_grid:
    bg = run(at)
    if bg is None:
        print(f"  {at:10.2e}  FAILED")
        results_scan.append(None)
        continue

    # Ψ at matter-Λ equality z ≈ 0.3
    psi_zeq = float(np.interp(0.3, bg['z'], bg['psi']))
    th      = theta_star(bg)
    c2      = chi2_DESI(bg)
    ten     = np.sqrt(c2/len(DESI_z))
    c2t     = ((th-theta_obs)/theta_obs_err)**2 if not np.isnan(th) else np.nan

    print(f"  {at:10.2e}  {psi_zeq:9.5f}  {bg['psi0']:8.5f}  "
          f"{bg['Feff0']:8.5f}  {th:8.5f}  {ten:8.3f}σ")

    results_scan.append(dict(alpha=at, bg=bg, psi_zeq=psi_zeq,
                             theta=th, chi2=c2, tension=ten, chi2_theta=c2t))

valid = [r for r in results_scan if r is not None and not np.isnan(r.get('theta',np.nan))]
if valid:
    best = min(valid, key=lambda r: r['chi2'] + r['chi2_theta'])
    print(f"\n  Best: α̃={best['alpha']:.2e}, tension={best['tension']:.3f}σ, "
          f"100θ*={best['theta']:.5f}, Δχ²={best['chi2']-SIM121C_chi2:+.3f}")
    bg_best = best['bg']

    print(f"\n  H(z) at best:")
    print(f"    {'z':>6}  {'H_obs':>7}  {'H_mod':>7}  {'pull':>6}")
    for z_d,Ho,sd in zip(DESI_z,DESI_H,DESI_s):
        Hm = H_interp(z_d, bg_best)
        print(f"    {z_d:6.3f}  {Ho:7.1f}  {Hm:7.2f}  {(Hm-Ho)/sd:6.2f}")
else:
    best = None; bg_best = None

if best:
    pass_desi  = best['tension'] < 2.0
    pass_theta = best['chi2_theta'] < 4.0
    beats_ref  = best['chi2'] < SIM121C_chi2
    if pass_desi and pass_theta:   verdict = "PASS"
    elif beats_ref and pass_theta: verdict = "PARTIAL"
    elif beats_ref:                verdict = "PARTIAL"
    else:                          verdict = "FAIL"
else:
    verdict = "FAIL"; pass_desi=pass_theta=beats_ref=False

print()
print("="*72)
print(f"SIM133 RESULT:  {verdict}")
if best:
    print(f"  Best α̃={best['alpha']:.2e}:  tension={best['tension']:.3f}σ  "
          f"(Phase1 floor: {SIM121C_tension:.2f}σ),  Δχ²={best['chi2']-SIM121C_chi2:+.3f}")
print("="*72)

# ─── Figures ─────────────────────────────────────────────────────────────────
# Fig 1: Ψ(z) for all α̃
if len(valid) > 1:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(valid)))

    for i, r in enumerate(valid):
        z_p  = r['bg']['z']; mask = z_p < 5
        axes[0].plot(z_p[mask], r['bg']['psi'][mask], lw=2, color=colors[i],
                     label=f"α̃={r['alpha']:.0e}")
    axes[0].axhline(0, color='gray', ls='--', lw=1)
    axes[0].axvline(0.3, color='orange', ls=':', lw=1.5, label='z≈0.3 (Λ-eq)')
    axes[0].set_xlabel(r'Redshift $z$'); axes[0].set_ylabel(r'$\Psi(z)$ [$M_{\rm Pl}$]')
    axes[0].set_title(r'SIM133: GB-sourced $\Psi(z)$ — sign flip at $z\approx0.3$')
    axes[0].legend(fontsize=8); axes[0].grid(alpha=0.3)

    for i, r in enumerate(valid):
        axes[1].semilogy([r['alpha']], [r['tension']], 'o', ms=8, color=colors[i])
    axes[1].axhline(SIM121C_tension, color='#2166ac', ls='--', lw=1.5,
                    label=f'Phase1 {SIM121C_tension:.2f}σ')
    axes[1].axhline(2.0, color='green', ls=':', lw=1.5, label='2σ PASS')
    alpha_arr = [r['alpha'] for r in valid]; tens_arr = [r['tension'] for r in valid]
    axes[1].semilogx(alpha_arr, tens_arr, '-', color='#d73027', lw=2)
    axes[1].set_xlabel(r'$\tilde{\alpha}$ (GB coupling)'); axes[1].set_ylabel('DESI tension [σ]')
    axes[1].set_title('SIM133: DESI tension vs GB coupling')
    axes[1].legend(fontsize=9); axes[1].grid(alpha=0.3)

    plt.tight_layout()
    for ext in ('pdf','png'):
        fig.savefig(os.path.join(OUT,f'sim133_psi_tension.{ext}'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("  Saved sim133_psi_tension")

# Fig 2: H(z) for best
if bg_best is not None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    z_fine = np.linspace(0.01, 2.5, 300)
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(valid)))

    for i, r in enumerate(valid):
        H_fine = [H_interp(z, r['bg']) for z in z_fine]
        axes[0].plot(z_fine, H_fine, lw=2, color=colors[i],
                     label=f"α̃={r['alpha']:.0e} ({r['tension']:.2f}σ)")
    axes[0].errorbar(DESI_z, DESI_H, yerr=DESI_s, fmt='ko', ms=7, capsize=4,
                     label='DESI', zorder=5)
    axes[0].set_xlabel(r'Redshift $z$'); axes[0].set_ylabel(r'$H(z)$ [km/s/Mpc]')
    axes[0].set_title('SIM133: H(z) for all GB couplings')
    axes[0].legend(fontsize=8); axes[0].grid(alpha=0.3)

    # Ψ at z<5 for best
    z_p = bg_best['z']; mask = z_p < 5
    axes[1].plot(z_p[mask], bg_best['psi'][mask], color='#d73027', lw=2,
                 label=f"α̃={best['alpha']:.1e} (best)")
    axes[1].axhline(0, color='gray', ls='--', lw=1)
    axes[1].axvline(0.3, color='orange', ls=':', lw=1.5, label='z≈0.3 sign flip')
    axes[1].set_xlabel(r'Redshift $z$'); axes[1].set_ylabel(r'$\Psi(z)$')
    axes[1].set_title('GB-sourced Ψ: matter→Λ sign flip')
    axes[1].legend(fontsize=9); axes[1].grid(alpha=0.3)

    plt.tight_layout()
    for ext in ('pdf','png'):
        fig.savefig(os.path.join(OUT,f'sim133_Hz_best.{ext}'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("  Saved sim133_Hz_best")

# ─── Save JSON ───────────────────────────────────────────────────────────────
out_json = {
    "sim": "SIM133", "phase": "Phase 3", "verdict": verdict,
    "description": "Gauss-Bonnet sourcing: α̃·Ψ·𝒢 added, Phase 1 backbone unchanged",
    "action": "S=∫d⁴x√g[½(1+2Λ₀Ψ²)R - ½(∂Ψ)² + α̃·Ψ·𝒢] + S_SM",
    "Lambda0": Lambda0,
    "scan": [{"alpha": r['alpha'], "tension": r['tension'], "chi2_DESI": r['chi2'],
               "theta_100": r['theta'], "psi0": r['bg']['psi0'],
               "Feff0": r['bg']['Feff0'], "psi_z03": r['psi_zeq']}
              for r in valid],
    "best": {"alpha": best['alpha'], "tension": best['tension'],
             "chi2_DESI": best['chi2'], "delta_chi2": best['chi2']-SIM121C_chi2,
             "theta_100": best['theta']} if best else None,
    "baseline": {"tension": SIM121C_tension, "chi2": SIM121C_chi2},
    "pass_desi": bool(pass_desi), "pass_theta": bool(pass_theta), "beats_ref": bool(beats_ref),
}
with open(os.path.join(OUT,'sim133_results.json'),'w') as f:
    json.dump(out_json, f, indent=2)
print("  Saved sim133_results.json")
print(f"\nAll outputs: {OUT}")
print("SIM133 complete.")
