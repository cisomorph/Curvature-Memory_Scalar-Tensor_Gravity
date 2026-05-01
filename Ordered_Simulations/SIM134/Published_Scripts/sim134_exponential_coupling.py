"""
SIM134 — CMSTG Phase 3: Exponential Coupling (Dilaton-type)
===========================================================
Replace the quadratic backbone with an exponential:
  F_eff = ½·exp(−2Λ₀Ψ²) + ξΨ

Properties:
  - No gravity flip possible: exp term always positive, bounded in (0, ½]
  - Stable attractor at Ψ* where dF_eff/dΨ = 0:
      −2Λ₀Ψ*·exp(−2Λ₀Ψ*²) + ξ = 0
  - F_eff(Ψ*) = (½ + 2Λ₀Ψ*²)·exp(−2Λ₀Ψ*²) < Phase1 F₀=0.521 when Λ₀ ≳ 0.30
  - CMB naturally protected: Ψ≈0 at z=1090 (hasn't grown yet) → F_eff≈½ (GR)
  - First-principles: dilaton coupling exp(−Λ₀Ψ²)·R appears in string eff. actions

Action:
  S = ∫d⁴x√(-g) [ F_eff(Ψ)·R − ½(∂Ψ)² ] + S_SM
  F_eff(Ψ) = ½·exp(−2Λ₀Ψ²) + ξΨ

Field EOM (FRW, N = ln a):
  Ψ'' + (3−ε_H)Ψ' = F_eff′(Ψ)·R/H² = [−2Λ₀Ψ·exp(−2Λ₀Ψ²) + ξ]·6(2−ε_H)

  At Ψ=0: source = ξ·R/H² > 0 → Ψ grows from zero (no tuning needed)
  As Ψ → Ψ*: source → 0 → field locks at attractor

Modified Friedmann:
  H²[3F_eff + 6·F_eff′(Ψ)·Ψ' − ½Ψ'²] = ω_m a⁻³ + ω_r a⁻⁴ + Λ_bare

Attractor properties (analytic):
  Let y* = 2Λ₀Ψ*²: F_eff(Ψ*) = (½+y*)·exp(−y*)
  y* → 0 as Λ₀ → ∞ → F_eff → ½ (GR)
  F_eff(Ψ*) < 0.521 when Λ₀ ≳ 0.30 (estimated from analytic formula)

Scan: Λ₀ ∈ [0.20, 0.33, 0.50, 1.00], ξ = 1/6 (fixed), Ψ_ini = 0
Baseline: SIM121C DESI tension = 2.77σ
"""

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq
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

XI             = 1.0/6.0
F0_phase1      = 0.5 + 0.003*2.62**2
Lbare_phase1   = 3.0*F0_phase1*h_target**2 - omh2_m - omh2_r

# ─── Physics functions ───────────────────────────────────────────────────────

def F_eff(Psi, L0):
    return 0.5*np.exp(-2.0*L0*Psi**2) + XI*Psi

def F_eff_prime(Psi, L0):
    return -2.0*L0*Psi*np.exp(-2.0*L0*Psi**2) + XI

def Psi_attractor(L0):
    """Solve F_eff′(Ψ*) = 0: −2Λ₀Ψ*·exp(−2Λ₀Ψ*²) + ξ = 0"""
    def eq(x):
        return -2.0*L0*x*np.exp(-2.0*L0*x**2) + XI
    # Attractor exists for x > 0 where exp term can balance ξ
    # For large L0, attractor is near 0; for small L0 it's larger
    try:
        return brentq(eq, 1e-6, 10.0, xtol=1e-8)
    except ValueError:
        return np.nan

def F_eff_at_attractor(L0):
    Ps = Psi_attractor(L0)
    if np.isnan(Ps):
        return np.nan
    return F_eff(Ps, L0)

def get_H2(Psi, y, N, Lbare, L0):
    a    = np.exp(N)
    rhs  = omh2_m/a**3 + omh2_r/a**4 + Lbare
    Fv   = F_eff(Psi, L0)
    Fp   = F_eff_prime(Psi, L0)
    coef = 3.0*Fv + 6.0*Fp*y - 0.5*y**2
    if coef < 0.1:
        coef = max(3.0*Fv, 0.05)
    return rhs/coef

def get_eps_H(Psi, y, N, Lbare, L0):
    eps = 5e-4
    H2p = get_H2(Psi, y, N+eps, Lbare, L0)
    H2m = get_H2(Psi, y, N-eps, Lbare, L0)
    H2  = get_H2(Psi, y, N, Lbare, L0)
    if H2 < 1e-40:
        return 0.0
    return -0.5*(H2p-H2m)/(2.0*eps*H2)

def ode(N, state, Lbare, L0):
    Psi, y = state
    H2     = get_H2(Psi, y, N, Lbare, L0)
    eps_H  = get_eps_H(Psi, y, N, Lbare, L0)
    source = F_eff_prime(Psi, L0) * 6.0*(2.0-eps_H)
    dy_dN  = source - (3.0-eps_H)*y
    return [y, dy_dN]

def integrate(psi_ini, Lbare, L0, z_ini=1e5, n_pts=4000):
    N_ini  = np.log(1.0/(1.0+z_ini))
    N_eval = np.linspace(N_ini, 0.0, n_pts)
    sol = solve_ivp(ode, (N_ini, 0.0), [psi_ini, 0.0],
                    args=(Lbare, L0), t_eval=N_eval,
                    method='RK45', rtol=1e-8, atol=1e-12, max_step=0.05)
    if not sol.success:
        return None
    N_arr  = sol.t; z_arr = np.exp(-N_arr)-1.0
    ps_arr = sol.y[0]; y_arr = sol.y[1]
    H_arr, Fv_arr = [], []
    for i in range(len(N_arr)):
        H2 = get_H2(ps_arr[i], y_arr[i], N_arr[i], Lbare, L0)
        H_arr.append(H100*np.sqrt(max(H2,0.0)))
        Fv_arr.append(F_eff(ps_arr[i], L0))
    H_arr = np.array(H_arr); Fv_arr = np.array(Fv_arr)
    idx = np.argsort(z_arr)
    return dict(z=z_arr[idx], psi=ps_arr[idx], y=y_arr[idx],
                H=H_arr[idx], Feff=Fv_arr[idx],
                psi0=float(ps_arr[-1]), y0=float(y_arr[-1]),
                Feff0=float(Fv_arr[-1]), H0=float(H_arr[-1]),
                Lbare=Lbare, L0=L0)

def calibrate(psi0, y0, L0):
    Fv   = F_eff(psi0, L0)
    Fp   = F_eff_prime(psi0, L0)
    coef = 3.0*Fv + 6.0*Fp*y0 - 0.5*y0**2
    return h_target**2*coef - omh2_m - omh2_r

def run(L0, n_iter=6):
    Lbare = Lbare_phase1
    bg = None
    for _ in range(n_iter):
        bg = integrate(0.0, Lbare, L0)
        if bg is None:
            return None
        Lbare_new = calibrate(bg['psi0'], bg['y0'], L0)
        if abs(Lbare_new-Lbare) < 1e-10:
            break
        Lbare = Lbare_new
    return integrate(0.0, Lbare, L0)

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

# ─── Analytic attractor survey ────────────────────────────────────────────────
print("="*72)
print("SIM134 — Phase 3: Exponential Coupling F_eff = ½exp(−2Λ₀Ψ²) + ξΨ")
print("="*72)
print(f"  ξ = 1/6 (fixed), Ψ_ini = 0")
print(f"  No gravity flip possible; attractor at Ψ* where F_eff′=0")
print()
print("  Analytic attractor preview:")
print(f"  {'Λ₀':>6}  {'Ψ*':>7}  {'F_eff(Ψ*)':>10}  {'vs Phase1':>10}")
for L0 in [0.10, 0.20, 0.33, 0.50, 1.00]:
    Ps  = Psi_attractor(L0)
    Fv  = F_eff_at_attractor(L0) if not np.isnan(Ps) else np.nan
    cmp = f"{'better' if not np.isnan(Fv) and Fv < F0_phase1 else 'worse':>10}"
    print(f"  {L0:6.3f}  {Ps:7.4f}  {Fv:10.6f}  {cmp}")
print()

# ─── Scan ────────────────────────────────────────────────────────────────────
Lambda0_grid = [0.20, 0.33, 0.50, 1.00]

results_scan = []

print(f"  {'Λ₀':>6}  {'Ψ*':>7}  {'F*(th)':>8}  {'Ψ₀':>8}  {'F₀':>8}  "
      f"{'Feff(z*)':>9}  {'100θ*':>8}  {'tension':>8}")

for L0 in Lambda0_grid:
    Ps_th  = Psi_attractor(L0)
    Fv_th  = F_eff_at_attractor(L0)

    bg = run(L0)
    if bg is None:
        print(f"  {L0:6.3f}  {Ps_th:7.4f}  {Fv_th:8.5f}  FAILED")
        results_scan.append(None)
        continue

    # F_eff at CMB epoch
    Feff_CMB = float(np.interp(z_star, bg['z'], bg['Feff']))
    th  = theta_star(bg)
    c2  = chi2_DESI(bg)
    ten = np.sqrt(c2/len(DESI_z))
    c2t = ((th-theta_obs)/theta_obs_err)**2 if not np.isnan(th) else np.nan

    print(f"  {L0:6.3f}  {Ps_th:7.4f}  {Fv_th:8.5f}  {bg['psi0']:8.5f}  "
          f"{bg['Feff0']:8.5f}  {Feff_CMB:9.6f}  {th:8.5f}  {ten:8.3f}σ")

    results_scan.append(dict(L0=L0, bg=bg, Ps_th=Ps_th, Fv_th=Fv_th,
                             Feff_CMB=Feff_CMB, theta=th, chi2=c2,
                             tension=ten, chi2_theta=c2t))

valid = [r for r in results_scan if r is not None and not np.isnan(r.get('theta',np.nan))]
if valid:
    best = min(valid, key=lambda r: r['chi2'] + r['chi2_theta'])
    print(f"\n  Best: Λ₀={best['L0']:.3f}, tension={best['tension']:.3f}σ, "
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
print(f"SIM134 RESULT:  {verdict}")
if best:
    print(f"  Best Λ₀={best['L0']:.3f}:  tension={best['tension']:.3f}σ  "
          f"(Phase1 floor: {SIM121C_tension:.2f}σ),  Δχ²={best['chi2']-SIM121C_chi2:+.3f}")
print("="*72)

# ─── Figures ─────────────────────────────────────────────────────────────────
# Fig 1: F_eff shape for each Λ₀
Psi_range = np.linspace(-0.5, 2.0, 300)
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
colors = plt.cm.plasma(np.linspace(0.1, 0.85, len(Lambda0_grid)))

for i, L0 in enumerate(Lambda0_grid):
    Fv_range = [F_eff(p, L0) for p in Psi_range]
    Ps  = Psi_attractor(L0)
    axes[0].plot(Psi_range, Fv_range, lw=2, color=colors[i], label=f'Λ₀={L0:.2f}')
    if not np.isnan(Ps):
        axes[0].axvline(Ps, color=colors[i], ls=':', lw=1.2)
axes[0].axhline(F0_phase1, color='#2166ac', ls='--', lw=1.5, label=f'Phase1 F₀={F0_phase1:.4f}')
axes[0].axhline(0.5, color='gray', ls=':', lw=1.5, label='GR')
axes[0].set_xlabel(r'$\Psi$ [$M_{\rm Pl}$]')
axes[0].set_ylabel(r'$F_{\rm eff}(\Psi) = \frac{1}{2}e^{-2\Lambda_0\Psi^2} + \xi\Psi$')
axes[0].set_title('SIM134: F_eff shape and attractor (dashed)')
axes[0].legend(fontsize=9); axes[0].grid(alpha=0.3)
axes[0].set_ylim(0.3, 0.8)

# Tension vs Λ₀
if len(valid) > 1:
    L0s  = [r['L0'] for r in valid]
    tens = [r['tension'] for r in valid]
    Fvs  = [r['bg']['Feff0'] for r in valid]

    ax2 = axes[1]
    ax2.plot(L0s, tens, 'o-', color='#d73027', lw=2, ms=7, label='DESI tension')
    ax2.axhline(SIM121C_tension, color='#2166ac', ls='--', lw=1.5,
                label=f'Phase1 {SIM121C_tension:.2f}σ')
    ax2.axhline(2.0, color='green', ls=':', lw=1.5, label='2σ PASS')
    ax2.set_xlabel(r'$\Lambda_0$'); ax2.set_ylabel('DESI tension [σ]')
    ax2.set_title(r'SIM134: DESI tension vs $\Lambda_0$')
    ax2.legend(fontsize=9); ax2.grid(alpha=0.3)

plt.tight_layout()
for ext in ('pdf','png'):
    fig.savefig(os.path.join(OUT,f'sim134_Feff_tension.{ext}'), dpi=150, bbox_inches='tight')
plt.close(fig)
print("  Saved sim134_Feff_tension")

# Fig 2: Ψ(z) and H(z)
if bg_best is not None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    for i, r in enumerate(valid):
        z_p = r['bg']['z']; mask = z_p < 3000
        axes[0].semilogx(z_p[mask]+1, r['bg']['psi'][mask], lw=2,
                         color=colors[i], label=f"Λ₀={r['L0']:.2f}")
    axes[0].axhline(0, color='gray', ls='--', lw=1)
    axes[0].axvline(z_star+1, color='orange', ls=':', lw=1.5, label=f'z*')
    axes[0].set_xlabel(r'$1+z$'); axes[0].set_ylabel(r'$\Psi(z)$ [$M_{\rm Pl}$]')
    axes[0].set_title('SIM134: Ψ growing to stable attractor')
    axes[0].legend(fontsize=9); axes[0].grid(alpha=0.3)

    z_fine = np.linspace(0.01, 2.5, 300)
    for i, r in enumerate(valid):
        H_fine = [H_interp(z, r['bg']) for z in z_fine]
        axes[1].plot(z_fine, H_fine, lw=2, color=colors[i],
                     label=f"Λ₀={r['L0']:.2f} ({r['tension']:.2f}σ)")
    axes[1].errorbar(DESI_z, DESI_H, yerr=DESI_s, fmt='ko', ms=7, capsize=4,
                     label='DESI', zorder=5)
    axes[1].set_xlabel(r'Redshift $z$'); axes[1].set_ylabel(r'$H(z)$ [km/s/Mpc]')
    axes[1].set_title('SIM134: H(z) vs DESI')
    axes[1].legend(fontsize=9); axes[1].grid(alpha=0.3)

    plt.tight_layout()
    for ext in ('pdf','png'):
        fig.savefig(os.path.join(OUT,f'sim134_psi_Hz.{ext}'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("  Saved sim134_psi_Hz")

# ─── Save JSON ───────────────────────────────────────────────────────────────
out_json = {
    "sim": "SIM134", "phase": "Phase 3", "verdict": verdict,
    "description": "Exponential coupling: F_eff=½exp(-2Λ₀Ψ²)+ξΨ, dilaton-type, no gravity flip",
    "action": "S=∫d⁴x√g[F_eff(Ψ)R - ½(∂Ψ)²]+S_SM, F_eff=½exp(-2Λ₀Ψ²)+ξΨ",
    "xi": XI,
    "scan": [{"L0": r['L0'], "Psi_attractor": r['Ps_th'], "Fatt_theory": r['Fv_th'],
               "tension": r['tension'], "chi2_DESI": r['chi2'],
               "theta_100": r['theta'], "psi0": r['bg']['psi0'],
               "Feff0": r['bg']['Feff0'], "Feff_CMB": r['Feff_CMB']}
              for r in valid],
    "best": {"L0": best['L0'], "tension": best['tension'],
             "chi2_DESI": best['chi2'], "delta_chi2": best['chi2']-SIM121C_chi2,
             "theta_100": best['theta']} if best else None,
    "baseline": {"tension": SIM121C_tension, "chi2": SIM121C_chi2},
    "pass_desi": bool(pass_desi), "pass_theta": bool(pass_theta), "beats_ref": bool(beats_ref),
}
with open(os.path.join(OUT,'sim134_results.json'),'w') as f:
    json.dump(out_json, f, indent=2)
print("  Saved sim134_results.json")
print(f"\nAll outputs: {OUT}")
print("SIM134 complete.")
