"""
SIM135 — CMSTG Phase 3: Bi-scalar (Ψ frozen at Phase 1 canonical, φ sourced by curvature)
===========================================================================================
SIM131–134 all FAIL because Ψ_ini=0 gives F_eff≈½ at z_CMB; Phase 1 canonical
requires F_eff=0.521. Curvature sourcing drives Ψ away from zero before z_CMB, so
CMB acoustic scale (100θ*≈0.92 instead of 1.041) is broken in every variant.

Bi-scalar fix:
  - LOCK Ψ at Phase 1 canonical (Ψ̄=2.62 M_Pl, F₀=0.5206).
    F₀ is a constant background contribution to G_eff.
  - Add a SECOND scalar φ driven by curvature for late-time DESI dynamics.
    φ_ini=0  →  F_total = F₀ = 0.521 at z_CMB  (CMB preserved by construction).
    φ grows at late times  →  F_total changes  →  H(z) shifts.

Action:
  S = ∫d⁴x√g [ (F₀ + ξ_φ·φ)·R  −  ½(∂φ)² ] + S_SM
  F₀ = 0.52059 (Phase 1 canonical), φ_ini = 0

Field EOM (FRW, N = ln a):
  φ'' + (3−ε_H)φ' = ξ_φ · R/H² = ξ_φ · 6(2−ε_H)

Modified Friedmann:
  H² [ 3(F₀+ξ_φφ) + 6ξ_φφ' − ½φ'² ] = ω_m a⁻³ + ω_r a⁻⁴ + Λ_bare

ε_H self-consistent formula (derived analytically from d(H²)/dN):
  With A = 12ξ−3y, B = y−6ξ (where y=φ'):
  ε_H = [ −Q̇/(2Q) + (3ξy − A·B)/(2D) ] / [ 1 + B²/(2D) ]
  where Q=ρ_tot, D=3F+6ξy−½y², Q̇=dQ/dN

Scan: ξ_φ ∈ {0.01, 0.03, 0.05, 0.08, 0.10, 0.15}
Baseline: SIM121C DESI tension = 2.77σ, chi²_DESI = 18.26
PASS: tension < 2.0σ  AND  |100θ*−1.04101| < 2×0.00029
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

# ─── Constants ────────────────────────────────────────────────────────────────
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

# Phase 1 canonical frozen contribution
F0 = 0.5 + 0.003 * 2.62**2        # = 0.520593
Lbare_phase1 = 3.0 * F0 * h_target**2 - omh2_m - omh2_r

# ─── Physics functions ────────────────────────────────────────────────────────

def F_eff(phi, xi):
    return F0 + xi * phi

def F_eff_prime(xi):
    return xi   # dF/dphi = xi (constant)

def get_eps_H_analytical(phi, y, N, Lbare, xi):
    """
    Self-consistent ε_H from d(H²)/dN.
    Eliminates circular dependency on y' by solving algebraically.
    """
    a   = np.exp(N)
    Q   = omh2_m/a**3 + omh2_r/a**4 + Lbare
    dQ  = -3.0*omh2_m/a**3 - 4.0*omh2_r/a**4
    F   = F0 + xi*phi
    D   = 3.0*F + 6.0*xi*y - 0.5*y**2
    if abs(D) < 1e-10:
        return 1.5   # fallback: matter domination
    # Coefficients in y' = A + eps_H * B
    A   = 12.0*xi - 3.0*y
    B   = y - 6.0*xi
    # Self-consistent solution
    num = -dQ/(2.0*Q) + (3.0*xi*y - A*B)/(2.0*D)
    den = 1.0 + B**2/(2.0*D)
    if abs(den) < 1e-10:
        return 1.5
    return num / den

def get_H2(phi, y, N, Lbare, xi):
    a   = np.exp(N)
    rho = omh2_m/a**3 + omh2_r/a**4 + Lbare
    F   = F0 + xi*phi
    D   = 3.0*F + 6.0*xi*y - 0.5*y**2
    if D < 0.1:
        D = max(3.0*F, 0.01)
    return rho / D

def ode(N, state, Lbare, xi):
    phi, y = state
    eps_H  = get_eps_H_analytical(phi, y, N, Lbare, xi)
    eps_H  = np.clip(eps_H, -1.0, 4.0)
    source = xi * 6.0*(2.0 - eps_H)
    dy_dN  = source - (3.0 - eps_H)*y
    return [y, dy_dN]

def integrate(Lbare, xi, z_ini=1e5, n_pts=4000):
    N_ini  = np.log(1.0/(1.0+z_ini))
    N_eval = np.linspace(N_ini, 0.0, n_pts)
    sol = solve_ivp(ode, (N_ini, 0.0), [0.0, 0.0],
                    args=(Lbare, xi), t_eval=N_eval,
                    method='RK45', rtol=1e-8, atol=1e-12, max_step=0.05)
    if not sol.success:
        return None
    N_arr  = sol.t
    z_arr  = np.exp(-N_arr) - 1.0
    ph_arr = sol.y[0]
    y_arr  = sol.y[1]
    H_arr  = np.array([H100*np.sqrt(max(get_H2(ph_arr[i],y_arr[i],N_arr[i],Lbare,xi),0.0))
                       for i in range(len(N_arr))])
    Fv_arr = np.array([F_eff(ph_arr[i],xi) for i in range(len(N_arr))])
    idx    = np.argsort(z_arr)
    return dict(z=z_arr[idx], N=N_arr[idx], phi=ph_arr[idx], y=y_arr[idx],
                H=H_arr[idx], Feff=Fv_arr[idx],
                phi0=float(ph_arr[-1]), y0=float(y_arr[-1]),
                Feff0=float(Fv_arr[-1]), H0=float(H_arr[-1]),
                Lbare=Lbare, xi=xi)

def calibrate(phi0, y0, xi):
    F   = F0 + xi*phi0
    D   = 3.0*F + 6.0*xi*y0 - 0.5*y0**2
    return h_target**2 * D - omh2_m - omh2_r

def run(xi, n_iter=8):
    Lbare = Lbare_phase1
    bg    = None
    for _ in range(n_iter):
        bg = integrate(Lbare, xi)
        if bg is None:
            return None
        Lbare_new = calibrate(bg['phi0'], bg['y0'], xi)
        if abs(Lbare_new - Lbare) < 1e-10:
            break
        Lbare = Lbare_new
    return integrate(Lbare, xi)

def H_interp(z, bg):
    return float(np.interp(z, bg['z'], bg['H']))

def theta_star(bg):
    H0 = bg['H0']; h = H0/100.0
    Ogam = 2.469e-5/h**2
    def rs_int(z):
        R  = (3.0*omh2_b/h**2) / (4.0*Ogam*(1+z))
        cs = c_kms / np.sqrt(3.0*(1.0+R))
        Hz = H_interp(z, bg)
        return cs/Hz if Hz > 0 else 0.0
    rs, _ = quad(rs_int, z_drag, 1e4, limit=200, epsrel=1e-5)
    def DC_int(z):
        Hz = H_interp(z, bg)
        return c_kms/Hz if Hz > 0 else 0.0
    DC, _ = quad(DC_int, 0, z_star, limit=200, epsrel=1e-5)
    return 100.0*rs/DC if DC > 0 else np.nan

def chi2_DESI(bg):
    return float(np.sum(((np.array([H_interp(z,bg) for z in DESI_z])-DESI_H)/DESI_s)**2))

# ─── Main scan ────────────────────────────────────────────────────────────────
print("="*72)
print("SIM135 — Phase 3: Bi-scalar (Ψ frozen, φ curvature-sourced)")
print("="*72)
print(f"  Action: (F₀ + ξ_φ·φ)·R − ½(∂φ)²,  F₀={F0:.5f},  φ_ini=0")
print(f"  EOM:  φ''+( 3−ε_H)φ' = ξ_φ·6(2−ε_H)")
print(f"  Scan: ξ_φ ∈ [0.01, 0.03, 0.05, 0.08, 0.10, 0.15]")
print(f"  Baseline: SIM121C {SIM121C_tension:.2f}σ, χ²={SIM121C_chi2:.2f}")
print()

xi_grid = [0.01, 0.03, 0.05, 0.08, 0.10, 0.15]

results = []

print(f"  {'ξ_φ':>6}  {'φ₀':>8}  {'F_eff0':>8}  {'F_CMB':>8}  {'100θ*':>8}  {'tension':>8}")
print(f"  {'-'*62}")

for xi in xi_grid:
    bg = run(xi)
    if bg is None:
        print(f"  {xi:6.3f}  FAILED")
        results.append(None)
        continue

    # Get F_eff at z_CMB
    z_arr = bg['z']
    F_arr = bg['Feff']
    F_cmb = float(np.interp(z_star, z_arr, F_arr))

    th  = theta_star(bg)
    c2  = chi2_DESI(bg)
    ten = np.sqrt(c2/len(DESI_z))
    c2t = ((th-theta_obs)/theta_obs_err)**2 if not np.isnan(th) else np.nan

    print(f"  {xi:6.3f}  {bg['phi0']:8.4f}  {bg['Feff0']:8.5f}  "
          f"{F_cmb:8.5f}  {th:8.5f}  {ten:8.3f}σ")

    results.append(dict(xi=xi, bg=bg, theta=th, chi2=c2, tension=ten,
                        chi2_theta=c2t, F_cmb=F_cmb))

# Best by joint DESI + CMB chi²
valid = [r for r in results if r is not None and not np.isnan(r['theta'])]

if valid:
    best = min(valid, key=lambda r: r['chi2'] + r['chi2_theta'])
    print(f"\n  Best: ξ_φ={best['xi']:.3f},  tension={best['tension']:.3f}σ  "
          f"(Phase1 floor: {SIM121C_tension:.2f}σ),  Δχ²={best['chi2']-SIM121C_chi2:+.3f}")
    print(f"        100θ*={best['theta']:.5f},  F_eff(z_CMB)={best['F_cmb']:.5f}")
    bg_best = best['bg']

    print(f"\n  H(z) pulls at best ξ_φ={best['xi']:.3f}:")
    print(f"    {'z':>6}  {'H_obs':>7}  {'H_mod':>7}  {'pull':>6}")
    for z_d, Ho, sd in zip(DESI_z, DESI_H, DESI_s):
        Hm = H_interp(z_d, bg_best)
        print(f"    {z_d:6.3f}  {Ho:7.1f}  {Hm:7.2f}  {(Hm-Ho)/sd:6.2f}")
else:
    best     = None
    bg_best  = None

# ─── Verdict ──────────────────────────────────────────────────────────────────
if best:
    pass_desi  = best['tension'] < 2.0
    pass_theta = best['chi2_theta'] < 4.0
    beats_ref  = best['chi2'] < SIM121C_chi2
    if pass_desi and pass_theta:    verdict = "PASS"
    elif pass_theta and beats_ref:  verdict = "PARTIAL"
    elif beats_ref:                 verdict = "PARTIAL"
    else:                           verdict = "FAIL"
else:
    verdict = "FAIL"
    pass_desi = pass_theta = beats_ref = False

print()
print("="*72)
print(f"SIM135 RESULT:  {verdict}")
if best:
    print(f"  Best ξ_φ={best['xi']:.3f}:  tension={best['tension']:.3f}σ  "
          f"(Phase1: {SIM121C_tension:.2f}σ),  Δχ²={best['chi2']-SIM121C_chi2:+.3f}")
    print(f"  100θ*={best['theta']:.5f}  (obs 1.04101±0.00029)")
print("="*72)

# ─── Figures ──────────────────────────────────────────────────────────────────
xis    = [r['xi']      for r in valid]
tens   = [r['tension'] for r in valid]
thetas = [r['theta']   for r in valid]
F_cmbs = [r['F_cmb']   for r in valid]
phi0s  = [r['bg']['phi0'] for r in valid]
F0s    = [r['bg']['Feff0'] for r in valid]

if len(valid) > 1:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    ax = axes[0]
    ax.plot(xis, tens, 's-', color='#d73027', lw=2, ms=7, label='DESI tension')
    ax.axhline(SIM121C_tension, color='#2166ac', ls='--', lw=1.5,
               label=f'Phase1 {SIM121C_tension:.2f}σ')
    ax.axhline(2.0, color='green', ls=':', lw=1.5, label='2σ PASS')
    ax2 = ax.twinx()
    ax2.plot(xis, thetas, '^-', color='#4575b4', lw=2, ms=7, label=r'$100\theta_*$')
    ax2.axhline(theta_obs, color='#4575b4', ls='--', lw=1, alpha=0.6)
    ax2.set_ylabel(r'$100\theta_*$', color='#4575b4')
    ax.set_xlabel(r'$\xi_\phi$'); ax.set_ylabel('DESI tension [σ]')
    ax.set_title('SIM135: Tension and CMB scale vs coupling')
    lines1, labs1 = ax.get_legend_handles_labels()
    lines2, labs2 = ax2.get_legend_handles_labels()
    ax.legend(lines1+lines2, labs1+labs2, fontsize=9)
    ax.grid(alpha=0.3)

    ax = axes[1]
    ax.plot(xis, F_cmbs, 'o-', color='#d73027', lw=2, ms=7, label=r'$F_{\rm eff}(z_{\rm CMB})$')
    ax.plot(xis, F0s,    'v-', color='#1a9850', lw=2, ms=7, label=r'$F_{\rm eff}(z=0)$')
    ax.axhline(F0, color='#2166ac', ls='--', lw=1.5, label=f'Phase1 F₀={F0:.4f}')
    ax.set_xlabel(r'$\xi_\phi$'); ax.set_ylabel(r'$F_{\rm eff}$')
    ax.set_title('SIM135: F_eff evolution extent vs coupling')
    ax.legend(fontsize=9); ax.grid(alpha=0.3)

    plt.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT, f'sim135_tension_Feff.{ext}'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("  Saved sim135_tension_Feff")

# Figure 2: H(z) for all xi values vs DESI
if valid:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    ax = axes[0]
    colors = plt.cm.plasma(np.linspace(0.1, 0.85, len(valid)))
    z_fine = np.linspace(0.01, 2.6, 300)
    for r, col in zip(valid, colors):
        H_fine = [H_interp(z, r['bg']) for z in z_fine]
        ax.plot(z_fine, H_fine, lw=2, color=col,
                label=f'ξ={r["xi"]:.2f} ({r["tension"]:.2f}σ)')
    ax.errorbar(DESI_z, DESI_H, yerr=DESI_s, fmt='ko', ms=7, capsize=4,
                label='DESI', zorder=5)
    ax.set_xlabel(r'Redshift $z$'); ax.set_ylabel(r'$H(z)$ [km/s/Mpc]')
    ax.set_title('SIM135: H(z) vs DESI')
    ax.legend(fontsize=8); ax.grid(alpha=0.3)

    ax = axes[1]
    z_plot = bg_best['z'] if bg_best else valid[0]['bg']['z']
    bg_p   = bg_best        if bg_best else valid[0]['bg']
    mask   = (z_plot > 1e-3) & (z_plot < 1200)
    ax.semilogx(z_plot[mask]+1, bg_p['phi'][mask], color='#d73027', lw=2,
                label=f'φ(z) best ξ={best["xi"]:.3f}')
    ax_r = ax.twinx()
    ax_r.semilogx(z_plot[mask]+1, bg_p['Feff'][mask], color='#4575b4', lw=2, ls='--',
                   label=r'$F_{\rm eff}(z)$ best')
    ax_r.axhline(F0, color='#4575b4', ls=':', lw=1, alpha=0.5)
    ax_r.set_ylabel(r'$F_{\rm eff}$', color='#4575b4')
    ax.axvline(z_star+1, color='gray', ls=':', lw=1)
    ax.set_xlabel(r'$1+z$'); ax.set_ylabel(r'$\phi(z)$ [$M_{\rm Pl}$]')
    ax.set_title('SIM135: Field and F_eff evolution (best)')
    lines1, labs1 = ax.get_legend_handles_labels()
    lines2, labs2 = ax_r.get_legend_handles_labels()
    ax.legend(lines1+lines2, labs1+labs2, fontsize=9)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT, f'sim135_Hz_phi.{ext}'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("  Saved sim135_Hz_phi")

# ─── Save JSON ────────────────────────────────────────────────────────────────
out_json = {
    "sim": "SIM135",
    "phase": "Phase 3",
    "verdict": verdict,
    "description": "Bi-scalar: Ψ frozen at Phase 1 canonical (F₀=0.521); φ sourced by curvature",
    "action": "S=∫d⁴x√g[(F₀+ξ_φ·φ)R − ½(∂φ)²]+S_SM, F₀=0.521 (frozen), φ_ini=0",
    "F0": F0,
    "scan": [
        {
            "xi": r['xi'],
            "phi0": r['bg']['phi0'],
            "Feff0": r['bg']['Feff0'],
            "Feff_CMB": r['F_cmb'],
            "theta_100": r['theta'],
            "chi2_DESI": r['chi2'],
            "tension": r['tension'],
        }
        for r in valid
    ],
    "best": {
        "xi": best['xi'],
        "tension": best['tension'],
        "chi2_DESI": best['chi2'],
        "delta_chi2": best['chi2'] - SIM121C_chi2,
        "theta_100": best['theta'],
        "F_CMB": best['F_cmb'],
    } if best else None,
    "baseline": {"tension": SIM121C_tension, "chi2": SIM121C_chi2},
    "pass_desi":  bool(pass_desi),
    "pass_theta": bool(pass_theta),
    "beats_ref":  bool(beats_ref),
}
with open(os.path.join(OUT, 'sim135_results.json'), 'w') as f:
    json.dump(out_json, f, indent=2)
print("  Saved sim135_results.json")
print(f"\nAll outputs: {OUT}")
print("SIM135 complete.")
