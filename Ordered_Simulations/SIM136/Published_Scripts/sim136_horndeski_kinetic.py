"""
SIM136 — CMSTG Phase 3: Horndeski non-minimal kinetic coupling
=============================================================
Motivation: SIM123 showed Ψ frozen at Ψ≈0.001 throughout because curvature
sourcing is too weak to overcome Hubble friction for small Ψ.  The Horndeski
non-minimal kinetic coupling G^{μν}∂_μΨ∂_νΨ sources kinetic energy from
curvature directly — at large H the field acquires more kinetic energy.

Action:
  S = ∫d⁴x√g [ (1+2Λ₀Ψ²)/2 · R  −  ½(∂Ψ)²  +  α·G^{μν}∂_μΨ∂_νΨ ] + S_SM

In FRW, G^{μν}∂_μΨ∂_νΨ = 3H²Ψ̇² (G^{00}=3H² in (-,+,+,+)), so the kinetic
energy density is enhanced: ½Ψ̇²→½(1+6αH²)Ψ̇².  The Friedmann becomes:
  H²[ 3F_eff + 6F_eff'Ψ' − ½(1+6αH²)Ψ'² ] = ω_m a⁻³ + ω_r a⁻⁴ + Λ_bare

The Ψ EOM is unchanged at leading order:
  Ψ'' + (3−ε_H)Ψ' = F_eff'·6(2−ε_H)

(Kinetic coupling modifies Friedmann through enhanced kinetic energy; the
friction/source ratio in the EOM is unaffected to leading order in α.)

Stability: requires (1+6αH²) > 0; for α>0 always stable; for α<0 requires
|α| < 1/(6H²_CMB) ≈ 10⁻⁹ to avoid ghost at recombination.

Test A — Ψ_ini=0 (curvature-roll from zero):
  Does the modified Friedmann change H(z) enough to help DESI while Ψ evolves?
  α ∈ {1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3}, Λ₀=0.003, Ψ_ini=0

Test B — Ψ_ini=Ψ̄=2.62 (Phase 1 canonical):
  Does Horndeski modify H(z) when field is near-frozen at attractor?
  Same α scan.

Baseline: SIM121C DESI tension = 2.77σ, chi²=18.26
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

Lambda0  = 0.003
Psi_bar  = 2.62   # Phase 1 canonical field value
F0       = 0.5 + Lambda0 * Psi_bar**2   # = 0.52059
Lbare_phase1 = 3.0 * F0 * h_target**2 - omh2_m - omh2_r

# ─── Physics functions ────────────────────────────────────────────────────────

def F_eff(Psi):
    return 0.5 + Lambda0 * Psi**2

def F_prime(Psi):
    return 2.0 * Lambda0 * Psi

def get_H2_horndeski(Psi, y, N, Lbare, alpha, H2_guess=None):
    """
    Solve H² from modified Friedmann:
      H²D₀ - 3α·H⁴·y² = Q
    i.e. 3α·y²·H⁴ - D₀·H² + Q = 0
    taking the root that → Q/D₀ as α→0.
    Uses iterative solution: H²_new = (Q + 3α·H²_old²·y²) / D₀
    """
    a   = np.exp(N)
    Q   = omh2_m/a**3 + omh2_r/a**4 + Lbare
    F   = F_eff(Psi)
    Fp  = F_prime(Psi)
    D0  = 3.0*F + 6.0*Fp*y - 0.5*y**2
    if D0 < 0.01:
        D0 = max(3.0*F, 0.01)
    if H2_guess is None:
        H2 = Q / D0
    else:
        H2 = H2_guess

    if abs(alpha) < 1e-30 or abs(y) < 1e-30:
        return Q / D0

    # Iterate to convergence
    for _ in range(20):
        H2_new = (Q + 3.0*alpha*H2**2*y**2) / D0
        if abs(H2_new - H2) < 1e-14 * max(abs(H2), 1e-30):
            break
        H2 = 0.5*(H2 + H2_new)   # damped iteration for stability
    return max(H2_new, 1e-40)

def get_eps_H(Psi, y, N, Lbare, alpha):
    eps = 5e-4
    H2p = get_H2_horndeski(Psi, y, N+eps, Lbare, alpha)
    H2m = get_H2_horndeski(Psi, y, N-eps, Lbare, alpha)
    H2  = get_H2_horndeski(Psi, y, N,     Lbare, alpha)
    if H2 < 1e-40:
        return 0.0
    return -0.5*(H2p - H2m) / (2.0*eps*H2)

def ode(N, state, Lbare, alpha):
    Psi, y = state
    H2     = get_H2_horndeski(Psi, y, N, Lbare, alpha)
    eps_H  = get_eps_H(Psi, y, N, Lbare, alpha)
    eps_H  = np.clip(eps_H, -1.0, 4.0)
    # Check ghost stability: (1+6αH²) > 0
    kin_factor = 1.0 + 6.0*alpha*H2
    if kin_factor < 0.01:
        return [0.0, 0.0]   # ghost regime: abort field evolution
    # Standard EOM (kinetic coupling modifies Friedmann but not EOM at leading order)
    source = F_prime(Psi) * 6.0*(2.0 - eps_H)
    dy_dN  = source - (3.0 - eps_H)*y
    return [y, dy_dN]

def integrate(Psi_ini, Lbare, alpha, z_ini=1e5, n_pts=4000):
    N_ini  = np.log(1.0/(1.0+z_ini))
    N_eval = np.linspace(N_ini, 0.0, n_pts)
    sol = solve_ivp(ode, (N_ini, 0.0), [Psi_ini, 0.0],
                    args=(Lbare, alpha), t_eval=N_eval,
                    method='RK45', rtol=1e-8, atol=1e-12, max_step=0.05)
    if not sol.success:
        return None
    N_arr  = sol.t
    z_arr  = np.exp(-N_arr) - 1.0
    ps_arr = sol.y[0]
    y_arr  = sol.y[1]
    H_arr  = np.array([H100*np.sqrt(max(get_H2_horndeski(ps_arr[i],y_arr[i],N_arr[i],Lbare,alpha),0.0))
                       for i in range(len(N_arr))])
    Fv_arr = np.array([F_eff(ps_arr[i]) for i in range(len(N_arr))])
    idx    = np.argsort(z_arr)
    return dict(z=z_arr[idx], N=N_arr[idx], psi=ps_arr[idx], y=y_arr[idx],
                H=H_arr[idx], Feff=Fv_arr[idx],
                psi0=float(ps_arr[-1]), y0=float(y_arr[-1]),
                Feff0=float(Fv_arr[-1]), H0=float(H_arr[-1]),
                Lbare=Lbare, alpha=alpha)

def calibrate(psi0, y0, alpha, bg):
    H0   = bg['H0']/H100   # h value
    H2_0 = H0**2
    F    = F_eff(psi0)
    Fp   = F_prime(psi0)
    D0   = 3.0*F + 6.0*Fp*y0 - 0.5*y0**2
    Lbare_new = h_target**2 * D0 - omh2_m - omh2_r
    return Lbare_new

def run(Psi_ini, alpha, n_iter=8):
    Lbare = Lbare_phase1
    bg    = None
    for _ in range(n_iter):
        bg = integrate(Psi_ini, Lbare, alpha)
        if bg is None:
            return None
        Lbare_new = calibrate(bg['psi0'], bg['y0'], alpha, bg)
        if abs(Lbare_new - Lbare) < 1e-10:
            break
        Lbare = Lbare_new
    return integrate(Psi_ini, Lbare, alpha)

def H_interp(z, bg):
    return float(np.interp(z, bg['z'], bg['H']))

def theta_star(bg):
    H0   = bg['H0']; h = H0/100.0
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

def run_test(label, Psi_ini, alpha_grid):
    print(f"\n  {'─'*64}")
    print(f"  Test {label}: Ψ_ini = {Psi_ini:.3f} M_Pl")
    print(f"  {'α':>10}  {'Ψ₀':>8}  {'F_eff0':>8}  {'F_CMB':>8}  {'100θ*':>8}  {'tension':>8}  {'kin_z0':>8}")
    print(f"  {'─'*64}")

    results = []
    for alpha in alpha_grid:
        bg = run(Psi_ini, alpha)
        if bg is None:
            print(f"  {alpha:10.2e}  FAILED")
            results.append(None)
            continue

        F_cmb = float(np.interp(z_star, bg['z'], bg['Feff']))
        H2_0  = (bg['H0']/H100)**2
        kin0  = 1.0 + 6.0*alpha*H2_0   # kinetic factor at z=0

        th  = theta_star(bg)
        c2  = chi2_DESI(bg)
        ten = np.sqrt(c2/len(DESI_z))
        c2t = ((th-theta_obs)/theta_obs_err)**2 if not np.isnan(th) else np.nan

        print(f"  {alpha:10.2e}  {bg['psi0']:8.4f}  {bg['Feff0']:8.5f}  "
              f"{F_cmb:8.5f}  {th:8.5f}  {ten:8.3f}σ  {kin0:8.4f}")

        results.append(dict(alpha=alpha, bg=bg, theta=th, chi2=c2,
                            tension=ten, chi2_theta=c2t, F_cmb=F_cmb, kin0=kin0))
    return results

# ─── Main ─────────────────────────────────────────────────────────────────────
print("="*72)
print("SIM136 — Phase 3: Horndeski non-minimal kinetic coupling")
print("="*72)
print(f"  Action: F(Ψ)R − ½(∂Ψ)² + α·G^{{μν}}∂_μΨ∂_νΨ")
print(f"  F(Ψ) = (1+2Λ₀Ψ²)/2,  Λ₀={Lambda0}")
print(f"  Friedmann: H²[3F+6F'Ψ'−½(1+6αH²)Ψ'²] = ρ_tot")
print(f"  Baseline: {SIM121C_tension:.2f}σ, χ²={SIM121C_chi2:.2f}")

alpha_grid = [1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3]

# Test A: Ψ_ini = 0
res_A = run_test("A", 0.0, alpha_grid)

# Test B: Ψ_ini = Phase 1 canonical
res_B = run_test("B", Psi_bar, alpha_grid)

# ─── Combined verdict ─────────────────────────────────────────────────────────
all_valid = [r for r in (res_A + res_B) if r is not None and not np.isnan(r['theta'])]

if all_valid:
    best = min(all_valid, key=lambda r: r['chi2'] + r['chi2_theta'])
    print(f"\n  Best overall: α={best['alpha']:.2e},  tension={best['tension']:.3f}σ,  "
          f"Δχ²={best['chi2']-SIM121C_chi2:+.3f},  100θ*={best['theta']:.5f}")
    bg_best = best['bg']
else:
    best     = None
    bg_best  = None

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
print(f"SIM136 RESULT:  {verdict}")
if best:
    print(f"  Best: α={best['alpha']:.2e},  tension={best['tension']:.3f}σ  "
          f"(Phase1: {SIM121C_tension:.2f}σ),  Δχ²={best['chi2']-SIM121C_chi2:+.3f}")
    print(f"  100θ*={best['theta']:.5f}  (obs 1.04101±0.00029)")
print("="*72)

# ─── Figures ──────────────────────────────────────────────────────────────────
valid_A = [r for r in res_A if r is not None and not np.isnan(r['theta'])]
valid_B = [r for r in res_B if r is not None and not np.isnan(r['theta'])]

if valid_A or valid_B:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    ax = axes[0]
    if valid_A:
        alphas_A = [r['alpha'] for r in valid_A]
        tens_A   = [r['tension'] for r in valid_A]
        th_A     = [r['theta'] for r in valid_A]
        ax.semilogx(alphas_A, tens_A, 's-', color='#d73027', lw=2, ms=7, label='Test A (Ψ_ini=0)')
    if valid_B:
        alphas_B = [r['alpha'] for r in valid_B]
        tens_B   = [r['tension'] for r in valid_B]
        ax.semilogx(alphas_B, tens_B, 'o-', color='#1a9850', lw=2, ms=7,
                    label=f'Test B (Ψ_ini={Psi_bar})')
    ax.axhline(SIM121C_tension, color='#2166ac', ls='--', lw=1.5,
               label=f'SIM121C {SIM121C_tension:.2f}σ')
    ax.axhline(2.0, color='green', ls=':', lw=1.5, label='2σ PASS')
    ax.set_xlabel(r'Horndeski coupling $\alpha$'); ax.set_ylabel('DESI tension [σ]')
    ax.set_title('SIM136: DESI tension vs Horndeski α')
    ax.legend(fontsize=9); ax.grid(alpha=0.3)

    ax = axes[1]
    if valid_A:
        thetas_A = [r['theta'] for r in valid_A]
        Fcmbs_A  = [r['F_cmb'] for r in valid_A]
        ax.semilogx(alphas_A, thetas_A, 's-', color='#d73027', lw=2, ms=7,
                    label=r'$100\theta_*$ A (Ψ_ini=0)')
    if valid_B:
        thetas_B = [r['theta'] for r in valid_B]
        Fcmbs_B  = [r['F_cmb'] for r in valid_B]
        ax.semilogx(alphas_B, thetas_B, 'o-', color='#1a9850', lw=2, ms=7,
                    label=r'$100\theta_*$ B (Ψ_ini=2.62)')
    ax.axhline(theta_obs, color='#2166ac', ls='--', lw=1.5, label=f'Planck {theta_obs}')
    ax.axhline(theta_obs - 2*theta_obs_err, color='#2166ac', ls=':', lw=1)
    ax.axhline(theta_obs + 2*theta_obs_err, color='#2166ac', ls=':', lw=1)
    ax.set_xlabel(r'Horndeski coupling $\alpha$'); ax.set_ylabel(r'$100\theta_*$')
    ax.set_title('SIM136: CMB acoustic scale vs Horndeski α')
    ax.legend(fontsize=9); ax.grid(alpha=0.3)

    plt.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT, f'sim136_tension_theta.{ext}'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("  Saved sim136_tension_theta")

# H(z) figure for best
if bg_best is not None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    ax = axes[0]
    z_fine = np.linspace(0.01, 2.6, 300)
    colors_A = plt.cm.Reds(np.linspace(0.4, 0.9, len(valid_A)))
    colors_B = plt.cm.Greens(np.linspace(0.4, 0.9, len(valid_B)))
    for r, col in zip(valid_A, colors_A):
        H_fine = [H_interp(z, r['bg']) for z in z_fine]
        ax.plot(z_fine, H_fine, lw=1.5, color=col, ls='-',
                label=f'A α={r["alpha"]:.1e} ({r["tension"]:.2f}σ)')
    for r, col in zip(valid_B[:3], colors_B[:3]):
        H_fine = [H_interp(z, r['bg']) for z in z_fine]
        ax.plot(z_fine, H_fine, lw=1.5, color=col, ls='--',
                label=f'B α={r["alpha"]:.1e} ({r["tension"]:.2f}σ)')
    ax.errorbar(DESI_z, DESI_H, yerr=DESI_s, fmt='ko', ms=7, capsize=4,
                label='DESI', zorder=5)
    ax.set_xlabel(r'Redshift $z$'); ax.set_ylabel(r'$H(z)$ [km/s/Mpc]')
    ax.set_title('SIM136: H(z) comparison')
    ax.legend(fontsize=7, ncol=2); ax.grid(alpha=0.3)

    ax = axes[1]
    z_plot = bg_best['z']
    mask   = (z_plot > 1e-3) & (z_plot < 1200)
    ax.semilogx(z_plot[mask]+1, bg_best['psi'][mask], color='#d73027', lw=2,
                label=f'Ψ(z), best α={best["alpha"]:.2e}')
    ax_r = ax.twinx()
    ax_r.semilogx(z_plot[mask]+1, bg_best['Feff'][mask], color='#4575b4', lw=2, ls='--',
                   label=r'$F_{\rm eff}(z)$ best')
    ax_r.axhline(F0, color='#4575b4', ls=':', lw=1, alpha=0.5, label=f'F₀={F0:.4f}')
    ax_r.set_ylabel(r'$F_{\rm eff}$', color='#4575b4')
    ax.axvline(z_star+1, color='gray', ls=':', lw=1)
    ax.set_xlabel(r'$1+z$'); ax.set_ylabel(r'$\Psi(z)$ [$M_{\rm Pl}$]')
    ax.set_title('SIM136: Best field trajectory')
    lines1, labs1 = ax.get_legend_handles_labels()
    lines2, labs2 = ax_r.get_legend_handles_labels()
    ax.legend(lines1+lines2, labs1+labs2, fontsize=9)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT, f'sim136_Hz_psi.{ext}'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("  Saved sim136_Hz_psi")

# ─── Save JSON ────────────────────────────────────────────────────────────────
def fmt_r(r):
    if r is None: return None
    return {"alpha": r['alpha'], "psi0": r['bg']['psi0'],
            "Feff0": r['bg']['Feff0'], "Feff_CMB": r['F_cmb'],
            "theta_100": r['theta'], "chi2_DESI": r['chi2'],
            "tension": r['tension'], "kin_factor_z0": r['kin0']}

out_json = {
    "sim": "SIM136",
    "phase": "Phase 3",
    "verdict": verdict,
    "description": "Horndeski kinetic: G^{μν}∂_μΨ∂_νΨ coupling; Phase 1 F backbone",
    "action": "S=∫d⁴x√g[(1+2Λ₀Ψ²)/2·R − ½(∂Ψ)² + α·G^{μν}∂_μΨ∂_νΨ]+S_SM",
    "Lambda0": Lambda0,
    "scan_A_psi_ini0": [fmt_r(r) for r in res_A],
    "scan_B_psi_ini_canonical": [fmt_r(r) for r in res_B],
    "best": {"alpha": best['alpha'], "tension": best['tension'],
             "chi2_DESI": best['chi2'], "delta_chi2": best['chi2']-SIM121C_chi2,
             "theta_100": best['theta']} if best else None,
    "baseline": {"tension": SIM121C_tension, "chi2": SIM121C_chi2},
    "pass_desi":  bool(pass_desi),
    "pass_theta": bool(pass_theta),
    "beats_ref":  bool(beats_ref),
}
with open(os.path.join(OUT, 'sim136_results.json'), 'w') as f:
    json.dump(out_json, f, indent=2)
print("  Saved sim136_results.json")
print(f"\nAll outputs: {OUT}")
print("SIM136 complete.")
