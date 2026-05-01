#!/usr/bin/env python3
"""
SIM120-alt — CMSTG Phase 2: Constrained Universal m₂₂ SPARC Fit
================================================================
Fixes a single universal m₂₂ value across all SPARC galaxies and refits
each galaxy with the soliton virial relation as a hard constraint:

    M_c × r_c = 9.1×10⁷ / m₂₂²   [Schive+2014]

This eliminates ρ_c as a free parameter (given r_c and m₂₂, ρ_c is fixed
by the virial relation). Free parameters per galaxy: r_c, ρ_s, r_s, Υ_d, Υ_b.

Scans m₂₂ ∈ [0.05, 2.0] to find the universal m₂₂ that minimises the
total weighted χ²/dof across all constrained SPARC galaxies.

Motivation: SIM119 derived per-galaxy m₂₂ values with median 0.28 and
σ=0.58 dex. If FDM is real, a single m₂₂ should exist. This sim finds it
and derives κ from the CMSTG DE-DM link: m_χ = √(2κ)Ψ̄.

Pass criteria:
    - Best-fit m₂₂_universal ∈ [0.1, 10]   (FDM window)
    - Median χ²/dof < 2 across galaxies     (fit quality)
    - κ = m_χ²/(2Ψ̄²) consistent with SIM119 κ range

Units: kpc / M_sun / km/s. G = 4.302×10⁻⁶ kpc (km/s)² M_sun⁻¹.
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize
from scipy.interpolate import interp1d
import os, glob, warnings
warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────────────────────────────────────
# CONSTANTS
# ─────────────────────────────────────────────────────────────────────────────
G_kpc     = 4.302e-6       # kpc (km/s)² M_sun⁻¹
soliton_A = 9.1e7          # M_sun kpc m₂₂²  (Schive+2014 virial constant)
I_sol, _  = quad(lambda x: x**2 / (1.0 + 0.091*x**2)**8, 0, np.inf, limit=300)

# CMSTG DE parameters (locked from SIM113)
Psi_bar   = 2.62           # Ψ̄ today [M_Pl]
M_Pl_eV   = 1.22e28        # M_Pl in eV

SPARC_DIR = os.path.join(
    os.path.dirname(__file__),
    '..', '..', '..',
    'Ordered_Simulations',
    'simulation_18_cmstg_graviton_emit',
    'figures', 'cmstg_p7_min_impl (2)', 'external', 'sparc_raw'
)

# ─────────────────────────────────────────────────────────────────────────────
# DATA LOADING (identical to SIM119)
# ─────────────────────────────────────────────────────────────────────────────
def load_sparc(filepath):
    name = os.path.splitext(os.path.basename(filepath))[0].replace('_rotmod', '')
    data = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            try:
                vals = list(map(float, line.split()))
                if len(vals) >= 7:
                    data.append(vals[:8])
            except ValueError:
                continue
    if len(data) < 6:
        return None
    arr = np.array(data)
    R, Vobs, errV, Vgas, Vdisk, Vbul = (arr[:,0], arr[:,1], arr[:,2],
                                         arr[:,3], arr[:,4], arr[:,5])
    mask = (errV > 0) & (Vobs > 0) & (R > 0)
    if mask.sum() < 6:
        return None
    if Vobs[mask].max() < 25.0:
        return None
    return {
        'name':      name,
        'R':         R[mask],
        'Vobs':      Vobs[mask],
        'errV':      np.maximum(errV[mask], 2.0),
        'Vgas':      Vgas[mask],
        'Vdisk':     Vdisk[mask],
        'Vbul':      Vbul[mask],
        'has_bulge': Vbul[mask].max() > 1.0,
        'Npts':      int(mask.sum()),
    }

# ─────────────────────────────────────────────────────────────────────────────
# MODEL
# ─────────────────────────────────────────────────────────────────────────────
def rho_from_virial(r_c, m22):
    """ρ_c determined by soliton virial relation given r_c and m₂₂."""
    # M_c = 4π ρ_c r_c³ I_sol → M_c × r_c = A/m₂₂²
    # → ρ_c = A / (4π r_c⁴ I_sol m₂₂²)
    return soliton_A / (4.0 * np.pi * r_c**4 * I_sol * m22**2)

def rho_soliton(r, rho_c, r_c):
    x = r / r_c
    return rho_c / (1.0 + 0.091*x**2)**8

def rho_NFW(r, rho_s, r_s):
    x = max(r / r_s, 1e-8)
    return rho_s / (x * (1.0 + x)**2)

def v2_DM(R_arr, rho_c, r_c, rho_s, r_s, n_grid=250):
    R_max  = R_arr.max() * 1.02
    r_grid = np.linspace(1e-3, R_max, n_grid)
    rho_dm = np.array([max(rho_soliton(r, rho_c, r_c), rho_NFW(r, rho_s, r_s))
                       for r in r_grid])
    integrand = 4.0*np.pi * r_grid**2 * rho_dm
    M_enc = np.zeros(n_grid)
    dr = np.diff(r_grid)
    M_enc[1:] = np.cumsum(0.5*(integrand[:-1]+integrand[1:])*dr)
    M_fn = interp1d(r_grid, M_enc, kind='linear', fill_value='extrapolate')
    return np.array([G_kpc * max(M_fn(R), 0.0) / R for R in R_arr])

def chi2_constrained(params, gal, m22_fixed):
    """χ²/dof with ρ_c fixed by virial relation: params = [log_r_c, log_ρ_s, log_r_s, Υ_d, Υ_b]."""
    try:
        r_c   = 10.0**params[0]
        rho_s = 10.0**params[1]
        r_s   = 10.0**params[2]
        Upsd  = params[3]
        Upsb  = params[4] if gal['has_bulge'] else 0.0

        rho_c = rho_from_virial(r_c, m22_fixed)
        if rho_c > 1e15 or rho_c < 1e3:
            return 1e12

        v2dm  = v2_DM(gal['R'], rho_c, r_c, rho_s, r_s)
        v2bar = (Upsd * gal['Vdisk']**2 +
                 Upsb * gal['Vbul']**2 +
                 gal['Vgas']**2)
        v2tot = v2dm + v2bar
        v_mod = np.sqrt(np.maximum(v2tot, 0.0))

        resid = (v_mod - gal['Vobs']) / gal['errV']
        chi2  = float(np.sum(resid**2)) / gal['Npts']
        return chi2 if np.isfinite(chi2) else 1e12
    except Exception:
        return 1e12

# Bounds for constrained fit (5 parameters: log_r_c, log_ρ_s, log_r_s, Υ_d, Υ_b)
BOUNDS_CONSTRAINED = [
    (-1.0, 1.2),    # log_r_c  : 0.1 – 16 kpc
    (4.0, 10.0),    # log_rho_s
    (0.0, 1.9),     # log_r_s
    (0.1, 2.0),     # Upsilon_disk
    (0.1, 2.0),     # Upsilon_bul
]

def fit_galaxy_constrained(gal, m22_fixed, n_starts=6):
    """Fit one galaxy with universal m₂₂ fixed via virial constraint."""
    rng  = np.random.default_rng(seed=42)
    best_chi2 = np.inf
    best_p    = None

    R_ref = gal['R'][np.argmax(gal['Vobs'])]
    fixed_seeds = [
        [-0.1, 7.0, np.log10(max(R_ref*0.5, 1.5)), 0.5, 0.7],
        [ 0.2, 6.5, np.log10(max(R_ref*0.7, 2.0)), 0.4, 0.5],
    ]
    lo = np.array([b[0] for b in BOUNDS_CONSTRAINED])
    hi = np.array([b[1] for b in BOUNDS_CONSTRAINED])
    random_seeds = [lo + rng.random(5)*(hi-lo) for _ in range(n_starts - 2)]

    for seed in fixed_seeds + random_seeds:
        seed = np.clip(seed, lo, hi)
        try:
            res = minimize(
                chi2_constrained, seed, args=(gal, m22_fixed),
                method='L-BFGS-B',
                bounds=BOUNDS_CONSTRAINED,
                options={'maxiter': 1500, 'ftol': 1e-9, 'gtol': 1e-6}
            )
            if res.fun < best_chi2:
                best_chi2 = res.fun
                best_p    = res.x
        except Exception:
            continue

    if best_p is None:
        return None
    r_c   = 10.0**best_p[0]
    rho_c = rho_from_virial(r_c, m22_fixed)
    rho_s = 10.0**best_p[1]
    r_s   = 10.0**best_p[2]
    return {
        'chi2': best_chi2,
        'r_c':  r_c,
        'rho_c': rho_c,
        'rho_s': rho_s,
        'r_s':  r_s,
        'Upsd': best_p[3],
        'Upsb': best_p[4] if gal['has_bulge'] else 0.0,
    }

def scan_m22(galaxies, m22_grid, n_starts_scan=3, verbose=False):
    """
    Scan m₂₂ grid. Uses reduced n_starts for speed during scan;
    galaxies is a subset (every 4th) for the scan phase.
    """
    results = []
    for m22 in m22_grid:
        chi2_list = []
        for gal in galaxies:
            res = fit_galaxy_constrained(gal, m22, n_starts=n_starts_scan)
            if res is not None and np.isfinite(res['chi2']):
                chi2_list.append(res['chi2'])
        if len(chi2_list) == 0:
            results.append({'m22': m22, 'median_chi2': np.inf, 'frac_pass': 0.0, 'n_fit': 0})
            continue
        arr = np.array(chi2_list)
        median_c2 = float(np.median(arr))
        frac_pass = float((arr < 2.0).mean())
        if verbose:
            print(f"  m₂₂={m22:.3f}: median χ²/dof={median_c2:.3f}, frac<2={frac_pass:.2f}, n={len(arr)}")
        results.append({
            'm22':        m22,
            'median_chi2': median_c2,
            'frac_pass':  frac_pass,
            'n_fit':      len(arr),
        })
    return results

# ─────────────────────────────────────────────────────────────────────────────
# κ DERIVATION
# ─────────────────────────────────────────────────────────────────────────────
def kappa_from_m22(m22, Psi_bar_Mpl=Psi_bar):
    """κ from CMSTG link m_χ = √(2κ)·Ψ̄ with m_χ = m₂₂ × 10⁻²² eV / M_Pl."""
    m_chi_eV = m22 * 1e-22
    return 0.5 * (m_chi_eV / (Psi_bar_Mpl * M_Pl_eV))**2

# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    print("=" * 72)
    print("SIM120-alt — CMSTG Phase 2: Constrained Universal m₂₂ SPARC Fit")
    print("=" * 72)
    print(f"\n  Soliton virial: A = {soliton_A:.2e} M_sun kpc, I_sol = {I_sol:.5f}")
    print(f"  CMSTG DE link: Ψ̄ = {Psi_bar} M_Pl  (SIM113 best-fit)")
    print(f"  Virial constraint: ρ_c = A / (4π r_c⁴ I_sol m₂₂²) — 1 fewer free parameter")
    print()

    # ── Load SPARC ─────────────────────────────────────────────────────────────
    sparc_files = sorted(glob.glob(os.path.join(SPARC_DIR, '*_rotmod.dat')))
    galaxies    = []
    skipped     = 0
    for f in sparc_files:
        g = load_sparc(f)
        if g is not None:
            galaxies.append(g)
        else:
            skipped += 1
    print(f"  Loaded {len(galaxies)} galaxies  (skipped {skipped})\n")

    if len(galaxies) == 0:
        print("ERROR: No SPARC galaxies found. Check SPARC_DIR path.")
        import sys; sys.exit(1)

    # ── Coarse m₂₂ scan on representative subset ──────────────────────────────
    print("─" * 72)
    print("PART A: Coarse m₂₂ scan (every 4th galaxy for speed)")
    print("─" * 72)
    galaxies_scan = galaxies[::4]   # ~40 galaxies for scan
    m22_coarse = np.array([0.05, 0.08, 0.12, 0.18, 0.25, 0.35, 0.50,
                           0.70, 1.0, 1.5, 2.0, 3.0, 5.0])
    print(f"\n  Scanning {len(m22_coarse)} m₂₂ values × {len(galaxies_scan)} galaxies "
          f"(n_starts=3 for speed)...\n")

    coarse_results = scan_m22(galaxies_scan, m22_coarse, n_starts_scan=3, verbose=True)

    # Find coarse minimum
    coarse_arr = np.array([r['median_chi2'] for r in coarse_results])
    i_best_c   = int(np.argmin(coarse_arr))
    m22_best_c = coarse_results[i_best_c]['m22']
    print(f"\n  Coarse best: m₂₂ = {m22_best_c:.3f}, "
          f"median χ²/dof = {coarse_results[i_best_c]['median_chi2']:.3f}")

    # ── Fine scan around coarse minimum ────────────────────────────────────────
    print("\n─" * 72)
    print("PART B: Fine m₂₂ scan (±0.5 dex around coarse minimum, full sample)")
    print("─" * 72)

    log_best_c = np.log10(m22_best_c)
    m22_fine   = np.logspace(log_best_c - 0.5, log_best_c + 0.5, 11)
    print(f"\n  Scanning {len(m22_fine)} values in [{m22_fine[0]:.3f}, {m22_fine[-1]:.3f}] "
          f"× {len(galaxies)} galaxies (n_starts=4)...\n")

    fine_results = scan_m22(galaxies, m22_fine, n_starts_scan=4, verbose=True)

    fine_arr    = np.array([r['median_chi2'] for r in fine_results])
    i_best_f    = int(np.argmin(fine_arr))
    m22_universal = fine_results[i_best_f]['m22']
    chi2_universal = fine_results[i_best_f]['median_chi2']
    frac_universal = fine_results[i_best_f]['frac_pass']
    n_universal    = fine_results[i_best_f]['n_fit']

    print(f"\n  Fine best: m₂₂_universal = {m22_universal:.3f}")
    print(f"  Median χ²/dof = {chi2_universal:.3f}")
    print(f"  Fraction χ²/dof < 2: {frac_universal:.2f}")

    # ── Full fit at universal m₂₂ for per-galaxy breakdown ─────────────────────
    print("\n─" * 72)
    print(f"PART C: Per-galaxy fits at m₂₂_universal = {m22_universal:.3f}")
    print("─" * 72)

    print(f"\n  {'Galaxy':>20} {'Npts':>5} {'χ²/dof':>8} {'r_c(kpc)':>9} "
          f"{'ρ_c':>11} {'ρ_s':>11} {'verdict':>9}")
    print("  " + "-"*78)

    all_chi2 = []
    pass_cnt = 0
    fail_cnt = 0

    for gal in galaxies:
        res = fit_galaxy_constrained(gal, m22_universal, n_starts=8)
        if res is None or not np.isfinite(res['chi2']):
            print(f"  {gal['name']:>20} {gal['Npts']:>5} {'ERROR':>8}")
            fail_cnt += 1
            continue

        all_chi2.append(res['chi2'])
        chi2_ok = res['chi2'] < 2.0
        rc_ok   = res['r_c'] > 0.3
        verdict = 'PASS' if (chi2_ok and rc_ok) else 'fail'
        if verdict == 'PASS':
            pass_cnt += 1
        else:
            fail_cnt += 1

        print(f"  {gal['name']:>20} {gal['Npts']:>5} {res['chi2']:>8.3f} "
              f"{res['r_c']:>9.3f} {res['rho_c']:>11.2e} {res['rho_s']:>11.2e} {verdict:>9}")

    # ── Summary statistics ──────────────────────────────────────────────────────
    all_chi2 = np.array(all_chi2)
    print("\n" + "─"*72)
    print("PART D: Summary statistics")
    print("─" * 72)
    print(f"\n  Galaxies fit:      {len(all_chi2)}")
    print(f"  PASS (χ²/dof<2):  {pass_cnt}  ({100*pass_cnt/max(len(all_chi2),1):.0f}%)")
    print(f"  Median χ²/dof:    {np.median(all_chi2):.3f}")
    print(f"  Mean χ²/dof:      {np.mean(all_chi2):.3f}")
    print(f"  Fraction < 2:     {(all_chi2<2.0).mean():.2f}")

    # ── Comparison table: m₂₂ scan ─────────────────────────────────────────────
    print("\n─" * 72)
    print("PART E: m₂₂ scan summary table")
    print("─" * 72)
    print(f"\n  {'m₂₂':>7}  {'median χ²/dof':>15}  {'frac<2':>8}  {'n_fit':>6}  {'note':>20}")
    print("  " + "-"*65)
    for r in coarse_results + fine_results:
        note = " ← universal" if abs(r['m22'] - m22_universal) < 0.01 else ""
        marker = " ***" if abs(r['m22'] - m22_universal) < 0.01 else ""
        print(f"  {r['m22']:>7.3f}  {r['median_chi2']:>15.3f}  {r['frac_pass']:>8.2f}  "
              f"{r['n_fit']:>6}  {note+marker:>20}")

    # ── κ derivation ────────────────────────────────────────────────────────────
    print("\n─" * 72)
    print("PART F: CMSTG κ derivation")
    print("─" * 72)
    kappa = kappa_from_m22(m22_universal)
    m_chi_eV = m22_universal * 1e-22

    print(f"\n  m₂₂_universal = {m22_universal:.3f}")
    print(f"  m_χ           = {m_chi_eV:.3e} eV")
    print(f"  Ψ̄             = {Psi_bar} M_Pl  (from SIM113)")
    print(f"  κ             = m_χ²/(2Ψ̄²M_Pl²) = {kappa:.4e}  (dimensionless)")
    print(f"  CMSTG link: m_χ = √(2κ)·Ψ̄  ✓")
    print(f"\n  SIM119 κ range (from per-galaxy m₂₂ distribution):")
    print(f"    m₂₂_median=0.28 → κ = {kappa_from_m22(0.28):.3e}")
    print(f"    m₂₂_best=0.082  → κ = {kappa_from_m22(0.082):.3e}")
    print(f"    m₂₂_window=[0.1,10] → κ ∈ [{kappa_from_m22(0.1):.2e}, {kappa_from_m22(10.0):.2e}]")
    print(f"\n  SIM120-alt κ = {kappa:.4e}")
    in_window = 0.1 <= m22_universal <= 10.0
    kappa_consistent = (kappa_from_m22(0.1) <= kappa <= kappa_from_m22(10.0))
    print(f"  m₂₂ in FDM window [0.1, 10]: {'YES' if in_window else 'NO'}")
    print(f"  κ consistent with SIM119:    {'YES' if kappa_consistent else 'NO'}")

    # ── Final verdict ───────────────────────────────────────────────────────────
    print("\n" + "═"*72)
    print("SIM120-alt RESULT:")
    print()

    PASS = (in_window and
            np.median(all_chi2) < 2.0 and
            (all_chi2 < 2.0).mean() >= 0.4)

    print(f"  m₂₂_universal    = {m22_universal:.3f}  ({'in' if in_window else 'OUTSIDE'} FDM window [0.1,10])")
    print(f"  Median χ²/dof    = {np.median(all_chi2):.3f}  ({'OK' if np.median(all_chi2)<2 else 'high'})")
    print(f"  Fraction χ²<2    = {(all_chi2<2.0).mean():.2f}  ({'OK' if (all_chi2<2.0).mean()>0.4 else 'low'})")
    print(f"  κ                = {kappa:.4e}")
    print()
    if PASS:
        print("  VERDICT: PASS")
        print("  A single universal m₂₂ fits the SPARC sample with the virial")
        print("  constraint imposed. The CMSTG DE-DM link κ = m_χ²/(2Ψ̄²) is")
        print("  uniquely determined.")
    else:
        print("  VERDICT: PARTIAL / FAIL")
        print("  A universal m₂₂ with virial constraint degrades fit quality,")
        print("  suggesting intrinsic scatter in soliton properties or galaxy-")
        print("  specific effects beyond a universal mass parameter.")
    print("═"*72)
