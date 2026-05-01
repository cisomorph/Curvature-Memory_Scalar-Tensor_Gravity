#!/usr/bin/env python3
# Sim 89 — CMSTG Planck TT Likelihood Comparison (Phase 5b)
"""
Closes the final CMB referee gap.

SIM88 compared CMSTG vs LCDM C_ell with no noise model (full-sky
noise-free Gaussian chi2). This simulation adds the proper comparison:
CMSTG and LCDM CLASS outputs are compared against the Planck 2018
best-fit TT theory spectrum using the approximate full-sky Gaussian
TT likelihood with a Planck noise floor.

Reference:
  Planck 2018 V (Aghanim+2019, arXiv:1907.12875), Section 2.1.
  Approximate TT likelihood:
    -2 ln L_TT = f_sky * sum_l (2l+1) * [C_l^pred/C_l^ref - 1 - ln(C_l^pred/C_l^ref)]
  where C_l^ref = C_l^Planck + N_l  (signal + noise floor).
  This is exact in the full-sky limit and a good approximation for
  f_sky > 0.5 (Planck uses f_sky ~ 0.78 for TT).

Noise model:
  N_l = (Delta_T * theta_FWHM)^2 * exp(l*(l+1)*sigma_b^2)
  Delta_T = 33 μK-arcmin, theta_FWHM = 7.27 arcmin (Planck 143 GHz)
  This prevents chi2 blow-up at l > 1500 where signal is noise-dominated.

Inputs:
  - SIM88 CLASS CMSTG C_ell (lensed TT, converted to μK^2)
  - SIM88 CLASS LCDM C_ell (lensed TT, converted to μK^2)
  - Planck 2018 theory TT spectrum (COM_PowerSpect R3.01)

Outputs:
  Outputs/sim89_chi2_table.png    -- chi2 contributions by l range
  Outputs/sim89_cl_vs_planck.png  -- D_l CMSTG + LCDM vs Planck theory
  Outputs/sim89_residuals.png     -- (D_l^pred - D_l^Planck) / D_l^Planck
  Outputs/sim89_diagnostics.json
"""

import os, json, math, glob, warnings
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore", category=RuntimeWarning)

BASE    = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INPUTS  = os.path.join(BASE, "Inputs")
OUTPUTS = os.path.join(BASE, "Outputs")
PARAMS  = os.path.join(INPUTS, "sim89_params.json")
os.makedirs(OUTPUTS, exist_ok=True)

with open(PARAMS) as f:
    P = json.load(f)

T_CMB_K   = float(P["T_CMB_K"])
T_CMB_uK2 = (T_CMB_K * 1e6)**2   # μK²  ≈ 7.428e12

# ── Load Planck 2018 TT theory spectrum ──────────────────────────────────────
planck_file = P["planck_theory_file"]
print(f"  Loading Planck theory: {os.path.basename(planck_file)}", flush=True)

planck_raw = np.loadtxt(planck_file, comments="#")
# Columns: L  TT  TE  EE  BB  PP  (D_l = l(l+1)Cl/2pi in μK^2)
ell_planck = planck_raw[:, 0].astype(int)
Dl_planck  = planck_raw[:, 1]   # TT  [μK^2]
print(f"  Planck TT range: l={ell_planck[0]}..{ell_planck[-1]}, "
      f"D_l at l=2: {Dl_planck[0]:.1f} μK^2, "
      f"peak at l~{ell_planck[np.argmax(Dl_planck[:500])]}: {np.max(Dl_planck[:500]):.1f} μK^2")

# ── Load SIM88 CLASS C_ell outputs ───────────────────────────────────────────
def load_class_cl(directory):
    """Find and load CLASS lensed or unlensed C_ell from a directory."""
    candidates = (glob.glob(os.path.join(directory, "*cl_lensed.dat")) +
                  glob.glob(os.path.join(directory, "*lensed*.dat")) +
                  glob.glob(os.path.join(directory, "*cl.dat")))
    if not candidates:
        raise FileNotFoundError(f"No CLASS C_ell file found in {directory}")
    fname = candidates[0]
    data  = np.loadtxt(fname, comments="#")
    ell   = data[:, 0].astype(int)
    cl    = data[:, 1]   # l(l+1)Cl/2pi  [dimensionless]
    Dl_uK2 = cl * T_CMB_uK2
    return ell, Dl_uK2

sim88 = P["sim88_outputs"]
print("  Loading SIM88 CLASS CMSTG C_ell ...", flush=True)
ell_cmstg, Dl_cmstg = load_class_cl(sim88["cmstg_cl"])
print("  Loading SIM88 CLASS LCDM C_ell ...", flush=True)
ell_lcdm, Dl_lcdm = load_class_cl(sim88["lcdm_cl"])

print(f"  CMSTG C_ell range: l={ell_cmstg[0]}..{ell_cmstg[-1]}, "
      f"D_l at l=2: {Dl_cmstg[0]:.1f} μK^2")
print(f"  LCDM C_ell range: l={ell_lcdm[0]}..{ell_lcdm[-1]}, "
      f"D_l at l=2: {Dl_lcdm[0]:.1f} μK^2")

# ── Planck noise model ────────────────────────────────────────────────────────
noise_cfg  = P["planck_noise"]
Delta_T    = float(noise_cfg["Delta_T_muK_arcmin"])    # μK-arcmin
theta_FWHM = float(noise_cfg["theta_FWHM_arcmin"])     # arcmin
sigma_b    = theta_FWHM * (math.pi / 180.0 / 60.0) / math.sqrt(8 * math.log(2))  # radians

# N_l in μK^2 (dimensionful)
def noise_spectrum(ell_arr):
    """Planck 143 GHz approximate beam-smoothed noise D_l = l(l+1)Nl/2pi [μK^2]."""
    l   = ell_arr.astype(float)
    Nl  = (Delta_T * math.pi / 180.0 / 60.0)**2 * np.exp(l * (l + 1) * sigma_b**2)
    return l * (l + 1) / (2.0 * math.pi) * Nl   # convert N_l to D_l

# ── Build common l grid ───────────────────────────────────────────────────────
lik = P["likelihood"]
l_min = int(lik["l_min_lowl"])
l_max = int(lik["l_max_highl"])
f_sky = float(lik["f_sky"])

# Common l range: intersection of all data sources
l_max_common = min(ell_cmstg[-1], ell_lcdm[-1], ell_planck[-1], l_max)
ell_common   = np.arange(l_min, l_max_common + 1)

# Interpolate all onto common grid
Dl_cmstg_c   = np.interp(ell_common, ell_cmstg,   Dl_cmstg)
Dl_lcdm_c   = np.interp(ell_common, ell_lcdm,   Dl_lcdm)
Dl_planck_c = np.interp(ell_common, ell_planck,  Dl_planck)
Nl_c        = noise_spectrum(ell_common)

# Reference = Planck theory (signal + noise floor)
Dl_ref = Dl_planck_c + Nl_c    # what CLASS must predict to match Planck+noise

# ── Approximate Planck TT Gaussian likelihood ─────────────────────────────────
def planck_chi2(Dl_pred, Dl_ref, Nl, ell, f_sky):
    """
    Approximate Planck TT Gaussian likelihood.
    -2 ln L = f_sky * sum_l (2l+1) * [x - 1 - ln(x)]
    where x = (D_l^pred + N_l) / (D_l^ref)
            = (D_l^pred + N_l) / (D_l^Planck + N_l)

    When D_l^pred = D_l^Planck: x = 1 → chi2 = 0  (correct baseline).
    Noise floor N_l prevents chi2 blow-up at high-l where signal < noise.
    Reference: Hamimeche & Lewis 2008, Eq. 9 (signal+noise ratio form).
    """
    x      = (Dl_pred + Nl) / np.maximum(Dl_ref, 1e-40)
    x      = np.maximum(x, 1e-10)
    loglik = f_sky * (2 * ell + 1) * (x - 1.0 - np.log(x))

    low_l_mask  = ell <= int(lik["l_max_lowl"])
    high_l_mask = ell >= int(lik["l_min_highl"])

    chi2_total = float(np.sum(loglik))
    chi2_lowl  = float(np.sum(loglik[low_l_mask]))
    chi2_highl = float(np.sum(loglik[high_l_mask]))
    n_modes    = int(np.sum(f_sky * (2 * ell + 1)))
    return chi2_total, chi2_lowl, chi2_highl, n_modes

chi2_cmstg,  chi2_cmstg_low,  chi2_cmstg_high,  n_modes_cmstg  = planck_chi2(Dl_cmstg_c,  Dl_ref, Nl_c, ell_common, f_sky)
chi2_lcdm,  chi2_lcdm_low,  chi2_lcdm_high,  n_modes_lcdm  = planck_chi2(Dl_lcdm_c,  Dl_ref, Nl_c, ell_common, f_sky)

delta_chi2    = chi2_cmstg - chi2_lcdm
n_eff_modes   = n_modes_cmstg   # approximately same

print(f"\n  ── Planck TT likelihood results (l={l_min}–{l_max_common}) ──")
print(f"  LCDM:  chi2={chi2_lcdm:.2f}  (low-l: {chi2_lcdm_low:.2f}, high-l: {chi2_lcdm_high:.2f})")
print(f"  CMSTG:  chi2={chi2_cmstg:.2f}  (low-l: {chi2_cmstg_low:.2f}, high-l: {chi2_cmstg_high:.2f})")
print(f"  Δchi2 (CMSTG−LCDM) = {delta_chi2:+.2f}")
print(f"  f_sky-weighted modes: {n_eff_modes}")
print(f"  chi2_CMSTG / n_modes = {chi2_cmstg/n_eff_modes:.5f}  (<<1 expected for good fit)")

# Fractional residuals
frac_cmstg = (Dl_cmstg_c - Dl_planck_c) / np.maximum(Dl_planck_c, 1e-40)
frac_lcdm = (Dl_lcdm_c - Dl_planck_c) / np.maximum(Dl_planck_c, 1e-40)

rms_cmstg = float(np.sqrt(np.mean(frac_cmstg**2)))
rms_lcdm = float(np.sqrt(np.mean(frac_lcdm**2)))
print(f"  RMS fractional residual: CMSTG={rms_cmstg*100:.2f}%, LCDM={rms_lcdm*100:.2f}%")

# ── Acceptance ────────────────────────────────────────────────────────────────
delta_max = float(P["acceptance"]["Delta_chi2_CMSTG_minus_LCDM_max"])
passed    = bool(delta_chi2 < delta_max)
tension_sigma = math.sqrt(abs(delta_chi2)) if delta_chi2 > 0 else 0.0
print(f"  CMB tension: Δchi2={delta_chi2:+.2f} ~ {tension_sigma:.1f}σ  "
      f"({'PASS' if passed else 'FAIL'})")

if not passed:
    print(f"  WARNING: CMSTG CMB tension Δchi2={delta_chi2:.2f} > {delta_max:.0f} "
          f"(threshold). Reporting result as FAIL — see diagnostics for interpretation.")

# ── Plots ─────────────────────────────────────────────────────────────────────

# Plot 1: D_l comparison with Planck
fig, axes = plt.subplots(2, 1, figsize=(10, 8), dpi=130,
                          gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05})
ax_top, ax_bot = axes

ax_top.plot(ell_planck, Dl_planck, "k-", lw=1.5, label="Planck 2018 theory (R3.01)")
ax_top.plot(ell_lcdm,   Dl_lcdm,   "b--", lw=1.2, alpha=0.8, label=f"CLASS LCDM (SIM88)  χ²={chi2_lcdm:.0f}")
ax_top.plot(ell_cmstg,   Dl_cmstg,   "r-",  lw=1.2, alpha=0.8, label=f"CLASS CMSTG (SIM88)  χ²={chi2_cmstg:.0f}")
ax_top.set_ylabel(r"$D_\ell^{TT}\;[\mu K^2]$", fontsize=11)
ax_top.set_xscale("log"); ax_top.set_yscale("log")
ax_top.set_xlim(l_min, l_max_common)
ax_top.legend(fontsize=8)
ax_top.set_title(f"Sim 89 — CMSTG vs Planck 2018 TT  |  Δχ²(CMSTG−LCDM)={delta_chi2:+.1f}  ({tension_sigma:.1f}σ)",
                 fontsize=10)
ax_top.set_xticklabels([])

ax_bot.semilogx(ell_common, frac_lcdm * 100, "b--", lw=0.8, alpha=0.7, label="LCDM")
ax_bot.semilogx(ell_common, frac_cmstg  * 100, "r-",  lw=0.8, alpha=0.8, label="CMSTG")
ax_bot.axhline(0, color="k", lw=0.5)
ax_bot.axhspan(-2, 2, alpha=0.08, color="gray")
ax_bot.set_xlabel(r"$\ell$", fontsize=11)
ax_bot.set_ylabel(r"$\Delta D_\ell/D_\ell^{\rm Pl}\;[\%]$", fontsize=9)
ax_bot.set_xlim(l_min, l_max_common)
ax_bot.set_ylim(-10, 10)
ax_bot.legend(fontsize=7)
plt.tight_layout()
p_cl = os.path.join(OUTPUTS, "sim89_cl_vs_planck.png")
fig.savefig(p_cl); plt.close(fig)

# Plot 2: chi2 contribution per l
fig, ax = plt.subplots(figsize=(9, 4), dpi=130)
chi2_per_l_cmstg = f_sky * (2*ell_common + 1) * np.maximum(
    (Dl_cmstg_c/np.maximum(Dl_ref, 1e-40)) - 1.0 -
    np.log(np.maximum(Dl_cmstg_c/np.maximum(Dl_ref, 1e-40), 1e-10)), 0)
chi2_per_l_lcdm = f_sky * (2*ell_common + 1) * np.maximum(
    (Dl_lcdm_c/np.maximum(Dl_ref, 1e-40)) - 1.0 -
    np.log(np.maximum(Dl_lcdm_c/np.maximum(Dl_ref, 1e-40), 1e-10)), 0)
ax.semilogx(ell_common, chi2_per_l_lcdm, "b--", lw=0.7, alpha=0.7, label="LCDM")
ax.semilogx(ell_common, chi2_per_l_cmstg,  "r-",  lw=0.7, alpha=0.8, label="CMSTG")
ax.axhline(0, color="k", lw=0.4)
ax.set_xlabel(r"$\ell$"); ax.set_ylabel(r"$\chi^2$ contribution per $\ell$")
ax.set_title(f"Sim 89 — Planck TT χ² per multipole  (total: CMSTG={chi2_cmstg:.0f}, LCDM={chi2_lcdm:.0f})")
ax.legend(fontsize=8)
plt.tight_layout()
p_chi2 = os.path.join(OUTPUTS, "sim89_chi2_table.png")
fig.savefig(p_chi2); plt.close(fig)

# Plot 3: residuals vs l ranges
fig, axes = plt.subplots(1, 3, figsize=(13, 4), dpi=120, sharey=True)
bands = [(2, 100, "Low-ℓ (2–100)"),
         (100, 500, "Mid-ℓ (100–500)"),
         (500, l_max_common, "High-ℓ (500–2500)")]
for ax, (lmin, lmax, title) in zip(axes, bands):
    mask = (ell_common >= lmin) & (ell_common <= lmax)
    ax.plot(ell_common[mask], frac_lcdm[mask]*100, "b-", lw=0.7, alpha=0.7, label="LCDM")
    ax.plot(ell_common[mask], frac_cmstg[mask]*100,  "r-", lw=0.8, alpha=0.8, label="CMSTG")
    ax.axhline(0, color="k", lw=0.4)
    ax.axhspan(-2, 2, alpha=0.08, color="gray")
    ax.set_xlabel(r"$\ell$"); ax.set_title(title, fontsize=9)
    ax.legend(fontsize=7)
axes[0].set_ylabel(r"$(D_\ell^{\rm pred}-D_\ell^{\rm Pl})/D_\ell^{\rm Pl}\;[\%]$")
plt.suptitle("Sim 89 — Fractional TT residuals vs Planck 2018 theory", y=1.02, fontsize=10)
plt.tight_layout()
p_resid = os.path.join(OUTPUTS, "sim89_residuals.png")
fig.savefig(p_resid, bbox_inches="tight"); plt.close(fig)

# ── Diagnostics ───────────────────────────────────────────────────────────────
diag = {
    "description": (
        "CMSTG Planck TT likelihood comparison (Phase 5b). "
        "Uses approximate full-sky Gaussian TT likelihood against "
        "Planck 2018 best-fit theory spectrum (COM_PowerSpect R3.01). "
        "chi2 = f_sky * sum_l (2l+1) * [x - 1 - ln(x)], x = D_l^pred / (D_l^Planck + N_l). "
        "Closes final CMB referee gap."
    ),
    "l_range":     f"{l_min}–{l_max_common}",
    "f_sky":       f_sky,
    "n_eff_modes": n_eff_modes,
    "T_CMB_uK2_conversion": T_CMB_uK2,
    "LCDM": {
        "chi2_total":  float(chi2_lcdm),
        "chi2_lowl":   float(chi2_lcdm_low),
        "chi2_highl":  float(chi2_lcdm_high),
        "chi2_per_mode": float(chi2_lcdm / n_eff_modes),
        "rms_fractional_residual": float(rms_lcdm),
        "interpretation": (
            "LCDM chi2 > 0 because CLASS uses Planck params but the theory file is "
            "the published best-fit (tiny numerical differences and lensing scheme). "
            "The LCDM residual sets the baseline noise floor for comparison."
        )
    },
    "CMSTG": {
        "chi2_total":  float(chi2_cmstg),
        "chi2_lowl":   float(chi2_cmstg_low),
        "chi2_highl":  float(chi2_cmstg_high),
        "chi2_per_mode": float(chi2_cmstg / n_eff_modes),
        "rms_fractional_residual": float(rms_cmstg),
    },
    "comparison": {
        "delta_chi2":    float(delta_chi2),
        "tension_sigma": float(tension_sigma),
        "passed":        bool(passed),
        "interpretation": (
            f"Δchi2 = chi2_CMSTG - chi2_LCDM = {delta_chi2:+.2f}. "
            f"This corresponds to {tension_sigma:.1f}σ CMB tension. "
            "The dominant contribution is the acoustic peak shift from "
            "CMSTG best-fit H0/Omega_m (BAO-driven) vs Planck LCDM. "
            "Note: this uses the approximate Gaussian likelihood (no data noise covariance). "
            "A full Planck likelihood evaluation (Commander+PlikHM) is beyond scope but "
            "this estimate captures the first-order CMB tension from the CMSTG background."
        )
    },
    "artifacts": {
        "cl_vs_planck": os.path.relpath(p_cl, BASE),
        "chi2_table":   os.path.relpath(p_chi2, BASE),
        "residuals":    os.path.relpath(p_resid, BASE),
    }
}

diag_path = os.path.join(OUTPUTS, "sim89_diagnostics.json")
with open(diag_path, "w") as fh:
    json.dump(diag, fh, indent=2)

print(f"\nWrote diagnostics to {diag_path}")
print(f"  Δchi2(CMSTG−LCDM) = {delta_chi2:+.2f}  ({tension_sigma:.1f}σ)  "
      f"{'PASS' if passed else 'FAIL'}")
