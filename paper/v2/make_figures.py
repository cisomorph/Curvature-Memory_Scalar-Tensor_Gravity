"""
make_figures.py — Publication figures for rift_paper_v2.tex

Generates three PDF figures:
  fig1_bao_residuals.pdf   — BAO data vs RIFT best-fit model (SIM87)
  fig2_cmb_Cl.pdf          — CMB Cℓ comparison RIFT vs ΛCDM (SIM88)
  fig3_lambda0_scan.pdf    — Λ₀ sensitivity scan: growth factor & G_eff (SIM91)

Reads data from simulation outputs; does not re-run any simulation.
"""

import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

# ── Paths ────────────────────────────────────────────────────────────────────
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..", "Ordered_Simulations")

SIM87_DIAG = os.path.join(ROOT, "SIM87/Outputs/sim87_diagnostics.json")
SIM88_RIFT = os.path.join(ROOT, "SIM88/Outputs/class_rift/00_00_cl.dat")
SIM88_LCDM = os.path.join(ROOT, "SIM88/Outputs/class_lcdm/00_00_cl.dat")
SIM91_DIAG = os.path.join(ROOT, "SIM91/Outputs/sim91_diagnostics.json")

# ── Style ────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":      "serif",
    "font.size":        11,
    "axes.labelsize":   12,
    "axes.titlesize":   12,
    "legend.fontsize":  10,
    "xtick.direction":  "in",
    "ytick.direction":  "in",
    "xtick.top":        True,
    "ytick.right":      True,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "figure.dpi":       150,
    "savefig.dpi":      300,
    "savefig.bbox":     "tight",
})

COLORS = {
    "rift":  "#1f77b4",   # blue
    "lcdm":  "#d62728",   # red
    "data":  "k",
}

# ═══════════════════════════════════════════════════════════════════════════
# Figure 1 — BAO residuals
# ═══════════════════════════════════════════════════════════════════════════

def compute_rift_predictions(H0, Omega_m, rd, Lambda0=0.003,
                              Omega_r=9.2e-5, Omega_b=0.049):
    """
    Compute BAO observables for the 12-point data vector.
    D_H = c/H(z), D_M = comoving angular diameter distance,
    D_V = (z * D_H * D_M^2)^{1/3}.
    All divided by r_d.
    """
    from scipy.integrate import quad
    c_kms = 299792.458  # km/s

    Psi0   = 0.01
    lam    = Lambda0 * Psi0**2
    Geff   = 1.0 / (1.0 + 16.0 * np.pi * lam)

    def E(z):
        a = 1.0 / (1.0 + z)
        OmL = 1.0 - Omega_m - Omega_r
        E2  = Omega_m / a**3 + Omega_r / a**4 + OmL
        return np.sqrt(max(E2, 1e-30))

    def DH(z):
        return c_kms / (H0 * E(z))

    def DM(z):
        integrand = lambda zp: 1.0 / (H0 * E(zp))
        val, _ = quad(integrand, 0, z)
        return c_kms * val

    def DV(z):
        return (z * DH(z) * DM(z)**2)**(1.0/3.0)

    BAO_DATA = [
        (0.15, "DV_over_rd",  4.47,    0.17,    "6dFGS+MGS"),
        (0.38, "DH_over_rd", 25.00,    0.76,    "BOSS DR12"),
        (0.38, "DM_over_rd", 10.23,    0.17,    "BOSS DR12"),
        (0.51, "DH_over_rd", 22.33,    0.58,    "BOSS DR12"),
        (0.51, "DM_over_rd", 13.36,    0.21,    "BOSS DR12"),
        (0.70, "DH_over_rd", 19.33,    0.53,    "eBOSS DR16 LRG"),
        (0.70, "DM_over_rd", 17.86,    0.33,    "eBOSS DR16 LRG"),
        (0.85, "DV_over_rd", 18.33,    0.595,   "eBOSS DR16 QSO"),
        (1.48, "DH_over_rd", 13.26,    0.55,    "eBOSS DR16 QSO"),
        (1.48, "DM_over_rd", 30.69,    0.80,    "eBOSS DR16 QSO"),
        (2.33, "DH_over_rd",  8.9906,  0.2161,  "Ly-α DR16"),
        (2.33, "DM_over_rd", 37.4334,  1.2669,  "Ly-α DR16"),
    ]

    preds = []
    for z, qty, obs, sig, survey in BAO_DATA:
        if qty == "DH_over_rd":
            pred = DH(z) / rd
        elif qty == "DM_over_rd":
            pred = DM(z) / rd
        else:
            pred = DV(z) / rd
        preds.append((z, qty, obs, sig, survey, pred))
    return preds


def make_fig1():
    diag = json.load(open(SIM87_DIAG))
    bf   = diag["best_fit_full_cov"]["params"]
    H0   = bf["H0"]
    Om   = bf["Omega_m"]
    rd   = bf["r_d"]
    L0   = bf["Lambda0"]

    preds = compute_rift_predictions(H0, Om, rd, L0)

    # Separate by type
    dh = [(z, obs, sig, pred) for z, qty, obs, sig, _, pred in preds if "DH" in qty]
    dm = [(z, obs, sig, pred) for z, qty, obs, sig, _, pred in preds if "DM" in qty]
    dv = [(z, obs, sig, pred) for z, qty, obs, sig, _, pred in preds if "DV" in qty]

    fig, axes = plt.subplots(3, 1, figsize=(6, 7), sharex=False)
    fig.subplots_adjust(hspace=0.45)

    def plot_residuals(ax, data, ylabel, title):
        zs   = np.array([d[0] for d in data])
        obs  = np.array([d[1] for d in data])
        sigs = np.array([d[2] for d in data])
        pred = np.array([d[3] for d in data])
        pull = (obs - pred) / sigs

        ax.axhline(0, color="gray", lw=0.8, ls="--")
        ax.axhspan(-1, 1, color="gray", alpha=0.12, label=r"$\pm1\sigma$")
        ax.axhspan(-2, 2, color="gray", alpha=0.06, label=r"$\pm2\sigma$")
        ax.errorbar(zs, pull, yerr=1.0, fmt="o", color=COLORS["data"],
                    ms=5, capsize=3, lw=1.2, label="Data pull")
        ax.set_ylabel(r"$(d_\mathrm{obs} - d_\mathrm{pred})/\sigma$", fontsize=10)
        ax.set_title(title, fontsize=10, pad=4)
        ax.set_ylim(-3.5, 3.5)
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))

    plot_residuals(axes[0], dh, r"$D_H/r_d$", r"$D_H/r_d$ residuals")
    plot_residuals(axes[1], dm, r"$D_M/r_d$", r"$D_M/r_d$ residuals")
    plot_residuals(axes[2], dv, r"$D_V/r_d$", r"$D_V/r_d$ residuals")

    for ax, data in zip(axes, [dh, dm, dv]):
        zs = [d[0] for d in data]
        ax.set_xlim(min(zs) - 0.15, max(zs) + 0.15)
        ax.set_xlabel(r"Redshift $z$", fontsize=10)

    # Add legend only to top panel
    axes[0].legend(loc="upper right", framealpha=0.8, fontsize=9)

    # Annotate chi2
    chi2 = diag["best_fit_full_cov"]["chi2"]
    dof  = diag["n_dof"]
    fig.text(0.98, 0.01,
             rf"$\chi^2/\mathrm{{dof}} = {chi2:.2f}/{dof} = {chi2/dof:.2f}$"
             rf"  ($H_0={H0:.1f}$, $\Omega_m={Om:.3f}$, $\Lambda_0={L0:.3f}$)",
             ha="right", va="bottom", fontsize=8.5, color="0.4")

    fig.suptitle("RIFT BAO Full-Covariance Fit (Sim 87)", fontsize=12, y=1.01)
    out = os.path.join(HERE, "fig1_bao_residuals.pdf")
    fig.savefig(out)
    plt.close(fig)
    print(f"Saved {out}")


# ═══════════════════════════════════════════════════════════════════════════
# Figure 2 — CMB Cℓ comparison
# ═══════════════════════════════════════════════════════════════════════════

def make_fig2():
    # CLASS output columns: ell, TT (dimensionless C_ell), EE, TE
    # Convert to D_ell = ell*(ell+1)*C_ell/(2*pi)  [dimensionless, then
    # multiply by T_CMB^2 in μK^2 for conventional units]
    T_CMB_uK = 2.7255e6   # μK

    def load_cl(path):
        data = np.loadtxt(path)
        ell  = data[:, 0].astype(int)
        cl_tt = data[:, 1]   # dimensionless C_ell^TT from CLASS
        dl_tt = ell * (ell + 1) * cl_tt / (2.0 * np.pi) * T_CMB_uK**2
        return ell, dl_tt

    ell_r, dl_r = load_cl(SIM88_RIFT)
    ell_l, dl_l = load_cl(SIM88_LCDM)

    # Align on common ell range
    l_min = max(ell_r[0], ell_l[0])
    l_max = min(ell_r[-1], ell_l[-1])
    mask_r = (ell_r >= l_min) & (ell_r <= l_max)
    mask_l = (ell_l >= l_min) & (ell_l <= l_max)
    ell   = ell_r[mask_r]
    dl_rift = dl_r[mask_r]
    dl_lcdm = dl_l[mask_l]

    # Fractional difference
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = (dl_rift - dl_lcdm) / np.abs(dl_lcdm)

    # Smooth for display (running mean, window=20)
    def running_mean(x, w=20):
        return np.convolve(x, np.ones(w)/w, mode="same")

    frac_smooth = running_mean(frac, 20)

    diag = json.load(open(os.path.join(ROOT, "SIM88/Outputs/sim88_diagnostics.json")))
    rms  = diag["class_results"]["rms_dCl_over_Cl"]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 6),
                                   gridspec_kw={"height_ratios": [3, 1.4]})
    fig.subplots_adjust(hspace=0.1)

    # Top: power spectra
    ax1.plot(ell, dl_lcdm, color=COLORS["lcdm"],  lw=1.2, label=r"$\Lambda$CDM (Planck 2018)")
    ax1.plot(ell, dl_rift,  color=COLORS["rift"], lw=1.2, ls="--",
             label=rf"RIFT ($\Lambda_0=0.003$, $H_0=68.1$, $\Omega_m=0.294$)")
    ax1.set_xscale("log")
    ax1.set_ylabel(r"$D_\ell^{TT}\ [\mu\mathrm{K}^2]$")
    ax1.set_xlim(2, 1500)
    ax1.set_ylim(0, None)
    ax1.xaxis.set_ticklabels([])
    ax1.legend(loc="upper right", framealpha=0.85)
    ax1.set_title("CMB TT Power Spectrum: RIFT vs $\\Lambda$CDM (Sim 88)", pad=6)

    # Bottom: fractional difference
    ax2.axhline(0, color="gray", lw=0.8, ls="--")
    ax2.fill_between(ell, frac * 100, 0, color=COLORS["rift"], alpha=0.25)
    ax2.plot(ell, frac_smooth * 100, color=COLORS["rift"], lw=1.0)
    ax2.axhline( rms * 100, color="0.5", lw=0.8, ls=":", label=rf"RMS $= {rms*100:.2f}\%$")
    ax2.axhline(-rms * 100, color="0.5", lw=0.8, ls=":")
    ax2.set_xscale("log")
    ax2.set_xlabel(r"Multipole $\ell$")
    ax2.set_ylabel(r"$\Delta D_\ell/D_\ell^{\Lambda\mathrm{CDM}}\ [\%]$", fontsize=10)
    ax2.set_xlim(2, 1500)
    ax2.legend(loc="upper right", framealpha=0.85, fontsize=9)
    ax2.yaxis.set_minor_locator(AutoMinorLocator(2))

    out = os.path.join(HERE, "fig2_cmb_Cl.pdf")
    fig.savefig(out)
    plt.close(fig)
    print(f"Saved {out}")


# ═══════════════════════════════════════════════════════════════════════════
# Figure 3 — Λ₀ sensitivity scan
# ═══════════════════════════════════════════════════════════════════════════

def make_fig3():
    diag    = json.load(open(SIM91_DIAG))
    scan    = diag["scan_results"]
    thresh  = diag["detection_thresholds"]

    L0     = np.array([s["Lambda0"]        for s in scan])
    Gdev   = np.array([s["Geff_deviation"] for s in scan]) * 1e6   # → ppm
    Ddev   = np.array([s["D_deviation"]    for s in scan]) * 100   # → %
    chi2   = np.array([s["chi2_BAO"]       for s in scan])

    # Thresholds (keys use descriptive names in SIM91 output)
    G_thresh_ppm = 1000.0   # 0.1% = 1000 ppm
    D_thresh_pct = 0.1      # 0.1%

    fig, axes = plt.subplots(1, 3, figsize=(10, 3.8))
    fig.subplots_adjust(wspace=0.38)

    # Panel 1: G_eff deviation
    ax = axes[0]
    ax.plot(L0, Gdev, "o-", color=COLORS["rift"], ms=4, lw=1.5)
    ax.axhline(G_thresh_ppm, color="orange", lw=1.0, ls="--",
               label=r"$0.1\%$ threshold (1000 ppm)")
    ax.axvline(0.003, color="gray", lw=0.8, ls=":", label=r"BAO best-fit $\Lambda_0$")
    ax.set_xlabel(r"$\Lambda_0$")
    ax.set_ylabel(r"$|G_\mathrm{eff}/G - 1|\ [\mathrm{ppm}]$")
    ax.set_title(r"Effective Newton Constant", fontsize=10)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend(fontsize=8, framealpha=0.85)
    ax.xaxis.set_minor_locator(AutoMinorLocator(3))

    # Panel 2: growth factor deviation
    ax = axes[1]
    ax.plot(L0, Ddev, "s-", color="#2ca02c", ms=4, lw=1.5)
    ax.axhline(D_thresh_pct, color="orange", lw=1.0, ls="--",
               label=rf"$0.1\%$ detection threshold")
    ax.axvline(0.003, color="gray", lw=0.8, ls=":", label=r"BAO best-fit $\Lambda_0$")
    ax.axvline(0.05,  color="purple", lw=0.8, ls="-.", label=r"$\Lambda_0=0.05$ (structure)")
    ax.set_xlabel(r"$\Lambda_0$")
    ax.set_ylabel(r"$|D_\mathrm{RIFT}/D_{\Lambda\mathrm{CDM}} - 1|\ [\%]$")
    ax.set_title(r"Linear Growth Factor", fontsize=10)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend(fontsize=8, framealpha=0.85)

    # Panel 3: BAO chi2
    ax = axes[2]
    chi2_ref = chi2[0]   # Lambda0=0 value
    ax.plot(L0, chi2 - chi2_ref, "^-", color=COLORS["lcdm"], ms=4, lw=1.5)
    ax.axhline(0, color="gray", lw=0.8, ls="--")
    ax.axvline(0.003, color="gray", lw=0.8, ls=":", label=r"BAO best-fit $\Lambda_0$")
    ax.set_xlabel(r"$\Lambda_0$")
    ax.set_ylabel(r"$\Delta\chi^2_\mathrm{BAO}$ (relative to $\Lambda_0=0$)")
    ax.set_title(r"BAO $\chi^2$ Sensitivity", fontsize=10)
    ax.set_xscale("log")
    ax.legend(fontsize=8, framealpha=0.85)
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    fig.suptitle(r"RIFT Coupling-Strength Sensitivity Scan (Sim 91)"
                 r"  $[H_0=67.6,\,\Omega_m=0.312]$",
                 fontsize=11, y=1.02)

    out = os.path.join(HERE, "fig3_lambda0_scan.pdf")
    fig.savefig(out)
    plt.close(fig)
    print(f"Saved {out}")


# ═══════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    make_fig1()
    make_fig2()
    make_fig3()
    print("All figures written to", HERE)
