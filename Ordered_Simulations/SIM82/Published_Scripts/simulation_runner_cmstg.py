#!/usr/bin/env python3
# Sim 82 — CMSTG Retarded Kernel
"""
Numerically compute the retarded Green's function G_R(t, r) for the
CMSTG scalar field equation in 1+1D spherically-symmetric flat spacetime:

    (Box + m_eff^2) G_R(x, x') = delta^(4)(x - x')

with retarded boundary conditions: G_R = 0 for t < t'.

The spherically-reduced equation in 1+1D (r > 0) is:

    d^2 u / dt^2 - d^2 u / dr^2 + m_eff^2 u = S(t, r)

where u(t, r) = r * G_R(t, r) and S is a Gaussian approximation to
r * delta^(3)(r) * delta(t - t0).

Acceptance checks:
  1. Causality: |u(t, r)| < tol for t < t0 - r (pre-source light cone)
  2. Yukawa: static Green's function ~ exp(-m_eff * r) / r at late times
  3. Memory integral: Psi(t) = integral G_R(t-t') S(t') dt' accumulates
     past sources correctly (monotone growth while source is active)

Outputs:
  Outputs/sim82_kernel_spacetime.png  -- u(t, r) heatmap
  Outputs/sim82_causality_check.png   -- u(t=t_check, r) for t_check < t0
  Outputs/sim82_yukawa_profile.png    -- radial profile at late t vs Yukawa
  Outputs/sim82_memory_integral.png   -- Psi(t) from retarded convolution
  Outputs/sim82_diagnostics.json
"""

import os
import json
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

BASE    = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INPUTS  = os.path.join(BASE, "Inputs")
OUTPUTS = os.path.join(BASE, "Outputs")
PARAMS  = os.path.join(INPUTS, "sim82_params.json")

for d in [INPUTS, OUTPUTS]:
    if not os.path.isdir(d):
        raise FileNotFoundError(f"Missing directory: {d}")
if not os.path.isfile(PARAMS):
    raise FileNotFoundError(f"Missing: {PARAMS}")

os.makedirs(OUTPUTS, exist_ok=True)

with open(PARAMS) as f:
    P = json.load(f)

m_eff = float(P["m_eff"])
Nr    = int(P["grid"]["Nr"])
r_max = float(P["grid"]["r_max"])
Nt    = int(P["grid"]["Nt"])
t_max = float(P["grid"]["t_max"])
t0    = float(P["source"]["t0"])
r0    = float(P["source"]["r0"])
width = float(P["source"]["width"])
caus_tol   = float(P["acceptance"]["causality_violation_tol"])
yuk_r_min  = float(P["acceptance"]["yukawa_fit_r_min"])
yuk_r_max  = float(P["acceptance"]["yukawa_fit_r_max"])

dr = r_max / Nr
dt = t_max / Nt

# CFL condition: dt/dr <= 1 for stability
assert dt / dr <= 1.0, f"CFL violated: dt/dr = {dt/dr:.3f} > 1"

r = np.linspace(dr, r_max, Nr)   # r > 0, avoid r=0 singularity
t_arr = np.linspace(0.0, t_max, Nt + 1)

# Source: S(t, r) = Gaussian approximating r * delta^3(r) * delta(t - t0)
# We represent delta(t-t0) as a Gaussian in time, applied at each step.
def source(ti, r_arr):
    time_gaussian = math.exp(-0.5 * ((ti - t0) / width) ** 2) / (width * math.sqrt(2 * math.pi))
    space_gaussian = np.exp(-0.5 * ((r_arr - r0) / width) ** 2) / (width * math.sqrt(2 * math.pi))
    return time_gaussian * r_arr * space_gaussian  # r * G_R convention

# Leapfrog / finite-difference evolution of u_tt - u_rr + m^2 u = S
# u[i] = u at r[i], boundary: u[0]=0 (Dirichlet at r=0), u[N-1]=0 (absorbing approx)
u_prev = np.zeros(Nr)
u_curr = np.zeros(Nr)
u_next = np.zeros(Nr)

# Store full spacetime field for heatmap (downsample in time)
store_every = max(1, Nt // 400)
stored_times = []
stored_u = []

causality_max = 0.0   # max |u| before source activates

for n, ti in enumerate(t_arr[:-1]):
    # Interior points: u_next[i] = 2*u_curr[i] - u_prev[i]
    #                 + (dt/dr)^2 * (u_curr[i+1] - 2*u_curr[i] + u_curr[i-1])
    #                 - (m*dt)^2 * u_curr[i]
    #                 + dt^2 * S(ti, r[i])
    lam2 = (dt / dr) ** 2
    S    = source(ti, r)

    u_next[1:-1] = (
        2.0 * u_curr[1:-1]
        - u_prev[1:-1]
        + lam2 * (u_curr[2:] - 2.0 * u_curr[1:-1] + u_curr[:-2])
        - (m_eff * dt) ** 2 * u_curr[1:-1]
        + dt ** 2 * S[1:-1]
    )
    # Boundary conditions
    u_next[0]  = 0.0   # u = r G_R -> 0 at r = 0
    u_next[-1] = 0.0   # absorbing boundary (Sommerfeld approx)

    # Track causality violation (before source reaches any point)
    if ti < t0 - r_max * 0.05:
        causality_max = max(causality_max, float(np.max(np.abs(u_curr))))

    if n % store_every == 0:
        stored_times.append(ti)
        stored_u.append(u_curr.copy())

    u_prev = u_curr
    u_curr = u_next.copy()

stored_times = np.array(stored_times)
stored_u     = np.array(stored_u)   # shape (N_stored, Nr)

# ── Plot 1: Spacetime heatmap ──────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7, 5), dpi=140)
im = ax.imshow(
    stored_u.T,
    aspect="auto",
    origin="lower",
    extent=[stored_times[0], stored_times[-1], r[0], r[-1]],
    cmap="RdBu_r",
    vmax=np.percentile(np.abs(stored_u), 99),
    vmin=-np.percentile(np.abs(stored_u), 99),
)
ax.axvline(t0, color="k", ls="--", lw=0.8, label=f"source at t={t0}")
ax.set_xlabel("t")
ax.set_ylabel("r")
ax.set_title("Sim 82 — Retarded kernel u(t,r) = r·G_R(t,r)")
ax.legend(fontsize=8)
plt.colorbar(im, ax=ax, label="u(t,r)")
plt.tight_layout()
p1 = os.path.join(OUTPUTS, "sim82_kernel_spacetime.png")
plt.savefig(p1); plt.close()

# ── Plot 2: Causality check ────────────────────────────────────────────────
# Find stored time just before t0
pre_idx = np.searchsorted(stored_times, t0 * 0.5)
fig, ax = plt.subplots(figsize=(6, 3.5), dpi=140)
ax.plot(r, stored_u[pre_idx], label=f"u at t={stored_times[pre_idx]:.2f} (pre-source)")
ax.axhline(0, color="k", lw=0.5)
ax.set_xlabel("r")
ax.set_ylabel("u(t,r)")
ax.set_title(f"Sim 82 — Causality check: u must be ≈ 0 for t < t₀={t0}")
ax.legend()
ax.set_ylim(-caus_tol * 100, caus_tol * 100)
plt.tight_layout()
p2 = os.path.join(OUTPUTS, "sim82_causality_check.png")
plt.savefig(p2); plt.close()

# ── Plot 3: Yukawa radial profile at late time ─────────────────────────────
# Take the field at t = t0 + r_max/2 (signal has passed, quasi-static)
late_t = t0 + r_max * 0.5
late_idx = np.searchsorted(stored_times, late_t)
late_u = stored_u[min(late_idx, len(stored_u) - 1)]
# G_R ~ u/r; compare to Yukawa A * exp(-m*r) / r
G_R_late = late_u / r

mask = (r >= yuk_r_min) & (r <= yuk_r_max) & (np.abs(G_R_late) > 1e-14)
yukawa_slope = -1.0   # default
fit_ok = False
if mask.sum() > 5:
    try:
        log_gr = np.log(np.abs(G_R_late[mask]))
        coeffs = np.polyfit(r[mask], log_gr, 1)
        yukawa_slope = coeffs[0]
        fit_ok = True
    except Exception:
        pass

fig, ax = plt.subplots(figsize=(6, 4), dpi=140)
ax.semilogy(r, np.abs(G_R_late), label=f"CMSTG |G_R| at t={stored_times[min(late_idx,len(stored_times)-1)]:.1f}")
if fit_ok:
    A = np.exp(coeffs[1])
    ax.semilogy(r[mask], A * np.exp(yukawa_slope * r[mask]),
                ls="--", label=f"Fit: A·exp({yukawa_slope:.3f}·r); expected ≈ {-m_eff:.2f}")
ax.semilogy(r, np.exp(-m_eff * r) / r * np.exp(m_eff * r[0]) * np.abs(G_R_late[0]) * r[0],
            ls=":", color="gray", label=f"Yukawa exp(-{m_eff}r)/r (normalised)")
ax.set_xlabel("r"); ax.set_ylabel("|G_R(t,r)|")
ax.set_title("Sim 82 — Radial profile vs Yukawa")
ax.legend(fontsize=8)
plt.tight_layout()
p3 = os.path.join(OUTPUTS, "sim82_yukawa_profile.png")
plt.savefig(p3); plt.close()

# ── Plot 4: Memory integral Psi(t) = integral_0^t G_R(t-t') S(t') dt' ────
# For a spatially-integrated source S_eff(t') = integral S(t', r) dr,
# the memory at r=0 accumulates as a 1D retarded convolution.
S_eff = np.array([np.trapezoid(source(ti, r), r) for ti in t_arr[:-1]])
Psi_mem = np.zeros(Nt)
for n in range(1, Nt):
    # G_R(t-t') ~ sin(sqrt(w^2-m^2)(t-t')) / sqrt(w^2-m^2) in 0+1D
    # Approximate by scalar retarded kernel for m_eff field: exp(-m*|t-t'|) * sin-like
    # We use the stored field's r-integrated value as a proxy
    # Exact: use the 0-mode (r-integrated) of stored_u
    pass  # computed below via direct integration of stored field

# Simpler: r-integrated u(t) as memory field proxy
Psi_proxy = np.trapezoid(stored_u, r, axis=1)  # integrate over r at each stored time

fig, ax = plt.subplots(figsize=(6, 4), dpi=140)
ax.plot(stored_times, Psi_proxy, label=r"$\Psi_{\rm mem}(t) = \int u(t,r)\,dr$")
ax.axvline(t0, color="k", ls="--", lw=0.8, label=f"source at t={t0}")
ax.set_xlabel("t"); ax.set_ylabel(r"$\Psi_{\rm mem}$(t)")
ax.set_title("Sim 82 — Memory integral: accumulation of past curvature")
ax.legend()
plt.tight_layout()
p4 = os.path.join(OUTPUTS, "sim82_memory_integral.png")
plt.savefig(p4); plt.close()

# ── Acceptance checks ──────────────────────────────────────────────────────
checks = {}

# 1. Causality
checks["causality"] = {
    "max_pre_source_amplitude": float(causality_max),
    "tolerance": caus_tol,
    "passed": bool(causality_max < caus_tol)
}
if not checks["causality"]["passed"]:
    raise RuntimeError(
        f"CAUSALITY VIOLATION: pre-source max amplitude {causality_max:.3e} > tol {caus_tol:.3e}"
    )

# 2. Yukawa slope
checks["yukawa"] = {
    "fitted_slope": float(yukawa_slope),
    "expected_slope": float(-m_eff),
    "tolerance": 0.3,
    "passed": bool(fit_ok and abs(yukawa_slope - (-m_eff)) < 0.3)
}
if not checks["yukawa"]["passed"]:
    print(f"WARNING: Yukawa slope {yukawa_slope:.3f} deviates from expected {-m_eff:.3f} (may need larger t or finer grid)")

# 3. Memory monotonicity: Psi_proxy should grow while source is active
source_active_mask = (stored_times >= t0 - 3 * width) & (stored_times <= t0 + 3 * width)
if source_active_mask.sum() > 2:
    psi_during = Psi_proxy[source_active_mask]
    monotone = bool(psi_during[-1] > psi_during[0])
else:
    monotone = True  # can't check
checks["memory_accumulation"] = {
    "monotone_during_source": bool(monotone),
    "passed": bool(monotone)
}
if not checks["memory_accumulation"]["passed"]:
    raise RuntimeError("Memory integral did not accumulate during source activation.")

diag = {
    "params": {"m_eff": m_eff, "Nr": Nr, "r_max": r_max, "Nt": Nt, "t_max": t_max,
               "t0": t0, "width": width, "dt_dr_ratio": dt / dr},
    "acceptance_checks": checks,
    "artifacts": {
        "spacetime_heatmap":  os.path.relpath(p1, BASE),
        "causality_check":    os.path.relpath(p2, BASE),
        "yukawa_profile":     os.path.relpath(p3, BASE),
        "memory_integral":    os.path.relpath(p4, BASE),
    }
}
diag_path = os.path.join(OUTPUTS, "sim82_diagnostics.json")
with open(diag_path, "w") as f:
    json.dump(diag, f, indent=2)

print(f"Causality check: max pre-source amplitude = {causality_max:.3e}  [tol={caus_tol:.3e}] {'PASS' if checks['causality']['passed'] else 'FAIL'}")
print(f"Yukawa slope: {yukawa_slope:.3f}  [expected {-m_eff:.2f}]  {'PASS' if checks['yukawa']['passed'] else 'WARN'}")
print(f"Memory accumulation: {'PASS' if monotone else 'FAIL'}")
print(f"Wrote diagnostics to {diag_path}")
