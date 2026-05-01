#!/usr/bin/env python3
# Sim 80 v2 — Resummed recursive damping re-run
"""
Usage
-----
    python simulation_runner_cmstg.py

Reads:
    Inputs/sim80_params.json
Writes:
    Outputs/sim80_uv_falloff.png
    Outputs/sim80_bubble_cutoff.png
    Outputs/sim80_diagnostics.json  (UV slopes, bubble slopes, params)
Fails loudly if required inputs or outputs are missing.
"""

import os, json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")  # headless-safe
import matplotlib.pyplot as plt

# ---------- Paths & input validation ----------
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INPUTS = os.path.join(BASE, "Inputs")
OUTPUTS = os.path.join(BASE, "Outputs")
PARAMS = os.path.join(INPUTS, "sim80_params.json")

if not os.path.isdir(INPUTS):
    raise FileNotFoundError("Missing Inputs/ directory")
if not os.path.isfile(PARAMS):
    raise FileNotFoundError("Missing Inputs/sim80_params.json")

os.makedirs(OUTPUTS, exist_ok=True)

# ---------- Load parameters ----------
with open(PARAMS, "r") as f:
    P = json.load(f)

m          = float(P["m"])
g          = float(P["g"])
Lambda_mem = float(P["Lambda_mem"])
gamma      = float(P["gamma"])

p_min = float(P["p_grid"]["p_min"])
p_max = float(P["p_grid"]["p_max"])
Np    = int(P["p_grid"]["Np"])

p_ext      = float(P["bubble"]["p_ext"])
Lambda_min = float(P["bubble"]["Lambda_min"])
Lambda_max = float(P["bubble"]["Lambda_max"])
N_L        = int(P["bubble"]["N_L"])
Ntheta     = int(P["bubble"]["Ntheta"])
Nx_log     = int(P["bubble"]["Nx_log"])
q_min      = float(P["bubble"]["q_min"])

# ---------- Propagators ----------
def G0(p):
    return 1.0 / (p**2 + m**2)

def G_resum(p):
    # CMSTG (resummed) damping: for p >> Λ_mem, G ~ p^{-2(1+γ)}
    return 1.0/(p**2 + m**2) * (1.0 + (p**2)/(Lambda_mem**2))**(-gamma)

def tail_slope(x, y, frac_tail=0.33):
    """Fit log y = a + b log x on the last frac_tail of the grid; return b."""
    n0 = int((1.0 - frac_tail) * len(x))
    lx = np.log(x[n0:])
    ly = np.log(y[n0:])
    b, a = np.polyfit(lx, ly, 1)
    return b

# ---------- 1) UV falloff ----------
p = np.logspace(np.log10(p_min), np.log10(p_max), Np)
G0_abs = np.abs(G0(p))
Gr_abs = np.abs(G_resum(p))

plt.figure(figsize=(6, 4.2), dpi=140)
plt.loglog(p, G0_abs, label="Baseline $G_0(p)$")
plt.loglog(p, Gr_abs, label=fr"CMSTG $G(p)$, $\gamma={gamma}$")
plt.xlabel("p (units of m)")
plt.ylabel("|G(p)|")
plt.title("Sim 80 v2 — UV falloff (log–log)")
plt.legend()
plt.tight_layout()
uv_path = os.path.join(OUTPUTS, "sim80_uv_falloff.png")
plt.savefig(uv_path)
plt.close()

s0 = tail_slope(p, G0_abs)
sr = tail_slope(p, Gr_abs)

# ---------- 2) Bubble integral on log-q grid ----------
# Angular GL on θ ∈ [0, π]
t_nodes, t_w = np.polynomial.legendre.leggauss(Ntheta)
theta  = 0.5*(t_nodes + 1.0) * math.pi
w_th   = 0.5 * math.pi * t_w
sin2   = np.sin(theta)**2
cos_th = np.cos(theta)

def bubble_I(Lc, use_cmstg: bool):
    # q ∈ (q_min, Lc], uniform in x = ln q
    x = np.linspace(np.log(q_min), np.log(Lc), Nx_log)
    q = np.exp(x)
    # Measure: q^3 dq = q^4 dx, integrated in x with uniform spacing
    Gq = G_resum(q) if use_cmstg else G0(q)
    q_col = q[:, None]
    k = np.sqrt(q_col**2 + p_ext**2 + 2.0*q_col*p_ext*cos_th[None, :])
    Gk = G_resum(k) if use_cmstg else G0(k)
    ang = (w_th[None, :] * sin2[None, :] * Gk).sum(axis=1)  # shape (Nx_log,)
    integrand_x = (q**4) * Gq * ang
    return (1.0/(8.0*math.pi**2)) * np.trapz(integrand_x, x)

L = np.logspace(np.log10(Lambda_min), np.log10(Lambda_max), N_L)
I0 = np.array([bubble_I(Lc, False) for Lc in L])
Ir = np.array([bubble_I(Lc, True)  for Lc in L])

plt.figure(figsize=(6, 4.2), dpi=140)
plt.plot(L, I0, label=r"Baseline $I_0(\Lambda)$")
plt.plot(L, Ir, label=rf"CMSTG $I_{{\rm CMSTG}}(\Lambda)$, $\gamma={gamma}$")
plt.xscale("log")
plt.xlabel(r"Cutoff $\Lambda$")
plt.ylabel("Bubble integral I(Λ)")
plt.title("Sim 80 v2 — One-loop bubble vs cutoff (log-q)")
plt.legend()
plt.tight_layout()
bubble_path = os.path.join(OUTPUTS, "sim80_bubble_cutoff.png")
plt.savefig(bubble_path)
plt.close()

# ---------- 3) Console summary ----------
print(f"Tail slopes: baseline ≈ {s0:.3f}, CMSTG ≈ {sr:.3f}")
print("Wrote figures to Outputs/.")

# ---------- 4) Acceptance checks + diagnostics ----------
def slope_vs_logL(L, I, tail_frac=0.20):
    """Linear slope dI / d(ln Λ) over the last tail_frac of points."""
    n0 = max(0, int((1.0 - tail_frac) * len(L)))
    x = np.log(L[n0:])
    y = I[n0:]
    b, a = np.polyfit(x, y, 1)
    return b

# UV acceptance: CMSTG slope must be < -2 by a margin
uv_margin = 0.05
if not (sr < -2.0 - uv_margin):
    raise RuntimeError(
        f"UV acceptance failed: CMSTG slope {sr:.3f} is not steeper than -2 by > {uv_margin}"
    )

# Cutoff acceptance: CMSTG bubble growth must be flatter than baseline
s_base = slope_vs_logL(L, I0)
s_cmstg = slope_vs_logL(L, Ir)
if not (abs(s_cmstg) < 0.5 * abs(s_base)):
    raise RuntimeError(
        f"Cutoff acceptance failed: CMSTG bubble slope {s_cmstg:.4g} not sufficiently flatter than baseline {s_base:.4g}"
    )

diag = {
    "uv_tail_slopes": {"baseline": float(s0), "cmstg": float(sr)},
    "bubble_slope_logLambda": {"baseline": float(s_base), "cmstg": float(s_cmstg)},
    "params": {
        "m": m, "g": g, "Lambda_mem": Lambda_mem, "gamma": gamma,
        "p_ext": p_ext, "Lambda_min": Lambda_min, "Lambda_max": Lambda_max,
        "Ntheta": Ntheta, "Nx_log": Nx_log, "q_min": q_min
    },
    "artifacts": {
        "uv_plot": os.path.relpath(uv_path, BASE),
        "bubble_plot": os.path.relpath(bubble_path, BASE)
    }
}
with open(os.path.join(OUTPUTS, "sim80_diagnostics.json"), "w") as f:
    json.dump(diag, f, indent=2)

# Final file existence checks
for _p in [uv_path, bubble_path]:
    if not os.path.isfile(_p):
        raise FileNotFoundError(f"Expected output missing: {_p}")

print("Acceptance checks passed. Wrote Outputs/sim80_diagnostics.json")

