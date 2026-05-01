#!/usr/bin/env python3
# Sim 84 — CMSTG Coupling Survey
"""
Survey Lambda(Psi) functional forms for the non-minimal curvature coupling.

The gravitational equation derived from the CMSTG action gives an effective
Newton constant:

    G_eff(Psi) = G / (1 + 16*pi*G * Lambda(Psi))

Two requirements must hold for all Psi:
  1. G_eff > 0: requires Lambda(Psi) > -1/(16*pi*G)
  2. GR recovery: Lambda(Psi) -> 0 as Psi -> 0

Additionally, for stability (no ghost graviton):
  3. 1/(8*pi*G) + 2*Lambda(Psi) > 0, i.e. Lambda(Psi) > -1/(16*pi*G)
     (same as condition 1)

For each functional form, this simulation:
  - Plots Lambda(Psi) and G_eff(Psi) vs Psi
  - Checks conditions 1-2
  - Computes the derivative Lambda'(Psi) (controls source term in scalar eq)
  - Identifies the memory feedback coefficient at the dark matter scale

Outputs:
  Outputs/sim84_coupling_survey.png   -- Lambda(Psi) for all models
  Outputs/sim84_Geff_survey.png       -- G_eff(Psi) for all models
  Outputs/sim84_derivatives.png       -- Lambda'(Psi) for all models
  Outputs/sim84_diagnostics.json
"""

import os, json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE    = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INPUTS  = os.path.join(BASE, "Inputs")
OUTPUTS = os.path.join(BASE, "Outputs")
PARAMS  = os.path.join(INPUTS, "sim84_params.json")
os.makedirs(OUTPUTS, exist_ok=True)

with open(PARAMS) as f:
    P = json.load(f)

G    = float(P["G_Newton"])
pmin = float(P["Psi_range"]["min"])
pmax = float(P["Psi_range"]["max"])
N    = int(P["Psi_range"]["N"])
tol_gr     = float(P["acceptance"]["GR_recovery_Lambda_tol"])
psi_gr_thr = float(P["acceptance"]["GR_recovery_Psi_threshold"])

Psi = np.linspace(pmin, pmax, N)

def make_Lambda(model):
    name = model["name"]
    L0   = float(model.get("Lambda0", 0.01))
    beta = float(model.get("beta", 1.0))
    if name == "linear":
        return L0 * Psi
    elif name == "quadratic":
        return L0 * Psi**2
    elif name == "gaussian":
        return L0 * (1.0 - np.exp(-beta * Psi**2))
    elif name == "tanh":
        return L0 * np.tanh(Psi**2)
    elif name == "saturating":
        return L0 * Psi**2 / (1.0 + Psi**2)
    else:
        raise ValueError(f"Unknown model: {name}")

def G_eff(Lambda_arr):
    return G / (1.0 + 16.0 * math.pi * G * Lambda_arr)

models     = P["coupling_models"]
results    = []

fig1, ax1 = plt.subplots(figsize=(7, 4.5), dpi=140)
fig2, ax2 = plt.subplots(figsize=(7, 4.5), dpi=140)
fig3, ax3 = plt.subplots(figsize=(7, 4.5), dpi=140)

for model in models:
    name     = model["name"]
    Lam      = make_Lambda(model)
    Ge       = G_eff(Lam)
    dPsi     = Psi[1] - Psi[0]
    dLam_dPsi = np.gradient(Lam, dPsi)

    # Acceptance checks
    Ge_positive  = bool(np.all(Ge > 0))
    idx_gr       = Psi < psi_gr_thr
    gr_recovery  = bool(np.max(np.abs(Lam[idx_gr])) < tol_gr) if idx_gr.any() else False

    # Memory feedback scale at Psi ~ 1 (representative DM scale)
    idx_dm    = np.argmin(np.abs(Psi - 1.0))
    feedback  = float(dLam_dPsi[idx_dm])

    results.append({
        "model": name,
        "form": model["form"],
        "G_eff_positive": Ge_positive,
        "GR_recovery": gr_recovery,
        "Lambda_at_Psi1": float(Lam[idx_dm]),
        "dLambda_dPsi_at_Psi1": feedback,
        "G_eff_at_Psi1": float(Ge[idx_dm]),
        "passed": Ge_positive and gr_recovery
    })

    ax1.plot(Psi, Lam,       label=name)
    ax2.plot(Psi, Ge / G,    label=name)
    ax3.plot(Psi, dLam_dPsi, label=name)

# ── Finalize plots ─────────────────────────────────────────────────────────
ax1.axhline(0, color="k", lw=0.5)
ax1.set_xlabel(r"$\Psi$"); ax1.set_ylabel(r"$\Lambda(\Psi)$")
ax1.set_title("Sim 84 — Coupling function survey"); ax1.legend(fontsize=8)
plt.figure(fig1.number); plt.tight_layout()
p1 = os.path.join(OUTPUTS, "sim84_coupling_survey.png")
fig1.savefig(p1); plt.close(fig1)

ax2.axhline(1.0, color="k", ls="--", lw=0.7, label="GR limit $G_{eff}=G$")
ax2.set_xlabel(r"$\Psi$"); ax2.set_ylabel(r"$G_{\rm eff}/G$")
ax2.set_title("Sim 84 — Effective Newton constant vs field"); ax2.legend(fontsize=8)
plt.figure(fig2.number); plt.tight_layout()
p2 = os.path.join(OUTPUTS, "sim84_Geff_survey.png")
fig2.savefig(p2); plt.close(fig2)

ax3.axhline(0, color="k", lw=0.5)
ax3.set_xlabel(r"$\Psi$"); ax3.set_ylabel(r"$\Lambda'(\Psi) = d\Lambda/d\Psi$")
ax3.set_title("Sim 84 — Coupling derivative (memory source amplitude)"); ax3.legend(fontsize=8)
plt.figure(fig3.number); plt.tight_layout()
p3 = os.path.join(OUTPUTS, "sim84_derivatives.png")
fig3.savefig(p3); plt.close(fig3)

# ── Acceptance summary ─────────────────────────────────────────────────────
n_failed = sum(1 for r in results if not r["passed"])
if n_failed > 0:
    failed_names = [r["model"] for r in results if not r["passed"]]
    raise RuntimeError(f"Coupling models failed constraints: {failed_names}")

diag = {
    "params": P,
    "results": results,
    "acceptance_summary": {"n_models": len(results), "n_failed": n_failed},
    "artifacts": {
        "coupling_survey": os.path.relpath(p1, BASE),
        "Geff_survey":     os.path.relpath(p2, BASE),
        "derivatives":     os.path.relpath(p3, BASE),
    }
}
diag_path = os.path.join(OUTPUTS, "sim84_diagnostics.json")
with open(diag_path, "w") as f:
    json.dump(diag, f, indent=2)

for r in results:
    status = "PASS" if r["passed"] else "FAIL"
    print(f"  {r['model']:12s}  G_eff>0={r['G_eff_positive']}  GR_recovery={r['GR_recovery']}  "
          f"Λ'(1)={r['dLambda_dPsi_at_Psi1']:.4f}  {status}")
print(f"Wrote diagnostics to {diag_path}")
