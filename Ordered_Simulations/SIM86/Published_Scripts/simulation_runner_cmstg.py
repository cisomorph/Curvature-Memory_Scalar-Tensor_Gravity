#!/usr/bin/env python3
# Sim 86 — CMSTG Observer Field (Self-Consistent Back-Reaction)
"""
Tests the closed CMSTG loop: Psi -> R -> Psi.

In the full CMSTG system, Psi is sourced by R through the retarded kernel
(eq:memory), while R itself depends on Psi through the gravitational
equation (eq:gravity) via G_eff(Psi) and T^Psi_{mu nu}.

This creates a self-referential (observer-field) structure: the field's
own curvature contribution feeds back into its source. The question is
whether this loop is:
  - Stable: Psi relaxes to a fixed point (attractor)
  - Oscillatory: Psi undergoes bounded oscillations
  - Runaway: Psi diverges (theory is ill-posed for these parameters)

We test this in the FLRW background using the full coupled system:

  FLRW scalar equation:
    Psi'' + 3H Psi' + m_eff^2 Psi = Lambda'(Psi) * R + V'(Psi)
    R = 6*(dH/dt + 2H^2)  [Ricci scalar, computed self-consistently]

  Modified Friedmann:
    3H^2 = 8pi G_eff [rho_m + rho_Psi - 6H Lambda'(Psi) Psi']

The "observer field" interpretation: in CMSTG, an observer is a region
of spacetime that has accumulated curvature memory (high Psi). Back-reaction
tests whether observers are self-sustaining (attractor fixed point) or
whether memory self-amplifies without bound (runaway).

Acceptance checks:
  1. Zero IC -> Psi stays zero (vacuum stability, GR recovery)
  2. Small IC -> Psi remains bounded and relaxes (no runaway)
  3. Large IC -> field either converges or explicitly flags divergence
  4. Friedmann constraint satisfied throughout

Outputs:
  Outputs/sim86_Psi_vs_a.png         -- Psi(a) for all ICs
  Outputs/sim86_Geff_vs_a.png        -- G_eff(a) for all ICs
  Outputs/sim86_phase_portrait.png   -- Psi vs dPsi/d(ln a)
  Outputs/sim86_diagnostics.json
"""

import os, json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE    = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INPUTS  = os.path.join(BASE, "Inputs")
OUTPUTS = os.path.join(BASE, "Outputs")
PARAMS  = os.path.join(INPUTS, "sim86_params.json")
os.makedirs(OUTPUTS, exist_ok=True)

with open(PARAMS) as f:
    P = json.load(f)

G        = float(P["G_Newton"])
m0       = float(P["field"]["m0"])
alpha    = float(P["field"]["alpha"])
beta_s   = float(P["field"]["beta"])
Lambda0  = float(P["field"]["Lambda0"])
a_ini    = float(P["grid"]["a_ini"])
a_fin    = float(P["grid"]["a_fin"])
N_steps  = int(P["grid"]["N_steps"])
gr_lim   = float(P["acceptance"]["GR_limit_max_Psi"])
stab_max = float(P["acceptance"]["stability_max_Psi"])
cstr_tol = float(P["acceptance"]["constraint_tol"])

rho_crit0 = 3.0 / (8.0 * math.pi * G)
omega_m   = 0.3
omega_r   = 9e-5
rho_m0    = omega_m * rho_crit0
rho_r0    = omega_r * rho_crit0

d_ln_a = (math.log(a_fin) - math.log(a_ini)) / N_steps

def Lambda_fn(Psi):
    return Lambda0 * Psi**2

def dLambda_fn(Psi):
    return 2.0 * Lambda0 * Psi

def m_eff_sq(Psi):
    return m0**2 * (1.0 + alpha * Psi**2 * math.exp(-beta_s * Psi**2))

def V(Psi):       return 0.0
def dV_dPsi(Psi): return 0.0

def H_from_constraint(a, Psi, Pi, H_prev):
    """Solve 3H^2 = 8pi G_eff [rho_total - 6H Lambda' Pi H] for H."""
    Lam  = Lambda_fn(Psi)
    dLam = dLambda_fn(Psi)
    rho_m_a = rho_m0 / a**3
    rho_r_a = rho_r0 / a**4
    H = H_prev
    for _ in range(8):
        dPsi_dt = Pi * H
        rho_Psi = 0.5 * dPsi_dt**2 + 0.5 * m_eff_sq(Psi) * Psi**2 + V(Psi)
        rho_tot = rho_m_a + rho_r_a + rho_Psi
        coeff   = 3.0 * (1.0 + 16.0*math.pi*G*Lam) + 48.0*math.pi*G*dLam*Pi
        if coeff <= 0:
            H = 1e-30; break
        H2 = 8.0 * math.pi * G * rho_tot / coeff
        H = math.sqrt(max(H2, 0.0))
    return H

def R_scalar(H, dH_dlna):
    """R = 6*(dH/dt + 2H^2) = 6*(H * dH/d(lna) + 2H^2) = 6H^2*(dH/Hd(lna) + 2)."""
    return 6.0 * (H * dH_dlna + 2.0 * H**2)

def derivs(ln_a, Psi, Pi, H, H_prev, d_ln_a_step):
    """RHS in terms of ln(a): returns (dPsi/d_lna, dPi/d_lna)."""
    a    = math.exp(ln_a)
    Lam  = Lambda_fn(Psi)
    dLam = dLambda_fn(Psi)
    m2   = m_eff_sq(Psi)
    # Estimate dH/d(ln a) from finite difference
    dH_dlna = (H - H_prev) / d_ln_a_step if d_ln_a_step > 0 else 0.0
    R = R_scalar(H, dH_dlna)
    if H**2 < 1e-60:
        return Pi, 0.0
    dPsi_dlna  = Pi
    d2Psi_dlna = (
        -3.0 * Pi
        - m2 * Psi / H**2
        + dLam * R / H**2
        + dV_dPsi(Psi) / H**2
    )
    return dPsi_dlna, d2Psi_dlna

def integrate_IC(Psi0, dPsi0, label):
    ln_a = math.log(a_ini)
    a    = a_ini
    H    = math.sqrt(8.0*math.pi*G*rho_r0/3.0) / a**2
    Psi  = Psi0
    Pi   = dPsi0 / H if H > 0 else 0.0
    H_prev = H

    store_every = max(1, N_steps // 3000)
    a_arr, Psi_arr, Pi_arr, H_arr, Geff_arr, cstr_arr = [], [], [], [], [], []

    runaway = False

    for step in range(N_steps):
        a = math.exp(ln_a)
        H = H_from_constraint(a, Psi, Pi, H)

        if step % store_every == 0:
            Lam = Lambda_fn(Psi)
            Ge  = G / (1.0 + 16.0*math.pi*G*Lam)
            dPsi_dt = Pi * H
            rho_Psi = 0.5*dPsi_dt**2 + 0.5*m_eff_sq(Psi)*Psi**2
            dLam = dLambda_fn(Psi)
            rho_m_a = rho_m0/a**3 + rho_r0/a**4
            lhs = 3.0*H**2
            rhs = 8.0*math.pi*Ge*(rho_m_a + rho_Psi - 6.0*H*dLam*dPsi_dt)
            cstr = abs(lhs - rhs)/(lhs + 1e-30)
            a_arr.append(a); Psi_arr.append(Psi)
            Pi_arr.append(Pi); H_arr.append(H)
            Geff_arr.append(Ge); cstr_arr.append(cstr)

        if abs(Psi) > stab_max:
            runaway = True
            print(f"  [{label}] RUNAWAY at a={a:.4e}, Psi={Psi:.3e}")
            break

        # RK4
        H2 = H_from_constraint(math.exp(ln_a + 0.5*d_ln_a), Psi, Pi, H)
        k1p, k1Pi = derivs(ln_a, Psi, Pi, H, H_prev, d_ln_a)
        k2p, k2Pi = derivs(ln_a + 0.5*d_ln_a,
                            Psi + 0.5*d_ln_a*k1p, Pi + 0.5*d_ln_a*k1Pi, H2, H, 0.5*d_ln_a)
        k3p, k3Pi = derivs(ln_a + 0.5*d_ln_a,
                            Psi + 0.5*d_ln_a*k2p, Pi + 0.5*d_ln_a*k2Pi, H2, H, 0.5*d_ln_a)
        H3 = H_from_constraint(math.exp(ln_a + d_ln_a), Psi + d_ln_a*k3p, Pi + d_ln_a*k3Pi, H)
        k4p, k4Pi = derivs(ln_a + d_ln_a,
                            Psi + d_ln_a*k3p, Pi + d_ln_a*k3Pi, H3, H, d_ln_a)

        H_prev = H
        Psi   += d_ln_a * (k1p  + 2*k2p  + 2*k3p  + k4p)  / 6.0
        Pi    += d_ln_a * (k1Pi + 2*k2Pi + 2*k3Pi + k4Pi) / 6.0
        ln_a  += d_ln_a

    return (np.array(a_arr), np.array(Psi_arr), np.array(Pi_arr),
            np.array(H_arr), np.array(Geff_arr), np.array(cstr_arr), runaway)

# ── Run all ICs ─────────────────────────────────────────────────────────────
fig1, ax1 = plt.subplots(figsize=(7, 5), dpi=140)
fig2, ax2 = plt.subplots(figsize=(7, 5), dpi=140)
fig3, ax3 = plt.subplots(figsize=(7, 5), dpi=140)

all_results = []

for ic in P["initial_conditions"]:
    name     = ic["name"]
    Psi0     = float(ic["Psi_ini"])
    dPsi0    = float(ic["dPsi_ini"])
    print(f"  Running IC: {name} (Psi0={Psi0})", flush=True)
    a_a, Psi_a, Pi_a, H_a, Ge_a, cstr_a, runaway = integrate_IC(Psi0, dPsi0, name)

    ax1.semilogx(a_a, Psi_a, label=f"{name}: $\\Psi_0={Psi0}$")
    ax2.semilogx(a_a, Ge_a / G, label=f"{name}: $\\Psi_0={Psi0}$")
    ax3.plot(Psi_a, Pi_a, label=f"{name}: $\\Psi_0={Psi0}$", alpha=0.7)

    max_cstr = float(np.max(cstr_a)) if len(cstr_a) > 0 else float("nan")
    final_Psi = float(Psi_a[-1]) if len(Psi_a) > 0 else float("nan")

    # Acceptance
    if name == "zero":
        passed = bool(np.all(np.abs(Psi_a) < gr_lim))
        note = "GR vacuum stability"
    else:
        passed = bool(not runaway and max_cstr < cstr_tol)
        note = "bounded evolution" if not runaway else "RUNAWAY"

    all_results.append({
        "IC": name, "Psi0": Psi0,
        "final_Psi": final_Psi,
        "max_constraint_violation": max_cstr,
        "runaway": bool(runaway),
        "passed": passed,
        "note": note
    })

# ── Plots ───────────────────────────────────────────────────────────────────
ax1.axhline(0, color="k", lw=0.4)
ax1.set_xlabel("$a$"); ax1.set_ylabel(r"$\Psi(a)$")
ax1.set_title("Sim 86 — Observer field: $\\Psi(a)$ for varying IC")
ax1.legend(fontsize=8)
plt.figure(fig1.number); plt.tight_layout()
p1 = os.path.join(OUTPUTS, "sim86_Psi_vs_a.png"); fig1.savefig(p1); plt.close(fig1)

ax2.axhline(1.0, color="k", ls="--", lw=0.7, label="GR: $G_{eff}=G$")
ax2.set_xlabel("$a$"); ax2.set_ylabel(r"$G_{\rm eff}(a)/G$")
ax2.set_title("Sim 86 — Effective Newton constant evolution")
ax2.legend(fontsize=8)
plt.figure(fig2.number); plt.tight_layout()
p2 = os.path.join(OUTPUTS, "sim86_Geff_vs_a.png"); fig2.savefig(p2); plt.close(fig2)

ax3.axhline(0, color="k", lw=0.4); ax3.axvline(0, color="k", lw=0.4)
ax3.set_xlabel(r"$\Psi$"); ax3.set_ylabel(r"$d\Psi/d(\ln a)$")
ax3.set_title("Sim 86 — Phase portrait: attractor structure")
ax3.legend(fontsize=8)
plt.figure(fig3.number); plt.tight_layout()
p3 = os.path.join(OUTPUTS, "sim86_phase_portrait.png"); fig3.savefig(p3); plt.close(fig3)

# ── Acceptance ──────────────────────────────────────────────────────────────
failed = [r for r in all_results if not r["passed"]]
if failed:
    names = [r["IC"] for r in failed]
    raise RuntimeError(f"Observer field test failed for ICs: {names}")

# ── Diagnostics ─────────────────────────────────────────────────────────────
diag = {
    "description": (
        "Self-consistent back-reaction test. Psi sources R through G_eff and T^Psi, "
        "which feeds back into the scalar equation via Lambda'(Psi)*R. "
        "Tests attractor stability of the closed CMSTG loop."
    ),
    "params": P,
    "results": all_results,
    "interpretation": {
        "zero_IC": "GR recovery: Psi=0 is a fixed point",
        "nonzero_IC": "Non-zero Psi evolves self-consistently; attractor behaviour indicates memory stability",
        "phase_portrait": "Trajectories converging to Psi=0 confirm GR attractor; nonzero attractors would indicate persistent memory",
    },
    "artifacts": {
        "Psi_vs_a":       os.path.relpath(p1, BASE),
        "Geff_vs_a":      os.path.relpath(p2, BASE),
        "phase_portrait": os.path.relpath(p3, BASE),
    }
}
diag_path = os.path.join(OUTPUTS, "sim86_diagnostics.json")
with open(diag_path, "w") as f:
    json.dump(diag, f, indent=2)

for r in all_results:
    status = "PASS" if r["passed"] else "FAIL"
    print(f"  {r['IC']:10s}  Psi_final={r['final_Psi']:.3e}  "
          f"cstr={r['max_constraint_violation']:.2e}  "
          f"runaway={r['runaway']}  {status}  [{r['note']}]")
print(f"\nWrote diagnostics to {diag_path}")
