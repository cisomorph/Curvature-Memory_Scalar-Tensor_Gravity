# Simulation Programme

70+ simulations across four phases (SIM80–SIM153). Each sim folder lives at `Ordered_Simulations/SIM###/` in the repository and contains `Inputs/`, `Outputs/`, `Published_Scripts/`, and `results.md`.

**Verdict legend:** PASS · FAIL · PARTIAL · STRUCTURAL · DEFERRED

---

## Phase 1 — Locked Action (SIM80–SIM111)

Phase 1 establishes the canonical action and runs 22 observational tests. The 2.77σ DESI Y1 BAO residual is first identified here as a structural floor.

| SIM | Topic | Verdict |
|-----|-------|---------|
| SIM80 | Resummed memory propagator / UV falloff | PASS |
| SIM82 | Retarded kernel causality and Yukawa profile | PARTIAL |
| SIM83–86 | Memory cutoff, coupling survey, Friedmann stability | PASS |
| SIM87 | BAO full-covariance refit (BOSS + eBOSS + Lyα) | PASS |
| SIM88 | GW speed c_T = c; Solar System G_eff | PASS |
| SIM89 | Planck acoustic-peak shift diagnostic | PASS |
| SIM90 | Joint CMB+BAO parameter fit (Δχ² ≈ 0) | PASS |
| SIM91–96 | RSD, σ₈, lensing, BBN, neutrinos, Ly-α | PASS |
| SIM97 | Full Planck plikHM TTTEEE (BAO-only params) | FAIL — Δχ² = +7.0 |
| SIM98 | Joint plikHM + BAO refit | PASS — Δ(−2lnL) = −0.013 |
| SIM99–100 | Ψ-sector rotation curves (two approaches) | FAIL → Paper IV |
| SIM101–106 | BAO rₛ, two-loop UV, Ward identity, RG flow | PASS |
| SIM107–108 | Solar System, GW speed | PASS |
| SIM109–110 | DESI w(z), inflation | FAIL |
| SIM111 | m₀ scan — DESI tension floor 3.44σ | STRUCTURAL |

**Phase 1 outcome:** Canonical action declared final. The 2.77σ DESI Y1 BAO residual is identified as a structural floor, not a tuning failure.

---

## Phase 2 — Matter Sector (SIM112–SIM130)

Extends the framework with a matter-sector χ-DM field. The χ-DM oscillating field passes the SPARC 161-galaxy rotation curve fit. The joint DE+DM fit converges on the same 2.77σ floor.

| SIM | Topic | Verdict |
|-----|-------|---------|
| SIM112–113 | λΨ⁴ and hilltop quintessence | PARTIAL |
| SIM114–115 | χ-DM condensate and gradient soliton | FAIL |
| SIM116–119 | χ-DM oscillating field, halos, SPARC 161-galaxy fit | PASS |
| SIM120–121C | Joint DE+DM, Ly-α, H₀, DESI+Planck MCMC | PARTIAL — 2.77σ floor |
| SIM122–130 | Phase 3 probes (matter sector) | FAIL |

**Phase 2 outcome:** χ-DM is the surviving dark matter route. The DE tension floor is unchanged by matter-sector extensions.

---

## Phase 3 — Curvature-Memory Scalars (SIM131–SIM136)

Six curvature-sourced scalar mechanisms tested. All fail. Theorem 1 (Paper I) is established from this phase.

| SIM | Mechanism | Tension | Verdict |
|-----|-----------|---------|---------|
| SIM131 | Additive ξΨR (conformal coupling) | 5.97σ | FAIL |
| SIM132 | Deser–Woodard analog | 3.59–3.83σ | FAIL |
| SIM133 | Gauss–Bonnet sourcing | 9.68σ | FAIL |
| SIM134 | Dilaton / exponential coupling | 3.59–3.84σ | FAIL |
| SIM135 | Bi-scalar (Ψ frozen, φ sourced by R) | 3.73–5.60σ | FAIL |
| SIM136 | Horndeski G_μν∂Ψ² kinetic | 3.53–3.75σ | FAIL |

**Phase 3 outcome → Theorem 1:** Any scalar sourced by R > 0 grows monotonically, increasing F_eff and suppressing H(z) at DESI redshifts. All such mechanisms worsen the tension.

---

## Phase 4 — Late-Time Extensions (SIM137–SIM144)

Tier 1 diagnostics establish the failure pattern; Tier 2 exhausts the remaining mechanism space. Theorem 2 (Paper I) is established from this phase.

| SIM | Tier | Mechanism | Verdict |
|-----|------|-----------|---------|
| SIM137 | T1 | SPARC failure-mode analysis | STRUCTURAL PATTERN |
| SIM138 | T1 | DESI per-bin sensitivity | DISTRIBUTED |
| SIM139 | T1 | RSD shape diagnostic | SHAPE FIXABLE |
| SIM140 | T2 | Step potential (P4-A) | PREDICTED FAIL (not run) |
| SIM141 | T2 | Running Λ₀(a) / BD analog (P4-B) | FAIL — θ* 63σ |
| SIM142 | T2 | Galileon G₃(Ψ)□Ψ (P4-C) | FAIL |
| SIM143 | T2 | Bi-scalar Ψ+φ quintessence (P4-D) | FAIL |
| SIM144 | T2 | Mixed-source φ (ξ_R + β_m) — completeness probe | FAIL |
| SIM145–146 | T3 | UV recheck / predictions | DEFERRED |

**Phase 4 outcome → Theorem 2:** Any mechanism raising H(z) at late times violates the Planck CMB acoustic-peak bound θ*. Tier 2 mechanism space exhaustively closed by SIM144.

---

## Phase 5 — Loophole Probes (SIM147–SIM150)

| SIM | Topic | Verdict |
|-----|-------|---------|
| SIM147 | Ψ_pre loophole — probe 1 | FAIL |
| SIM148 | Ψ_pre loophole — probe 2 | FAIL |
| SIM149 | Ψ_local loophole | FAIL |
| SIM150 | Λ₀ sweep: G_eff ceiling universality | STRUCTURAL |

**Phase 5 outcome:** Both Ψ_pre and Ψ_local loopholes triply closed. G_eff ceiling is universal across the Λ₀ parameter space.

---

## Phase A — Pre-Recombination Options (SIM151–SIM153)

| SIM | Topic | Verdict |
|-----|-------|---------|
| SIM151 | Option A baseline | PASS (reference) |
| SIM152 | Option B — early dark energy coupling | FAIL — χ²/dof = 91.9 |
| SIM153 | Option C — extended mechanism | HALT-PC2 |

**Phase A outcome:** Option B fails badly. Option C halted at Phase C2 checkpoint. Companion no-go note written. Paper VII contingency roadmap active.

---

## Simulation Folder Structure

Each sim at `Ordered_Simulations/SIM###/`:

```
SIM###/
├── Inputs/
│   └── sim###_params.json       # All run parameters
├── Outputs/
│   ├── sim###_diagnostics.json  # χ², σ, convergence flags
│   └── *.pdf / *.png            # Plots and figures
├── Published_Scripts/
│   └── simulation_runner_cmstg.py
└── results.md                   # Human-readable verdict and notes
```

Simulation index: [`Ordered_Simulations/CMSTG_Simulation_Folders.csv`](https://github.com/cisomorph/Curvature-Memory_Scalar-Tensor_Gravity/blob/main/Ordered_Simulations/CMSTG_Simulation_Folders.csv)
