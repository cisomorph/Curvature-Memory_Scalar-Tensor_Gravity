# CMSTG: Curvature-Memory Scalar-Tensor Gravity

**A Lagrangian-based framework for curvature memory and modified gravity**

CMSTG is a covariant scalar-tensor extension of GR in which spacetime curvature couples to a non-local scalar field Ψ through a retarded Green's function, producing *curvature memory*. Dark energy emerges from the Ψ–geometry coupling without a bare cosmological constant. The theory is UV-finite through two loops in the resummed propagator, with a negative beta function driving Λ₀ to its GR-limit fixed point (Λ₀ = 0) and no infrared Landau pole. The non-local memory kernel is defined in the cosmological rest frame as a feature of the effective theory; the frame dependence and EFT scope are addressed in Paper III §2.2(d). The dark matter sector (Phase 2) introduces an ultra-light χ field seeded by Ψ.

**Status:** Phase 1 canonical action is locked and declared final. Two structural no-go theorems (Paper I) prove that no covariant scalar-tensor extension within the CMSTG class can reduce the 2.77σ DESI Y1 residual without modifying pre-recombination physics. A completeness probe (SIM144) exhausts the mixed-source scalar sector, confirming the no-go covers all 14 Tier 2 mechanism tests (SIM131–SIM144). Next decisive test: DESI Y3 (~2027).

---

## The Action

**Phase 1 canonical (locked):**

$$S = \int d^4x\sqrt{-g}\left[\frac{1+2\Lambda_0\Psi^2}{2}R - \frac{1}{2}(\nabla\Psi)^2\right] + S_{\rm SM}$$

| Parameter | Value | Status |
|-----------|-------|--------|
| Λ₀ | 0.003 M_Pl⁻² | locked |
| Ψ̄ | 2.62 M_Pl | attractor |
| F_eff = (1+2Λ₀Ψ̄²)/2 | 0.521 | locked |
| H₀ | 67.59 km/s/Mpc | best-fit |
| DESI tension | 2.77σ | structural floor |
| ln B vs ΛCDM | −0.71 | inconclusive |

---

## Papers

Four papers with LaTeX source document the full research programme, plus a master compilation. All source lives in `Papers/`.

```
Papers/
├── paper1_nogo/             Paper I   — Two Structural No-Go Theorems
├── paper2_framework/        Paper II  — Framework, Field Equations, Observational Tests
├── paper3_uv/               Paper III — UV Finiteness and Λ₀ Fixed Point
├── paper4_galactic/         Paper IV  — Galactic-Scale Constraints
└── master/                  Master compilation (Papers I–IV unified)
```

### Paper summary

| Paper | Title | Sims | Verdict |
|-------|-------|------|---------|
| I | Two structural no-go theorems for late-time modifications | SIM131–144 | FAIL/NO-GO |
| II | Framework, field equations, and 22 observational tests | SIM80–111 | PASS (2.77σ floor) |
| III | UV finiteness, two-loop; Λ₀ fixed point (negative β function); §2.2(d) kernel frame scope | SIM102–106 | PASS |
| IV | Galactic constraints; Ψ cannot replace dark matter | SIM99–100, 103 | FAIL |

---

## No-Go Theorems

Two structural no-go theorems together prove that the 2.77σ DESI residual is irreducible within CMSTG:

**Theorem 1 (Paper I — Phase 3):**
Any scalar field φ sourced by curvature R > 0 and initialized at zero grows monotonically, increasing F_eff and thereby *suppressing* H(z) at DESI redshifts. This is the wrong direction.

**Theorem 2 (Paper I — Phase 4):**
Any mechanism that raises H(z) in the interval z ∈ [0, z_rec] with fixed sound horizon r_s necessarily increases θ* = r_s/D_C* above the Planck bound 100θ* = 1.04101 ± 0.00029. The DESI fix and CMB preservation are mutually exclusive for any pure late-time modification.

**Loophole (Phase 5, if warranted):** Simultaneously increase r_s via pre-recombination physics (early dark energy, extra relativistic species). This is separately constrained by Planck CMB power spectra and constitutes a distinct theoretical programme.

---

## Simulations

All simulations live in `Ordered_Simulations/`. Every folder contains `results.md` summarising the outcome.

```
Ordered_Simulations/
  SIM###/
    SPEC.md              pre-run: action, success criteria
    run.py               executable simulation
    output.json          numerical results
    RESULT.md            post-run: verdict and structural diagnosis
    results.md           standardised summary
    figures/             plots (.pdf/.png)
```

### Phase 1 — Locked Action (SIM80–SIM111)

| Sim | Topic | Verdict |
|-----|-------|---------|
| SIM80 | Resummed memory propagator / UV falloff | PASS |
| SIM82 | Retarded kernel causality and Yukawa profile | PARTIAL |
| SIM83–86 | Memory cutoff, coupling survey, Friedmann stability | PASS |
| SIM87 | BAO full-covariance refit (BOSS+eBOSS+Lyα) | PASS |
| SIM88 | GW speed c_T = c; Solar System G_eff | PASS |
| SIM89 | Planck acoustic-peak shift diagnostic | PASS |
| SIM90 | Joint CMB+BAO parameter fit (Δχ² ≈ 0) | PASS |
| SIM91–96 | RSD, σ₈, lensing, BBN, neutrinos, Ly-α | PASS |
| SIM97 | Full Planck plikHM TTTEEE (BAO-only params) | FAIL (+7.0) |
| SIM98 | Joint plikHM+BAO refit (Δ(−2lnL) = −0.013) | PASS |
| SIM99–100 | Ψ-sector rotation curves (two approaches) | FAIL |
| SIM101–106 | BAO r_d, two-loop UV, Ward identity, RG flow | PASS |
| SIM107–108 | Solar System, GW speed | PASS |
| SIM109–110 | DESI w(z), inflation | FAIL |
| SIM111 | m₀ scan — DESI tension floor 3.44σ | STRUCTURAL |

### Phase 2 — Matter Sector (SIM112–SIM130)

| Sim | Topic | Verdict |
|-----|-------|---------|
| SIM112–113 | λΨ⁴ and hilltop quintessence | PARTIAL |
| SIM114–115 | χ-DM condensate and gradient soliton | FAIL |
| SIM116–119 | χ-DM oscillating field, halos, SPARC 161-galaxy fit | PASS |
| SIM120–121C | Joint DE+DM, Ly-α, H₀, DESI+Planck MCMC | PARTIAL (2.77σ floor) |
| SIM122–130 | Phase 3 probes (matter sector) | FAIL |

### Phase 3 — Curvature-Memory Scalars (SIM131–SIM136)

All six probes fail. Structural no-go established (Paper VII).

| Sim | Mechanism | Tension | Verdict |
|-----|-----------|---------|---------|
| SIM131 | Additive ξΨR (conformal) | 5.97σ | FAIL |
| SIM132 | Deser–Woodard analog | 3.59–3.83σ | FAIL |
| SIM133 | Gauss–Bonnet sourcing | 9.68σ | FAIL |
| SIM134 | Dilaton/exponential coupling | 3.59–3.84σ | FAIL |
| SIM135 | Bi-scalar (Ψ frozen, φ sourced by R) | 3.73–5.60σ | FAIL |
| SIM136 | Horndeski G^{μν}∂Ψ² kinetic | 3.53–3.75σ | FAIL |

### Phase 4 — Late-Time Extensions (SIM137–SIM144)

Second structural no-go established. Tier 2 mechanism space exhaustively covered by SIM144.

| Sim | Tier | Mechanism | Verdict |
|-----|------|-----------|---------|
| SIM137 | T1 | SPARC failure-mode analysis | STRUCTURAL\_PATTERN |
| SIM138 | T1 | DESI per-bin sensitivity | DISTRIBUTED |
| SIM139 | T1 | RSD shape diagnostic | SHAPE\_FIXABLE |
| SIM140 | T2 | Step potential (P4-A) | PREDICTED FAIL (not run) |
| SIM141 | T2 | Running Λ₀(a) / BD analog (P4-B) | FAIL/STRUCTURAL — θ* 63σ |
| SIM142 | T2 | Galileon G₃(Ψ)□Ψ (P4-C) | FAIL/STRUCTURAL |
| SIM143 | T2 | Bi-scalar Ψ+φ quintessence (P4-D) | FAIL/STRUCTURAL |
| SIM144 | T2 | Mixed-source φ (ξ_R + β_m) — completeness probe (P4-E) | FAIL/STRUCTURAL |
| SIM145–146 | T3 | UV recheck / predictions | DEFERRED (no Tier 2 PASS) |

---

## Observational Predictions

Phase 1 canonical CMSTG makes the following predictions, distinguishable from ΛCDM by DESI Y3, Euclid, or CMB-S4:

1. **DESI Y3 BAO tension:** If the 2.77σ floor persists at ≥3σ with more data, CMSTG requires pre-recombination new physics (Phase 5). If tension falls below 2σ, Phase 1 is fully validated.
2. **χ-DM mass:** m₂₂ ≈ 0.28–0.34 for low-mass gas-dominated dwarfs and LSBs (v_flat ≲ 120 km/s). Testable by next-generation stellar kinematics.
3. **G_eff/G = 1 − 15 ppm** at z = 0 — negligible at current precision but in principle testable by future GW standard-siren surveys.

---

## Repository Structure

```
Curvature-Memory_Scalar-Tensor_Gravity/
├── README.md                         ← this file
├── RESEARCH_RULES.md                 standing rules for all simulation work
├── CITATION.cff                      citation metadata
├── Papers/                           all papers (LaTeX source)
└── Ordered_Simulations/              all simulations (SIM80–SIM146)
```

---

## Running a Simulation

Requirements: Python 3.10+, NumPy, SciPy, Matplotlib.

```bash
cd Ordered_Simulations/SIM98
python run.py
# output.json and figures/ are written in-place
```

SIM97/98 additionally require [CLASS](https://github.com/lesgourg/class_public) and the [Planck 2018 plikHM likelihood](https://pla.esac.esa.int/pla/#cosmology) with [clipy](https://github.com/benabed/clipy).

---

## License

CC BY-NC 4.0 — Attribution-NonCommercial. Free to share and adapt with credit; commercial use prohibited.

**Author:** Christopher Robert Barrick Wilson  
**Contact:** C.Isomorph@gmail.com  
**Last updated:** 2026-05-14
