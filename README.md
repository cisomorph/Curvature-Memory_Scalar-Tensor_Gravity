# RIFT: Recursive Intelligence-Field Theory

**A Lagrangian-based framework for curvature memory and emergent gravitation**

This repository contains the simulation code and papers for RIFT — a modified gravity theory in which spacetime curvature couples to a non-local scalar field Ψ through a retarded Green's function, producing *curvature memory*. The theory derives dark energy from the scalar field without a bare cosmological constant; dark matter is addressed in Phase 2 via a direct matter coupling.

---

## The Action

**Phase 1 (locked):**

$$S = \int d^4x\sqrt{-g}\left[\frac{M_{\rm Pl}^2 + 2\Lambda_0\Psi^2}{2}R - \frac{1}{2}(\nabla\Psi)^2 - \frac{1}{2}m_0^2\Psi^2\right] + S_{\rm SM}$$

Parameters locked at: **Λ₀ = 0.003**, **m₀ = 0.001 H₀**.

**Phase 2 (unlocked):** bare Λ removed; dark energy from scalar potential V(Ψ); dark matter from direct coupling βΨ²L_m.

---

## Papers

Research is organised by phase. Full details in [`simulations/SIM_MAP.md`](simulations/SIM_MAP.md).

```
paper/
├── phase1/                  # Phase 1 — Locked Action
│   ├── I_original/          # Paper I: Theory and foundations (SIM80–SIM86)
│   ├── II_observational/    # Paper II: Full observational programme (SIM87–SIM111)
│   ├── III_sound_horizon/   # Paper III: BAO sound horizon (SIM101)
│   ├── IV_uv_structure/     # Paper IV: UV structure and RG flow (SIM102–SIM106)
│   ├── V_dark_matter/       # Paper V: Dark matter failure analysis (SIM99–SIM103)
│   ├── master/              # Master compilation — Papers I–V unified, 46pp
│   └── v2/                  # v2 draft (precursor to five-paper split)
└── phase2/                  # Phase 2 — Unlocked Action
    └── P2_I_framework/      # P2 Paper I: Framework derivation (SIM112–SIM115)
```

---

## Simulations

Each simulation lives in `simulations/SIMxx/` with three subdirectories:

```
simulations/
  SIMxx/
    Inputs/              # sim_params.json — all parameters
    Published_Scripts/   # simulation_runner_rift.py — self-contained runner
    Outputs/             # diagnostics.json, figures (.pdf/.png), data
  SIM_MAP.md             # Full simulation index with results and notes
```

### Phase 1 — Locked Action

#### Paper I: Theory Foundations (SIM80–SIM86)

| Sim | Topic | Result |
|-----|-------|--------|
| SIM80 | Resummed memory propagator / UV falloff validation | PASS |
| SIM82 | Retarded kernel / causality | PARTIAL |
| SIM83 | Memory cutoff oscillation period | PASS |
| SIM84 | Coupling form survey (G_eff > 0 throughout) | PASS |
| SIM85 | Friedmann memory term validation | PASS |
| SIM86 | Observer field / back-reaction stability | PASS |

#### Paper II: Observational Tests (SIM87–SIM111)

| Sim | Topic | Result |
|-----|-------|--------|
| SIM87 | BAO full-covariance refit (BOSS+eBOSS+Lyα) | PASS |
| SIM88 | CMB TT via CLASS full Boltzmann | PASS |
| SIM89 | Planck TT likelihood baseline | PASS |
| SIM90 | Joint CMB+BAO parameter fit | PASS |
| SIM91 | Λ₀ sensitivity scan | PASS |
| SIM92 | Non-linear structure growth (σ₈, HALOFIT) | PASS |
| SIM93 | Bayesian model comparison (Savage-Dickey) | PASS |
| SIM94 | DESI Year 1 BAO refit | PASS |
| SIM95 | CMB polarization EE/TE via CLASS | PASS |
| SIM96 | RSD / f·σ₈ growth rate (BOSS+eBOSS) | PASS |
| SIM97 | Real Planck 2018 plikHM TTTEEE likelihood | PASS* |
| SIM98 | Joint plikHM+BAO re-optimisation | PASS |
| SIM107 | Solar System viability / fifth-force bounds | PASS |
| SIM108 | Gravitational wave speed and GW constraints | PASS |
| SIM109 | Dark energy w(z) vs Planck and DESI | PASS / FAIL† |
| SIM110 | Slow-roll inflation — locked action | FAIL |
| SIM111 | m₀ scan for DESI tension — full parameter sweep | FAIL |

#### Paper III: BAO Sound Horizon (SIM101)

| Sim | Topic | Result |
|-----|-------|--------|
| SIM101 | Sound horizon r_d from RIFT perturbation theory | PASS‡ |

#### Paper IV: UV Structure (SIM102–SIM106)

| Sim | Topic | Result |
|-----|-------|--------|
| SIM102 | One-loop Ψ self-energy / UV finiteness via memory damping | PASS |
| SIM104 | Complete one-loop UV structure (Ψ, graviton WFR, vertex) | PASS |
| SIM105 | Beta function / RG flow: Λ₀ asymptotic freedom | PASS |
| SIM106 | Two-loop graviton sector UV-finiteness | PASS |

#### Paper V: Dark Matter (SIM99–SIM103)

| Sim | Topic | Result |
|-----|-------|--------|
| SIM99 | Galaxy rotation curves — Yukawa halo scan | FAIL |
| SIM100 | 1+1D time-domain galactic field / tachyonic condensation | FAIL |
| SIM103 | Vainshtein mechanism for dark matter halos | FAIL |

---

### Phase 2 — Unlocked Action

#### P2 Paper I: Framework (SIM112–SIM115)

| Sim | Topic | Result |
|-----|-------|--------|
| SIM112 | λΨ⁴ quintessence — no bare Λ | FAIL |
| SIM113 | SSB hilltop V = λ(Ψ²−v²)² quintessence | PARTIAL |
| SIM114 | βΨ²ρ_m condensate: galactic rotation curves | FAIL |
| SIM115 | Gradient soliton dark matter | FAIL |
| SIM116 | Ψ-switch ξ condensate halo mechanism | PLANNED |
| SIM117 | Memory-of-memory ξ field (Level-2 recursion) | PLANNED |

---

### Notes

*SIM97: Δ(−2lnL) = +7.0 at SIM90 joint best-fit. Driven by parameter offset between SIM90 and Planck CMB optimum, not G_eff (42 ppm). RIFT at the Planck parameter point gives Δ ≈ 0.

†SIM109: PASS vs Planck (w₀ = −0.992, 1.3σ); FAIL vs DESI (3.6σ). SIM111 confirms the DESI tension floor is structural: minimum achievable tension across the full m₀ scan is 3.44σ. Unlocking the action (Phase 2) is required.

‡SIM101: RIFT shift in r_d is < 22 ppm — degenerate with ΛCDM within current BAO precision.

SIM81 was an internal scratch test and is not included.

---

## Best-fit Parameters (SIM90 Joint CMB+BAO)

| Parameter | RIFT | ΛCDM (Planck 2018) |
|-----------|------|---------------------|
| H₀ [km/s/Mpc] | 67.59 | 67.36 |
| Ωm | 0.312 | 0.3153 |
| rd [Mpc] | 147.56 | — |
| Λ₀ | 0.008 | 0 (by definition) |
| σ₈ | 0.824 | 0.811 |

---

## Running a Simulation

Requirements: Python 3.10+, NumPy, SciPy, Matplotlib.

SIM88, SIM95, and SIM97 additionally require [CLASS](https://github.com/lesgourg/class_public) compiled at the path set in `class_executable` in the params JSON. SIM97 also requires the [Planck 2018 likelihood package](https://pla.esac.esa.int/pla/#cosmology) and [clipy](https://github.com/benabed/clipy) (v0.15).

```bash
cd simulations/SIM96
python Published_Scripts/simulation_runner_rift.py
# outputs go to Outputs/
```

---

## License

CC BY-NC 4.0 — Attribution-NonCommercial. Free to share and adapt with credit; commercial use prohibited.

Author: Christopher Robert Barrick Wilson
Date: 2026-04-18
