# RIFT: Recursive Intelligence-Field Theory

**A Lagrangian-based framework for curvature memory and emergent gravitation**

This repository contains the simulation code and paper draft for RIFT — a
modified gravity theory in which spacetime curvature couples to a non-local
scalar field Ψ through a retarded Green's function, producing *curvature memory*.
The theory replaces dark matter and dark energy with a single recursive field
derived from a modified Einstein-Hilbert action.

---

## The Action

$$S = \int d^4x\sqrt{-g}\left[\left(\frac{1}{16\pi G} + \Lambda(\Psi)\right)R + \frac{1}{2}(\nabla\Psi)^2 - \frac{1}{2}m(\Psi)^2\Psi^2 - V(\Psi)\right] + S_{\rm matter}$$

Λ(Ψ) = Λ₀Ψ² couples the recursive field to the Ricci scalar, modifying the
effective gravitational constant: G_eff = G / (1 + 16πΛ(Ψ)).

---

## Paper

`paper/rift_paper_v2.tex` / `rift_paper_v2.pdf` — 22 pages, compiles with pdflatex.

**Key results (Sims 87–96):**

| Test | Result |
|------|--------|
| BAO full-covariance refit (SIM87) | χ²/dof = 1.22, Δχ²(RIFT−ΛCDM) = 0 |
| CMB TT via CLASS (SIM88) | RMS(ΔC_ℓ/C_ℓ) = 2.93% |
| Planck TT likelihood (SIM89) | ΛCDM baseline χ² = 6.20 ✓ |
| Joint CMB+BAO fit (SIM90) | Δχ² = 0; H₀=67.59, Ωm=0.312, Λ₀=0.008 |
| Λ₀ sensitivity scan (SIM91) | BAO flat; growth dev >0.1% at Λ₀≥0.05 |
| Non-linear σ₈ / HALOFIT (SIM92) | Δσ₈ < 0.007% at Λ₀=0.003 |
| Bayesian model comparison (SIM93) | ln B = −0.71 (Inconclusive — not excluded) |
| DESI Y1 BAO refit (SIM94) | χ²/dof = 1.36, Δχ² = 0 |
| CMB polarization EE/TE (SIM95) | RMS(ΔEE/EE) = 0.19%, RMS(ΔTE/TE) = 0.43% |
| RSD / f·σ₈ growth rate (SIM96) | χ²/dof = 0.86, Δχ²(RIFT−ΛCDM) = −0.14 |

---

## Simulations

Each simulation lives in `simulations/SIMxx/` with three subdirectories:

```
simulations/
  SIMxx/
    Inputs/         # sim_params.json — all parameters
    Published_Scripts/  # simulation_runner_rift.py — self-contained runner
    Outputs/        # diagnostics.json, figures (.pdf/.png), data
```

| Sim | Topic | Status |
|-----|-------|--------|
| SIM82 | Retarded kernel / causality | PASS |
| SIM83 | Memory cutoff oscillation period | PASS |
| SIM84 | Coupling form survey (G_eff > 0) | PASS |
| SIM85 | Friedmann memory term validation | PASS |
| SIM86 | Observer field / back-reaction stability | PASS |
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

---

## Running a simulation

Requirements: Python 3.10+, NumPy, SciPy, Matplotlib. SIM88/SIM95 also require
[CLASS](https://github.com/lesgourg/class_public) compiled at the path set in
the `class_executable` field of the params JSON.

```bash
cd simulations/SIM96
python Published_Scripts/simulation_runner_rift.py
# outputs go to Outputs/
```

---

## Best-fit parameters (SIM90 joint CMB+BAO)

| Parameter | RIFT | ΛCDM (Planck 2018) |
|-----------|------|---------------------|
| H₀ [km/s/Mpc] | 67.59 | 67.36 |
| Ωm | 0.312 | 0.3153 |
| rd [Mpc] | 147.56 | — |
| Λ₀ | 0.008 | 0 (by definition) |
| σ₈ | 0.824 | 0.811 |

---

## License

CC BY-NC 4.0 — Attribution-NonCommercial. Free to share and adapt with credit; commercial use prohibited.

Author: cisomorph  
Date: 2026-04-09
