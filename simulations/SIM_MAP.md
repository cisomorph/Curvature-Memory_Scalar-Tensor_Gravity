# RIFT Simulation Map

All published simulations for the Recursive Intelligence-Field Theory (RIFT) cosmology programme.
Organised by phase and paper. Results: **PASS** / **FAIL** / **PARTIAL** (see notes).

---

## Phase 1 — Locked Action

Action: `S = ∫d⁴x√(-g) [(M_Pl²+2Λ₀Ψ²)/2 R - ½(∇Ψ)² - ½m₀²Ψ² - ρ_Λ] + S_SM`  
Parameters locked at: **Λ₀ = 0.003**, **m₀ = 0.001 H₀**, bare **ρ_Λ** retained.  
Papers: `paper/phase1/`

---

### Paper I — Theoretical Foundation (`I_original/`)

| Sim | Topic | Result | Notes |
|-----|-------|--------|-------|
| SIM80 | UV Finiteness — recursive memory self-energy, bubble integral softening | PASS | Bubble slope softened vs baseline; Λ_mem fixed |
| SIM82 | Retarded Kernel — G_R(t,r) in 1+1D flat spacetime, causality check | PARTIAL | Causality & memory accumulation pass; Yukawa slope fail (0.11 vs −1.0) |
| SIM83 | Memory Cutoff — L_mem = 1/m_eff as function of field amplitude | PASS | |
| SIM84 | Coupling Survey — Λ(Ψ) functional forms; verify G_eff > 0 and GR recovery | PASS | |
| SIM85 | FLRW Friedmann Test — modified equations reduce to GR in limits | PASS | |
| SIM86 | Back-Reaction — Ψ sources R via G_eff; attractor stability of closed RIFT loop | PASS | |

---

### Paper II — Observational Tests (`II_observational/`)

| Sim | Topic | Result | Notes |
|-----|-------|--------|-------|
| SIM87 | BAO Full-Covariance Refit — off-diagonal BOSS DR12 + eBOSS DR16 correlations | PASS | Fixes diagonal-only vulnerability |
| SIM88 | CMB Full CLASS External_Pk — RIFT P(k) from G_eff(Ψ), full Boltzmann transfer | PASS | |
| SIM89 | Planck TT Likelihood — approximate Gaussian TT vs Planck 2018 best-fit | PASS | Baseline check for SIM90 |
| SIM90 | Joint CMB+BAO Parameter Fit — Λ₀ = 0.003, m₀ = 0.001 H₀ | PASS | Resolves Phase 5b tension |
| SIM91 | Λ₀ Sensitivity Scan — G_eff, growth factor, σ₈, H(z), BAO χ² vs Λ₀ | PASS | Detectability threshold mapped |
| SIM92 | Non-linear Structure Growth — CLASS linear P(k) + HALOFIT Takahashi+2012 | PASS | |
| SIM93 | Bayesian Model Comparison — Savage-Dickey density ratio, RIFT vs ΛCDM | PASS | ln B ≈ −0.71 (inconclusive); LCDM nested at Λ₀=0 |
| SIM94 | DESI Year 1 BAO Refit — DESI Y1 BAO data refit | PASS | |
| SIM95 | CMB Polarization EE/TE — CLASS full Boltzmann | PASS | |
| SIM96 | RSD / f·σ₈ Growth Rate — redshift-space distortions | PASS | |
| SIM97 | Full Planck 2018 plikHM TTTEEE Likelihood | PASS | Real likelihood, not approximation |
| SIM98 | Joint plikHM CMB + BAO Parameter Fit | PASS | Δ(−2 ln L) = −0.013 vs ΛCDM |
| SIM107 | Solar System Viability — effective matter-Ψ coupling, fifth-force bounds | PASS | Cassini constraint satisfied |
| SIM108 | Gravitational Wave Speed — GW170817 c_T constraint | PASS | 6/6 tests pass |
| SIM109 | Dark Energy w(z) — full RIFT background w₀, wₐ | PASS / FAIL | PASS vs Planck; 3.6σ tension vs DESI |
| SIM110 | Slow-Roll Inflation — (nₛ, r) from locked action | FAIL | nₛ = 0.932 (need 0.9649); requires ξ ~ 10⁴ Starobinsky |

---

### Paper III — BAO Sound Horizon (`III_sound_horizon/`)

| Sim | Topic | Result | Notes |
|-----|-------|--------|-------|
| SIM101 | BAO Sound Horizon r_d from First Principles — CLASS Boltzmann integral | PASS | r_d = 147.2 Mpc matches Planck; RIFT δr_d/r_d < 0.3% |

---

### Paper IV — UV Structure (`IV_uv_structure/`)

| Sim | Topic | Result | Notes |
|-----|-------|--------|-------|
| SIM102 | One-Loop Scalar Self-Energy — UV structure and memory regulation | PASS | Memory damping regulates loop integrals |
| SIM104 | Full One-Loop UV Structure — all diagrams (Ψ self-energy, graviton, vertex) | PASS | UV finiteness at full one-loop confirmed |
| SIM105 | RG Flow of Λ₀ — beta function and running coupling | PASS | Λ₀ asymptotically free; dΛ₀/d ln μ < 0 |
| SIM106 | Two-Loop Graviton Sector — graviton self-energy and mixed diagrams | PASS | Final UV caveat closed |

---

### Paper V — Dark Matter (`V_dark_matter/`)

| Sim | Topic | Result | Notes |
|-----|-------|--------|-------|
| SIM99 | Galactic Rotation Curves — shooting method, quartic potential | FAIL | No (Ψ_bc, λ_gal) gives flat curves within cosmological constraints |
| SIM100 | 1+1D Time-Domain Galactic Evolution — tachyonic condensation test | FAIL | Condensate does not reproduce flat rotation curves |
| SIM103 | Vainshtein Screening — fifth-force screening mechanism at galactic scales | PASS | Screening works; but DM mechanism not viable in locked action |

---

### Phase 1 DESI Extension (post-master)

| Sim | Topic | Result | Notes |
|-----|-------|--------|-------|
| SIM111 | m₀ Scan for DESI w₀-wₐ Tension — full scan m₀ ∈ [0.001, 10] H₀ | FAIL | Min achievable DESI tension = **3.44σ** (structural obstruction) |

---

## Phase 2 — Unlocked Action

Action: bare Λ removed; DE from scalar potential; DM from β coupling.  
Current form: `S = ∫d⁴x√(-g) [(M_Pl²+2Λ₀Ψ²)/2 R - ½(∇Ψ)² + V(Ψ) + βΨ²L_m] + S_SM`  
Parameters: **Λ₀ = 0.003** (locked), λ and β derived from first principles.  
Papers: `paper/phase2/`

---

### P2 Paper I — Framework (`P2_I_framework/`)

| Sim | Topic | Result | Notes |
|-----|-------|--------|-------|
| SIM112 | λΨ⁴ Quintessence — no bare Λ, DE from plateau | FAIL | Always freezing quintessence (wₐ > 0); wrong sign for DESI |
| SIM113 | SSB Hilltop V_J = λ(Ψ²−v²)² — no bare Λ | PARTIAL | Thawing confirmed (wₐ < 0); DESI tension 3.6σ → 2.7σ; gap persists |
| SIM114 | βΨ²ρ_m condensate: galactic Ψ profile and rotation curves | FAIL | Trilemma: condensate needs β > 8πΛ₀=0.075 but then G_eff→0; β_needed(~10⁻⁴) << threshold; v=13.16 >> Ψ_max=1.29 (10×) structural |
| SIM115 | Gradient soliton: Ψ dips below Ψ_cosmo in dense regions (most-RIFT approach) | FAIL | G_eff safe (Ψ varies only 4×10⁻⁷!); ρ_DM ~10⁻⁵–10⁻¹⁴ × needed; H₀² suppression of λ_kpc structural — not tunable |
| SIM116 | Ψ-switch ξ condensate: separate ξ field with V(ξ,Ψ) = λ_ξ(ξ²−Ψ²/λ_ξ)²; ξ condenses in halo (Ψ≈Ψ_cosmo), absent in disk (Ψ≈0) | PLANNED | Energy scale set by λ_ξ independently of H₀²; tests Ψ→0 switching as halo geometry mechanism |
| SIM117 | Memory-of-memory ξ: □ξ + m_ξ²ξ = g∫G_R(x−x')Ψ²(x')d⁴x'; ξ sourced by time-integrated Ψ² history | PLANNED | Most RIFT philosophically; ρ_DM ∝ integrated baryon history → naturally traces galaxy age/mass; RIFT recursion at level 2 |

---

## Summary Statistics

| Phase | Sims | PASS | FAIL | PARTIAL | PLANNED |
|-------|------|------|------|---------|---------|
| Phase 1 (P-I)    | 6  | 5 | 0 | 1 | — |
| Phase 1 (P-II)   | 16 | 14 | 2 | 0 | — |
| Phase 1 (P-III)  | 1  | 1 | 0 | 0 | — |
| Phase 1 (P-IV)   | 4  | 4 | 0 | 0 | — |
| Phase 1 (P-V)    | 3  | 1 | 2 | 0 | — |
| Phase 1 extension| 1  | 0 | 1 | 0 | — |
| Phase 2          | 6  | 0 | 3 | 1 | 2 |
| **Total**        | **37** | **25** | **8** | **2** | **2** |

---

## Simulation Numbering Notes

- **SIM81** was skipped (reserved / not published).
- **SIM1–SIM79**: Pre-series exploratory runs (old `simulation_XX` format, pre-locked-action framework). Superseded by SIM80+ series. Key issues: SIM04/SIM11 READMEs swapped, SIM09 Ψ→−40 unphysical, SIM14 sweep broken.
- **SIM114+**: Phase 2 programme pending (DM sector, full cosmological tests, Solar System screening with β coupling).
