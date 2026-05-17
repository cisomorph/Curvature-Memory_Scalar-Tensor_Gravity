# Framework

CMSTG is a covariant scalar-tensor theory in which a scalar field Ψ couples non-minimally to the Ricci scalar via a *curvature-memory* kernel — a retarded non-local operator that gives Ψ an effective mass set by the integrated curvature history of the universe.

---

## The Action

The Phase 1 canonical action (dimension-4, covariant, locked after SIM98):

```
S = ∫ d⁴x √−g  { (M_Pl²/2)(1 + 2Λ₀Ψ²)R
                 − (1/2)(∂Ψ)²
                 − V(Ψ)
                 − (Λ₀/M²) Ψ ∫ K(x,x') R(x') d⁴x' √−g'  }
    + S_matter[g_μν, ψ_m]
```

**Key parameters:**
- `Λ₀ = 0.003 M_Pl⁻²` — non-minimal coupling (locked by joint CMB+BAO fit)
- `K(x,x')` — retarded memory kernel (Yukawa profile; causal by construction)
- `V(Ψ)` — scalar potential; flat near the attractor

---

## Effective Planck Mass

The non-minimal coupling generates a running effective Newton constant:

```
F_eff(Ψ) = (1 + 2Λ₀Ψ²) / 2
```

At the cosmological attractor Ψ̄ = 2.62 M_Pl:

```
F_eff = 0.521   (locked)
G_eff / G = 1 − 15 ppm  at z = 0
```

---

## Field Equations

**Scalar equation** (Ψ sourced by curvature and memory):
```
□Ψ = −∂V/∂Ψ + 2Λ₀ΨR + (Λ₀/M²) ∫ K(x,x') R(x') d⁴x'
```
The curvature-sourcing term drives Ψ toward its attractor. This is the origin of Theorem 1: any scalar sourced by R > 0 grows monotonically, increasing F_eff and suppressing H(z) at DESI redshifts — the wrong direction.

**Modified Friedmann equation:**
```
3 F_eff H² = ρ_m + ρ_r + ρ_Ψ + ρ_memory
```

**GW speed:**
```
c_T² = 1   (exact, by construction)
```

---

## Memory Kernel

| Property | Value / result |
|----------|---------------|
| Profile | Yukawa-type: K ~ exp(−r/λ) / r |
| Causality | Strictly retarded — no advanced contributions (SIM82 PASS) |
| UV behaviour | Resummed propagator; two-loop UV finite (Paper III / SIM102–106) |
| RG flow | Negative β-function → Λ₀ → 0 in UV (asymptotic-freedom-like) |
| IR cutoff | Set by memory scale M; constrained by BAO rₛ (SIM101 PASS) |

---

## Canonical Parameters

| Symbol | Value | Unit | Status | Note |
|--------|-------|------|--------|------|
| Λ₀ | 0.003 | M_Pl⁻² | locked | Non-minimal coupling |
| Ψ̄ | 2.62 | M_Pl | attractor | Cosmological VEV |
| F_eff | 0.521 | — | locked | (1+2Λ₀Ψ̄²)/2 |
| H₀ | 67.59 | km/s/Mpc | best-fit | Joint CMB+BAO |
| σ_DESI | 2.77 | σ | structural floor | DESI Y1 BAO residual |
| ln B | −0.71 | — | inconclusive | Bayes factor vs ΛCDM |

---

## Phase 1 Observational Tests (SIM80–SIM111)

| Sector | Test | SIM | Result |
|--------|------|-----|--------|
| BAO | BOSS + eBOSS + Lyα full covariance | SIM87 | PASS |
| BAO | DESI Y1 w(z) | SIM109–110 | FAIL — 2.77σ floor |
| CMB | Planck acoustic-peak shift θ* | SIM89 | PASS |
| CMB | plikHM TTTEEE alone | SIM97 | FAIL Δχ²=+7.0 |
| CMB | Joint plikHM + BAO refit | SIM98 | PASS Δ(−2lnL)=−0.013 |
| GW | Speed c_T = c | SIM88 | PASS |
| Solar System | G_eff / G PPN | SIM88, 107 | PASS |
| RSD | f σ₈ shape | SIM91 | PASS |
| σ₈ | Lensing + clustering | SIM92 | PASS |
| BBN | Neff, Yp | SIM93 | PASS |
| Neutrinos | Σmν, Neff | SIM94 | PASS |
| Ly-α | Power spectrum | SIM95 | PASS |
| UV | Two-loop renormalizability | SIM102–106 | PASS |
| UV | Ward identity | SIM103 | PASS |
| DM | Ψ-sector rotation curves | SIM99–100 | FAIL → Paper IV |
| DESI | m₀ scan, structural floor | SIM111 | 3.44σ minimum |

---

## Development Phases

| Phase | Sims | Focus | Outcome |
|-------|------|-------|---------|
| 1 | SIM80–111 | Locked canonical action; 22 obs. tests | Canonical declared; 2.77σ floor identified |
| 2 | SIM112–130 | Matter sector; χ-DM dark matter | χ-DM passes SPARC 161-gal; floor persists |
| 3 | SIM131–136 | Curvature-sourced scalar extensions | All 6 FAIL → Theorem 1 |
| 4 | SIM137–144 | Late-time mechanism space | All FAIL → Theorem 2; Tier 2 closed |
| 5 | SIM147–150 | Ψ_pre and Ψ_local loophole probes | Both loopholes triply closed |
| A | SIM151–153 | Pre-recombination options B & C | Option B FAIL; Option C HALT |

See [Simulation Programme](Simulation-Programme) for full sim tables.
