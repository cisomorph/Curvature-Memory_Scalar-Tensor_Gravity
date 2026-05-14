# SIM147 RESULT — Curvature-Memory Kernel Decay Characterization

**Date:** 2026-05-14  
**Phase:** 5 (Phase 5 gate calculation)  
**Gate for:** SIM148 (Ψ_pre conjecture), SIM149 (Ψ_local conjecture)

---

## Gate Decision: RED

> **τ_mem ≪ 1 Gyr and w(13.4 Gyr) < 10⁻⁴.**  
> Both Ψ_pre (SIM148) and Ψ_local (SIM149) conjectures are effectively dead at the
> locked k_m = 10 Mpc⁻¹. SIM148 and SIM149 should run for completeness but results
> must be flagged as conditional on unphysically large initial Ψ amplitudes.

---

## Summary of Results

### Analytic kernel derivation

The CMSTG retarded Green's function in k-space (massless limit, m₀ ≪ k_m):

    G̃_R(Δt; k) = θ(Δt) · sin(ω_k Δt)/ω_k · e^{−k²/k_m²}

The position-space kernel at coincident points (x = x') is:

    K(Δt) = (1/2π²) ∫₀^∞ dk k sin(kΔt) e^{−k²/k_m²}
           = (k_m³ / 8π^{3/2}) · Δt · exp(−(k_m Δt)²/4)

This is a **Gaussian envelope in Δt**, not exponential or power-law. The Fourier
transform of a Gaussian in k is a Gaussian in t.

**The mass correction is negligible:** m₀/k_m ≈ 2×10⁻⁸, so the massless approximation
is exact to O((m₀/k_m)²) ≈ 10⁻¹⁶.

### Locked parameters

| Parameter | Value |
|---|---|
| Λ₀ | 0.003 |
| m₀ | 2.245×10⁻⁷ Mpc⁻¹ = 10⁻³ H₀ |
| k_m (locked, SIM102) | 10 Mpc⁻¹ |
| τ_mem = 2/k_m | 0.2 Mpc = **2.05×10⁻⁴ Gyr ≈ 205,000 years** |
| Δt_peak = √2/k_m | 0.141 Mpc = 1.45×10⁻⁴ Gyr |

### Fit comparison — locked k_m = 10 Mpc⁻¹

| Form | τ_mem [Gyr] | α | R² | Note |
|---|---|---|---|---|
| Exponential | 1.02×10⁻⁴ | — | −0.003 | Best among three candidates |
| Exp × power law | 1.02×10⁻⁴ | 0.5 | −0.143 | |
| Power law | — | 1212 | −4.9×10⁵ | |

> **None of the three candidate forms fits well (all R² < 0).**  
> This is expected: the kernel has a Gaussian shape, which is not among the
> three candidates. The correct form is derived analytically above.

### Memory weights at reference horizons

All weights are zero to all floating-point precision (double underflow).

At Δt = 10 Gyr in Mpc: 10 × 977.8 = 9778 Mpc.  
Gaussian argument: (k_m × Δt)²/4 = (10 × 9778)²/4 = 2.39×10⁹.  
→ K(10 Gyr) ~ exp(−2.39×10⁹) ≈ 0 (not representable in IEEE 754 double).

| k_m [Mpc⁻¹] | τ_mem [Gyr] | w(10 Gyr) | w(13.4 Gyr) | w(13.8 Gyr) | Gate |
|---|---|---|---|---|---|
| 1 | 2.05×10⁻³ | 0 | 0 | 0 | RED |
| 3 | 6.82×10⁻⁴ | 0 | 0 | 0 | RED |
| **10** (locked) | **2.05×10⁻⁴** | **0** | **0** | **0** | **RED** |
| 30 | 6.82×10⁻⁵ | 0 | 0 | 0 | RED |
| 100 | 2.05×10⁻⁵ | 0 | 0 | 0 | RED |

**Entire sweep is RED.** τ_mem = 2/k_m; all five k_m values give τ_mem ≪ 1 Gyr.

### Parameter dominance

k_m dominates completely: m₀/k_m ≈ 2×10⁻⁸ for the locked k_m = 10 Mpc⁻¹.
The decay timescale is set entirely by k_m. To reach τ_mem = 10 Gyr would require
k_m ≈ 2×10⁻⁴ Mpc⁻¹ — five orders of magnitude below the SIM102 value.

### Spatial coherence

| Quantity | Value |
|---|---|
| ℓ_coh = 1/k_m | 0.100 Mpc |
| R_vir (L* galaxy) | ~0.100 Mpc |
| ℓ_coh / R_vir | 1.00 |
| Kernel uniform across L*? | **Marginally yes** (ℓ_coh = R_vir) |

SIM149 can treat the kernel as approximately spatially uniform across an L* galaxy
at the locked k_m. For smaller k_m the kernel is clearly uniform; for larger k_m it
is not. This finding is less relevant given the temporal decay rules out SIM149.

---

## Interpretation

The CMSTG memory kernel decays with a characteristic Gaussian timescale
τ_mem = 2/k_m determined entirely by the memory scale k_m. For the locked
k_m = 10 Mpc⁻¹ (established by SIM102), τ_mem ≈ 205,000 years — six
orders of magnitude shorter than the galactic formation timescale (~1 Gyr)
and eight orders of magnitude shorter than the recombination horizon (13.4 Gyr).

The memory integral Ψ(x) = ∫d⁴x' G_R(x,x') S[R(x')] accumulates curvature
contributions only within a time window ~τ_mem ≈ 200,000 years around each
event. The Ψ_pre conjecture (memory of curvature at the Big Bang surviving
to recombination) requires temporal support across 13.4 Gyr — a factor
~6.5×10⁴ beyond τ_mem. The Ψ_local conjecture (10 Gyr of halo curvature
accumulating) requires temporal support ~4.9×10⁴ × τ_mem.

Neither conjecture can hold at the locked parameters unless Ψ_pre or Ψ_local
are fine-tuned to unphysically large amplitudes to compensate.

**Physical origin:** The Gaussian memory damping e^{−k²/k_m²} is a UV
regulator that suppresses high-k modes. Its Fourier transform in time is
a Gaussian in Δt, which is the fastest possible decay (super-exponential).
This is the CMSTG design choice: memory is *local in time* at the scale
1/k_m. The theory was designed to give instantaneous-looking cosmological
sources at Hubble scales; it was not designed to retain memory on Gyr
timescales.

---

## Handoff to SIM148 and SIM149

### Artifacts saved to `sims/sim147_output/`

| File | Contents |
|---|---|
| `kernel_sim147.pkl` | Dict with `dt_Mpc`, `K_norm` arrays; formula string; conversion factors |
| `sim147_metadata.json` | All key outputs: τ_mem, gate, l_coh, sweep table |
| `sim147_kernel_decay.pdf` | K(Δt) plot and τ_mem vs k_m plot |

### Loading the kernel (Python)

```python
import pickle, numpy as np
from scipy.interpolate import interp1d

with open('sims/sim147_output/kernel_sim147.pkl', 'rb') as f:
    kd = pickle.load(f)

# Analytic formula (exact): K in Mpc units
km  = kd['km_locked_Mpcinv']  # 10.0 Mpc⁻¹
def K_func(dt_Mpc):
    return dt_Mpc * np.exp(-0.25 * (km * dt_Mpc)**2)

# Or use interpolation table:
K_interp = interp1d(kd['dt_Mpc'], kd['K_norm'], bounds_error=False, fill_value=0.0)

# Convert seconds to Mpc:
dt_Mpc = dt_seconds * kd['s_to_Mpc']
K_val  = K_interp(dt_Mpc)
```

### JSON metadata keys

```json
{
  "tau_mem_Gyr": 2.045e-4,
  "tau_mem_Mpc": 0.2,
  "l_coh_Mpc": 0.1,
  "k_m_locked": 10.0,
  "w_13.4Gyr": 0.0,
  "w_10Gyr": 0.0,
  "gate": "RED",
  "km_for_tau_1Gyr_Mpcinv": 2.044e-3,
  "km_for_tau_10Gyr_Mpcinv": 2.044e-4
}
```

---

## Gate Summary for SIM148 / SIM149

**GATE: RED**

- τ_mem ≈ 2.05×10⁻⁴ Gyr at locked k_m = 10 Mpc⁻¹
- w(10 Gyr) = w(13.4 Gyr) = w(13.8 Gyr) = 0.000 (double underflow)
- Kernel decays as Gaussian, not power law — no long-range temporal tail
- Both Ψ_pre and Ψ_local require Ψ amplitudes ~exp(10⁹) × physical values
- Conjectures are effectively dead at locked parameters
- SIM148 and SIM149 should proceed for formal completeness with this flag

To resurrect conjectures would require k_m ≤ 2×10⁻³ Mpc⁻¹ (for YELLOW gate).
This is ~5000× below the SIM102 value and inconsistent with UV finiteness
(SIM104–106). Phase 5 should not treat this as a free parameter.
