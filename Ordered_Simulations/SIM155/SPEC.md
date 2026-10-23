# SIM155 — SPEC (pre-registered 2026-10-23, before any fit was run)

Full-CMB test of the SIM154 survivors (ROADMAP Phase 3.4). The data, models, method and pass criteria are
those fixed in `ROADMAP.md` (commit `566a9a2`), repeated here:

- **Data:**
  - D1 = Planck 2018 (low-ℓ TT + low-ℓ EE + plik-lite TTTEEE + lensing; native Cobaya) + DESI DR2 BAO + Pantheon+
  - D2 = D1 with Union3 instead of Pantheon+
  - D3 = D1 with DES-Y5 instead of Pantheon+
- **Models** (hi_class v3.4 + `hi_class_model/cmstg_attractor.patch`, via `common/cobaya_hiclass.py`):
  - ΛCDM: ω_b, ω_cdm, H0, log A_s, n_s, τ (+ A_planck)
  - w0waCDM: + w0, wa (CLASS fluid with PPF)
  - CMSTG attractor family A: + Λ0 ∈ [0.02, 1], ξ ∈ [0, 0.3]; φ = 0 on the radiation-era attractor at a = 1e-10; bare Λ by closure
  - Shared: 1 × 0.06 eV neutrino (N_eff = 3.046); halofit; YHe = 0.2454 + 0.438(ω_b − 0.02238) + 0.16(S − 1)
- **Method:** Cobaya `minimize` (BOBYQA), 3 independent starting points per model and dataset. Report the best of the three.
- **CMSTG passes** if:
  - (a) Δχ² vs ΛCDM ≥ 6 on D1;
  - (b) at the best fit, Cassini/LLR (C2) and BBN (C3) from SIM154 still hold;
  - (c) Δχ² ≥ 6 on D2 or D3.
- **Always reported:** Δχ² and ΔAIC vs ΛCDM and vs w0waCDM. "Competitive with w0waCDM" means χ² within 3.
