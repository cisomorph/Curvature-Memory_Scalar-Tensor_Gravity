# CMSTG attractor model in hi_class (ROADMAP Phase 3)

`cmstg_attractor.patch` adds the CMSTG attractor models (SIM132/SIM134/SIM154) to
[hi_class](https://github.com/miguelzuma/hi_class_public) v3.4.0.1 (commit `ece81a57`) as a covariant
Horndeski model:

    G2 = X − V0,    V0 = 3 p0 H0²   (bare cosmological constant, tuned for closure)
    G4 = F(φ):  family 0  F = ½(1 − 2Λ0φ²) + ξφ        (subtractive, SIM132)
                family 1  F = ½e^{−2Λ0φ²} + ξφ         (exponential, SIM134)

    gravity_model  = cmstg_attractor
    parameters_smg = p0, Lambda0, xi, family, phi_ini, phi_prime_ini

`phi_prime_ini = 0` puts φ' on its radiation-era attractor, φ' = 3a²ρ_m F′ / [4F aH (1 + 3F′²/F)],
refined by hi_class's secant step. Starting at rest makes the stiff background solver crawl, the same
issue hi_class documents for Brans–Dicke.

## Build
    git clone https://github.com/miguelzuma/hi_class_public.git hi_class
    cd hi_class && git checkout ece81a57 && git am ../cmstg_attractor.patch
    make class

## Required settings (see `A_k6.ini`)
    a_ini_over_a_today_default = 1.e-10
    a_ini_test_qs_smg = 1.e-10
    pert_ic_ini_z_ref_smg = 1.e9
    kineticity_safe_smg = 1e-6
    cs2_safe_smg = 1e-5

The model's scalar has no kinetic energy at the start and turns around near its attractor (z ≈ 7.4 for
the SIM154 survivor). At those points the code's kinetic matrix and c_s² formula become 0/0, so a small
kineticity regulator is needed, and a tiny negative c_s² noise (≈ −1e-7) must be tolerated.

## Validation (2026-10-23)
| test | result |
|---|---|
| GR limit (Λ0 = ξ = 0): background vs ΛCDM | identical 100θ_s (1.043001); z_drag to 6e-9 |
| GR limit: lensed C_ℓ vs ΛCDM | 2e-4 (TT), 4e-4 (EE), 9e-4 (TE) of peak; independent of the regulator, i.e. hi_class's scalar-sector numerics |
| regulator convergence (survivor): κ = 1e-5 vs 1e-6 | 4e-6 (TT), 8e-6 (EE), 1e-5 (TE) |
| regulator convergence: κ = 1e-4 vs 1e-6 | 4e-5 … 8e-5 |
| survivor vs ΛCDM at the same H0 (the physical effect) | 2.5e-3 (TT), 6.6e-3 (EE), 1.1e-2 (TE) |
| background vs the project's own `st_solver` (survivor, z_ini = 1e8) | H/H_ΛCDM to 2e-5; F(z) to 1e-5 (z = 0 … 1090) |

Runtime is about 5 s per full C_ℓ evaluation on one core.
