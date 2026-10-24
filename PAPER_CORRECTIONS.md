# Paper and README corrections implied by the audit (ROADMAP Phase 6, draft 2026-10-23)

These are the factual statements the audited reruns contradict. Each line gives the claim, what the reruns show, and
the evidence. Wording of the papers themselves is the author's call; nothing in `Papers/` has been edited.

## README / Paper II (framework, Phase 1)

| claim | what the audited reruns show | evidence |
|---|---|---|
| Phase 1 passes Solar-System tests (SIM107 PASS) | Fails Cassini by 41× at the canonical Ψ̄ = 2.62 (\|γ−1\| = 9.5e-4); no screening. SIM150 also found 9.5e-4. | SIM107, SIM150 docstring |
| "2.77σ DESI structural floor" | Not a property of CMSTG. It came from a non-DESI H(z) table, √(χ²/N), a truncated r_s and a mislabelled H0. With Ψ frozen, F0 is exactly degenerate with H0 for BAO + θ*. ΛCDM itself is 2.1σ (DR1) / 3.1σ (DR2) from DESI in this setup. | SIM121C, SIM122, B2–B6, B11 |
| "Dark energy emerges … without a bare cosmological constant" | Every audited model uses Λ_bare (or a potential that acts as one). | THEORY_NOTES 4.3 |
| "G_eff/G = 1 − 15 ppm at z = 0" | For F_eff = 0.521, G_eff/G = 0.5/F = 0.96 (−4%). | SIM121C notes |
| "Next decisive test: DESI Y3 (~2027)" | DESI DR2 (three years of data) was published in March 2025; every audited rerun includes it. | `common/data/` |
| H0 tension | CMSTG lowers the true H0 (65.9–66.6), worsening the SH0ES tension (5.6–6.2σ vs 5.0σ). | SIM121B |
| Phase 1 unlocked Ψ "negligible" (SIM130) | It improves on ΛCDM (Δχ² = +5.4 DR1 / +13.9 DR2), but Ψ0 ≈ 2.88 fails Cassini (49×). | SIM130 |

## Paper I (no-go theorems)

| claim | what the audited reruns show | evidence |
|---|---|---|
| Theorem 2: the DESI fix and θ* are mutually exclusive for late-time modifications ("θ* 63σ") | The 63σ came from holding H0 and r_s fixed. With H0 solved from θ*, w0waCDM and several CMSTG variants improve DESI + θ* (Δχ² ≈ +5 … +9.5 DR1). The theorem's premise (raise H(z) everywhere in [0, z_rec]) is not what the DESI preference requires. | SIM141, SIM121C |
| Theorem 1: curvature-sourced growth of F is "the wrong direction" | Also assumed a fixed H0. Growing F improves the fit once H0 floats (SIM130, SIM123, SIM135). | SIM130/123/135 |
| SIM131–144 all FAIL; SIM144 establishes completeness | Several flip to cosmology PASS (SIM130, 132, 134, 135, 141–144). What excludes most of them is local gravity (Cassini, LLR, GW170817), not DESI/CMB. The attractor models SIM132/134 survive local gravity. | AUDIT_NOTES §0 |
| Any statement that CMSTG explains the DESI preference (added 2026-10-24) | The background-level gains above do not survive a full fit. With hi_class + Planck 2018 + DESI DR2 + SNe (SIM155): **(1)** the best local-gravity-safe attractor model (family A) improves on ΛCDM by only Δχ² = +4.65 (D1), +4.50 (Union3) and +3.19 (DES-Y5). **(2)** Even that is only at a change of G since BBN (8%) that BBN excludes (Y_p +3.3σ); where BBN allows, the gain is +3.25 for two extra parameters (ΔAIC +0.75). **(3)** w0waCDM gains +9.2 … +19.6 on the same data. **(4)** The model's change of G is complete before recombination, so it cannot produce the late-time signal DESI + SNe prefer. **(5)** It raises S8 by 2.5–7.5%. A paper should report this test. | SIM155 RESULT.md |
| Specific model equations | Friedmann without a GR limit (Phase 3); doubled F'Ψ' and kinetic terms; Gauss–Bonnet term missing (SIM133); Horndeski energy and EOM wrong, c_T unchecked (SIM136); G3 term (SIM142); force/energy factor 2 (SIM143); matter-coupling sign (SIM144). | AUDIT_NOTES B9–B10, sections 6 |

## Paper IV (galactic)
Not audited. SIM150's Phase-5 INCOMPATIBLE verdict stands on independent grounds.

## Paper III (UV)
Not audited (THEORY_NOTES 4.4).
