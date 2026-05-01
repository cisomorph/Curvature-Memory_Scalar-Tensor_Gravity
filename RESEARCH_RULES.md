# Standing Rules for the CMSTG Repository

These rules govern all simulation, paper, and commit work. They derive from Phase 1, Phase 2, and Phase 3 experience and are not negotiable without explicit discussion.

---

## 1. Project context

This is the CMSTG (Curvature-Memory Scalar-Tensor Gravity) simulation repository, maintained by Christopher R. B. Wilson. CMSTG is a scalar-tensor extension of GR with a non-minimal coupling Λ(Ψ)R and a retarded memory kernel. The Phase 1 canonical action is locked:

- Action: S = ∫d⁴x √-g [(1/16πG + Λ(Ψ))R + ½(∇Ψ)² − ½m(Ψ)²Ψ² − V(Ψ)] + S_matter
- Λ₀ = 0.003 (locked)
- Ψ̄ = 2.62 M_Pl
- F_eff = (1 + 2Λ₀Ψ̄²)/2 = 0.521
- Phase 1 passes: BAO, CMB TT/EE/TE, RSD, GW speed, Solar System, UV finiteness through two loops
- Phase 1 residual: 2.77σ tension with DESI Y1 BAO (structural)

Phase 2 (SIM112–130) exhausted matter-sector extensions. Phase 3 (SIM131–136) exhausted curvature-sourced scalar extensions and produced a structural no-go theorem. Phase 4 (SIM137 onward) is the current programme.

---

## 2. Core theoretical discipline

### 2.1 Foundations before phases
Do not advance to later simulation phases before the earlier phase is closed. If a diagnostic sim (Tier 1) changes what the mechanism sim (Tier 2) should target, update the spec before running.

### 2.2 Derived vs. phenomenological
Always distinguish:
- **Derived results:** follow from variation of the Lagrangian
- **Phenomenological parametrizations:** functional forms imposed by hand (e.g., f_ψ(a), running Λ₀(a))

Flag phenomenological inputs explicitly in sim output JSON and in any paper text. Do not let phenomenological fits be cited later as derivations.

### 2.3 Every new Lagrangian term must be checked against
Before any mechanism simulation reports a PASS, verify:

1. **GR recovery:** Λ(0) = 0, V(0) = 0, G_eff > 0
2. **Graviton sector:** c_T = c (GW170817 bound); m_g = 0 (Ward identity Π_hh(0) = 0)
3. **No tachyon:** ω² = k² + m_eff² > 0 for all physical modes
4. **UV behaviour:** no new quartic divergences that the memory kernel does not suppress
5. **Causality:** retarded Green's function structure preserved

A sim that passes observational criteria but fails any of the above is a FAIL, not a PASS.

### 2.4 The Sim 8 lesson
V(Ψ) = αΨ² exp(−βΨ²) did not produce monotone decaying halo solutions without a matter source term. The lesson: potential shape alone does not guarantee the desired field profile. Check that the EOM actually admits the solution class you're looking for before scanning parameters.

---

## 3. Simulation conventions

### 3.1 Directory structure
```
cmstg/
├── RESEARCH_RULES.md            # this file
├── PHASE4_ROADMAP.md            # current roadmap
├── sims/
│   ├── SIM137/
│   │   ├── SPEC.md              # pre-run: action, criteria
│   │   ├── run.py               # executable simulation
│   │   ├── output.json          # numerical results
│   │   ├── RESULT.md            # post-run: verdict, diagnosis
│   │   └── figures/
│   ├── SIM138/
│   └── ...
├── papers/
└── shared/                      # common utilities, solvers
```

### 3.2 Output JSON schema
Every sim must produce `output.json` with at minimum:
```json
{
  "sim_id": "SIM137",
  "timestamp": "ISO8601",
  "action_spec": "Phase 1 canonical" | "<full Lagrangian description>",
  "parameters": { ... },
  "observational_targets": { "dataset": "...", "chi2": ..., "dof": ... },
  "theoretical_checks": {
    "gr_recovery": true,
    "c_T_eq_c": true,
    "no_tachyon": true,
    "ward_identity": true,
    "uv_finite": true
  },
  "verdict": "PASS" | "PARTIAL" | "FAIL",
  "failure_mode": "<structural diagnosis if FAIL or PARTIAL>",
  "derived_vs_phenom": { "<param>": "derived" | "phenomenological" }
}
```

### 3.3 Success criteria discipline
Each sim spec states explicit, verifiable success criteria. Loop until all are verified or until a structural obstruction is diagnosed. Do not report PASS on partial evidence.

### 3.4 Commit after every sim
After each sim completes (PASS, PARTIAL, or FAIL), immediately:
1. Commit the `SIM<N>/` directory with output.json, RESULT.md, figures
2. Update `PHASE4_ROADMAP.md` with the verdict and any scope changes
3. Push to GitHub before starting the next sim

This preserves the paper trail that made Phase 2 and Phase 3 documents strong.

---

## 4. Coding style

### 4.1 Simplicity first
- Minimum code that solves the problem
- No speculative features, no abstractions for single-use code
- Match existing style in `shared/` rather than refactoring
- If a sim script exceeds 300 lines, factor the common numerical machinery into `shared/`

### 4.2 Numerical verification
- Every new solver: compare against a known analytic limit
- Every new likelihood: verify against the Phase 1 canonical parameters reproduce SIM90/SIM98 χ² within 0.01
- Every new ODE integrator: convergence test with step size halving

### 4.3 External dependencies
- Prefer CLASS for Boltzmann work (used in Phase 1)
- Prefer `clipy` for Planck plikHM likelihood (used in SIM97/98)
- Use full covariance matrices for BAO (BOSS+eBOSS+Lyα, DESI Y1)
- When in doubt about a citation or numerical value, verify with a web search before using it

---

## 5. Writing and papers

### 5.1 When a phase closes, write the conclusions paper
Follow the structure of Papers VI and VII:
- Introduction + motivation
- Per-sim results with actions, parameters, outcomes
- Summary table with verdicts and Δχ²
- Structural diagnosis (no-go theorem if applicable)
- Directions for the next phase
- Canonical parameter set if changed

### 5.2 No overreach
- "CMSTG" is a scalar-tensor extension of GR with a retarded memory kernel. Not "modified gravity," not a TOE candidate in the paper text.
- Speculative material (e.g., SM coupling appendix in the master paper) must be explicitly labeled speculative.
- The Bayesian evidence is ln B = -0.71 (Inconclusive). Do not claim CMSTG is preferred over ΛCDM; claim it is not excluded.

### 5.3 Publication strategy
- Zenodo with CC BY 4.0 for timestamping and priority protection
- Peer review venues: Phys. Rev. D, JCAP for cosmology sims; Class. Quantum Grav. for the canonical action paper
