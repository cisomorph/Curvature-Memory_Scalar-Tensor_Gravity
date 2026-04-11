# RIFT Paper v2 — Submission Guide
Generated: 2026-04-11 (updated from 2026-04-09)

---

## STEP 1: arXiv (submit first, before JCAP)

**URL:** https://arxiv.org/submit  
**File to upload:** `rift_paper_v2_arxiv.tar.gz` (400K)

### Metadata to paste in:

**Title:**
```
Recursive Intelligence-Field Theory (RIFT): A Lagrangian-Based Framework for Curvature Memory and Emergent Gravitation
```

**Authors:**
```
Christopher Robert Barrick Wilson
```

**Abstract** (plain text, no LaTeX math — paste this):
```
Current cosmological theory requires dark matter, dark energy, and inflation — none
confirmed in laboratory physics. This paper presents Recursive Intelligence-Field
Theory (RIFT), a covariant scalar-tensor field theory in which spacetime curvature
couples dynamically to a scalar field Psi that recursively encodes the history of
past curvature via a retarded Green's function (curvature memory). The theory is
derived from a single covariant action; all field equations follow by standard
variation without additional postulates.

We subject RIFT to fourteen independent tests (Sims 87--100). Key results:
(i) BAO full-covariance refit (BOSS+eBOSS+Ly-alpha): chi2/dof=1.22,
    Delta_chi2(RIFT-LCDM)=0;
(ii) CMB TT via full CLASS Boltzmann solver: RMS deviation=2.93% (l=2-1500),
    G_eff/G = 1 - 16 ppm at best-fit coupling;
(iii) joint CMB+BAO parameter fit: Delta_chi2=0, H0=67.59 km/s/Mpc,
    Omega_m=0.312, Lambda0~0.013 (numerically unconstrained; chi2 flat across [0,0.1]);
(iv) non-linear structure growth: |Delta_sigma8| < 0.007% at Lambda0=0.003;
(v) Bayesian model comparison (Savage-Dickey): ln B = -0.71 (Inconclusive; RIFT not excluded);
(vi) DESI Year 1 BAO refit: chi2/dof=1.36, Delta_chi2=0;
(vii) CMB polarization EE/TE via CLASS: RMS < 0.2% (EE), < 0.5% (TE);
(viii) RSD/f*sigma8 growth rate (9 surveys, z=0.07-1.48): chi2/dof=0.86, Delta_chi2=-0.14;
(ix) official Planck 2018 plikHM TTTEEE likelihood: Delta(-2lnL)=-0.013~0 at CMB optimum;
(x) galaxy rotation curves, static (Sim 99): locked action cannot replace dark matter
    -- Yukawa halo exponentially steeper than 1/r^2 isothermal needed for flat curves;
(xi) galaxy rotation curves, time-domain PDE (Sim 100): galactic curvature coupling
    ~5e-9 kpc^-2 is 5e7x weaker than m0^2; curvature condensation ruled out.
In all cosmological tests (Sims 87--98), RIFT at the joint best-fit is statistically
indistinguishable from LCDM at current precision. The coupling Lambda0 < 0.095 (95% CL)
is cosmologically inert. Galactic-scale predictions await a screening mechanism.
Falsifiable predictions are given for DESI Y5, Euclid, and CMB-S4.
```

**Primary category:** `astro-ph.CO`  
**Cross-list:** `gr-qc`

**Comments field:**
```
32 pages, 7 figures. Simulation code and data: https://github.com/cisomorph/rift-cosmology
```

**License:** Creative Commons Attribution 4.0 International (CC BY 4.0)

---

## STEP 2: JCAP Submission

**Journal:** Journal of Cosmology and Astroparticle Physics (JCAP)  
**Publisher:** IOP Publishing / SISSA  
**Submission portal:** https://www.sissa.it/app/jcap/ (or via https://iopscience.iop.org/journal/1475-7516)

**Files to upload:**
- `rift_paper_v2.tex` (main manuscript)
- `fig1_bao_residuals.pdf`
- `fig2_cmb_Cl.pdf`
- `fig3_lambda0_scan.pdf`
- `fig4_sigma8_scan.pdf`
- `fig5_desi_bao.pdf`
- `fig6_polarization.pdf`
- `fig7_rsd_fsigma8.pdf`
- Cover letter (paste text below)

**OR** upload the single `rift_paper_v2_arxiv.tar.gz` if the portal accepts archives.

**Article type:** Regular article  
**Section:** Modified gravity / Dark energy / Cosmological perturbations

### Cover Letter (paste into portal):

---

Dear JCAP Editors,

We submit for your consideration "Recursive Intelligence-Field Theory (RIFT):
A Lagrangian-Based Framework for Curvature Memory and Emergent Gravitation"
by Christopher Robert Barrick Wilson (Independent Researcher).

**Summary.** RIFT is a covariant scalar-tensor modified gravity theory in which
a non-minimally coupled scalar field Psi produces *curvature memory*: the geometry
at each spacetime event retains a causal, retarded integral over past curvature.
The complete Lagrangian is derived from first principles with no postulates beyond
the modified Einstein-Hilbert action; all field equations follow by standard
variation. The theory recovers GR exactly in the Psi -> 0 limit.

**Observational programme.** The paper presents fourteen independent numerical
tests against precision cosmological data, implemented as fully reproducible Python
simulations (code released at https://github.com/cisomorph/rift-cosmology):

1. BAO full-covariance refit (BOSS DR12 + eBOSS DR16 + Ly-alpha, 12x12 covariance):
   chi2/dof = 1.22, Delta_chi2(RIFT-LCDM) = 0.
2. CMB TT power spectrum via CLASS full Boltzmann code: RMS = 2.93%, G_eff/G = 1-16 ppm.
3. Planck TT likelihood: LCDM baseline chi2 = 6.20 (validated).
4. Joint CMB+BAO parameter fit: Delta_chi2 = 0; H0=67.59 km/s/Mpc, Omega_m=0.312.
5. Lambda0 sensitivity scan: BAO chi2 flat across Lambda0 in [0, 0.1].
6. Non-linear structure growth (HALOFIT): |Delta_sigma8| < 0.007% at Lambda0=0.003.
7. Bayesian model comparison (Savage-Dickey): ln B = -0.71 (Inconclusive).
8. DESI Year 1 BAO refit: chi2/dof = 1.36, Delta_chi2 = 0.
9. CMB polarization EE/TE via CLASS: RMS < 0.2% (EE), < 0.5% (TE).
10. RSD/f*sigma8 growth rate (6dFGRS, BOSS, eBOSS, 9 data points): chi2/dof = 0.86.
11. Official Planck 2018 plikHM TTTEEE likelihood (via clipy): Delta(-2lnL) = -0.013~0
    at joint CMB+BAO optimum.
12. Galaxy rotation curves, static test (Sim 99, NGC 3198 + NGC 6503): the locked
    action at Lambda0=0.003 cannot replace dark matter. The Yukawa halo profile is
    exponentially steeper than the 1/r^2 isothermal needed for flat curves.
13. Galaxy rotation curves, time-domain PDE (Sim 100): galactic curvature coupling
    2*Lambda0*|R_gal| ~ 5e-9 kpc^-2 is 5e7x weaker than m0^2; curvature condensation
    at galactic densities is ruled out; lambda<0 self-coupling causes tachyonic blow-up.

In all cosmological tests (Sims 87--98), RIFT at the joint best-fit is statistically
indistinguishable from LCDM. The coupling parameter Lambda0 is bounded to
Lambda0 < 0.095 (95% CL) and is cosmologically inert at current observational
precision. The paper explicitly documents where the theory currently fails
(rotation curves) and what mechanism is needed to address galactic scales.

**Why JCAP.** This paper sits at the intersection of modified gravity theory and
precision cosmological data analysis — exactly the scope of JCAP. All results are
fully reproducible, the code is open-source, and the paper includes falsifiable
predictions for DESI Y5, Euclid, and CMB-S4.

The manuscript has not been submitted elsewhere. There are no conflicts of interest.

Sincerely,  
Christopher Robert Barrick Wilson  
C.Isomorph@gmail.com  
https://github.com/cisomorph/rift-cosmology

---

## STEP 3: After arXiv appears (~1 business day)

1. Note your arXiv ID (e.g., 2604.XXXXX)
2. Add it to the JCAP submission if not yet submitted
3. Update the paper date: change `\date{July 2025 (v2, revised \today)}` to
   `\date{April 2026 (v2)}` and push to GitHub

## STEP 4: Add arXiv DOI to GitHub repo description

Once you have the arXiv ID, run:
```bash
gh repo edit cisomorph/rift-cosmology --description "RIFT paper (arXiv:2604.XXXXX) — curvature memory modified gravity. 14 cosmological tests."
```

---

## Follow-on Research Programme (post-submission)

Three open problems identified during v2 preparation, in order of scientific leverage:

### Paper III — r_d from first principles
Derive the BAO sound horizon from RIFT field equations. Currently r_d is a free
parameter; a RIFT-specific prediction of r_d would be the first genuinely distinct
result over LCDM. Requires: RIFT perturbation equations around Psi(t) background,
effective sound speed c_s(Psi), integrate r_d = integral(c_s dt/a).

### Paper IV — UV loop calculation
Compute the one-loop correction to the RIFT propagator with the Lambda(Psi)R*Psi
interaction vertex. UV finiteness is currently a conjecture. This either proves or
kills the "no renormalization" claim. Hard calculation — may require dimensional
regularization + heat kernel methods.

### Paper V — Dark matter mechanism
SIM99/100 ruled out G_eff and curvature condensation. Derive a screening or
higher-order stabilization mechanism from the action that could sustain a
galactic-scale halo. Candidates: (a) non-minimal kinetic term X(Psi)(dPsi)^2,
(b) Vainshtein-type screening, (c) topological soliton solutions.
