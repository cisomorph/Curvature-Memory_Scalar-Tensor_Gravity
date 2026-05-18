// app.jsx — main CMSTG site; depends on data.js and tweaks-panel.jsx

const D   = window.CMSTG_DATA;
const Ctx = window.TweaksContext;

const REPO_BASE = D.meta.repo + "/tree/main/";

function simRepoUrl(id) {
  const match = id.match(/SIM\d+[a-zA-Z]*/);
  return match ? REPO_BASE + "Ordered_Simulations/" + match[0] : null;
}
function paperRepoUrl(p) {
  return p.repoPath ? REPO_BASE + p.repoPath : null;
}

// ── Tooltip content dictionaries ──────────────────────────────
const PARAM_TIPS = {
  "\\Lambda_0":
    "Non-minimal coupling constant (0.003 M_Pl⁻²). Locked by the UV fixed point — the unique value at which the β-function vanishes and two-loop finiteness holds.",
  "\\bar\\Psi":
    "Cosmological VEV of the curvature-memory scalar Ψ. Any FRW background evolves to this attractor regardless of initial conditions, setting the amplitude of all memory effects.",
  "F_\\mathrm{eff}":
    "Effective gravitational coupling = (1 + 2Λ₀Ψ̄²)/2 = 0.521. Determines G_eff/G — the local gravitational constant differs from Newton's G by only 15 parts per million at z = 0.",
  "H_0":
    "Hubble constant from joint CMB + BAO fit (67.59 km/s/Mpc). Consistent with Planck 2018. The Hubble tension with local distance-ladder measurements (~73 km/s/Mpc) persists independently of CMSTG.",
  "\\sigma_\\mathrm{DESI}":
    "Irreducible 2.77σ tension with DESI Y1 BAO data. Theorems 1 & 2 prove this cannot be reduced by any late-time modification of the CMSTG action — it is a structural floor, not a fit parameter.",
  "\\ln B":
    "Log Bayes factor vs ΛCDM (−0.71). On the Jeffreys scale, |ln B| < 1 is 'inconclusive' — CMSTG is statistically indistinguishable from ΛCDM at current data precision.",
};

const STATUS_TIPS = {
  "locked":
    "Value fixed by UV fixed-point analysis. Cannot be varied without violating two-loop finiteness or the negative β-function that drives Λ₀ → 0 in the UV.",
  "attractor":
    "Not a free parameter — FRW dynamics drive Ψ to this value regardless of initial conditions. Changing it would require initial conditions far outside the attractor basin.",
  "best-fit":
    "Best-fit value from joint CMB + BAO parameter estimation. Has measurement uncertainty; the value is well-constrained by current data.",
  "structural floor":
    "Not a fitted parameter but the minimum achievable tension given the locked action. Theorems 1 & 2 prove this cannot be reduced by late-time modifications.",
  "inconclusive":
    "Bayes factor magnitude below the threshold for decisive model selection. The Jeffreys scale requires |ln B| > 5 for 'strong' evidence.",
};

const VERDICT_TIPS = {
  "PASS":
    "Meets all observational constraints within acceptance thresholds.",
  "FAIL":
    "Statistically significant conflict with one or more observational datasets.",
  "PARTIAL":
    "Satisfies some but not all constraints; the failure mode is potentially mitigatable.",
  "STRUCTURAL":
    "Failure is a consequence of theoretical architecture — persists across all parameter values.",
  "STRUCTURAL_PATTERN":
    "Multiple simulations share the same failure mode, confirming a systematic structural origin rather than an isolated numerical result.",
  "DISTRIBUTED":
    "The observational residual is spread smoothly across redshift bins, ruling out a localised dataset error and confirming a structural floor.",
  "SHAPE_FIXABLE":
    "The failure is in overall normalisation/shape and could in principle be corrected by parameter adjustment — but Theorems 1 & 2 prevent this in practice.",
  "PREDICTED FAIL":
    "Outcome predicted by Theorems 1 & 2 before running. Not executed because the result is mathematically predetermined.",
  "DEFERRED":
    "Not run. Upstream structural results make the outcome predetermined; the simulation would add no new information.",
};

const PAPER_TIPS = {
  fail:
    "FAIL — the simulation programme found a fundamental conflict with observational data that cannot be resolved within the CMSTG action class.",
  pass:
    "PASS — consistent with all observational constraints. A '2.77σ floor' qualifier means a structural minimum tension persists but is not a violation of the theory.",
  partial:
    "PARTIAL — satisfies some but not all constraints; further work is ongoing.",
};

const PHASE_TIPS = {
  phase1:
    "Phase 1 locked the canonical CMSTG action after 22 independent observational tests. This action is the foundation for all subsequent work — no parameters were changed after Phase 1.",
  phase2:
    "Phase 2 extended the matter sector by adding the χ dark matter field. The χ-DM oscillating field passes the SPARC 161-galaxy rotation curve fit with m₂₂ ≈ 0.31 (ultralight axion mass ≈ 0.31 × 10⁻²² eV).",
  phase3:
    "Phase 3 exhaustively tested curvature-sourced scalar mechanisms for reducing the DESI tension. All six failed, establishing Theorem 1: curvature-sourced scalars always grow in the wrong direction.",
  phase4:
    "Phase 4 ran Tier 1 diagnostics and Tier 2 mechanism tests covering the full late-time extension space. All eight Tier 2 mechanisms failed, establishing Theorem 2 and closing the late-time possibility space.",
};

const SIM_TIPS = {
  "SIM80":
    "Derives and validates the resummed curvature-memory propagator — the Green's function of the memory kernel K(x,x'), summed to all loop orders. Without resummation the bare propagator diverges at two loops.",
  "SIM82":
    "Verifies that K(x,x') is causal (zero outside the light cone) and has the spatial Yukawa profile K(r) ∝ e^{−m₀r}/r — the expected form of a massive scalar propagator, required for unitarity.",
  "SIM83–86":
    "Four-simulation sweep over memory cutoff scale and coupling parameter space. Confirms Friedmann equation stability — the scalar background does not destabilise FRW expansion under perturbations.",
  "SIM87":
    "Full-covariance BAO refit combining BOSS galaxy clustering, eBOSS quasar correlations, and Lyman-α forest data simultaneously, accounting for cross-dataset correlations. Establishes the 2.77σ DESI floor.",
  "SIM88":
    "Gravitational wave speed test (|cᴛ/c − 1| < 10⁻¹⁵ from GW170817 + GRB 170817A) and Solar System effective gravitational constant G_eff. Both pass well within observational bounds.",
  "SIM89":
    "Tests whether CMSTG shifts the CMB acoustic peak scale θ* = rₛ/Dₐ outside Planck's measurement 100θ* = 1.04101 ± 0.00029. Any shift here would immediately rule out the theory — this passes.",
  "SIM90":
    "Joint parameter fit to Planck CMB and BAO data simultaneously. Δχ² ≈ 0 vs ΛCDM: the two models are statistically indistinguishable. Λ₀ is largely unconstrained (flat landscape around Λ₀ ≈ 0.013).",
  "SIM91–96":
    "Six independent probes: Redshift-Space Distortions (growth rate f·σ₈), matter fluctuation amplitude σ₈, CMB lensing, Big Bang Nucleosynthesis light-element yields, neutrino mass bounds, Lyman-α forest power spectrum. All pass.",
  "SIM97":
    "Full Planck plikHM TTTEEE likelihood evaluated at parameters fitted from BAO only. The Δχ² = +7.0 penalty shows BAO-only parameters are suboptimal for the full CMB — motivates the joint refit in SIM98.",
  "SIM98":
    "Joint optimisation of Planck plikHM + BAO simultaneously. Δ(−2lnL) = −0.013 is negligible: a single parameter set satisfies both datasets with essentially no CMSTG–ΛCDM discrimination.",
  "SIM99–100":
    "Two independent tests of whether the CMSTG scalar Ψ can replace dark matter in galaxy rotation curves. Both fail: Ψ produces a profile too steep at sub-kpc scales. Phase 2 χ-DM is the surviving route.",
  "SIM101–106":
    "Six precision calculations: sound horizon rᴅ, two-loop UV divergence cancellation, Ward identity (gauge invariance), and renormalisation group flow of Λ₀ confirming the UV fixed point at Λ₀ → 0.",
  "SIM107–108":
    "Post-Newtonian Solar System tests (PPN parameters γ and β) and GW propagation-speed constraint from GW170817. Both pass well within current observational bounds.",
  "SIM109–110":
    "DESI dark energy equation-of-state w(z) fit and inflationary initial conditions. Neither reduces the 2.77σ DESI floor; the tension is confirmed as robust to w(z) parameterisation and inflationary priors.",
  "SIM111":
    "Full scan of the memory-kernel mass m₀ across its prior. Minimum achievable DESI tension is 3.44σ at m₀ ≈ 0.18 Mpc⁻¹ — establishes the structural floor and locks the Phase 1 action.",
  "SIM112–113":
    "λΨ⁴ self-interaction and hilltop quintessence potential added to the scalar sector. Both modify the dark energy equation of state but cannot push the BAO tension below the structural floor.",
  "SIM114–115":
    "χ dark matter as a BEC condensate (single ground state) and as a gradient soliton (spatially varying). Both fail rotation curve fitting before oscillating-field physics is included — the condensate profile is too extended.",
  "SIM116–119":
    "Four-simulation validation of the χ-DM oscillating field against the SPARC 161-galaxy dataset. Reproduces flat rotation curves across v_flat ~ 30–300 km/s with m₂₂ ≈ 0.31. Key Phase 2 success.",
  "SIM120–121C":
    "Joint DE+DM optimisation: both Ψ (dark energy) and χ (dark matter) sectors free simultaneously. Even with maximum freedom, converges on the 2.77σ DESI floor — confirms it is structural.",
  "SIM122–130":
    "Nine matter-sector preparatory probes for Phase 3. All fail to reduce the DESI tension, ruling out further matter-sector modifications and motivating the curvature-memory scalar tests in Phase 3.",
  "SIM131":
    "Additive conformal coupling ξΨR. The curvature source R > 0 drives Ψ to larger values, increasing F_eff and suppressing H(z) — the wrong direction. DESI tension: 5.97σ.",
  "SIM132":
    "Non-local Deser–Woodard analog: ψ sourced by □⁻¹R. Same directional failure despite different non-local structure. Tension: 3.59–3.83σ. Confirms the curvature-sourced pattern.",
  "SIM133":
    "Gauss–Bonnet sourcing: φ driven by G = R² − 4R_μν² + R_μνρσ². The topological combination amplifies the curvature signal, producing extreme 9.68σ DESI tension.",
  "SIM134":
    "Dilaton / exponential coupling e^{αΨ}R. Despite a different functional form, still curvature-sourced and fails in the same direction: tension 3.59–3.84σ depending on α.",
  "SIM135":
    "Bi-scalar with Ψ frozen at its attractor and φ separately sourced by R. Removing Ψ freedom does not help — φ's curvature sourcing still fails: tension 3.73–5.60σ.",
  "SIM136":
    "Horndeski kinetic coupling G^{μν}∂_μΨ∂_νΨ — couples scalar kinetic energy to the Einstein tensor. Tension 3.53–3.75σ. Closes the curvature-sourced mechanism space; Theorem 1 established.",
  "SIM137":
    "Tier 1 diagnostic: analyses the pattern of SPARC rotation-curve failures. Residuals show a systematic radial bias consistent with a structural Ψ-sector problem, not random scatter.",
  "SIM138":
    "Tier 1 diagnostic: DESI per-bin sensitivity analysis. The 2.77σ residual is spread across bins z = 0.3–1.5 with no dominant outlier, ruling out a localised dataset error.",
  "SIM139":
    "Tier 1 diagnostic: RSD growth-rate shape. The f(z) residual could in principle be fixed by adjusting F_eff, but Theorems 1 & 2 prevent this adjustment — classified SHAPE_FIXABLE, not PASS.",
  "SIM140":
    "Tier 2 P4-A: step-function potential V(Ψ) with a transition at z ≈ 0.5. Predicted to fail by Theorem 2 without running — any late-time H(z) boost shifts θ* outside the Planck bound.",
  "SIM141":
    "Tier 2 P4-B: running Λ₀(a) / Brans–Dicke analog, allowing the coupling to evolve with the scale factor. Fails catastrophically: shifts the CMB acoustic scale by 63σ.",
  "SIM142":
    "Tier 2 P4-C: Galileon cubic G₃(Ψ)□Ψ from the Horndeski class. Provides kinetic screening at small scales but cannot raise H(z) enough to close the DESI gap without violating GW speed bounds.",
  "SIM143":
    "Tier 2 P4-D: bi-scalar Ψ + φ quintessence with independent potentials. Despite the additional freedom, still converges on the 2.77σ DESI floor — both scalars are late-time only.",
  "SIM144":
    "Tier 2 completeness probe P4-E: mixed-source φ with both curvature (ξ_R) and matter (β_m) couplings, designed to span the full remaining mechanism space. Fails, closing Tier 2 and establishing Theorem 2.",
  "SIM145–146":
    "Tier 3 UV recheck and predictions update — deferred because no Tier 2 mechanism passed. Would have validated UV properties of a successful Tier 2 mechanism; redundant given Theorem 2.",
};

const THEOREM_TIPS = {
  "1":
    "Proved by SIM131–136. The constraint is directional: in an expanding FRW background R > 0, any curvature-sourced scalar necessarily grows monotonically, increasing F_eff and suppressing H(z) — the opposite of what is needed to reduce DESI tension.",
  "2":
    "Proved by SIM137–144. The angular scale θ* = rₛ/D_C* couples the sound horizon (set before recombination) to the angular diameter distance (set by late-time expansion). Raising H(z) to fix DESI contracts D_C* and shifts θ* outside Planck's bound.",
};

const PRED_TIPS = {
  "01":
    "The decisive test. DESI Y3 (full five-year dataset, expected ≈2027) will either confirm or relieve the tension. ≥3σ → Phase 5 pre-recombination new physics required; <2σ → Phase 1 validated without modification.",
  "02":
    "m₂₂ ≡ m_χ/(10⁻²² eV), the dimensionless ultralight axion mass. The predicted range 0.28–0.34 for gas-dominated dwarfs is testable by next-generation IFU spectroscopy (JWST NIRSpec, ELT-HARMONI) measuring stellar velocity dispersions.",
  "03":
    "G_eff/G = 1 − 15 ppm at z = 0. Current Solar System precision is ~10 ppm — already close. Third-generation GW standard-siren networks (Einstein Telescope, Cosmic Explorer) may reach sub-ppm sensitivity.",
};

// ── Tooltip component (portal-based, follows cursor) ──────────
function Tip({ text, children, block = false }) {
  const [pos, setPos] = React.useState(null);
  if (!text) return children;
  const Tag = block ? 'div' : 'span';
  return (
    <Tag
      className="has-tip"
      onMouseEnter={e => setPos({ x: e.clientX, y: e.clientY })}
      onMouseMove={e => setPos({ x: e.clientX, y: e.clientY })}
      onMouseLeave={() => setPos(null)}
    >
      {children}
      {pos && ReactDOM.createPortal(
        <div className="tooltip-bubble" style={{ left: pos.x + 14, top: pos.y + 14 }}>
          {text}
        </div>,
        document.body
      )}
    </Tag>
  );
}

// ── Modal overlay ─────────────────────────────────────────────
function Modal({ url, onClose }) {
  React.useEffect(() => {
    const handler = e => { if (e.key === 'Escape') onClose(); };
    document.addEventListener('keydown', handler);
    return () => document.removeEventListener('keydown', handler);
  }, [onClose]);

  return (
    <div className="modal-backdrop" onClick={onClose}>
      <div className="modal-frame" onClick={e => e.stopPropagation()}>
        <button className="modal-close" onClick={onClose} aria-label="Close">✕</button>
        <iframe src={url} className="modal-iframe" title="CMSTG Explainer" />
      </div>
    </div>
  );
}

// ── KaTeX helper ─────────────────────────────────────────────
function K({ s }) {
  const ref = React.useRef(null);
  React.useEffect(() => {
    if (ref.current && window.katex) {
      katex.render(s, ref.current, { throwOnError: false, displayMode: false });
    }
  }, [s]);
  return <span ref={ref} />;
}

// ── Verdict badge ─────────────────────────────────────────────
function VBadge({ kind, label }) {
  return <span className={`vbadge vbadge--${kind}`}>{label}</span>;
}

// ── Sim verdict display (strip underscores) ───────────────────
function formatVerdict(v) { return v.replace(/_/g, ' '); }

// ── Sim verdict CSS class ─────────────────────────────────────
function simClass(verdict = '') {
  const v = verdict.toUpperCase();
  if (v === 'PASS')              return 'sim--pass';
  if (v.startsWith('FAIL') || v === 'PREDICTED FAIL') return 'sim--fail';
  if (v.startsWith('PARTIAL'))   return 'sim--partial';
  if (v === 'STRUCTURAL' || v === 'STRUCTURAL_PATTERN' || v === 'DISTRIBUTED' || v === 'SHAPE_FIXABLE') return 'sim--structural';
  if (v === 'DEFERRED')          return 'sim--deferred';
  return 'sim--other';
}

// ── NavBar ────────────────────────────────────────────────────
function NavBar() {
  const sections = ['parameters', 'papers', 'theorems', 'phases', 'predictions'];
  return (
    <nav className="navbar">
      <span className="nav-brand">CMSTG</span>
      <div className="nav-links">
        {sections.map(id => (
          <a key={id} href={`#${id}`} className="nav-link">{id}</a>
        ))}
        <a
          href={D.meta.repo}
          target="_blank"
          rel="noreferrer"
          className="nav-link nav-link--repo"
        >
          GitHub
        </a>
      </div>
    </nav>
  );
}

// ── Hero ──────────────────────────────────────────────────────
function Hero() {
  const m = D.meta;
  return (
    <header className="hero">
      <div className="hero-eyebrow">
        <span className="mono">CMSTG · v{m.version} · {m.released}</span>
        <span className="mono">{m.license}</span>
      </div>
      <h1 className="hero-title">{m.title}</h1>
      <p className="hero-tagline">{m.tagline}</p>
      <div className="hero-meta">
        <span>{m.author}</span>
        <a href={m.repo} target="_blank" rel="noreferrer" className="hero-link">
          GitHub repository →
        </a>
      </div>
    </header>
  );
}

// ── Parameters ────────────────────────────────────────────────
function Parameters() {
  return (
    <section className="section" id="parameters">
      <h2 className="sec-title">Canonical Parameters</h2>
      <div className="table-scroll">
      <table className="param-table">
        <thead>
          <tr>
            <th>Symbol</th>
            <th>Value</th>
            <th>Unit</th>
            <th>Status</th>
            <th>Note</th>
          </tr>
        </thead>
        <tbody>
          {D.parameters.map(p => (
            <tr key={p.sym}>
              <td className="mono">
                <Tip text={PARAM_TIPS[p.sym]}><K s={p.sym} /></Tip>
              </td>
              <td className="mono">{p.val}</td>
              <td className="mono">{p.unit ? <K s={p.unit} /> : '—'}</td>
              <td>
                <Tip text={STATUS_TIPS[p.status]}>
                  <span className={`pstatus pstatus--${p.status.replace(/[\s/]+/g, '-')}`}>
                    {p.status}
                  </span>
                </Tip>
              </td>
              <td className="td-note">{p.note ? (p.note.includes('\\') ? <K s={p.note} /> : p.note) : ''}</td>
            </tr>
          ))}
        </tbody>
      </table>
      </div>
    </section>
  );
}

// ── Papers ────────────────────────────────────────────────────
function Papers() {
  const [activeExplainer, setActiveExplainer] = React.useState(null);

  return (
    <section className="section" id="papers">
      <h2 className="sec-title">Publication Record</h2>
      <div className="papers-grid">
        {D.papers.map(p => (
          <article
            key={p.id}
            className={`paper-card paper-card--${p.verdictKind}`}
            onClick={() => p.explainer && setActiveExplainer(p.explainer)}
            title={p.explainer ? 'Click to open explainer' : undefined}
          >
            <div className="paper-head">
              <a
                href={p.explainer || paperRepoUrl(p)}
                target="_blank"
                rel="noreferrer"
                className="paper-id paper-id--link"
                onClick={e => e.stopPropagation()}
              >
                Paper {p.id} ↗
              </a>
              <Tip text={PAPER_TIPS[p.verdictKind]}>
                <VBadge kind={p.verdictKind} label={p.verdict} />
              </Tip>
            </div>
            <h3 className="paper-title">{p.title}</h3>
            <p className="paper-summary">{p.summary}</p>
            <div className="paper-foot">
              <span className="paper-sims">{p.sims}</span>
              <a
                href={paperRepoUrl(p)}
                target="_blank"
                rel="noreferrer"
                className="paper-link"
                onClick={e => e.stopPropagation()}
              >
                Repo →
              </a>
            </div>
          </article>
        ))}
      </div>
      {activeExplainer && (
        <Modal url={activeExplainer} onClose={() => setActiveExplainer(null)} />
      )}
    </section>
  );
}

// ── Theorems ──────────────────────────────────────────────────
function Theorems() {
  return (
    <section className="section" id="theorems">
      <h2 className="sec-title">Structural No-Go Theorems</h2>
      <div className="theorems-stack">
        {D.theorems.map(t => (
          <div key={t.n} className="theorem-box">
            <span className="theorem-label">Theorem {t.n} · {t.phase}</span>
            <Tip text={THEOREM_TIPS[t.n]} block>
              <div className="theorem-title">{t.title}</div>
            </Tip>
            <div className="theorem-body">{t.body}</div>
          </div>
        ))}
      </div>
      <div className="loophole-box">
        <span className="loophole-label">Surviving route</span>
        <p className="loophole-text">{D.loophole}</p>
      </div>
    </section>
  );
}

// ── Simulation phases ─────────────────────────────────────────
function PhaseBlock({ phase }) {
  const { compactPhases, showAllSims } = React.useContext(Ctx);
  const [open, setOpen] = React.useState(true);

  React.useEffect(() => {
    setOpen(!compactPhases);
  }, [compactPhases]);

  const hasTiers = phase.sims.some(s => s.tier);

  return (
    <div className="phase-block">
      <button className="phase-header" onClick={() => setOpen(o => !o)}>
        <div className="phase-left">
          <Tip text={PHASE_TIPS[phase.key]}>
            <span className="phase-label">{phase.label} · {phase.range}</span>
          </Tip>
          <span className="phase-title">{phase.title}</span>
        </div>
        <span className="phase-chevron">{open ? '▾' : '▸'}</span>
      </button>
      <p className="phase-summary">{phase.summary}</p>
      {open && (
        <div className="phase-sims">
          <table className="sim-table">
            <thead>
              <tr>
                <th>ID</th>
                <th>Topic</th>
                {hasTiers && <th>Tier</th>}
                <th>Verdict</th>
                {showAllSims && <th>Detail</th>}
              </tr>
            </thead>
            <tbody>
              {phase.sims.map(sim => {
                const url = simRepoUrl(sim.id);
                return (
                  <tr
                    key={sim.id}
                    className={simClass(sim.verdict) + (url ? ' sim--linked' : '')}
                    onClick={url ? () => window.open(url, '_blank', 'noreferrer') : undefined}
                    title={url ? `Open ${sim.id} in repository` : undefined}
                  >
                    <td className="sim-id">{sim.id}</td>
                    <td className="sim-topic">
                      <Tip text={SIM_TIPS[sim.id]}>{sim.topic}</Tip>
                    </td>
                    {hasTiers && <td className="sim-tier">{sim.tier || ''}</td>}
                    <td className="sim-verdict">
                      <Tip text={VERDICT_TIPS[sim.verdict]}>{formatVerdict(sim.verdict)}</Tip>
                    </td>
                    {showAllSims && <td className="sim-detail">{sim.detail || ''}</td>}
                  </tr>
                );
              })}
            </tbody>
          </table>
        </div>
      )}
    </div>
  );
}

function Phases() {
  return (
    <section className="section" id="phases">
      <h2 className="sec-title">Simulation Programme</h2>
      <div className="phases-stack">
        {D.phases.map(p => <PhaseBlock key={p.key} phase={p} />)}
      </div>
    </section>
  );
}

// ── Predictions ───────────────────────────────────────────────
function Predictions() {
  return (
    <section className="section" id="predictions">
      <h2 className="sec-title">Falsifiable Predictions</h2>
      <div className="pred-grid">
        {D.predictions.map(p => (
          <div key={p.n} className="pred-card">
            <div className="pred-head">
              <span className="pred-n">P{p.n}</span>
              <span className="pred-milestone">{p.milestone}</span>
            </div>
            <Tip text={PRED_TIPS[p.n]} block>
              <div className="pred-title">{p.title}</div>
            </Tip>
            <p className="pred-body">{p.body}</p>
          </div>
        ))}
      </div>
    </section>
  );
}

// ── Footer ────────────────────────────────────────────────────
function Footer() {
  const m = D.meta;
  return (
    <footer className="footer">
      <span className="mono">{m.title} · v{m.version} · {m.license}</span>
      <a href={m.repo} target="_blank" rel="noreferrer" className="footer-link mono">
        cisomorph/Curvature-Memory_Scalar-Tensor_Gravity
      </a>
    </footer>
  );
}

// ── Root app ──────────────────────────────────────────────────
function App() {
  return (
    <window.TweaksProvider>
      <div className="app-wrap">
        <NavBar />
        <main className="main">
          <Hero />
          <Parameters />
          <Papers />
          <Theorems />
          <Phases />
          <Predictions />
        </main>
        <Footer />
        <window.TweaksPanel />
      </div>
    </window.TweaksProvider>
  );
}

const root = ReactDOM.createRoot(document.getElementById('root'));
root.render(<App />);
