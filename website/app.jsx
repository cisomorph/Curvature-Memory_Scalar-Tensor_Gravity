// app.jsx — main CMSTG site; depends on data.js and tweaks-panel.jsx

const D   = window.CMSTG_DATA;
const Ctx = window.TweaksContext;

const REPO_BASE = D.meta.repo + "/tree/main/";

// Build a GitHub tree URL for a sim ID string (handles ranges like "SIM83–86")
function simRepoUrl(id) {
  const match = id.match(/SIM\d+[a-zA-Z]*/);
  return match ? REPO_BASE + "Ordered_Simulations/" + match[0] : null;
}

function paperRepoUrl(p) {
  return p.repoPath ? REPO_BASE + p.repoPath : null;
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
              <td className="mono"><K s={p.sym} /></td>
              <td className="mono">{p.val}</td>
              <td className="mono">{p.unit ? <K s={p.unit} /> : '—'}</td>
              <td>
                <span className={`pstatus pstatus--${p.status.replace(/[\s/]+/g, '-')}`}>
                  {p.status}
                </span>
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
              <VBadge kind={p.verdictKind} label={p.verdict} />
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
            <div className="theorem-title">{t.title}</div>
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
          <span className="phase-label">{phase.label} · {phase.range}</span>
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
                    <td className="sim-topic">{sim.topic}</td>
                    {hasTiers && <td className="sim-tier">{sim.tier || ''}</td>}
                    <td className="sim-verdict">{formatVerdict(sim.verdict)}</td>
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
            <div className="pred-title">{p.title}</div>
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
