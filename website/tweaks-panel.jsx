// tweaks-panel.jsx — display-options panel; loaded before app.jsx
// Exports: window.TweaksContext, window.TweaksProvider, window.TweaksPanel

const TweaksContext = React.createContext({
  showAllSims: false,
  setShowAllSims: () => {},
  compactPhases: false,
  setCompactPhases: () => {},
});
window.TweaksContext = TweaksContext;

function TweaksProvider({ children }) {
  const [showAllSims,    setShowAllSims]    = React.useState(false);
  const [compactPhases,  setCompactPhases]  = React.useState(false);

  return (
    <TweaksContext.Provider value={{ showAllSims, setShowAllSims, compactPhases, setCompactPhases }}>
      {children}
    </TweaksContext.Provider>
  );
}
window.TweaksProvider = TweaksProvider;

function TweaksPanel() {
  const { showAllSims, setShowAllSims, compactPhases, setCompactPhases } = React.useContext(TweaksContext);
  const [open, setOpen] = React.useState(false);

  return (
    <div className="tweaks-panel">
      <button
        className="tweaks-toggle"
        onClick={() => setOpen(o => !o)}
        title="Display options"
        aria-label="Toggle display options"
      >
        ⚙
      </button>
      {open && (
        <div className="tweaks-menu" role="menu">
          <div className="tweaks-title">Display</div>
          <label className="tweaks-row">
            <input
              type="checkbox"
              checked={showAllSims}
              onChange={e => setShowAllSims(e.target.checked)}
            />
            Show sim details
          </label>
          <label className="tweaks-row">
            <input
              type="checkbox"
              checked={compactPhases}
              onChange={e => setCompactPhases(e.target.checked)}
            />
            Collapse phases
          </label>
        </div>
      )}
    </div>
  );
}
window.TweaksPanel = TweaksPanel;
