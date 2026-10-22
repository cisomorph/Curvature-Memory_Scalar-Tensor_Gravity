#!/usr/bin/env python3
"""Summary figure of the audit: best θ*-matched Δχ² vs ΛCDM per model (DESI DR1 and DR2),
coloured by local-gravity viability. Numbers are taken from the audited output files."""
import json, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R = '/home/aion/CMSTG-fixes/Ordered_Simulations'
def best_matched(path, key='scan'):
    d = json.load(open(path))
    rows = d.get(key) or d.get('scan_results') or []
    best = None
    for r in rows:
        o = r.get('theta_matched')
        if o and (best is None or o['DR1']['delta_chi2_tot_vs_LCDM'] > best['DR1']['delta_chi2_tot_vs_LCDM']):
            best = o
    return best

entries = []   # (label, dDR1, dDR2, status)
def add(label, path, status, key='scan'):
    b = best_matched(os.path.join(R, path), key)
    if b: entries.append((label, b['DR1']['delta_chi2_tot_vs_LCDM'], b['DR2']['delta_chi2_tot_vs_LCDM'], status))

c = json.load(open(f'{R}/SIM121C/Outputs/sim121c_results.json'))['fits_same_pipeline']
c2 = json.load(open(f'{R}/SIM121C/Outputs/DR2/sim121c_results.json'))['fits_same_pipeline']
entries.append(('w0waCDM (SIM121C)', c['w0waCDM']['delta_chi2_vs_LCDM'], c2['w0waCDM']['delta_chi2_vs_LCDM'], 'reference'))
entries.append(('Phase 1 frozen F0 (SIM121C)', 0.0, 0.0, 'Cassini fail'))
add('SIM130 Phase 1 unlocked', 'SIM130/Outputs/sim130_results.json', 'Cassini fail')
add('SIM123 hilltop DE', 'SIM123/Outputs/sim123_results.json', 'Cassini fail')
add('SIM132 subtractive F', 'SIM132/Outputs/sim132_results.json', 'passes Cassini+LLR')
add('SIM134 exponential F', 'SIM134/Outputs/sim134_results.json', 'passes Cassini+LLR')
add('SIM135 bi-scalar ξφR', 'SIM135/Outputs/sim135_results.json', 'Cassini fail')
add('SIM136 Horndeski', 'SIM136/Outputs/sim136_results.json', 'c_T fail (GW170817)')
add('SIM141 running Λ0(a)', 'SIM141/output.json', 'LLR / Cassini fail', key='scan_results')
add('SIM142 G3 kinetic', 'SIM142/output.json', 'Cassini fail', key='scan_results')
add('SIM143 + quintessence', 'SIM143/output.json', 'Cassini fail', key='scan_results')
add('SIM144 mixed source', 'SIM144/SIM144_results.json', 'Cassini fail', key='scan_results')
r124 = json.load(open(f'{R}/SIM124/Outputs/sim124_results.json')).get('results_C', [])
if r124:
    b = max(r124, key=lambda r: r['obs']['DR1']['delta_chi2_tot_vs_LCDM'])['obs']
    entries.append(('SIM124 K(Ψ) hilltop', b['DR1']['delta_chi2_tot_vs_LCDM'], b['DR2']['delta_chi2_tot_vs_LCDM'], 'Cassini fail'))
add('SIM131 ξΨR', 'SIM131/Outputs/sim131_results.json', 'cosmology fail', key='scan_xi')

col = {'reference': '#7f7f7f', 'Cassini fail': '#d6604d', 'passes Cassini+LLR': '#1a9850',
       'c_T fail (GW170817)': '#b2abd2', 'LLR / Cassini fail': '#f4a582', 'cosmology fail': '#404040'}
entries.sort(key=lambda e: e[1])
fig, ax = plt.subplots(figsize=(10, 6.5))
y = np.arange(len(entries))
ax.barh(y + 0.2, [e[1] for e in entries], 0.4, color=[col[e[3]] for e in entries], label='DESI DR1')
ax.barh(y - 0.2, [e[2] for e in entries], 0.4, color=[col[e[3]] for e in entries], alpha=0.45, hatch='//', label='DESI DR2')
ax.set_yticks(y); ax.set_yticklabels([e[0] for e in entries], fontsize=9)
ax.axvline(0, color='k', lw=1)
ax.set_xlabel('best θ*-matched Δχ² vs ΛCDM  (> 0 = better than ΛCDM)')
ax.set_xlim(left=max(-60, min(e[1] for e in entries) - 5))
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=v, label=k) for k, v in col.items()] +
          [Patch(facecolor='w', edgecolor='k', label='solid: DR1'), Patch(facecolor='w', edgecolor='k', hatch='//', label='hatched: DR2')],
          fontsize=8, loc='lower right')
ax.set_title('CMSTG audit (2026-10-22): DESI + θ* fit vs ΛCDM, coloured by local-gravity status')
fig.tight_layout()
out = '/home/aion/CMSTG-fixes/audit_logs/audit_summary.png'
fig.savefig(out, dpi=150)
print('saved', out)
for e in entries: print(f"  {e[0]:30s} DR1 {e[1]:+8.2f}  DR2 {e[2]:+8.2f}  {e[3]}")
