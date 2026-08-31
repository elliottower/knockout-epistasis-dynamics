"""Generate Boolean vs ODE scatter plot for all 28 networks."""
import json
import glob
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

results = {}
for f in glob.glob('results/ode_full/*_ode.json'):
    with open(f) as fh:
        d = json.load(fh)
    results[d['model']] = d
for f in glob.glob('results/grn_v2/ode_full/*_ode.json'):
    if 'summary' in f:
        continue
    with open(f) as fh:
        d = json.load(fh)
    results[d['model']] = d

with open('results/grn_v2/merged_all_27_analysis.json') as f:
    merged = json.load(f)
bool_deltas = {m['model']: m.get('delta_o3plus') for m in merged['all_models']}
bool_deltas['grieco_bladder'] = 0.038

xs, ys, names = [], [], []
for name, r in sorted(results.items()):
    od = r.get('ode_delta_o3plus', 0) * 100
    bd = bool_deltas.get(name)
    if bd is None:
        continue
    xs.append(bd * 100)
    ys.append(od)
    names.append(name)

xs = np.array(xs)
ys = np.array(ys)

creation = (xs > 0.5) & (ys > 0.5)
destruction = (xs < -0.5) & (ys < -0.5)
null_agree = (~creation) & (~destruction) & (np.abs(xs) < 0.5) & (np.abs(ys) < 0.5)

fig, ax = plt.subplots(figsize=(4.5, 4.5))

c_create = '#2166ac'
c_destroy = '#b2182b'
c_null = '#878787'

ax.scatter(xs[creation], ys[creation], s=45, c=c_create,
           edgecolors='white', linewidths=0.4, zorder=3, alpha=0.85,
           label='Creation')
ax.scatter(xs[destruction], ys[destruction], s=45, c=c_destroy,
           edgecolors='white', linewidths=0.4, zorder=3, alpha=0.85,
           label='Destruction')
remaining = ~creation & ~destruction
ax.scatter(xs[remaining], ys[remaining], s=45, c=c_null,
           edgecolors='white', linewidths=0.4, zorder=3, alpha=0.85,
           label='Near-null')

lo = min(xs.min(), ys.min()) - 5
hi = max(xs.max(), ys.max()) + 5
ax.plot([lo, hi], [lo, hi], '--', color='#999999', linewidth=0.8, zorder=1)

r_val = np.corrcoef(xs, ys)[0, 1]
ax.text(0.04, 0.96, f'$r = {r_val:.3f}$\n$n = {len(xs)}$',
        transform=ax.transAxes, fontsize=10, va='top', ha='left',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                  edgecolor='#cccccc', alpha=0.9))

ax.set_xlabel(r'Boolean $\Delta_{3+}$ (pp)', fontsize=11)
ax.set_ylabel(r'Hill-function ODE $\Delta_{3+}$ (pp)', fontsize=11)
ax.set_xlim(lo, hi)
ax.set_ylim(lo, hi)
ax.set_aspect('equal')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend(fontsize=8, loc='lower right', framealpha=0.9,
          edgecolor='#cccccc')

plt.tight_layout()
out = Path('paper/figures/fig_ode_scatter.pdf')
plt.savefig(out, dpi=300, bbox_inches='tight')
plt.savefig(out.with_suffix('.png'), dpi=300, bbox_inches='tight')
print(f'Saved {out} and {out.with_suffix(".png")}')
print(f'  {len(xs)} networks plotted')
print(f'  Pearson r = {r_val:.3f}')
sign_ok = sum(1 for x, y in zip(xs, ys)
              if (x > 0 and y > 0) or (x < 0 and y < 0)
              or (abs(x) < 0.5 and abs(y) < 0.5))
print(f'  Sign preserved: {sign_ok}/{len(xs)}')
