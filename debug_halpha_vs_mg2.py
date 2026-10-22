"""Does H-alpha/H-beta follow Mg II/H-beta in Figure 6 (bottom row)? Compares
the two ratio shifts cell by cell, separately on lit faces and in shadows.
Reads a Figure 6 map cache; under a second.

Run:  python3 debug_halpha_vs_mg2.py [cache.npz]
"""
import sys

import numpy as np

c = np.load(sys.argv[1] if len(sys.argv) > 1 else 'plots/vdm_frame_maps.npz')
d = {}
for k in ('C4', 'Mg2', 'Halpha', 'HI'):
    w, n = c['with_' + k], c['no_' + k]
    with np.errstate(divide='ignore', invalid='ignore'):
        d[k] = np.where((w > 0) & (n > 0), np.log10(w / n), np.nan)
mask = np.ones_like(d['HI'], bool)
for k in d:
    mask &= c['no_' + k] > 1e-3 * c['no_' + k].max()
ha, mg = d['Halpha'] - d['HI'], d['Mg2'] - d['HI']
for name, sel in (('lit (H-beta up > 0.3 dex)', d['HI'] > 0.3),
                  ('shadow (H-beta down > 0.3 dex)', d['HI'] < -0.3),
                  ('all responding cells', np.ones_like(mask))):
    m = mask & sel & np.isfinite(ha) & np.isfinite(mg)
    print('%-32s %5d cells: median Mg II/Hb %+5.2f, H-alpha/Hb %+5.2f; '
          'correlation %.2f; same sign in %2.0f%% of cells'
          % (name, m.sum(), np.median(mg[m]), np.median(ha[m]),
             np.corrcoef(mg[m], ha[m])[0, 1], 100 * np.mean(np.sign(mg[m]) == np.sign(ha[m]))))
