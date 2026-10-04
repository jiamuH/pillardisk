"""Measure the residual amplitudes of Figures 7 (one pillar) and 8 (100
pillars) from the flux-map caches, and check that Figure 8 uses the same
pillars as Figure 6. Reads caches only; under a second.

Run:  python3 debug_fig78_amplitudes.py
"""
import numpy as np

for suffix in ('1pillar', 'N100'):
    c = np.load(f'plots/time_evolving_maps_{suffix}.npz')
    print(f'--- {suffix} ---')
    for k in ('C4', 'Mg2', 'Halpha'):
        f = c['with_' + k]
        res = (f - f.mean(axis=0)[None, :]) / f.mean()
        # Line-integrated flux: fractional change over time
        L = f.sum(axis=1)
        print('%-6s residual / <F>: min %+5.2f  max %+5.2f  1%% %+5.2f  99%% %+5.2f | '
              'line-integrated flux varies by +-%.1f per cent (half range)'
              % (k, res.min(), res.max(), *np.percentile(res, [1, 99]),
                 50 * (L.max() - L.min()) / L.mean()))

a = np.load('plots/time_evolving_maps_N100.npz')['pillar_rphi']
b = np.load('plots/vdm_frame_maps.npz')['pillar_rphi']
print('\nFigure 8 pillars identical to Figure 6:', a.shape == b.shape and np.allclose(a, b))
