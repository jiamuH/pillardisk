"""diag_moc_diffuse.py - MOCASSIN (3D, diffuse field included) vs a direct-only straight-ray
march on its own grid: which of its ionised cells / line light lie beyond
the lamp's direct front? (docs/NOTE_diffuse_field.md)
Run: python3 -m pillardisk.transient.debug.diag_moc_diffuse <out.npz>  (~40 s)"""
import numpy as np, time, sys
t0 = time.time()
from pillardisk.transient.moc_maps import load_run
from pillardisk.transient.scripts.plot_cmi_compare import march_mask
run = '/data2/jhuang/runs/mocassin/moc_edd10'
d = load_run(run)
print(f'loaded {time.time() - t0:.0f} s', flush=True)
g = np.load(f'{run}/grid.npz')
x, y, z = (a / d['r0_cm'] for a in d['axes'])
nH = d['nH']
ion = march_mask(x, y, z, nH, float(g['q_ion']), d['r0_cm'],
                 n_phi=360, n_th=180, n_s=600)
print(f'march {time.time() - t0:.0f} s', flush=True)
act = d['active'] & (z[None, None, :] > 0)
xh = d['ions'][1][..., 1]
L = d.get('lines', {})
lines = {k: L[k] for k in ('Ha',)} if 'Ha' in L else {}
if 'MgII2804' in L:
    lines['MgII'] = L['MgII2804'] + L['MgII2796']
    lines['CIV'] = L['CIV1551'] + L['CIV1548']
np.savez(sys.argv[1], direct=ion, active=act, xh=xh, nH=nH, Te=d['Te'], x=x, y=y, z=z,
         **{f'L_{k}': v for k, v in lines.items()})
beyond = act & ~ion
lit = act & ion
print(f'upper-half active cells: {act.sum()}, inside direct front {lit.sum()}, beyond {beyond.sum()}')
for lab, m in (('inside', lit), ('beyond', beyond)):
    w = nH[m]
    print(f'{lab:7s} x(H+) median {np.median(xh[m]):.3f}, frac x>0.5: {np.mean(xh[m] > 0.5):.0%}, '
          f'nH-weighted mean x {np.sum(w * xh[m]) / np.sum(w):.3f}')
for k, v in lines.items():
    tot = v[act].sum()
    print(f'{k:5s} share of upper-half light beyond the direct front: {v[beyond].sum() / tot:.1%}')

# ---- where the light beyond the direct front sits (cells from the front)
from scipy.ndimage import distance_transform_edt
from pillardisk.transient.cloudy_rays import N_SURF
dist = distance_transform_edt(~ion)
print('cells beyond   n      x(H+) median  x>0.1   x>0.5   ' + '  '.join(lines))
for lo, hi, lab in ((0.5, 1.5, '1'), (1.5, 2.5, '2'), (2.5, 5.5, '3-5'),
                    (5.5, 1e9, '>5')):
    m = beyond & (dist >= lo) & (dist < hi)
    print(f'{lab:12s} {m.sum():6d}  {np.median(xh[m]):10.3g}  '
          f'{np.mean(xh[m] > 0.1):6.1%}  {np.mean(xh[m] > 0.5):6.1%}  ' +
          '  '.join(f'{v[m].sum() / v[act].sum():6.1%}'
                    for v in lines.values()))
up = z > 0
for lab, want in (('top dense cell reached directly', True),
                  ('top dense cell beyond the front', False)):
    xs = []
    for i in range(nH.shape[0]):
        for j in range(nH.shape[1]):
            k = np.nonzero(act[i, j] & (nH[i, j] > N_SURF) & up)[0]
            if len(k) and bool(ion[i, j, k[-1]]) == want:
                xs.append(xh[i, j, k[-1]])
    print(f'{lab}: n = {len(xs)}, x(H+) median {np.median(xs):.3g}, '
          f'16-84% {np.percentile(xs, 16):.2g}-{np.percentile(xs, 84):.2g}')
