"""Measure what Figure 6 (100-pillar three-line velocity-delay maps) shows, to
check its description in Section 6.2. Reads the map cache written by
make_vdm_frame; computes nothing new.

Run:  python3 debug_fig6_text.py
"""
import numpy as np

C_KMS = 2.998e5
LAMBDA0 = {'C4': 1549.0, 'Mg2': 2798.0, 'Halpha': 6562.8, 'HI': 4861.32}

c = np.load('plots/vdm_frame_maps.npz')
tau = c['tau']
dlog, v = {}, {}
for k in LAMBDA0:
    w, n = c['with_' + k], c['no_' + k]
    v[k] = C_KMS * (c['lam_' + k] - LAMBDA0[k]) / LAMBDA0[k]
    with np.errstate(divide='ignore', invalid='ignore'):
        dlog[k] = np.where((w > 0) & (n > 0), np.log10(w / n), np.nan)
    # Shape of the with-pillar response (top row)
    tot = w.sum()
    tmean = (w.sum(axis=1) * tau).sum() / tot
    vrms = np.sqrt((w.sum(axis=0) * v[k] ** 2).sum() / tot)
    early = w[tau < 5].sum() / tot
    print('%-6s top row: mean delay %5.1f d, rms velocity %5.0f km/s, '
          'fraction at delay < 5 d %.2f' % (k, tmean, vrms, early))

# Only pixels carrying a meaningful response (1e-3 of each map's peak)
mask = {k: c['no_' + k] > 1e-3 * c['no_' + k].max() for k in LAMBDA0}
print('\nMiddle row, Delta log Psi (own response), over responding pixels:')
for k in ('C4', 'Mg2', 'Halpha', 'HI'):
    d = dlog[k][mask[k]]
    d = d[np.isfinite(d)]
    print('%-6s min %+5.2f  1%% %+5.2f  99%% %+5.2f  max %+5.2f   '
          'fraction below -0.3 dex %.2f, above +0.3 dex %.2f'
          % (k, d.min(), *np.percentile(d, [1, 99]), d.max(),
             (d < -0.3).mean(), (d > 0.3).mean()))

print('\nPixel-by-pixel correlation of Delta log Psi between lines:')
for a, b in (('C4', 'Mg2'), ('C4', 'Halpha'), ('Mg2', 'Halpha'), ('Halpha', 'HI')):
    m = mask[a] & mask[b] & np.isfinite(dlog[a]) & np.isfinite(dlog[b])
    print('  %-6s vs %-6s r = %.2f' % (a, b, np.corrcoef(dlog[a][m], dlog[b][m])[0, 1]))

# Bottom row: Delta log(Psi/Psi_Hbeta), split by where H-beta itself drops
# (shadowed cells) or rises (lit pillar faces). Same pixel grid in velocity
# is assumed: the maps share the velocity grid if nlambda and the velocity
# range are the same, which is checked here.
print('\nBottom row, Delta log(Psi/Psi_Hbeta):')
for k in ('C4', 'Mg2', 'Halpha'):
    assert np.allclose(v[k], v['HI'], atol=50), 'velocity grids differ'
    r = dlog[k] - dlog['HI']
    m = mask[k] & mask['HI'] & np.isfinite(r)
    shadow = m & (dlog['HI'] < -0.3)
    lit = m & (dlog['HI'] > 0.3)
    print('%-6s all: 1%% %+5.2f 99%% %+5.2f | H-beta shadowed (%4d px): median %+5.2f'
          ' | H-beta lit (%4d px): median %+5.2f'
          % (k, *np.percentile(r[m], [1, 99]), shadow.sum(), np.median(r[shadow]),
             lit.sum(), np.median(r[lit])))

print('\nSpread of the bottom row within H-beta-shadowed and H-beta-lit cells:')
for k in ('C4', 'Mg2', 'Halpha'):
    r = dlog[k] - dlog['HI']
    m = mask[k] & mask['HI'] & np.isfinite(r)
    for name, sel in (('shadowed', m & (dlog['HI'] < -0.3)), ('lit', m & (dlog['HI'] > 0.3))):
        x = r[sel]
        print('%-6s %-8s 10%% %+5.2f  90%% %+5.2f  fraction > +0.1 %.2f, < -0.1 %.2f'
              % (k, name, *np.percentile(x, [10, 90]), (x > 0.1).mean(), (x < -0.1).mean()))

print('\nDelay range of cells changed by more than 0.3 dex (middle row):')
for k in ('C4', 'Mg2', 'Halpha'):
    it, il = np.nonzero(mask[k] & (np.abs(np.nan_to_num(dlog[k])) > 0.3))
    print('%-6s delay 5-95%%: %4.1f to %4.1f d; |v| 95%%: %5.0f km/s'
          % (k, *np.percentile(tau[it], [5, 95]), np.percentile(np.abs(v[k][il]), 95)))
