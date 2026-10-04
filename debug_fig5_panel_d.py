"""Measure what panel (d) of Figure 5 (two-pillar H-alpha map, log-ratio of
with/without pillars) actually shows, to rewrite its description in Section 6.3.

Run:  python3 debug_fig5_panel_d.py
"""
import numpy as np

C_KMS = 2.998e5
d = np.load('plots/vdm_manual_pillars_maps.npz')
lam, tau, psi, psi_no = d['lam'], d['tau'], d['psi'], d['psi_no']
l0 = float(d['lambda0'])
v = C_KMS * (lam - l0) / l0
with np.errstate(divide='ignore', invalid='ignore'):
    r = np.where((psi > 0) & (psi_no > 0), np.log10(psi / psi_no), np.nan)

for label, sel in (('deficit below -0.3 dex', r < -0.3),
                   ('deficit below -0.1 dex', r < -0.1),
                   ('excess above +0.1 dex', r > 0.1)):
    it, il = np.nonzero(sel)
    if it.size == 0:
        print(label, ': none')
        continue
    print('%-24s %5d pixels; velocity %+6.0f to %+6.0f km/s; delay %4.1f to %4.1f d'
          % (label, it.size, v[il].min(), v[il].max(), tau[it].min(), tau[it].max()))

it, il = np.unravel_index(np.nanargmin(r), r.shape)
print('strongest deficit %.2f dex (%+.0f%%) at v = %+.0f km/s, tau = %.1f d'
      % (r[it, il], 100 * (10 ** r[it, il] - 1), v[il], tau[it]))
it, il = np.unravel_index(np.nanargmax(r), r.shape)
print('strongest excess  %.2f dex (%+.0f%%) at v = %+.0f km/s, tau = %.1f d'
      % (r[it, il], 100 * (10 ** r[it, il] - 1), v[il], tau[it]))

# Response-weighted change near each pillar's hotspot
for name, vc, tc in (('pillar phi = 3pi/4', 4700., 6.0), ('pillar phi = 0', 0., 1.2)):
    box = (np.abs(v[None, :] - vc) < 1000) & (np.abs(tau[:, None] - tc) < 1.5)
    print('%-20s summed response in a +-1000 km/s, +-1.5 d box: change %+.0f%%'
          % (name, 100 * (psi[box].sum() / psi_no[box].sum() - 1)))
print('total response change over the whole map: %+.1f%%'
      % (100 * (psi.sum() / psi_no.sum() - 1)))
