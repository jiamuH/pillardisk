"""How heavily is the pillar belt shadowed, and what is the real ambient
temperature there?

The whole-grid shadow covering fraction is diluted by the inner disc, which has
no pillars. This prints the radial profile of the covering fraction and the
azimuthal temperature distribution at the belt, for the wide and narrow pillar
geometries and for the bowl, so the ambient temperature that the embedded-star
heating has to compete against is measured rather than assumed.
"""
import numpy as np

import test_lag_spectrum as tls
from fluxflux_lib import DiscGrid, eta_lp

HLAMP = 50 * 0.00399
tls.config['disk']['nphi'] = 720


def build(sigma_phi):
    dk = tls.build_disk(hlamp=HLAMP)
    if sigma_phi > 0:
        tls.add_pillars(dk, 100, 0.5, 0.25, sigma_phi, r_min=5.0)
    dk.d = 287.0 * 1e6 * 3.086e18 / (2.998e10 * 86400.0)
    g = DiscGrid(dk)
    dk.tx_base *= (50.0 / eta_lp(g)) ** 0.25
    dk.tv_base *= 1.5
    g.tv2 *= 1.5
    return dk, g


for sp in [0.30, 0.05, 0.0]:
    dk, g = build(sp)
    dk.tx_base *= 1.35                       # bright state
    T = dk.get_temperature(g.r2, g.p2, compute_shadows=True)
    if len(dk.pillars) > 0:
        smask = dk._compute_shadow_mask(g.r2, g.p2, g.h2)
    else:
        smask = np.ones_like(g.r2)
    label = f"sigma_phi = {sp}" if sp > 0 else "bowl (no pillars)"
    print(f"\n=== {label} ===")
    print(f"whole-grid shadow fraction: {1.0 - smask.mean():.3f}")
    print(f"{'r [ld]':>8}{'f_shadow':>10}{'T mean':>9}{'T 10%':>9}"
          f"{'T 50%':>9}{'T 90%':>9}{'T_visc':>9}")
    for rr in [3.0, 7.0, 10.0, 15.0, 19.0]:
        i = int(np.argmin(np.abs(dk.r - rr)))
        Ti = T[i]
        print(f"{dk.r[i]:>8.2f}{1.0 - smask[i].mean():>10.3f}{Ti.mean():>9.0f}"
              f"{np.percentile(Ti, 10):>9.0f}{np.percentile(Ti, 50):>9.0f}"
              f"{np.percentile(Ti, 90):>9.0f}"
              f"{np.interp(dk.r[i], dk.r, dk.tv_base):>9.0f}")
