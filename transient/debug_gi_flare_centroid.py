"""Diagnose what drives the flare-window g vs i ICCF centroid.

Prints, for the detrended flare-window g-i CCF: the peak, the extent of
the contiguous r > 0.8 r_max segment used for the centroid, and where the
FR/RSS Monte-Carlo centroids actually land (fraction near the peak vs at
the |tau| ~ 100 d edges).

Run: python3 transient/debug_gi_flare_centroid.py
"""

import numpy as np

import lag_iccf as L


def segment_bounds(tau_grid, r):
    ipk = np.nanargmax(r)
    above = np.isfinite(r) & (r > 0.8 * r[ipk])
    lo = ipk
    while lo > 0 and above[lo - 1]:
        lo -= 1
    hi = ipk
    while hi < len(r) - 1 and above[hi + 1]:
        hi += 1
    return tau_grid[lo], tau_grid[hi], tau_grid[ipk], r[ipk]


def main():
    t1, f1, e1 = L.load_lc('g', 58200.0, 59300.0)
    t2, f2, e2 = L.load_lc('i', 58200.0, 59300.0)
    f1d = L.detrend(t1, f1, L.DETREND_SIGMA)
    f2d = L.detrend(t2, f2, L.DETREND_SIGMA)
    print(f'points: g {len(t1)}, i {len(t2)} (flare window)')

    r = L.iccf(t1, f1d, t2, f2d)
    lo, hi, tpk, rmax = segment_bounds(L.TAU_GRID, r)
    _, cent, _ = L.peak_and_centroid(L.TAU_GRID, r)
    print(f'real CCF: rmax = {rmax:.3f} at tau = {tpk:+.2f} d')
    print(f'centroid segment (contiguous r > 0.8 rmax): '
          f'[{lo:+.2f}, {hi:+.2f}] d -> centroid {cent:+.2f} d')

    nmc = 400
    cents = []
    for _ in range(nmc):
        rr = L.iccf(*L.fr_rss(t1, f1d, e1)[:2], *L.fr_rss(t2, f2d, e2)[:2])
        _, ct, _ = L.peak_and_centroid(L.TAU_GRID, rr)
        if np.isfinite(ct):
            cents.append(ct)
    cents = np.array(cents)
    print(f'\nMC centroids ({len(cents)} realizations):')
    for name, m in [('|tau| < 25 d      ', np.abs(cents) < 25),
                    ('25 < |tau| < 60 d ', (np.abs(cents) >= 25)
                     & (np.abs(cents) < 60)),
                    ('|tau| > 60 d (edge)', np.abs(cents) >= 60)]:
        print(f'  {name}: {m.mean() * 100:5.1f} %')
    print(f'  median {np.median(cents):+.2f}, '
          f'16/84 pct {np.percentile(cents, 15.87):+.2f} / '
          f'{np.percentile(cents, 84.13):+.2f}')


if __name__ == '__main__':
    main()
