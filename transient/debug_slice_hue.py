#!/usr/bin/env python3
"""Why is the RGB slice gray? Print the global-rank hue triplets of the
rays at the star azimuth and the luminance distribution of the deposit."""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from transient.sim_spectrum import load_sim, to_physical, star_position, CONFIG
from transient.sim_cloudy_spectrum import load_grid, compute_rays
from transient.sim_cloudy_linemaps import LINES as LINE_WINS


def main():
    sim = load_sim(os.path.join(HERE, 'data', 'sim', 'disk.out1.00012.athdf'))
    phys = to_physical(sim, CONFIG)
    grid = load_grid()
    rays = compute_rays(sim, phys, CONFIG, grid)
    wave = grid[0]
    r, th, ph = sim['r'], sim['th'], sim['ph']
    nph, nth = len(ph), len(th)
    xs, ys = star_position(sim)
    ip = np.argmin(np.abs(ph - np.mod(np.arctan2(ys, xs), 2 * np.pi)))

    F_all, dA_all = rays['F_line'], rays['dA']
    Lg = []
    for key in ('halpha', 'mgii', 'civ'):
        _, w1, w2 = LINE_WINS[key]
        m = (wave >= w1) & (wave <= w2)
        Lg.append(dA_all * np.trapezoid(F_all[:, m] / wave[None, m],
                                        wave[m], axis=1))
    Lg = np.stack(Lg).reshape(3, nph, nth)
    litg = Lg.sum(axis=0) > 0
    hue_g = np.zeros_like(Lg)
    for k in range(3):
        vals = np.log10(Lg[k][litg] + 1e-30)
        order = np.argsort(np.argsort(vals))
        hue_g[k][litg] = (order + 1) / vals.size
    hue = hue_g[:, ip, :]
    up = th <= np.pi / 2
    matter = rays['matter'].reshape(nph, nth)[ip]
    print("theta_idx  matter  rank(Ha)  rank(MgII)  rank(CIV)  spread")
    for j in np.where(up)[0][::4]:
        h = hue[:, j]
        print(f"  {j:3d}      {str(matter[j])[0]}     "
              f"{h[0]:.2f}      {h[1]:.2f}       {h[2]:.2f}     "
              f"{h.max()-h.min():.2f}")
    sp = (hue[:, up].max(axis=0) - hue[:, up].min(axis=0))
    print(f"hue spread (max-min of ranks), upper rays: "
          f"median {np.median(sp):.2f}, 90th {np.percentile(sp, 90):.2f}")

    # luminance: how much of the deposit sits within DR of the max?
    w = rays['w'][ip].copy()
    ws = w.sum(axis=1)
    w[ws <= 0, 0] = 1.0
    frac = w / w.sum(axis=1)[:, None]
    dtot = (Lg[:, ip, :].sum(axis=0)[:, None] * frac)
    logd = np.log10(dtot[dtot > 0])
    print(f"deposit log range: max {logd.max():.1f}, "
          f"p50 {np.median(logd):.1f}, p10 {np.percentile(logd, 10):.1f} "
          f"(so {100*(logd > logd.max()-2.5).mean():.0f}% of lit cells "
          f"within 2.5 dex of max)")


if __name__ == '__main__':
    main()
