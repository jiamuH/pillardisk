#!/usr/bin/env python3
"""Diagnose the RGB line-map channels: after the per-channel log stretch,
which line 'wins' where, and how far apart the channels sit. Decides how
to fix the muddy composite."""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from transient.sim_spectrum import load_sim, to_physical, CONFIG
from transient.sim_cloudy_spectrum import load_grid, compute_rays

LINES = {'civ': (1515., 1585.), 'mgii': (2755., 2845.),
         'halpha': (6510., 6620.)}


def main():
    path = os.path.join(HERE, 'data', 'sim', 'disk.out1.00012.athdf')
    cfg = CONFIG
    sim = load_sim(path)
    phys = to_physical(sim, cfg)
    grid = load_grid()
    rays = compute_rays(sim, phys, cfg, grid)
    wave = grid[0]
    r, th, ph = sim['r'], sim['th'], sim['ph']
    F_line = rays['F_line'].reshape(len(ph), len(th), -1)
    dA = rays['dA'].reshape(len(ph), len(th))
    w = rays['w'].copy()
    wsum = w.sum(axis=2)
    zw = wsum <= 0
    w[zw, :] = 0.0
    w[zw, 0] = 1.0
    frac = w / (w.sum(axis=2)[..., None])

    maps = {}
    for key, (w1, w2) in LINES.items():
        m = (wave >= w1) & (wave <= w2)
        Flam = F_line[:, :, m] / wave[None, None, m]
        L_ray = dA * np.trapz(Flam, wave[m], axis=2)
        maps[key] = (L_ray[..., None] * frac).sum(axis=1)

    DR = 2.5
    v = {}
    for key in LINES:
        lm = maps[key]
        v[key] = np.clip((np.log10(lm + 1e-30)
                          - (np.log10(lm.max()) - DR)) / DR, 0, 1)
        print(f"{key:7s}: map max {lm.max():.2e}, "
              f"stretched percentiles 50/90/99 = "
              f"{np.percentile(v[key], 50):.2f} "
              f"{np.percentile(v[key], 90):.2f} "
              f"{np.percentile(v[key], 99):.2f}")

    stack = np.stack([v['halpha'], v['mgii'], v['civ']])
    lum = stack.mean(axis=0)
    lit = lum > 0.05
    winner = stack.argmax(axis=0)
    names = ['halpha', 'mgii', 'civ']
    print(f"\nOf {lit.sum()} lit pixels (lum > 0.05):")
    for k, name in enumerate(names):
        fracwin = (winner[lit] == k).mean()
        print(f"  {name:7s} wins {100*fracwin:5.1f}%")
    # how big is the winning margin typically?
    srt = np.sort(stack, axis=0)
    margin = srt[2] - srt[1]
    print(f"win margin over 2nd channel, lit pixels: "
          f"median {np.median(margin[lit]):.3f}, "
          f"90th {np.percentile(margin[lit], 90):.3f}")
    # correlation between stretched channels on lit pixels
    for a, b in [('halpha', 'mgii'), ('halpha', 'civ'), ('mgii', 'civ')]:
        c = np.corrcoef(v[a][lit], v[b][lit])[0, 1]
        print(f"corr({a}, {b}) on lit pixels = {c:.3f}")


if __name__ == '__main__':
    main()
