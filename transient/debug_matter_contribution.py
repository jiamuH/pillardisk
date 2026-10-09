#!/usr/bin/env python3
"""How much of the reprocessed light comes from matter-bounded
(transparent) rays vs radiation-bounded (front-forming) rays, in the
current upper-hemisphere sum at i = 45 deg."""
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
    sim = load_sim(path)
    phys = to_physical(sim, CONFIG)
    grid = load_grid()
    rays = compute_rays(sim, phys, CONFIG, grid)
    wave = grid[0]
    matter = rays['matter'].reshape(-1)
    dA = rays['dA']                       # upper hemisphere only (dA=0 below)
    active = dA > 0
    print(f"active (upper) rays: {active.sum()}, of which "
          f"{(matter & active).sum()} matter-bounded "
          f"({100*(matter & active).sum()/active.sum():.0f}%)")

    for name, F in [('line', rays['F_line']), ('continuum', rays['F_cont'])]:
        Ltot = (dA[:, None] * F).sum(axis=0)
        Lmat = (dA[matter, None] * F[matter]).sum(axis=0)
        tot = np.trapezoid(Ltot / wave, wave)
        mat = np.trapezoid(Lmat / wave, wave)
        print(f"{name:9s}: total {tot:.2e} erg/s, "
              f"matter-bounded share {100*mat/tot:.0f}%")

    Flam = rays['F_line'] / wave[None, :]
    for key, (w1, w2) in LINES.items():
        m = (wave >= w1) & (wave <= w2)
        tot = np.trapezoid(dA[:, None] * Flam[:, m], wave[m], axis=1).sum()
        mat = np.trapezoid(dA[matter, None] * Flam[matter][:, m],
                           wave[m], axis=1).sum()
        print(f"{key:7s}: matter-bounded share {100*mat/tot:.0f}%")


if __name__ == '__main__':
    main()
