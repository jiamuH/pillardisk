#!/usr/bin/env python3
"""Diagnose the blocky rectangular tabs at the outer edge of the
bottom-left arm in the face-on maps (phi ~ 230-270 deg): is the sharp
outer boundary of the lit region jumping because of the DATA (front
radius / density) or because of the Cloudy lookup (clipping)?

Prints, per azimuth column: the outermost radius that receives Halpha
deposit, the max sub-cell front radius among upper rays, the
emission-weighted log n, and the clipped fraction."""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from transient.sim_spectrum import load_sim, to_physical, CONFIG
from transient.sim_cloudy_spectrum import load_grid, compute_rays


def main():
    sim = load_sim(os.path.join(HERE, 'data', 'sim', 'disk.out1.00012.athdf'))
    phys = to_physical(sim, CONFIG)
    grid = load_grid()
    rays = compute_rays(sim, phys, CONFIG, grid)
    r, th, ph = sim['r'], sim['th'], sim['ph']
    nph, nth = len(ph), len(th)
    dA = rays['dA'].reshape(nph, nth)
    up = dA > 0
    r_if = rays['r_if'].reshape(nph, nth)
    matter = rays['matter'].reshape(nph, nth)
    pv, hv, cv = grid[2]
    logn = rays['logn'].reshape(nph, nth)

    sel = (np.degrees(ph) >= 228) & (np.degrees(ph) <= 272)
    print("deg    max_rIF(rad-bnd,upper)  n_matter  n_clipfloor")
    for i in np.where(sel)[0][::2]:
        rb = up[i] & ~matter[i]
        rmax = r_if[i][rb].max() if rb.any() else np.nan
        ncl = int(((logn[i] <= hv[0] + 1e-6) & up[i]).sum())
        print(f"{np.degrees(ph[i]):6.1f}     {rmax:5.3f}              "
              f"{int((matter[i] & up[i]).sum()):2d}        {ncl:2d}")


if __name__ == '__main__':
    main()
