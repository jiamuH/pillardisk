#!/usr/bin/env python3
"""Diagnose the terraced 'layer' structure in the RGB line map.

Hypothesis: the face-on map at fixed azimuth is a sum over the ~64 polar
rays, each depositing its luminosity along its own ionized radial segment
(ending at that ray's ionization front). Every time r crosses one ray's
front radius, a ray drops out of the sum -> a discrete brightness step.
Check: do the radial step locations in the H alpha map line up with the
distribution of per-ray ionization-front radii r_IF(theta)?
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from transient.sim_spectrum import load_sim, to_physical, CONFIG
from transient.sim_cloudy_spectrum import load_grid, compute_rays


def main():
    path = os.path.join(HERE, 'data', 'sim', 'disk.out1.00012.athdf')
    cfg = CONFIG
    sim = load_sim(path)
    phys = to_physical(sim, cfg)
    grid = load_grid()
    rays = compute_rays(sim, phys, cfg, grid)
    wave = grid[0]
    r, th, ph = sim['r'], sim['th'], sim['ph']
    nph, nth = len(ph), len(th)
    F_line = rays['F_line'].reshape(nph, nth, -1)
    dA = rays['dA'].reshape(nph, nth)
    r_if = rays['r_if'].reshape(nph, nth)

    w = rays['w'].copy()
    wsum = w.sum(axis=2)
    zw = wsum <= 0
    w[zw, :] = 0.0
    w[zw, 0] = 1.0
    frac = w / (w.sum(axis=2)[..., None])

    m = (wave >= 6510) & (wave <= 6620)
    Flam = F_line[:, :, m] / wave[None, None, m]
    L_ray = dA * np.trapz(Flam, wave[m], axis=2)
    lmap = (L_ray[..., None] * frac).sum(axis=1)     # (ph, r)

    # pick an azimuth on the visibly layered left side (phi ~ 180 deg)
    ipick = np.argmin(np.abs(np.degrees(ph) - 180.0))
    prof = lmap[ipick]
    # radial locations of the biggest downward log steps
    lp = np.log10(prof + 1e-30)
    dstep = np.diff(lp)
    big = np.argsort(dstep)[:12]                     # 12 largest drops
    r_steps = np.sort(0.5 * (r[big] + r[big + 1]))

    # per-theta front radii at that azimuth (radiation-bounded rays only)
    matter = rays['matter'].reshape(nph, nth)
    rif_here = np.sort(r_if[ipick][~matter[ipick]])
    print(f"azimuth {np.degrees(ph[ipick]):.1f} deg")
    print(f"largest radial brightness drops at r = "
          f"{np.array2string(r_steps, precision=3)}")
    print(f"per-theta ionization-front radii (radiation-bounded rays):\n"
          f"{np.array2string(rif_here, precision=3)}")
    # match: fraction of big steps within one radial cell of some r_IF
    dr = np.median(np.diff(r))
    match = [np.min(np.abs(rif_here - rs)) < 1.5 * dr for rs in r_steps]
    print(f"steps within 1.5 radial cells of a front radius: "
          f"{sum(match)}/{len(match)}")

    # how many rays still deposit at each radius (the 'layer count')
    nofront = (~matter[ipick]).sum()
    print(f"\nrays at this azimuth: {nth} total, {nofront} radiation-"
          f"bounded (finite front), {matter[ipick].sum()} matter-bounded")
    print(f"distinct front radii: {np.unique(rif_here).size}")

    # also: how much do the Cloudy lookup params cluster (flat patches)?
    logn = rays['logn'].reshape(nph, nth)
    logN = rays['logN'].reshape(nph, nth)
    logphi = rays['logphi'].reshape(nph, nth)
    for name, arr, lo, hi in [('logn', logn, 9.0, 12.0),
                              ('logN', logN, 21.5, 23.5),
                              ('logphi', logphi, 17.0, 21.0)]:
        clip = ((arr <= lo + 1e-6) | (arr >= hi - 1e-6)).mean()
        print(f"{name:7s}: {100*clip:4.1f}% of rays clipped at grid edge")


if __name__ == '__main__':
    main()
