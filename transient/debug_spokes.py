#!/usr/bin/env python3
"""Diagnose the radial 'spoke' structure in the face-on line maps: do
the dark radial lanes coincide with azimuths whose ionization fronts sit
at small radius across all theta (shadow cones of dense arm segments),
or with azimuths whose rays clip at the Cloudy grid edge?"""
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
    wave = grid[0]
    r, th, ph = sim['r'], sim['th'], sim['ph']
    nph, nth = len(ph), len(th)

    dA = rays['dA'].reshape(nph, nth)
    up = dA > 0
    r_if = rays['r_if'].reshape(nph, nth)
    matter = rays['matter'].reshape(nph, nth)
    F_line = rays['F_line'].reshape(nph, nth, -1)
    m = (wave >= 6350) & (wave <= 6800)
    L_ray = dA * np.trapezoid(F_line[:, :, m] / wave[None, None, m],
                              wave[m], axis=2)

    # per-azimuth: how far out does the lit region reach, and how bright?
    Lphi = L_ray.sum(axis=1)
    rif_max = np.where(up, r_if, 0).max(axis=1)     # outermost lit radius
    fmat = (matter & up).sum(axis=1) / up.sum(axis=1)
    # clipping fraction per azimuth
    pv, hv, cv = grid[2]
    logn = rays['logn'].reshape(nph, nth)
    clip = ((logn <= hv[0] + 1e-6) | (logn >= hv[-1] - 1e-6)) & up
    fclip = clip.sum(axis=1) / up.sum(axis=1)

    # the darkest azimuths (candidate spokes) in per-azimuth luminosity
    order = np.argsort(Lphi)
    print("deg     L_phi/med  outermost_lit_r  matter_frac  clip_frac")
    med = np.median(Lphi)
    for i in order[:8]:
        print(f"{np.degrees(ph[i]):6.1f}   {Lphi[i]/med:6.2f}     "
              f"{rif_max[i]:5.2f}          {fmat[i]:.2f}        {fclip[i]:.2f}")
    print("...brightest for contrast:")
    for i in order[-3:]:
        print(f"{np.degrees(ph[i]):6.1f}   {Lphi[i]/med:6.2f}     "
              f"{rif_max[i]:5.2f}          {fmat[i]:.2f}        {fclip[i]:.2f}")
    c1 = np.corrcoef(Lphi, rif_max)[0, 1]
    c2 = np.corrcoef(Lphi, fclip)[0, 1]
    c3 = np.corrcoef(Lphi, fmat)[0, 1]
    print(f"\ncorr(L_phi, outermost lit radius) = {c1:.2f}")
    print(f"corr(L_phi, clip fraction)        = {c2:.2f}")
    print(f"corr(L_phi, matter fraction)      = {c3:.2f}")


if __name__ == '__main__':
    main()
