#!/usr/bin/env python3
"""Are the nested yellow rings in the emissivity-placed maps the
individual theta-ray fronts (discretization artifact) or the arm walls
(physical)? Reproduce the map's radial profile at one azimuth with the
CURRENT smooth-ingredient interpolation, and compare its peak radii
with the per-theta front radii."""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from transient.sim_spectrum import load_sim, to_physical, CONFIG
from transient.sim_cloudy_spectrum import load_grid, compute_rays
from transient.sim_cloudy_linemaps import (LINES, ems_xi_curves, geval)


def main():
    sim = load_sim(os.path.join(HERE, 'data', 'sim', 'disk.out1.00012.athdf'))
    phys = to_physical(sim, CONFIG)
    grid = load_grid()
    rays = compute_rays(sim, phys, CONFIG, grid)
    wave = grid[0]
    r, th, ph = sim['r'], sim['th'], sim['ph']
    nph, nth, nr = len(ph), len(th), len(r)
    dA = rays['dA'].reshape(nph, nth)
    F_line = rays['F_line'].reshape(nph, nth, -1)
    r_if = rays['r_if'].reshape(nph, nth)
    matter = rays['matter'].reshape(nph, nth)
    rf_edges = sim['rf']
    drc = np.diff(rf_edges)
    rif_eff = np.where(matter, rf_edges[-1], r_if)

    i = np.argmin(np.abs(np.degrees(ph) - 170.0))   # ringy lower-left
    NSUP = 64
    th_fine = np.linspace(th[0], th[-1], NSUP * nth)
    tpos = np.interp(th_fine, th, np.arange(nth))
    j0 = np.clip(np.floor(tpos).astype(int), 0, nth - 2)
    tw = (tpos - j0)[:, None]
    _, _, curves, ionc, nxi, S, dS = ems_xi_curves(sim, phys, rays, CONFIG)

    key = 'halpha'
    tex, w1, w2 = LINES[key]
    m = (wave >= w1) & (wave <= w2)
    L = dA * np.trapezoid(F_line[:, :, m] / wave[None, None, m],
                          wave[m], axis=2)
    rif_f = np.interp(th_fine, th, rif_eff[i])
    trunc = np.clip((rif_f[:, None] - rf_edges[None, :-1])
                    / drc[None, :], 0.0, 1.0)
    Sf = (1 - tw) * S[i, j0, :] + tw * S[i, j0 + 1, :]
    dSf = (1 - tw) * dS[i, j0, :] + tw * dS[i, j0 + 1, :]
    idxr = np.clip(np.searchsorted(r, rif_f) - 1, 0, nr - 2)
    trr = np.clip((rif_f - r[idxr]) / (r[idxr + 1] - r[idxr]), 0, 1)
    Send_f = (np.take_along_axis(Sf, idxr[:, None], 1)[:, 0] * (1 - trr)
              + np.take_along_axis(Sf, (idxr + 1)[:, None], 1)[:, 0] * trr)
    Send_f = np.maximum(Send_f, 1e-300)[:, None]
    xh = np.clip(Sf / Send_f, 0.0, 1.0)
    xl = np.clip((Sf - dSf) / Send_f, 0.0, 1.0)
    cf = (1 - tw) * curves[key][i, j0, :] + tw * curves[key][i, j0 + 1, :]
    prof = np.maximum(geval(cf, xh, nxi) - geval(cf, xl, nxi), 0) * trunc
    s = prof.sum(axis=1)
    prof[s <= 0, 0] = 1.0
    prof /= prof.sum(axis=1)[:, None]
    Lf = np.interp(th_fine, th, L[i]) / NSUP
    row = Lf @ prof                                  # map radial profile

    lr = np.log10(row + 1e-30)
    pk = [k for k in range(1, nr - 1)
          if lr[k] > lr[k - 1] and lr[k] >= lr[k + 1]
          and lr[k] > lr.max() - 3]
    print(f"azimuth {np.degrees(ph[i]):.1f} deg, line {key}")
    print(f"map profile peaks at r = "
          f"{np.array2string(r[pk], precision=3)}")
    up = dA[i] > 0
    rifs = np.sort(np.unique(np.round(rif_eff[i][up & ~matter[i]], 3)))
    print(f"per-theta front radii (radiation-bounded, upper): \n"
          f"{np.array2string(rifs, precision=3)}")
    match = [np.min(np.abs(rifs - rp)) < 0.02 for rp in r[pk]]
    print(f"peaks within 0.02 r0 of an individual front: "
          f"{sum(match)}/{len(match)}")
    # sharpness: biggest single-cell log10 jumps in the profile
    dlr = np.diff(lr)
    big = np.argsort(np.abs(dlr))[-8:][::-1]
    print("largest cell-to-cell log10 jumps:")
    for b in sorted(big):
        print(f"  r = {0.5*(r[b]+r[b+1]):.3f}: dlog = {dlr[b]:+.2f}")


if __name__ == '__main__':
    main()
