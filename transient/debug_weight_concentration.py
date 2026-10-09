#!/usr/bin/env python3
"""How concentrated is the n^2 r^2 dr deposition weight along each
radiation-bounded ray? If most of the weight sits just inside the
ionization front, the 'uniform recombination-weight' painting is closer
to front-weighted than to uniform."""
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
    rays = compute_rays(sim, phys, CONFIG, load_grid())
    r = sim['r']
    nph, nth = len(sim['ph']), len(sim['th'])
    w = rays['w']                                  # (ph, th, r)
    matter = rays['matter']
    dA = rays['dA'].reshape(nph, nth)
    rb = (~matter) & (dA > 0) & (w.sum(axis=2) > 0)  # radiation-bounded, upper

    wsel = w[rb]                                   # (nray_rb, nr)
    r_if = rays['r_if'][rb]
    frac_last = []
    for q in (0.1, 0.2):
        # weight fraction inside the outermost q of the ionized segment
        seg_lo = r[0] + (1 - q) * (r_if - r[0])
        inlast = r[None, :] >= seg_lo[:, None]
        f = (wsel * inlast).sum(axis=1) / wsel.sum(axis=1)
        frac_last.append(f)
        print(f"weight in outermost {int(q*100)}% of the ionized segment: "
              f"median {np.median(f):.2f}, 25th {np.percentile(f, 25):.2f}, "
              f"75th {np.percentile(f, 75):.2f}")
    print(f"rays considered: {rb.sum()} radiation-bounded upper rays")


if __name__ == '__main__':
    main()
