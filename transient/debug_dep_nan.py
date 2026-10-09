#!/usr/bin/env python3
"""Find the NaN source in the demo slice deposition."""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from transient.sim_spectrum import load_sim, to_physical, star_position, CONFIG
from transient.sim_cloudy_spectrum import load_grid, compute_rays
from transient.sim_cloudy_linemaps import LINES, ems_deposition

sim = load_sim(os.path.join(HERE, 'data', 'sim', 'disk.out1.00012.athdf'))
phys = to_physical(sim, CONFIG)
grid = load_grid()
rays = compute_rays(sim, phys, CONFIG, grid)
wave = grid[0]
nph, nth = len(sim['ph']), len(sim['th'])
xs, ys = star_position(sim)
ip = np.argmin(np.abs(sim['ph'] - np.mod(np.arctan2(ys, xs), 2 * np.pi)))
F_line = rays['F_line'].reshape(nph, nth, -1)[ip]
dAp = rays['dA'].reshape(nph, nth)[ip]
depos = ems_deposition(sim, phys, rays, CONFIG)
for key, (tex, w1, w2) in LINES.items():
    m = (wave >= w1) & (wave <= w2)
    L_ray = dAp * np.trapezoid(F_line[:, m] / wave[None, m], wave[m], axis=1)
    dep = L_ray[:, None] * depos[key][ip]
    logl = np.log10(dep + 1e-30)
    print(f"{key:7s}: NaN in dep {np.isnan(dep).sum()}, "
          f"min dep {dep.min():.3e}, max dep {np.nanmax(dep):.3e}, "
          f"NaN in log {np.isnan(logl).sum()}, "
          f"logl.max {np.nanmax(logl):.2f} plainmax {logl.max():.2f}")
