#!/usr/bin/env python3
"""
plot_sim_ifront_physical.py - the ionization front rendered as a physical
surface: side view in cylindrical (R, z), lamp at the origin. For each
azimuth phi, the front along every lamp ray (polar angle theta) sits at
radius r_IF(theta, phi); that point is drawn at
    R = r_IF sin(theta),  z = r_IF cos(theta).
Each azimuth gives one wall-profile curve (colored by phi, cyclic map);
only the radiation-bounded part is drawn - directions where the light
escapes the box entirely are left blank. The star's azimuth is
highlighted. This is the sim's version of the toy-model bowl + arm-wall
cartoon: everything outside a curve (at its azimuth) is in shadow.

Run:  python3 transient/plot_sim_ifront_physical.py [dump.athdf]
"""

import os
import sys

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from transient.sim_spectrum import (  # noqa: E402
    load_sim, to_physical, add_colorbar, star_position, CONFIG, PLOTDIR)
from transient.sim_cloudy_spectrum import load_grid, compute_rays  # noqa: E402

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else \
        os.path.join(HERE, 'data', 'sim', 'disk.out1.00012.athdf')
    cfg = CONFIG
    sim = load_sim(path)
    phys = to_physical(sim, cfg)
    rays = compute_rays(sim, phys, cfg, load_grid())
    r_if = rays['r_if']                       # (ph, th)
    matter = rays['ion'][:, :, -1]            # light escapes -> no wall
    th, ph = sim['th'], sim['ph']

    R = np.where(matter, np.nan, r_if) * np.sin(th)[None, :]
    Z = np.where(matter, np.nan, r_if) * np.cos(th)[None, :]

    fig, ax = plt.subplots(figsize=(12, 6.5))
    cmap = mpl.cm.twilight
    for k in range(len(ph)):
        ax.plot(R[k], Z[k], '-', color=cmap(ph[k] / (2 * np.pi)),
                lw=1.0, alpha=0.45)
    # star's azimuth highlighted
    xs, ys = star_position(sim)
    phs = np.arctan2(ys, xs) % (2 * np.pi)
    ks = np.argmin(np.abs(ph - phs))
    ax.plot(R[ks], Z[ks], '-', color='black', lw=2.8, alpha=0.95,
            label=r'$\rm star~azimuth$')
    # geometry guides: wedge boundaries, inner/outer box radii, lamp
    tt = np.linspace(th.min(), th.max(), 200)
    for rr, style in [(sim['rf'][0], '--'), (sim['rf'][-1], '--')]:
        ax.plot(rr * np.sin(tt), rr * np.cos(tt), style, color='gray',
                lw=1.2, alpha=0.8)
    for tb in (th.min(), th.max()):
        ax.plot([0, 2.15 * np.sin(tb)], [0, 2.15 * np.cos(tb)], ':',
                color='gray', lw=1.2)
    ax.plot(0, 0, marker='o', ms=13, mfc='gold', mec='black', mew=1.5)
    ax.text(0.045, 0.02, r'$\rm lamp$', fontsize=14)
    ax.text(0.53, -0.28, r'$\rm inner~rim$', fontsize=12, color='gray',
            rotation=90)
    ax.text(1.97, -0.28, r'$\rm outer~box$', fontsize=12, color='gray',
            rotation=75)
    ax.set_xlabel(r'$R~[r_0]$', fontsize=18)
    ax.set_ylabel(r'$z~[r_0]$', fontsize=18)
    ax.set_xlim(0, 2.2)
    ax.set_ylim(-0.72, 0.72)
    ax.set_aspect('equal')
    ax.legend(fontsize=13, frameon=False, loc='upper left')
    sm = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(0, 360),
                               cmap=cmap)
    cb = add_colorbar(sm, ax, r'$\rm azimuth~\phi~[deg]$')
    cb.set_ticks([0, 90, 180, 270, 360])
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=14)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()
    os.makedirs(PLOTDIR, exist_ok=True)
    out = os.path.join(PLOTDIR, 'sim_ifront_physical.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")


if __name__ == '__main__':
    main()
