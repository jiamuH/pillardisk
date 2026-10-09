#!/usr/bin/env python3
"""
plot_pillar_xy.py - top-down (x-y plane) view of the transient pillar
geometry: surface temperature map (the TDE-heated spiral stands out), height
contours of the raised arm, and the region of the disk whose sight line to
the observer is occulted (opaque mode) or absorbed (balmer_abs mode).

The temperature map and the height contours are inclination-independent
(they are the physical disk surface); the occulted region depends on the
viewing angle, so one figure is rendered per inclination in INCLINATIONS.

Run:  python3 transient/plot_pillar_xy.py [config_transient.yaml]
"""

import copy
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from transient.transient_disk import disks_from_config  # noqa: E402

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

HERE = os.path.dirname(os.path.abspath(__file__))
INCLINATIONS = [0.0, 45.0, 70.0, 80.0]


def occulted_region(cfg, inc_deg):
    """Visibility field (1 = visible, 0 = fully occulted/absorbed) at the
    given inclination, on the disk grid."""
    cfg_i = copy.deepcopy(cfg)
    cfg_i['observation']['inclination'] = inc_deg
    disk, _, _ = disks_from_config(cfg_i, nr=400, nphi=720)
    if disk.occult_mode == 'balmer_abs':
        col_ref = 2.0 * max(pp['sigma_r'] for pp in disk.pillars)
        vis = np.exp(-disk.tau_edge * disk.compute_los_column() / col_ref)
    else:
        vis = disk.compute_observer_occultation()
    return vis


def main():
    cfg_path = sys.argv[1] if len(sys.argv) > 1 else \
        os.path.join(HERE, 'config_transient.yaml')
    with open(cfg_path) as f:
        cfg = yaml.safe_load(f)

    # inclination-independent fields: temperature and arm height
    disk, _, _ = disks_from_config(cfg, nr=400, nphi=720)
    r_2d, phi_2d = np.meshgrid(disk.r, disk.phi, indexing='ij')
    x = r_2d * np.cos(phi_2d)
    y = r_2d * np.sin(phi_2d)
    T = disk.surface_temperature_map(r_2d, phi_2d)  # includes TDE heating
    h = disk.get_height(r_2d, phi_2d)
    hb = np.interp(r_2d.ravel(), disk.r, disk.h_base).reshape(r_2d.shape)
    h_arm = h - hb
    hmax = float(h_arm.max())
    levels = [round(hmax * f, 2) for f in (0.3, 0.55, 0.8)]

    def close(a):
        return np.concatenate([a, a[:, :1]], axis=1)
    xc, yc, Tc, hac = close(x), close(y), close(T), close(h_arm)

    absorbing = disk.occult_mode == 'balmer_abs'
    cyan_note = (r'$\rm cyan~dashed:~sight~lines~absorbed~by~the~arm$'
                 if absorbing else
                 r'$\rm cyan~dashed:~inner~disk~occulted~by~the~arm$')

    for inc in INCLINATIONS:
        vis = close(occulted_region(cfg, inc))

        from matplotlib.colors import LogNorm
        fig, ax = plt.subplots(figsize=(10, 9))
        im = ax.pcolormesh(xc, yc, Tc, cmap='inferno', shading='auto',
                           norm=LogNorm(vmin=3e3, vmax=3e4))
        cb = plt.colorbar(im, ax=ax, shrink=0.8, aspect=20)
        cb.set_label(r'$\rm surface~temperature~[K]$', fontsize=16)
        cb.ax.tick_params(direction='in', labelsize=13)

        ax.contour(xc, yc, hac, levels=levels, colors='white', linewidths=1.5)
        # occulted region (visibility below 0.5): shade + dashed boundary
        if vis.min() < 0.5:
            ax.contourf(xc, yc, vis, levels=[0.0, 0.5], colors='cyan',
                        alpha=0.20)
            ax.contour(xc, yc, vis, levels=[0.5], colors='cyan',
                       linewidths=2, linestyles='dashed')

        lim = 8.0
        ax.annotate('', xy=(lim - 0.4, 0), xytext=(lim - 2.4, 0),
                    arrowprops=dict(color='limegreen', width=2.5, headwidth=10))
        ax.text(lim - 2.4, 0.45, r'$\rm to~observer$', color='limegreen',
                fontsize=15)
        ax.text(0.03, 0.97, rf'$i = {inc:.0f}^\circ$',
                transform=ax.transAxes, va='top', fontsize=20, color='white')
        ax.text(0.03, 0.90,
                rf'$\rm white~contours:~arm~height~'
                rf'({levels[0]},~{levels[1]},~{levels[2]}~ld)$',
                transform=ax.transAxes, va='top', fontsize=13, color='white')
        ax.text(0.03, 0.85, cyan_note,
                transform=ax.transAxes, va='top', fontsize=13, color='cyan')

        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect('equal')
        ax.set_xlabel(r'$x~[\rm light~days]$', fontsize=18)
        ax.set_ylabel(r'$y~[\rm light~days]$', fontsize=18)
        ax.tick_params(which='major', direction='in', length=8, width=1.5,
                       top=True, right=True, labelsize=14)
        ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                       top=True, right=True)
        ax.minorticks_on()

        out = os.path.join(HERE, 'plots', f'transient_pillar_xy_i{inc:.0f}.png')
        os.makedirs(os.path.dirname(out), exist_ok=True)
        plt.savefig(out, dpi=200, bbox_inches='tight')
        plt.close()
        print(f"Saved {out}")


if __name__ == '__main__':
    main()
