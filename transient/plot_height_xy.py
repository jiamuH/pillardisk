#!/usr/bin/env python3
"""
plot_height_xy.py - top-down (x-y plane) map of the disk + pillar HEIGHT
above the midplane, showing the two-armed spiral shape of the sheared
TDE debris.

Run:  python3 transient/plot_height_xy.py [config_transient.yaml]
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from transient.transient_disk import TransientPillarDisk  # noqa: E402

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    cfg_path = sys.argv[1] if len(sys.argv) > 1 else \
        os.path.join(HERE, 'config_transient.yaml')
    with open(cfg_path) as f:
        cfg = yaml.safe_load(f)

    from transient.transient_disk import disks_from_config
    disk, _, _ = disks_from_config(cfg, nr=400, nphi=720)
    p = cfg['transient']['pillar']

    r_2d, phi_2d = np.meshgrid(disk.r, disk.phi, indexing='ij')
    x = r_2d * np.cos(phi_2d)
    y = r_2d * np.sin(phi_2d)
    h = disk.get_height(r_2d, phi_2d)

    def close(a):
        return np.concatenate([a, a[:, :1]], axis=1)
    x, y, h = close(x), close(y), close(h)

    fig, ax = plt.subplots(figsize=(10, 9))
    im = ax.pcolormesh(x, y, h, cmap='viridis', shading='auto',
                       vmin=0.0, vmax=p['height'])
    cb = plt.colorbar(im, ax=ax, shrink=0.8, aspect=20)
    cb.set_label(r'$\rm height~above~midplane~[light~days]$', fontsize=16)
    cb.ax.tick_params(direction='in', labelsize=13)

    lim = 10.0
    ax.plot(p['r_pillar'] * np.cos(p['phi_pillar']),
            p['r_pillar'] * np.sin(p['phi_pillar']), '*', color='red',
            markersize=18, markeredgecolor='white', markeredgewidth=0.8)
    ax.annotate('', xy=(lim - 0.5, 0), xytext=(lim - 3.0, 0),
                arrowprops=dict(color='limegreen', width=2.5, headwidth=10))
    ax.text(lim - 3.0, 0.55, r'$\rm to~observer$', color='limegreen',
            fontsize=15)
    ax.text(0.03, 0.97, r'$\rm star:~embedded~star~position$',
            transform=ax.transAxes, va='top', fontsize=13, color='white')

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

    out = os.path.join(HERE, 'plots', 'transient_height_xy.png')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")


if __name__ == '__main__':
    main()
