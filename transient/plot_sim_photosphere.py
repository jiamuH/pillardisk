#!/usr/bin/env python3
"""
plot_sim_photosphere.py - column-density surface height map: the height z
(from the midplane) at which a VERTICAL column from the upper wedge
boundary into the disk reaches a fixed TOTAL hydrogen column,

    N_H(>z) = int n_H ds  (ds ~ r dtheta)  =  N_thresh.

Pure N_H, no opacity assumed: the region is well inside the dust
sublimation radius (no dust), and electron scattering depends on the
ionization state - so the opacity interpretation is left to the
ionization machinery (compute_rays). The mapped fiducial is
N_H = 1e23 cm^-2; stats for 1e22/1e23/1e24 are printed. These
constant-column surfaces are the sim analog of the toy model's height
surface get_height(r, phi); the Gaussian sigma_theta (plot_sim_hr.py)
remains the thermodynamic diagnostic. Depends linearly on NMID_CM3.

Run:  python3 transient/plot_sim_photosphere.py [dump.athdf]
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from transient.sim_spectrum import (  # noqa: E402
    load_sim, to_physical, add_colorbar, mark_star, CONFIG, LD_CM, PLOTDIR)

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

N_LEVELS = [1e22, 1e23, 1e24]   # cm^-2: fiducial column thresholds
N_MAP = 1e23                    # the mapped surface


def photosphere(sim, phys, cfg, thresh):
    """height z/r of the tau=1 surface per (phi, r), from the top down."""
    r, th = sim['r'], sim['th']
    r0_cm = cfg['R0_LD'] * LD_CM
    dth = np.gradient(th)
    # vertical path element ds = r dtheta; integrate from the wedge top
    ds = (r[None, None, :] * r0_cm) * dth[None, :, None]
    up = th <= np.pi / 2
    Ncum = np.cumsum(phys['nH'][:, up, :] * ds[:, up, :], axis=1)
    above = Ncum >= thresh                          # first True = surface
    j = np.argmax(above, axis=1)                    # (ph, r)
    zr = (np.pi / 2 - th[up][j])                    # z/r in rad
    zr[~above.any(axis=1)] = 0.0                    # never opaque
    return zr


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else \
        os.path.join(HERE, 'data', 'sim', 'disk.out1.00012.athdf')
    cfg = CONFIG
    sim = load_sim(path)
    phys = to_physical(sim, cfg)
    r, ph = sim['r'], sim['ph']

    for nthr in N_LEVELS:
        zr = photosphere(sim, phys, cfg, nthr)
        med, p95 = np.percentile(zr[zr > 0], [50, 95])
        print(f"N_H = {nthr:.0e}: surface z/r median={med:.3f}, "
              f"95%={p95:.3f} (x sigma_theta ~ {med/0.055:.1f} / "
              f"{p95/0.055:.1f})")
    zr_map = photosphere(sim, phys, cfg, N_MAP)

    R, P = np.meshgrid(r, ph)
    X, Y = R * np.cos(P), R * np.sin(P)
    fig, ax = plt.subplots(figsize=(9.5, 8))
    pc = ax.pcolormesh(X, Y, zr_map, cmap='cividis', shading='auto')
    ax.set_aspect('equal')
    add_colorbar(pc, ax, r'$z(N_{\rm H}=10^{23}~\rm cm^{-2})/r$')
    mark_star(ax, sim)
    ax.set_xlabel(r'$x~[r_0]$', fontsize=16)
    ax.set_ylabel(r'$y~[r_0]$', fontsize=16)
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=13)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()
    os.makedirs(PLOTDIR, exist_ok=True)
    out = os.path.join(PLOTDIR, 'sim_photosphere_map.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")


if __name__ == '__main__':
    main()
