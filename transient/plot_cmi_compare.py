#!/usr/bin/env python3
"""plot_cmi_compare.py - compare the CMacIonize Monte Carlo
photoionization result with our photon-counting march, both on the SAME
Cartesian density grid (transient/data/cmi_test/grid.npz), so that any
difference comes from the physics (diffuse field, Monte Carlo transfer)
and not from resolution.

Our march: straight rays from the lamp at the origin, nearest-cell
density sampling, cumulative recombinations per steradian
n^2 alpha_B s^2 ds (case B, alpha_B = 2.59e-13, as in the pipeline)
compared with the photon budget Q/4pi. A cell is "ionized" by the march
if the budget is not yet exhausted when the ray reaches it. CMacIonize:
a cell is "ionized" if its neutral hydrogen fraction is below 0.5.

Outputs (in transient/plots):
  cmi_slice_compare.png    vertical slice at the star azimuth: CMacIonize
                           neutral fraction with both ionization fronts
  cmi_front_height.png     face-on map of the CMacIonize front height
                           (top of the neutral gas in each column, upper
                           half)
  cmi_front_height_diff.png  face-on map of front height, CMacIonize
                           minus our march

Run:  python3 -m transient.plot_cmi_compare [--snapshot disk_020.hdf5]
"""

import argparse
import os
import sys

import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from transient.sim_spectrum import (  # noqa: E402
    load_sim, star_position, add_colorbar, CONFIG, PLOTDIR, LD_CM)
from transient.sim_cloudy_spectrum import ALPHA_B  # noqa: E402

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

RUNDIR = os.path.join(HERE, 'data', 'cmi_test')
R_IN = 0.5          # inner edge of the simulated volume [r0]


def load_cmi(path, x, y, z):
    """CMacIonize neutral fraction reshaped onto the (x, y, z) grid of
    grid.npz, using the stored cell coordinates (order-independent)."""
    r0_m = CONFIG['R0_LD'] * LD_CM / 100.0
    with h5py.File(path, 'r') as f:
        c = f['PartType0/Coordinates'][:] / r0_m
        xh = f['PartType0/NeutralFractionH'][:]
        box = f['Header'].attrs['BoxSize'] / r0_m
    # Gadget output coordinates run from 0 to BoxSize; shift to centred
    c = c - 0.5 * box[None, :]
    idx = [np.clip(np.round((c[:, k] - g[0]) / (g[1] - g[0])).astype(int),
                   0, g.size - 1) for k, g in enumerate((x, y, z))]
    out = np.full((x.size, y.size, z.size), np.nan)
    out[idx[0], idx[1], idx[2]] = xh
    return out


def nearest(g, v):
    return np.clip(np.round((v - g[0]) / (g[1] - g[0])).astype(int),
                   0, g.size - 1)


def march_mask(x, y, z, nH, q_ion, r0_cm, n_phi=720, n_th=360, n_s=900):
    """3D ionized mask (upper half) from our straight-ray march on the
    Cartesian grid: each direction marks the cells it passes while its
    photon budget lasts."""
    budget = q_ion / (4 * np.pi)
    smax = np.sqrt(x[-1] ** 2 + y[-1] ** 2 + z[-1] ** 2)
    s = np.linspace(1e-3, smax, n_s)
    ds_cm = (s[1] - s[0]) * r0_cm
    th = np.linspace(0.02, np.pi / 2, n_th)           # upper half only
    ion = np.zeros(nH.shape, dtype=bool)
    for ph in np.linspace(0, 2 * np.pi, n_phi, endpoint=False):
        X = s[None, :] * np.sin(th)[:, None] * np.cos(ph)
        Y = s[None, :] * np.sin(th)[:, None] * np.sin(ph)
        Z = s[None, :] * np.cos(th)[:, None]
        inb = ((np.abs(X) <= x[-1]) & (np.abs(Y) <= y[-1])
               & (Z <= z[-1]))
        ix, iy, iz = nearest(x, X), nearest(y, Y), nearest(z, Z)
        n = np.where(inb, nH[ix, iy, iz], 0.0)
        S = np.cumsum(n ** 2 * ALPHA_B * (s[None, :] * r0_cm) ** 2 * ds_cm,
                      axis=1)
        lit = inb & (S < budget)
        ion[ix[lit], iy[lit], iz[lit]] = True
    return ion


def front_height(ionized, z):
    """Top of the neutral gas in each (x, y) column, upper half: the
    largest z of a non-ionized cell (NaN where the column is fully
    ionized)."""
    up = z > 0
    neutral = ~ionized[:, :, up]
    zu = z[up]
    has = neutral.any(axis=2)
    ktop = zu.size - 1 - np.argmax(neutral[:, :, ::-1], axis=2)
    return np.where(has, zu[ktop], np.nan)


def style(ax):
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=13)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--snapshot', default='disk_020.hdf5')
    ap.add_argument('--dump', default=os.path.join(
        HERE, 'data', 'sim', 'disk.out1.00012.athdf'))
    a = ap.parse_args()

    g = np.load(os.path.join(RUNDIR, 'grid.npz'))
    x, y, z, nH = g['x'], g['y'], g['z'], g['nH']
    q_ion = float(g['q_ion'])
    r0_cm = CONFIG['R0_LD'] * LD_CM
    xh = load_cmi(os.path.join(RUNDIR, a.snapshot), x, y, z)
    print(f"CMacIonize: {np.isnan(xh).sum()} unfilled cells; neutral "
          f"fraction median {np.nanmedian(xh):.3e}")
    ion_cmi = xh < 0.5
    ion_ours = march_mask(x, y, z, nH, q_ion, r0_cm)
    R3 = np.sqrt(x[:, None, None] ** 2 + y[None, :, None] ** 2)
    gas = (R3 >= R_IN) & (z[None, None, :] > 0)
    agree = (ion_cmi == ion_ours)[gas].mean()
    print(f"cells in the upper simulated volume classified the same "
          f"(ionized vs neutral): {100 * agree:.1f}%")
    print(f"ionized fraction of those cells: CMacIonize "
          f"{100 * ion_cmi[gas].mean():.1f}%, our march "
          f"{100 * ion_ours[gas].mean():.1f}%")

    sim = load_sim(a.dump)
    xs, ys = star_position(sim)
    phs = np.arctan2(ys, xs)
    rstar = np.hypot(xs, ys)
    os.makedirs(PLOTDIR, exist_ok=True)

    # ---- 1. vertical slice at the star azimuth ----
    Rg = np.linspace(0.0, x[-1], 600)
    Zg = np.linspace(0.0, z[-1], 260)
    RR, ZZ = np.meshgrid(Rg, Zg)
    ix = nearest(x, RR * np.cos(phs))
    iy = nearest(y, RR * np.sin(phs))
    iz = nearest(z, ZZ)
    sl = xh[ix, iy, iz]
    sl_ours = ion_ours[ix, iy, iz].astype(float)
    sl_cmi = ion_cmi[ix, iy, iz].astype(float)
    fig, ax = plt.subplots(figsize=(12, 5.2))
    pc = ax.pcolormesh(RR, ZZ, np.clip(sl, 1e-6, 1.0),
                       norm=LogNorm(1e-6, 1.0), cmap='viridis',
                       shading='auto')
    add_colorbar(pc, ax, r'$x_{\rm H\,I}~\rm (CMacIonize)$', fontsize=14)
    ax.contour(RR, ZZ, sl_cmi, levels=[0.5], colors='white',
               linewidths=2.5)
    ax.contour(RR, ZZ, sl_ours, levels=[0.5], colors='crimson',
               linewidths=2.5, linestyles='--')
    ax.plot([], [], color='white', lw=2.5,
            label=r'$\rm CMacIonize~front~(diffuse~field)$')
    ax.plot([], [], color='crimson', lw=2.5, ls='--',
            label=r'$\rm our~march~(case~B)$')
    ax.plot(0, 0, marker='*', ms=26, mfc='gold', mec='black', mew=1.5,
            ls='none', zorder=6, clip_on=False)
    ax.plot(rstar, 0, marker='*', ms=18, mfc='white', mec='black',
            mew=1.5, ls='none', zorder=6, clip_on=False)
    ax.set_xlim(0, x[-1])
    ax.set_ylim(0, z[-1])
    ax.set_xlabel(r'$R~[r_0]$', fontsize=17)
    ax.set_ylabel(r'$z~[r_0]$', fontsize=17)
    ax.legend(fontsize=13, loc='upper left', frameon=False,
              labelcolor='white')
    style(ax)
    out = os.path.join(PLOTDIR, 'cmi_slice_compare.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")

    # ---- 2. face-on front height (CMacIonize) and difference ----
    zf_cmi = front_height(ion_cmi, z)
    zf_ours = front_height(ion_ours, z)
    Rxy = np.sqrt(x[:, None] ** 2 + y[None, :] ** 2)
    outside = (Rxy < R_IN) | (Rxy > x[-1])
    zf_cmi[outside] = np.nan
    zf_ours[outside] = np.nan
    X2, Y2 = np.meshgrid(x, y, indexing='ij')
    for name, arr, cmap, norm, lab in [
            ('cmi_front_height.png', zf_cmi, 'magma', None,
             r'$z_{\rm front}~[r_0]~\rm (CMacIonize)$'),
            ('cmi_front_height_diff.png', zf_cmi - zf_ours, 'RdBu_r',
             TwoSlopeNorm(0.0, -0.1, 0.1),
             r'$z_{\rm front,\,CMacIonize} - z_{\rm front,\,march}~[r_0]$')]:
        fig, ax = plt.subplots(figsize=(9.5, 8))
        pc = ax.pcolormesh(X2, Y2, arr, cmap=cmap, norm=norm,
                           shading='auto')
        ax.set_aspect('equal')
        add_colorbar(pc, ax, lab, fontsize=14)
        ax.plot(xs, ys, marker='*', ms=22, mfc='white', mec='black',
                mew=1.5, ls='none', zorder=5)
        ax.set_xlabel(r'$x~[r_0]$', fontsize=16)
        ax.set_ylabel(r'$y~[r_0]$', fontsize=16)
        style(ax)
        out = os.path.join(PLOTDIR, name)
        plt.savefig(out, dpi=200, bbox_inches='tight')
        plt.close()
        print(f"Saved {out}")
    d = (zf_cmi - zf_ours)[~np.isnan(zf_cmi - zf_ours)]
    print(f"front height difference (CMacIonize - march): median "
          f"{np.median(d):+.4f} r0, 10th/90th percentiles "
          f"{np.percentile(d, 10):+.4f} / {np.percentile(d, 90):+.4f} r0")


if __name__ == '__main__':
    main()
