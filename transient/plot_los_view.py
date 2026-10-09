#!/usr/bin/env python3
"""
plot_los_view.py - observer's line-of-sight view of the disk: the surface is
projected onto the sky plane as seen from the inclined viewing direction
e = (sin i, 0, cos i), depth-sorted (painter's algorithm) so the raised
spiral arms hide the surface behind them, and colored by surface
temperature. Cells whose sight line passes through the absorbing arm are
dimmed by the transmission exp(-tau) at the Balmer edge.

Renders one figure per inclination in INCLINATIONS (the last one is the
adopted viewing angle from the config).

Run:  python3 transient/plot_los_view.py [config_transient.yaml]
"""

import copy
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from transient.transient_disk import disks_from_config  # noqa: E402

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

HERE = os.path.dirname(os.path.abspath(__file__))
INCLINATIONS = [0.0, 45.0, 70.0, 80.0]


def render(cfg, inc_deg):
    cfg_i = copy.deepcopy(cfg)
    cfg_i['observation']['inclination'] = inc_deg
    disk, _, _ = disks_from_config(cfg_i, nr=500, nphi=1000)

    r_2d, phi_2d = np.meshgrid(disk.r, disk.phi, indexing='ij')
    h_2d = disk.get_height(r_2d, phi_2d)
    T_2d = disk.surface_temperature_map(r_2d, phi_2d)  # includes TDE heating
    x = r_2d * np.cos(phi_2d)
    y = r_2d * np.sin(phi_2d)
    z = h_2d

    # facing-the-observer cut (full surface normal, both slopes)
    nx, ny, nz = disk._surface_normal(r_2d, phi_2d, h_2d)
    dot = nx * disk.ex + ny * disk.ey + nz * disk.ez

    # dimming of each cell as seen by the observer: geometric visibility
    # (opaque mode) or Balmer-edge transmission (balmer_abs mode)
    if disk.occult_mode == 'balmer_abs':
        col_ref = 2.0 * max(p['sigma_r'] for p in disk.pillars)
        trans = np.exp(-disk.tau_edge * disk.compute_los_column() / col_ref)
        dim_note = r'$\rm dimmed:~disk~seen~through~the~arm$' \
                   r'$\rm ~(Balmer\mbox{-}edge~transmission)$'
    else:
        trans = disk.compute_observer_occultation()

    # sky-plane coordinates: X along +y, Y along (-cos i, 0, sin i)
    x_sky = y
    y_sky = -x * disk.cosi + z * disk.sini
    depth = x * disk.sini + z * disk.cosi  # larger = closer to observer

    # simple painter's render: draw front-facing cells far-to-near, so the
    # near-side arm (drawn last) covers the inner disk it hides.
    lim = 6.0
    yfac = float(np.clip(disk.cosi + 0.28, 0.4, 1.0))
    vis = (dot > 0.0)
    order = np.argsort(depth[vis].ravel())          # far first, near on top
    xs = x_sky[vis].ravel()[order]
    ys = y_sky[vis].ravel()[order]
    Ts = T_2d[vis].ravel()[order]
    tr = trans[vis].ravel()[order]
    rr = r_2d[vis].ravel()[order]

    from matplotlib.colors import LogNorm
    norm = LogNorm(vmin=3e3, vmax=3e4)
    rgba = cm.inferno(norm(Ts))
    if disk.occult_mode == 'balmer_abs':
        # translucent screen: disk seen dimmed through the arm
        rgba[:, :3] *= (0.15 + 0.85 * tr)[:, None]
    # opaque mode: no dimming; the near-side arm covers the disk it hides.

    # marker area matched to the (log-spaced) local cell size so the surface
    # tiles without gaps
    sizes = np.clip((1.2 * rr) ** 2, 4.0, 300.0)
    fig, ax = plt.subplots(figsize=(11, 3.0 + 8.0 * yfac))
    ax.set_facecolor('black')
    ax.scatter(xs, ys, c=rgba, s=sizes, marker='s', linewidths=0,
               rasterized=True)
    sm = cm.ScalarMappable(norm=norm, cmap='inferno')
    cb = plt.colorbar(sm, ax=ax, shrink=0.8, aspect=20)
    cb.set_label(r'$\rm surface~temperature~[K]$', fontsize=16)
    cb.ax.tick_params(direction='in', labelsize=13)

    p = cfg['transient']['pillar']
    T_tde = p.get('pillar_temp', 0.0)
    T_amb = float(np.interp(p['r_pillar'], disk.r, disk.t_base))
    occ_uv = float(disk.compute_occulted_fraction(np.array([2500.0]))[0])
    ax.text(0.02, 0.97, rf'$i = {inc_deg:.0f}^\circ$',
            transform=ax.transAxes, va='top', fontsize=20, color='white')
    ax.text(0.02, 0.87,
            rf'$r_p = {p["r_pillar"]:.0f}~{{\rm ld}},~'
            rf'h = {p["height"]:.1f}~{{\rm ld}},~'
            rf'\sigma_\phi = {p["sigma_phi"]:.1f},~'
            rf'A = {p["spiral_shear"]:.0f}$',
            transform=ax.transAxes, va='top', fontsize=13, color='white')
    ax.text(0.02, 0.80,
            rf'$T_{{\rm pillar}}^{{\rm peak}} = {T_tde / 1e3:.0f}~{{\rm kK}}$'
            rf'$~~(T_{{\rm disk}}(r_p) = {T_amb / 1e3:.0f}~{{\rm kK}})$',
            transform=ax.transAxes, va='top', fontsize=13, color='white')
    ax.text(0.02, 0.73,
            rf'$\rm occulted~UV~(2500~\AA) = {occ_uv * 100:.0f}\%$',
            transform=ax.transAxes, va='top', fontsize=13, color='white')
    if inc_deg > 0:
        ax.text(0.98, 0.03, r'$\rm near~side$',
                transform=ax.transAxes, ha='right', fontsize=13,
                color='lightgray')
        ax.text(0.98, 0.97, r'$\rm far~side$',
                transform=ax.transAxes, ha='right', va='top', fontsize=13,
                color='lightgray')

    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim * yfac, lim * yfac)
    ax.set_aspect('equal')
    ax.set_xlabel(r'$\rm sky~X~[light~days]$', fontsize=18)
    ax.set_ylabel(r'$\rm sky~Y~[light~days]$', fontsize=18)
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=14, color='gray')
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True, color='gray')
    ax.minorticks_on()

    out = os.path.join(HERE, 'plots',
                       f'transient_los_view_i{inc_deg:.0f}.png')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")


def main():
    cfg_path = sys.argv[1] if len(sys.argv) > 1 else \
        os.path.join(HERE, 'config_transient.yaml')
    with open(cfg_path) as f:
        cfg = yaml.safe_load(f)
    for inc in INCLINATIONS:
        render(cfg, inc)


if __name__ == '__main__':
    main()
