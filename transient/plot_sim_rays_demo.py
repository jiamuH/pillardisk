#!/usr/bin/env python3
"""plot_sim_rays_demo.py - demonstration of the ionization ray march.

A vertical (R, z) slice of the simulation at the azimuth of the
disrupted star, showing ONLY the observer-side (upper) hemisphere, the
one that enters the observed sum. Log density is the background, the
lamp-post ionizing source sits at the origin, and a fan of radial rays
is drawn: solid gold over the ionized segment, dotted gray through the
shadowed gas behind the front, royal blue for matter-bounded rays that
stay transparent to the edge of the box. The crimson curve is the
ionization front r_IF(theta) of this azimuth; the white star marks the
disrupted star (midplane).

Output: plots/sim_rays_demo.png
Run:  python3 transient/plot_sim_rays_demo.py [dump.athdf]
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from transient.sim_spectrum import (  # noqa: E402
    load_sim, to_physical, star_position, add_colorbar, CONFIG, PLOTDIR)
from transient.sim_cloudy_spectrum import load_grid, compute_rays  # noqa: E402

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

N_RAYS = 7           # rays drawn across the upper wedge


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else \
        os.path.join(HERE, 'data', 'sim', 'disk.out1.00012.athdf')
    cfg = CONFIG
    sim = load_sim(path)
    phys = to_physical(sim, cfg)
    grid = load_grid()
    rays = compute_rays(sim, phys, cfg, grid)

    r, th, ph = sim['r'], sim['th'], sim['ph']
    nph, nth = len(ph), len(th)
    r_if = rays['r_if'].reshape(nph, nth)
    matter = rays['matter'].reshape(nph, nth)

    # observer-side (upper) hemisphere only: the lower half is the hidden
    # face of the opaque disk and does not enter the observed sum
    up = th <= np.pi / 2
    thu = th[up]

    # slice at the star's azimuth
    xs, ys = star_position(sim)
    ps = np.mod(np.arctan2(ys, xs), 2 * np.pi)
    ip = np.argmin(np.abs(ph - ps))
    rstar = np.hypot(xs, ys)

    # background: log n_H in the (R, z) plane at this azimuth
    TH, RR = np.meshgrid(thu, r, indexing='ij')
    Rc, Zc = RR * np.sin(TH), RR * np.cos(TH)
    fig, ax = plt.subplots(figsize=(12, 5.2))
    pc = ax.pcolormesh(Rc, Zc, np.log10(phys['nH'][ip][up]),
                       cmap='Greys', shading='auto')
    add_colorbar(pc, ax, r'$\log n_{\rm H}~[\rm cm^{-3}]$', fontsize=14)

    # the ionization front of this azimuth
    ax.plot(r_if[ip, up] * np.sin(thu), r_if[ip, up] * np.cos(thu),
            color='crimson', lw=2.5, alpha=0.9, zorder=4,
            label=r'${\rm ionization~front}~r_{\rm IF}(\theta)$')

    # the fan of demonstration rays
    ju = np.where(up)[0]
    lab_lit = lab_shad = lab_mat = False
    for j in ju[np.unique(np.linspace(0, ju.size - 1,
                                      N_RAYS).astype(int))]:
        s, c = np.sin(th[j]), np.cos(th[j])
        if matter[ip, j]:
            ax.plot([0, r[-1] * s], [0, r[-1] * c], color='royalblue',
                    lw=2.2, alpha=0.85, zorder=3,
                    label=None if lab_mat else
                    r'$\rm matter{-}bounded~ray$')
            lab_mat = True
        else:
            rf = r_if[ip, j]
            ax.plot([0, rf * s], [0, rf * c], color='goldenrod',
                    lw=2.2, alpha=0.9, zorder=3,
                    label=None if lab_lit else
                    r'$\rm ionized~segment$')
            ax.plot([rf * s, r[-1] * s], [rf * c, r[-1] * c],
                    color='dimgray', lw=1.5, ls=':', alpha=0.8, zorder=2,
                    label=None if lab_shad else r'$\rm shadowed$')
            ax.plot(rf * s, rf * c, 'o', ms=7, mfc='crimson',
                    mec='black', mew=0.8, zorder=5)
            lab_lit = lab_shad = True

    # observer direction (inclination from the +z disk axis)
    i_obs = np.radians(CONFIG['INCL_DEG'])
    ox, oz = np.sin(i_obs), np.cos(i_obs)
    ax.annotate('', xytext=(1.42, 0.36), xy=(1.42 + 0.24 * ox,
                                             0.36 + 0.24 * oz),
                arrowprops=dict(arrowstyle='-|>', lw=2.5, color='black'))
    ax.text(1.41, 0.40, r'$\rm to~observer~(i=45^\circ)$',
            fontsize=14, ha='right', va='bottom')

    # the lamp and the star
    ax.plot(0, 0, marker='*', ms=30, mfc='gold', mec='black', mew=1.5,
            ls='none', zorder=6, clip_on=False)
    ax.annotate(r'$\rm lamp$', (0, 0), xytext=(0.05, 0.05),
                textcoords='data', fontsize=16)
    ax.plot(rstar, 0, marker='*', ms=20, mfc='white', mec='black',
            mew=1.5, ls='none', zorder=6, clip_on=False)
    ax.annotate(r'$\rm star$', (rstar, 0), xytext=(rstar + 0.04, 0.04),
                textcoords='data', fontsize=16, color='white')

    ax.set_xlim(-0.08, r[-1] * 1.02)
    zmax = (r[-1] * np.cos(thu)).max()
    ax.set_ylim(-0.04, zmax * 1.30)
    ax.set_xlabel(r'$R~[r_0]$', fontsize=17)
    ax.set_ylabel(r'$z~[r_0]$', fontsize=17)
    ax.legend(fontsize=13, loc='upper left', frameon=False)
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=13)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()
    os.makedirs(PLOTDIR, exist_ok=True)
    out = os.path.join(PLOTDIR, 'sim_rays_demo.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")

    # ---- second figure: the same slice, but the emitting gas painted
    # with the three-line composite (R = Halpha, G = MgII, B = CIV).
    # Each ray's line luminosities come from its ONE Cloudy slab, and the
    # placement along the ray follows Cloudy's OWN depth-resolved
    # emissivities (the arm_column2 'save lines emissivity' output,
    # mapped by cumulative recombination fraction; see ems_deposition in
    # sim_cloudy_linemaps.py). Per-line totals are conserved. ----
    from transient.sim_cloudy_linemaps import (LINES as LINE_WINS,
                                               ems_deposition)
    nw = len(grid[0])
    wave = grid[0]
    F_line = rays['F_line'].reshape(nph, nth, nw)[ip]     # (nth, nw)
    dAp = rays['dA'].reshape(nph, nth)[ip]                # 0 below midplane
    depos = ems_deposition(sim, phys, rays, cfg)
    dep = {}
    for key, (tex, w1, w2) in LINE_WINS.items():
        m = (wave >= w1) & (wave <= w2)
        L_ray = dAp * np.trapezoid(F_line[:, m] / wave[None, m],
                                   wave[m], axis=1)
        dep[key] = L_ray[:, None] * depos[key][ip]

    # sample onto a Cartesian (R, z) grid of the upper wedge
    ngx, ngz = 1000, 420
    xg = np.linspace(0.0, r[-1] * 1.02, ngx)
    zg = np.linspace(0.0, zmax * 1.30, ngz)
    Xg, Zg = np.meshgrid(xg, zg)
    rr = np.hypot(Xg, Zg)
    tt = np.arctan2(Xg, Zg)                               # polar angle
    inside = (rr >= r[0]) & (rr <= r[-1]) & (tt >= th[0]) & (tt <= np.pi/2)
    it = np.clip(np.searchsorted(th, tt), 0, nth - 1)
    ir = np.clip(np.searchsorted(r, rr), 0, len(r) - 1)
    # per-cell composite: hue differences along a ray are now REAL (the
    # U-stratified deposition above), so color per cell. Same recipe as
    # the face-on map: coarse-knot (13 percentile anchors) smooth
    # equalization per channel over a WIDE dynamic range, so the faint
    # emission stays visible, then a chroma boost.
    # wide range (8 dex) and no darkening gamma, so the faint
    # matter-bounded upper layers stay visible next to the front skin
    DR, SAT = 8.0, 1.8
    chan = np.zeros((3, nth, len(r)))
    for k, key in enumerate(['halpha', 'mgii', 'civ']):
        lm = dep[key]
        logl = np.log10(lm + 1e-30)
        lit = logl > logl.max() - DR
        qs = np.linspace(0.0, 100.0, 21)
        anchors = np.percentile(logl[lit], qs)
        anchors = np.maximum.accumulate(anchors)
        anchors += 1e-9 * np.arange(anchors.size)
        chan[k] = np.interp(logl, anchors, qs / 100.0)
        chan[k][~lit] = 0.0
    rgb = chan[:, it, ir].transpose(1, 2, 0)
    rgb[~inside] = 0.0
    lum = rgb.mean(axis=2, keepdims=True)
    rgb = np.clip(lum + SAT * (rgb - lum), 0, 1)

    fig, ax = plt.subplots(figsize=(12, 5.2))
    ax.imshow(rgb, extent=[xg[0], xg[-1], zg[0], zg[-1]],
              origin='lower', aspect='auto', interpolation='bilinear')
    ax.plot(r_if[ip, up] * np.sin(thu), r_if[ip, up] * np.cos(thu),
            color='white', lw=1.5, alpha=0.7, zorder=4)
    ax.plot(0, 0, marker='*', ms=30, mfc='gold', mec='white', mew=1.5,
            ls='none', zorder=6, clip_on=False)
    ax.annotate(r'$\rm lamp$', (0, 0), xytext=(0.05, 0.05),
                textcoords='data', fontsize=16, color='white')
    ax.plot(rstar, 0, marker='*', ms=20, mfc='white', mec='black',
            mew=1.5, ls='none', zorder=6, clip_on=False)
    ax.annotate(r'$\rm star$', (rstar, 0), xytext=(rstar + 0.04, 0.04),
                textcoords='data', fontsize=16, color='white')
    ax.annotate('', xytext=(1.42, 0.36), xy=(1.42 + 0.24 * ox,
                                             0.36 + 0.24 * oz),
                arrowprops=dict(arrowstyle='-|>', lw=2.5, color='white'))
    ax.text(1.41, 0.40, r'$\rm to~observer~(i=45^\circ)$',
            fontsize=14, ha='right', va='bottom', color='white')
    for kk, (txt, col) in enumerate([
            (r'$\rm H\alpha$', 'crimson'),
            (r'$\rm Mg\,II$', 'mediumseagreen'),
            (r'$\rm C\,IV$', 'cornflowerblue')]):
        ax.text(0.02, 0.95 - 0.09 * kk, txt, transform=ax.transAxes,
                va='top', fontsize=17, color=col)
    ax.set_xlim(-0.08, r[-1] * 1.02)
    ax.set_ylim(-0.04, zmax * 1.30)
    ax.set_facecolor('black')
    ax.set_xlabel(r'$R~[r_0]$', fontsize=17)
    ax.set_ylabel(r'$z~[r_0]$', fontsize=17)
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=13, color='white')
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True, color='white')
    ax.minorticks_on()
    out = os.path.join(PLOTDIR, 'sim_rays_demo_rgb.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")


if __name__ == '__main__':
    main()
