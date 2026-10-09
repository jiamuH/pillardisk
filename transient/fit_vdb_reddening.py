#!/usr/bin/env python3
"""
fit_vdb_reddening.py - Ma+26 reddening test with the Vanden Berk (2001)
quasar composite as the intrinsic template (as they use a composite):

    f_obs(lambda) = N (lambda/3000)^alpha VdB(lambda) 10^(-0.4 A_V xi(lambda))

VdB already contains the small blue bump (Fe II + Balmer continuum) and the
broad lines. We fit N, tilt alpha, A_V, and small-grain steepness c2 to the
2019 spectrum AND (separately) to the 2000 spectrum, on line-masked
continuum points, and see whether reddening a normal quasar reproduces the
2019 bump.

Run:  python3 transient/fit_vdb_reddening.py
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from transient.transient_disk import TransientPillarDisk as TPD  # noqa: E402
from transient.make_transient_figure import (  # noqa: E402
    load_epoch, rebin_R, LINE_WINDOWS, PLOTDIR)

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

Z = 0.494
VDB = '/Users/jiamuh/python/lrd_theory/data/vanden_berk_2001.txt'


def main():
    v = np.loadtxt(VDB)
    vw, vf = v[:, 0], v[:, 1]

    reb = {}
    for n in ('sdss2000', 'eboss2019'):
        ep = load_epoch(n)
        wb, fb, eb = rebin_R(ep['wave_obs'], ep['flam'], ep['ivar'], R=500.0)
        reb[n] = (wb / (1.0 + Z), fb, eb)
    grid = reb['eboss2019'][0]
    sel = (grid > 2450) & (grid < 6200)
    wl = grid[sel]
    vdb = np.interp(wl, vw, vf)
    good = np.ones_like(wl, dtype=bool)
    for w1, w2 in LINE_WINDOWS:
        good &= ~((wl >= w1) & (wl <= w2))

    def fit_epoch(f, e):
        def model(p):
            lnN, al, av, c2 = p
            xi = np.array([TPD.smallgrain_extinction_shape(w, 27., c2, 5.5, 0.)
                           for w in wl])
            return np.exp(lnN) * (wl / 3000.0) ** al * vdb \
                * 10.0 ** (-0.4 * av * xi)
        sol = least_squares(lambda p: ((model(p) - f) / e)[good],
                            [np.log(np.median(f)), -0.5, 0.5, 4.0],
                            bounds=([-10, -3, -1, 1], [10, 3, 5, 12]))
        chi2 = np.sum(((model(sol.x) - f) / e)[good] ** 2) / (good.sum() - 4)
        return chi2, sol.x, model(sol.x)

    f19, e19 = reb['eboss2019'][1][sel], reb['eboss2019'][2][sel]
    f00, e00 = np.interp(wl, reb['sdss2000'][0], reb['sdss2000'][1]), \
        np.interp(wl, reb['sdss2000'][0], reb['sdss2000'][2])
    c19, p19, m19 = fit_epoch(f19, e19)
    c00, p00, m00 = fit_epoch(f00, e00)
    print(f"2019 = reddened VdB:  chi2/dof={c19:.1f}  "
          f"alpha={p19[1]:+.2f} A_V={p19[2]:+.2f} c2={p19[3]:.1f}")
    print(f"2000 = reddened VdB:  chi2/dof={c00:.1f}  "
          f"alpha={p00[1]:+.2f} A_V={p00[2]:+.2f} c2={p00[3]:.1f}")

    fig, ax = plt.subplots(figsize=(12, 7))
    ax.fill_between(wl, f19 - e19, f19 + e19, step='mid', color='black',
                    alpha=0.15, lw=0)
    ax.plot(wl, f19, drawstyle='steps-mid', color='black', lw=1.5, alpha=0.85,
            label=r'$\rm SDSS\mbox{-}2019~(data)$')
    ax.plot(wl, f00, drawstyle='steps-mid', color='crimson', lw=1.1,
            alpha=0.6, label=r'$\rm SDSS\mbox{-}2000~(data)$')
    ax.plot(wl, m19, '-', color='royalblue', lw=3, alpha=0.95,
            label=rf'$\rm reddened~VdB~\to~2019~(A_V={p19[2]:.2f},~'
                  rf'\chi^2_\nu={c19:.0f})$')
    ax.plot(wl, m00, '--', color='seagreen', lw=2, alpha=0.9,
            label=rf'$\rm reddened~VdB~\to~2000~(A_V={p00[2]:.2f},~'
                  rf'\chi^2_\nu={c00:.0f})$')
    ax.axvline(3646, color='gray', ls=':', lw=1.5)
    for w1, w2 in LINE_WINDOWS:
        ax.axvspan(w1, w2, color='gray', alpha=0.12, lw=0)
    ax.set_xlabel(r'$\rm rest\mbox{-}frame~wavelength~[\AA]$', fontsize=18)
    ax.set_ylabel(r'$f_\lambda~[\rm 10^{-17}~erg~s^{-1}~cm^{-2}~\AA^{-1}]$',
                  fontsize=17)
    ax.legend(fontsize=12, frameon=False)
    ax.set_xlim(2450, 6200)
    ax.set_ylim(0, None)
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=14)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()
    out = os.path.join(PLOTDIR, 'transient_vdb_reddening.png')
    os.makedirs(PLOTDIR, exist_ok=True)
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")


if __name__ == '__main__':
    main()
