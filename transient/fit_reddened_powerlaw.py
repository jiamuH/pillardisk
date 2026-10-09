#!/usr/bin/env python3
"""
fit_reddened_powerlaw.py - test the Ma et al. (2026) red-quasar picture on
the RAW 2019 spectrum (not the difference): is the 2019 'bump' just the
reddening peak of an intrinsic blue power-law?

    f_obs(lambda) = N (lambda/3000)^alpha * 10^(-0.4 A_V xi(lambda))

with xi the SMC or small-grain (Ma+26) extinction curve. Fit N, alpha, A_V
(and the small-grain steepness c2) to the line-masked 2019 continuum. The
2000 spectrum is fit the same way for comparison (a bluer / less-reddened
power-law).

Run:  python3 transient/fit_reddened_powerlaw.py
"""

import itertools
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from transient.transient_disk import TransientPillarDisk as TPD  # noqa: E402
from transient.make_transient_figure import (  # noqa: E402
    load_epoch, rebin_R, LINE_WINDOWS, PLOTDIR)

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

Z = 0.494
AV_GRID = np.arange(0.0, 3.01, 0.15)
ALPHA_GRID = np.arange(-2.6, 0.61, 0.2)   # f_lambda slope
LAWS = {'SMC': None, 'small-grain': (27., 6., 5.5, 0.)}


def xi_of(law_c, wl):
    if law_c is None:
        x = np.array([TPD.smc_extinction_shape(w) for w in wl])
        return x / TPD.smc_extinction_shape(5500.0)
    return np.array([TPD.smallgrain_extinction_shape(w, *law_c) for w in wl])


def fit_epoch(wl, f, e, good):
    """best reddened power-law over (law, A_V, alpha); N solved linearly."""
    best = None
    for lname, law_c in LAWS.items():
        xi = xi_of(law_c, wl)
        for av, al in itertools.product(AV_GRID, ALPHA_GRID):
            shape = (wl / 3000.0) ** al * 10.0 ** (-0.4 * av * xi)
            w = 1.0 / e[good] ** 2
            N = np.sum(f[good] * shape[good] * w) / np.sum(shape[good] ** 2 * w)
            model = N * shape
            chi2 = np.sum(((model[good] - f[good]) / e[good]) ** 2) \
                / (good.sum() - 3)
            if best is None or chi2 < best[0]:
                best = (chi2, lname, av, al, N, model.copy())
    return best


def main():
    b = {}
    for n in ('sdss2000', 'eboss2019'):
        ep = load_epoch(n)
        wb, fb, eb = rebin_R(ep['wave_obs'], ep['flam'], ep['ivar'], R=500.0)
        b[n] = (wb / (1.0 + Z), fb, eb)
    grid = b['sdss2000'][0]
    sel = (grid > 2450) & (grid < 6200)
    wl = grid[sel]
    f00, e00 = b['sdss2000'][1][sel], b['sdss2000'][2][sel]
    f19 = np.interp(wl, b['eboss2019'][0], b['eboss2019'][1])
    e19 = np.interp(wl, b['eboss2019'][0], b['eboss2019'][2])
    good = np.ones_like(wl, dtype=bool)
    for w1, w2 in LINE_WINDOWS:
        good &= ~((wl >= w1) & (wl <= w2))

    best19 = fit_epoch(wl, f19, e19, good)
    best00 = fit_epoch(wl, f00, e00, good)
    print(f"2019: chi2/dof={best19[0]:.1f}  law={best19[1]}  "
          f"A_V={best19[2]:.2f}  alpha={best19[3]:.1f}")
    print(f"2000: chi2/dof={best00[0]:.1f}  law={best00[1]}  "
          f"A_V={best00[2]:.2f}  alpha={best00[3]:.1f}")

    fig, ax = plt.subplots(figsize=(12, 7))
    ax.fill_between(wl, f19 - e19, f19 + e19, step='mid', color='black',
                    alpha=0.15, lw=0)
    ax.plot(wl, f19, drawstyle='steps-mid', color='black', lw=1.4, alpha=0.85,
            label=r'$\rm SDSS\mbox{-}2019~(data)$')
    ax.plot(wl, f00, drawstyle='steps-mid', color='crimson', lw=1.2,
            alpha=0.7, label=r'$\rm SDSS\mbox{-}2000~(data)$')
    ax.plot(wl, best19[5], '-', color='royalblue', lw=3, alpha=0.95,
            label=rf'$\rm reddened~PL~fit~to~2019~({best19[1]},~'
                  rf'A_V={best19[2]:.1f},~\chi^2_\nu={best19[0]:.0f})$')
    # intrinsic (dereddened) power-law behind the 2019 fit
    xi = xi_of(LAWS[best19[1]], wl)
    intrinsic = best19[4] * (wl / 3000.0) ** best19[3]
    ax.plot(wl, intrinsic, '--', color='seagreen', lw=2, alpha=0.9,
            label=r'$\rm intrinsic~power\mbox{-}law~(dereddened)$')
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
    out = os.path.join(PLOTDIR, 'transient_reddened_powerlaw.png')
    os.makedirs(PLOTDIR, exist_ok=True)
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")


if __name__ == '__main__':
    main()
