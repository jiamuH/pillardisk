#!/usr/bin/env python3
"""
fit_redden_quiescent.py - Ma+26-style reddening test, but using OUR OWN
quiescent state as the intrinsic quasar template (no generic composite):

    f_2019(lambda) = B * f_quiescent(lambda) * 10^(-0.4 A_V xi(lambda))

The intrinsic is the quiescent (2021 DESI or 2000 SDSS) spectrum, allowed to
brighten by a free factor B (optionally a power-law tilt), then reddened by
the SMC or small-grain (Ma+26) law. Fit B, A_V (and tilt) to the 2019
spectrum. Tests whether the 2019 'bump' is the reddening peak of a
brightened intrinsic continuum.

Run:  python3 transient/fit_redden_quiescent.py
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
AV_GRID = np.arange(-1.0, 2.01, 0.1)     # <0 = de-reddening (dust cleared)
TILT_GRID = np.arange(-1.0, 1.01, 0.25)  # extra power-law tilt on intrinsic
LAWS = {'SMC': None, 'small-grain': (27., 6., 5.5, 0.)}
INTRINSIC = 'sdss2000'   # our quiescent template (2000 or desi2021)


def xi_of(law_c, wl):
    if law_c is None:
        x = np.array([TPD.smc_extinction_shape(w) for w in wl])
        return x / TPD.smc_extinction_shape(5500.0)
    return np.array([TPD.smallgrain_extinction_shape(w, *law_c) for w in wl])


def main():
    reb = {}
    for n in ('sdss2000', 'eboss2019', 'desi2021'):
        ep = load_epoch(n)
        if ep is None:
            continue
        wb, fb, eb = rebin_R(ep['wave_obs'], ep['flam'], ep['ivar'], R=500.0)
        reb[n] = (wb / (1.0 + Z), fb, eb)
    grid = reb['eboss2019'][0]
    sel = (grid > 2450) & (grid < 6200)
    wl = grid[sel]
    f19 = reb['eboss2019'][1][sel]
    e19 = reb['eboss2019'][2][sel]
    fq = np.interp(wl, reb[INTRINSIC][0], reb[INTRINSIC][1])   # intrinsic
    good = np.ones_like(wl, dtype=bool)
    for w1, w2 in LINE_WINDOWS:
        good &= ~((wl >= w1) & (wl <= w2))

    best = None
    for lname, law_c in LAWS.items():
        xi = xi_of(law_c, wl)
        for av, tilt in itertools.product(AV_GRID, TILT_GRID):
            shape = fq * (wl / 3000.0) ** tilt * 10.0 ** (-0.4 * av * xi)
            w = 1.0 / e19[good] ** 2
            B = np.sum(f19[good] * shape[good] * w) / np.sum(shape[good]**2 * w)
            model = B * shape
            chi2 = np.sum(((model[good] - f19[good]) / e19[good]) ** 2) \
                / (good.sum() - 3)
            if best is None or chi2 < best[0]:
                best = (chi2, lname, av, tilt, B, model.copy())
    chi2, lname, av, tilt, B, model = best
    print(f"intrinsic = {INTRINSIC}")
    print(f"BEST reddened-brightened-quiescent fit to 2019:")
    print(f"  chi2/dof = {chi2:.1f}  law={lname}  A_V={av:+.2f}  "
          f"tilt={tilt:+.2f}  brighten B={B:.2f}")

    fig, ax = plt.subplots(figsize=(12, 7))
    ax.fill_between(wl, f19 - e19, f19 + e19, step='mid', color='black',
                    alpha=0.15, lw=0)
    ax.plot(wl, f19, drawstyle='steps-mid', color='black', lw=1.5, alpha=0.85,
            label=r'$\rm SDSS\mbox{-}2019~(data,~the~bump)$')
    ax.plot(wl, fq, drawstyle='steps-mid', color='crimson', lw=1.2, alpha=0.7,
            label=rf'$\rm intrinsic~=~{INTRINSIC}~(quiescent)$')
    ax.plot(wl, B * fq * (wl / 3000.0) ** tilt, '--', color='seagreen', lw=2,
            alpha=0.9, label=rf'$\rm brightened~intrinsic~(B={B:.1f})$')
    ax.plot(wl, model, '-', color='royalblue', lw=3, alpha=0.95,
            label=rf'$\rm reddened~({lname},~A_V={av:+.1f},~'
                  rf'\chi^2_\nu={chi2:.0f})$')
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
    out = os.path.join(PLOTDIR, 'transient_redden_quiescent.png')
    os.makedirs(PLOTDIR, exist_ok=True)
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")


if __name__ == '__main__':
    main()
