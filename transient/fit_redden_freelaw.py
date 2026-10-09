#!/usr/bin/env python3
"""
fit_redden_freelaw.py - the Ma+26 reddening test done faithfully:
intrinsic = brightened quiescent (2000), reddened by a FREE-PARAMETER
small-grain law whose c-parameters are FIT (not fixed):

    f_2019 = B * f_2000 * 10^(-0.4 A_V xi(lambda; c1,c2,c3,c4))

Fit B, A_V, c2, c4 (c1, c3 fixed; A_V and c1 are degenerate). Then show the
model-independent diagnostic: the attenuation actually required by the data,
  Delta m(lambda) = -2.5 log10(f_2019 / f_2000),
which ANY reddening model must equal (up to the constant -2.5 log10 B). A
real extinction curve A_V*xi is monotonic (rises to the blue); the required
Delta m is NOT - it has an interior minimum (a 'transparency window') at the
bump, which no dust law of any parameters can produce.

Run:  python3 transient/fit_redden_freelaw.py
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


def main():
    reb = {}
    for n in ('sdss2000', 'eboss2019'):
        ep = load_epoch(n)
        wb, fb, eb = rebin_R(ep['wave_obs'], ep['flam'], ep['ivar'], R=500.0)
        reb[n] = (wb / (1.0 + Z), fb, eb)
    grid = reb['eboss2019'][0]
    sel = (grid > 2450) & (grid < 6200)
    wl = grid[sel]
    f19 = reb['eboss2019'][1][sel]
    e19 = reb['eboss2019'][2][sel]
    f00 = np.interp(wl, reb['sdss2000'][0], reb['sdss2000'][1])
    good = np.ones_like(wl, dtype=bool)
    for w1, w2 in LINE_WINDOWS:
        good &= ~((wl >= w1) & (wl <= w2))

    # ---- fit free-parameter small-grain reddening of brightened 2000 ----
    def model(p):
        B, av, c2, c4 = p
        xi = np.array([TPD.smallgrain_extinction_shape(w, 27., c2, 5.5, c4)
                       for w in wl])
        return B * f00 * 10.0 ** (-0.4 * av * xi)

    def resid(p):
        return ((model(p) - f19) / e19)[good]

    p0 = [1.5, 0.5, 4.0, 0.0]
    sol = least_squares(resid, p0, bounds=([0.1, -2, 1, 0],
                                           [10, 5, 12, 2]))
    B, av, c2, c4 = sol.x
    chi2 = np.sum(resid(sol.x) ** 2) / (good.sum() - 4)
    print(f"free-law reddening fit: chi2/dof = {chi2:.1f}  "
          f"(B={B:.2f}, A_V={av:+.2f}, c2={c2:.2f}, c4={c4:.2f})")

    # ---- model-independent required attenuation ----
    ratio = f19 / f00
    dm_req = -2.5 * np.log10(ratio)               # = -2.5log B + A_V xi
    xi_best = np.array([TPD.smallgrain_extinction_shape(w, 27., c2, 5.5, c4)
                        for w in wl])
    dm_model = -2.5 * np.log10(B) + av * xi_best
    imin = wl[good][np.argmin(dm_req[good])]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=True,
                                   gridspec_kw={'height_ratios': [1.4, 1],
                                                'hspace': 0.06})
    ax1.fill_between(wl, f19 - e19, f19 + e19, step='mid', color='black',
                     alpha=0.15, lw=0)
    ax1.plot(wl, f19, drawstyle='steps-mid', color='black', lw=1.5,
             alpha=0.85, label=r'$\rm SDSS\mbox{-}2019~(data)$')
    ax1.plot(wl, model(sol.x), '-', color='royalblue', lw=3, alpha=0.95,
             label=rf'$\rm best~free\mbox{{-}}law~reddening~'
                   rf'(\chi^2_\nu={chi2:.0f})$')
    ax1.set_ylabel(r'$f_\lambda~[\rm 10^{-17}~cgs]$', fontsize=17)
    ax1.legend(fontsize=13, frameon=False)
    ax1.set_ylim(0, None)

    ax2.plot(wl[good], dm_req[good], 'o', color='darkgreen', ms=3.5,
             label=r'$\rm required:~-2.5\log(f_{2019}/f_{2000})$')
    ax2.plot(wl, dm_model, '-', color='royalblue', lw=2.5,
             label=r'$\rm any~dust~law:~{\rm const}+A_V\,\xi(\lambda)$')
    ax2.axvline(imin, color='crimson', ls='--', lw=1.5)
    ax2.text(imin + 60, ax2.get_ylim()[1], r'$\rm transparency~window$',
             color='crimson', fontsize=13, rotation=90, va='top')
    ax2.invert_yaxis()
    ax2.set_ylabel(r'$\Delta m~[\rm mag]$', fontsize=17)
    ax2.set_xlabel(r'$\rm rest\mbox{-}frame~wavelength~[\AA]$', fontsize=18)
    ax2.legend(fontsize=12, frameon=False, loc='lower right')

    for ax in (ax1, ax2):
        ax.axvline(3646, color='gray', ls=':', lw=1.2)
        for w1, w2 in LINE_WINDOWS:
            ax.axvspan(w1, w2, color='gray', alpha=0.12, lw=0)
        ax.tick_params(which='major', direction='in', length=8, width=1.5,
                       top=True, right=True, labelsize=14)
        ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                       top=True, right=True)
        ax.minorticks_on()
    ax1.set_xlim(2450, 6200)
    out = os.path.join(PLOTDIR, 'transient_redden_freelaw.png')
    os.makedirs(PLOTDIR, exist_ok=True)
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"required attenuation minimum (transparency window) at "
          f"{imin:.0f} A")
    print(f"Saved {out}")


if __name__ == '__main__':
    main()
