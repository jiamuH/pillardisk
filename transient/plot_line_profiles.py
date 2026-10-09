#!/usr/bin/env python3
"""
plot_line_profiles.py - emission-line profiles in velocity space for the
three epochs (SDSS 2000, eBOSS 2019 flare, DESI 2021), one figure per
line, with the 2019 - 2000 difference overlaid. Looks for line-profile
changes associated with the flare (arm kinematics).

Lines: Mg II 2798, Hgamma 4340, Hbeta 4861, [O III] 5007 (narrow-line
flux-calibration check).

Run:  python3 transient/plot_line_profiles.py
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from transient.make_transient_figure import (  # noqa: E402
    load_epoch, rebin_R, PLOTDIR)

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

Z = 0.494
C_KMS = 2.99792458e5
REBIN_R = 1500.0          # keep enough resolution for profile shapes
VMAX = 25000.0            # velocity half-range (km/s)

LINES = [
    ('mg2', 2798.75, r'$\rm Mg\,II~\lambda 2799$'),
    ('hgamma', 4340.47, r'$\rm H\gamma~\lambda 4340$'),
    ('hbeta', 4861.32, r'$\rm H\beta~\lambda 4861$'),
    ('oiii', 5006.84, r'$\rm [O\,III]~\lambda 5007$'),
]

COLORS = {'sdss2000': 'crimson', 'eboss2019': 'black',
          'desi2021': 'royalblue'}
LABELS = {'sdss2000': r'$\rm SDSS\mbox{-}2000$',
          'eboss2019': r'$\rm SDSS\mbox{-}2019~(flare)$',
          'desi2021': r'$\rm DESI\mbox{-}2021$'}


def main():
    binned = {}
    for name in ('sdss2000', 'eboss2019', 'desi2021'):
        ep = load_epoch(name)
        if ep is None:
            continue
        wb, fb, eb = rebin_R(ep['wave_obs'], ep['flam'], ep['ivar'],
                             R=REBIN_R)
        binned[name] = (wb / (1.0 + Z), fb, eb)

    os.makedirs(PLOTDIR, exist_ok=True)
    for key, lam0, label in LINES:
        fig, ax = plt.subplots(figsize=(10, 7))
        for name in ('sdss2000', 'eboss2019', 'desi2021'):
            if name not in binned:
                continue
            wr, fr, er = binned[name]
            v = (wr / lam0 - 1.0) * C_KMS
            m = np.abs(v) < VMAX
            ax.fill_between(v[m] / 1e3, fr[m] - er[m], fr[m] + er[m],
                            step='mid', color=COLORS[name], alpha=0.18,
                            lw=0)
            ax.plot(v[m] / 1e3, fr[m], drawstyle='steps-mid',
                    color=COLORS[name], lw=2, alpha=0.9, label=LABELS[name])

        # 2019 - 2000 difference on the 2000 grid
        w0, f0, e0 = binned['sdss2000']
        f19 = np.interp(w0, binned['eboss2019'][0], binned['eboss2019'][1])
        e19 = np.interp(w0, binned['eboss2019'][0], binned['eboss2019'][2])
        v = (w0 / lam0 - 1.0) * C_KMS
        m = np.abs(v) < VMAX
        derr = np.sqrt(e0 ** 2 + e19 ** 2)
        ax.fill_between(v[m] / 1e3, (f19 - f0 - derr)[m],
                        (f19 - f0 + derr)[m], step='mid', color='seagreen',
                        alpha=0.18, lw=0)
        ax.plot(v[m] / 1e3, (f19 - f0)[m], drawstyle='steps-mid',
                color='seagreen', lw=2, alpha=0.9,
                label=r'$\rm 2019-2000$')

        ax.axhline(0, color='gray', lw=1)
        ax.axvline(0, color='gray', ls=':', lw=1.5)
        ax.text(0.03, 0.95, label, transform=ax.transAxes, va='top',
                fontsize=17)
        ax.set_xlabel(r'$\rm velocity~[10^3~km~s^{-1}]$', fontsize=18)
        ax.set_ylabel(
            r'$f_\lambda~[\rm 10^{-17}~erg~s^{-1}~cm^{-2}~\AA^{-1}]$',
            fontsize=17)
        ax.set_xticks([-20, -10, 0, 10, 20])
        ax.legend(fontsize=13, frameon=False, loc='upper right')
        ax.tick_params(which='major', direction='in', length=8, width=1.5,
                       top=True, right=True, labelsize=14)
        ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                       top=True, right=True)
        ax.minorticks_on()
        out = os.path.join(PLOTDIR, f'line_profile_{key}.png')
        plt.savefig(out, dpi=200, bbox_inches='tight')
        plt.close()
        print(f"Saved {out}")


if __name__ == '__main__':
    main()
