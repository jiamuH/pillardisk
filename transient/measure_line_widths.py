#!/usr/bin/env python3
"""
measure_line_widths.py - quantitative line-profile comparison between
epochs for Hbeta and Mg II:

  1. subtract a local linear continuum (fitted in velocity side windows)
     from each epoch's spectrum,
  2. measure the FWHM of the continuum-subtracted line per epoch,
  3. measure the FWHM of the 2019 - 2000 difference profile (the
     RESPONDING gas), and
  4. plot (a) the continuum-subtracted profiles + difference and (b) the
     ratio 2019/2000 of continuum-subtracted profiles vs velocity.

Run:  python3 transient/measure_line_widths.py
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
REBIN_R = 1500.0

# line: (key, lam0, label, continuum side windows [km/s], plot range)
# Hbeta windows avoid He II 4686 (-10,800) and [O III] 4959/5007
# (+6,000/+9,000); Mg II windows sit on the smooth ramp.
LINES = [
    ('hbeta', 4861.32, r'$\rm H\beta$',
     [(-16000., -11500.), (11500., 17000.)], 12000.),
    ('mg2', 2798.75, r'$\rm Mg\,II$',
     [(-20000., -12000.), (12000., 20000.)], 12000.),
]
COLORS = {'sdss2000': 'crimson', 'eboss2019': 'black',
          'desi2021': 'royalblue'}
LABELS = {'sdss2000': r'$\rm SDSS\mbox{-}2000$',
          'eboss2019': r'$\rm SDSS\mbox{-}2019~(flare)$',
          'desi2021': r'$\rm DESI\mbox{-}2021$'}


def fwhm_of(v, prof):
    """FWHM (km/s) from interpolated half-maximum crossings around the
    profile peak; also returns the peak velocity."""
    ipk = np.argmax(prof)
    half = prof[ipk] / 2.0
    # walk outward from the peak to the first half crossings
    lo = ipk
    while lo > 0 and prof[lo] > half:
        lo -= 1
    hi = ipk
    while hi < len(prof) - 1 and prof[hi] > half:
        hi += 1
    if lo == 0 or hi == len(prof) - 1:
        return np.nan, v[ipk]
    vlo = np.interp(half, [prof[lo], prof[lo + 1]], [v[lo], v[lo + 1]])
    vhi = np.interp(half, [prof[hi], prof[hi - 1]], [v[hi], v[hi - 1]])
    return vhi - vlo, v[ipk]


def subtracted_profiles(lam0, cwins, vmax=25000.):
    """Continuum-subtracted line profiles of all epochs on the 2000-epoch
    velocity grid. Returns (v [km/s], sub {epoch: profile},
    err {epoch: 1-sigma error})."""
    binned = {}
    for name in ('sdss2000', 'eboss2019', 'desi2021'):
        ep = load_epoch(name)
        if ep is None:
            continue
        wb, fb, eb = rebin_R(ep['wave_obs'], ep['flam'], ep['ivar'],
                             R=REBIN_R)
        binned[name] = (wb / (1.0 + Z), fb, eb)
    w0 = binned['sdss2000'][0]
    v = (w0 / lam0 - 1.0) * C_KMS
    m_all = np.abs(v) < vmax
    v = v[m_all]
    sub, err = {}, {}
    for name, (wr, fr, er) in binned.items():
        f = np.interp(w0, wr, fr)[m_all]
        csel = np.zeros_like(v, dtype=bool)
        for v1, v2 in cwins:
            csel |= (v >= v1) & (v <= v2)
        coef = np.polyfit(v[csel], f[csel], 1)
        sub[name] = f - np.polyval(coef, v)
        err[name] = np.interp(w0, wr, er)[m_all]
    return v, sub, err


def main():
    os.makedirs(PLOTDIR, exist_ok=True)
    for key, lam0, label, cwins, vplot in LINES:
        v, sub, err = subtracted_profiles(lam0, cwins)
        print(f"\n{key}:  FWHM of continuum-subtracted line")
        for name in ('sdss2000', 'eboss2019', 'desi2021'):
            if name not in sub:
                continue
            fw, vpk = fwhm_of(v, sub[name])
            print(f"  {name:>9}: FWHM = {fw:7.0f} km/s   "
                  f"(peak at {vpk:+5.0f} km/s)")
        diff = sub['eboss2019'] - sub['sdss2000']
        fw, vpk = fwhm_of(v, diff)
        print(f"  {'2019-2000':>9}: FWHM = {fw:7.0f} km/s   "
              f"(peak at {vpk:+5.0f} km/s)   [responding gas]")

        # ---- continuum-subtracted profiles ----
        fig, ax = plt.subplots(figsize=(10, 7))
        for name in ('sdss2000', 'eboss2019', 'desi2021'):
            ax.plot(v / 1e3, sub[name], drawstyle='steps-mid',
                    color=COLORS[name], lw=2, alpha=0.9, label=LABELS[name])
        ax.plot(v / 1e3, diff, drawstyle='steps-mid', color='seagreen',
                lw=2, alpha=0.9, label=r'$\rm 2019-2000$')
        ax.axhline(0, color='gray', lw=1)
        ax.axvline(0, color='gray', ls=':', lw=1.5)
        ax.text(0.03, 0.95, label + r'$\rm ~(continuum~subtracted)$',
                transform=ax.transAxes, va='top', fontsize=17)
        ax.set_xlabel(r'$\rm velocity~[10^3~km~s^{-1}]$', fontsize=18)
        ax.set_ylabel(
            r'$f_\lambda~[\rm 10^{-17}~erg~s^{-1}~cm^{-2}~\AA^{-1}]$',
            fontsize=17)
        ax.set_xlim(-vplot / 1e3, vplot / 1e3)
        ax.set_xticks([-9, -6, -3, 0, 3, 6, 9])
        ax.legend(fontsize=13, frameon=False, loc='upper right')
        ax.tick_params(which='major', direction='in', length=8, width=1.5,
                       top=True, right=True, labelsize=14)
        ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                       top=True, right=True)
        ax.minorticks_on()
        out = os.path.join(PLOTDIR, f'line_profile_sub_{key}.png')
        plt.savefig(out, dpi=200, bbox_inches='tight')
        plt.close()
        print(f"  Saved {out}")

        # ---- ratio of continuum-subtracted profiles ----
        fig, ax = plt.subplots(figsize=(10, 7))
        base = sub['sdss2000']
        ok = base > 0.15 * np.max(base)      # avoid dividing by noise
        ax.plot(v[ok] / 1e3, (sub['eboss2019'] / base)[ok],
                drawstyle='steps-mid', color='black', lw=2, alpha=0.9,
                label=r'$\rm 2019/2000$')
        if 'desi2021' in sub:
            ax.plot(v[ok] / 1e3, (sub['desi2021'] / base)[ok],
                    drawstyle='steps-mid', color='royalblue', lw=2,
                    alpha=0.9, label=r'$\rm 2021/2000$')
        ax.axhline(1, color='gray', lw=1)
        ax.axvline(0, color='gray', ls=':', lw=1.5)
        ax.text(0.03, 0.95, label + r'$\rm ~line~ratio$',
                transform=ax.transAxes, va='top', fontsize=17)
        ax.set_xlabel(r'$\rm velocity~[10^3~km~s^{-1}]$', fontsize=18)
        ax.set_ylabel(r'$\rm ratio~(continuum~subtracted)$', fontsize=17)
        ax.set_xlim(-vplot / 1e3, vplot / 1e3)
        ax.set_xticks([-9, -6, -3, 0, 3, 6, 9])
        ax.legend(fontsize=13, frameon=False, loc='upper right')
        ax.tick_params(which='major', direction='in', length=8, width=1.5,
                       top=True, right=True, labelsize=14)
        ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                       top=True, right=True)
        ax.minorticks_on()
        out = os.path.join(PLOTDIR, f'line_ratio_{key}.png')
        plt.savefig(out, dpi=200, bbox_inches='tight')
        plt.close()
        print(f"  Saved {out}")


if __name__ == '__main__':
    main()
