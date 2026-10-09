#!/usr/bin/env python3
"""
fit_selfconsistent.py - fit the fully self-consistent disk model to the
observed 2019 - 2000 difference. Everything comes from compute_sed and the
geometry:
  - the disk is a sum of blackbodies from the temperature profile,
  - the near-side arm OCCULTS the inner disk (geometric ray-march),
  - the same arm EMITS hydrogen recombination continuum (heat_mode balmer),
    which is red-truncated at the Balmer edge.
There are NO free additive components: the occultation is a geometric mask
on the disk blackbodies, not a separate curve.

Fitted: inclination i, electron temperature T_e, thermalized fraction f_bb,
and the emission amplitude (linear, solved analytically per grid point).

Run:  python3 transient/fit_selfconsistent.py
"""

import copy
import itertools
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from transient.transient_disk import disks_from_config  # noqa: E402
from transient.make_transient_figure import (  # noqa: E402
    load_epoch, rebin_R, LINE_WINDOWS, PLOTDIR)

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

HERE = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(HERE, 'config_transient.yaml')) as fh:
    CFG = yaml.safe_load(fh)
CFG['transient']['occultation']['n_ray'] = 150
Z = CFG['observation']['redshift']
TO_FLAM = 1e-17 / (1.0 + Z)
NR, NPHI = 200, 400
REF_T = 15000.0   # reference pillar_temp; emission scales linearly with amp

I_GRID = [46., 54., 62., 70., 78.]
TE_GRID = [8000., 11000., 14000., 18000., 22000.]
BB_GRID = [0.0, 0.2, 0.4]


def base_cfg(inc):
    c = copy.deepcopy(CFG)
    c['observation']['inclination'] = inc
    c['transient']['pillar']['phi_pillar'] = 0.0        # near side
    c['transient']['pillar']['heat_mode'] = 'balmer'    # recombination arm
    c['transient']['pillar'].pop('heat_r', None)        # emit over the arm
    c['transient']['pillar'].pop('heat_sigma_r', None)
    c['transient']['occultation']['mode'] = 'opaque'    # geometric occultation
    return c


def main():
    # ---- data difference on a line-free rest grid ----
    b = {}
    for n in ('sdss2000', 'eboss2019'):
        ep = load_epoch(n)
        wb, fb, eb = rebin_R(ep['wave_obs'], ep['flam'], ep['ivar'], R=500.0)
        b[n] = (wb / (1.0 + Z), fb, eb)
    grid = b['sdss2000'][0]
    sel = (grid > 2450) & (grid < 6000)
    WL = grid[sel]
    f2000 = b['sdss2000'][1][sel]
    f2019 = np.interp(WL, b['eboss2019'][0], b['eboss2019'][1])
    e2000 = b['sdss2000'][2][sel]
    e2019 = np.interp(WL, b['eboss2019'][0], b['eboss2019'][2])
    ddiff = f2019 - f2000
    derr = np.sqrt(e2000 ** 2 + e2019 ** 2)
    good = np.ones_like(WL, dtype=bool)
    for w1, w2 in LINE_WINDOWS:
        good &= ~((WL >= w1) & (WL <= w2))
    nw = CFG['transient']['normalization']['rest_window']
    normsel = (WL >= nw[0]) & (WL <= nw[1])

    results = []
    for inc in I_GRID:
        ci = base_cfg(inc)
        c0 = copy.deepcopy(ci)
        c0['transient']['pillar']['pillar_temp'] = 0.0   # occult, no emission
        disk0, quiet0, _ = disks_from_config(c0, nr=NR, nphi=NPHI)
        q = quiet0.compute_sed(WL) * TO_FLAM
        occ0 = disk0.compute_sed(WL) * TO_FLAM
        scale = np.median(f2000[normsel]) / np.median(q[normsel])
        q *= scale
        occ0 *= scale
        occ_def = q - occ0                               # >0, UV-weighted
        for te, bb in itertools.product(TE_GRID, BB_GRID):
            ce = copy.deepcopy(ci)
            ce['transient']['pillar']['pillar_temp'] = REF_T
            ce['transient']['pillar']['balmer_te'] = te
            ce['transient']['pillar']['balmer_bb_frac'] = bb
            diskE, _, _ = disks_from_config(ce, nr=NR, nphi=NPHI)
            Eref = diskE.compute_sed(WL) * TO_FLAM * scale - occ0  # arm emission
            # linear amplitude: model = s*Eref - occ_def  ~  ddiff
            w = 1.0 / derr[good] ** 2
            num = np.sum((ddiff[good] + occ_def[good]) * Eref[good] * w)
            den = np.sum(Eref[good] ** 2 * w)
            s = max(0.0, num / den)
            model = s * Eref - occ_def
            chi2 = np.sum(((model[good] - ddiff[good]) / derr[good]) ** 2)
            results.append((chi2 / (good.sum() - 4), inc, te, bb, s,
                            model.copy(), s * Eref, occ_def.copy()))
        print(f"  i = {inc:.0f} done")

    results.sort(key=lambda r: r[0])
    print(f"\n  {'chi2/dof':>9}  {'i':>4}  {'T_e[kK]':>7}  {'f_bb':>5}  {'amp':>6}")
    for r in results[:8]:
        print(f"  {r[0]:>9.2f}  {r[1]:>4.0f}  {r[2]/1e3:>7.1f}  {r[3]:>5.2f}"
              f"  {r[4]:>6.2f}")
    best = results[0]
    print(f"\n  BEST: i = {best[1]:.0f} deg, T_e = {best[2]/1e3:.1f} kK, "
          f"f_bb = {best[3]:.2f}  (chi2/dof = {best[0]:.2f})")

    # ---- plot the best self-consistent model ----
    _, inc, te, bb, s, model, emis, occ_def = best

    # decompose the arm emission into its recombination (free-bound) and
    # thermalized-blackbody parts: shape = (1-fbb)*recomb + fbb*BB
    ci = base_cfg(inc)
    c0 = copy.deepcopy(ci)
    c0['transient']['pillar']['pillar_temp'] = 0.0
    disk0, quiet0, _ = disks_from_config(c0, nr=NR, nphi=NPHI)
    q = quiet0.compute_sed(WL) * TO_FLAM
    occ0 = disk0.compute_sed(WL) * TO_FLAM
    scale = np.median(f2000[normsel]) / np.median(q[normsel])
    occ0 *= scale

    def emis_at(fbb_val):
        ce = copy.deepcopy(ci)
        ce['transient']['pillar']['pillar_temp'] = REF_T
        ce['transient']['pillar']['balmer_te'] = te
        ce['transient']['pillar']['balmer_bb_frac'] = fbb_val
        d, _, _ = disks_from_config(ce, nr=NR, nphi=NPHI)
        return d.compute_sed(WL) * TO_FLAM * scale - occ0

    recomb_comp = s * (1.0 - bb) * emis_at(0.0)     # free-bound continuum
    bb_comp = s * bb * emis_at(1.0)                  # thermalized blackbody

    fig, ax = plt.subplots(figsize=(12, 7))
    ax.fill_between(WL, ddiff - derr, ddiff + derr, step='mid',
                    color='darkgreen', alpha=0.2, lw=0)
    ax.plot(WL, ddiff, drawstyle='steps-mid', color='darkgreen', lw=1.5,
            alpha=0.9, label=r'$\rm 2019-2000~(data)$')
    ax.plot(WL, model, '-', color='crimson', lw=3, alpha=0.95,
            label=r'$\rm self\mbox{-}consistent~model$')
    ax.plot(WL, recomb_comp, '--', color='purple', lw=1.8, alpha=0.85,
            label=r'$\rm recombination~(free\mbox{-}bound)$')
    ax.plot(WL, bb_comp, '--', color='darkorange', lw=1.8, alpha=0.85,
            label=rf'$\rm blackbody~(f_{{bb}}={bb:.1f},~T_e)$')
    ax.plot(WL, -occ_def, ':', color='olive', lw=1.8, alpha=0.85,
            label=r'$\rm occultation~deficit~(from~geometry)$')
    ax.axhline(0, color='gray', lw=1)
    ax.axvline(3646, color='gray', ls=':', lw=1.5)
    for w1, w2 in LINE_WINDOWS:
        ax.axvspan(w1, w2, color='gray', alpha=0.12, lw=0)
    ax.text(0.03, 0.95,
            rf'$i={inc:.0f}^\circ,~T_e={te/1e3:.0f}~{{\rm kK}},~'
            rf'f_{{\rm bb}}={bb:.1f}$',
            transform=ax.transAxes, va='top', fontsize=15)
    ax.set_xlabel(r'$\rm rest\mbox{-}frame~wavelength~[\AA]$', fontsize=18)
    ax.set_ylabel(r'$\Delta f_\lambda~[\rm 10^{-17}~erg~s^{-1}~cm^{-2}~\AA^{-1}]$',
                  fontsize=17)
    ax.legend(fontsize=13, frameon=False)
    ax.set_xlim(2450, 6000)
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=14)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()
    out = os.path.join(PLOTDIR, 'transient_selfconsistent_fit.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")


if __name__ == '__main__':
    main()
