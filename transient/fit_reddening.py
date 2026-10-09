#!/usr/bin/env python3
"""
fit_reddening.py - can dust reddening produce the observed UV break?

Two physically distinct configurations are tested against the 2019 - 2000
difference spectrum (SMC-like extinction, Pei 1992, no 2175 A bump):

  B) FULL SCREEN: dust covers the whole nucleus during the transient;
     flare = (quiescent + arm emission) x 10^(-0.4 A(lambda)).
     Occultation disabled. Free: A_V, T_e, f_bb, emission amplitude.

  D) DUSTY ARM: the same spiral arm, but its sight lines REDDEN the inner
     disk (occult_mode='dust_abs', tau ~ tau_V * SMC(lambda)) instead of
     blocking gray. Free: tau_V, i, T_e, f_bb, emission amplitude.

Both are compared to the adopted gray geometric occultation model.

Run:  python3 transient/fit_reddening.py
"""

import copy
import itertools
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from transient.transient_disk import (  # noqa: E402
    disks_from_config, TransientPillarDisk)
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
REF_T = 15000.0


def get_data():
    b = {}
    for n in ('sdss2000', 'eboss2019'):
        ep = load_epoch(n)
        wb, fb, eb = rebin_R(ep['wave_obs'], ep['flam'], ep['ivar'], R=500.0)
        b[n] = (wb / (1.0 + Z), fb, eb)
    grid = b['sdss2000'][0]
    sel = (grid > 2450) & (grid < 6000)
    wl = grid[sel]
    f2000 = b['sdss2000'][1][sel]
    f2019 = np.interp(wl, b['eboss2019'][0], b['eboss2019'][1])
    e2000 = b['sdss2000'][2][sel]
    e2019 = np.interp(wl, b['eboss2019'][0], b['eboss2019'][2])
    good = np.ones_like(wl, dtype=bool)
    for w1, w2 in LINE_WINDOWS:
        good &= ~((wl >= w1) & (wl <= w2))
    nw = CFG['transient']['normalization']['rest_window']
    normsel = (wl >= nw[0]) & (wl <= nw[1])
    return wl, f2019 - f2000, np.sqrt(e2000**2 + e2019**2), f2000, good, normsel


def fit_amp(target, basis, extra, err, good):
    """weighted LSQ amplitude s>=0 for model = s*basis + extra."""
    w = 1.0 / err[good] ** 2
    num = np.sum((target[good] - extra[good]) * basis[good] * w)
    den = np.sum(basis[good] ** 2 * w)
    return max(0.0, num / den)


def main():
    WL, ddiff, derr, f2000, good, normsel = get_data()
    ndof_base = good.sum()
    xi = np.array([TransientPillarDisk.smc_extinction_shape(w) for w in WL])
    xi_norm = xi / TransientPillarDisk.smc_extinction_shape(5500.0)

    def emission_and_quiet(cfg_mod):
        flare, quiet, _ = disks_from_config(cfg_mod, nr=NR, nphi=NPHI)
        return flare.compute_sed(WL) * TO_FLAM, quiet.compute_sed(WL) * TO_FLAM

    # ---------------- B: full-screen dust ----------------
    bestB = None
    for te, fbb in itertools.product([8000., 11000., 14000., 18000.],
                                     [0.0, 0.2]):
        c = copy.deepcopy(CFG)
        c['observation']['inclination'] = 46.0
        c['transient']['occultation']['enabled'] = False
        c['transient']['pillar'].update(pillar_temp=REF_T, balmer_te=te,
                                        balmer_bb_frac=fbb)
        fl, q = emission_and_quiet(c)
        scale = np.median(f2000[normsel]) / np.median(q[normsel])
        q *= scale
        emis = fl * scale - q                       # pure arm emission
        for av in [0.05, 0.1, 0.15, 0.2, 0.3, 0.45]:
            dust = 10.0 ** (-0.4 * av * xi_norm)
            # diff = dust*(q + s*emis) - q  ->  s*(dust*emis) + (dust-1)*q
            s = fit_amp(ddiff, dust * emis, (dust - 1.0) * q, derr, good)
            model = s * dust * emis + (dust - 1.0) * q
            chi2 = np.sum(((model[good] - ddiff[good]) / derr[good]) ** 2) \
                / (ndof_base - 4)
            if bestB is None or chi2 < bestB[0]:
                bestB = (chi2, te, fbb, av, s, model.copy())
    print(f"B  full-screen dust: chi2/dof = {bestB[0]:.1f}  "
          f"(T_e={bestB[1]/1e3:.0f} kK, f_bb={bestB[2]:.1f}, "
          f"A_V={bestB[3]:.2f})")

    # ---------------- D: dusty arm ----------------
    bestD = None
    DUST = [('smc', None), ('smallgrain', (27., 3., 5.5, 0.)),
            ('smallgrain', (27., 6., 5.5, 0.))]
    for (dlaw, dc), inc, tauv in itertools.product(DUST, [46., 62.],
                                                   [2., 5., 10., 20.]):
        c0 = copy.deepcopy(CFG)
        c0['observation']['inclination'] = inc
        c0['transient']['occultation']['mode'] = 'dust_abs'
        c0['transient']['occultation']['tau_edge'] = tauv
        c0['transient']['occultation']['dust_law'] = dlaw
        if dc is not None:
            c0['transient']['occultation']['dust_c'] = list(dc)
        c0['transient']['pillar']['pillar_temp'] = 0.0
        d0, q0 = emission_and_quiet(c0)
        scale = np.median(f2000[normsel]) / np.median(q0[normsel])
        occ_def = (q0 - d0) * scale                 # reddening deficit
        for te, fbb in itertools.product([11000., 14000.], [0.2]):
            ce = copy.deepcopy(c0)
            ce['transient']['pillar'].update(pillar_temp=REF_T, balmer_te=te,
                                             balmer_bb_frac=fbb)
            fe, _ = emission_and_quiet(ce)
            emis = fe * scale - (d0 * scale)        # arm emission on top
            s = fit_amp(ddiff, emis, -occ_def, derr, good)
            model = s * emis - occ_def
            chi2 = np.sum(((model[good] - ddiff[good]) / derr[good]) ** 2) \
                / (ndof_base - 4)
            if bestD is None or chi2 < bestD[0]:
                bestD = (chi2, inc, tauv, te, fbb, s, model.copy(), dlaw,
                         2 if dc is None else dc[1])
        print(f"  D grid: {dlaw} i={inc:.0f}, tau_V={tauv:.0f} done")
    print(f"D  dusty arm: chi2/dof = {bestD[0]:.1f}  "
          f"(law={bestD[7]}, i={bestD[1]:.0f}, tau_V={bestD[2]:.0f}, "
          f"T_e={bestD[3]/1e3:.0f} kK)")

    # ---------------- adopted gray occultation (reference) ----------------
    cg = copy.deepcopy(CFG)
    cg['observation']['inclination'] = 46.0
    cg0 = copy.deepcopy(cg)
    cg0['transient']['pillar']['pillar_temp'] = 0.0
    dg0, qg = emission_and_quiet(cg0)
    scale = np.median(f2000[normsel]) / np.median(qg[normsel])
    occ_def_g = (qg - dg0) * scale
    cge = copy.deepcopy(cg)
    cge['transient']['pillar']['pillar_temp'] = REF_T
    fg, _ = emission_and_quiet(cge)
    emis_g = fg * scale - dg0 * scale
    s = fit_amp(ddiff, emis_g, -occ_def_g, derr, good)
    model_g = s * emis_g - occ_def_g
    chi2g = np.sum(((model_g[good] - ddiff[good]) / derr[good]) ** 2) \
        / (ndof_base - 4)
    print(f"G  gray occultation (adopted): chi2/dof = {chi2g:.1f}")

    # ---------------- A: dust-only, no emission (de-reddening) ----------------
    # transient = pure change in the line-of-sight dust column:
    #   f_2019 = f_2000 x 10^(0.4 dA_V xi),  diff = f_2000 (10^(0.4 dA_V xi)-1)
    # dA_V > 0 means dust CLEARED between 2000 and 2019 (UV brightens).
    bestA = None
    laws = {'SMC': xi_norm,
            'small-grain': np.array([TransientPillarDisk.smallgrain_extinction_shape(
                w, 27., 3., 5.5, 0.) for w in WL])}
    for lname, xic in laws.items():
        for dav in np.linspace(-0.4, 1.4, 46):
            model = f2000 * (10.0 ** (0.4 * dav * xic) - 1.0)
            chi2 = np.sum(((model[good] - ddiff[good]) / derr[good]) ** 2) \
                / (ndof_base - 1)
            if bestA is None or chi2 < bestA[0]:
                bestA = (chi2, lname, dav, model.copy())
    print(f"A  dust-only (de-redden): chi2/dof = {bestA[0]:.1f}  "
          f"(law={bestA[1]}, dA_V={bestA[2]:.2f})")

    # ---------------- plot ----------------
    os.makedirs(PLOTDIR, exist_ok=True)
    fig, ax = plt.subplots(figsize=(12, 7))
    ax.fill_between(WL, ddiff - derr, ddiff + derr, step='mid',
                    color='darkgreen', alpha=0.2, lw=0)
    ax.plot(WL, ddiff, drawstyle='steps-mid', color='darkgreen', lw=1.5,
            alpha=0.9, label=r'$\rm 2019-2000~(data)$')
    ax.plot(WL, model_g, '-', color='crimson', lw=3, alpha=0.95,
            label=rf'$\rm gray~occultation~(\chi^2_\nu={chi2g:.0f})$')
    ax.plot(WL, bestD[6], '-', color='royalblue', lw=2.5, alpha=0.95,
            label=rf'$\rm dusty~arm~SMC~(\chi^2_\nu={bestD[0]:.0f},~'
                  rf'\tau_V={bestD[2]:.0f})$')
    ax.plot(WL, bestB[5], '--', color='saddlebrown', lw=2.5, alpha=0.95,
            label=rf'$\rm full~dust~screen~(\chi^2_\nu={bestB[0]:.0f},~'
                  rf'A_V={bestB[3]:.2f})$')
    ax.plot(WL, bestA[3], ':', color='purple', lw=2.8, alpha=0.95,
            label=rf'$\rm dust~only,~no~emission~(\chi^2_\nu={bestA[0]:.0f},~'
                  rf'\Delta A_V={bestA[2]:.2f})$')
    ax.axhline(0, color='gray', lw=1)
    ax.axvline(3646, color='gray', ls=':', lw=1.5)
    for w1, w2 in LINE_WINDOWS:
        ax.axvspan(w1, w2, color='gray', alpha=0.12, lw=0)
    ax.set_xlabel(r'$\rm rest\mbox{-}frame~wavelength~[\AA]$', fontsize=18)
    ax.set_ylabel(r'$\Delta f_\lambda~[\rm 10^{-17}~erg~s^{-1}~cm^{-2}~\AA^{-1}]$',
                  fontsize=17)
    ax.legend(fontsize=12, frameon=False)
    ax.set_xlim(2450, 6000)
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=14)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()
    out = os.path.join(PLOTDIR, 'transient_reddening_test.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")


if __name__ == '__main__':
    main()
