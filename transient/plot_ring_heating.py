#!/usr/bin/env python3
"""
plot_ring_heating.py - three-panel comparison of the multi-ring
temperature-change (Planck derivative) flare model and the adopted
Balmer-edge arm model:

  top    : observed spectra (2000 / 2019 / 2021) with the model
           quiescent disk and BOTH flare models overlaid
  middle : 2019 - 2000 difference, data vs both models
  bottom : 2019 / 2000 ratio, data vs both models

Ring heating: T(r) -> T(r) * [1 + a exp(-(r-r0)^2 / 2 sr^2)] on the
quiet disk (best case of fit_ring_heating.py), plus the same raised-arm
occultation as the adopted model.

Run:  python3 transient/plot_ring_heating.py
"""

import copy
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import transient.tune_redside as tr  # noqa: E402
from transient.transient_disk import disks_from_config  # noqa: E402
from transient.make_transient_figure import (  # noqa: E402
    load_epoch, rebin_R, LINE_WINDOWS, PLOTDIR)
from pillar_disk import C, DAY, FNU_TO_MJY  # noqa: E402

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

RING = dict(r0=2.0, sr=1.0, a=0.50)              # best ring-heating case
ANA = dict(te=12000., fbb=0.10, vs=55000., pas=0.0)   # adopted (i = 60)
NR, NPHI = 200, 400


def main():
    cfg = tr.CFG
    z = cfg['observation']['redshift']
    to_flam = 1e-17 / (1.0 + z)
    comp = cfg['computation']
    wrest = np.logspace(np.log10(comp['wobs_min']),
                        np.log10(comp['wobs_max']),
                        comp['nwavelengths']) / (1.0 + z)

    # adopted-model amplitude from the fit machinery
    WL, ddiff, derr, good, occ_def, A, flare_ref = tr.setup()
    shape_fit = np.array([flare_ref._arm_emission_shape(
        w, ANA['te'], ANA['fbb'], ANA['vs'], ANA['pas']) for w in WL])
    E = A * shape_fit
    w_fit = 1.0 / derr[good] ** 2
    num = np.sum((ddiff[good] + occ_def[good]) * E[good] * w_fit)
    den = np.sum(E[good] ** 2 * w_fit)
    s = max(0.0, num / den)

    # quiet disk + occulted (arm raised, unheated) on the full grid
    c0 = copy.deepcopy(cfg)
    c0['transient']['pillar']['pillar_temp'] = 0.0
    flare0, quiet, _ = disks_from_config(c0, nr=NR, nphi=NPHI)
    q_full = quiet.compute_sed(wrest) * to_flam
    occ0_full = flare0.compute_sed(wrest) * to_flam

    # data and normalization
    binned = {}
    for name in ('sdss2000', 'eboss2019', 'desi2021'):
        ep = load_epoch(name)
        if ep is None:
            continue
        wb, fb, eb = rebin_R(ep['wave_obs'], ep['flam'], ep['ivar'], R=500.0)
        binned[name] = (wb / (1.0 + z), fb, eb)
    w1, w2 = cfg['transient']['normalization']['rest_window']
    wr, fr, _ = binned[cfg['transient']['normalization']['ref_epoch']]
    sel_d = (wr >= w1) & (wr <= w2)
    sel_m = (wrest >= w1) & (wrest <= w2)
    scale = np.median(fr[sel_d]) / np.median(q_full[sel_m])
    q_full *= scale
    occ0_full *= scale

    # adopted Balmer-edge flare model on the full grid
    p0 = cfg['transient']['pillar']
    E_ref_full = flare_ref.compute_sed(wrest) * to_flam * scale - occ0_full
    shape_ref_full = np.array([flare_ref._arm_emission_shape(
        w, p0['balmer_te'], p0['balmer_bb_frac'], p0['balmer_vsmear'],
        p0.get('balmer_paschen', 0.12)) for w in wrest])
    A_full = np.median(E_ref_full / shape_ref_full)
    shape_full = np.array([flare_ref._arm_emission_shape(
        w, ANA['te'], ANA['fbb'], ANA['vs'], ANA['pas']) for w in wrest])
    flam_edge = occ0_full + s * A_full * shape_full

    # ring-heating flare model: heated quiet disk + same occultation
    weight, _ = quiet._sed_weights(False)
    W_r = weight.sum(axis=1)
    r_2d, phi_2d = np.meshgrid(quiet.r, quiet.phi, indexing='ij')
    T_r = quiet.get_temperature(r_2d, phi_2d)[:-1, :].mean(axis=1)
    rr = quiet.r[:-1]
    T2 = T_r * (1.0 + RING['a']
                * np.exp(-0.5 * ((rr - RING['r0']) / RING['sr']) ** 2))
    ld_to_cm = C * DAY
    norm = ld_to_cm ** 2 / (quiet.d * ld_to_cm) ** 2 * FNU_TO_MJY
    dF = np.array([np.sum(W_r * (quiet.planck_function(w, T2)
                                 - quiet.planck_function(w, T_r)))
                   for w in wrest]) * norm * to_flam * scale
    flam_ring = occ0_full + dF
    flam_quiet = q_full

    # ---- figure ----
    os.makedirs(PLOTDIR, exist_ok=True)
    fig, axes = plt.subplots(3, 1, figsize=(12, 13), sharex=True,
                             gridspec_kw={'height_ratios': [2.2, 1, 1],
                                          'hspace': 0.06})
    ax1, ax2, ax3 = axes
    colors = {'sdss2000': 'crimson', 'eboss2019': 'black',
              'desi2021': 'royalblue'}
    labels = {'sdss2000': r'$\rm SDSS\mbox{-}2000$',
              'eboss2019': r'$\rm SDSS\mbox{-}2019$',
              'desi2021': r'$\rm DESI\mbox{-}2021$'}
    for name in ('sdss2000', 'eboss2019', 'desi2021'):
        if name in binned:
            wrb, frb, erb = binned[name]
            ax1.fill_between(wrb, frb - erb, frb + erb, step='mid',
                             color=colors[name], alpha=0.22, lw=0)
            ax1.plot(wrb, frb, drawstyle='steps-mid', color=colors[name],
                     lw=1.3, alpha=0.9, label=labels[name])
    ax1.plot(wrest, flam_quiet, '--', color='darkorange', lw=3, alpha=0.9,
             label=r'$\rm model~quiescent~(no~arm)$')
    ax1.plot(wrest, flam_edge, '-', color='seagreen', lw=3, alpha=0.95,
             label=r'$\rm adopted~flare~(Balmer\mbox{-}edge~arm)$')
    ax1.plot(wrest, flam_ring, '-', color='mediumvioletred', lw=3,
             alpha=0.9, label=r'$\rm flare~(ring~heating)$')

    grid_rest = binned['sdss2000'][0]
    f2000 = binned['sdss2000'][1]
    e2000 = binned['sdss2000'][2]
    f2019 = np.interp(grid_rest, binned['eboss2019'][0],
                      binned['eboss2019'][1])
    e2019 = np.interp(grid_rest, binned['eboss2019'][0],
                      binned['eboss2019'][2])
    d_err = np.sqrt(e2019 ** 2 + e2000 ** 2)
    ax2.fill_between(grid_rest, f2019 - f2000 - d_err,
                     f2019 - f2000 + d_err, step='mid', color='darkgreen',
                     alpha=0.2, lw=0)
    ax2.plot(grid_rest, f2019 - f2000, drawstyle='steps-mid',
             color='darkgreen', lw=1.5, alpha=0.9,
             label=r'$\rm 2019-2000~(data)$')
    ax2.plot(wrest, flam_edge - flam_quiet, '-', color='seagreen', lw=3,
             alpha=0.95, label=r'$\rm adopted$')
    ax2.plot(wrest, flam_ring - flam_quiet, '-', color='mediumvioletred',
             lw=3, alpha=0.9, label=r'$\rm ring~heating$')

    r_val = f2019 / f2000
    r_err = r_val * np.sqrt((e2019 / f2019) ** 2 + (e2000 / f2000) ** 2)
    ax3.fill_between(grid_rest, r_val - r_err, r_val + r_err, step='mid',
                     color='darkgreen', alpha=0.2, lw=0)
    ax3.plot(grid_rest, r_val, drawstyle='steps-mid', color='darkgreen',
             lw=1.5, alpha=0.9, label=r'$\rm 2019/2000~(data)$')
    ax3.plot(wrest, flam_edge / flam_quiet, '-', color='seagreen', lw=3,
             alpha=0.95, label=r'$\rm adopted$')
    ax3.plot(wrest, flam_ring / flam_quiet, '-', color='mediumvioletred',
             lw=3, alpha=0.9, label=r'$\rm ring~heating$')

    ax1.set_ylabel(r'$f_\lambda~[\rm 10^{-17}~erg~s^{-1}~cm^{-2}~\AA^{-1}]$',
                   fontsize=18)
    ax1.legend(fontsize=12, frameon=False, ncol=2, loc='upper right')
    ax1.text(0.02, 0.95,
             rf'$\rm ring~heating:~r_0 = {RING["r0"]:.0f}~ld,~'
             rf'\sigma_r = {RING["sr"]:.0f}~ld,~a = {RING["a"]:.1f};~'
             rf'T(r_0) = {np.interp(RING["r0"], rr, T_r)/1e3:.0f}~kK$',
             transform=ax1.transAxes, va='top', fontsize=15)
    ax2.axhline(0, color='gray', lw=1)
    ax2.set_ylabel(r'$\Delta f_\lambda$', fontsize=18)
    ax2.legend(fontsize=12, frameon=False, loc='upper right')
    ax3.axhline(1, color='gray', lw=1)
    ax3.set_ylabel(r'$\rm ratio$', fontsize=18)
    ax3.set_xlabel(r'$\rm rest\mbox{-}frame~wavelength~[\AA]$', fontsize=18)
    ax3.legend(fontsize=12, frameon=False, loc='upper left')

    for ax in axes:
        for w1_, w2_ in LINE_WINDOWS:
            ax.axvspan(w1_, w2_, color='gray', alpha=0.12, lw=0)
        ax.tick_params(which='major', direction='in', length=8, width=1.5,
                       top=True, right=True, labelsize=14)
        ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                       top=True, right=True)
        ax.minorticks_on()
    ax1.set_xlim(2300, 6900)
    ax1.set_ylim(0, None)

    out = os.path.join(PLOTDIR, 'transient_ring_heating.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")

    # ---- the corresponding temperature profiles ----
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.plot(rr, T_r, '--', color='darkorange', lw=3, alpha=0.9,
            label=r'$\rm quiescent~T(r)$')
    ax.plot(rr, T2, '-', color='mediumvioletred', lw=3, alpha=0.9,
            label=r'$\rm ring\mbox{-}heated~T(r)$')
    ax.axvline(RING['r0'], color='gray', ls=':', lw=1.5)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(0.3, 30)
    ax.set_ylim(1e3, 1e5)
    ax.set_xlabel(r'$\rm radius~[light~days]$', fontsize=18)
    ax.set_ylabel(r'$\rm T~[K]$', fontsize=18)
    ax.legend(fontsize=14, frameon=False)
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=14)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()
    out2 = os.path.join(PLOTDIR, 'transient_ring_heating_tprofile.png')
    plt.savefig(out2, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out2}")


if __name__ == '__main__':
    main()
