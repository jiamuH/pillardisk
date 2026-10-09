#!/usr/bin/env python3
"""
plot_flare_model.py - full-spectrum view of the tuned flare model:

  top    : observed SDSS-2000 / SDSS-2019 / DESI-2021 spectra with the
           model quiescent disk (bare Keplerian bowl, no arm) and the
           model flare state (arm occultation + tuned arm emission)
  middle : 2019 - 2000 difference, data vs model
  bottom : 2019 / 2000 ratio, data vs model

Arm emission: the tuned analytic recombination edge (tune_redside best:
T_e = 8 kK, f_bb = 0.025, no Paschen term, 55,000 km/s Gaussian smear);
amplitude from the same linear solve as the fits.

Run:  python3 transient/plot_flare_model.py
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

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

TE, FBB, VS, PAS = 12000., 0.10, 55000., 0.0     # tuned arm emission (i=60)
NR, NPHI = 200, 400


def main():
    cfg = tr.CFG
    z = cfg['observation']['redshift']
    to_flam = 1e-17 / (1.0 + z)
    comp = cfg['computation']
    wrest = np.logspace(np.log10(comp['wobs_min']),
                        np.log10(comp['wobs_max']),
                        comp['nwavelengths']) / (1.0 + z)

    # ---- amplitude from the same machinery as the fits ----
    WL, ddiff, derr, good, occ_def, A, flare_ref = tr.setup()
    shape_fit = np.array([flare_ref._arm_emission_shape(w, TE, FBB, VS, PAS)
                          for w in WL])
    E = A * shape_fit
    w_fit = 1.0 / derr[good] ** 2
    num = np.sum((ddiff[good] + occ_def[good]) * E[good] * w_fit)
    den = np.sum(E[good] ** 2 * w_fit)
    s = max(0.0, num / den)
    chi2 = float(np.sum(((s * E - occ_def)[good] - ddiff[good]) ** 2
                        * w_fit))
    print(f"amplitude s = {s:.3f} (chi2/dof = {chi2/(good.sum()-5):.1f})")

    # ---- full-wavelength model spectra ----
    c0 = copy.deepcopy(cfg)
    c0['transient']['pillar']['pillar_temp'] = 0.0
    flare0, quiet, _ = disks_from_config(c0, nr=NR, nphi=NPHI)
    q_full = quiet.compute_sed(wrest) * to_flam
    occ0_full = flare0.compute_sed(wrest) * to_flam

    # ---- data ----
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

    # arm-emission prefactor on the full grid (same scalar as A on the
    # fit grid; recomputed here as a cross-check)
    p0 = cfg['transient']['pillar']
    E_ref_full = flare_ref.compute_sed(wrest) * to_flam * scale - occ0_full
    shape_ref_full = np.array([flare_ref._arm_emission_shape(
        w, p0['balmer_te'], p0['balmer_bb_frac'], p0['balmer_vsmear'],
        p0.get('balmer_paschen', 0.12)) for w in wrest])
    A_full = np.median(E_ref_full / shape_ref_full)
    assert abs(A_full / A - 1.0) < 0.02, (A_full, A)

    shape_full = np.array([flare_ref._arm_emission_shape(w, TE, FBB, VS,
                                                         PAS)
                           for w in wrest])
    flam_quiet = q_full
    flam_flare = occ0_full + s * A_full * shape_full

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
    ax1.plot(wrest, flam_flare, '-', color='seagreen', lw=3, alpha=0.95,
             label=r'$\rm model~flare~(occultation+arm)$')

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
    ax2.plot(wrest, flam_flare - flam_quiet, '-', color='seagreen', lw=3,
             alpha=0.95, label=r'$\rm model$')

    r_val = f2019 / f2000
    r_err = r_val * np.sqrt((e2019 / f2019) ** 2 + (e2000 / f2000) ** 2)
    ax3.fill_between(grid_rest, r_val - r_err, r_val + r_err, step='mid',
                     color='darkgreen', alpha=0.2, lw=0)
    ax3.plot(grid_rest, r_val, drawstyle='steps-mid', color='darkgreen',
             lw=1.5, alpha=0.9, label=r'$\rm 2019/2000~(data)$')
    ax3.plot(wrest, flam_flare / flam_quiet, '-', color='seagreen', lw=3,
             alpha=0.95, label=r'$\rm model$')

    ax1.set_ylabel(r'$f_\lambda~[\rm 10^{-17}~erg~s^{-1}~cm^{-2}~\AA^{-1}]$',
                   fontsize=18)
    ax1.legend(fontsize=12, frameon=False, ncol=2, loc='upper right')
    ax1.text(0.02, 0.95,
             rf'$\rm i = {cfg["observation"]["inclination"]:.0f}^\circ,~'
             rf'T_e = {TE/1e3:.0f}~kK,~f_{{bb}} = {FBB:.3f},~'
             rf'v_{{smear}} = {VS/1e3:.0f}{{,}}000~km~s^{{-1}}$',
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

    out = os.path.join(PLOTDIR, 'transient_flare_model.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")


if __name__ == '__main__':
    main()
