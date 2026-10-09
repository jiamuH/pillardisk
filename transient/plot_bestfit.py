#!/usr/bin/env python3
"""
plot_bestfit.py - three-panel figure of the single best-fit first-principles
model (from scan_ipT.py) overlaid on the observed spectra: model spectra,
difference (flare - quiescent) vs 2019 - 2000, and ratio. A single blackbody
peaking at the observed bump peak is shown for reference in all panels.

Best-fit parameters (edit BEST if the scan is rerun):
    inclination, arm azimuth phi, pillar peak temperature.

Run:  python3 transient/plot_bestfit.py [config_transient.yaml]
"""

import copy
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

BEST = dict(inclination=58.0, phi_deg=90.0, pillar_temp=9500.0)


def main():
    cfg_path = sys.argv[1] if len(sys.argv) > 1 else \
        os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     'config_transient.yaml')
    with open(cfg_path) as f:
        cfg = yaml.safe_load(f)
    z = cfg['observation']['redshift']
    to_flam = 1e-17 / (1.0 + z)

    comp = cfg['computation']
    wobs = np.logspace(np.log10(comp['wobs_min']), np.log10(comp['wobs_max']),
                       comp['nwavelengths'])
    wrest = wobs / (1.0 + z)

    cfgb = copy.deepcopy(cfg)
    cfgb['observation']['inclination'] = BEST['inclination']
    cfgb['transient']['pillar']['phi_pillar'] = float(np.radians(BEST['phi_deg']))
    cfgb['transient']['pillar']['pillar_temp'] = BEST['pillar_temp']
    flare, quiet, _ = disks_from_config(cfgb)
    flam_flare = flare.compute_sed(wrest) * to_flam
    flam_quiet = quiet.compute_sed(wrest) * to_flam

    # ---------------- data ----------------
    binned = {}
    for name in ('sdss2000', 'eboss2019', 'desi2021'):
        ep = load_epoch(name)
        if ep is None:
            continue
        wb, fb, eb = rebin_R(ep['wave_obs'], ep['flam'], ep['ivar'], R=500.0)
        binned[name] = (wb / (1.0 + z), fb, eb)
    ref = cfg['transient']['normalization']['ref_epoch']
    grid_rest = binned[ref][0]
    w1, w2 = cfg['transient']['normalization']['rest_window']
    wr, fr, _ = binned[ref]
    sel_d = (wr >= w1) & (wr <= w2)
    sel_m = (wrest >= w1) & (wrest <= w2)
    scale = np.median(fr[sel_d]) / np.median(flam_quiet[sel_m])
    flam_quiet *= scale
    flam_flare *= scale

    # ---------------- figure ----------------
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
             label=r'$\rm model~quiescent$')
    ax1.plot(wrest, flam_flare, '-', color='seagreen', lw=3, alpha=0.95,
             label=r'$\rm best\mbox{-}fit~transient$')

    f2000 = np.interp(grid_rest, binned['sdss2000'][0], binned['sdss2000'][1])
    f2019 = np.interp(grid_rest, binned['eboss2019'][0], binned['eboss2019'][1])
    e2000 = np.interp(grid_rest, binned['sdss2000'][0], binned['sdss2000'][2])
    e2019 = np.interp(grid_rest, binned['eboss2019'][0], binned['eboss2019'][2])
    d_err = np.sqrt(e2019 ** 2 + e2000 ** 2)
    r_val = f2019 / f2000
    r_err = r_val * np.sqrt((e2019 / f2019) ** 2 + (e2000 / f2000) ** 2)
    ax2.fill_between(grid_rest, f2019 - f2000 - d_err, f2019 - f2000 + d_err,
                     step='mid', color='darkgreen', alpha=0.2, lw=0)
    ax2.plot(grid_rest, f2019 - f2000, drawstyle='steps-mid', color='darkgreen',
             lw=1.5, alpha=0.9, label=r'$\rm 2019-2000~(data)$')
    ax2.plot(wrest, flam_flare - flam_quiet, '-', color='seagreen', lw=3,
             alpha=0.95, label=r'$\rm best\mbox{-}fit$')
    ax3.fill_between(grid_rest, r_val - r_err, r_val + r_err, step='mid',
                     color='darkgreen', alpha=0.2, lw=0)
    ax3.plot(grid_rest, r_val, drawstyle='steps-mid', color='darkgreen',
             lw=1.5, alpha=0.9, label=r'$\rm 2019/2000~(data)$')
    ax3.plot(wrest, flam_flare / flam_quiet, '-', color='seagreen', lw=3,
             alpha=0.95, label=r'$\rm best\mbox{-}fit$')

    # single blackbody peaking at the observed bump peak (reference)
    diff_data = f2019 - f2000
    m = (grid_rest > 2900) & (grid_rest < 4200)
    lam_pk = grid_rest[m][np.argmax(diff_data[m])]
    amp_pk = np.max(diff_data[m])
    T_bb = 2.8978e7 / lam_pk

    def b_lambda(lam, T):
        x = 1.4388e8 / (lam * T)
        return lam ** -5 / (np.exp(np.clip(x, 0, 700)) - 1.0)

    bb = b_lambda(wrest, T_bb) / b_lambda(lam_pk, T_bb) * amp_pk
    ax1.plot(wrest, flam_quiet + bb, '--', color='black', lw=2, alpha=0.8,
             label=r'$\rm quiescent + single~BB$')
    ax2.plot(wrest, bb, '--', color='black', lw=2, alpha=0.8,
             label=rf'$\rm single~BB~(T={T_bb:.0f}~K)$')
    ax3.plot(wrest, (flam_quiet + bb) / flam_quiet, '--', color='black',
             lw=2, alpha=0.8, label=r'$\rm single~BB$')

    ax1.set_ylabel(r'$f_\lambda~[\rm 10^{-17}~erg~s^{-1}~cm^{-2}~\AA^{-1}]$',
                   fontsize=18)
    ax1.legend(fontsize=12, frameon=False, ncol=2, loc='upper right')
    ax1.text(0.02, 0.95,
             rf'$\rm best\mbox{{-}}fit:~i = {BEST["inclination"]:.0f}^\circ,~'
             rf'\phi_{{\rm arm}} = {BEST["phi_deg"]:.0f}^\circ,~'
             rf'T_{{\rm pillar}}^{{\rm peak}} = '
             rf'{BEST["pillar_temp"]/1e3:.1f}~{{\rm kK}}$',
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

    out = os.path.join(PLOTDIR, 'transient_sed_bestfit.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")


if __name__ == '__main__':
    main()
