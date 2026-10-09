#!/usr/bin/env python3
"""
make_transient_figure.py - build the transient-pillar model, overlay the real
SDSS/eBOSS/DESI spectra, and render a three-panel figure in the style of
Liu et al. Figure 1 (spectra / difference / ratio vs rest wavelength).

Run:  python3 transient/make_transient_figure.py [config_transient.yaml]
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from transient.transient_disk import TransientPillarDisk  # noqa: E402

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

HERE = os.path.dirname(os.path.abspath(__file__))
DATADIR = os.path.join(HERE, 'data')
PLOTDIR = os.path.join(HERE, 'plots')
C_A = 2.998e18  # speed of light in Angstrom/s

# Rest-frame emission-line windows (masked in normalization, shaded in panels)
LINE_WINDOWS = [(2740., 2860.), (3717., 3737.), (3860., 3880.),
                (4090., 4110.), (4320., 4360.), (4810., 4910.),
                (4949., 5017.), (6480., 6650.)]


def load_epoch(name):
    path = os.path.join(DATADIR, f'epoch_{name}.npz')
    if not os.path.exists(path):
        print(f"  epoch '{name}' not found ({path}); skipping")
        return None
    d = np.load(path)
    return dict(wave_obs=d['wave_obs'], flam=d['flam'], ivar=d['ivar'],
                mjd=float(d['mjd']), label=str(d['label']))


def rebin_R(wave, flam, ivar, R=500.0):
    """Inverse-variance-weighted rebin onto a log-wavelength grid of
    resolution R. Returns bin-center wavelengths, binned flux, and the
    1-sigma error on the binned flux (1/sqrt(sum ivar) per bin)."""
    lnw = np.log(wave)
    edges = np.arange(lnw.min(), lnw.max(), 1.0 / R)
    idx = np.digitize(lnw, edges)
    w = np.where(ivar > 0, ivar, 0.0)
    nb = len(edges) + 1
    sum_w = np.bincount(idx, weights=w, minlength=nb)
    sum_wf = np.bincount(idx, weights=w * flam, minlength=nb)
    sum_wl = np.bincount(idx, weights=w * wave, minlength=nb)
    good = sum_w > 0
    return (sum_wl[good] / sum_w[good], sum_wf[good] / sum_w[good],
            1.0 / np.sqrt(sum_w[good]))


def line_free(rest_wave):
    mask = np.ones_like(rest_wave, dtype=bool)
    for w1, w2 in LINE_WINDOWS:
        mask &= ~((rest_wave >= w1) & (rest_wave <= w2))
    return mask


def main():
    cfg_path = sys.argv[1] if len(sys.argv) > 1 else \
        os.path.join(HERE, 'config_transient.yaml')
    with open(cfg_path) as f:
        cfg = yaml.safe_load(f)

    z = cfg['observation']['redshift']
    from transient.transient_disk import disks_from_config
    disk, quiet_disk, dmpc = disks_from_config(cfg)
    p = cfg['transient']['pillar']

    comp = cfg['computation']
    wobs = np.logspace(np.log10(comp['wobs_min']), np.log10(comp['wobs_max']),
                       comp['nwavelengths'])
    wrest = wobs / (1.0 + z)

    print(f"d_L = {dmpc:.0f} Mpc, i = {cfg['observation']['inclination']} deg, "
          f"h_p*tan(i) = {p['height'] * np.tan(np.radians(cfg['observation']['inclination'])):.1f} ld "
          f"vs r_p = {p['r_pillar']} ld")
    print("Computing flare SED (pillar + occultation)...")
    fnu_flare = disk.compute_sed(wrest)
    print("Computing quiescent SED (no pillar)...")
    fnu_quiet = quiet_disk.compute_sed(wrest)

    # Parent compute_sed actually returns rest-frame f_lambda (its Planck
    # function is B_lambda per cm despite the B_nu docstring) times 1e26.
    # Convert to observed-frame f_lambda in 1e-17 erg/s/cm2/A:
    #   f_lam,obs(lam_obs) = f_lam,rest(lam_rest) / (1+z)  at d = d_L
    to_flam = 1e-17 / (1.0 + z)
    flam_flare = fnu_flare * to_flam
    flam_quiet = fnu_quiet * to_flam

    # inclination family of transient SEDs: the adopted model (from the base
    # config) at each viewing angle -- only the inclination is varied; the
    # occultation, emission mode and amplitude all come from the base config.
    inc_curves = {}
    for inc in cfg['transient'].get('inclination_curves', []):
        import copy as _copy
        cfg_i = _copy.deepcopy(cfg)
        cfg_i['observation']['inclination'] = float(inc)
        from transient.transient_disk import disks_from_config as _dfc
        flare_i, quiet_i, _ = _dfc(cfg_i)
        print(f"Computing inclination curve i = {inc:.0f} deg...")
        inc_curves[float(inc)] = (flare_i.compute_sed(wrest) * to_flam,
                                  quiet_i.compute_sed(wrest) * to_flam)

    # occultation diagnostic
    occ_frac = disk.compute_occulted_fraction(np.array([2000., 2500., 3000., 5000.]))
    print("Occulted flux fraction at rest 2000/2500/3000/5000 A: "
          + "/".join(f"{f:.2f}" for f in occ_frac))
    dpk = (flam_flare - flam_quiet)
    print(f"Model bump peak (difference spectrum): rest "
          f"{wrest[np.argmax(dpk)]:.0f} A")

    # ---------------- data ----------------
    epochs = {name: load_epoch(name)
              for name in ('sdss2000', 'eboss2019', 'desi2021')}
    grid_rest = None
    binned = {}
    for name, ep in epochs.items():
        if ep is None:
            continue
        wb, fb, eb = rebin_R(ep['wave_obs'], ep['flam'], ep['ivar'], R=500.0)
        binned[name] = (wb / (1.0 + z), fb, eb)
    if binned:
        # common rest grid for difference/ratio = grid of the reference epoch
        ref = cfg['transient']['normalization']['ref_epoch']
        grid_rest = binned.get(ref, list(binned.values())[0])[0]

    # normalization of the model to the reference epoch in a line-free window
    scale = 1.0
    ref = cfg['transient']['normalization']['ref_epoch']
    w1, w2 = cfg['transient']['normalization']['rest_window']
    if ref in binned:
        wr, fr, _er = binned[ref]
        sel_d = (wr >= w1) & (wr <= w2)
        sel_m = (wrest >= w1) & (wrest <= w2)
        scale = np.median(fr[sel_d]) / np.median(flam_quiet[sel_m])
        print(f"Model flux scale factor (match {ref} at rest "
              f"{w1:.0f}-{w2:.0f} A): {scale:.3f}")
    flam_flare *= scale
    flam_quiet *= scale
    for inc, (ff, fq) in list(inc_curves.items()):
        s_i = np.median(fr[sel_d]) / np.median(fq[sel_m]) \
            if ref in binned else 1.0
        inc_curves[inc] = (ff * s_i, fq * s_i)

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
            wr, fr, er = binned[name]
            ax1.fill_between(wr, fr - er, fr + er, step='mid',
                             color=colors[name], alpha=0.22, lw=0)
            ax1.plot(wr, fr, drawstyle='steps-mid', color=colors[name],
                     lw=1.3, alpha=0.9, label=labels[name])
    ax1.plot(wrest, flam_quiet, '--', color='darkorange', lw=3, alpha=0.9,
             label=r'$\rm model~quiescent$')
    show_balmer = cfg['transient'].get('plot_balmer_model', True)
    if show_balmer:
        ax1.plot(wrest, flam_flare, '-', color='darkorange', lw=3,
                 alpha=0.9, label=r'$\rm model~transient~(Balmer)$')

    # inclination family (geometric occultation), encoded by a colorbar
    import matplotlib as mpl
    inc_cmap = mpl.cm.plasma
    inc_norm = mpl.colors.Normalize(vmin=0.0, vmax=90.0)
    inc_vals = sorted(inc_curves)
    for inc in inc_vals:
        ff, fq = inc_curves[inc]
        ax1.plot(wrest, ff, '-', color=inc_cmap(inc_norm(inc)), lw=2.2,
                 alpha=0.95)
    ax1.set_ylabel(r'$f_\lambda~[\rm 10^{-17}~erg~s^{-1}~cm^{-2}~\AA^{-1}]$',
                   fontsize=18)
    ax1.legend(fontsize=13, frameon=False, ncol=2, loc='upper right')

    # difference and ratio panels (data: 2019 - 2000 on the common grid)
    f2000 = f2019 = None
    if 'sdss2000' in binned and 'eboss2019' in binned:
        f2000 = np.interp(grid_rest, binned['sdss2000'][0], binned['sdss2000'][1])
        f2019 = np.interp(grid_rest, binned['eboss2019'][0], binned['eboss2019'][1])
        e2000 = np.interp(grid_rest, binned['sdss2000'][0], binned['sdss2000'][2])
        e2019 = np.interp(grid_rest, binned['eboss2019'][0], binned['eboss2019'][2])
        d_err = np.sqrt(e2019 ** 2 + e2000 ** 2)
        r_val = f2019 / f2000
        r_err = r_val * np.sqrt((e2019 / f2019) ** 2 + (e2000 / f2000) ** 2)
        ax2.fill_between(grid_rest, f2019 - f2000 - d_err, f2019 - f2000 + d_err,
                         step='mid', color='darkgreen', alpha=0.2, lw=0)
        ax2.plot(grid_rest, f2019 - f2000, drawstyle='steps-mid',
                 color='darkgreen', lw=1.5, alpha=0.9,
                 label=r'$\rm 2019-2000~(data)$')
        ax3.fill_between(grid_rest, r_val - r_err, r_val + r_err, step='mid',
                         color='darkgreen', alpha=0.2, lw=0)
        ax3.plot(grid_rest, r_val, drawstyle='steps-mid', color='darkgreen',
                 lw=1.5, alpha=0.9, label=r'$\rm 2019/2000~(data)$')
    if show_balmer:
        ax2.plot(wrest, flam_flare - flam_quiet, '-', color='darkorange',
                 lw=3, alpha=0.9, label=r'$\rm model~(Balmer~opacity)$')
    for inc in inc_vals:
        ff, fq = inc_curves[inc]
        c = inc_cmap(inc_norm(inc))
        ax2.plot(wrest, ff - fq, '-', color=c, lw=2.2, alpha=0.95)
        ax3.plot(wrest, ff / fq, '-', color=c, lw=2.2, alpha=0.95)
    ax2.axhline(0, color='gray', lw=1)
    ax2.set_ylabel(r'$\Delta f_\lambda$', fontsize=18)
    ax2.legend(fontsize=12, frameon=False, loc='upper right')

    if show_balmer:
        ax3.plot(wrest, flam_flare / flam_quiet, '-', color='darkorange',
                 lw=3, alpha=0.9, label=r'$\rm model~(Balmer~opacity)$')
    ax3.axhline(1, color='gray', lw=1)
    ax3.set_ylabel(r'$\rm ratio$', fontsize=18)
    ax3.set_xlabel(r'$\rm rest\mbox{-}frame~wavelength~[\AA]$', fontsize=18)
    ax3.legend(fontsize=12, frameon=False, loc='upper left')

    # single colorbar encoding inclination, spanning the three panels
    sm = mpl.cm.ScalarMappable(cmap=inc_cmap, norm=inc_norm)
    sm.set_array([])
    cb = fig.colorbar(sm, ax=list(axes), fraction=0.03, pad=0.015,
                      aspect=40)
    cb.set_label(r'$\rm viewing~inclination~i~[deg]$', fontsize=16)
    cb.set_ticks([0, 15, 30, 45, 60, 75, 90])
    cb.ax.tick_params(direction='in', labelsize=13)

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

    out = os.path.join(PLOTDIR, 'transient_sed_figure.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")

    # ---------------- diagnostics ----------------
    mask = disk.compute_observer_occultation()
    fig, ax = plt.subplots(figsize=(9, 7))
    im = ax.pcolormesh(np.degrees(disk.phi), disk.r, mask, cmap='inferno',
                       shading='auto', vmin=0, vmax=1)
    ax.set_yscale('log')
    ax.set_xlabel(r'$\phi~[\rm deg]$', fontsize=18)
    ax.set_ylabel(r'$r~[\rm light~days]$', fontsize=18)
    cb = plt.colorbar(im, ax=ax, shrink=0.8, aspect=20)
    cb.set_label(r'$\rm visibility$', fontsize=16)
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=14)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()
    out2 = os.path.join(PLOTDIR, 'transient_occultation_map.png')
    plt.savefig(out2, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out2}")


if __name__ == '__main__':
    main()
