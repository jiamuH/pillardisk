#!/usr/bin/env python3
"""
plot_line_model_selfcon.py - SELF-CONSISTENT three-panel broad-line
figures: both line-profile shapes come from the disk geometry, with one
shared velocity scale (M_BH) per line.

  quiescent line = full axisymmetric Keplerian ring at r_pillar (the arm
                   gas before shearing/heating; ring_velocity_field)
                   + narrow NLR Gaussian + constant baseline
  flare line     = quiescent + heated-arm kernel (arm_velocity_field,
                   with occultation weighting) as the response

Fitted per line: log M_BH (shared by ring and arm), rotation sense,
narrow width, and linear amplitudes; local broadening fixed at 500 km/s.
Panels: data with model quiescent/flare; difference; ratio.

Lines: Hbeta, Hgamma, Mg II. (Halpha has no quiescent baseline: the
SDSS-2000 spectrum ends at rest 6158 A.)

Run:  python3 transient/plot_line_model_selfcon.py
"""

import copy
import itertools
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import transient.tune_redside as tr  # noqa: E402
from transient.transient_disk import disks_from_config  # noqa: E402
from transient.measure_line_widths import (  # noqa: E402
    subtracted_profiles, fwhm_of)
from transient.make_transient_figure import PLOTDIR  # noqa: E402

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

C_KMS = 2.99792458e5
VFIT = 11000.0
SIG_LOCAL = 500.0
MREF = 1e8
LOGM_GRID = np.arange(5.0, 7.01, 0.125)
SIGN_GRID = [200., 300., 400., 500.]           # narrow NLR width
NR, NPHI = 200, 400

# key: (lam0, TeX label, continuum windows [km/s], fit masks [km/s])
LINES = {
    'hbeta': (4861.32, r'$\rm H\beta$',
              [(-16000., -11500.), (11500., 17000.)],
              [(5200., 7000.), (8000., 10500.)]),      # [O III] 4959/5007
    'hgamma': (4340.47, r'$\rm H\gamma$',
               [(-14500., -10500.), (7000., 12000.)],
               [(800., 2400.)]),                       # [O III] 4363
    'mg2': (2798.75, r'$\rm Mg\,II$',
            [(-20000., -12000.), (12000., 20000.)],
            []),
}


def profile_from_field(v_grid, v0, w0, mfac, orient):
    """Histogram of scaled cell velocities + local Gaussian broadening."""
    vc = orient * v0 * mfac
    dv = np.median(np.diff(v_grid))
    edges = np.concatenate([v_grid - dv / 2, [v_grid[-1] + dv / 2]])
    prof, _ = np.histogram(vc, bins=edges, weights=w0)
    nk = int(4 * SIG_LOCAL / dv) + 1
    k = np.exp(-0.5 * (np.arange(-nk, nk + 1) * dv / SIG_LOCAL) ** 2)
    return np.convolve(prof, k / k.sum(), mode='same')


def fit_line(key, lam0, cwins, masks, v_arm, w_arm, v_ring, w_ring):
    v, sub, err = subtracted_profiles(lam0, cwins)
    diff = sub['eboss2019'] - sub['sdss2000']
    derr = np.sqrt(err['eboss2019'] ** 2 + err['sdss2000'] ** 2)
    good = np.abs(v) < VFIT
    for v1, v2 in masks:
        good &= ~((v >= v1) & (v <= v2))
    w_q = 1.0 / err['sdss2000'][good] ** 2
    w_r = 1.0 / derr[good] ** 2

    best = None
    for logm, orient, sig_n in itertools.product(LOGM_GRID, (+1, -1),
                                                 SIGN_GRID):
        mfac = np.sqrt(10.0 ** logm / MREF)
        ring = profile_from_field(v, v_ring, w_ring, mfac, orient)
        arm = profile_from_field(v, v_arm, w_arm, mfac, orient)
        ngau = np.exp(-0.5 * (v / sig_n) ** 2)
        comps = [ring, ngau, np.ones_like(v)]
        M = np.array([[np.sum(ci[good] * cj[good] * w_q) for cj in comps]
                      for ci in comps])
        b = np.array([np.sum(sub['sdss2000'][good] * ci[good] * w_q)
                      for ci in comps])
        try:
            a = np.linalg.solve(M, b)
        except np.linalg.LinAlgError:
            continue
        a[:2] = np.clip(a[:2], 0.0, None)
        model_q = a[0] * ring + a[1] * ngau + a[2]
        chi2_q = float(np.sum((model_q[good] - sub['sdss2000'][good]) ** 2
                              * w_q))
        num = np.sum(diff[good] * arm[good] * w_r)
        den = np.sum(arm[good] ** 2 * w_r)
        s = max(0.0, num / den)
        chi2_r = float(np.sum((s * arm[good] - diff[good]) ** 2 * w_r))
        tot = chi2_q + chi2_r
        if best is None or tot < best[0]:
            best = (tot, chi2_q, chi2_r, logm, orient, sig_n, a.copy(), s,
                    model_q, s * arm)

    (tot, chi2_q, chi2_r, logm, orient, sig_n, a, s, model_q,
     model_r) = best
    ndof_q = good.sum() - 5
    ndof_r = good.sum() - 3
    fw_q, _ = fwhm_of(v, model_q - a[2])
    fw_r, _ = fwhm_of(v, model_r) if s > 0 else (np.nan, 0.0)
    fw_dq, _ = fwhm_of(v, sub['sdss2000'])
    fw_dr, _ = fwhm_of(v, diff)
    print(f"\n{key}: logM = {logm:.3f}, rot = "
          f"{'+' if orient > 0 else '-'}, narrow sig = {sig_n:.0f} km/s")
    print(f"  quiescent chi2/dof = {chi2_q/ndof_q:.2f}, "
          f"FWHM model/data = {fw_q:.0f}/{fw_dq:.0f} km/s")
    print(f"  response  chi2/dof = {chi2_r/ndof_r:.2f}, "
          f"FWHM model/data = {fw_r:.0f}/{fw_dr:.0f} km/s")
    return dict(v=v, sub=sub, err=err, diff=diff, derr=derr, masks=masks,
                model_q=model_q, model_r=model_r, logm=logm,
                orient=orient)


def plot_line(key, label, r):
    v, sub, err = r['v'], r['sub'], r['err']
    model_q, model_r = r['model_q'], r['model_r']
    model_f = model_q + model_r
    fig, axes = plt.subplots(3, 1, figsize=(11, 13), sharex=True,
                             gridspec_kw={'height_ratios': [2.2, 1, 1],
                                          'hspace': 0.06})
    ax1, ax2, ax3 = axes
    colors = {'sdss2000': 'crimson', 'eboss2019': 'black',
              'desi2021': 'royalblue'}
    labels = {'sdss2000': r'$\rm SDSS\mbox{-}2000$',
              'eboss2019': r'$\rm SDSS\mbox{-}2019$',
              'desi2021': r'$\rm DESI\mbox{-}2021$'}
    for name in ('sdss2000', 'eboss2019', 'desi2021'):
        if name not in sub:
            continue
        ax1.fill_between(v / 1e3, sub[name] - err[name],
                         sub[name] + err[name], step='mid',
                         color=colors[name], alpha=0.18, lw=0)
        ax1.plot(v / 1e3, sub[name], drawstyle='steps-mid',
                 color=colors[name], lw=1.5, alpha=0.9, label=labels[name])
    ax1.plot(v / 1e3, model_q, '--', color='darkorange', lw=3, alpha=0.9,
             label=r'$\rm model~quiescent~(Keplerian~ring)$')
    ax1.plot(v / 1e3, model_f, '-', color='seagreen', lw=3, alpha=0.95,
             label=r'$\rm model~flare~(+heated~arm)$')
    ax1.text(0.03, 0.95,
             label + rf'$;~\log M_{{BH}} = {r["logm"]:.2f},~'
             rf'\rm rot~{"+" if r["orient"] > 0 else "-"}$',
             transform=ax1.transAxes, va='top', fontsize=15)

    ax2.fill_between(v / 1e3, r['diff'] - r['derr'],
                     r['diff'] + r['derr'], step='mid', color='darkgreen',
                     alpha=0.2, lw=0)
    ax2.plot(v / 1e3, r['diff'], drawstyle='steps-mid', color='darkgreen',
             lw=1.5, alpha=0.9, label=r'$\rm 2019-2000~(data)$')
    ax2.plot(v / 1e3, model_r, '-', color='seagreen', lw=3, alpha=0.95,
             label=r'$\rm model~(heated~arm)$')

    base = sub['sdss2000']
    ok = base > 0.15 * np.max(base)
    r_val = np.where(ok, sub['eboss2019'] / np.where(ok, base, 1.0),
                     np.nan)
    r_err = np.where(ok, np.sqrt(err['eboss2019'] ** 2
                                 + (r_val * err['sdss2000']) ** 2)
                     / np.where(ok, base, 1.0), np.nan)
    ax3.fill_between(v / 1e3, r_val - r_err, r_val + r_err, step='mid',
                     color='darkgreen', alpha=0.2, lw=0)
    ax3.plot(v / 1e3, r_val, drawstyle='steps-mid', color='darkgreen',
             lw=1.5, alpha=0.9, label=r'$\rm 2019/2000~(data)$')
    mask_m = model_q > 0.15 * np.max(model_q)
    ax3.plot(v[mask_m] / 1e3, model_f[mask_m] / model_q[mask_m], '-',
             color='seagreen', lw=3, alpha=0.95, label=r'$\rm model$')

    ax1.set_ylabel(
        r'$f_\lambda~[\rm 10^{-17}~erg~s^{-1}~cm^{-2}~\AA^{-1}]$',
        fontsize=17)
    ax1.legend(fontsize=12, frameon=False, ncol=2, loc='upper right')
    ax2.axhline(0, color='gray', lw=1)
    ax2.set_ylabel(r'$\Delta f_\lambda$', fontsize=17)
    ax2.legend(fontsize=12, frameon=False, loc='upper right')
    ax3.axhline(1, color='gray', lw=1)
    ax3.set_ylabel(r'$\rm ratio$', fontsize=17)
    ax3.set_xlabel(r'$\rm velocity~[10^3~km~s^{-1}]$', fontsize=18)
    ax3.legend(fontsize=12, frameon=False, loc='upper left')
    for ax in axes:
        for v1, v2 in r['masks']:
            ax.axvspan(v1 / 1e3, v2 / 1e3, color='gray', alpha=0.12, lw=0)
        ax.axvline(0, color='gray', ls=':', lw=1.5)
        ax.tick_params(which='major', direction='in', length=8, width=1.5,
                       top=True, right=True, labelsize=14)
        ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                       top=True, right=True)
        ax.minorticks_on()
    ax1.set_xlim(-12, 12)
    ax3.set_xticks([-12, -9, -6, -3, 0, 3, 6, 9, 12])

    out = os.path.join(PLOTDIR, f'transient_line_model_selfcon_{key}.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved {out}")


def main():
    cfg = copy.deepcopy(tr.CFG)
    flare, quiet, _ = disks_from_config(cfg, nr=NR, nphi=NPHI)
    p = cfg['transient']['pillar']
    v_arm, w_arm = flare.arm_velocity_field(MREF)
    v_ring, w_ring = quiet.ring_velocity_field(MREF, p['r_pillar'],
                                               p['sigma_r'])
    os.makedirs(PLOTDIR, exist_ok=True)
    for key, (lam0, label, cwins, masks) in LINES.items():
        res = fit_line(key, lam0, cwins, masks, v_arm, w_arm,
                       v_ring, w_ring)
        plot_line(key, label, res)


if __name__ == '__main__':
    main()
