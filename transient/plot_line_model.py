#!/usr/bin/env python3
"""
plot_line_model.py - three-panel Hbeta line-profile figure:

  top    : observed continuum-subtracted Hbeta profiles (2000, 2019,
           2021) with the model quiescent line (Keplerian annulus) and
           the model flare line (quiescent + responding annulus)
  middle : 2019 - 2000 difference, data vs model
  bottom : 2019 / 2000 ratio, data vs model

Model: gas on circular Keplerian orbits at inclination i (config), in an
annulus with log-normal radial emissivity. A ring at projected orbital
speed v_p has the classic double-horned profile dP/dv = 1/(pi
sqrt(v_p^2 - v^2)); the annulus integrates over v_p(r) with local
(thermal/turbulent) Gaussian broadening sig_loc. The quiescent line and
the responding component are fitted independently (v_p0, radial width,
sig_loc, amplitude); the flare profile is their sum.

Run:  python3 transient/plot_line_model.py
"""

import itertools
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from transient.measure_line_widths import (  # noqa: E402
    subtracted_profiles, fwhm_of)
from transient.make_transient_figure import PLOTDIR  # noqa: E402

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

C_KMS = 2.99792458e5
LAM_HB = 4861.32
HB_WINS = [(-16000., -11500.), (11500., 17000.)]
MASKS = [(5200., 7000.), (8000., 10500.)]      # [O III] 4959 / 5007
VFIT = 11000.0

VP_GRID = np.arange(400., 4001., 100.)
SLNR_GRID = [0.05, 0.2, 0.4, 0.6, 0.8, 1.0]
SLOC_GRID = [200., 400., 600., 800., 1000., 1300., 1700., 2200., 2800.]

# implied radius bookkeeping (degenerate: r/r_g = (c sin i / v_p)^2)
INCL = 46.0
RG_LD_1E8 = 0.0057      # r_g in light days at M_BH = 1e8 Msun


def annulus_profile(v_grid, vp0, sig_lnr, sig_loc, nr=31):
    """Keplerian annulus line profile on v_grid (unit peak arbitrary):
    log-normal radial emissivity -> v_p = vp0 * exp(-x/2) for ln-radius
    offset x; per-ring bin masses from the arcsin cumulative; then local
    Gaussian broadening."""
    dv = np.median(np.diff(v_grid))
    edges = np.concatenate([v_grid - dv / 2, [v_grid[-1] + dv / 2]])
    xs = np.linspace(-3 * sig_lnr, 3 * sig_lnr, nr)
    wts = np.exp(-0.5 * (xs / sig_lnr) ** 2)
    prof = np.zeros_like(v_grid)
    for x, wt in zip(xs, wts / wts.sum()):
        vp = vp0 * np.exp(-x / 2)
        cum = np.arcsin(np.clip(edges / vp, -1.0, 1.0)) / np.pi
        prof += wt * np.diff(cum)
    nk = int(4 * sig_loc / dv) + 1
    k = np.exp(-0.5 * (np.arange(-nk, nk + 1) * dv / sig_loc) ** 2)
    return np.convolve(prof, k / k.sum(), mode='same')


def fit_annulus(v, target, err, good, narrow=False):
    """Grid fit of the annulus family; with narrow=True a narrow Gaussian
    (the non-varying NLR component) is added, both amplitudes solved
    jointly. Returns (params, broad model, narrow model, chi2)."""
    w = 1.0 / err[good] ** 2
    dv = np.median(np.diff(v))
    sig_n_grid = [200., 300., 400., 500.] if narrow else [None]
    best = None
    for vp0, slnr, sloc, sig_n in itertools.product(
            VP_GRID, SLNR_GRID, SLOC_GRID, sig_n_grid):
        prof = annulus_profile(v, vp0, slnr, sloc)
        if sig_n is None:
            num = np.sum(target[good] * prof[good] * w)
            den = np.sum(prof[good] ** 2 * w)
            s = max(0.0, num / den)
            model_b, model_n = s * prof, np.zeros_like(prof)
            pars = (vp0, slnr, sloc, s, None)
        else:
            ngau = np.exp(-0.5 * (v / sig_n) ** 2)
            ones = np.ones_like(v)
            # joint linear solve for (broad, narrow, constant baseline);
            # the constant absorbs the continuum-subtraction residual and
            # may be negative
            comps = [prof, ngau, ones]
            M = np.array([[np.sum(ci[good] * cj[good] * w)
                           for cj in comps] for ci in comps])
            b = np.array([np.sum(target[good] * ci[good] * w)
                          for ci in comps])
            try:
                a = np.linalg.solve(M, b)
            except np.linalg.LinAlgError:
                continue
            a[:2] = np.clip(a[:2], 0.0, None)
            model_b = a[0] * prof + a[2] * ones
            model_n = a[1] * ngau
            pars = (vp0, slnr, sloc, a[0], sig_n)
        resid = (model_b + model_n - target)[good]
        chi2 = float(np.sum(resid ** 2 * w))
        if best is None or chi2 < best[-1]:
            best = (pars, model_b, model_n, chi2)
    return best


def main():
    v, sub, err = subtracted_profiles(LAM_HB, HB_WINS)
    diff = sub['eboss2019'] - sub['sdss2000']
    derr = np.sqrt(err['eboss2019'] ** 2 + err['sdss2000'] ** 2)
    good = np.abs(v) < VFIT
    for v1, v2 in MASKS:
        good &= ~((v >= v1) & (v <= v2))

    ndof = good.sum() - 5
    (p_q, broad_q, narrow_q, chi2_q) = fit_annulus(
        v, sub['sdss2000'], err['sdss2000'], good, narrow=True)
    model_q = broad_q + narrow_q
    (p_r, model_r, _, chi2_r) = fit_annulus(v, diff, derr, good)
    model_f = model_q + model_r

    sini = np.sin(np.radians(INCL))
    for tag, p, chi2 in (('quiescent', p_q, chi2_q),
                         ('response', p_r, chi2_r)):
        vp0, slnr, sloc, s, sig_n = p
        r_ld = RG_LD_1E8 * (C_KMS * sini / vp0) ** 2
        ntag = (f", narrow sig = {sig_n:.0f} km/s"
                if sig_n is not None else "")
        print(f"{tag:>10}: v_p = {vp0:.0f} km/s, sig_lnr = {slnr:.2f}, "
              f"sig_loc = {sloc:.0f} km/s{ntag}, "
              f"chi2/dof = {chi2/ndof:.2f}; "
              f"r = {r_ld:.0f} ld (M=1e8; r scales with M)")
    fw_q, _ = fwhm_of(v, model_q)
    fw_f, _ = fwhm_of(v, model_f)
    print(f"model FWHM: quiescent = {fw_q:.0f} km/s, "
          f"flare = {fw_f:.0f} km/s")

    # ---- figure ----
    os.makedirs(PLOTDIR, exist_ok=True)
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
             label=r'$\rm model~quiescent~(Keplerian)$')
    ax1.plot(v / 1e3, model_f, '-', color='seagreen', lw=3, alpha=0.95,
             label=r'$\rm model~flare~(+response)$')

    ax2.fill_between(v / 1e3, diff - derr, diff + derr, step='mid',
                     color='darkgreen', alpha=0.2, lw=0)
    ax2.plot(v / 1e3, diff, drawstyle='steps-mid', color='darkgreen',
             lw=1.5, alpha=0.9, label=r'$\rm 2019-2000~(data)$')
    ax2.plot(v / 1e3, model_r, '-', color='seagreen', lw=3, alpha=0.95,
             label=r'$\rm model$')

    base = sub['sdss2000']
    ok = base > 0.15 * np.max(base)
    r_val = np.where(ok, sub['eboss2019'] / np.where(ok, base, 1.0), np.nan)
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
    ax1.text(0.03, 0.95, r'$\rm H\beta~(continuum~subtracted)$',
             transform=ax1.transAxes, va='top', fontsize=16)
    ax2.axhline(0, color='gray', lw=1)
    ax2.set_ylabel(r'$\Delta f_\lambda$', fontsize=17)
    ax2.legend(fontsize=12, frameon=False, loc='upper right')
    ax3.axhline(1, color='gray', lw=1)
    ax3.set_ylabel(r'$\rm ratio$', fontsize=17)
    ax3.set_xlabel(r'$\rm velocity~[10^3~km~s^{-1}]$', fontsize=18)
    ax3.legend(fontsize=12, frameon=False, loc='upper left')
    for ax in axes:
        for v1, v2 in MASKS:
            ax.axvspan(v1 / 1e3, v2 / 1e3, color='gray', alpha=0.12, lw=0)
        ax.axvline(0, color='gray', ls=':', lw=1.5)
        ax.tick_params(which='major', direction='in', length=8, width=1.5,
                       top=True, right=True, labelsize=14)
        ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                       top=True, right=True)
        ax.minorticks_on()
    ax1.set_xlim(-12, 12)
    ax3.set_xticks([-12, -9, -6, -3, 0, 3, 6, 9, 12])

    out = os.path.join(PLOTDIR, 'transient_line_model.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")


if __name__ == '__main__':
    main()
