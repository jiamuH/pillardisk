#!/usr/bin/env python3
"""
fit_line_budget.py - recombination energy-budget test of the bump.

If the 2019 bump is hydrogen recombination continuum (Balmer continuum),
Case B ties its luminosity to the hydrogen recombination LINES:
   L(Balmer continuum) ~ 8 * L(Hbeta)   (T_e ~ 1e4 K)
so a bump of luminosity L_bump demands the Balmer lines brighten by
~L_bump/8. We measure the Hbeta (and Hgamma) broad-line flux CHANGE from
2000 to 2019 and compare it to what the bump requires. If the lines did not
brighten by that much, the bump cannot be recombination continuum.

Run:  python3 transient/fit_line_budget.py
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
DL_CM = 2792.0 * 3.0857e24          # luminosity distance in cm (z=0.494)
FOURPI_DL2 = 4 * np.pi * DL_CM ** 2
RATIO_BAC_HB = 8.0                  # L(Balmer continuum)/L(Hbeta), Case B 1e4 K


def line_flux(wl, f, e, line, blue, red):
    """Integrated flux above a linear continuum between the blue/red anchor
    windows, in erg/s/cm2 (input f in 1e-17 cgs). Returns (flux, err)."""
    def med(win):
        m = (wl >= win[0]) & (wl <= win[1])
        return np.median(wl[m]), np.median(f[m])
    xb, yb = med(blue)
    xr, yr = med(red)
    m = (wl >= line[0]) & (wl <= line[1])
    cont = yb + (yr - yb) * (wl[m] - xb) / (xr - xb)
    dl = np.gradient(wl)[m]
    flux = np.sum((f[m] - cont) * dl) * 1e-17
    err = np.sqrt(np.sum((e[m] * dl) ** 2)) * 1e-17
    return flux, err


def main():
    reb = {}
    for n in ('sdss2000', 'eboss2019'):
        ep = load_epoch(n)
        wb, fb, eb = rebin_R(ep['wave_obs'], ep['flam'], ep['ivar'], R=1500.0)
        reb[n] = (wb / (1.0 + Z), fb, eb)

    # ---- Hbeta and Hgamma broad-line fluxes (hydrogen recombination lines) ----
    lines = {'Hbeta': (dict(line=(4760, 4950), blue=(4600, 4680),
                            red=(5090, 5160))),
             'Hgamma': (dict(line=(4290, 4400), blue=(4180, 4240),
                             red=(4440, 4500)))}
    print(f"{'line':8} {'F_2000':>10} {'F_2019':>10} {'dF':>10}  [erg/s/cm2]")
    dL_line = {}
    for name, w in lines.items():
        f00, e00 = line_flux(*reb['sdss2000'], **w)
        f19, e19 = line_flux(*reb['eboss2019'], **w)
        dL_line[name] = (f19 - f00) * FOURPI_DL2
        print(f"{name:8} {f00:10.2e} {f19:10.2e} {f19-f00:10.2e}"
              f"   dL={ (f19-f00)*FOURPI_DL2:.2e} erg/s")

    # ---- bump continuum luminosity (integrate the difference continuum) ----
    wl = reb['sdss2000'][0]
    f00 = reb['sdss2000'][1]
    e00 = reb['sdss2000'][2]
    f19 = np.interp(wl, reb['eboss2019'][0], reb['eboss2019'][1])
    e19 = np.interp(wl, reb['eboss2019'][0], reb['eboss2019'][2])
    diff = f19 - f00
    ediff = np.sqrt(e00 ** 2 + e19 ** 2)

    # robust Hbeta CHANGE: measure the line on the DIFFERENCE spectrum, with
    # anchors just outside the broad line (avoids the bump-tilted continuum)
    dF_hb, dF_hb_e = line_flux(wl, diff, ediff, line=(4760, 4930),
                               blue=(4700, 4745), red=(5030, 5075))
    dL_line['Hbeta'] = dF_hb * FOURPI_DL2
    print(f"\nHbeta CHANGE (measured on the difference): "
          f"dF={dF_hb:.2e} +/- {dF_hb_e:.1e} cgs, "
          f"dL={dF_hb*FOURPI_DL2:.2e} erg/s")
    # line-free continuum points across the bump, interpolate across lines
    cont_anchors = np.array([2650., 2900., 3050., 3200., 3400., 3600., 3900.,
                             4000., 4250., 4550., 4700.])
    da = np.interp(cont_anchors, wl, diff)
    # Balmer-continuum region 2700-3646, and the whole bump 2700-5000
    def integ(w1, w2):
        m = (wl > w1) & (wl < w2)
        dc = np.interp(wl[m], cont_anchors, da)   # smooth continuum diff
        return np.trapz(np.clip(dc, 0, None), wl[m]) * 1e-17

    F_bac = integ(2700., 3646.)          # Balmer-continuum region
    F_bump = integ(2700., 5000.)         # whole recombination bump
    L_bac = F_bac * FOURPI_DL2
    L_bump = F_bump * FOURPI_DL2
    print(f"\nBalmer-continuum region (2700-3646 A): "
          f"F={F_bac:.2e} cgs, L={L_bac:.2e} erg/s")
    print(f"whole bump (2700-5000 A):              "
          f"F={F_bump:.2e} cgs, L={L_bump:.2e} erg/s")

    # ---- the test ----
    L_hb_required = L_bac / RATIO_BAC_HB
    L_hb_observed = dL_line['Hbeta']
    print(f"\n--- recombination energy budget ---")
    print(f"Hbeta luminosity REQUIRED if bump = Balmer continuum:")
    print(f"   L(Hbeta) = L(BaC)/{RATIO_BAC_HB:.0f} = {L_hb_required:.2e} erg/s")
    print(f"Hbeta luminosity CHANGE observed 2000->2019:")
    print(f"   dL(Hbeta) = {L_hb_observed:.2e} erg/s")
    ratio = L_hb_observed / L_hb_required if L_hb_required else float('nan')
    print(f"observed / required = {ratio:.2f}")
    if ratio < 0.3:
        print("VERDICT: lines did NOT brighten enough -> bump is NOT "
              "recombination continuum (tension with the adopted model).")
    elif ratio > 0.5:
        print("VERDICT: line brightening CONSISTENT with recombination "
              "continuum (supports the adopted model).")
    else:
        print("VERDICT: marginal.")

    # ---------------- figure: observed vs recombination-expected ----------
    # smooth bump continuum from the difference (interpolate across lines)
    anchors = np.array([2550., 2650., 3050., 3200., 3400., 3600., 3900.,
                        4000., 4250., 4550., 4700., 5200., 5600., 6000.])
    bump_cont = np.interp(wl, anchors, np.interp(anchors, wl, diff))

    # Balmer lines Case B demands if the bump is Balmer continuum:
    # L(Hbeta) = L(BaC)/8, and the decrement fixes the others.
    F_hb = L_hb_required / FOURPI_DL2 / 1e-17            # integrated, 1e-17 units
    caseB = {4861.: 1.0, 4340.: 0.468, 4102.: 0.259, 3970.: 0.159,
             3889.: 0.105}                              # Hb, Hg, Hd, He, H8
    vfwhm = 3000.0                                       # broad-line FWHM km/s
    exp_lines = np.zeros_like(wl)
    for lam0, r in caseB.items():
        sig = vfwhm / 2.355 / 2.998e5 * lam0
        exp_lines += (F_hb * r) / (sig * np.sqrt(2 * np.pi)) \
            * np.exp(-0.5 * ((wl - lam0) / sig) ** 2)
    expected = bump_cont + exp_lines

    expected_2019 = f00 + bump_cont + exp_lines   # 2000 + bump + Case B lines

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 10), sharex=True,
                                   gridspec_kw={'height_ratios': [1.8, 1],
                                                'hspace': 0.05})
    m = (wl > 2500) & (wl < 6050)
    # TOP: raw spectra
    ax1.plot(wl[m], f00[m], drawstyle='steps-mid', color='crimson', lw=1.2,
             alpha=0.7, label=r'$\rm 2000~(quiescent)$')
    ax1.plot(wl[m], f19[m], drawstyle='steps-mid', color='black', lw=1.4,
             alpha=0.85, label=r'$\rm 2019~(observed)$')
    ax1.plot(wl[m], expected_2019[m], '-', color='royalblue', lw=2.3,
             alpha=0.9,
             label=r'$\rm expected~2019~if~bump=Balmer~continuum$')
    for lam0, nm in [(4861, r'$\rm H\beta$'), (4340, r'$\rm H\gamma$'),
                     (4102, r'$\rm H\delta$')]:
        ax1.text(lam0, expected_2019[np.argmin(np.abs(wl - lam0))] + 3, nm,
                 fontsize=13, ha='center', color='royalblue')
    ax1.set_ylabel(r'$f_\lambda~[\rm 10^{-17}~erg~s^{-1}~cm^{-2}~\AA^{-1}]$',
                   fontsize=16)
    ax1.legend(fontsize=13, frameon=False, loc='upper right')
    ax1.set_ylim(0, None)
    # BOTTOM: difference
    ax2.fill_between(wl[m], (diff - ediff)[m], (diff + ediff)[m], step='mid',
                     color='black', alpha=0.15, lw=0)
    ax2.plot(wl[m], diff[m], drawstyle='steps-mid', color='black', lw=1.5,
             alpha=0.9, label=r'$\rm observed~(2019-2000)$')
    ax2.plot(wl[m], expected[m], '-', color='royalblue', lw=2.3, alpha=0.9,
             label=r'$\rm expected~(bump + Case~B~lines)$')
    ax2.plot(wl[m], bump_cont[m], '--', color='seagreen', lw=1.6, alpha=0.8,
             label=r'$\rm bump~continuum~L(BaC)$')
    ax2.axhline(0, color='gray', lw=1)
    ax2.set_ylabel(r'$\Delta f_\lambda$', fontsize=17)
    ax2.set_xlabel(r'$\rm rest\mbox{-}frame~wavelength~[\AA]$', fontsize=18)
    ax2.legend(fontsize=12, frameon=False, loc='upper right')
    for ax in (ax1, ax2):
        ax.axvline(3646, color='gray', ls=':', lw=1.2)
        ax.set_xlim(2500, 6050)
        ax.tick_params(which='major', direction='in', length=8, width=1.5,
                       top=True, right=True, labelsize=14)
        ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                       top=True, right=True)
        ax.minorticks_on()
    out = os.path.join(PLOTDIR, 'transient_line_budget.png')
    os.makedirs(PLOTDIR, exist_ok=True)
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")


if __name__ == '__main__':
    main()
