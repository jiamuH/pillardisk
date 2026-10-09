#!/usr/bin/env python3
"""
fit_components.py - decompose the observed 2019 - 2000 difference spectrum
into physical continuum components with a non-negative least-squares fit:

    Delta f_lambda = a_bf * Balmer_bf(T_e)
                     + sum_k a_k * BB(T_k)
                     - a_occ * occultation_deficit

The Balmer bound-free continuum and blackbodies are positive emission; the
occultation deficit (computed from the disk model) is a negative term that
suppresses the far-UV. The fit shows whether a smooth (line/Fe II-free)
combination of these pieces reproduces the observed bump + UV turnover.

Run:  python3 transient/fit_components.py
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import nnls

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from transient.transient_disk import TransientPillarDisk  # noqa: E402
from transient.make_transient_figure import (  # noqa: E402
    load_epoch, rebin_R, LINE_WINDOWS, PLOTDIR)

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

Z = 0.494
TO_FLAM = 1e-17 / (1.0 + Z)
BB_TEMPS = [6000., 9000., 13000., 20000., 30000.]   # blackbody grid
TE_BF = 1.0e4          # Balmer bound-free electron temperature (K)
VSMEAR = 15000.        # Doppler smearing of the Balmer edge (km/s)
OCC_INC = 70.0         # inclination used for the occultation-deficit shape


def b_lambda(lam, T):
    x = 1.4388e8 / (lam * T)
    return lam ** -5 / (np.exp(np.clip(x, 0, 700)) - 1.0)


def main():
    # ---------- data difference on a rest grid ----------
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
    ddiff = f2019 - f2000
    derr = np.sqrt(e2000 ** 2 + e2019 ** 2)

    # fit only line-free pixels
    good = np.ones_like(wl, dtype=bool)
    for w1, w2 in LINE_WINDOWS:
        good &= ~((wl >= w1) & (wl <= w2))

    # ---------- basis functions on wl ----------
    disk = TransientPillarDisk(
        rin=0.05, rout=100.0, nr=200, nphi=400, h1=2.0, r0=100.0, beta=10.0,
        tv1=650.0, tirrad_tvisc_ratio=0.15, fcol=1.0, hlamp=0.05,
        dmpc=2792.0, cosi=np.cos(np.radians(OCC_INC)), occultation=True,
        occult_mode='opaque', n_ray=150)
    # Balmer bound-free continuum shape (smeared recombination, no thermal part)
    bf = np.array([disk._arm_emission_shape(float(w), TE_BF, 0.0, VSMEAR)
                   for w in wl])

    # occultation deficit shape: a pillar that only occults (no heating)
    disk.add_pillar(2.0, 0.0, height=1.2, sigma_r=1.5, sigma_phi=0.9,
                    spiral_shear=2.0, arm_length=2.5, pillar_temp=0.0,
                    heat_mode='blackbody')
    f_occ = disk.compute_sed(wl) * TO_FLAM
    f_noocc = disk.compute_sed_no_pillars(wl) * TO_FLAM
    occ_deficit = f_occ - f_noocc            # <= 0, strongest in the UV

    # ---------- build design matrix (normalize columns) ----------
    cols, names = [], []
    cols.append(bf / np.max(bf));                names.append('Balmer b-f')
    for T in BB_TEMPS:
        c = b_lambda(wl, T)
        cols.append(c / np.max(c));              names.append(f'BB {T/1e3:.0f} kK')
    cols.append(occ_deficit / np.max(np.abs(occ_deficit)))
    names.append('occultation')
    A = np.vstack(cols).T                        # (nwl, ncomp)

    # weighted NNLS over line-free pixels
    w = 1.0 / derr
    coef, _ = nnls((A[good] * w[good, None]), ddiff[good] * w[good])
    model = A @ coef
    chi2 = np.sum(((model[good] - ddiff[good]) / derr[good]) ** 2)
    ndof = good.sum() - np.sum(coef > 0)
    print(f"reduced chi^2 = {chi2/ndof:.2f}  ({int(good.sum())} pts, "
          f"{int(np.sum(coef>0))} active components)")
    for nm, c in zip(names, coef):
        if c > 0:
            print(f"  {nm:<14} amplitude = {c:.2f}")

    # ---------- plot ----------
    os.makedirs(PLOTDIR, exist_ok=True)
    fig, ax = plt.subplots(figsize=(12, 7))
    ax.fill_between(wl, ddiff - derr, ddiff + derr, step='mid',
                    color='darkgreen', alpha=0.2, lw=0)
    ax.plot(wl, ddiff, drawstyle='steps-mid', color='darkgreen', lw=1.5,
            alpha=0.9, label=r'$\rm 2019-2000~(data)$')
    ax.plot(wl, model, '-', color='crimson', lw=3, alpha=0.95,
            label=r'$\rm NNLS~combined~fit$')
    comp_colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(names)))
    for nm, c, col in zip(names, coef, comp_colors):
        if c > 0:
            j = names.index(nm)
            ax.plot(wl, A[:, j] * c, '--', color=col, lw=1.8, alpha=0.85,
                    label=nm)
    ax.axhline(0, color='gray', lw=1)
    ax.axvline(3646, color='gray', ls=':', lw=1.5)
    ax.text(3646, ax.get_ylim()[1] * 0.9, r'$\rm Balmer~edge$', fontsize=12,
            rotation=90, va='top', ha='right', color='gray')
    for w1, w2 in LINE_WINDOWS:
        ax.axvspan(w1, w2, color='gray', alpha=0.12, lw=0)
    ax.set_xlabel(r'$\rm rest\mbox{-}frame~wavelength~[\AA]$', fontsize=18)
    ax.set_ylabel(r'$\Delta f_\lambda~[\rm 10^{-17}~erg~s^{-1}~cm^{-2}~\AA^{-1}]$',
                  fontsize=17)
    ax.legend(fontsize=12, frameon=False, ncol=2)
    ax.set_xlim(2450, 6000)
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=14)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()
    out = os.path.join(PLOTDIR, 'transient_component_fit.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")


if __name__ == '__main__':
    main()
