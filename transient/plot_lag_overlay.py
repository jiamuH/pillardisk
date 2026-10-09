"""Overlay measured ZTF inter-band lags on the model lag spectrum.

Model curves (flare disk with the arm, quiet smooth bowl) come from
transient/data/model_lag_spectrum.npz (compute_model_lags.py); measured
points are the PyROA lag posteriors (data/pyroa_lags.npz) and the
detrended full-baseline ICCF g-r centroid (data/iccf_cache_gr_detrended.npz).
Everything is observed-frame lag relative to the g band.

Run: python3 transient/plot_lag_overlay.py
"""

import os

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

HERE = os.path.dirname(os.path.abspath(__file__))
DATADIR = os.path.join(HERE, 'data')
PLOTDIR = os.path.join(HERE, 'plots')

# ZTF effective wavelengths, observed frame [Angstrom]
WEFF = {'g': 4783.0, 'r': 6417.0, 'i': 7867.0}


def pct(dist):
    lo, med, hi = np.percentile(dist, [15.87, 50.0, 84.13])
    return med, med - lo, hi - med


def main():
    os.makedirs(PLOTDIR, exist_ok=True)
    model = np.load(os.path.join(DATADIR, 'model_lag_spectrum.npz'))
    pyroa = np.load(os.path.join(DATADIR, 'pyroa_lags.npz'))
    iccf = np.load(os.path.join(DATADIR, 'iccf_cache_gr_detrended.npz'))

    fig, ax = plt.subplots(figsize=(9, 6.5))

    wobs = model['wobs']
    for name, color, label in [
            ('flare', 'orangered', r'$\rm model,~flare~disk~(arm)$'),
            ('quiet', 'royalblue', r'$\rm model,~quiet~disk$')]:
        tau = model[f'tau_{name}_obs']
        dtau = tau - np.interp(WEFF['g'], wobs, tau)
        ax.plot(wobs, dtau, lw=3, alpha=0.8, color=color, label=label)

    for band, marker in [('r', 'o'), ('i', 's')]:
        med, elo, ehi = pct(pyroa[f'tau_{band}'])
        ax.errorbar(WEFF[band], med, yerr=[[elo], [ehi]], fmt=marker,
                    ms=10, color='black', capsize=4, lw=2,
                    label=rf'$\rm PyROA~{band}$' if band == 'r' else
                          rf'$\rm PyROA~{band}$')
    med, elo, ehi = pct(iccf['full_cents'])
    ax.errorbar(WEFF['r'] + 120.0, med, yerr=[[elo], [ehi]], fmt='D', ms=9,
                markerfacecolor='none', color='gray', capsize=4, lw=2,
                label=r'$\rm ICCF~g\mbox{-}r~(detrended)$')

    ax.axhline(0.0, color='gray', lw=1, ls=':')
    ax.axvline(WEFF['g'], color='gray', lw=1, ls=':')
    ax.set_xlabel(r'$\rm Observed~wavelength~[\AA]$', fontsize=18)
    ax.set_ylabel(r'$\rm Lag~relative~to~g~[days]$', fontsize=18)
    ax.legend(fontsize=13, frameon=False, loc='upper left')
    ax.minorticks_on()
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=14)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)

    out = os.path.join(PLOTDIR, 'transient_lag_overlay.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'Saved {out}')


if __name__ == '__main__':
    main()
