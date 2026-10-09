"""Plot the cleaned ZTF multi-band light curves of the transient quasar.

Single panel, flux [mJy] vs MJD, with the eBOSS (flare) and DESI
(quiescent) spectral epochs marked.

Run: python3 transient/plot_lightcurves.py
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

BANDS = [('g', 'royalblue'), ('r', 'orangered'), ('i', 'darkgoldenrod')]
EPOCHS = [(58523.0, r'$\rm eBOSS~2019$'), (59518.0, r'$\rm DESI~2021$')]


def load_lc(band):
    path = os.path.join(DATADIR, f'lc_ztf_{band}.npz')
    if not os.path.exists(path):
        print(f'missing {path}, skipped')
        return None
    return np.load(path)


def main():
    os.makedirs(PLOTDIR, exist_ok=True)
    fig, ax = plt.subplots(figsize=(12, 6))

    for band, color in BANDS:
        lc = load_lc(band)
        if lc is None:
            continue
        ax.errorbar(lc['mjd'], lc['flux_mjy'], yerr=lc['fluxerr_mjy'],
                    fmt='o', ms=4, lw=1, elinewidth=1, color=color,
                    alpha=0.8, label=rf'$\rm ZTF~{band}$')

    for mjd, label in EPOCHS:
        ax.axvline(mjd, color='gray', ls='--', lw=1.5)
        ax.text(mjd + 15, 0.985, label, rotation=90, va='top', ha='left',
                fontsize=13, color='gray',
                transform=ax.get_xaxis_transform())

    ax.set_xlabel(r'$\rm MJD~[days]$', fontsize=18)
    ax.set_ylabel(r'$\rm Flux~[mJy]$', fontsize=18)
    ax.legend(fontsize=14, frameon=False, loc='upper right')
    ax.minorticks_on()
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=14)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)

    out = os.path.join(PLOTDIR, 'transient_lightcurves.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'Saved {out}')


if __name__ == '__main__':
    main()
