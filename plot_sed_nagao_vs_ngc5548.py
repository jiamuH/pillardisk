"""Referee reply, Section 4 bullet 5: compare the Nagao (2006) SED adopted in
the paper with the NGC5548.sed table shipped with Cloudy.

(1) Plots the incident continuum of the two single-zone runs at
    log Phi_H = 19, log n_H = 10, solar metallicity. Both runs are normalised
    to the same hydrogen-ionizing photon flux by phi(H) 19.
(2) Prints log(line / H-beta) for both SEDs at log Phi_H = 19 (illuminated)
    and 17.5 (shadow), and the illuminated-to-shadow change.

Run:  python3 plot_sed_nagao_vs_ngc5548.py
"""
import os
import numpy as np
import matplotlib.pyplot as plt

from compare_cloudy_sed_test import D, read_linelist, find_key

PLOT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'plots')

SEDS = [('Nagao', 'nagao', r'$\rm Nagao~et~al.~(2006)$', 'orangered'),
        ('NGC 5548', 'ngc5548', r'$\rm NGC~5548~(Cloudy~table)$', 'royalblue')]
LIT, SHADOW = '19', '17.5'
LINES = [('Ly-alpha', 'h  1 1215.67A'), ('C IV', 'blnd 1549.00A'),
         ('Mg II', 'blnd 2798.00A'), ('H-alpha', 'h  1 6562.80A')]
HBETA = 'h  1 4861.32A'

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy', 'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'


def incident_continuum(prefix):
    """wavelength [A] and incident nu F_nu [erg s^-1 cm^-2] from save continuum."""
    a = np.loadtxt(os.path.join(D, prefix + '_SED.conA'), comments='#', usecols=(0, 1))
    return a[:, 0], a[:, 1]


def plot_sed():
    fig, ax = plt.subplots(figsize=(9, 6.5))
    ymax = 0.
    for _, stem, label, color in SEDS:
        lam, nufnu = incident_continuum('%s_n10_phi%s' % (stem, LIT))
        good = nufnu > 0
        ax.plot(lam[good], nufnu[good], color=color, lw=3, alpha=0.8, label=label)
        ymax = max(ymax, nufnu.max())
    ax.axvline(912., color='gray', ls=':', lw=2)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1., 1e6)
    ax.set_ylim(ymax * 1e-5, ymax * 3)
    ax.set_xlabel(r'$\rm Wavelength~[\mbox{\AA}]$')
    ax.set_ylabel(r'$\lambda F_\lambda~[\rm erg~s^{-1}~cm^{-2}]$')
    ax.minorticks_on()
    ax.tick_params(which='both', direction='in', top=True, right=True)
    ax.tick_params(which='major', length=8, width=1.5)
    ax.tick_params(which='minor', length=4, width=1)
    ax.legend(fontsize=14, frameon=False, loc='lower left')
    fig.tight_layout()
    out = os.path.join(PLOT_DIR, 'sed_nagao_vs_ngc5548.png')
    fig.savefig(out, dpi=200)
    print('saved', out)


def ratio_table():
    print('\nlog(line / H-beta) at log n_H = 10, solar metallicity')
    print('%-10s %12s %12s %12s %12s %12s %12s' % (
        'line', 'Nagao lit', 'Nagao shadow', 'Nagao change',
        '5548 lit', '5548 shadow', '5548 change'))
    for name, key in LINES:
        row = []
        for _, stem, _, _ in SEDS:
            values = []
            for phi in (LIT, SHADOW):
                t = read_linelist('%s_n10_phi%s' % (stem, phi))
                values.append(np.log10(t[find_key(t, key)] / t[find_key(t, HBETA)]))
            row += [values[0], values[1], values[1] - values[0]]
        print('%-10s' % name + ''.join(' %12.2f' % v for v in row))


if __name__ == '__main__':
    plot_sed()
    ratio_table()
