#!/usr/bin/env python3
"""
plot_cloudy_arm_fit.py - compare the observed 2019-2000 difference
spectrum with (a) the tuned analytic smeared-Balmer-edge model and
(b) the best-fit self-consistent Cloudy arm spectra (continuum + lines
from the same cloud), all through the same geometry (occultation
deficit + linear emission amplitude).

Run:  python3 transient/plot_cloudy_arm_fit.py
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from transient.tune_redside import setup  # noqa: E402
from transient.fit_cloudy_arm import (  # noqa: E402
    smear_loggrid, load_cache, WNORM)
from transient.make_transient_figure import LINE_WINDOWS, PLOTDIR  # noqa: E402

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

# analytic tuned best (tune_redside scan)
ANA = dict(te=8000., fbb=0.025, vs=55000., pas=0.0)
# grid to draw and Cloudy cases: (face, phi, n_H, par3, vsmear)
GRID_NAME = 'arm_column'
CLOUDY_CASES = [
    ('shield', 21.0, 12.0, 22.5, 20000.),  # best total chi2
    ('shield', 18.0, 10.0, 21.5, 15000.),  # best red side
    ('shield', 21.0, 12.0, 23.5, 20000.),  # thickest column, same corner
]
PAR3_LABEL = r'\log N'


def main():
    wave, faces, params = load_cache(GRID_NAME)
    phi, hden = params['phi'], params['n_H']
    p3 = list(params.values())[2]

    WL, ddiff, derr, good, occ_def, A, flare_ref = setup()
    w_fit = 1.0 / derr[good] ** 2

    def solve_amp(E):
        num = np.sum((ddiff[good] + occ_def[good]) * E[good] * w_fit)
        den = np.sum(E[good] ** 2 * w_fit)
        return max(0.0, num / den)

    fig, ax = plt.subplots(figsize=(12, 7))
    ax.fill_between(WL, ddiff - derr, ddiff + derr, step='mid',
                    color='darkgreen', alpha=0.2, lw=0)
    ax.plot(WL, ddiff, drawstyle='steps-mid', color='darkgreen', lw=1.5,
            alpha=0.9, label=r'$\rm 2019-2000~(data)$')

    # analytic tuned best
    shape = np.array([flare_ref._arm_emission_shape(
        w, ANA['te'], ANA['fbb'], ANA['vs'], ANA['pas']) for w in WL])
    E = A * shape
    s = solve_amp(E)
    ax.plot(WL, s * E - occ_def, '-', color='crimson', lw=3, alpha=0.95,
            label=(r'$\rm analytic~edge~(T_e=8~kK,~'
                   r'v=55{,}000~km~s^{-1})$'))

    colors = ['royalblue', 'darkorange', 'purple']
    for (face, p, n, v3, vs), col in zip(CLOUDY_CASES, colors):
        ipt = np.argmin(np.abs(phi - p) + np.abs(hden - n) + np.abs(p3 - v3))
        wg, fg = smear_loggrid(wave, faces[face][ipt] / wave, vs)
        sh = np.interp(WL, wg, fg) / np.interp(WNORM, wg, fg)
        E = A * sh
        s = solve_amp(E)
        lab = (rf'$\rm Cloudy~{face},~\log\varphi={p:.0f},~\log n={n:.0f},~'
               rf'{PAR3_LABEL}={v3:.1f},~v={vs/1e3:.0f}{{,}}000$')
        ax.plot(WL, s * E - occ_def, '-', color=col, lw=2.2, alpha=0.85,
                label=lab)

    ax.axhline(0, color='gray', lw=1)
    ax.axvline(3646, color='gray', ls=':', lw=1.5)
    for w1, w2 in LINE_WINDOWS:
        ax.axvspan(w1, w2, color='gray', alpha=0.12, lw=0)
    ax.set_xlabel(r'$\rm rest\mbox{-}frame~wavelength~[\AA]$', fontsize=18)
    ax.set_ylabel(r'$\Delta f_\lambda~[\rm 10^{-17}~erg~s^{-1}~cm^{-2}~\AA^{-1}]$',
                  fontsize=17)
    ax.legend(fontsize=12, frameon=False)
    ax.set_xlim(2450, 6000)
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=14)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()
    os.makedirs(PLOTDIR, exist_ok=True)
    out = os.path.join(PLOTDIR, 'transient_cloudy_arm_fit.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")


if __name__ == '__main__':
    main()
