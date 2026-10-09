#!/usr/bin/env python3
"""
plot_extinction_laws.py - compare the SMC and small-grain (Ma et al. 2026,
Eq. 1) extinction curves over a wide wavelength range, showing that they are
nearly identical across our transient's rest-frame data window (2500-6000 A)
and diverge only in the far-UV (< ~2000 A), which our data does not reach.

Run:  python3 transient/plot_extinction_laws.py
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from transient.transient_disk import TransientPillarDisk as TPD  # noqa: E402
from transient.make_transient_figure import PLOTDIR  # noqa: E402

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

lam = np.logspace(np.log10(1000.), np.log10(8000.), 400)
smc = np.array([TPD.smc_extinction_shape(w) for w in lam])
smc /= TPD.smc_extinction_shape(5500.0)

fig, ax = plt.subplots(figsize=(11, 7))
ax.axvspan(2500, 6000, color='gold', alpha=0.15, lw=0)
ax.text(3700, 0.5, r'$\rm our~data$', color='goldenrod', fontsize=15, ha='center')
ax.plot(lam, smc, '-', color='black', lw=3, label=r'$\rm SMC~(Pei~1992)$')
for c2, col in [(3., 'royalblue'), (6., 'crimson'), (10., 'darkorange')]:
    sg = np.array([TPD.smallgrain_extinction_shape(w, 27., c2, 5.5, 0.) for w in lam])
    ax.plot(lam, sg, '--', color=col, lw=2.2,
            label=rf'$\rm small\mbox{{-}}grain~(c_2={c2:.0f})$')
ax.axvline(3646, color='gray', ls=':', lw=1.5)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$\rm rest\mbox{-}frame~wavelength~[\AA]$', fontsize=18)
ax.set_ylabel(r'$A_\lambda / A_V$', fontsize=18)
ax.set_xlim(1000, 8000)
ax.set_ylim(0.3, 12)
ax.legend(fontsize=14, frameon=False, loc='upper right')
from matplotlib.ticker import ScalarFormatter
for axis in (ax.xaxis, ax.yaxis):
    axis.set_major_formatter(ScalarFormatter())
ax.tick_params(which='major', direction='in', length=8, width=1.5,
               top=True, right=True, labelsize=14)
ax.tick_params(which='minor', direction='in', length=4, width=1.0,
               top=True, right=True)
os.makedirs(PLOTDIR, exist_ok=True)
out = os.path.join(PLOTDIR, 'transient_extinction_laws.png')
plt.savefig(out, dpi=200, bbox_inches='tight')
plt.close()
print(f"Saved {out}")
