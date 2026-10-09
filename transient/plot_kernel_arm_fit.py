#!/usr/bin/env python3
"""
plot_kernel_arm_fit.py - show the geometric-kernel arm models against the
observed 2019-2000 difference spectrum:

  - best kernel fit overall (arm azimuth 30 deg, logM = 9.5, edge
    T_e = 8 kK, f_bb = 0.05, rotation '-'),
  - best kernel fit at a PHYSICAL black-hole mass (logM = 9.0, same
    azimuth, Cloudy shield slab phi=21 n=12 N=23.0),
  - the ad hoc Gaussian model (T_e = 8 kK, v = 55,000 km/s) for
    reference.

Also saves a second figure with the kernel shape itself (at 1e8 Msun
and scaled to the best-fit logM = 9.5).

Run:  python3 transient/plot_kernel_arm_fit.py
"""

import copy
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import transient.tune_redside as tr  # noqa: E402
from transient.fit_cloudy_arm import load_cache  # noqa: E402
from transient.fit_kernel_arm import (  # noqa: E402
    kernel_on_grid, WNORM, DLN, MREF, LOCAL_VS)
from transient.make_transient_figure import LINE_WINDOWS, PLOTDIR  # noqa: E402

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

PHI_ARM = 30.0
ANA = dict(te=8000., fbb=0.025, vs=55000., pas=0.0)   # Gaussian reference


def solve_amp(E, ddiff, occ_def, derr, good):
    w = 1.0 / derr[good] ** 2
    num = np.sum((ddiff[good] + occ_def[good]) * E[good] * w)
    den = np.sum(E[good] ** 2 * w)
    return max(0.0, num / den)


def main():
    g = np.arange(np.log(1500.0), np.log(11000.0), DLN)
    wg = np.exp(g)

    # ---- reference geometry (config, phi_arm = 0): Gaussian model ----
    WL, ddiff, derr, good, occ_def0, A0, flare0 = tr.setup()
    shape_g = np.array([flare0._arm_emission_shape(
        w, ANA['te'], ANA['fbb'], ANA['vs'], ANA['pas']) for w in WL])
    E = A0 * shape_g
    s = solve_amp(E, ddiff, occ_def0, derr, good)
    model_gauss = s * E - occ_def0

    # ---- kernel geometry (phi_arm = 30 deg) ----
    cfg = copy.deepcopy(tr.CFG)
    cfg['transient']['pillar']['phi_pillar'] = float(np.radians(PHI_ARM))
    _, _, _, _, occ_def, A, flare_ref = tr.setup(cfg)
    v0, w0 = flare_ref.arm_velocity_field(MREF)

    S_edge = np.array([flare_ref._arm_emission_shape(
        w, 8000., 0.05, LOCAL_VS, 0.0) for w in wg])
    cwave, cfaces, cparams = load_cache('arm_column')
    cp = list(cparams.values())
    ipt = int(np.argmin(np.abs(cp[0] - 21.) + np.abs(cp[1] - 12.)
                        + np.abs(cp[2] - 23.)))
    S_cl = np.interp(wg, cwave, cfaces['shield'][ipt] / cwave)

    cases = [
        (S_edge, 9.5, -1, 'darkorange',
         r'$\rm kernel,~\log M=9.5,~edge~T_e=8~kK$'),
        (S_cl, 9.0, -1, 'royalblue',
         r'$\rm kernel,~\log M=9.0,~Cloudy~slab$'),
    ]

    fig, ax = plt.subplots(figsize=(12, 7))
    ax.fill_between(WL, ddiff - derr, ddiff + derr, step='mid',
                    color='darkgreen', alpha=0.2, lw=0)
    ax.plot(WL, ddiff, drawstyle='steps-mid', color='darkgreen', lw=1.5,
            alpha=0.9, label=r'$\rm 2019-2000~(data)$')
    ax.plot(WL, model_gauss, '-', color='crimson', lw=3, alpha=0.95,
            label=r'$\rm ad~hoc~Gaussian~(v=55{,}000~km~s^{-1})$')
    for S, logm, orient, col, lab in cases:
        kern = kernel_on_grid(v0, w0, np.sqrt(10.0 ** logm / MREF), orient)
        F = fftconvolve(S, kern, mode='same')
        shape = np.interp(WL, wg, F) / np.interp(WNORM, wg, F)
        E = A * shape
        s = solve_amp(E, ddiff, occ_def, derr, good)
        ax.plot(WL, s * E - occ_def, '-', color=col, lw=2.2, alpha=0.9,
                label=lab)

    ax.axhline(0, color='gray', lw=1)
    ax.axvline(3646, color='gray', ls=':', lw=1.5)
    for w1, w2 in LINE_WINDOWS:
        ax.axvspan(w1, w2, color='gray', alpha=0.12, lw=0)
    ax.set_xlabel(r'$\rm rest\mbox{-}frame~wavelength~[\AA]$', fontsize=18)
    ax.set_ylabel(r'$\Delta f_\lambda~[\rm 10^{-17}~erg~s^{-1}~cm^{-2}~\AA^{-1}]$',
                  fontsize=17)
    ax.legend(fontsize=13, frameon=False)
    ax.set_xlim(2450, 6000)
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=14)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()
    os.makedirs(PLOTDIR, exist_ok=True)
    out = os.path.join(PLOTDIR, 'transient_kernel_arm_fit.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")

    # ---- the kernel itself (binned to the data pixel scale, R ~ 500) ----
    dln_show = 2e-3
    fig, ax = plt.subplots(figsize=(10, 6.5))
    for logm, col in ((8.0, 'royalblue'), (9.5, 'darkorange')):
        kern = kernel_on_grid(v0, w0, np.sqrt(10.0 ** logm / MREF), -1,
                              dln=dln_show)
        n = len(kern)
        sgrid = (np.arange(n) - n // 2) * dln_show
        vk = (np.expm1(sgrid)) * 2.99792458e5 / 1e3
        ax.plot(vk, kern / kern.max(), drawstyle='steps-mid', color=col,
                lw=2.5, alpha=0.9,
                label=rf'$\rm \log M_{{BH}} = {logm:.1f}$')
    ax.axvline(0, color='gray', ls=':', lw=1.5)
    ax.set_xlabel(r'$\rm line\mbox{-}of\mbox{-}sight~velocity~[10^3~km~s^{-1}]$',
                  fontsize=18)
    ax.set_ylabel(r'$\rm kernel~(peak~normalized)$', fontsize=17)
    ax.legend(fontsize=14, frameon=False)
    ax.set_xlim(-60, 60)
    ax.set_xticks([-60, -40, -20, 0, 20, 40, 60])
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=14)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()
    out2 = os.path.join(PLOTDIR, 'arm_velocity_kernel.png')
    plt.savefig(out2, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out2}")


if __name__ == '__main__':
    main()
