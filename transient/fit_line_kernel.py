#!/usr/bin/env python3
"""
fit_line_kernel.py - fit the ARM's geometric velocity kernel to the
observed RESPONDING Hbeta profile (the continuum-subtracted 2019 - 2000
difference), with M_BH free. This asks: if the flare's Hbeta response
came from the arm itself, what M_BH (velocity scale) does its width
imply, and does the kernel's geometric shape (asymmetry) match?

The kernel shape is fixed by the arm geometry + inclination per arm
azimuth; the velocity scale is sqrt(M_BH); the line is assumed
intrinsically narrow (local sigma 500 km/s) so the profile IS the
kernel. Amplitude is solved linearly. [O III] 4959/5007 windows are
masked on the red side.

Run:  python3 transient/fit_line_kernel.py
"""

import copy
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
MREF = 1e8
LAM_HB = 4861.32
HB_WINS = [(-16000., -11500.), (11500., 17000.)]
MASKS = [(5200., 7000.), (8000., 10500.)]     # [O III] 4959 / 5007
VFIT = 11000.0                                 # fit window half-range
SIG_LOCAL = 500.0                              # intrinsic line width (km/s)

PHI_ARM_DEG = [0., 30., 60., 90., 120., 150.]
LOGM_GRID = np.arange(5.5, 9.01, 0.125)
NR, NPHI = 200, 400


def model_profile(v_grid, v0, w0, mfac, orient):
    """Arm line profile on the data velocity grid: histogram of the
    scaled cell velocities, then a small Gaussian local broadening."""
    vc = orient * v0 * mfac
    dv = np.median(np.diff(v_grid))
    edges = np.concatenate([v_grid - dv / 2, [v_grid[-1] + dv / 2]])
    prof, _ = np.histogram(vc, bins=edges, weights=w0)
    nk = int(4 * SIG_LOCAL / dv) + 1
    k = np.exp(-0.5 * (np.arange(-nk, nk + 1) * dv / SIG_LOCAL) ** 2)
    return np.convolve(prof, k / k.sum(), mode='same')


def main():
    v, sub, err = subtracted_profiles(LAM_HB, HB_WINS)
    diff = sub['eboss2019'] - sub['sdss2000']
    derr = np.sqrt(err['eboss2019'] ** 2 + err['sdss2000'] ** 2)
    good = np.abs(v) < VFIT
    for v1, v2 in MASKS:
        good &= ~((v >= v1) & (v <= v2))
    w_fit = 1.0 / derr[good] ** 2
    fw_obs, vpk_obs = fwhm_of(v, diff)
    print(f"observed responding Hbeta: FWHM = {fw_obs:.0f} km/s, "
          f"peak at {vpk_obs:+.0f} km/s")

    fields = {}
    results = []
    for phi_deg in PHI_ARM_DEG:
        cfg = copy.deepcopy(tr.CFG)
        cfg['transient']['pillar']['phi_pillar'] = float(np.radians(phi_deg))
        flare, _, _ = disks_from_config(cfg, nr=NR, nphi=NPHI)
        v0, w0 = flare.arm_velocity_field(MREF)
        fields[phi_deg] = (v0, w0)
        for logm in LOGM_GRID:
            mfac = np.sqrt(10.0 ** logm / MREF)
            for orient in (+1, -1):
                prof = model_profile(v, v0, w0, mfac, orient)
                num = np.sum(diff[good] * prof[good] * w_fit)
                den = np.sum(prof[good] ** 2 * w_fit)
                s = max(0.0, num / den)
                chi2 = float(np.sum(
                    ((s * prof[good] - diff[good]) / derr[good]) ** 2))
                results.append((chi2, phi_deg, logm, orient, s))
        print(f"  phi_arm = {phi_deg:.0f} done", flush=True)

    ndof = good.sum() - 4
    results.sort(key=lambda r: r[0])
    print(f"\n  ndata = {good.sum()}")
    print(f"  {'chi2/dof':>9} {'phi':>5} {'logM':>6} {'rot':>4} {'amp':>9}")
    for r in results[:10]:
        print(f"  {r[0]/ndof:>9.2f} {r[1]:>5.0f} {r[2]:>6.2f} "
              f"{'+' if r[3] > 0 else '-':>4} {r[4]:>9.2e}")

    best = results[0]
    _, bphi, blogm, borient, bs = best
    v0, w0 = fields[bphi]
    prof_b = bs * model_profile(v, v0, w0,
                                np.sqrt(10.0 ** blogm / MREF), borient)
    fw_mod, vpk_mod = fwhm_of(v, prof_b)
    print(f"\n  best: phi_arm = {bphi:.0f} deg, logM = {blogm:.2f}, "
          f"rot = {'+' if borient > 0 else '-'} "
          f"(chi2/dof = {best[0]/ndof:.2f})")
    print(f"  model FWHM = {fw_mod:.0f} km/s, peak at {vpk_mod:+.0f} km/s")

    # a contrasting mass for the figure (10x heavier)
    prof_h = model_profile(v, v0, w0,
                           np.sqrt(10.0 ** (blogm + 1) / MREF), borient)
    num = np.sum(diff[good] * prof_h[good] * w_fit)
    den = np.sum(prof_h[good] ** 2 * w_fit)
    prof_h *= max(0.0, num / den)

    fig, ax = plt.subplots(figsize=(10, 7))
    ax.fill_between(v / 1e3, diff - derr, diff + derr, step='mid',
                    color='seagreen', alpha=0.2, lw=0)
    ax.plot(v / 1e3, diff, drawstyle='steps-mid', color='seagreen', lw=2,
            alpha=0.9, label=r'$\rm H\beta~response~(2019-2000)$')
    ax.plot(v / 1e3, prof_b, '-', color='crimson', lw=3, alpha=0.95,
            label=(rf'$\rm arm~kernel,~\log M={blogm:.1f},~'
                   rf'\phi_{{\rm arm}}={bphi:.0f}^\circ$'))
    ax.plot(v / 1e3, prof_h, '--', color='darkorange', lw=2.2, alpha=0.9,
            label=rf'$\rm arm~kernel,~\log M={blogm + 1:.1f}$')
    for v1, v2 in MASKS:
        ax.axvspan(v1 / 1e3, v2 / 1e3, color='gray', alpha=0.12, lw=0)
    ax.axhline(0, color='gray', lw=1)
    ax.axvline(0, color='gray', ls=':', lw=1.5)
    ax.set_xlabel(r'$\rm velocity~[10^3~km~s^{-1}]$', fontsize=18)
    ax.set_ylabel(
        r'$\Delta f_\lambda~[\rm 10^{-17}~erg~s^{-1}~cm^{-2}~\AA^{-1}]$',
        fontsize=17)
    ax.set_xlim(-12, 12)
    ax.set_xticks([-12, -9, -6, -3, 0, 3, 6, 9, 12])
    ax.legend(fontsize=13, frameon=False, loc='upper right')
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=14)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()
    os.makedirs(PLOTDIR, exist_ok=True)
    out = os.path.join(PLOTDIR, 'line_profile_arm_model.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved {out}")


if __name__ == '__main__':
    main()
