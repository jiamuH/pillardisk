#!/usr/bin/env python3
"""
plot_sim_hr.py - disk aspect ratio H/r versus radius from the Athena dump,
measured two independent ways:

  1. thermal:  H/r = c_s / v_K = sqrt(cs^2 * r)   (code units, v_K = r^-1/2)
     using the phi-median midplane cs^2 = P/rho;
  2. density:  H/r = sigma_theta, the density-weighted standard deviation
     of (theta - pi/2) per (phi, r) column (exact for a Gaussian vertical
     profile), phi-median with a 5-95 percentile band showing the
     azimuthal variation (arm thickening).

Run:  python3 transient/plot_sim_hr.py [dump.athdf]
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from transient.sim_spectrum import load_sim, PLOTDIR  # noqa: E402

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else \
        os.path.join(HERE, 'data', 'sim', 'disk.out1.00012.athdf')
    sim = load_sim(path)
    r, th, ph = sim['r'], sim['th'], sim['ph']
    rho, press = sim['rho'], sim['press']            # (ph, th, r)
    jm = np.argmin(np.abs(th - np.pi / 2))

    # thermal H/r from the midplane sound speed
    cs2_mid = np.median(press[:, jm, :] / rho[:, jm, :], axis=0)   # (r,)
    hr_cs = np.sqrt(cs2_mid * r)

    # density H/r: sigma_theta per (phi, r) column
    dth = np.gradient(th)
    wgt = rho * dth[None, :, None]
    z = (th - np.pi / 2)[None, :, None]
    sig = np.sqrt((wgt * z ** 2).sum(axis=1) / wgt.sum(axis=1))    # (ph, r)
    hr_med = np.median(sig, axis=0)
    hr_lo, hr_hi = np.percentile(sig, [5, 95], axis=0)

    fig, ax = plt.subplots(figsize=(11, 7))
    ax.fill_between(r, hr_lo, hr_hi, color='royalblue', alpha=0.18, lw=0,
                    label=r'$\rm azimuthal~5\mbox{-}95\%~spread$')
    ax.plot(r, hr_med, '-', color='royalblue', lw=3, alpha=0.9,
            label=r'$\rm density~scale~height~\sigma_\theta$')
    ax.plot(r, hr_cs, '--', color='orangered', lw=3, alpha=0.9,
            label=r'$c_s/v_K~\rm (midplane)$')
    ax.set_xlabel(r'$r~[r_0]$', fontsize=18)
    ax.set_ylabel(r'$H/r$', fontsize=18)
    ax.set_xlim(r.min(), r.max())
    ax.legend(fontsize=14, frameon=False)
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=14)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()
    os.makedirs(PLOTDIR, exist_ok=True)
    out = os.path.join(PLOTDIR, 'sim_hr_profile.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"H/r at r=1: cs-based {np.interp(1.0, r, hr_cs):.3f}, "
          f"density-based {np.interp(1.0, r, hr_med):.3f}")
    print(f"Saved {out}")

    # ---- 2D face-on map of the local scale height H/r(phi, r) ----
    from transient.sim_spectrum import add_colorbar, mark_star
    R, P = np.meshgrid(r, ph)
    X, Y = R * np.cos(P), R * np.sin(P)
    fig, ax = plt.subplots(figsize=(9.5, 8))
    pc = ax.pcolormesh(X, Y, sig, cmap='viridis', shading='auto')
    ax.set_aspect('equal')
    add_colorbar(pc, ax, r'$H/r~(\sigma_\theta)$')
    mark_star(ax, sim)
    ax.set_xlabel(r'$x~[r_0]$', fontsize=16)
    ax.set_ylabel(r'$y~[r_0]$', fontsize=16)
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=13)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()
    out = os.path.join(PLOTDIR, 'sim_hr_map.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")


if __name__ == '__main__':
    main()
