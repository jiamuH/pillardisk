#!/usr/bin/env python3
"""
fit_ring_heating.py - test the MULTI-RING TEMPERATURE-CHANGE idea for the
bump: instead of arm Balmer emission, the flare heats the disk's own
rings around some radius, T(r) -> T(r) * [1 + a * exp(-(r-r0)^2/2sr^2)].
The difference of the two multi-ring configurations is a sum of Planck
DERIVATIVES, which is narrower than any single blackbody. The raised-arm
occultation (adopted geometry, i = 60) is kept.

Model difference = [F_heated - F_quiet](scaled) - occ_def, with T(r) and
the per-ring emitting weights taken from the actual disk. No linear
amplitude: the heating amplitude a IS the amplitude.

Run:  python3 transient/fit_ring_heating.py
"""

import copy
import itertools
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from transient.tune_redside import data_on_grid, CFG, TO_FLAM, NR, NPHI  # noqa: E402
from transient.transient_disk import disks_from_config  # noqa: E402
from pillar_disk import C, DAY, FNU_TO_MJY  # noqa: E402

R0_GRID = [1.0, 1.5, 2.0, 3.0, 4.0, 6.0]
SR_GRID = [0.5, 1.0, 1.5, 2.5]
A_GRID = [0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.8, 1.2]


def main():
    WL, ddiff, derr, f2000, good = data_on_grid()
    nw = CFG['transient']['normalization']['rest_window']
    normsel = (WL >= nw[0]) & (WL <= nw[1])
    red = WL > 4500

    # adopted geometry: quiet disk + raised-arm occultation deficit
    c0 = copy.deepcopy(CFG)
    c0['transient']['pillar']['pillar_temp'] = 0.0
    flare0, quiet, _ = disks_from_config(c0, nr=NR, nphi=NPHI)
    q = quiet.compute_sed(WL) * TO_FLAM
    scale = np.median(f2000[normsel]) / np.median(q[normsel])
    occ0 = flare0.compute_sed(WL) * TO_FLAM * scale
    occ_def = q * scale - occ0

    # per-ring weights and temperatures of the QUIET disk (axisymmetric)
    weight, _ = quiet._sed_weights(False)          # (nr-1, nphi)
    W_r = weight.sum(axis=1)                       # (nr-1,)
    r_2d, phi_2d = np.meshgrid(quiet.r, quiet.phi, indexing='ij')
    T_r = quiet.get_temperature(r_2d, phi_2d)[:-1, :].mean(axis=1)
    rr = quiet.r[:-1]
    ld_to_cm = C * DAY
    norm = ld_to_cm ** 2 / (quiet.d * ld_to_cm) ** 2 * FNU_TO_MJY
    print(f"disk T(r): {np.interp(2.0, rr, T_r):.0f} K at 2 ld, "
          f"{np.interp(4.0, rr, T_r):.0f} K at 4 ld")

    def dflux(r0, sr, a):
        T2 = T_r * (1.0 + a * np.exp(-0.5 * ((rr - r0) / sr) ** 2))
        out = np.empty(WL.shape)
        for iw, w in enumerate(WL):
            out[iw] = np.sum(W_r * (quiet.planck_function(w, T2)
                                    - quiet.planck_function(w, T_r)))
        return out * norm * TO_FLAM * scale

    results = []
    for r0, sr, a in itertools.product(R0_GRID, SR_GRID, A_GRID):
        dF = dflux(r0, sr, a)
        model = dF - occ_def
        chi2 = float(np.sum(((model[good] - ddiff[good]) / derr[good]) ** 2))
        chi2_red = float(np.sum(((model[good & red] - ddiff[good & red])
                                 / derr[good & red]) ** 2))
        ipk = np.argmax(dF)
        tail = dF[np.argmin(np.abs(WL - 5500.))] / dF[ipk]
        results.append((chi2, chi2_red, r0, sr, a, WL[ipk], tail))

    ndof = good.sum() - 4
    results.sort(key=lambda r: r[0])
    print(f"\n  top 12 (heating shape red tail = dF(5500)/dF(peak); "
          f"single BB reference ~0.35-0.40):")
    print(f"  {'chi2/dof':>9} {'chi2_red':>9} {'r0[ld]':>7} {'sr':>5} "
          f"{'a':>5} {'peak[A]':>8} {'tail':>6}")
    for r in results[:12]:
        print(f"  {r[0]/ndof:>9.1f} {r[1]:>9.0f} {r[2]:>7.1f} {r[3]:>5.1f} "
              f"{r[4]:>5.2f} {r[5]:>8.0f} {r[6]:>6.3f}")

    best = results[0]
    dF = dflux(best[2], best[3], best[4])
    model = dF - occ_def
    print(f"\n  {'wl[A]':>7} {'data':>7} {'model':>7} {'-occ':>7}")
    for wtest in (2550., 3050., 3400., 3800., 4200., 4600., 5000., 5400.,
                  5800.):
        i = np.argmin(np.abs(WL - wtest))
        print(f"  {WL[i]:>7.0f} {ddiff[i]:>7.1f} {model[i]:>7.1f} "
              f"{-occ_def[i]:>7.1f}")


if __name__ == '__main__':
    main()
