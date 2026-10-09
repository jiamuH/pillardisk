#!/usr/bin/env python3
"""
fit_kernel_phiarm.py - joint scan of the ARM AZIMUTH and M_BH for the
geometric-kernel arm model (fit_kernel_arm.py). Moving the arm in
azimuth changes BOTH the Doppler kernel (an arm at phi = +-90 deg sits
at maximum line-of-sight velocity) and the occultation deficit (deepest
for the near-side arm at phi = 0), so the geometry is rebuilt per
azimuth.

Intrinsic shapes kept small: the best analytic-edge family and the best
matter-bounded Cloudy slabs from the previous scans.

Run:  python3 transient/fit_kernel_phiarm.py
"""

import copy
import os
import sys

import numpy as np
from scipy.signal import fftconvolve

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import transient.tune_redside as tr  # noqa: E402
from transient.fit_cloudy_arm import load_cache  # noqa: E402
from transient.fit_kernel_arm import (  # noqa: E402
    kernel_on_grid, C_KMS, WNORM, DLN, MREF, LOCAL_VS)

PHI_ARM_DEG = [0., 30., 60., 90., 120., 150.]
LOGM_GRID = np.arange(7.5, 10.0, 0.25)
EDGE_SET = [(8000., 0.0), (8000., 0.025), (8000., 0.05),
            (10000., 0.025), (12000., 0.025)]        # (T_e, f_bb), pas = 0
CLOUDY_SET = [('shield', 21., 12., 22.5), ('shield', 21., 12., 23.0),
              ('illum', 17., 9., 21.5), ('shield', 18., 10., 22.5)]


def main():
    # intrinsic shapes on the fine grid (geometry-independent)
    g = np.arange(np.log(1500.0), np.log(11000.0), DLN)
    wg = np.exp(g)
    cwave, cfaces, cparams = load_cache('arm_column')
    cp = list(cparams.values())

    results = []
    for phi_deg in PHI_ARM_DEG:
        cfg = copy.deepcopy(tr.CFG)
        cfg['transient']['pillar']['phi_pillar'] = float(np.radians(phi_deg))
        WL, ddiff, derr, good, occ_def, A, flare_ref = tr.setup(cfg)
        w_fit = 1.0 / derr[good] ** 2
        red = WL > 4500

        shapes = {}
        for te, fbb in EDGE_SET:
            shapes[('edge', te, fbb)] = np.array(
                [flare_ref._arm_emission_shape(w, te, fbb, LOCAL_VS, 0.0)
                 for w in wg])
        for face, p, n, N in CLOUDY_SET:
            ipt = int(np.argmin(np.abs(cp[0] - p) + np.abs(cp[1] - n)
                                + np.abs(cp[2] - N)))
            shapes[('cloudy', face, p, n, N)] = np.interp(
                wg, cwave, cfaces[face][ipt] / cwave)

        v0, w0 = flare_ref.arm_velocity_field(MREF)
        vrms = np.sqrt(np.sum(v0 ** 2 * w0) / np.sum(w0))
        occ_uv = float(np.interp(2550.0, WL, occ_def))
        print(f"phi_arm = {phi_deg:5.0f} deg: kernel rms(1e8) = "
              f"{vrms:6.0f} km/s, occultation deficit at 2550 A = "
              f"{occ_uv:5.1f}", flush=True)

        for logm in LOGM_GRID:
            mfac = np.sqrt(10.0 ** logm / MREF)
            for orient in (+1, -1):
                kern = kernel_on_grid(v0, w0, mfac, orient)
                for key, S in shapes.items():
                    F = fftconvolve(S, kern, mode='same')
                    fnorm = np.interp(WNORM, wg, F)
                    if fnorm <= 0:
                        continue
                    shape = np.interp(WL, wg, F) / fnorm
                    E = A * shape
                    num = np.sum((ddiff[good] + occ_def[good])
                                 * E[good] * w_fit)
                    den = np.sum(E[good] ** 2 * w_fit)
                    s = max(0.0, num / den)
                    model = s * E - occ_def
                    chi2 = float(np.sum(
                        ((model[good] - ddiff[good]) / derr[good]) ** 2))
                    chi2_red = float(np.sum(
                        ((model[good & red] - ddiff[good & red])
                         / derr[good & red]) ** 2))
                    results.append((chi2, chi2_red, phi_deg, logm, orient,
                                    key, s))
        ndata = good.sum()

    ndof = ndata - 6

    def show(rows):
        print(f"  {'chi2/dof':>9} {'chi2_red':>9} {'phi':>5} {'logM':>5} "
              f"{'rot':>4} {'amp':>6}  intrinsic shape")
        for r in rows:
            key = r[5]
            if key[0] == 'edge':
                desc = f"edge T_e={key[1]/1e3:.0f}kK f_bb={key[2]:.3f}"
            else:
                desc = (f"cloudy {key[1]} phi={key[2]:.0f} n={key[3]:.0f} "
                        f"N={key[4]:.1f}")
            print(f"  {r[0]/ndof:>9.1f} {r[1]:>9.0f} {r[2]:>5.0f} "
                  f"{r[3]:>5.2f} {'+' if r[4] > 0 else '-':>4} "
                  f"{r[6]:>6.2f}  {desc}")

    results.sort(key=lambda r: r[0])
    print(f"\n  top 15 by TOTAL chi2:")
    show(results[:15])
    print(f"\n  top 8 by RED-SIDE chi2:")
    show(sorted(results, key=lambda r: r[1])[:8])


if __name__ == '__main__':
    main()
