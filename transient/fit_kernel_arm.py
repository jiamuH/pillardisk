#!/usr/bin/env python3
"""
fit_kernel_arm.py - replace the ad hoc Gaussian Doppler smearing of the
arm emission with the GEOMETRIC velocity kernel of the arm itself:
every arm cell is Doppler-shifted by its circular Keplerian
line-of-sight velocity (transient_disk.arm_velocity_field), so the
smearing is set by the arm geometry + inclination, with M_BH the only
new (free) parameter. v_los scales exactly as sqrt(M_BH), so one
velocity field serves every M_BH.

Intrinsic (rest-frame) arm emission shapes tried:
  - analytic recombination edge (_arm_emission_shape) with only a small
    LOCAL broadening (1,000 km/s; thermal + turbulence), scanning
    (T_e, f_bb, paschen);
  - the matter-bounded Cloudy slabs (arm_column grid), both faces.

Both kernel orientations (disk rotation sense) are fitted. The model
difference spectrum is s * A * (shape conv kernel) - occ_def as in
tune_redside / fit_cloudy_arm.

Run:  python3 transient/fit_kernel_arm.py
"""

import os
import sys

import numpy as np
from scipy.signal import fftconvolve

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from transient.tune_redside import setup  # noqa: E402
from transient.fit_cloudy_arm import load_cache  # noqa: E402

C_KMS = 2.99792458e5
WNORM = 3640.0
DLN = 2.5e-4                       # ln-lambda working grid (R = 4000)
MREF = 1e8                         # M_BH at which the field is computed
LOCAL_VS = 1000.0                  # local (thermal/turbulent) broadening

LOGM_GRID = np.arange(7.5, 10.0, 0.25)
TE_GRID = [7000., 8000., 9000., 10000., 11000., 12000., 14000.]
FBB_GRID = [0.0, 0.025, 0.05]
PAS_GRID = [0.0, 0.04]


def kernel_on_grid(v0, w0, mfac, orient, dln=DLN):
    """Doppler kernel on the uniform ln-lambda grid: histogram of
    s = ln(1 + beta) over arm cells, beta = orient * v0 * mfac / c.
    Returns an odd-length array centered on zero shift, summing to 1.
    dln sets the bin width (default: the fine convolution grid; pass a
    coarser value, e.g. matched to the data pixels, for display)."""
    beta = orient * v0 * mfac / C_KMS
    ok = np.abs(beta) < 0.7
    s = np.log1p(beta[ok])
    w = w0[ok]
    off = np.rint(s / dln).astype(int)
    nk = max(abs(off.min()), abs(off.max())) + 1
    kern = np.bincount(off + nk, weights=w, minlength=2 * nk + 1)
    return kern / kern.sum()


def main():
    WL, ddiff, derr, good, occ_def, A, flare_ref = setup()
    w_fit = 1.0 / derr[good] ** 2
    red = WL > 4500

    # geometric velocity field (scales as sqrt(M_BH / MREF))
    v0, w0 = flare_ref.arm_velocity_field(MREF)
    vbar = np.sum(v0 * w0) / np.sum(w0)
    vrms = np.sqrt(np.sum(v0 ** 2 * w0) / np.sum(w0))
    print(f"arm velocity field at M_BH = {MREF:.0e} Msun: "
          f"mean = {vbar:+.0f} km/s, rms = {vrms:.0f} km/s "
          f"(both scale as sqrt(M))")

    # fine ln-lambda grid
    g = np.arange(np.log(1500.0), np.log(11000.0), DLN)
    wg = np.exp(g)

    # intrinsic shapes: analytic edge family
    shapes = {}
    for te in TE_GRID:
        for fbb in FBB_GRID:
            for pas in PAS_GRID:
                key = ('edge', te, fbb, pas)
                shapes[key] = np.array([flare_ref._arm_emission_shape(
                    w, te, fbb, LOCAL_VS, pas) for w in wg])
    # intrinsic shapes: matter-bounded Cloudy slabs
    cwave, cfaces, cparams = load_cache('arm_column')
    cp = list(cparams.values())
    for face in cfaces:
        flam = cfaces[face] / cwave[None, :]
        for ipt in range(len(cp[0])):
            key = ('cloudy', face, cp[0][ipt], cp[1][ipt], cp[2][ipt])
            shapes[key] = np.interp(wg, cwave, flam[ipt])

    results = []
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
                num = np.sum((ddiff[good] + occ_def[good]) * E[good] * w_fit)
                den = np.sum(E[good] ** 2 * w_fit)
                s = max(0.0, num / den)
                model = s * E - occ_def
                chi2 = float(np.sum(
                    ((model[good] - ddiff[good]) / derr[good]) ** 2))
                chi2_red = float(np.sum(
                    ((model[good & red] - ddiff[good & red])
                     / derr[good & red]) ** 2))
                results.append((chi2, chi2_red, logm, orient, key, s))

    ndof = good.sum() - 5

    def show(rows):
        print(f"  {'chi2/dof':>9} {'chi2_red':>9} {'logM':>5} {'rot':>4} "
              f"{'amp':>6}  intrinsic shape")
        for r in rows:
            key = r[4]
            if key[0] == 'edge':
                desc = (f"edge T_e={key[1]/1e3:.0f}kK f_bb={key[2]:.3f} "
                        f"pas={key[3]:.2f}")
            else:
                desc = (f"cloudy {key[1]} phi={key[2]:.0f} n={key[3]:.0f} "
                        f"N={key[4]:.1f}")
            print(f"  {r[0]/ndof:>9.1f} {r[1]:>9.0f} {r[2]:>5.2f} "
                  f"{'+' if r[3] > 0 else '-':>4} {r[5]:>6.2f}  {desc}")

    results.sort(key=lambda r: r[0])
    print(f"\n  top 12 by TOTAL chi2 (all shapes):")
    show(results[:12])
    edge_res = [r for r in results if r[4][0] == 'edge']
    cl_res = [r for r in results if r[4][0] == 'cloudy']
    print(f"\n  top 6 analytic-edge:")
    show(edge_res[:6])
    print(f"\n  top 6 Cloudy:")
    show(cl_res[:6])
    print(f"\n  top 6 by RED-SIDE chi2:")
    show(sorted(results, key=lambda r: r[1])[:6])

    # residual profile of the overall best
    chi2, chi2_red, logm, orient, key, s = results[0][:6]
    kern = kernel_on_grid(v0, w0, np.sqrt(10.0 ** logm / MREF), orient)
    F = fftconvolve(shapes[key], kern, mode='same')
    shape = np.interp(WL, wg, F) / np.interp(WNORM, wg, F)
    model = s * A * shape - occ_def
    print(f"\n  best: logM = {logm:.2f}, rot = {'+' if orient > 0 else '-'},"
          f" {key}  (chi2/dof = {chi2/ndof:.1f}; "
          f"kernel rms = {vrms * np.sqrt(10**logm/MREF):.0f} km/s)")
    print(f"  {'wl[A]':>7} {'data':>7} {'model':>7} {'-occ':>7}")
    for wtest in (2550., 3050., 3400., 3800., 4200., 4600., 5000., 5400.,
                  5800.):
        i = np.argmin(np.abs(WL - wtest))
        print(f"  {WL[i]:>7.0f} {ddiff[i]:>7.1f} {model[i]:>7.1f} "
              f"{-occ_def[i]:>7.1f}")


if __name__ == '__main__':
    main()
