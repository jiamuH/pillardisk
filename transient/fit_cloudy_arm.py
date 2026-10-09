#!/usr/bin/env python3
"""
fit_cloudy_arm.py - fit the transient arm emission with SELF-CONSISTENT
Cloudy spectra (continuum + lines from the same photoionized cloud)
instead of the analytic smeared-Balmer-edge shape.

For every grid point and both cloud faces (illuminated = reflected,
shielded = outward diffuse), the emitted nu*F_nu spectrum is converted
to an f_lambda shape, Doppler-smeared by the arm's orbital motion, and
fitted to the observed 2019-2000 difference spectrum with the
linear-amplitude machinery of tune_redside.setup():
    model = s * A * shape(w) - occ_def.

Run:  python3 transient/fit_cloudy_arm.py [loc_metal_flux|arm_column]
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from transient.tune_redside import setup  # noqa: E402
from transient.extract_cloudy_arm import cache_path, GRIDS  # noqa: E402

VS_GRID = [5000., 10000., 15000., 20000., 25000.]   # km/s (physical range)
C_KMS = 2.99792458e5
WNORM = 3640.0        # shape normalization wavelength (A)


def smear_loggrid(wave, flam, vsmear_kms):
    """Gaussian Doppler smear (sigma = v/c in ln-lambda) on a uniform
    ln-lambda grid; returns (wave_out, flam_out) on that grid."""
    lnw = np.log(wave)
    dln = 2.5e-4                                   # R = 4000 working grid
    g = np.arange(lnw[0], lnw[-1], dln)
    f = np.interp(g, lnw, flam)
    sig = vsmear_kms / C_KMS
    if sig > 0:
        nk = int(5 * sig / dln)
        k = np.exp(-0.5 * (np.arange(-nk, nk + 1) * dln / sig) ** 2)
        k /= k.sum()
        f = np.convolve(f, k, mode='same')
    return np.exp(g), f


def load_cache(grid_name):
    """Returns (wave, faces dict, params dict) for a cached grid."""
    d = np.load(cache_path(grid_name))
    par3 = GRIDS[grid_name]['par3']
    faces = {'illum': d['refl_cont'] + d['refl_line'],
             'shield': d['out_cont'] + d['out_line']}
    params = {'phi': d['phi'], 'n_H': d['hden'], par3: d[par3]}
    return d['wave'], faces, params


def main():
    grid_name = sys.argv[1] if len(sys.argv) > 1 else 'loc_metal_flux'
    wave, faces, params = load_cache(grid_name)
    pnames = list(params)           # ['phi', 'n_H', <par3>]
    pvals = list(params.values())
    npts = len(pvals[0])

    WL, ddiff, derr, good, occ_def, A, _ = setup()
    w_fit = 1.0 / derr[good] ** 2
    red = WL > 4500

    results = []
    for face, spec in faces.items():
        # nu*F_nu -> f_lambda shape (arbitrary norm before division at WNORM)
        flam_all = spec / wave[None, :]
        for ipt in range(npts):
            for vs in VS_GRID:
                wg, fg = smear_loggrid(wave, flam_all[ipt], vs)
                fnorm = np.interp(WNORM, wg, fg)
                if fnorm <= 0:
                    continue
                shape = np.interp(WL, wg, fg) / fnorm
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
                results.append((chi2, chi2_red, face,
                                pvals[0][ipt], pvals[1][ipt], pvals[2][ipt],
                                vs, s, ipt))

    ndof = good.sum() - 4
    header = (f"  {'chi2/dof':>9} {'chi2_red':>9} {'face':>7} "
              f"{pnames[0]:>5} {pnames[1]:>6} {pnames[2]:>7} "
              f"{'v[1e3km/s]':>11} {'amp':>8}")

    def show(rows):
        print(header)
        for r in rows:
            print(f"  {r[0]/ndof:>9.1f} {r[1]:>9.0f} {r[2]:>7} {r[3]:>5.1f} "
                  f"{r[4]:>6.1f} {r[5]:>7.1f} {r[6]/1e3:>11.0f} {r[7]:>8.2f}")

    results.sort(key=lambda r: r[0])
    print(f"grid '{grid_name}', ndata = {good.sum()}, "
          f"red-side ndata = {(good & red).sum()}")
    print(f"\n  top 15 by TOTAL chi2:")
    show(results[:15])
    print(f"\n  top 8 by RED-SIDE chi2:")
    show(sorted(results, key=lambda r: r[1])[:8])

    # residual profile of the overall best
    chi2, chi2_red, face, b1, b2, b3, bvs, bs, ipt = results[0]
    wg, fg = smear_loggrid(wave, faces[face][ipt] / wave, bvs)
    shape = np.interp(WL, wg, fg) / np.interp(WNORM, wg, fg)
    model = bs * A * shape - occ_def
    print(f"\n  best: {face} face, {pnames[0]}={b1:.1f}, {pnames[1]}={b2:.1f},"
          f" {pnames[2]}={b3:.1f}, v={bvs/1e3:.0f}e3 km/s"
          f" (chi2/dof={chi2/ndof:.1f})")
    print(f"  {'wl[A]':>7} {'data':>7} {'model':>7} {'-occ':>7}")
    for wtest in (2550., 3050., 3400., 3800., 4200., 4600., 5000., 5400.,
                  5800.):
        i = np.argmin(np.abs(WL - wtest))
        print(f"  {WL[i]:>7.0f} {ddiff[i]:>7.1f} {model[i]:>7.1f} "
              f"{-occ_def[i]:>7.1f}")


if __name__ == '__main__':
    main()
