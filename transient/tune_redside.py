#!/usr/bin/env python3
"""
tune_redside.py - scan the arm-emission shape parameters (balmer_te,
balmer_bb_frac, balmer_vsmear, balmer_paschen) against the observed
2019 - 2000 difference spectrum, targeting the red-side residual
(model must fall to ~0 redward of 4500 A) without losing the fit
at 3800-4600 A.

Uses the separability of the balmer-mode arm emission in compute_sed:
    E(w) = area_w * amp * shape(w),
where area_w (geometry) and amp (pillar_temp) are scalars. The disk is
built ONCE (occultation deficit + reference emission); each grid point
only re-evaluates shape(w) and solves the emission amplitude linearly.

Run:  python3 transient/tune_redside.py
"""

import copy
import itertools
import os
import sys

import numpy as np
import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from transient.transient_disk import disks_from_config  # noqa: E402
from transient.make_transient_figure import (  # noqa: E402
    load_epoch, rebin_R, LINE_WINDOWS)

HERE = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(HERE, 'config_transient.yaml')) as fh:
    CFG = yaml.safe_load(fh)
CFG['transient']['occultation']['n_ray'] = 150
Z = CFG['observation']['redshift']
TO_FLAM = 1e-17 / (1.0 + Z)
NR, NPHI = 200, 400

TE_GRID = [7000., 8000., 9000., 10000., 11000., 12000., 13000.]
FBB_GRID = [0.0, 0.025, 0.05, 0.075, 0.10]
VS_GRID = [45000., 55000., 65000., 75000., 90000.]
PAS_GRID = [0.0, 0.02, 0.04]


def data_on_grid():
    b = {}
    for n in ('sdss2000', 'eboss2019'):
        ep = load_epoch(n)
        wb, fb, eb = rebin_R(ep['wave_obs'], ep['flam'], ep['ivar'], R=500.0)
        b[n] = (wb / (1.0 + Z), fb, eb)
    grid = b['sdss2000'][0]
    sel = (grid > 2450) & (grid < 6000)
    WL = grid[sel]
    f2000 = b['sdss2000'][1][sel]
    f2019 = np.interp(WL, b['eboss2019'][0], b['eboss2019'][1])
    e2000 = b['sdss2000'][2][sel]
    e2019 = np.interp(WL, b['eboss2019'][0], b['eboss2019'][2])
    good = np.ones_like(WL, dtype=bool)
    for w1, w2 in LINE_WINDOWS:
        good &= ~((WL >= w1) & (WL <= w2))
    return WL, f2019 - f2000, np.sqrt(e2000 ** 2 + e2019 ** 2), f2000, good


def setup(cfg=None):
    """Data difference + one-time geometry: returns (WL, ddiff, derr, good,
    occ_def, A, flare_ref) where the arm emission for any shape parameters
    is E(w) = A * flare_ref._arm_emission_shape(w, te, fbb, vs, pas) and
    the model difference spectrum is s * E - occ_def. Pass cfg to override
    the module-level config (e.g. a different arm azimuth or inclination)."""
    cfg = CFG if cfg is None else cfg
    WL, ddiff, derr, f2000, good = data_on_grid()
    nw = cfg['transient']['normalization']['rest_window']
    normsel = (WL >= nw[0]) & (WL <= nw[1])

    # geometry (built once): quiet, occultation deficit, prefactor
    c0 = copy.deepcopy(cfg)
    c0['transient']['pillar']['pillar_temp'] = 0.0
    flare0, quiet, _ = disks_from_config(c0, nr=NR, nphi=NPHI)
    q = quiet.compute_sed(WL) * TO_FLAM
    scale = np.median(f2000[normsel]) / np.median(q[normsel])
    q *= scale
    occ0 = flare0.compute_sed(WL) * TO_FLAM * scale
    occ_def = q - occ0

    p0 = cfg['transient']['pillar']
    cr = copy.deepcopy(cfg)
    flare_ref, _, _ = disks_from_config(cr, nr=NR, nphi=NPHI)
    E_ref = flare_ref.compute_sed(WL) * TO_FLAM * scale - occ0
    shape_ref = np.array([flare_ref._arm_emission_shape(
        w, p0['balmer_te'], p0['balmer_bb_frac'], p0['balmer_vsmear'],
        p0.get('balmer_paschen', 0.12)) for w in WL])
    ratio = E_ref / shape_ref
    A = np.median(ratio)
    assert np.std(ratio) / A < 1e-6, "emission not separable?!"
    return WL, ddiff, derr, good, occ_def, A, flare_ref


def main():
    WL, ddiff, derr, good, occ_def, A, flare_ref = setup()

    # ---- scan ----
    w_fit = 1.0 / derr[good] ** 2
    red = WL > 4500
    results = []
    for te, fbb, vs, pas in itertools.product(TE_GRID, FBB_GRID,
                                              VS_GRID, PAS_GRID):
        shape = np.array([flare_ref._arm_emission_shape(w, te, fbb, vs, pas)
                          for w in WL])
        E = A * shape
        num = np.sum((ddiff[good] + occ_def[good]) * E[good] * w_fit)
        den = np.sum(E[good] ** 2 * w_fit)
        s = max(0.0, num / den)
        model = s * E - occ_def
        resid = (model[good] - ddiff[good]) / derr[good]
        chi2 = float(np.sum(resid ** 2))
        chi2_red = float(np.sum(((model[good & red] - ddiff[good & red])
                                 / derr[good & red]) ** 2))
        results.append((chi2, chi2_red, te, fbb, vs, pas, s, model))

    ndof = good.sum() - 5
    results.sort(key=lambda r: r[0])
    print(f"ndata = {good.sum()}, red-side ndata = {(good & red).sum()}")
    print(f"\n  top 12 by TOTAL chi2 (chi2/dof, red-side chi2):")
    print(f"  {'chi2/dof':>9} {'chi2_red':>9} {'T_e[kK]':>8} {'f_bb':>5} "
          f"{'v[1e3km/s]':>11} {'paschen':>8} {'amp':>6}")
    for r in results[:12]:
        print(f"  {r[0]/ndof:>9.1f} {r[1]:>9.0f} {r[2]/1e3:>8.0f} "
              f"{r[3]:>5.2f} {r[4]/1e3:>11.0f} {r[5]:>8.2f} {r[6]:>6.2f}")

    by_red = sorted(results, key=lambda r: r[1])
    print(f"\n  top 8 by RED-SIDE chi2:")
    for r in by_red[:8]:
        print(f"  {r[0]/ndof:>9.1f} {r[1]:>9.0f} {r[2]/1e3:>8.0f} "
              f"{r[3]:>5.2f} {r[4]/1e3:>11.0f} {r[5]:>8.2f} {r[6]:>6.2f}")

    # residual profile of the overall best
    best = results[0]
    model = best[7]
    print(f"\n  best model vs data (occultation deficit shown separately):")
    print(f"  {'wl[A]':>7} {'data':>7} {'model':>7} {'-occ':>7}")
    for wtest in (2550., 3050., 3400., 3800., 4200., 4600., 5000., 5400.,
                  5800.):
        i = np.argmin(np.abs(WL - wtest))
        print(f"  {WL[i]:>7.0f} {ddiff[i]:>7.1f} {model[i]:>7.1f} "
              f"{-occ_def[i]:>7.1f}")


if __name__ == '__main__':
    main()
