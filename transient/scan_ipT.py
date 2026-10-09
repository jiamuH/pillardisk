#!/usr/bin/env python3
"""
scan_ipT.py - coarse chi^2 grid fit over (inclination, arm azimuth,
pillar peak temperature) for the first-principles blackbody model.
Fits the model difference spectrum (flare - quiescent) to the observed
2019 - 2000 difference at line-free rest wavelengths, with the quiescent
normalized to the 2000 spectrum at rest 5050-5300 A (per inclination).

Run:  python3 transient/scan_ipT.py
"""

import copy
import itertools
import os
import sys
import time

import numpy as np
import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from transient.transient_disk import disks_from_config  # noqa: E402
from transient.make_transient_figure import load_epoch, rebin_R  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(HERE, 'config_transient.yaml')) as fh:
    CFG = yaml.safe_load(fh)
CFG['transient']['occultation']['n_ray'] = 150  # coarser for speed
Z = CFG['observation']['redshift']
TO_FLAM = 1e-17 / (1.0 + Z)
NR, NPHI = 200, 400

# line-free rest wavelengths for the fit (avoid MgII, [OII], Hb, [OIII], ...)
WL_FIT = np.array([2500., 2600., 2650., 3000., 3050., 3200., 3400., 3550.,
                   3600., 3950., 4000., 4200., 4250., 4450., 4550., 4650.,
                   5100., 5200., 5300., 5500., 5700., 5900.])
NORM = (WL_FIT >= 5500) & (WL_FIT <= 6050)

I_GRID = [40., 46., 52., 58., 64., 70.]
PHI_GRID = [30., 45., 60., 75., 90.]
T_GRID = [8500., 9500., 10500., 12000.]


def data_diff():
    d = {}
    for name in ('sdss2000', 'eboss2019'):
        ep = load_epoch(name)
        wb, fb, eb = rebin_R(ep['wave_obs'], ep['flam'], ep['ivar'], R=400.0)
        rest = wb / (1.0 + Z)
        d[name] = (np.interp(WL_FIT, rest, fb),
                   np.interp(WL_FIT, rest, eb))
    diff = d['eboss2019'][0] - d['sdss2000'][0]
    err = np.sqrt(d['eboss2019'][1] ** 2 + d['sdss2000'][1] ** 2)
    return diff, err, d['sdss2000'][0]


def main():
    ddiff, derr, d2000 = data_diff()
    ndof = len(WL_FIT) - 3

    results = []
    t0 = time.time()
    for k, inc in enumerate(I_GRID):
        cfg_i = copy.deepcopy(CFG)
        cfg_i['observation']['inclination'] = inc
        _, quiet, _ = disks_from_config(cfg_i, nr=NR, nphi=NPHI)
        fq = quiet.compute_sed(WL_FIT) * TO_FLAM
        scale = np.median(d2000[NORM]) / np.median(fq[NORM])
        fq_s = fq * scale
        for phi, T in itertools.product(PHI_GRID, T_GRID):
            cfg_t = copy.deepcopy(cfg_i)
            cfg_t['transient']['pillar']['phi_pillar'] = float(np.radians(phi))
            cfg_t['transient']['pillar']['pillar_temp'] = T
            flare, _, _ = disks_from_config(cfg_t, nr=NR, nphi=NPHI)
            ff = flare.compute_sed(WL_FIT) * TO_FLAM * scale
            mdiff = ff - fq_s
            chi2 = float(np.sum(((mdiff - ddiff) / derr) ** 2))
            results.append((chi2 / ndof, inc, phi, T))
        print(f"  i={inc:.0f} done ({(k+1)}/{len(I_GRID)}, "
              f"{time.time()-t0:.0f} s)")

    results.sort(key=lambda r: r[0])
    print("\n  best (reduced chi^2):")
    print(f"  {'chi2/dof':>9}  {'i[deg]':>7}  {'phi[deg]':>8}  {'T[kK]':>6}")
    for chi2r, inc, phi, T in results[:8]:
        print(f"  {chi2r:>9.2f}  {inc:>7.0f}  {phi:>8.0f}  {T/1e3:>6.1f}")
    best = results[0]
    print(f"\n  BEST: i = {best[1]:.0f} deg, phi = {best[2]:.0f} deg, "
          f"T_peak = {best[3]/1e3:.1f} kK  (chi2/dof = {best[0]:.2f})")


if __name__ == '__main__':
    main()
