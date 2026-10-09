#!/usr/bin/env python3
"""
diag_redside.py - decompose the RED SIDE (> 4500 A rest) of the model
2019-2000 difference spectrum at the current config_transient.yaml
parameters, to identify which component produces the flat ~10 (1e-17 cgs)
residual where the data go to zero.

Components (all from compute_sed, no free additive curves):
  - occultation deficit  : quiet - flare(pillar_temp=0)   [negative in diff]
  - arm recombination    : (1-fbb) * [flare(fbb=0)  - flare(pillar_temp=0)]
  - arm thermalized BB   :   fbb   * [flare(fbb=1)  - flare(pillar_temp=0)]
Within the recombination part, the Paschen term is isolated by evaluating
_arm_emission_shape with the module-level Paschen ratio zeroed.

Run:  python3 transient/diag_redside.py
"""

import copy
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


def sed_with(fbb=None, pillar_temp=None):
    c = copy.deepcopy(CFG)
    if fbb is not None:
        c['transient']['pillar']['balmer_bb_frac'] = fbb
    if pillar_temp is not None:
        c['transient']['pillar']['pillar_temp'] = pillar_temp
    flare, quiet, _ = disks_from_config(c, nr=NR, nphi=NPHI)
    return flare, quiet


def main():
    WL, ddiff, derr, f2000, good = data_on_grid()
    nw = CFG['transient']['normalization']['rest_window']
    normsel = (WL >= nw[0]) & (WL <= nw[1])

    fbb = CFG['transient']['pillar']['balmer_bb_frac']

    # quiet + occultation-only (arm raised, not heated)
    flare0, quiet = sed_with(pillar_temp=0.0)
    q = quiet.compute_sed(WL) * TO_FLAM
    scale = np.median(f2000[normsel]) / np.median(q[normsel])
    q *= scale
    occ0 = flare0.compute_sed(WL) * TO_FLAM * scale
    occ_def = q - occ0                       # >0; enters diff as -occ_def

    # full model and fbb end-members
    flare_full, _ = sed_with()
    full = flare_full.compute_sed(WL) * TO_FLAM * scale
    flare_rec, _ = sed_with(fbb=0.0)
    E_rec = flare_rec.compute_sed(WL) * TO_FLAM * scale - occ0
    flare_bb, _ = sed_with(fbb=1.0)
    E_bb = flare_bb.compute_sed(WL) * TO_FLAM * scale - occ0

    mdiff = full - q
    rec_part = (1.0 - fbb) * E_rec
    bb_part = fbb * E_bb

    # Paschen share of the recombination part: re-evaluate the emission
    # shape with the Paschen ratio zeroed (module constant is hardcoded;
    # compute the shape ratio directly from the disk's own methods)
    p = flare_rec.pillars[0]
    te, vs = p['balmer_te'], p['balmer_vsmear']
    shape_full = np.array([flare_rec._arm_emission_shape(w, te, 0.0, vs)
                           for w in WL])
    import transient.transient_disk as td
    from pillar_disk import H as HP, K as KB, ANGSTROM as ANG, C as CC
    kte = KB * te

    def shape_no_paschen(w):
        e_photon = HP * CC / (w * ANG)
        sig = e_photon * vs * 1e5 / CC
        e_b = HP * CC / (3646.0 * ANG)
        e_p = HP * CC / (8204.0 * ANG)
        s_edge = 1.0 + 0.12 * np.exp(-(e_b - e_p) / kte)
        s = td.TransientPillarDisk._smeared_edge(e_photon, e_b, kte, sig)
        return s / s_edge

    shape_nop = np.array([shape_no_paschen(w) for w in WL])
    with np.errstate(divide='ignore', invalid='ignore'):
        pas_frac_of_rec = np.where(shape_full > 0,
                                   1.0 - shape_nop / shape_full, 0.0)
    pas_part = rec_part * pas_frac_of_rec

    chi2 = np.sum(((mdiff[good] - ddiff[good]) / derr[good]) ** 2)
    print(f"config: i={CFG['observation']['inclination']:.0f} deg, "
          f"T_e={te/1e3:.1f} kK, f_bb={fbb:.2f}, vsmear={vs:.0f} km/s, "
          f"T_pillar={CFG['transient']['pillar']['pillar_temp']:.0f} K")
    print(f"chi2 (line-free 2450-6000) = {chi2:.0f}  "
          f"(ndata = {good.sum()})")

    red = WL > 4500
    chi2_red = np.sum(((mdiff[good & red] - ddiff[good & red])
                       / derr[good & red]) ** 2)
    print(f"chi2 redward of 4500 A     = {chi2_red:.0f}  "
          f"(ndata = {(good & red).sum()})\n")

    hdr = (f"{'wl[A]':>7} {'data':>7} {'model':>7} {'rec-Bal':>8} "
           f"{'rec-Pas':>8} {'therm-BB':>9} {'-occ':>7}")
    print(hdr)
    for wtest in (3800., 4200., 4600., 5000., 5400., 5800.):
        i = np.argmin(np.abs(WL - wtest))
        print(f"{WL[i]:>7.0f} {ddiff[i]:>7.1f} {mdiff[i]:>7.1f} "
              f"{rec_part[i] - pas_part[i]:>8.1f} {pas_part[i]:>8.1f} "
              f"{bb_part[i]:>9.1f} {-occ_def[i]:>7.1f}")


if __name__ == '__main__':
    main()
