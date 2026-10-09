#!/usr/bin/env python3
"""
extract_cloudy_arm.py - one-time extraction of the emitted (diffuse)
spectra from a Cloudy grid's `save continuum units Angstroms` output into
a compact npz cache for the transient-arm modeling.

Supported grids (argument = grid name, default loc_metal_flux):
  loc_metal_flux : the 1827-point BLR grid (phi, hden, Z)
  arm_column     : the 100-point matter-bounded arm grid (phi, hden, colden)

Output npz (transient/data/cloudy_arm_spectra_<grid>.npz):
    wave      (nw,)        wavelength [A], common to all blocks
    refl_cont (npts, nw)   reflected (illuminated-face) continuum
    refl_line (npts, nw)   reflected lines
    out_cont  (npts, nw)   outward (shielded-face) diffuse continuum
    out_line  (npts, nw)   outward lines
    phi, hden (npts,)      grid parameters per block
    Z or colden (npts,)    third grid parameter (name per grid)
All spectra are nu*F_nu [erg s^-1 cm^-2] at the cloud surface.

Run:  python3 transient/extract_cloudy_arm.py [loc_metal_flux|arm_column]
"""

import os
import sys

import numpy as np

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
GRIDS = {
    'loc_metal_flux': dict(
        base='/Users/jiamuh/c23.01/my_models/loc_metal_flux',
        prefix='strong_LOC_varym_N25_v100_lineflux',
        con_suffix='_SED.conA',
        par3='Z'),
    'arm_column': dict(
        base='/Users/jiamuh/c23.01/my_models/arm_column',
        prefix='arm_column',
        con_suffix='_SED.conA',
        par3='colden'),
    'arm_column2': dict(
        base='/Users/jiamuh/c23.01/my_models/arm_column2',
        prefix='arm_column2',
        con_suffix='_SED.conA',
        par3='colden'),
}
WMIN, WMAX = 300.0, 12000.0
# .conA columns: 0 wave, 1 incident, 2 trans, 3 DiffOut, 4 net trans,
#                5 reflc, 6 total, 7 reflin, 8 outlin, 9 lineID, 10 cont,
#                11 nLine


def cache_path(grid_name):
    tag = '' if grid_name == 'loc_metal_flux' else f'_{grid_name}'
    return os.path.join(DATA_DIR, f'cloudy_arm_spectra{tag}.npz')


def read_grid_params(grd_file):
    """Three varied grid parameters per block, in vary-command order
    (columns 6-8 of the .grd file)."""
    p1, p2, p3 = [], [], []
    with open(grd_file) as fh:
        header = fh.readline()
        assert header.startswith('#Index')
        for line in fh:
            f = line.split('\t')
            p1.append(float(f[6]))
            p2.append(float(f[7]))
            p3.append(float(f[8]))
    return np.array(p1), np.array(p2), np.array(p3)


def main():
    grid_name = sys.argv[1] if len(sys.argv) > 1 else 'loc_metal_flux'
    g = GRIDS[grid_name]
    con = os.path.join(g['base'], g['prefix'] + g['con_suffix'])
    grd = os.path.join(g['base'], g['prefix'] + '.grd')

    phi, hden, p3 = read_grid_params(grd)
    npts = len(phi)
    print(f"grid '{grid_name}': {npts} points "
          f"(phi {phi.min():.1f}-{phi.max():.1f}, "
          f"hden {hden.min():.1f}-{hden.max():.1f}, "
          f"{g['par3']} {p3.min():.1f}-{p3.max():.1f})")

    blocks = []          # (nrow, 5) arrays [wave, DiffOut, reflc,
    cur = []             #                   reflin, outlin]
    nline = 0
    with open(con) as fh:
        for line in fh:
            nline += 1
            if line.startswith('#'):
                if 'GRID_DELIMIT' in line:
                    blocks.append(np.array(cur))
                    cur = []
                    if len(blocks) % 200 == 0:
                        print(f"  {len(blocks)} blocks read "
                              f"({nline/1e6:.0f}M lines)", flush=True)
                continue
            f = line.split('\t')
            w = float(f[0])
            if WMIN <= w <= WMAX:
                cur.append((w, float(f[3]), float(f[5]),
                            float(f[7]), float(f[8])))
    if cur:
        blocks.append(np.array(cur))
    print(f"read {len(blocks)} blocks total")
    assert len(blocks) == npts, f"{len(blocks)} blocks vs {npts} grid points"

    nw = len(blocks[0])
    wave = blocks[0][:, 0][::-1]                # ascending wavelength
    for b in blocks:
        assert len(b) == nw and np.allclose(b[:, 0][::-1], wave), \
            "wavelength mesh differs between blocks"

    out_cont = np.array([b[:, 1][::-1] for b in blocks])
    refl_cont = np.array([b[:, 2][::-1] for b in blocks])
    refl_line = np.array([b[:, 3][::-1] for b in blocks])
    out_line = np.array([b[:, 4][::-1] for b in blocks])

    out = cache_path(grid_name)
    os.makedirs(DATA_DIR, exist_ok=True)
    np.savez_compressed(out, wave=wave, refl_cont=refl_cont,
                        refl_line=refl_line, out_cont=out_cont,
                        out_line=out_line, phi=phi, hden=hden,
                        **{g['par3']: p3})
    print(f"saved {out}  (wave: {nw} points, "
          f"{wave.min():.0f}-{wave.max():.0f} A)")


if __name__ == '__main__':
    main()
