#!/usr/bin/env python3
"""extract_cloudy_ems.py - parse the depth-resolved line emissivities
('save lines emissivity') of the arm_column2 Cloudy grid into cumulative
emission-fraction curves versus hydrogen column density.

For each of the 196 grid models and each mapped line (H alpha, the C IV
and Mg II doublets summed), the curve F(N) is the fraction of that
line's total slab emission produced within column N of the illuminated
face (N = n_H x depth; the slabs are constant density). The curves are
resampled onto a common log N grid and cached as
data/cloudy_arm_ems_arm_column2.npz together with the grid parameters,
so the maps can place each line where Cloudy actually forms it.

Run:  python3 transient/extract_cloudy_ems.py
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from transient.extract_cloudy_arm import read_grid_params  # noqa: E402

BASE = '/Users/jiamuh/c23.01/my_models/arm_column2'
OUT = os.path.join(HERE, 'data', 'cloudy_arm_ems_arm_column2.npz')
XI = np.linspace(0.0, 1.0, 201)          # cumulative recombination frac
# .ems columns after depth: Halpha, CIV 1548, CIV 1551, MgII 2796, MgII 2803
COMBINE = {'halpha': [0], 'civ': [1, 2], 'mgii': [3, 4]}


def main():
    phi, hden, colden = read_grid_params(os.path.join(BASE,
                                                      'arm_column2.grd'))
    blocks = []
    cur = []
    with open(os.path.join(BASE, 'arm_column2.ems')) as fh:
        for line in fh:
            if 'GRID_DELIMIT' in line:
                blocks.append(np.array(cur))
                cur = []
            elif line.startswith('#'):
                continue
            else:
                cur.append([float(x) for x in line.split()])
    if cur:
        blocks.append(np.array(cur))
    assert len(blocks) == len(phi), (len(blocks), len(phi))

    # each line's cumulative emission fraction G(xi) as a function of the
    # slab's own cumulative RECOMBINATION fraction xi (traced by the
    # cumulative H alpha emission, a recombination line). The along-ray
    # coordinate matched to xi is the cumulative absorbed-photon fraction
    # of the Strommgren march, so H alpha reduces exactly to the n^2
    # recombination weight and the other lines shift relative to it as
    # Cloudy dictates.
    curves = {k: np.zeros((len(blocks), XI.size)) for k in COMBINE}
    for i, b in enumerate(blocks):
        depth = b[:, 0]
        cums = {}
        for key, cols in COMBINE.items():
            em = b[:, 1:][:, cols].sum(axis=1)
            cum = np.concatenate([[0.0], np.cumsum(
                0.5 * (em[1:] + em[:-1]) * np.diff(depth))])
            cum += em[0] * depth[0]
            cums[key] = cum
        xi_native = cums['halpha']
        if xi_native[-1] <= 0:
            for key in COMBINE:
                curves[key][i] = XI
            continue
        xi_native = xi_native / xi_native[-1]
        xi_native = np.maximum.accumulate(xi_native)
        for key, cum in cums.items():
            tot = cum[-1]
            if tot <= 0:
                curves[key][i] = XI
                continue
            curves[key][i] = np.interp(XI, xi_native, cum / tot)
    np.savez_compressed(OUT, xi=XI, phi=phi, hden=hden, colden=colden,
                        **{f'G_{k}': v for k, v in curves.items()})
    print(f"saved {OUT}  ({len(blocks)} models, G(xi) with "
          f"{XI.size} xi points)")


if __name__ == '__main__':
    main()
