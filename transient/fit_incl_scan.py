#!/usr/bin/env python3
"""
fit_incl_scan.py - inclination scan of the tuned analytic-edge arm model
(the current best continuum configuration). For each inclination the
geometry (occultation deficit, emission prefactor) is rebuilt and a small
(T_e, f_bb, vsmear) grid is refitted with the linear amplitude solve.

Run:  python3 transient/fit_incl_scan.py
"""

import copy
import itertools
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import transient.tune_redside as tr  # noqa: E402

I_GRID = [30., 38., 46., 58., 70., 80.]
TE_GRID = [7000., 8000., 9000., 10000.]
FBB_GRID = [0.0, 0.025, 0.05]
VS_GRID = [45000., 55000., 65000.]


def main():
    print(f"  {'i[deg]':>7} {'chi2/dof':>9} {'chi2_red':>9} "
          f"{'occ(2550)':>10} {'T_e[kK]':>8} {'f_bb':>6} {'v[1e3]':>7} "
          f"{'amp':>6}")
    for inc in I_GRID:
        cfg = copy.deepcopy(tr.CFG)
        cfg['observation']['inclination'] = inc
        WL, ddiff, derr, good, occ_def, A, flare_ref = tr.setup(cfg)
        w_fit = 1.0 / derr[good] ** 2
        red = WL > 4500
        occ_uv = float(np.interp(2550.0, WL, occ_def))
        best = None
        for te, fbb, vs in itertools.product(TE_GRID, FBB_GRID, VS_GRID):
            shape = np.array([flare_ref._arm_emission_shape(
                w, te, fbb, vs, 0.0) for w in WL])
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
            if best is None or chi2 < best[0]:
                best = (chi2, chi2_red, te, fbb, vs, s)
        ndof = good.sum() - 6
        chi2, chi2_red, te, fbb, vs, s = best
        print(f"  {inc:>7.0f} {chi2/ndof:>9.1f} {chi2_red:>9.0f} "
              f"{occ_uv:>10.1f} {te/1e3:>8.0f} {fbb:>6.3f} "
              f"{vs/1e3:>7.0f} {s:>6.2f}", flush=True)


if __name__ == '__main__':
    main()
