#!/usr/bin/env python3
"""moc_convert.py - convert an Athena++ disk dump into a MOCASSIN
photoionization run (full physics: the pipeline's AGN spectrum including
X-rays, helium and metals, self-consistent temperature).

Writes a self-contained run directory (default transient/data/moc_test):

  input/input.in      MOCASSIN input file
  input/density.dat   x y z [cm]  n_H [cm^-3], loop order x, y, z (z
                      fastest), as MOCASSIN reads it
  input/agn_sed.dat   the Cloudy-deck incident AGN spectrum (the same
                      interpolate table as the pipeline), as wavelength
                      [A] vs f_lambda (shape only: MOCASSIN rescales it
                      to LPhot = the pipeline's Q_ION)
  input/abun.in       solar abundances (Asplund et al. 2009) for H, He,
                      C, N, O, Ne, Mg, Si, S (and Fe with --with-iron);
                      all others zero
  grid.npz            the same Cartesian density as a numpy array
  data, dustData      symbolic links to MOCASSIN's atomic data
  output/             empty, for MOCASSIN's results

The grid has an ODD number of cells per axis, so one cell is centred on
the lamp at the origin. Cells outside the simulated volume get zero
density (MOCASSIN treats them as empty).

Run:  python3 -m transient.moc_convert
"""

import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from transient.sim_spectrum import (  # noqa: E402
    load_sim, to_physical, CONFIG, LD_CM)
from transient.sim_cloudy_spectrum import Q_ION  # noqa: E402
from transient.sim_full_spectrum import SED_X, SED_Y  # noqa: E402
from transient.cmi_convert import resample  # noqa: E402

MOC = '/Users/jiamuh/codes/MOCASSIN-2.0'

# abundances by number relative to H (Asplund et al. 2009, solar);
# elements not listed are switched off (zero) to save memory and time
ABUN = {1: 1.0, 2: 0.0851, 6: 2.69e-4, 7: 6.76e-5, 8: 4.90e-4,
        10: 8.51e-5, 12: 3.98e-5, 14: 3.24e-5, 16: 1.32e-5, 26: 3.16e-5}
SYMBOLS = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na',
           'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti',
           'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']

INPUT = """densityFile "input/density.dat"
nx {nx}
ny {ny}
nz {nz}
contShape "input/agn_sed.dat"
TStellar 100000.
LPhot {lphot:.4e}
nebComposition "input/abun.in"
TeStart 10000.
nuMin 1.001e-5
nuMax {numax:.1f}
nbins {nbins}
nPhotons {nphot:d}
autoPackets 0.20 2. {nphotmax:d}
maxIterateMC {niter} 95.
convLimit 0.05
nstages 7
Rin {rin:.6e}
output
"""


def write_sed(path, numax):
    """Cloudy-deck AGN SED (log E [Ryd] vs log f_nu) as wavelength [A]
    vs f_lambda, ascending wavelength, spanning 1e-5 Ryd to numax."""
    logE = np.linspace(-5.0, np.log10(numax), 3000)
    fnu = 10.0 ** np.interp(logE, SED_X, SED_Y)
    lam = 911.76 / 10.0 ** logE                         # [A]
    flam = fnu / lam ** 2                               # shape only
    order = np.argsort(lam)
    np.savetxt(path, np.column_stack([lam[order], flam[order]]),
               fmt='%.6e')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dump', default=os.path.join(
        HERE, 'data', 'sim', 'disk.out1.00012.athdf'))
    ap.add_argument('--outdir', default=os.path.join(HERE, 'data',
                                                     'moc_test'))
    ap.add_argument('--ncell', type=int, nargs=3, default=[65, 65, 33])
    ap.add_argument('--nphot', type=int, default=10000000)
    ap.add_argument('--niter', type=int, default=20)
    ap.add_argument('--numax', type=float, default=300.0)  # Ryd, ~4 keV
    ap.add_argument('--nbins', type=int, default=1500)
    ap.add_argument('--with-iron', action='store_true',
                    help='include iron (needs nForLevels raised in the '
                         'MOCASSIN source and much more memory)')
    a = ap.parse_args()
    abun = dict(ABUN)
    if not a.with_iron:
        abun.pop(26)
    if any(n % 2 == 0 for n in a.ncell):
        sys.exit("use odd cell counts, so a cell is centred on the lamp")

    sim = load_sim(a.dump)
    phys = to_physical(sim, CONFIG)
    rmax = sim['rf'][-1]
    zmax = rmax * np.cos(sim['thf'][0])
    half = [rmax, rmax, zmax * 1.02]
    cent, grid = resample(sim, phys, a.ncell, half, empty=0.0)

    for sub in ('input', 'output'):
        os.makedirs(os.path.join(a.outdir, sub), exist_ok=True)
    for link in ('data', 'dustData'):
        dst = os.path.join(a.outdir, link)
        if not os.path.lexists(dst):
            os.symlink(os.path.join(MOC, link), dst)

    r0_cm = CONFIG['R0_LD'] * LD_CM
    X, Y, Z = np.meshgrid(*cent, indexing='ij')         # z fastest
    np.savetxt(os.path.join(a.outdir, 'input', 'density.dat'),
               np.column_stack([X.ravel() * r0_cm, Y.ravel() * r0_cm,
                                Z.ravel() * r0_cm, grid.ravel()]),
               fmt='%.6e')
    np.savez_compressed(os.path.join(a.outdir, 'grid.npz'),
                        x=cent[0], y=cent[1], z=cent[2], nH=grid,
                        half=half, q_ion=Q_ION)
    write_sed(os.path.join(a.outdir, 'input', 'agn_sed.dat'), a.numax)
    with open(os.path.join(a.outdir, 'input', 'abun.in'), 'w') as fh:
        for zat, sym in enumerate(SYMBOLS, start=1):
            fh.write(f"{abun.get(zat, 0.0):.4e}   ! {sym}\n")
    with open(os.path.join(a.outdir, 'input', 'input.in'), 'w') as fh:
        fh.write(INPUT.format(nx=a.ncell[0], ny=a.ncell[1], nz=a.ncell[2],
                              lphot=Q_ION / 1e36, numax=a.numax,
                              nbins=a.nbins, nphot=a.nphot,
                              nphotmax=a.nphot * 8, niter=a.niter,
                              rin=0.48 * r0_cm))
    gas = grid > 0
    cell = [2 * half[k] / a.ncell[k] for k in range(3)]
    print(f"wrote {a.outdir}: {np.prod(a.ncell)} cells "
          f"({gas.sum()} with gas), cell size {cell[0]:.4f} x "
          f"{cell[1]:.4f} x {cell[2]:.4f} r0, n_H max {grid.max():.2e}")


if __name__ == '__main__':
    main()
