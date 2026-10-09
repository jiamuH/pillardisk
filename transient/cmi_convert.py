#!/usr/bin/env python3
"""cmi_convert.py - convert an Athena++ disk dump into a CMacIonize
photoionization test (minimal, hydrogen-only scope).

Resamples the physical hydrogen density of the simulation (spherical
r, theta, phi grid) onto a Cartesian grid by volume-averaging over
SUB^3 sub-points per Cartesian cell, and writes:

  <outdir>/density.txt     CMacIonize ASCII density input
                           (x y z n_H per line; x, y, z in units of r0,
                           n_H in cm^-3)
  <outdir>/test.param      CMacIonize parameter file: hydrogen only,
                           fixed temperature 1e4 K, monochromatic
                           13.6 eV point lamp at the origin with the
                           pipeline's Q_ION, diffuse field ON
  <outdir>/grid.npz        the same Cartesian density as a numpy array,
                           for running our own march on the identical
                           grid (so a comparison isolates physics from
                           resolution)

Cells outside the simulated volume (the central hole r < r_in, beyond
r_out, and above the wedge) get a negligible density N_EMPTY.

Run:  python3 -m transient.cmi_convert
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

N_EMPTY = 1.0          # cm^-3, density outside the simulated volume
SUB = 3                # sub-points per Cartesian cell per dimension


def resample(sim, phys, ncell, half, empty=N_EMPTY):
    """Volume-averaged n_H on a Cartesian grid spanning [-half, half]
    in each direction (units of r0). Returns (centres, n_H grid)."""
    rf, thf, phf = sim['rf'], sim['thf'], sim['phf']
    nH = phys['nH']                                    # (ph, th, r)
    nx, ny, nz = ncell
    edges = [np.linspace(-half[k], half[k], ncell[k] + 1) for k in range(3)]
    cent = [0.5 * (e[1:] + e[:-1]) for e in edges]
    # sub-point offsets within a cell, as fractions of the cell size
    off = (np.arange(SUB) + 0.5) / SUB - 0.5
    grid = np.zeros((nx, ny, nz))
    dz = edges[2][1] - edges[2][0]
    dx = edges[0][1] - edges[0][0]
    dy = edges[1][1] - edges[1][0]
    X, Y = np.meshgrid(cent[0], cent[1], indexing='ij')
    for k in range(nz):
        acc = np.zeros((nx, ny))
        for ox in off:
            for oy in off:
                for oz in off:
                    x = X + ox * dx
                    y = Y + oy * dy
                    z = cent[2][k] + oz * dz
                    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
                    th = np.arccos(np.clip(z / np.maximum(r, 1e-12), -1, 1))
                    ph = np.mod(np.arctan2(y, x), 2 * np.pi)
                    inside = ((r >= rf[0]) & (r < rf[-1])
                              & (th >= thf[0]) & (th < thf[-1]))
                    ir = np.clip(np.searchsorted(rf, r) - 1, 0, len(rf) - 2)
                    it = np.clip(np.searchsorted(thf, th) - 1, 0,
                                 len(thf) - 2)
                    ip = np.clip(np.searchsorted(phf, ph) - 1, 0,
                                 len(phf) - 2)
                    acc += np.where(inside, nH[ip, it, ir], empty)
        grid[:, :, k] = acc / SUB ** 3
    return cent, grid


PARAM = """# CMacIonize minimal test: hydrogen-only photoionization of the
# Athena++ disk dump {dump}, written by transient/cmi_convert.py
SimulationBox:
  anchor: [{ax:.6e} m, {ay:.6e} m, {az:.6e} m]
  sides: [{sx:.6e} m, {sy:.6e} m, {sz:.6e} m]
  periodicity: [false, false, false]

DensityGrid:
  type: Cartesian
  number of cells: [{nx}, {ny}, {nz}]

DensityFunction:
  type: AsciiFile
  filename: density.txt
  number of cells: [{nx}, {ny}, {nz}]
  box anchor: [{ax:.6e} m, {ay:.6e} m, {az:.6e} m]
  box sides: [{sx:.6e} m, {sy:.6e} m, {sz:.6e} m]
  temperature: 10000. K
  length unit: {lu:.6e} m
  density unit: 1. cm^-3

TemperatureCalculator:
  do temperature calculation: false

PhotonSourceDistribution:
  type: SingleStar
  position: [0. m, 0. m, 0. m]
  luminosity: {q:.4e} s^-1

PhotonSourceSpectrum:
  type: Monochromatic
  frequency: 13.6 eV

TaskBasedIonizationSimulation:
  number of iterations: {niter}
  number of photons: {nphot:d}
  diffuse field: true

IonizationSimulation:
  number of iterations: {niter}
  number of photons: {nphot:d}

DensityGridWriter:
  type: Gadget
  prefix: disk_
  padding: 3

DiffuseReemissionHandler:
  type: Physical

CrossSections:
  type: FixedValue
  hydrogen_0: 6.3e-18 cm^2
  helium_0: 0. m^2
  carbon_1: 0. m^2
  carbon_2: 0. m^2
  nitrogen_0: 0. m^2
  nitrogen_1: 0. m^2
  nitrogen_2: 0. m^2
  oxygen_0: 0. m^2
  oxygen_1: 0. m^2
  neon_0: 0. m^2
  neon_1: 0. m^2
  sulphur_1: 0. m^2
  sulphur_2: 0. m^2
  sulphur_3: 0. m^2

# case A total recombination rate at 1e4 K: the diffuse field is
# transported explicitly, so on-the-spot (case B) must NOT be assumed
RecombinationRates:
  type: FixedValue
  hydrogen_1: 4.18e-13 cm^3 s^-1
  helium_1: 0. m^3 s^-1
  carbon_2: 0. m^3 s^-1
  carbon_3: 0. m^3 s^-1
  nitrogen_1: 0. m^3 s^-1
  nitrogen_2: 0. m^3 s^-1
  nitrogen_3: 0. m^3 s^-1
  oxygen_1: 0. m^3 s^-1
  oxygen_2: 0. m^3 s^-1
  neon_1: 0. m^3 s^-1
  neon_2: 0. m^3 s^-1
  sulphur_2: 0. m^3 s^-1
  sulphur_3: 0. m^3 s^-1
  sulphur_4: 0. m^3 s^-1
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dump', default=os.path.join(
        HERE, 'data', 'sim', 'disk.out1.00012.athdf'))
    ap.add_argument('--outdir', default=os.path.join(HERE, 'data',
                                                     'cmi_test'))
    ap.add_argument('--ncell', type=int, nargs=3, default=[128, 128, 64])
    ap.add_argument('--nphot', type=float, default=1e7)
    ap.add_argument('--niter', type=int, default=20)
    a = ap.parse_args()

    sim = load_sim(a.dump)
    phys = to_physical(sim, CONFIG)
    rmax = sim['rf'][-1]
    zmax = rmax * np.cos(sim['thf'][0])               # top of the wedge
    half = [rmax * 1.0, rmax * 1.0, zmax * 1.02]
    cent, grid = resample(sim, phys, a.ncell, half)
    os.makedirs(a.outdir, exist_ok=True)

    X, Y, Z = np.meshgrid(*cent, indexing='ij')
    np.savetxt(os.path.join(a.outdir, 'density.txt'),
               np.column_stack([X.ravel(), Y.ravel(), Z.ravel(),
                                grid.ravel()]),
               fmt='%.6e', header='x y z [r0]  n_H [cm^-3]')
    np.savez_compressed(os.path.join(a.outdir, 'grid.npz'),
                        x=cent[0], y=cent[1], z=cent[2], nH=grid,
                        half=half, q_ion=Q_ION)

    r0_m = CONFIG['R0_LD'] * LD_CM / 100.0
    with open(os.path.join(a.outdir, 'test.param'), 'w') as fh:
        fh.write(PARAM.format(
            dump=os.path.basename(a.dump),
            ax=-half[0] * r0_m, ay=-half[1] * r0_m, az=-half[2] * r0_m,
            sx=2 * half[0] * r0_m, sy=2 * half[1] * r0_m,
            sz=2 * half[2] * r0_m,
            nx=a.ncell[0], ny=a.ncell[1], nz=a.ncell[2],
            lu=r0_m, q=Q_ION, niter=a.niter, nphot=int(a.nphot)))

    cell = [2 * half[k] / a.ncell[k] for k in range(3)]
    print(f"wrote {a.outdir}: {np.prod(a.ncell)} cells, cell size "
          f"{cell[0]:.4f} x {cell[1]:.4f} x {cell[2]:.4f} r0, "
          f"n_H range {grid.min():.2e} - {grid.max():.2e} cm^-3")


if __name__ == '__main__':
    main()
