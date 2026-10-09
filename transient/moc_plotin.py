#!/usr/bin/env python3
"""moc_plotin.py - choose the lines mocassinPlot writes per cell, from a
finished run's volume-integrated output/lineFlux.out.

mocassinPlot (MOCASSIN's per-cell emissivity driver) only writes the lines
listed in input/plot.in, by MOCASSIN's internal line code. This script
selects every line with lambda in [--wmin, --wmax] and flux >= --min-ratio
x H beta from the "Formal Solution" block of lineFlux.out (which lists each
line with its code):
  H I    (upper, lower) -> vacuum lambda from the Rydberg formula
  He I   34 recombination lines, lambda from the "Analytical" block (rounded
         air wavelengths, +-2 A)
  He II  (upper, lower) -> vacuum lambda from the Rydberg formula (Z = 2)
  metals (element, ion, lower, upper) with MOCASSIN's vacuum lambda
The seven lines of the original plot.in come first under their old names
(Ha, Hb, CIV1551, CIV1548, MgII2804, MgII2796, OIII5007), so the map panels
are unchanged; the rest follow, strongest first.

Writes:
  <run>/input/plot.in        mocassinPlot input ('mono', one line per row)
  <run>/input/plot_lines.txt column table for plot.out: name, code,
                             lambda [A], air (1 = lambda is in air),
                             ion mass [amu], ratio to H beta
Then run mocassinPlot in <run> (overwrites output/plot.out).

Run:  python3 -m transient.moc_plotin [--run DIR] [--min-ratio 0.01]
"""

import argparse
import os

import numpy as np

DEFAULT_RUN = '/data2/jhuang/runs/mocassin/moc_coarse'
R_H, R_HE2 = 109677.58, 4 * 109722.27              # cm^-1
MASS = {1: 1.008, 2: 4.003, 6: 12.011, 7: 14.007, 8: 15.999, 10: 20.180,
        12: 24.305, 14: 28.086, 16: 32.06}
SYM = {1: 'H', 2: 'He', 6: 'C', 7: 'N', 8: 'O', 10: 'Ne', 12: 'Mg',
       14: 'Si', 16: 'S'}
ROMAN = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII']
LEGACY = [('Ha', 1), ('Hb', 2), ('CIV1551', 1393), ('CIV1548', 1394),
          ('MgII2804', 6884), ('MgII2796', 6885), ('OIII5007', 3741)]


def parse_lineflux(path):
    """All lines of the first component: list of
    (code, lambda [A], air, mass, ratio, name)."""
    txt = open(path).read().splitlines()
    i0 = next(i for i, s in enumerate(txt) if 'Formal Solution' in s)
    i1 = next(i for i, s in enumerate(txt) if 'Component:' in s)
    # He I wavelengths: the Analytical block lists the 34 lines in code order
    j = next(i for i, s in enumerate(txt) if s.strip() == 'HeI')
    he1_lam = [float(s.split()[0]) for s in txt[j + 1:j + 35]]
    out, block = [], None
    for s in txt[i0:i1]:
        t = s.split()
        if 'HI recombination' in s:
            block = 'HI'
        elif s.strip() == 'HeI':
            block = 'HeI'
        elif 'HeII recombination' in s:
            block = 'HeII'
        elif 'forbidden' in s:
            block = 'metal'
        try:
            nums = [float(x) for x in t]
        except ValueError:
            continue
        if block == 'HI' and len(nums) == 4:
            up, lo, ratio, code = nums
            lam = 1e8 / (R_H * (1 / lo ** 2 - 1 / up ** 2))
            out.append((int(code), lam, 0, MASS[1], ratio,
                        f'HI_{int(up)}-{int(lo)}'))
        elif block == 'HeI' and len(nums) == 3:
            l, ratio, code = nums
            out.append((int(code), he1_lam[int(l) - 1], 1, MASS[2], ratio,
                        f'HeI_{he1_lam[int(l) - 1]:.0f}'))
        elif block == 'HeII' and len(nums) == 4:
            up, lo, ratio, code = nums
            lam = 1e8 / (R_HE2 * (1 / lo ** 2 - 1 / up ** 2))
            out.append((int(code), lam, 0, MASS[2], ratio,
                        f'HeII_{int(up)}-{int(lo)}'))
        elif block == 'metal' and len(nums) == 7:
            el, ion, lo, up, lam, ratio, code = nums
            el, ion = int(el), int(ion)
            out.append((int(code), lam, 0, MASS[el], ratio,
                        f'{SYM[el]}{ROMAN[ion - 1]}_{lam:.0f}'))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', default=DEFAULT_RUN)
    ap.add_argument('--min-ratio', type=float, default=0.01)
    ap.add_argument('--wmin', type=float, default=1000.0)
    ap.add_argument('--wmax', type=float, default=11000.0)
    a = ap.parse_args()
    lines = parse_lineflux(os.path.join(a.run, 'output', 'lineFlux.out'))
    by_code = {c[0]: c for c in lines}
    rows = [(name,) + by_code[code][:5] for name, code in LEGACY]
    keep = sorted((c for c in lines
                   if a.wmin <= c[1] <= a.wmax and c[4] >= a.min_ratio
                   and c[0] not in dict(LEGACY).values()),
                  key=lambda c: -c[4])
    names = {r[0] for r in rows}
    for c in keep:
        name = c[5]
        while name in names:                    # same ion, same rounded lambda
            name += "'"
        names.add(name)
        rows.append((name,) + c[:5])
    inp = os.path.join(a.run, 'input')
    with open(os.path.join(inp, 'plot.in'), 'w') as fh:
        fh.write('mono\n')
        for r in rows:
            fh.write(f"line {r[1]:<8d} {r[2]:.2f}   {r[2]:.2f}\n")
    with open(os.path.join(inp, 'plot_lines.txt'), 'w') as fh:
        fh.write('# name code lambda_A air mass_amu ratio_Hbeta '
                 '(column order of output/plot.out)\n')
        for r in rows:
            fh.write(f"{r[0]:14s} {r[1]:6d} {r[2]:12.4f} {r[3]:d} "
                     f"{r[4]:8.3f} {r[5]:.4e}\n")
    kinds = {}
    for r in rows:
        k = r[0].split('_')[0] if '_' in r[0] else 'legacy'
        kinds[k] = kinds.get(k, 0) + 1
    print(f"{len(rows)} lines ({a.wmin:.0f}-{a.wmax:.0f} A, >= "
          f"{a.min_ratio:g} H beta; 7 legacy first): {kinds}")
    print(f"wrote {inp}/plot.in and plot_lines.txt")


if __name__ == '__main__':
    main()
