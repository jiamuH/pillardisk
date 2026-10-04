"""Referee reply, Section 4 bullet 7: where do the paper's Cloudy grid models
actually stop, and what hydrogen column do they reach?

Parses the grid .out file. Each model's output block ends with a line
'GRID_DELIMIT -- gridNNNNNNNNN' giving its index; the grid parameters for that
index come from the .grd table.

Run:  python3 analyze_grid_stopping.py
"""
import re
import numpy as np

GRID_DIR = '/Users/jiamuh/c23.01/my_models/loc_metal_flux/'
STEM = 'strong_LOC_varym_N25_v100_lineflux'

stop_pattern = re.compile(r'stopped because (.*?)\. Iteration')
column_pattern = re.compile(r'^ Hydrogen\s+([-\d.]+)')
index_pattern = re.compile(r'GRID_DELIMIT -- grid(\d+)')

stop_reason, hydrogen_column = {}, {}
current_stop, current_column = None, None
with open(GRID_DIR + STEM + '.out', errors='replace') as f:
    for line in f:
        m = stop_pattern.search(line)
        if m:
            current_stop = m.group(1)
        elif 'Log10 Column density' in line:
            c = column_pattern.match(line)
            if c:
                current_column = float(c.group(1))
        else:
            m = index_pattern.search(line)
            if m:
                i = int(m.group(1))
                stop_reason[i] = current_stop
                hydrogen_column[i] = current_column
                current_stop, current_column = None, None

grid = np.loadtxt(GRID_DIR + STEM + '.grd', skiprows=1, usecols=(6, 7, 8))
phi, hden, metals = grid[:, 0], grid[:, 1], grid[:, 2]
n = len(grid)
print('models in .grd: %d, parsed from .out: %d' % (n, len(stop_reason)))

reasons = np.array([stop_reason.get(i, 'missing') for i in range(n)])
columns = np.array([hydrogen_column.get(i, np.nan) for i in range(n)])

print('\nstopping reason, all models')
for r in sorted(set(reasons)):
    sel = reasons == r
    print('  %-28s %5d models (%4.1f per cent), log N_H reached: median %.2f, range %.2f to %.2f'
          % (r, sel.sum(), 100. * sel.sum() / n, np.nanmedian(columns[sel]),
             np.nanmin(columns[sel]), np.nanmax(columns[sel])))

print('\nstopping reason at solar metallicity only (Z = 1), by ionizing flux')
solar = metals == 1.
print('%8s %10s %10s %10s %14s' % ('log Phi_H', 'lowest Te', 'H column', 'zone limit', 'median log N_H'))
for p in np.unique(phi):
    sel = solar & (phi == p)
    counts = [np.sum(sel & (reasons == r)) for r in
              ('lowest Te reached', 'H column dens reached', 'NZONE reached')]
    print('%8.1f %10d %10d %10d %14.2f'
          % (p, counts[0], counts[1], counts[2], np.nanmedian(columns[sel])))

print('\nsolar metallicity, by gas density')
print('%8s %10s %10s %10s %14s' % ('log n_H', 'lowest Te', 'H column', 'zone limit', 'median log N_H'))
for d in np.unique(hden):
    sel = solar & (hden == d)
    counts = [np.sum(sel & (reasons == r)) for r in
              ('lowest Te reached', 'H column dens reached', 'NZONE reached')]
    print('%8.1f %10d %10d %10d %14.2f'
          % (d, counts[0], counts[1], counts[2], np.nanmedian(columns[sel])))
