"""Referee reply, Section 4 bullet 7: where in the cloud do the lines form?

Reads the grid's saved line emissivity file (emissivity per zone against depth)
for two grid points at solar metallicity and log n_H = 10:
  index 870 = log Phi_H 19   (illuminated)
  index 261 = log Phi_H 17.5 (shadowed)
and reports the fraction of each line's emission that has accumulated by a
given hydrogen column density N_H = n_H * depth.

Run:  python3 analyze_line_depth.py
"""
import numpy as np

GRID_DIR = '/Users/jiamuh/c23.01/my_models/loc_metal_flux/'
STEM = 'strong_LOC_varym_N25_v100_lineflux'
MODELS = {870: 'illuminated, log Phi_H = 19', 261: 'shadowed, log Phi_H = 17.5'}
DENSITY = 1e10  # cm^-3, log n_H = 10 for both grid points
THRESHOLDS = [22., 23., 24., 24.5, 25.]
LABELS = {'H  1 4861.32A': 'H-beta', 'blnd 1549.00A': 'C IV',
          'blnd 2798.00A': 'Mg II', 'blnd 1909.00A': 'C III]'}


def read_blocks(wanted):
    """yield (index, header, array) for the requested model indices."""
    out, rows, header = {}, [], None
    index = 0
    with open(GRID_DIR + STEM + 'blrrnfb.ems', errors='replace') as f:
        for line in f:
            if line.startswith('#depth'):
                header = line.rstrip('\n').split('\t')
            elif 'GRID_DELIMIT' in line:
                if index in wanted:
                    out[index] = np.array(rows, dtype=float)
                rows = []
                index += 1
            elif line.strip():
                rows.append(line.rstrip('\n').split('\t'))
    return header, out


header, blocks = read_blocks(set(MODELS))
for index, description in MODELS.items():
    a = blocks[index]
    depth = a[:, 0]
    column = DENSITY * depth
    print('\nmodel %d (%s): %d zones, final log N_H = %.2f'
          % (index, description, len(depth), np.log10(column[-1])))
    print('%-8s' % 'line' + ''.join('%18s' % ('by log N_H = %.1f' % t) for t in THRESHOLDS))
    widths = np.diff(np.concatenate([[0.], depth]))
    for j, name in enumerate(header[1:], start=1):
        if name.strip() not in LABELS:
            continue
        cumulative = np.cumsum(a[:, j] * widths)
        fraction = cumulative / cumulative[-1]
        row = '%-8s' % LABELS[name.strip()]
        for t in THRESHOLDS:
            k = np.searchsorted(column, 10. ** t)
            row += '%18s' % ('%.3f' % (fraction[min(k, len(fraction) - 1)]
                                       if k < len(fraction) else 1.0))
        print(row)
