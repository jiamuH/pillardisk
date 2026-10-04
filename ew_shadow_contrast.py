"""How much do the LOC equivalent widths change inside a pillar shadow?

Compares partially transparent pillars (f_trans = 0.3) with opaque pillars,
using the same LOC grid, spline and EW convention as test_loc_emissivity.py
(EW = F_line * 1215 A / incident nu F_nu at 1215 A, natural cubic spline in
log space).

Shadow depths used, measured with the paper's own shadow and flux code
(debug_shadow_rim_ftrans.py):
  f_trans = 0.3 : log Phi_H lower by 0.52 dex in the lane, 1.05 dex just
                  behind the pillar
  opaque        : flux clipped to the code's floor, log Phi_H = 17, or to the
                  paper's stated floor, log Phi_H = 15

Run:  python3 ew_shadow_contrast.py
"""
import numpy as np
from scipy.interpolate import CubicSpline

from pillardisk.pillar_line_cloudy import load_cloudy_models

CLOUDY_FILE = '/Users/jiamuh/c23.01/my_models/loc_metal_flux/strong_LOC_varym_N25_v100_lineflux_LineList_BLR_Fe2_flux.txt'
EXTENSION_FILE = '/Users/jiamuh/c23.01/my_models/loc_metal_flux/strong_LOC_varym_N25_v100_lineflux_extlow_LineList_BLR_Fe2_flux.txt'
LINES = [('Halpha', 'H-alpha'), ('Mg2', 'Mg II'), ('C4', 'C IV')]
LAM_REF = 1215.0
Z_TARGET = 1.0

interp_dict, phi_grid, _ = load_cloudy_models(
    CLOUDY_FILE, Z_target=Z_TARGET, extension_file=EXTENSION_FILE)


def grid_values(key):
    return np.array([float(interp_dict[key](p, Z_TARGET, grid=False)) for p in phi_grid])


inci = grid_values('inci1215')
log_ew = {}
for key, _ in LINES:
    ew = np.where(inci > 0, grid_values(key) * LAM_REF / inci, np.nan)
    log_ew[key] = CubicSpline(phi_grid, np.log10(np.maximum(ew, 1e-300)),
                              bc_type='natural', extrapolate=False)

print('LOC grid spans log Phi_H = %.1f to %.1f\n' % (phi_grid[0], phi_grid[-1]))

for lit in (21.0, 20.0, 19.0):
    cases = [('f_trans = 0.3, lane (-0.52 dex)', lit - 0.52),
             ('f_trans = 0.3, behind pillar (-1.05 dex)', lit - 1.05),
             ('Section 6.2 text, log Phi_H = 18.5', 18.5),
             ('opaque, code floor log Phi_H = 17', 17.0),
             ('opaque, paper floor log Phi_H = 15', 15.0)]
    print('Unshadowed log Phi_H = %.1f' % lit)
    header = '  %-42s' % 'shadowed case (log Phi_H in shadow)'
    print(header + ''.join('%22s' % ('change in log EW, ' + name) for _, name in LINES))
    for label, shadow in cases:
        row = '  %-34s (%5.2f)' % (label, shadow)
        for key, _ in LINES:
            row += '%22.2f' % (log_ew[key](shadow) - log_ew[key](lit))
        print(row)
    print()
