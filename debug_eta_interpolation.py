"""Does the responsivity eta = dlog F / dlog Phi used in the velocity-delay
maps have kinks at the Cloudy grid knots?

Compares, for C IV, Mg II and H-alpha:
  (a) the maps' method: quadratic spline in linear flux (RectBivariateSpline,
      kx=2, from load_cloudy_models) with a symmetric finite difference of
      +-0.05 dex, as in compute_velocity_delay_map_cloudy;
  (b) the paper's stated method (and Fig. 3): natural cubic spline in log flux,
      differentiated analytically.

Run:  python3 -m pillardisk.debug_eta_interpolation
"""
import os

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

from pillardisk.pillar_line_cloudy import load_cloudy_models

CLOUDY = '/Users/jiamuh/c23.01/my_models/loc_metal_flux/strong_LOC_varym_N25_v100_lineflux_LineList_BLR_Fe2_flux.txt'
EXT = '/Users/jiamuh/c23.01/my_models/loc_metal_flux/strong_LOC_varym_N25_v100_lineflux_extlow_LineList_BLR_Fe2_flux.txt'
LINES = [('C4', 'C IV', 'royalblue'), ('Mg2', 'Mg II', 'forestgreen'),
         ('Halpha', 'H-alpha', 'orangered'), ('HI', 'H-beta', 'black')]
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'plots',
                   'debug_eta_interpolation.png')

interp, phi_grid, _ = load_cloudy_models(CLOUDY, Z_target=1.0, extension_file=EXT)
lp = np.linspace(phi_grid[0] + 0.05, phi_grid[-1] - 0.05, 2000)

plt.rcParams.update({'text.usetex': True, 'font.family': 'serif', 'font.size': 18})
fig, ax = plt.subplots(figsize=(9, 6))
for key, name, color in LINES:
    f = lambda x: np.array([float(interp[key](v, 1.0, grid=False)) for v in np.atleast_1d(x)])
    eta_maps = (np.log10(f(lp + 0.05)) - np.log10(f(lp - 0.05))) / 0.1
    knots = np.array([float(interp[key](p, 1.0, grid=False)) for p in phi_grid])
    eta_paper = CubicSpline(phi_grid, np.log10(knots), bc_type='natural').derivative()(lp)
    jump = np.max(np.abs(np.diff(eta_maps)))
    print('%-8s largest step in eta between neighbouring points (maps method): %.3f;'
          ' largest |eta_maps - eta_paper|: %.2f' % (name, jump, np.max(np.abs(eta_maps - eta_paper))))
    ax.plot(lp, eta_maps, color=color, lw=2, label=r'$\rm %s,~maps$' % name.replace(' ', '~'))
    ax.plot(lp, eta_paper, color=color, lw=2, ls='--')
for p in phi_grid:
    ax.axvline(p, color='gray', lw=0.5, alpha=0.5)
ax.set_xlabel(r'$\log\Phi_{\rm H}~[\rm cm^{-2}~s^{-1}]$')
ax.set_ylabel(r'$\eta$')
ax.set_ylim(-1, 3)
ax.minorticks_on()
ax.tick_params(which='both', direction='in', top=True, right=True)
ax.legend(fontsize=13, loc='upper right')
fig.tight_layout()
fig.savefig(OUT, dpi=150)
print('saved', OUT)

# Responsivity of each line relative to H-beta, log(eta_X / eta_Hbeta), from
# the paper's method (log-space cubic spline), over the disc's flux range.
def eta_of(key):
    knots = np.array([float(interp[key](p, 1.0, grid=False)) for p in phi_grid])
    return CubicSpline(phi_grid, np.log10(knots), bc_type='natural').derivative()

eta_hb = eta_of('HI')
for lphi in (16, 17, 18, 18.5, 19, 19.5, 20, 20.5, 21):
    row = '  log Phi_H = %4.1f  eta_Hbeta = %.2f' % (lphi, eta_hb(lphi))
    for key, name, _ in LINES[:3]:
        e = eta_of(key)(lphi)
        row += '   %s %.2f (log ratio %+.2f)' % (name, e, np.log10(e / eta_hb(lphi)) if e > 0 else np.nan)
    print(row)
