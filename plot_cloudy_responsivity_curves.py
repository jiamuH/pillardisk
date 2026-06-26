#!/usr/bin/env python3
r"""
Diagnostic plot of the Cloudy F-L curves and their local responsivities.

For each of C\,IV, Mg\,II, H\alpha (and H\beta for reference) at Z = Z_sun:

  - Top row: log10 F_line vs log Phi_H, interpolated from the
    RectBivariateSpline at 401 fine grid points. The slope of this
    curve at any point is exactly the responsivity eta below.
  - Bottom row: eta(log Phi_H) = d log F / d log Phi, computed by
    finite-differencing the same spline.

If eta is everywhere positive, the responsivity-weighted Psi has the
same sign as the emissivity-weighted Psi (only the per-cell scaling
changes). Negative eta means the F-L curve is declining at that Phi,
which would flip the sign of the cell's contribution.

Also prints a quick summary of the eta range per line and reports
whether negative eta values occur over the spliced [15, 21] grid.

Usage:
    python3 -m pillardisk.plot_cloudy_responsivity_curves
"""

import argparse
import os

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

from pillardisk.pillar_line_time_cloudy import _load_cloudy_raw_data


def main(outfile='plots/cloudy_responsivity_curves.png'):
    cloudy_file = ('/Users/jiamuh/c23.01/my_models/loc_metal_flux/'
                   'strong_LOC_varym_N25_v100_lineflux_LineList_BLR_Fe2_flux.txt')
    cloudy_ext = ('/Users/jiamuh/c23.01/my_models/loc_metal_flux/'
                  'strong_LOC_varym_N25_v100_lineflux_extlow_'
                  'LineList_BLR_Fe2_flux.txt')

    cloudy_data, phi_grid, Z_grid = _load_cloudy_raw_data(
        cloudy_file, Z_target=1.0, extension_file=cloudy_ext)

    print(f'phi_grid: {phi_grid}')
    print(f'Z_grid range: {Z_grid[0]} to {Z_grid[-1]} ({len(Z_grid)} pts)')

    # For the diagnostic we work with a 1-D slice at Z = Z_sun. We
    # interpolate log10 F (rather than F itself) using a natural cubic
    # spline, which is C^2 continuous so its derivative (= eta) is
    # smooth (C^1), with no kinks at knots. log-space interpolation
    # avoids the wild overshoot the linear-space spline shows for CIV
    # (whose F spans ~7 dex).
    Z_target = 1.0
    Z_idx = int(np.argmin(np.abs(Z_grid - Z_target)))
    log_F_grid = {
        key: np.log10(np.maximum(cloudy_data[key][:, Z_idx], 1e-300))
        for key in cloudy_data
    }
    log_interps = {
        key: CubicSpline(phi_grid, log_F_grid[key], bc_type='natural',
                         extrapolate=False)
        for key in cloudy_data
    }

    # Fine phi grid for plotting
    phi_fine = np.linspace(phi_grid[0], phi_grid[-1], 401)

    line_order = ['C4', 'Mg2', 'Halpha', 'HI']
    line_labels = {'C4': r'$\rm C\,IV$', 'Mg2': r'$\rm Mg\,II$',
                   'Halpha': r'$\rm H\alpha$', 'HI': r'$\rm H\beta$'}
    line_colors = {'C4': 'royalblue', 'Mg2': 'forestgreen',
                   'Halpha': 'orangered', 'HI': 'black'}

    plt.rcParams.update({'font.family': 'serif', 'font.size': 14})

    fig, axes = plt.subplots(2, 1, figsize=(9, 9), sharex=True)
    ax_F, ax_eta = axes

    print('\nResponsivity summary:')
    print(f'  {"line":<8} {"eta_min":>12} {"eta_max":>12} {"any negative?":>18}')
    for line_key in line_order:
        if line_key not in log_interps:
            continue
        # log F on the fine grid from the PCHIP interpolator
        log_F = log_interps[line_key](phi_fine)
        # eta is the analytic derivative of the PCHIP spline. This is
        # smoother and more accurate than np.gradient on a finite grid.
        eta_vals = log_interps[line_key].derivative()(phi_fine)
        eta_min = float(np.nanmin(eta_vals))
        eta_max = float(np.nanmax(eta_vals))
        any_neg = bool(eta_min < 0)
        print(f'  {line_key:<8} {eta_min:>12.3f} {eta_max:>12.3f} '
              f'{str(any_neg):>18}')

        # Also overplot the actual Cloudy grid points so we can see how
        # well the interpolation tracks them.
        color = line_colors.get(line_key, 'k')
        ax_F.plot(phi_fine, log_F, color=color, lw=2.5, alpha=0.85,
                  label=line_labels.get(line_key, line_key))
        ax_F.plot(phi_grid, log_F_grid[line_key], 'o', color=color,
                  ms=5, mfc='white', mec=color, mew=1.5)
        ax_eta.plot(phi_fine, eta_vals, color=color, lw=2.5, alpha=0.85,
                    label=line_labels.get(line_key, line_key))

    ax_F.set_ylabel(r'$\log_{10} F_{\rm line}~[\rm photons~cm^{-2}~s^{-1}]$',
                    fontsize=15)
    ax_F.legend(fontsize=13, frameon=False, loc='best')
    ax_F.minorticks_on()
    ax_F.tick_params(top=True, right=True, axis='both', which='major',
                     length=8, width=1.5, direction='in', labelsize=12)
    ax_F.tick_params(top=True, right=True, axis='both', which='minor',
                     length=5, width=1, direction='in')

    ax_eta.axhline(0, color='k', lw=1.0, alpha=0.5)
    ax_eta.set_xlabel(r'$\log_{10}\Phi_{\rm H}~[\rm photons~cm^{-2}~s^{-1}]$',
                      fontsize=15)
    ax_eta.set_ylabel(r'$\eta = d\log F / d\log\Phi$', fontsize=15)
    ax_eta.minorticks_on()
    ax_eta.tick_params(top=True, right=True, axis='both', which='major',
                       length=8, width=1.5, direction='in', labelsize=12)
    ax_eta.tick_params(top=True, right=True, axis='both', which='minor',
                       length=5, width=1, direction='in')

    plt.tight_layout()
    os.makedirs(os.path.dirname(outfile) or '.', exist_ok=True)
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'\nSaved {outfile}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outfile',
                        default='plots/cloudy_responsivity_curves.png')
    args = parser.parse_args()
    main(outfile=args.outfile)
