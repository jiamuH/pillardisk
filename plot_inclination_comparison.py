#!/usr/bin/env python3
"""
Compare time-evolving C\\,IV velocity-time maps for the 1-pillar setup at
several inclination angles, demonstrating that the near-side (suppression)
and far-side (enhancement) asymmetry along the line of sight grows with
inclination.

Runs the same disk + 1-pillar configuration three times, varying only the
inclination (default i = 0, 45, 80 deg), and stacks the resulting C\\,IV
F(lambda, t) maps vertically. Residuals are not plotted -- the raw flux
maps already show the foreshortening modulation at high i.

Usage:
    python3 -m pillardisk.plot_inclination_comparison [config_line.yaml]
"""

import argparse
import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, NullFormatter

from pillardisk.pillar_disk import load_config
from pillardisk.pillar_line_time_cloudy import (
    C, KM_TO_CM,
    compute_time_evolving_cloudy_map,
    load_cloudy_models,
    v_virial_from_mbh,
)
from pillardisk.plot_pillar_snapshots import _build_disk_with_one_pillar


def _set_inclination_deg(disk, i_deg):
    """Update disk's inclination-related attributes in place."""
    cosi = float(np.cos(np.radians(i_deg)))
    sini = float(np.sin(np.radians(i_deg)))
    disk.cosi = cosi
    disk.sini = sini
    disk.ex = sini
    disk.ey = 0.0
    disk.ez = cosi


def main(config_file='config_line.yaml',
         inclinations_deg=(30.0, 45.0, 80.0),
         tmax_days=1200.0,
         outfile='plots/fig_inclination_civ_1pillar.png'):
    config = load_config(config_file)
    if config is None:
        raise SystemExit(f'Could not load config: {config_file}')

    disk, r_pillar, phi_pillar_0 = _build_disk_with_one_pillar(config)
    print(f'Pillar at r = {r_pillar:.2f} ld, phi_0 = '
          f'{np.degrees(phi_pillar_0):.1f} deg')

    # Only need C IV.
    lambda0_dict = {'C4': 1549.0}
    line_params = config.get('emission_line', {})
    if 'M_BH' in line_params:
        v_virial = v_virial_from_mbh(float(line_params['M_BH']), disk.r0)
    else:
        v_virial = float(line_params.get('v_virial', 1000.0))

    cloudy_file = ('/Users/jiamuh/c23.01/my_models/loc_metal_flux/'
                   'strong_LOC_varym_N25_v100_lineflux_LineList_BLR_Fe2_flux.txt')
    cloudy_extension = ('/Users/jiamuh/c23.01/my_models/loc_metal_flux/'
                        'strong_LOC_varym_N25_v100_lineflux_extlow_'
                        'LineList_BLR_Fe2_flux.txt')
    cloudy_interp_dict, phi_grid, _ = load_cloudy_models(
        cloudy_file, Z_target=1.0, extension_file=cloudy_extension)

    comp_params = config.get('computation', {})
    nlambda = int(comp_params.get('nlambda', 200))
    ntime = int(comp_params.get('ntime', 50))
    n_cores = comp_params.get('n_cores', None)
    if n_cores is not None:
        n_cores = int(n_cores)

    # Run for each inclination.
    results = {}
    lambda_grid = None
    time_grid = None
    for i_deg in inclinations_deg:
        _set_inclination_deg(disk, i_deg)
        print(f'\n=== inclination = {i_deg:.0f} deg ===')
        lg, tg, fm = compute_time_evolving_cloudy_map(
            disk, cloudy_interp_dict, phi_grid, lambda0_dict, v_virial,
            nlambda=nlambda, ntime=ntime, tmax=tmax_days,
            use_absolute_flux=True, n_cores=n_cores,
        )
        results[i_deg] = fm['C4']
        lambda_grid = lg['C4']
        time_grid = tg

    # ── Plot: vertically stacked C IV panels, one per inclination ──
    n_rows = len(inclinations_deg)
    fig, axes = plt.subplots(n_rows, 1, figsize=(10, 3.6 * n_rows + 0.6),
                              sharex=True)
    if n_rows == 1:
        axes = [axes]

    lambda0 = lambda0_dict['C4']
    v_ticks_major = np.array([-9000, -6000, -3000, 0, 3000, 6000, 9000])
    v_ticks_minor = np.array([-7500, -4500, -1500, 1500, 4500, 7500])

    # Fixed shared colour scale [0, 3] so all panels (and the height
    # sweep figure) are directly comparable.
    vmax_shared = 3.0

    for k, i_deg in enumerate(inclinations_deg):
        ax = axes[k]
        F = results[i_deg]
        norm = float(np.mean(F))
        if not np.isfinite(norm) or norm == 0:
            norm = 1.0
        F_norm = F / norm

        im = ax.pcolormesh(lambda_grid, time_grid, F_norm,
                           shading='auto', cmap='rainbow',
                           vmin=0.0, vmax=vmax_shared)
        ax.axvline(lambda0, color='white', ls='--', lw=1.5, alpha=0.7)
        cb = plt.colorbar(im, ax=ax, shrink=0.85, aspect=20, extend='neither')
        cb.set_label(r'$F\,/\,\langle F\rangle$', fontsize=15)
        cb.ax.tick_params(labelsize=12)

        # Inclination label as in-panel text.
        ax.text(0.03, 0.93, rf'$i={i_deg:.0f}^\circ$',
                transform=ax.transAxes, fontsize=18,
                color='white', ha='left', va='top',
                bbox=dict(boxstyle='round,pad=0.25',
                          facecolor='black', alpha=0.45,
                          edgecolor='none'))

        ax.set_ylabel(r'$\rm Time~[days]$', fontsize=18)
        ax.set_xlim(lambda_grid[0], lambda_grid[-1])
        ax.set_ylim(time_grid[0], time_grid[-1])
        ax.minorticks_on()
        ax.tick_params(top=False, right=True, axis='both', which='major',
                       length=8, width=1.5, direction='in', labelsize=14)
        ax.tick_params(top=False, right=True, axis='both', which='minor',
                       length=5, width=1, direction='in')

        if k == n_rows - 1:
            ax.set_xlabel(r'$\lambda~[\rm \AA]$', fontsize=18)
        else:
            ax.set_xlabel('')

        # Top velocity axis only on the top panel. We use ax.twiny()
        # rather than secondary_xaxis('top', ...): the latter has a
        # zero-height bbox that silently clips inward-pointing tick
        # marks (direction='in') to 0 px, so the ticks vanish.
        if k == 0:
            v_lo = (C / KM_TO_CM) * (lambda_grid[0] - lambda0) / lambda0
            v_hi = (C / KM_TO_CM) * (lambda_grid[-1] - lambda0) / lambda0
            ax2 = ax.twiny()
            ax2.set_xlim(v_lo, v_hi)
            major = v_ticks_major[(v_ticks_major >= v_lo)
                                  & (v_ticks_major <= v_hi)]
            minor = v_ticks_minor[(v_ticks_minor >= v_lo)
                                  & (v_ticks_minor <= v_hi)]
            ax2.xaxis.set_major_locator(FixedLocator(major))
            ax2.xaxis.set_minor_locator(FixedLocator(minor))
            ax2.xaxis.set_minor_formatter(NullFormatter())
            ax2.tick_params(axis='x', which='major',
                            length=8, width=1.5, direction='in',
                            labelsize=14, top=True, bottom=False,
                            labeltop=True, labelbottom=False)
            ax2.tick_params(axis='x', which='minor',
                            length=5, width=1.0, direction='in',
                            top=True, bottom=False)
            ax2.xaxis.set_label_position('top')
            ax2.set_xlabel(r'$v~[\rm km\,s^{-1}]$', fontsize=18, labelpad=8)

    plt.tight_layout()
    os.makedirs(os.path.dirname(outfile) or '.', exist_ok=True)
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'\nInclination comparison plot saved to {outfile}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config', nargs='?', default='config_line.yaml',
                        help='Config YAML file (default: config_line.yaml)')
    parser.add_argument('--incs', type=float, nargs='+',
                        default=[30.0, 45.0, 80.0],
                        help='Inclination angles in deg (default: 30 45 80)')
    parser.add_argument('--tmax', type=float, default=1200.0,
                        help='Total time span in days (default: 1200)')
    parser.add_argument('--outfile',
                        default='plots/fig_inclination_civ_1pillar.png',
                        help='Output PNG path')
    args = parser.parse_args()
    main(config_file=args.config,
         inclinations_deg=tuple(args.incs),
         tmax_days=args.tmax,
         outfile=args.outfile)
