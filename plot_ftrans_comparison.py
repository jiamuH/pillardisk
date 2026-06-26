#!/usr/bin/env python3
"""
Compare time-evolving CIV/MgII/Halpha velocity-time maps for the 1-pillar
configuration computed with f_trans = 0.0 (fully opaque pillar) vs
f_trans = 0.3 (partially transparent), and plot the difference per line.

Two values of f_trans are run back-to-back through the same disk + Cloudy
setup that ``make_time_evolving_figures.py`` uses; only ``disk.f_trans``
is varied between calls. Output is a (n_lines x 3) figure: F(f_trans=0),
F(f_trans=0.3), and the difference F_0 - F_0.3.

Usage:
    python3 -m pillardisk.plot_ftrans_comparison [config_line.yaml]
"""

import argparse
import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, NullFormatter

from pillardisk.pillar_disk import load_config
from pillardisk.pillar_line_time_cloudy import (
    C, KM_TO_CM,
    compute_orbital_period,
    compute_time_evolving_cloudy_map,
    load_cloudy_models,
    v_virial_from_mbh,
)
from pillardisk.plot_pillar_snapshots import _build_disk_with_one_pillar


def _run_one_f_trans(disk, f_trans, cloudy_interp_dict, phi_grid,
                     lambda0_dict, v_virial, nlambda, ntime, tmax,
                     n_cores=None):
    """Set ``disk.f_trans`` and run the time-evolving Cloudy map. Returns
    ``(lambda_grid_dict, time_grid, flux_map_dict)`` for the with-pillars
    case."""
    disk.f_trans = float(f_trans)
    print(f"\n=== f_trans = {disk.f_trans} ===")
    return compute_time_evolving_cloudy_map(
        disk, cloudy_interp_dict, phi_grid, lambda0_dict, v_virial,
        nlambda=nlambda, ntime=ntime, tmax=tmax,
        use_absolute_flux=True, n_cores=n_cores,
    )


def _add_velocity_axis(ax, lambda_grid, lambda0):
    lam2v = lambda lam, l0=lambda0: (C / KM_TO_CM) * (lam - l0) / l0
    v2lam = lambda v, l0=lambda0: l0 * (1.0 + v * KM_TO_CM / C)
    ax2 = ax.secondary_xaxis('top', functions=(lam2v, v2lam))
    ax2.set_xlabel(r'$v~[\rm km\,s^{-1}]$', fontsize=16, labelpad=8)
    v_lo, v_hi = lam2v(lambda_grid[0]), lam2v(lambda_grid[-1])
    major = np.array([-9000, -6000, -3000, 0, 3000, 6000, 9000])
    minor = np.array([-7500, -4500, -1500, 1500, 4500, 7500])
    ax2.xaxis.set_major_locator(FixedLocator(major[(major >= v_lo) & (major <= v_hi)]))
    ax2.xaxis.set_minor_locator(FixedLocator(minor[(minor >= v_lo) & (minor <= v_hi)]))
    ax2.xaxis.set_minor_formatter(NullFormatter())
    ax2.tick_params(axis='both', which='major', length=8, width=1.5,
                    direction='in', labelsize=13)
    ax2.tick_params(axis='both', which='minor', length=5, width=1, direction='in')
    return ax2


def _style_main_axis(ax, lambda_grid, time_grid):
    ax.set_xlim(lambda_grid[0], lambda_grid[-1])
    ax.set_ylim(time_grid[0], time_grid[-1])
    ax.set_xlabel(r'$\lambda~[\rm \AA]$', fontsize=16)
    ax.set_ylabel(r'$\rm Time~[days]$', fontsize=16)
    ax.minorticks_on()
    ax.tick_params(top=False, right=True, axis='both', which='major',
                   length=8, width=1.5, direction='in', labelsize=13)
    ax.tick_params(top=False, right=True, axis='both', which='minor',
                   length=5, width=1, direction='in')


def main(config_file='config_line.yaml',
         f_trans_values=(0.0, 0.3),
         tmax_days=1500.0,
         outfile='plots/fig_ftrans_comparison_1pillar.png'):
    config = load_config(config_file)
    if config is None:
        raise SystemExit(f'Could not load config: {config_file}')

    # Build the same 1-pillar disk that the snapshot script uses.
    disk, r_pillar, phi_pillar_0 = _build_disk_with_one_pillar(config)

    # Emission line setup, mirroring time_cloudy_main.
    lambda0_dict = {'Halpha': 6562.8, 'Mg2': 2798.0, 'C4': 1549.0}
    line_params = config.get('emission_line', {})
    if 'M_BH' in line_params:
        v_virial = v_virial_from_mbh(float(line_params['M_BH']), disk.r0)
    else:
        v_virial = float(line_params.get('v_virial', 1000.0))
    print(f'v_virial = {v_virial:.0f} km/s at r0 = {disk.r0:.1f} ld')
    print(f'Pillar at r = {r_pillar:.2f} ld, phi_0 = '
          f'{np.degrees(phi_pillar_0):.1f} deg')

    # Cloudy models (absolute-flux grid + low-phi extension, same defaults
    # as time_cloudy_main).
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

    # Run twice; keep the resulting flux maps in memory.
    results = {}
    lambda_grid_dict = None
    time_grid = None
    for ft in f_trans_values:
        lg, tg, fm = _run_one_f_trans(
            disk, ft, cloudy_interp_dict, phi_grid, lambda0_dict,
            v_virial, nlambda, ntime, tmax_days, n_cores=n_cores)
        results[ft] = fm
        lambda_grid_dict = lg
        time_grid = tg

    # ── Plot ─────────────────────────────────────────────────
    line_order = ['C4', 'Mg2', 'Halpha']
    line_labels = {'Halpha': r'$\rm H\alpha$',
                   'Mg2':    r'$\rm Mg\,II$',
                   'C4':     r'$\rm C\,IV$'}

    n_rows = len(line_order)
    fig, axes = plt.subplots(n_rows, 3, figsize=(20, 4.6 * n_rows))
    if n_rows == 1:
        axes = axes[np.newaxis, :]

    ft_lo, ft_hi = f_trans_values
    for i, line_key in enumerate(line_order):
        lam = lambda_grid_dict[line_key]
        lam0 = lambda0_dict[line_key]
        F_lo = results[ft_lo][line_key]
        F_hi = results[ft_hi][line_key]

        # Per-line normalisation by max of mean spectrum at f_trans=ft_hi
        # (so all three columns share the same denominator).
        ref_mean = np.mean(F_hi, axis=0)
        norm = float(np.nanmax(ref_mean))
        if not np.isfinite(norm) or norm == 0:
            norm = 1.0

        F_lo_n = F_lo / norm
        F_hi_n = F_hi / norm
        diff_n = (F_lo - F_hi) / norm

        # Column 1: F at ft_lo
        ax0 = axes[i, 0]
        im0 = ax0.pcolormesh(lam, time_grid, F_lo_n, shading='auto',
                             cmap='rainbow', vmin=0.0, vmax=1.0)
        ax0.axvline(lam0, color='white', ls='--', lw=1.5, alpha=0.7)
        cb0 = plt.colorbar(im0, ax=ax0, shrink=0.8, aspect=20, extend='neither')
        cb0.set_label(r'$F\,/\,\max\langle F\rangle_t$', fontsize=14)
        cb0.ax.tick_params(labelsize=12)
        ax0.set_title(rf'{line_labels[line_key]}, $f_{{\rm trans}}={ft_lo}$',
                      fontsize=17, pad=10)

        # Column 2: F at ft_hi
        ax1 = axes[i, 1]
        im1 = ax1.pcolormesh(lam, time_grid, F_hi_n, shading='auto',
                             cmap='rainbow', vmin=0.0, vmax=1.0)
        ax1.axvline(lam0, color='white', ls='--', lw=1.5, alpha=0.7)
        cb1 = plt.colorbar(im1, ax=ax1, shrink=0.8, aspect=20, extend='neither')
        cb1.set_label(r'$F\,/\,\max\langle F\rangle_t$', fontsize=14)
        cb1.ax.tick_params(labelsize=12)
        ax1.set_title(rf'{line_labels[line_key]}, $f_{{\rm trans}}={ft_hi}$',
                      fontsize=17, pad=10)

        # Column 3: difference (ft_lo - ft_hi)
        ax2_panel = axes[i, 2]
        vmax = float(np.nanmax(np.abs(diff_n)))
        if not np.isfinite(vmax) or vmax == 0:
            vmax = 1.0
        im2 = ax2_panel.pcolormesh(lam, time_grid, diff_n, shading='auto',
                                   cmap='seismic', vmin=-vmax, vmax=vmax)
        ax2_panel.axvline(lam0, color='black', ls='--', lw=1.5, alpha=0.7)
        cb2 = plt.colorbar(im2, ax=ax2_panel, shrink=0.8, aspect=20,
                           extend='neither')
        cb2.set_label(r'$(F_{0} - F_{0.3})\,/\,\max\langle F\rangle_t$',
                      fontsize=14)
        cb2.ax.tick_params(labelsize=12)
        ax2_panel.set_title(rf'{line_labels[line_key]}, $\Delta$ '
                            rf'$(f_{{\rm trans}}={ft_lo}-{ft_hi})$',
                            fontsize=17, pad=10)

        for ax in (ax0, ax1, ax2_panel):
            _style_main_axis(ax, lam, time_grid)
            _add_velocity_axis(ax, lam, lam0)

    plt.tight_layout()
    os.makedirs(os.path.dirname(outfile) or '.', exist_ok=True)
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'\nf_trans comparison plot saved to {outfile}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config', nargs='?', default='config_line.yaml',
                        help='Config YAML file (default: config_line.yaml)')
    parser.add_argument('--ftrans-lo', type=float, default=0.0,
                        help='Lower f_trans (default: 0.0)')
    parser.add_argument('--ftrans-hi', type=float, default=0.3,
                        help='Higher f_trans (default: 0.3)')
    parser.add_argument('--tmax', type=float, default=1500.0,
                        help='Total time span in days (default: 1500)')
    parser.add_argument('--outfile',
                        default='plots/fig_ftrans_comparison_1pillar.png',
                        help='Output PNG path')
    args = parser.parse_args()
    main(config_file=args.config,
         f_trans_values=(args.ftrans_lo, args.ftrans_hi),
         tmax_days=args.tmax,
         outfile=args.outfile)
