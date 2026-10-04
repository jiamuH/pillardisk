#!/usr/bin/env python3
r"""
Side-by-side comparison of emissivity- and responsivity-weighted
velocity-delay maps for the 100-pillar configuration.

Each row is one emission line (C\,IV, Mg\,II, H\alpha). Each row has
three panels:

    | psi_emiss / max(psi_emiss) | psi_resp / max(|psi_resp|) | diff |

where the third panel is

    diff(lambda, tau) = psi_resp / max(|psi_resp|)
                       - psi_emiss / max(psi_emiss).

If responsivity is exactly proportional to emissivity for every cell,
both normalised maps coincide and the diff panel is zero everywhere.
Any cell-by-cell departure (non-power-law F-L curve, sign flip, etc.)
shows up directly in the diff panel.

Usage:
    python3 -m pillardisk.make_vdm_weighting_compare [config_line.yaml]
"""

import argparse
import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import AutoMinorLocator, MaxNLocator

VDM_CMAP = LinearSegmentedColormap.from_list('vdm_wrb',
                                              ['white', 'red', 'blue'])

from pillardisk.pillar_disk import load_config
from pillardisk.pillar_line_time_cloudy import (
    C, KM_TO_CM,
    _load_cloudy_raw_data,
    compute_velocity_delay_map_cloudy,
    v_virial_from_mbh,
)
from pillardisk.plot_pillar_snapshots import _build_disk_with_one_pillar
from pillardisk.make_time_evolving_figures import make_many_pillars_config


def _build_disk_with_n100(config_file):
    """Build the disk with the random 100-pillar ensemble (make_many=True)."""
    tmp = make_many_pillars_config(config_file)
    try:
        cfg = load_config(tmp)
    finally:
        os.unlink(tmp)
    # _build_disk_with_one_pillar only handles the manual list, so rebuild
    # from cfg manually.
    from pillardisk.pillar_disk import PillarDisk, resolve_rin
    disk_params = cfg.get('disk', {}).copy()
    temp_params = cfg.get('temperature', {}).copy()
    lamp_params = cfg.get('lamp', {}).copy()
    obs_params = cfg.get('observation', {}).copy()

    float_params = ['rin', 'rout', 'h1', 'r0', 'beta', 'tv1', 'alpha',
                    'tx1', 'tirrad_tvisc_ratio', 'fcol', 'hlamp', 'dmpc',
                    'cosi', 'inclination', 'redshift', 'M_BH', 'r_isco_rg',
                    'f_trans', 'f_turb_line']
    int_params = ['nr', 'nphi']
    for p in [disk_params, temp_params, lamp_params, obs_params]:
        for key, value in list(p.items()):
            if value is None:
                continue
            if isinstance(value, str) and value.lower() == 'auto':
                continue
            if key in int_params:
                p[key] = int(float(value))
            elif key in float_params:
                p[key] = float(value)
    if 'inclination' in obs_params and obs_params['inclination'] is not None:
        obs_params['cosi'] = np.cos(np.radians(obs_params['inclination']))
        del obs_params['inclination']
    merged = {**disk_params, **temp_params, **lamp_params, **obs_params}
    resolve_rin(merged)
    disk = PillarDisk(**merged)

    # Sample N random pillars matching the config.
    pillars_cfg = cfg.get('pillars', {})
    N = int(pillars_cfg.get('N_pillar', 100))
    r_mean = float(pillars_cfg.get('r_mean', 10.0))
    sig_r = float(pillars_cfg.get('sig_r', 5.0))
    rmin = float(pillars_cfg.get('rmin', 2.0))
    h_p = float(pillars_cfg.get('h_pillar', 0.04))
    sigma_r_p = float(pillars_cfg.get('sigma_r_pillar', 0.1))
    sigma_phi_p = float(pillars_cfg.get('sigma_phi_pillar', 0.05))
    modify_h = bool(pillars_cfg.get('modify_height_pillar', True))
    modify_T = bool(pillars_cfg.get('modify_temp_pillar', True))
    temp_factor = float(pillars_cfg.get('temp_factor_pillar', 1.0))

    np.random.seed(42)  # reproducibility for the comparison plot
    r_pillars = np.clip(np.random.normal(r_mean, sig_r, N), rmin, disk.rout)
    phi_pillars = np.random.uniform(0, 2 * np.pi, N)
    for r_p, phi_p in zip(r_pillars, phi_pillars):
        disk.add_pillar(r_pillar=float(r_p), phi_pillar=float(phi_p),
                        height=h_p, sigma_r=sigma_r_p, sigma_phi=sigma_phi_p,
                        modify_height=modify_h, modify_temp=modify_T,
                        temp_factor=temp_factor)
    return disk, cfg


def _add_velocity_top_axis(ax, lambda_grid, lambda0):
    # Plain twiny() — a regular full Axes pinned on top of `ax`. Its xlim
    # is the velocity range that corresponds to ax's wavelength range, so
    # both axes line up. Using twiny() (instead of secondary_xaxis) makes
    # direction='in' ticks render reliably; secondary_xaxis has a zero-
    # height bbox that swallows inward ticks even with clip_on=False.
    v_lo = (C / KM_TO_CM) * (lambda_grid[0] - lambda0) / lambda0
    v_hi = (C / KM_TO_CM) * (lambda_grid[-1] - lambda0) / lambda0
    ax2 = ax.twiny()
    # Follow the parent's size when a colorbar divider has shrunk it.
    ax2.set_axes_locator(ax.get_axes_locator())
    ax2.set_xlim(v_lo, v_hi)
    # Few labelled velocities so the larger tick labels do not collide.
    from matplotlib.ticker import FixedLocator, MultipleLocator
    ax2.xaxis.set_major_locator(FixedLocator([-6000, 0, 6000]))
    ax2.xaxis.set_minor_locator(MultipleLocator(1500))
    ax2.tick_params(axis='x', which='major', direction='in',
                    length=8, width=1.5, labelsize=17,
                    top=True, bottom=False, labeltop=True, labelbottom=False)
    ax2.tick_params(axis='x', which='minor', direction='in',
                    length=5, width=1.0,
                    top=True, bottom=False)
    ax2.set_xlabel(r'$v~[\rm km\,s^{-1}]$', fontsize=17, labelpad=6)
    ax2.xaxis.set_label_position('top')
    return ax2


def _add_colorbar(im, ax, label):
    """Colorbar with exactly the height of its panel; ticks facing inward."""
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    cax = make_axes_locatable(ax).append_axes('right', size='4.5%', pad=0.12)
    cb = ax.figure.colorbar(im, cax=cax, extend='neither')
    cb.set_label(label, fontsize=17)
    cb.ax.minorticks_on()
    cb.ax.tick_params(which='major', direction='in', length=6, width=1.2,
                      labelsize=17)
    cb.ax.tick_params(which='minor', direction='in', length=3, width=0.8)
    return cb


def main(config_file='config_line.yaml',
         outfile='plots/fig_vdm_weighting_compare.png',
         cache_file='plots/vdm_weighting_compare_maps.npz',
         lowres=False, replot=False):
    lambda0_dict = {'C4': 1549.0, 'Mg2': 2798.0, 'Halpha': 6562.8}
    line_order = ['C4', 'Mg2', 'Halpha']
    line_labels = {'C4': r'$\rm C\,IV$', 'Mg2': r'$\rm Mg\,II$',
                   'Halpha': r'$\rm H\alpha$'}

    if replot:
        d = np.load(cache_file)
        tau_grid = d['tau_grid']
        lambda_grid_dict = {k: d[f'lam_{k}'] for k in line_order}
        results_emiss = {k: d[f'emiss_{k}'] for k in line_order}
        results_resp = {k: d[f'resp_{k}'] for k in line_order}
        _plot_figure(line_order, lambda_grid_dict, tau_grid, results_emiss,
                     results_resp, d['emiss_HI'], d['resp_HI'],
                     lambda0_dict, line_labels, outfile=outfile)
        return

    if lowres:
        # Coarse disc grid for quick layout tests (16 times fewer cells).
        import tempfile
        import yaml
        with open(config_file) as f:
            raw = yaml.safe_load(f)
        raw['disk']['nr'] = 150
        raw['disk']['nphi'] = 90
        tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.yaml',
                                          delete=False)
        yaml.dump(raw, tmp, default_flow_style=False)
        tmp.close()
        config_file = tmp.name

    print('Building disk with 100 pillars...')
    disk, cfg = _build_disk_with_n100(config_file)

    line_params = cfg.get('emission_line', {})
    if 'M_BH' in line_params:
        v_virial = v_virial_from_mbh(float(line_params['M_BH']), disk.r0)
    else:
        v_virial = float(line_params.get('v_virial', 1000.0))
    print(f'v_virial = {v_virial:.0f} km/s at r0 = {disk.r0:.1f} ld')

    cloudy_file = ('/Users/jiamuh/c23.01/my_models/loc_metal_flux/'
                   'strong_LOC_varym_N25_v100_lineflux_LineList_BLR_Fe2_flux.txt')
    cloudy_ext = ('/Users/jiamuh/c23.01/my_models/loc_metal_flux/'
                  'strong_LOC_varym_N25_v100_lineflux_extlow_LineList_BLR_Fe2_flux.txt')
    print('Loading Cloudy grids (with extension)...')
    cloudy_data_raw, phi_grid, Z_grid = _load_cloudy_raw_data(
        cloudy_file, Z_target=1.0, extension_file=cloudy_ext)

    # Build per-line interpolators in the main process (no multiprocessing
    # because we're only running a few maps and want clean error messages).
    from scipy.interpolate import RectBivariateSpline
    cloudy_interp_dict = {
        key: RectBivariateSpline(phi_grid, Z_grid, cloudy_data_raw[key],
                                 kx=2, ky=2)
        for key in cloudy_data_raw
    }

    nlambda = 200
    ntau = 100
    taumax = 50.0

    results_emiss = {}
    results_resp = {}
    lambda_grid_dict = {}
    tau_grid = None
    for line_key in line_order:
        l0 = lambda0_dict[line_key]
        print(f'\n=== {line_key} (lambda0={l0} A) ===')
        print('  computing emissivity-weighted...')
        lg, tg, psi_e = compute_velocity_delay_map_cloudy(
            disk, line_key, l0, v_virial, cloudy_interp_dict,
            nlambda=nlambda, ntau=ntau, taumax=taumax,
            weighting='emissivity')
        print('  computing responsivity-weighted...')
        _, _, psi_r = compute_velocity_delay_map_cloudy(
            disk, line_key, l0, v_virial, cloudy_interp_dict,
            nlambda=nlambda, ntau=ntau, taumax=taumax,
            weighting='responsivity')
        results_emiss[line_key] = psi_e
        results_resp[line_key] = psi_r
        lambda_grid_dict[line_key] = lg
        tau_grid = tg
        # Quick diagnostic: how different are they cell-by-cell?
        max_e = float(np.nanmax(np.abs(psi_e)))
        max_r = float(np.nanmax(np.abs(psi_r)))
        print(f'  max|psi_emiss|={max_e:.3e}, max|psi_resp|={max_r:.3e}')
        if max_e > 0:
            ratio = max_r / max_e
            print(f'  amplitude ratio resp/emiss = {ratio:.3e}')

    # H-beta maps for both weightings, the denominators of the ratio maps.
    # All lines share the same velocity grid, so maps divide pixel by pixel.
    print('\n=== HI (H-beta, lambda0=4861.32 A) ===')
    _, _, hb_e = compute_velocity_delay_map_cloudy(
        disk, 'HI', 4861.32, v_virial, cloudy_interp_dict,
        nlambda=nlambda, ntau=ntau, taumax=taumax, weighting='emissivity')
    _, _, hb_r = compute_velocity_delay_map_cloudy(
        disk, 'HI', 4861.32, v_virial, cloudy_interp_dict,
        nlambda=nlambda, ntau=ntau, taumax=taumax, weighting='responsivity')

    os.makedirs(os.path.dirname(cache_file) or '.', exist_ok=True)
    np.savez(cache_file, tau_grid=tau_grid, emiss_HI=hb_e, resp_HI=hb_r,
             **{f'lam_{k}': lambda_grid_dict[k] for k in line_order},
             **{f'emiss_{k}': results_emiss[k] for k in line_order},
             **{f'resp_{k}': results_resp[k] for k in line_order})
    print(f'Saved maps to {cache_file}')

    _plot_figure(line_order, lambda_grid_dict, tau_grid, results_emiss,
                 results_resp, hb_e, hb_r, lambda0_dict, line_labels,
                 outfile=outfile)


def _plot_figure(line_order, lambda_grid_dict, tau_grid, results_emiss,
                 results_resp, hb_e, hb_r, lambda0_dict, line_labels,
                 outfile=None):
    """Render one 3x3 figure (rows = lines, cols = emiss/resp/diff).

    If ``normalize=True`` each map is divided by ∬ Ψ dλ dτ so the panels
    are unit-integrated probability densities and we use the shared
    ``psi_vmax_fixed = 1e-3`` scale of the paper. If ``normalize=False``
    the raw Ψ values are plotted with a per-row vmax set to the joint
    max of the emiss/resp maps (each row gets its own scale because
    absolute amplitudes differ by orders of magnitude across lines).
    """
    n_rows = len(line_order)
    fig, axes = plt.subplots(n_rows, 3, figsize=(18, 4.6 * n_rows))
    if n_rows == 1:
        axes = axes[np.newaxis, :]

    for i, line_key in enumerate(line_order):
        lam = lambda_grid_dict[line_key]
        l0 = lambda0_dict[line_key]
        psi_e = results_emiss[line_key]
        psi_r = results_resp[line_key]

        # Line response relative to H-beta, under each weighting. Both maps
        # of a pair come from the same disc cells, so the shared geometry
        # (bright short-delay band, Keplerian envelope) cancels.
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio_e = np.where((psi_e > 0) & (hb_e > 0),
                               np.log10(psi_e / hb_e), np.nan)
            ratio_r = np.where((psi_r > 0) & (hb_r > 0),
                               np.log10(psi_r / hb_r), np.nan)
        both = np.concatenate([ratio_e[np.isfinite(ratio_e)],
                               ratio_r[np.isfinite(ratio_r)]])
        r_lo, r_hi = np.nanpercentile(both, [2, 98])
        print(f'  {line_key}: log(X/Hbeta) emissivity median '
              f'{np.nanmedian(ratio_e):.2f}, responsivity median '
              f'{np.nanmedian(ratio_r):.2f}; colour range {r_lo:.2f} to {r_hi:.2f}')

        # Mean delay at each wavelength, and the overall mean delay, for both
        # weightings (the ratio is independent of the map normalisation).
        def _mean_delay_spectrum(psi):
            col = psi.sum(axis=0)
            with np.errstate(divide='ignore', invalid='ignore'):
                return np.where(col > 0, (psi * tau_grid[:, None]).sum(axis=0) / col,
                                np.nan)

        def _mean_delay(psi):
            prof = psi.sum(axis=1)
            return float((prof * tau_grid).sum() / prof.sum())

        tau_e = _mean_delay_spectrum(psi_e)
        tau_r = _mean_delay_spectrum(psi_r)
        mean_e, mean_r = _mean_delay(psi_e), _mean_delay(psi_r)
        print(f'  {line_key}: mean delay emissivity {mean_e:.1f} d, '
              f'responsivity {mean_r:.1f} d (ratio {mean_r / mean_e:.2f})')

        # --- Columns 1 and 2: log(X / H-beta), same color scale ---
        for col, (rmap, name) in enumerate(((ratio_e, r'Emissivity'),
                                            (ratio_r, r'Responsivity'))):
            ax = axes[i, col]
            im = ax.pcolormesh(lam, tau_grid, rmap, shading='auto',
                               cmap='viridis', vmin=r_lo, vmax=r_hi)
            ax.axvline(l0, color='white', ls='--', lw=1.5, alpha=0.8)
            _add_colorbar(im, ax, r'$\log\,(\Psi/\Psi_{\rm H\beta})$')
            ax.text(0.03, 0.95,
                    line_labels[line_key] + r'$\rm /H\beta~' + name + r'$',
                    transform=ax.transAxes, fontsize=17, color='black',
                    ha='left', va='top',
                    bbox=dict(boxstyle='round,pad=0.25', facecolor='white',
                              alpha=0.7, edgecolor='none'))

        # --- Column 3: difference of the two ratio maps, in dex ---
        ax = axes[i, 2]
        dratio = ratio_r - ratio_e
        print(f'  {line_key}: Delta log(X/Hbeta), median |.| = '
              f'{np.nanmedian(np.abs(dratio)):.2f} dex, 98th percentile = '
              f'{np.nanpercentile(np.abs(dratio), 98):.2f} dex')
        im = ax.pcolormesh(lam, tau_grid, dratio, shading='auto',
                           cmap='RdBu', vmin=-0.5, vmax=0.5)
        ax.axvline(l0, color='black', ls='--', lw=1.5, alpha=0.7)
        _add_colorbar(im, ax, r'$\Delta\log\,(\Psi/\Psi_{\rm H\beta})$')
        ax.text(0.03, 0.95, line_labels[line_key] + r'$\rm /H\beta~Difference$',
                transform=ax.transAxes, fontsize=17, color='black',
                ha='left', va='top',
                bbox=dict(boxstyle='round,pad=0.25', facecolor='white',
                          alpha=0.7, edgecolor='none'))

        # Axis cosmetics for all three panels in this row.
        for ax in axes[i, :]:
            ax.set_xlabel(r'$\lambda~[\rm \AA]$', fontsize=17)
            ax.set_ylabel(r'$\tau~[\rm days]$', fontsize=17)
            ax.set_xlim(lam[0], lam[-1])
            ax.set_ylim(0, tau_grid[-1])
            ax.minorticks_on()
            ax.tick_params(axis='both', which='major', length=8, width=1.5,
                           direction='in', labelsize=17, right=True)
            ax.tick_params(axis='both', which='minor', length=5, width=1.0,
                           direction='in', right=True)
            ax.tick_params(axis='x', which='both', top=False,
                           labeltop=False)
            _add_velocity_top_axis(ax, lam, l0)

    plt.tight_layout()
    os.makedirs(os.path.dirname(outfile) or '.', exist_ok=True)
    plt.savefig(outfile, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'\nSaved {outfile}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config', nargs='?', default='config_line.yaml')
    parser.add_argument('--outfile',
                        default='plots/fig_vdm_weighting_compare.png')
    parser.add_argument('--cache',
                        default='plots/vdm_weighting_compare_maps.npz')
    parser.add_argument('--lowres', action='store_true',
                        help='coarse disc grid for a quick layout test')
    parser.add_argument('--replot', action='store_true',
                        help='replot from the cached maps, no recomputation')
    args = parser.parse_args()
    main(config_file=args.config, outfile=args.outfile,
         cache_file=args.cache, lowres=args.lowres, replot=args.replot)
