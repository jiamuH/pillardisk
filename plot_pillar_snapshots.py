#!/usr/bin/env python3
"""
Plot the single-pillar disc geometry at several orbital phases.

Builds the same disc + 1-pillar setup that ``make_time_evolving_figures.py``
uses, then sweeps the pillar's azimuth ``phi_p`` over one orbital period and
saves a face-on / 3D snapshot at each phase. Useful for visually verifying
that the pillar position, shadow, and disc tilt look correct at every step
of the time-evolving computation (independently of the spectral output).

Usage:
    python3 -m pillardisk.plot_pillar_snapshots [config_line.yaml]
"""

import argparse
import os

import numpy as np

from pillardisk.pillar_disk import PillarDisk, load_config, resolve_rin
from pillardisk.pillar_line import parse_math_expr
from pillardisk.pillar_line_time_cloudy import (
    compute_orbital_period, v_virial_from_mbh,
)


def _to_list(x):
    return x if isinstance(x, list) else [x]


def _build_disk_with_one_pillar(config):
    """Construct a ``PillarDisk`` with exactly one pillar, taken as the first
    entry of the manual pillar list in ``config``. Returns
    ``(disk, r_pillar, phi_pillar_0)``."""
    disk_params = config.get('disk', {}).copy()
    temp_params = config.get('temperature', {}).copy()
    lamp_params = config.get('lamp', {}).copy()
    obs_params = config.get('observation', {}).copy()

    float_params = ['rin', 'rout', 'h1', 'r0', 'beta', 'tv1', 'alpha',
                    'tx1', 'tirrad_tvisc_ratio', 'fcol', 'hlamp', 'dmpc',
                    'cosi', 'inclination', 'redshift', 'M_BH', 'r_isco_rg',
                    'f_trans', 'f_turb_line']
    int_params = ['nr', 'nphi']

    for params in [disk_params, temp_params, lamp_params, obs_params]:
        for key, value in list(params.items()):
            if value is None:
                continue
            if isinstance(value, str) and value.lower() == 'auto':
                continue
            if key in int_params:
                params[key] = int(float(value))
            elif key in float_params:
                params[key] = float(value)

    if 'inclination' in obs_params and obs_params['inclination'] is not None:
        obs_params['cosi'] = np.cos(np.radians(obs_params['inclination']))
        del obs_params['inclination']

    merged = {**disk_params, **temp_params, **lamp_params, **obs_params}
    resolve_rin(merged)
    disk = PillarDisk(**merged)

    pillars_cfg = config.get('pillars', {})
    r_pillar = float(_to_list(pillars_cfg.get('r_pillar', [4.0]))[0])
    phi_pillar_0 = parse_math_expr(_to_list(pillars_cfg.get('phi_pillar', [0.0]))[0])
    height = float(_to_list(pillars_cfg.get('height', [0.15]))[0])
    sigma_r = float(_to_list(pillars_cfg.get('sigma_r', [0.5]))[0])
    sigma_phi = parse_math_expr(_to_list(pillars_cfg.get('sigma_phi', [0.2]))[0])
    temp_factor = float(_to_list(pillars_cfg.get('temp_factor', [5]))[0])
    modify_height = bool(_to_list(pillars_cfg.get('modify_height', [True]))[0])
    modify_temp = bool(_to_list(pillars_cfg.get('modify_temp', [True]))[0])

    disk.add_pillar(r_pillar=r_pillar, phi_pillar=phi_pillar_0,
                    height=height, sigma_r=sigma_r, sigma_phi=sigma_phi,
                    temp_factor=temp_factor,
                    modify_height=modify_height, modify_temp=modify_temp)

    return disk, r_pillar, phi_pillar_0


def main(config_file='config_line.yaml', n_snaps=8, tmax_days=1500.0,
         outdir='plots/pillar_snapshots', show_light_rays=False):
    config = load_config(config_file)
    if config is None:
        raise SystemExit(f'Could not load config: {config_file}')

    disk, r_pillar, phi_pillar_0 = _build_disk_with_one_pillar(config)

    line_params = config.get('emission_line', {})
    if 'M_BH' in line_params:
        v_virial = v_virial_from_mbh(float(line_params['M_BH']), disk.r0)
    else:
        v_virial = float(line_params.get('v_virial', 1000.0))

    T_orbital = compute_orbital_period(r_pillar, v_virial, disk.r0)
    omega_p = 2.0 * np.pi / T_orbital
    print(f'Pillar at r = {r_pillar:.2f} ld, phi_0 = '
          f'{np.degrees(phi_pillar_0):.1f} deg')
    print(f'Orbital period T_p = {T_orbital:.1f} days')
    print(f'Snapshot range: 0 to {tmax_days:.0f} days '
          f'({tmax_days/T_orbital:.2f} orbits)')

    os.makedirs(outdir, exist_ok=True)

    times_days = np.linspace(0.0, tmax_days, n_snaps, endpoint=False)
    for i, t_days in enumerate(times_days):
        phi_now = np.mod(phi_pillar_0 + omega_p * t_days, 2.0 * np.pi)
        for p in disk.pillars:
            p['phi'] = phi_now

        n_orb = t_days / T_orbital
        title = (rf'$t = {t_days:.0f}~\rm d$  '
                 rf'$({n_orb:.2f}\,T_p)$,  '
                 rf'$\phi_p = {np.degrees(phi_now):.0f}^\circ$')
        outfile = os.path.join(outdir, f'pillar_snapshot_{i:02d}.png')
        disk.plot_3d_geometry(show_light_rays=show_light_rays,
                              color_by='flux',
                              filename=outfile,
                              title=title,
                              n_phi_plot=180, n_r_plot=100)
        print(f'  [{i+1}/{n_snaps}] t={t_days:7.1f} d  '
              f'({n_orb:5.2f} T_p)  phi_p={np.degrees(phi_now):6.1f} deg  '
              f'-> {outfile}')

    print(f'\nDone. {n_snaps} snapshots saved to {outdir}/')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config', nargs='?', default='config_line.yaml',
                        help='Config YAML file (default: config_line.yaml)')
    parser.add_argument('--n-snaps', type=int, default=8,
                        help='Number of snapshots over [0, tmax) (default: 8)')
    parser.add_argument('--tmax', type=float, default=1500.0,
                        help='Total time span in days (default: 1500, '
                             'matches make_time_evolving_figures.py)')
    parser.add_argument('--outdir', default='plots/pillar_snapshots',
                        help='Output directory (default: plots/pillar_snapshots)')
    parser.add_argument('--rays', action='store_true',
                        help='Include lamp-to-pillar light rays in 3D view')
    args = parser.parse_args()
    main(config_file=args.config, n_snaps=args.n_snaps,
         tmax_days=args.tmax, outdir=args.outdir,
         show_light_rays=args.rays)
