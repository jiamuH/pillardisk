"""Referee reply, Section 7 bullet 2: measure the shadowed fraction of the
pillar belt quoted in the Results subsection on ionizing flux maps
("Approximately 20--30% of the disc surface in the pillar belt (5 < r < 15 ld)
is shadowed at the level Phi_H/Phi_0 < 0.5").

Figure 1 was made by pillar_line_time_cloudy.main, which draws the 100 pillar
positions without a fixed seed, so its exact layout cannot be reproduced. This
script measures the fraction over several independent draws from the same
distribution, with the config_line.yaml parameters and no flux floor (as in
Figure 1). Two definitions, both weighted by disc area within 5 < r < 15 ld:
  mask:  shadow mask S < 0.5
  flux:  Phi_H(with pillars) / Phi_H(same disc, no pillars) < 0.5

Run:  python3 measure_shadow_fraction.py
"""
import os

import numpy as np

from pillardisk.pillar_disk import PillarDisk, load_config, resolve_rin

CONFIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config_line.yaml')
FLOAT_KEYS = {'rin', 'rout', 'h1', 'r0', 'beta', 'tv1', 'alpha', 'tx1',
              'tirrad_tvisc_ratio', 'fcol', 'hlamp', 'dmpc', 'cosi', 'redshift',
              'M_BH', 'r_isco_rg', 'f_trans', 'log_phi_inner', 'f_turb_line'}
INT_KEYS = {'nr', 'nphi'}
N_DRAWS = 10
BELT = (5.0, 15.0)


def build_disk(**overrides):
    config = load_config(CONFIG)
    params = {}
    for section in ('disk', 'temperature', 'lamp', 'observation'):
        for key, value in config.get(section, {}).items():
            if isinstance(value, str) and value.lower() == 'auto':
                params[key] = 'auto'
            elif key in INT_KEYS:
                params[key] = int(float(value))
            elif key in FLOAT_KEYS:
                params[key] = float(value)
            else:
                params[key] = value
    params['cosi'] = np.cos(np.radians(params.pop('inclination')))
    params.update(overrides)
    resolve_rin(params)
    return PillarDisk(**params), config['pillars']


def add_random_pillars(disk, pillars_cfg, rng):
    """Same distribution as pillar_line_time_cloudy.main (make_many branch)."""
    n = int(pillars_cfg['N_pillar'])
    r_p = rng.normal(float(pillars_cfg['r_mean']), float(pillars_cfg['sig_r']), n)
    r_p = np.clip(r_p, max(disk.rin, float(pillars_cfg['rmin'])), disk.rout)
    phi_p = rng.uniform(0., 2. * np.pi, n)
    for i in range(n):
        disk.add_pillar(r_pillar=float(r_p[i]), phi_pillar=float(phi_p[i]),
                        height=float(pillars_cfg['h_pillar']),
                        sigma_r=float(pillars_cfg['sigma_r_pillar']),
                        sigma_phi=float(pillars_cfg['sigma_phi_pillar']),
                        modify_height=True, modify_temp=True, temp_factor=1.0)


r1 = np.linspace(0.5, 20.0, 400)
phi1 = np.linspace(0., 2. * np.pi, 720, endpoint=False)
r2, phi2 = np.meshgrid(r1, phi1, indexing='ij')
in_belt = (r2 > BELT[0]) & (r2 < BELT[1])
area = r2 * in_belt                      # uniform dr, dphi: area element ~ r

reference, pillars_cfg = build_disk(no_fluxfloor=True)
log_phi_bare = np.asarray(reference.compute_log_ionizing_flux(r2, phi2))

mask_fractions, flux_fractions = [], []
for seed in range(N_DRAWS):
    disk, _ = build_disk(no_fluxfloor=True)
    add_random_pillars(disk, pillars_cfg, np.random.default_rng(seed))
    mask = disk._compute_shadow_mask(r2, phi2, disk.get_height(r2, phi2))
    log_phi = np.asarray(disk.compute_log_ionizing_flux(r2, phi2))
    ratio = 10. ** (log_phi - log_phi_bare)
    mask_fractions.append(np.sum(area * (mask < 0.5)) / np.sum(area))
    flux_fractions.append(np.sum(area * (ratio < 0.5)) / np.sum(area))
    print('draw %2d  mask S < 0.5: %5.1f per cent   flux ratio < 0.5: %5.1f per cent'
          % (seed, 100 * mask_fractions[-1], 100 * flux_fractions[-1]))

for label, values in (('mask S < 0.5', mask_fractions),
                      ('flux ratio < 0.5', flux_fractions)):
    values = 100. * np.array(values)
    print('%-18s mean %.1f per cent, range %.1f to %.1f per cent over %d draws'
          % (label, values.mean(), values.min(), values.max(), N_DRAWS))
