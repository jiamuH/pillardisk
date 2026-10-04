"""Referee reply, Section 7 bullet 2: check two statements about pillar
shadows before consolidating the Results subsection on disc geometry.

(1) Transparency: with the fiducial config (f_trans = 0.3), how deep is the
    shadow in the shadow mask and in log Phi_H?
(2) Rim re-illumination: does the rising bowl rim see the lamp again beyond
    the shadow lane inside r_out = 20 light days, or does the lane run to the
    disc edge?

Uses the real PillarDisk shadow mask and ionizing-flux code with the values in
config_line.yaml.

Run:  python3 debug_shadow_rim_ftrans.py
"""
import os

import numpy as np

from pillardisk.pillar_disk import PillarDisk, load_config, resolve_rin

CONFIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config_line.yaml')
FLOAT_KEYS = {'rin', 'rout', 'h1', 'r0', 'beta', 'tv1', 'alpha', 'tx1',
              'tirrad_tvisc_ratio', 'fcol', 'hlamp', 'dmpc', 'cosi', 'redshift',
              'M_BH', 'r_isco_rg', 'f_trans', 'log_phi_inner', 'f_turb_line'}
INT_KEYS = {'nr', 'nphi'}


def build_disk(f_trans_override=None):
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
    if f_trans_override is not None:
        params['f_trans'] = f_trans_override
    resolve_rin(params)
    return PillarDisk(**params)


def lane_profile(label, r_p, phi_p, height, sigma_r, sigma_phi, f_trans=None):
    disk = build_disk(f_trans)
    disk.add_pillar(r_pillar=r_p, phi_pillar=phi_p, height=height,
                    sigma_r=sigma_r, sigma_phi=sigma_phi,
                    modify_height=True, modify_temp=True, temp_factor=1.0)
    # Full 2D grid with radius along axis 0, as the surface-normal code expects.
    r1 = np.linspace(0.5, disk.rout, 400)
    phi1 = np.linspace(0., 2. * np.pi, 1440, endpoint=False)
    r2, phi2 = np.meshgrid(r1, phi1, indexing='ij')
    mask2 = disk._compute_shadow_mask(r2, phi2, disk.get_height(r2, phi2))
    log_phi2 = np.asarray(disk.compute_log_ionizing_flux(r2, phi2))

    j_in = np.argmin(np.abs(np.angle(np.exp(1j * (phi1 - phi_p)))))
    j_out = np.argmin(np.abs(np.angle(np.exp(1j * (phi1 - phi_p - np.pi / 2.)))))
    behind = r1 > r_p + 0.3
    radii = r1[behind]
    mask = mask2[behind, j_in]
    log_phi_in = log_phi2[behind, j_in]
    log_phi_out = log_phi2[behind, j_out]

    shadowed = mask < 0.99
    print('\n%s: r_p = %.1f, h_p = %.2f, f_trans = %.2f, h_lamp = %.2f, r_out = %.0f'
          % (label, r_p, height, disk.f_trans, disk.hlamp, disk.rout))
    print('  minimum shadow mask along the lane           %.3f' % mask.min())
    if shadowed.any():
        last = radii[shadowed].max()
        print('  shadow reaches out to r = %.2f light days (%s)'
              % (last, 'disc edge' if last >= radii[-2] else 're-illuminated beyond'))
    else:
        print('  no shadow along the pillar azimuth')
    for r_probe in (r_p + 1., 10., 15., 19.5):
        if r_probe > radii[-1] or r_probe < radii[0]:
            continue
        k = np.argmin(np.abs(radii - r_probe))
        print('  r = %5.2f  mask %.3f  log Phi_H behind pillar %.2f  unshadowed %.2f  drop %.2f dex'
              % (radii[k], mask[k], log_phi_in[k], log_phi_out[k],
                 log_phi_out[k] - log_phi_in[k]))


if __name__ == '__main__':
    lane_profile('single pillar, fiducial config', 4.0, 3 * np.pi / 4, 0.15, 0.5, 0.2)
    lane_profile('single pillar, opaque (f_trans = 0) for comparison',
                 4.0, 3 * np.pi / 4, 0.15, 0.5, 0.2, f_trans=0.0)
    for r_p in (2.0, 4.0, 10.0):
        lane_profile('N = 100 style pillar', r_p, 1.0, 0.04, 0.1, 0.05)
