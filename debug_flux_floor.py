"""How much of the disc falls below the ionizing-flux floor? The code clips
log Phi_H at 17 (pillar_disk.compute_log_ionizing_flux, direct method) while
the text says 15. Area fractions with the seed-1 pillars of Figures 6 and 8,
computed without any floor. A few seconds.

Run:  python3 debug_flux_floor.py
"""
import numpy as np

from pillardisk.pillar_disk import PillarDisk, load_config, resolve_rin

FLOAT_KEYS = {'rin', 'rout', 'h1', 'r0', 'beta', 'tv1', 'alpha', 'tx1',
              'tirrad_tvisc_ratio', 'fcol', 'hlamp', 'dmpc', 'cosi', 'redshift',
              'M_BH', 'r_isco_rg', 'f_trans', 'log_phi_inner', 'f_turb_line'}


def build_disk():
    """Disc from config_line.yaml with no flux floor (as measure_shadow_fraction.py)."""
    config = load_config('config_line.yaml')
    params = {}
    for section in ('disk', 'temperature', 'lamp', 'observation'):
        for key, value in config.get(section, {}).items():
            if isinstance(value, str) and value.lower() == 'auto':
                params[key] = 'auto'
            elif key in ('nr', 'nphi'):
                params[key] = int(float(value))
            elif key in FLOAT_KEYS:
                params[key] = float(value)
            else:
                params[key] = value
    params['cosi'] = np.cos(np.radians(params.pop('inclination')))
    params['no_fluxfloor'] = True
    resolve_rin(params)
    return PillarDisk(**params), config['pillars']


disk, pc = build_disk()
rng = np.random.default_rng(1)          # same draws as Figures 6 and 8
n = int(pc['N_pillar'])
r_p = np.clip(rng.normal(float(pc['r_mean']), float(pc['sig_r']), n),
              max(disk.rin, float(pc['rmin'])), disk.rout)
phi_p = rng.uniform(0., 2. * np.pi, n)
for i in range(n):
    disk.add_pillar(r_pillar=float(r_p[i]), phi_pillar=float(phi_p[i]),
                    height=float(pc['h_pillar']), sigma_r=float(pc['sigma_r_pillar']),
                    sigma_phi=float(pc['sigma_phi_pillar']),
                    modify_height=True, modify_temp=True, temp_factor=1.5)

r1 = np.linspace(0.2, 20.0, 400)
phi1 = np.linspace(0., 2. * np.pi, 720, endpoint=False)
r2, phi2 = np.meshgrid(r1, phi1, indexing='ij')
log_phi = np.asarray(disk.compute_log_ionizing_flux(r2, phi2))
w = r2                                     # area element for uniform dr, dphi
for lo, hi in ((0.2, 20.0), (2.0, 20.0), (5.0, 15.0)):
    sel = (r2 >= lo) & (r2 <= hi)
    print('%4.1f-%4.1f ld: below 17: %5.1f%%  below 15: %5.1f%%  '
          'median log Phi %.2f, 5th percentile %.2f'
          % (lo, hi, 100 * np.sum(w[sel] * (log_phi[sel] < 17)) / np.sum(w[sel]),
             100 * np.sum(w[sel] * (log_phi[sel] < 15)) / np.sum(w[sel]),
             np.median(log_phi[sel]), np.percentile(log_phi[sel], 5)))
for r in (2, 5, 10, 15, 20):
    row = log_phi[np.argmin(np.abs(r1 - r))]
    print('r = %2d ld: log Phi min %.2f, median %.2f, max %.2f'
          % (r, row.min(), np.median(row), row.max()))
