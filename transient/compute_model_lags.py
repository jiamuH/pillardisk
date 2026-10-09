"""Model lag spectrum tau(lambda) for the transient-quasar geometry.

Builds the flare (arm + heating) and quiet (smooth bowl) disks from
transient/config_transient.yaml via disks_from_config, computes the
response-function mean delay vs wavelength with PillarDisk.compute_lag_spectrum
(rest-frame wavelengths; the arm enters through the get_height /
get_temperature overrides), and saves observed-frame results to
transient/data/model_lag_spectrum.npz. Note the observer occultation acts
only in compute_sed, so it is NOT part of this lag prediction.

Run: python3 transient/compute_model_lags.py [config.yaml] [--quick]
(--quick: 3 wavelengths only, to time the computation)
"""

import os
import sys
import time

import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from transient.transient_disk import disks_from_config  # noqa: E402

DATADIR = os.path.join(HERE, 'data')

NTAU = 500
TAUMAX = 150.0  # days; generous so the cold outer-disk response is not cut
WOBS_MIN, WOBS_MAX, NW = 3000.0, 10500.0, 25


def main():
    quick = '--quick' in sys.argv
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    cfgpath = args[0] if args else os.path.join(HERE, 'config_transient.yaml')
    with open(cfgpath) as fh:
        cfg = yaml.safe_load(fh)

    z = cfg['observation']['redshift']
    incl = cfg['observation']['inclination']
    flare, quiet, _ = disks_from_config(cfg)

    wobs = np.linspace(WOBS_MIN, WOBS_MAX, 3 if quick else NW)
    wrest = wobs / (1.0 + z)

    out = {'wobs': wobs, 'wrest': wrest, 'redshift': z, 'inclination': incl}
    for name, disk in [('flare', flare), ('quiet', quiet)]:
        t0 = time.time()
        tau_mean, _, _ = disk.compute_lag_spectrum(
            wrest, ntau=NTAU, taumax=TAUMAX, parallel=True)
        out[f'tau_{name}_rest'] = tau_mean
        out[f'tau_{name}_obs'] = tau_mean * (1.0 + z)
        print(f'{name}: {len(wrest)} wavelengths in {time.time() - t0:.1f} s; '
              f'tau_obs {out[f"tau_{name}_obs"].min():.2f}-'
              f'{out[f"tau_{name}_obs"].max():.2f} d')

    if quick:
        print('quick mode: not saving')
        return
    path = os.path.join(DATADIR, 'model_lag_spectrum.npz')
    np.savez(path, **out)
    print(f'saved {path}')


if __name__ == '__main__':
    main()
