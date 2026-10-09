#!/usr/bin/env python3
"""
tune_transient.py - quick metrics for hand-tuning the transient pillar model
against the observed epochs. Prints data targets vs model values at a few
rest wavelengths; no plots.

Run:  python3 transient/tune_transient.py [config_transient.yaml]
"""

import os
import sys

import numpy as np
import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from transient.transient_disk import TransientPillarDisk  # noqa: E402
from transient.make_transient_figure import load_epoch, rebin_R  # noqa: E402

# probe wavelengths chosen to avoid emission-line windows
PROBES = np.array([2620., 3050., 3220., 3620., 4050., 4550., 5150., 6050.])


def data_targets(z):
    out = {}
    for name in ('sdss2000', 'eboss2019'):
        ep = load_epoch(name)
        if ep is None:
            return None
        wb, fb = rebin_R(ep['wave_obs'], ep['flam'], ep['ivar'], R=300.0)
        out[name] = np.interp(PROBES, wb / (1.0 + z), fb,
                              left=np.nan, right=np.nan)
    return out


def main():
    cfg_path = sys.argv[1] if len(sys.argv) > 1 else \
        os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     'config_transient.yaml')
    with open(cfg_path) as f:
        cfg = yaml.safe_load(f)

    z = cfg['observation']['redshift']
    from transient.transient_disk import disks_from_config
    flare_disk, quiet_disk, dmpc = disks_from_config(cfg)

    fnu_flare = flare_disk.compute_sed(PROBES)
    fnu_quiet = quiet_disk.compute_sed(PROBES)
    # parent SED output is rest-frame f_lambda x 1e26 (see transient_disk.py)
    to_flam = 1e-17 / (1.0 + z)
    m_flare = fnu_flare * to_flam
    m_quiet = fnu_quiet * to_flam

    d = data_targets(z)
    i5000 = np.argmin(np.abs(PROBES - 5000.0))
    scale = 1.0
    if d is not None:
        scale = d['sdss2000'][i5000] / m_quiet[i5000]
    print(f"raw model quiescent flam at rest 5000 A = {m_quiet[i5000]:.3e} "
          f"(1e-17 cgs);  scale to data = {scale:.3e}")
    m_flare_s = m_flare * scale
    m_quiet_s = m_quiet * scale

    hdr = "rest A   " + "".join(f"{int(w):>9d}" for w in PROBES)
    print(hdr)
    if d is not None:
        print("dat 2000 " + "".join(f"{v:>9.1f}" for v in d['sdss2000']))
        print("dat 2019 " + "".join(f"{v:>9.1f}" for v in d['eboss2019']))
        print("dat diff " + "".join(
            f"{v:>9.1f}" for v in d['eboss2019'] - d['sdss2000']))
        print("dat ratio" + "".join(
            f"{v:>9.2f}" for v in d['eboss2019'] / d['sdss2000']))
    print("mod quiet" + "".join(f"{v:>9.1f}" for v in m_quiet_s))
    print("mod flare" + "".join(f"{v:>9.1f}" for v in m_flare_s))
    print("mod diff " + "".join(f"{v:>9.1f}" for v in m_flare_s - m_quiet_s))
    print("mod ratio" + "".join(f"{v:>9.2f}" for v in m_flare_s / m_quiet_s))

    occ_frac = flare_disk.compute_occulted_fraction(PROBES)
    print("occ frac " + "".join(f"{v:>9.2f}" for v in occ_frac))


if __name__ == '__main__':
    main()
