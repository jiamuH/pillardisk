#!/usr/bin/env python3
"""
diag_occ_vs_inc.py - how much of the disk UV/optical flux the arm OCCULTS
as a function of inclination, isolated from the arm's emission and from the
flux-conserving self-dimming. Uses compute_occulted_fraction(), which
compares the SED with the occultation on vs off at the SAME temperatures.

Run:  python3 transient/diag_occ_vs_inc.py
"""

import copy
import os
import sys

import numpy as np
import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from transient.transient_disk import disks_from_config  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
PROBES = np.array([2400., 2800., 3400., 5000.])
INCS = [0., 30., 45., 60., 70., 80.]


def main():
    with open(os.path.join(HERE, 'config_transient.yaml')) as f:
        cfg = yaml.safe_load(f)
    var = cfg['transient'].get('variant_geometric', {})

    print("Occulted flux fraction (geometric/opaque mode, h = "
          f"{var.get('height', cfg['transient']['pillar']['height'])} ld, "
          f"r_p = {cfg['transient']['pillar']['r_pillar']} ld)")
    print("rest A:  " + "".join(f"{int(w):>8d}" for w in PROBES))
    for inc in INCS:
        cfg_i = copy.deepcopy(cfg)
        cfg_i['observation']['inclination'] = inc
        cfg_i['transient']['occultation']['mode'] = 'opaque'
        if 'height' in var:
            cfg_i['transient']['pillar']['height'] = var['height']
        if 'pillar_temp' in var:
            cfg_i['transient']['pillar']['pillar_temp'] = var['pillar_temp']
        flare, _, _ = disks_from_config(cfg_i)
        occ = flare.compute_occulted_fraction(PROBES)
        crit = np.degrees(np.arctan(
            cfg['transient']['pillar']['r_pillar']
            / cfg_i['transient']['pillar']['height']))
        tag = "  <- occultation turns on" if inc >= crit else ""
        print(f"i={inc:4.0f}:" + "".join(f"{v:>8.2f}" for v in occ) + tag)
    print(f"\ncritical i ~ arctan(r_p/h) = {crit:.0f} deg")


if __name__ == '__main__':
    main()
