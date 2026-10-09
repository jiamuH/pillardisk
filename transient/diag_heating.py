#!/usr/bin/env python3
"""
diag_heating.py - does the model SED actually respond to the pillar heating?
Compute the flare-minus-quiescent bump at several pillar_temp values and
print the peak amplitude, to check whether the arm emission is entering the
SED (and whether its own emission is being self-occulted).

Run:  python3 transient/diag_heating.py
"""

import copy
import os
import sys

import numpy as np
import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from transient.transient_disk import disks_from_config  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(HERE, 'config_transient.yaml')) as f:
    CFG = yaml.safe_load(f)
Z = CFG['observation']['redshift']
TO_FLAM = 1e-17 / (1.0 + Z)
PROBES = np.array([2600., 3050., 3400., 3600., 4000., 5000.])


def flare_minus_quiet(pillar_temp, occultation=True):
    cfg = copy.deepcopy(CFG)
    cfg['transient']['pillar']['pillar_temp'] = pillar_temp
    cfg['transient']['occultation']['enabled'] = occultation
    flare, quiet, _ = disks_from_config(cfg, nr=250, nphi=500)
    mf = flare.compute_sed(PROBES) * TO_FLAM
    mq = quiet.compute_sed(PROBES) * TO_FLAM
    return mf, mq


print("rest A:      " + "".join(f"{int(w):>8d}" for w in PROBES))
print("--- difference (flare - quiescent), occultation ON ---")
for Tp in [0.0, 14000.0, 28000.0, 40000.0, 60000.0]:
    mf, mq = flare_minus_quiet(Tp, occultation=True)
    print(f"Tp={Tp/1e3:5.0f}kK " + "".join(f"{v:>8.1f}" for v in (mf - mq)))

print("--- difference, occultation OFF (pure arm emission) ---")
for Tp in [0.0, 14000.0, 28000.0, 40000.0, 60000.0]:
    mf, mq = flare_minus_quiet(Tp, occultation=False)
    print(f"Tp={Tp/1e3:5.0f}kK " + "".join(f"{v:>8.1f}" for v in (mf - mq)))

# how much of the arm's heated footprint is occulted (weight-averaged)?
cfg = copy.deepcopy(CFG)
flare, _, _ = disks_from_config(cfg, nr=250, nphi=500)
r2d, phi2d = np.meshgrid(flare.r, flare.phi, indexing='ij')
p = flare.pillars[0]
g = flare._pillar_footprint(r2d, phi2d, p, heat=True)
occ = flare.compute_observer_occultation()
vis_freac = np.sum(g * occ) / np.sum(g)
print(f"\nheated-footprint mean visibility (1=fully visible): {vis_freac:.2f}")
print(f"heated-footprint mean radius: "
      f"{np.sum(g * r2d) / np.sum(g):.2f} ld  (r_p = {p['r']} ld, "
      f"heat_r = {p.get('heat_r')} ld)")
