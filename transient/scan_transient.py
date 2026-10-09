#!/usr/bin/env python3
"""
scan_transient.py - compare the two bump-source hypotheses under the
column-density (Balmer bound-free) absorption screen:

  Mode A: absorption + intrinsically brightened disk (no arm emission)
  Mode B: absorption + arm recombination emission (disk unchanged)

Run:  python3 transient/scan_transient.py
"""

import copy
import itertools
import os
import sys
import time

import numpy as np
import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from transient.transient_disk import disks_from_config  # noqa: E402
from transient.tune_transient import PROBES, data_targets  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(HERE, 'config_transient.yaml')) as fh:
    BASE_CFG = yaml.safe_load(fh)
Z = BASE_CFG['observation']['redshift']
TO_FLAM = 1e-17 / (1.0 + Z)

d = data_targets(Z)
target_diff = d['eboss2019'] - d['sdss2000']


def run_case(mods_transient, mods_pillar, mods_obs=None, mods_temp=None,
             nr=200, nphi=360):
    cfg = copy.deepcopy(BASE_CFG)
    cfg['transient'].update(mods_transient)
    cfg['transient']['pillar'].update(mods_pillar)
    if mods_obs:
        cfg['observation'].update(mods_obs)
    if mods_temp:
        cfg['temperature'].update(mods_temp)
    cfg['lamp']['hlamp'] = mods_temp.pop('_hlamp') if mods_temp and '_hlamp' in mods_temp else cfg['lamp']['hlamp']
    cfg['transient']['occultation']['n_ray'] = 200
    flare, quiet, _ = disks_from_config(cfg, nr=nr, nphi=nphi)
    m_flare = flare.compute_sed(PROBES) * TO_FLAM
    m_quiet = quiet.compute_sed(PROBES) * TO_FLAM
    scale = d['sdss2000'][np.argmin(np.abs(PROBES - 5150.0))] \
        / m_quiet[np.argmin(np.abs(PROBES - 5150.0))]
    diff = (m_flare - m_quiet) * scale
    return np.sum((diff - target_diff) ** 2), diff


results = []
t0 = time.time()

# Geometric-occultation variant (no Balmer opacity), TDE-heated arm kept,
# higher inclination for stronger occultation
for inc, h, temp in itertools.product(
        [75.0, 80.0, 84.0], [0.8, 1.2, 1.6],
        [20000.0, 28000.0, 36000.0]):
    chi2, diff = run_case(
        {'flare_tv1_scale': 1.0,
         'occultation': dict(BASE_CFG['transient']['occultation'],
                             mode='opaque')},
        {'height': h, 'pillar_temp': temp},
        mods_obs={'inclination': inc})
    results.append((chi2, f"i={inc:.0f} h={h:.1f} T={temp:.0f}", diff))

print(f"{len(results)} cases in {time.time() - t0:.0f} s\n")
results.sort(key=lambda x: x[0])
print("data diff " + "".join(f"{v:>8.1f}" for v in target_diff))
print("-" * 76)
for chi2, tag, diff in results[:6]:
    print(f"{tag:<26} chi2={chi2:8.0f}")
    print("   diff   " + "".join(f"{v:>8.1f}" for v in diff))
