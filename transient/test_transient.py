#!/usr/bin/env python3
"""
test_transient.py - quick sanity tests for TransientPillarDisk.

Run:  python3 transient/test_transient.py
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pillar_disk import PillarDisk  # noqa: E402
from transient.transient_disk import TransientPillarDisk  # noqa: E402

# Quasar-like base parameters (small grid for speed where the slow parent runs).
# Bowl rim at the outer edge (r0 = rout); tv1 = 790 K at 100 ld gives
# T_visc(1 ld) ~ 2.5e4 K with the r^-3/4 law.
BASE = dict(rin=0.04, rout=100.0, h1=2.0, r0=100.0, beta=10.0,
            tv1=790.0, tirrad_tvisc_ratio=0.5, hlamp=0.05,
            dmpc=2750.0, cosi=0.5, fcol=1.0)

PILLAR = dict(r_pillar=1.0, phi_pillar=0.0, height=1.5,
              sigma_r=0.5, sigma_phi=1.5, pillar_temp=9000.0)

npass = 0
nfail = 0


def check(name, ok, detail):
    global npass, nfail
    status = "PASS" if ok else "FAIL"
    if ok:
        npass += 1
    else:
        nfail += 1
    print(f"[{status}] {name}: {detail}")


# ----------------------------------------------------------------------
# Test 1: parent equivalence (no pillars, occultation off)
# ----------------------------------------------------------------------
waves = np.logspace(np.log10(1500.0), np.log10(9000.0), 8)
parent = PillarDisk(nr=80, nphi=60, **BASE)
child = TransientPillarDisk(nr=80, nphi=60, occultation=False, **BASE)
f_parent = parent.compute_sed(waves)
f_child = child.compute_sed(waves)
maxdev = np.max(np.abs(f_child / f_parent - 1.0))
check("1 parent equivalence", maxdev < 1e-3,
      f"max fractional deviation = {maxdev:.2e} (< 1e-3 required)")

# ----------------------------------------------------------------------
# Test 2: opaque tall wall suppresses far-UV; face-on does not
# ----------------------------------------------------------------------
waves_uv = np.array([1500.0, 2000.0, 2500.0])
disk = TransientPillarDisk(nr=200, nphi=360, occultation=True, **BASE)
disk.add_pillar(**PILLAR)
occ_frac = disk.compute_occulted_fraction(waves_uv)
check("2a inclined occultation",
      occ_frac[0] > 0.5 and np.all(np.diff(occ_frac) < 0),
      f"occulted flux fraction at 1500/2000/2500 A = "
      f"{occ_frac[0]:.2f}/{occ_frac[1]:.2f}/{occ_frac[2]:.2f} "
      f"(> 0.5 at 1500 A and decreasing with wavelength required)")

face_on = dict(BASE, cosi=0.99)
disk_f = TransientPillarDisk(nr=200, nphi=360, occultation=True, **face_on)
disk_f.add_pillar(**PILLAR)
occ_frac_f = disk_f.compute_occulted_fraction(waves_uv)
check("2b face-on limit", np.all(occ_frac_f < 0.05),
      f"face-on occulted fraction = {np.max(occ_frac_f):.3f} (< 0.05 required)")

# ----------------------------------------------------------------------
# Test 3: heated pillar alone -> bump peak near blackbody peak
# ----------------------------------------------------------------------
waves_b = np.logspace(np.log10(1500.0), np.log10(9000.0), 60)
disk_h = TransientPillarDisk(nr=200, nphi=360, occultation=False, **BASE)
disk_h.add_pillar(**PILLAR)
f_flare = disk_h.compute_sed(waves_b)
f_quiet = disk_h.compute_sed_no_pillars(waves_b)
# difference in f_lambda units (f_nu / lambda^2), amplitude arbitrary
diff_flam = (f_flare - f_quiet) / waves_b ** 2
lam_peak = waves_b[np.argmax(diff_flam)]
lam_wien = 2.8978e7 / PILLAR['pillar_temp']  # B_lambda peak for T_TDE
check("3 heated-pillar bump peak", abs(lam_peak / lam_wien - 1.0) < 0.35,
      f"difference-spectrum peak = {lam_peak:.0f} A vs blackbody peak "
      f"{lam_wien:.0f} A (within 35% required; geometry shifts it)")

# ----------------------------------------------------------------------
# Test 4: ray-march convergence and smoothness
# ----------------------------------------------------------------------
disk_c1 = TransientPillarDisk(nr=200, nphi=360, occultation=True, n_ray=300, **BASE)
disk_c1.add_pillar(**PILLAR)
disk_c2 = TransientPillarDisk(nr=200, nphi=360, occultation=True, n_ray=600, **BASE)
disk_c2.add_pillar(**PILLAR)
f1 = disk_c1.compute_sed(waves_b)
f2 = disk_c2.compute_sed(waves_b)
conv = np.max(np.abs(f1 / f2 - 1.0))
check("4a n_ray convergence", conv < 0.005,
      f"n_ray 300 (default) vs 600 max deviation = {conv:.2e} (< 0.5% required)")

# smoothness of the flare-quiescent difference: no small-scale wiggles
d = np.gradient(np.log(np.abs((f2 - disk_c2.compute_sed_no_pillars(waves_b))
                              / waves_b ** 2) + 1e-300))
sign_changes = np.sum(np.diff(np.sign(np.gradient(d))) != 0)
check("4b smoothness", sign_changes <= 6,
      f"curvature sign changes in log-difference spectrum = {sign_changes} "
      f"(few = smooth)")

# ----------------------------------------------------------------------
# Test 5: far-side pillar -> no UV occultation
# ----------------------------------------------------------------------
disk_far = TransientPillarDisk(nr=200, nphi=360, occultation=True, **BASE)
far = dict(PILLAR, phi_pillar=np.pi)
disk_far.add_pillar(**far)
occ_far = disk_far.compute_occulted_fraction(waves_uv)
check("5 far-side pillar", np.all(occ_far < 0.1),
      f"far-side occulted fraction at 1500/2000/2500 A = "
      f"{occ_far[0]:.3f}/{occ_far[1]:.3f}/{occ_far[2]:.3f} (< 0.1 required)")

print(f"\n{npass} passed, {nfail} failed")
sys.exit(1 if nfail else 0)
