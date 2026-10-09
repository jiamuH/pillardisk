#!/usr/bin/env python3
"""
diag_arm_height.py - effective occulting-wall height of the spiral arm,
for calibrating the toy model's `height` parameter.

For each azimuth phi and polar angle theta (upper half), compute the
radial (grazing) hydrogen column through the arm region,
    N(phi, theta) = int n_H dr  over r in [0.8, 1.6] r0,
and find the height z = r_arm (pi/2 - theta) at which N drops below an
opacity threshold:
    Thomson  tau_T  = 1 : N = 1.5e24 cm^-2  (gray electron scattering)
    dust     tau_V  = 1 : N = 1.9e21 cm^-2  (MW dust-to-gas; dusty arm)
The tallest wall section over azimuth is what matters for occultation.
Uses the standard anchors (r0 = 2 ld, n_mid = 1e10 cm^-3).

Run:  python3 transient/diag_arm_height.py
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from transient.sim_spectrum import load_sim, to_physical, CONFIG, LD_CM

R_ARM = 1.2          # representative arm radius [r0]
RRANGE = (0.8, 1.6)  # radial extent of the arm region for the column
THRESH = {'Thomson tau=1 (gray)': 1.5e24, 'dust tau_V=1 (dusty arm)': 1.9e21}


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else \
        os.path.join(HERE, 'data', 'sim', 'disk.out1.00012.athdf')
    cfg = CONFIG
    sim = load_sim(path)
    phys = to_physical(sim, cfg)
    r, th, ph = sim['r'], sim['th'], sim['ph']
    r0_cm = cfg['R0_LD'] * LD_CM
    dr = np.diff(sim['rf']) * r0_cm
    msel = (r >= RRANGE[0]) & (r <= RRANGE[1])
    # radial column at each (phi, theta) through the arm region [cm^-2]
    N = (phys['nH'][:, :, msel] * dr[None, None, msel]).sum(axis=2)

    up = th <= np.pi / 2                      # upper half (midplane at 90)
    z_of_th = R_ARM * (np.pi / 2 - th[up])    # height [r0]
    r_arm_ld = R_ARM * cfg['R0_LD']
    print(f"arm radius {R_ARM} r0 = {r_arm_ld:.1f} ld;  wedge top = "
          f"{R_ARM*(np.pi/2-th.min()):.2f} r0 = "
          f"{R_ARM*(np.pi/2-th.min())*cfg['R0_LD']:.2f} ld")
    for name, thr in THRESH.items():
        # per azimuth: highest theta (largest z) with N above threshold
        above = N[:, up] > thr                # (ph, nup); z decreasing? no:
        h_phi = np.where(above.any(axis=1),
                         z_of_th[np.argmax(above, axis=1)], 0.0)
        # argmax finds FIRST True going from the wedge top downward
        # (th ascending = z descending), i.e. the highest opaque point
        med, p95, hmax = np.percentile(h_phi, [50, 95, 100])
        capped = np.mean(above[:, 0]) * 100   # opaque already at wedge top
        print(f"{name}: wall height h = {med:.3f} (median) / "
              f"{p95:.3f} (95%) / {hmax:.3f} (max) r0 "
              f"= {med*cfg['R0_LD']:.2f} / {p95*cfg['R0_LD']:.2f} / "
              f"{hmax*cfg['R0_LD']:.2f} ld"
              f"   [{capped:.0f}% of azimuths saturate the wedge -> "
              f"lower limit]" if capped > 0 else "")
        hr = p95 / R_ARM
        icrit = np.degrees(np.arctan2(1.0, hr))
        print(f"   -> h/r = {hr:.3f} (95%), occultation needs "
              f"i > arctan(r/h) = {icrit:.0f} deg")


if __name__ == '__main__':
    main()
