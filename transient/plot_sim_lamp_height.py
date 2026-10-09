#!/usr/bin/env python3
"""plot_sim_lamp_height.py - lamp-height experiment for the ionization
front geometry.

For lamp heights h = 0, 0.1, 0.2, 0.3 r0 on the disk axis, march
straight rays from the lamp at (0, 0, h) through the simulation density
cube in the (R, z) plane at the star's azimuth (the plane contains the
z axis, so the rays stay in it), accumulate recombinations per steradian
n^2 alpha_B s^2 ds (s = distance from the LAMP), and find where the
photon budget Q/4pi is exhausted. Overlays the resulting ionization
front curves on the density slice, showing how an elevated lamp sees
over the inner structures and pushes the shadows.

Output: plots/sim_lamp_height_fronts.png
Run:  python3 transient/plot_sim_lamp_height.py [dump.athdf]
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from transient.sim_spectrum import (  # noqa: E402
    load_sim, to_physical, star_position, add_colorbar, CONFIG, PLOTDIR,
    LD_CM)
from transient.sim_cloudy_spectrum import Q_ION, ALPHA_B  # noqa: E402

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

HEIGHTS = [0.0, 0.1, 0.2, 0.3]         # lamp heights [r0]
COLORS = ['gold', 'orangered', 'crimson', 'mediumorchid']
N_ALPHA = 600                          # ray directions per height
N_S = 1200                             # steps along each ray
S_MAX = 3.2                            # max path length [r0]


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else \
        os.path.join(HERE, 'data', 'sim', 'disk.out1.00012.athdf')
    cfg = CONFIG
    sim = load_sim(path)
    phys = to_physical(sim, cfg)
    r, th, ph = sim['r'], sim['th'], sim['ph']
    rf, thf = sim['rf'], sim['thf']
    r0_cm = cfg['R0_LD'] * LD_CM
    budget = Q_ION / (4 * np.pi)

    xs, ys = star_position(sim)
    ip = np.argmin(np.abs(ph - np.mod(np.arctan2(ys, xs), 2 * np.pi)))
    nH = phys['nH'][ip]                                # (nth, nr)
    rstar = np.hypot(xs, ys)

    # straight-ray march from (0, 0, h): angle alpha from the +z axis
    s = np.linspace(1e-3, S_MAX, N_S)
    ds_cm = (s[1] - s[0]) * r0_cm
    alpha = np.linspace(0.35, np.pi / 2 + 0.35, N_ALPHA)
    fronts = {}
    for h in HEIGHTS:
        Rp = s[None, :] * np.sin(alpha)[:, None]       # (na, ns)
        Zp = h + s[None, :] * np.cos(alpha)[:, None]
        rp = np.hypot(Rp, Zp)
        tp = np.arctan2(Rp, Zp)
        inside = ((rp >= rf[0]) & (rp <= rf[-1])
                  & (tp >= thf[0]) & (tp <= thf[-1]))
        it_ = np.clip(np.searchsorted(th, tp), 0, len(th) - 1)
        ir_ = np.clip(np.searchsorted(r, rp), 0, len(r) - 1)
        n = np.where(inside, nH[it_, ir_], 0.0)
        dS = n ** 2 * ALPHA_B * (s[None, :] * r0_cm) ** 2 * ds_cm
        S = np.cumsum(dS, axis=1)
        hit = S >= budget
        has = hit.any(axis=1)
        k = np.argmax(hit, axis=1)
        # sub-step interpolation of the crossing
        S_hi = S[np.arange(N_ALPHA), k]
        S_lo = np.where(k > 0, S[np.arange(N_ALPHA), np.maximum(k - 1, 0)],
                        0.0)
        fc = np.clip((budget - S_lo) / (S_hi - S_lo + 1e-300), 0, 1)
        s_if = s[np.maximum(k - 1, 0)] + fc * (s[1] - s[0])
        Rf = s_if * np.sin(alpha)
        Zf = h + s_if * np.cos(alpha)
        keep = has & (Zf >= -0.02)
        fronts[h] = (Rf[keep], Zf[keep])
        print(f"h = {h:.1f}: {100 * (~has).mean():.0f}% of directions "
              f"stay matter-bounded")

    # ---- figure: density slice + front curves ----
    up = th <= np.pi / 2
    thu = th[up]
    TH, RR = np.meshgrid(thu, r, indexing='ij')
    Rc, Zc = RR * np.sin(TH), RR * np.cos(TH)
    fig, ax = plt.subplots(figsize=(12, 5.6))
    pc = ax.pcolormesh(Rc, Zc, np.log10(nH[up]), cmap='Greys',
                       shading='auto')
    add_colorbar(pc, ax, r'$\log n_{\rm H}~[\rm cm^{-3}]$', fontsize=14)
    for h, c in zip(HEIGHTS, COLORS):
        Rf, Zf = fronts[h]
        o = np.argsort(np.arctan2(Rf, Zf - h))
        ax.plot(Rf[o], Zf[o], color=c, lw=2.5, alpha=0.9, zorder=4,
                label=rf'$h_{{\rm lamp}} = {h:.1f}~r_0$')
        ax.plot(0, h, marker='*', ms=22, mfc=c, mec='black', mew=1.2,
                ls='none', zorder=6, clip_on=False)
    ax.plot(rstar, 0, marker='*', ms=18, mfc='white', mec='black',
            mew=1.5, ls='none', zorder=6, clip_on=False)
    ax.annotate(r'$\rm star$', (rstar, 0), xytext=(rstar + 0.04, 0.04),
                textcoords='data', fontsize=15, color='white')
    ax.set_xlim(-0.08, r[-1] * 1.02)
    zmax = (r[-1] * np.cos(thu)).max()
    ax.set_ylim(-0.04, max(zmax * 1.30, 0.45))
    ax.set_xlabel(r'$R~[r_0]$', fontsize=17)
    ax.set_ylabel(r'$z~[r_0]$', fontsize=17)
    ax.legend(fontsize=13, loc='upper left', frameon=False)
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=13)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()
    os.makedirs(PLOTDIR, exist_ok=True)
    out = os.path.join(PLOTDIR, 'sim_lamp_height_fronts.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")


if __name__ == '__main__':
    main()
