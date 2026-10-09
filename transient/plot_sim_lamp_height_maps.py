#!/usr/bin/env python3
"""plot_sim_lamp_height_maps.py - face-on H alpha maps for elevated
lamp heights h = 0, 0.1, 0.2, 0.3 r0.

Straight rays are marched from the lamp at (0, 0, h) through the
density cube (the lamp sits on the axis, so each ray stays in its
meridional plane). The absorbed ionizing photons per path step follow
from clipping the cumulative recombination integral at the photon
budget, and H alpha is painted where the photons are absorbed:
L(Halpha) = 0.45 photons per recombination (case B) x E(6563 A). This
is exact for a recombination line (the Cloudy G(xi) curve for H alpha
is the identity) and needs no slab lookup, so the lamp-centered angular
grid can be made fine enough that no supersampling artifacts arise.
Only the observer-side (z >= 0) half is painted, matching the pipeline.

Outputs: plots/sim_lamp_maps_h00.png ... h03.png (+ printed totals)
Run:  python3 transient/plot_sim_lamp_height_maps.py [dump.athdf]
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from transient.sim_spectrum import (  # noqa: E402
    load_sim, to_physical, add_colorbar, mark_star, CONFIG, PLOTDIR, LD_CM)
from transient.sim_cloudy_spectrum import Q_ION, ALPHA_B  # noqa: E402

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

HEIGHTS = [0.0, 0.1, 0.2, 0.3]
N_ALPHA = 1024                     # lamp-centered polar directions
N_S = 1200                         # steps along each ray
S_MAX = 3.2                        # [r0]
F_HA = 0.45                        # case-B H alpha photons per recomb.
E_HA = 6.626e-27 * 2.998e10 / 6.563e-5   # erg per H alpha photon


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else \
        os.path.join(HERE, 'data', 'sim', 'disk.out1.00012.athdf')
    cfg = CONFIG
    sim = load_sim(path)
    phys = to_physical(sim, cfg)
    r, th, ph = sim['r'], sim['th'], sim['ph']
    rf, thf, phf = sim['rf'], sim['thf'], sim['phf']
    nr, nph = len(r), len(ph)
    r0_cm = cfg['R0_LD'] * LD_CM
    budget = Q_ION / (4 * np.pi)            # photons/s/sr

    s = np.linspace(1e-3, S_MAX, N_S)
    ds_cm = (s[1] - s[0]) * r0_cm
    alpha = np.linspace(0.02, np.pi / 2 + 0.45, N_ALPHA)
    dalpha = alpha[1] - alpha[0]
    dphi = np.diff(phf)
    dOm = np.sin(alpha) * dalpha            # per unit dphi

    Rp0 = s[None, :] * np.sin(alpha)[:, None]        # (na, ns), h-indep.
    os.makedirs(PLOTDIR, exist_ok=True)
    for h in HEIGHTS:
        Zp = h + s[None, :] * np.cos(alpha)[:, None]
        rp = np.hypot(Rp0, Zp)
        tp = np.arctan2(Rp0, Zp)
        inside = ((rp >= rf[0]) & (rp <= rf[-1])
                  & (tp >= thf[0]) & (tp <= thf[-1]) & (Zp >= 0.0))
        it_ = np.clip(np.searchsorted(th, tp), 0, len(th) - 1)
        ir_ = np.clip(np.searchsorted(r, rp), 0, nr - 1)
        irR = np.clip(np.searchsorted(r, Rp0), 0, nr - 1)  # midplane bin
        lmap = np.zeros((nph, nr))
        for i in range(nph):
            n = np.where(inside, phys['nH'][i][it_, ir_], 0.0)
            dS = n ** 2 * ALPHA_B * (s[None, :] * r0_cm) ** 2 * ds_cm
            S = np.cumsum(dS, axis=1)
            dabs = np.diff(np.clip(S, 0.0, budget), axis=1,
                           prepend=0.0)                 # photons/s/sr
            Lha = dabs * dOm[:, None] * dphi[i] * F_HA * E_HA
            np.add.at(lmap[i], irR.ravel(), Lha.ravel())
        print(f"h = {h:.1f} r0: L(Halpha, recombination count) = "
              f"{lmap.sum():.2e} erg/s")

        R, P = np.meshgrid(r, ph)
        X, Y = R * np.cos(P), R * np.sin(P)
        vmax = lmap.max()
        fig, ax = plt.subplots(figsize=(9.5, 8))
        pc = ax.pcolormesh(X, Y, np.clip(lmap, vmax / 3e3, None),
                           norm=LogNorm(vmin=vmax / 3e3, vmax=vmax),
                           cmap='magma', shading='auto')
        ax.set_aspect('equal')
        add_colorbar(pc, ax, r'$L_{\rm H\alpha}~\rm per~cell~[erg~s^{-1}]$',
                     fontsize=14)
        mark_star(ax, sim)
        ax.set_xlabel(r'$x~[r_0]$', fontsize=16)
        ax.set_ylabel(r'$y~[r_0]$', fontsize=16)
        _m, _e = f'{lmap.sum():.1e}'.split('e')
        ax.text(0.02, 0.02,
                rf'$h_{{\rm lamp}} = {h:.1f}~r_0$' + '\n'
                + rf'$L_{{\rm H\alpha}} = {_m}\times10^{{{int(_e)}}}'
                + r'~\rm erg~s^{-1}$',
                transform=ax.transAxes, va='bottom', fontsize=15,
                color='white',
                bbox=dict(facecolor='black', alpha=0.55,
                          edgecolor='none', pad=4))
        ax.tick_params(which='major', direction='in', length=8, width=1.5,
                       top=True, right=True, labelsize=13)
        ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                       top=True, right=True)
        ax.minorticks_on()
        out = os.path.join(PLOTDIR,
                           f'sim_lamp_maps_h{h:.1f}'.replace('.', '')
                           + '.png')
        plt.savefig(out, dpi=200, bbox_inches='tight')
        plt.close()
        print(f"Saved {out}")


if __name__ == '__main__':
    main()
