#!/usr/bin/env python3
"""plot_sim_lamp_height_rgb.py - face-on three-line RGB composites
(R = H alpha, G = Mg II, B = C IV) for elevated lamp heights.

Full Cloudy treatment on lamp-centered straight rays: for each height,
each (alpha, phi) ray is marched from the lamp at (0, 0, h), its slab
parameters are built exactly as in the main pipeline (foreshortened
ionizing flux at the front from the front-surface tilt in lamp
coordinates, absorption-weighted density, true or grid-max column), the
per-line luminosities come from the arm_column2 grid (window-integrated
line channel, both faces summed - the total emitted, with no viewing
weights), and each line is painted along the ray with its Cloudy G(xi)
emissivity curve. The lamp-centered angular grid (N_ALPHA directions)
is much finer than the sim's theta grid, so no supersampling artifacts.

Outputs: plots/sim_lamp_rgb_h00.png ... h03.png (+ printed line totals)
Run:  python3 transient/plot_sim_lamp_height_rgb.py [dump.athdf]
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from transient.sim_spectrum import (  # noqa: E402
    load_sim, to_physical, mark_star, CONFIG, PLOTDIR, LD_CM)
from transient.sim_cloudy_spectrum import (  # noqa: E402
    load_grid, Q_ION, ALPHA_B)
from transient.sim_cloudy_linemaps import LINES, EMS_NPZ, geval  # noqa: E402

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

HEIGHTS = [0.0, 0.1, 0.2, 0.3]
N_ALPHA = 1024
N_S = 1200
S_MAX = 3.2


def line_tables(grid):
    """Window-integrated line luminosity per cm^2 of slab (both faces)
    for each grid model, as log-space trilinear interpolators."""
    wave, itp, (pv, hv, cv) = grid
    d = np.load(os.path.join(HERE, 'data',
                             'cloudy_arm_spectra_arm_column2.npz'))
    tot_line = d['refl_line'] + d['out_line']          # (npts, nw) nuF_nu
    p, hh, c = d['phi'], d['hden'], d['colden']
    tabs = {}
    for key, (tex, w1, w2) in LINES.items():
        m = (wave >= w1) & (wave <= w2)
        F = np.trapezoid(tot_line[:, m] / wave[None, m], wave[m], axis=1)
        arr = np.zeros((pv.size, hv.size, cv.size))
        arr[np.searchsorted(pv, p), np.searchsorted(hv, hh),
            np.searchsorted(cv, c)] = F
        tabs[key] = RegularGridInterpolator(
            (pv, hv, cv), np.log10(arr + 1e-30))
    return tabs, (pv, hv, cv)


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else \
        os.path.join(HERE, 'data', 'sim', 'disk.out1.00012.athdf')
    cfg = CONFIG
    sim = load_sim(path)
    phys = to_physical(sim, cfg)
    grid = load_grid()
    tabs, (pv, hv, cv) = line_tables(grid)
    ems = np.load(EMS_NPZ)
    pu, hu, cu = (np.unique(ems['phi']), np.unique(ems['hden']),
                  np.unique(ems['colden']))
    M = np.zeros((pu.size, hu.size, cu.size), dtype=int)
    M[np.searchsorted(pu, ems['phi']), np.searchsorted(hu, ems['hden']),
      np.searchsorted(cu, ems['colden'])] = np.arange(ems['phi'].size)
    nxi = ems['xi'].size
    Gc = {key: ems[f'G_{key}'] for key in LINES}

    r, th, ph = sim['r'], sim['th'], sim['ph']
    rf, thf, phf = sim['rf'], sim['thf'], sim['phf']
    nr, nph = len(r), len(ph)
    r0_cm = cfg['R0_LD'] * LD_CM
    budget = Q_ION / (4 * np.pi)

    s = np.linspace(1e-3, S_MAX, N_S)
    ds_cm = (s[1] - s[0]) * r0_cm
    alpha = np.linspace(0.02, np.pi / 2 + 0.45, N_ALPHA)
    dOm_a = np.sin(alpha) * (alpha[1] - alpha[0])
    dphi = np.diff(phf)
    Rp0 = s[None, :] * np.sin(alpha)[:, None]

    for h in HEIGHTS:
        Zp = h + s[None, :] * np.cos(alpha)[:, None]
        rp = np.hypot(Rp0, Zp)
        tp = np.arctan2(Rp0, Zp)
        inside = ((rp >= rf[0]) & (rp <= rf[-1])
                  & (tp >= thf[0]) & (tp <= thf[-1]) & (Zp >= 0.0))
        it_ = np.clip(np.searchsorted(th, tp), 0, len(th) - 1)
        ir_ = np.clip(np.searchsorted(r, rp), 0, nr - 1)
        irR = np.clip(np.searchsorted(r, Rp0), 0, nr - 1)

        # ---- pass 1: march every azimuth, store per-ray properties ----
        # stored NORMALIZED by the photon budget: raw cumulative
        # recombination values (~1e53) overflow float32
        S_all = np.zeros((nph, N_ALPHA, N_S), dtype=np.float32)
        s_eff = np.full((nph, N_ALPHA), s[-1])
        matter = np.ones((nph, N_ALPHA), dtype=bool)
        logn = np.full((nph, N_ALPHA), hv[0])
        logN = np.full((nph, N_ALPHA), cv[0])
        hasgas = np.zeros((nph, N_ALPHA), dtype=bool)
        for i in range(nph):
            n = np.where(inside, phys['nH'][i][it_, ir_], 0.0)
            dS = n ** 2 * ALPHA_B * (s[None, :] * r0_cm) ** 2 * ds_cm
            S = np.cumsum(dS, axis=1)
            S_all[i] = (S / budget).astype(np.float32)
            hasgas[i] = S[:, -1] > 0
            hit = S >= budget
            has = hit.any(axis=1)
            matter[i] = ~has
            k = np.argmax(hit, axis=1)
            S_hi = S[np.arange(N_ALPHA), k]
            S_lo = np.where(k > 0,
                            S[np.arange(N_ALPHA), np.maximum(k - 1, 0)], 0)
            fc = np.clip((budget - S_lo) / (S_hi - S_lo + 1e-300), 0, 1)
            sif = s[np.maximum(k - 1, 0)] + fc * (s[1] - s[0])
            # matter rays: last in-gas point
            klast = N_S - 1 - np.argmax((n > 0)[:, ::-1], axis=1)
            s_eff[i] = np.where(has, sif, s[klast])
            dabs = np.diff(np.clip(S, 0, budget), axis=1, prepend=0.0)
            wsum = dabs.sum(axis=1)
            ok = wsum > 0
            ln = np.where(n > 0, np.log10(n + 1e-30), 0.0)
            logn[i][ok] = ((dabs * ln).sum(axis=1)[ok] / wsum[ok])
            Ncol = (n * ds_cm).sum(axis=1)
            logN[i] = np.where(has, cv[-1], np.log10(Ncol + 1.0))

        # ---- foreshortening from the front surface in lamp coords ----
        lnf = np.log(s_eff)
        dfa = np.gradient(lnf, alpha, axis=1)
        dfp = np.gradient(lnf, ph, axis=0)
        mu = 1.0 / np.sqrt(1.0 + dfa ** 2
                           + (dfp / np.sin(alpha)[None, :]) ** 2)
        logphi = np.log10(Q_ION / (4 * np.pi)
                          / (s_eff * r0_cm) ** 2 * mu + 1e-30)
        pts = np.stack([np.clip(logphi, pv[0], pv[-1]),
                        np.clip(logn, hv[0], hv[-1]),
                        np.clip(logN, cv[0], cv[-1])], axis=-1)
        dA = (s_eff * r0_cm) ** 2 * dOm_a[None, :] * dphi[:, None] / mu
        midx = M[np.argmin(np.abs(pts[..., 0:1] - pu[None, None, :]), -1),
                 np.argmin(np.abs(pts[..., 1:2] - hu[None, None, :]), -1),
                 np.argmin(np.abs(pts[..., 2:3] - cu[None, None, :]), -1)]
        L_ray = {}
        for key in LINES:
            L = 10.0 ** tabs[key](pts.reshape(-1, 3)).reshape(nph, N_ALPHA)
            L_ray[key] = np.where(hasgas, L * dA, 0.0)

        # ---- pass 2: paint with the Cloudy G(xi) curves ----
        maps = {key: np.zeros((nph, nr)) for key in LINES}
        for i in range(nph):
            S = S_all[i].astype(np.float64)          # in budget units
            Send = np.where(matter[i], S[:, -1], 1.0)
            Send = np.maximum(Send, 1e-300)[:, None]
            xi = np.clip(S / Send, 0.0, 1.0)
            for key in LINES:
                cf = Gc[key][midx[i]]                  # (na, nxi)
                Gv = geval(cf, xi, nxi)
                dfrac = np.maximum(np.diff(Gv, axis=1, prepend=0.0), 0.0)
                np.add.at(maps[key][i], irR.ravel(),
                          (L_ray[key][i][:, None] * dfrac).ravel())
        for key in LINES:
            print(f"h = {h:.1f}: {key:7s} L = {maps[key].sum():.2e} erg/s")

        # ---- RGB composite (same recipe as sim_linemap_rgb) ----
        ngrid = 900
        xg = np.linspace(-r.max(), r.max(), ngrid)
        Xg, Yg = np.meshgrid(xg, xg)
        Rg = np.hypot(Xg, Yg)
        Pg = np.mod(np.arctan2(Yg, Xg), 2 * np.pi)
        ins = (Rg >= r.min()) & (Rg <= r.max())
        fr = np.clip(np.interp(Rg, r, np.arange(nr)), 0, nr - 1 - 1e-9)
        i0 = fr.astype(int)
        i1 = np.minimum(i0 + 1, nr - 1)
        fs = fr - i0
        phe = np.concatenate([[ph[-1] - 2 * np.pi], ph,
                              [ph[0] + 2 * np.pi]])
        fp = np.interp(Pg, phe, np.arange(-1, nph + 1))
        j0 = np.floor(fp).astype(int) % nph
        j1 = (j0 + 1) % nph
        ft = fp - np.floor(fp)
        DR, SAT = 2.5, 5.0
        rgb = np.zeros((ngrid, ngrid, 3))
        for k, key in enumerate(['halpha', 'mgii', 'civ']):
            logl = np.log10(maps[key] + 1e-30)
            lit = logl > logl.max() - DR
            qs = np.linspace(0.0, 100.0, 13)
            anchors = np.percentile(logl[lit], qs)
            anchors = np.maximum.accumulate(anchors)
            anchors += 1e-9 * np.arange(anchors.size)
            eq = np.interp(logl, anchors, qs / 100.0)
            eq[~lit] = 0.0
            samp = ((1 - ft) * (1 - fs) * eq[j0, i0]
                    + ft * (1 - fs) * eq[j1, i0]
                    + (1 - ft) * fs * eq[j0, i1]
                    + ft * fs * eq[j1, i1])
            samp[~ins] = 0.0
            rgb[..., k] = samp
        lum = rgb.mean(axis=2, keepdims=True)
        rgb = np.clip(lum + SAT * (rgb - lum), 0, 1)

        fig, ax = plt.subplots(figsize=(9.5, 9))
        ax.imshow(rgb, extent=[-r.max(), r.max(), -r.max(), r.max()],
                  origin='lower', interpolation='bilinear')
        ax.set_aspect('equal')
        mark_star(ax, sim)
        for kk, (txt, col) in enumerate([
                (r'$\rm H\alpha$', 'crimson'),
                (r'$\rm Mg\,II$', 'mediumseagreen'),
                (r'$\rm C\,IV$', 'cornflowerblue')]):
            ax.text(0.02, 0.97 - 0.055 * kk, txt, transform=ax.transAxes,
                    va='top', fontsize=17, color=col)
        ax.text(0.98, 0.03, rf'$h_{{\rm lamp}} = {h:.1f}~r_0$',
                transform=ax.transAxes, ha='right', va='bottom',
                fontsize=16, color='white')
        ax.set_xlabel(r'$x~[r_0]$', fontsize=16)
        ax.set_ylabel(r'$y~[r_0]$', fontsize=16)
        ax.tick_params(which='major', direction='in', length=8, width=1.5,
                       top=True, right=True, labelsize=13)
        ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                       top=True, right=True)
        ax.minorticks_on()
        out = os.path.join(PLOTDIR,
                           f'sim_lamp_rgb_h{h:.1f}'.replace('.', '')
                           + '.png')
        plt.savefig(out, dpi=200, bbox_inches='tight')
        plt.close()
        print(f"Saved {out}")


if __name__ == '__main__':
    main()
