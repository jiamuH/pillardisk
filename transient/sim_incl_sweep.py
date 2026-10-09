#!/usr/bin/env python3
"""
sim_incl_sweep.py - the observed spectrum versus viewing inclination.

The Strommgren march, shadows, foreshortening and Cloudy slab lookups are
set by the LAMP's geometry and are inclination-independent (computed
once); per inclination only the observer-dependent steps are redone:
face selection (lit vs shadow side of each slab) and the Doppler
projection of the slab velocities. Observer-side blocking of one slab by
another is NOT modeled, so i is capped at 80 deg (edge-on would need it).

Outputs:
  sim_incl_spectra.png    - spectra at all inclinations (plasma colors)
  sim_incl_linecurves.png - C IV / Mg II / H alpha flux vs inclination
  sim_incl_trailed.png    - H alpha profile vs inclination (horn
                            separation growing as v sin i)
  data/sim_incl_spectra.npz - cached arrays

Run:  python3 transient/sim_incl_sweep.py [dump.athdf]
"""

import os
import sys

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, NullFormatter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from transient.sim_spectrum import (  # noqa: E402
    load_sim, to_physical, vlos_cube, add_colorbar, CONFIG, C_KMS, PLOTDIR)
from transient.sim_cloudy_spectrum import load_grid, compute_rays  # noqa: E402

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

INCLS = np.arange(0.0, 80.1, 5.0)
LINES = {'civ': (r'$\rm C\,IV$', 1500., 1600., 'royalblue'),
         'mgii': (r'$\rm Mg\,II$', 2750., 2850., 'seagreen'),
         'halpha': (r'$\rm H\alpha$', 6470., 6660., 'crimson')}


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else \
        os.path.join(HERE, 'data', 'sim', 'disk.out1.00012.athdf')
    cfg = CONFIG
    sim = load_sim(path)
    phys = to_physical(sim, cfg)
    grid = load_grid()
    rays = compute_rays(sim, phys, cfg, grid)
    wave = grid[0]
    from transient.sim_full_spectrum import thermal_and_agn
    L_th_base, L_agn = thermal_and_agn(sim, phys, cfg, wave)

    F = rays['F_raw']
    nx, ny, nz = rays['normal']
    matter = rays['matter']
    dA = rays['dA']
    w = rays['w']
    wsum = rays['wsum']

    spectra = np.zeros((len(INCLS), wave.size))
    line_spec = np.zeros((len(INCLS), wave.size))
    for k, inc in enumerate(INCLS):
        ir = np.radians(inc)
        si, ci = np.sin(ir), np.cos(ir)
        cosv = nx * si + nz * ci                     # observer at phi_obs=0
        s_out = 0.5 * (1.0 + np.tanh(cosv / 0.3))
        # Lambert projection with the opacity blend: opaque faces get
        # face selection x 4|cos v|; transparent slabs blend toward
        # isotropic by f_used (fraction of the photon budget consumed),
        # so the matter/radiation-bounded flip is continuous.
        proj_opq = 4.0 * np.abs(cosv)
        blend = np.where(matter, rays['f_used'], 1.0)
        w_out = (blend * s_out * proj_opq + (1.0 - blend)).reshape(-1)
        w_ref = (blend * (1.0 - s_out) * proj_opq
                 + (1.0 - blend)).reshape(-1)
        Fc = w_ref[:, None] * F['refl_cont'] + w_out[:, None] * F['out_cont']
        Fl = w_ref[:, None] * F['refl_line'] + w_out[:, None] * F['out_line']
        _, (r_o, t_o, p_o) = vlos_cube(sim, inc, 0.0)
        vlos = phys['vr'] * r_o + phys['vt'] * t_o + phys['vp'] * p_o
        v_ray = np.where(wsum > 0,
                         (w * vlos).sum(axis=2) / (wsum + 1e-300),
                         vlos[:, :, 0]).reshape(-1)
        shift = 1.0 - v_ray / C_KMS
        Lt = np.zeros_like(wave)
        Ll = np.zeros_like(wave)
        for i in range(len(dA)):
            ws = wave * shift[i]
            Lt += dA[i] * np.interp(wave, ws, (Fc[i] + Fl[i]) / ws,
                                    left=0, right=0)
            Ll += dA[i] * np.interp(wave, ws, Fl[i] / ws, left=0, right=0)
        # total = reprocessed + thermal (Lambert cos i) + direct AGN
        spectra[k] = Lt + 4.0 * ci * L_th_base + L_agn
        line_spec[k] = Ll
        print(f"i = {inc:4.0f} deg done")

    np.savez(os.path.join(HERE, 'data', 'sim_incl_spectra.npz'),
             incls=INCLS, wave=wave, spectra=spectra, line_spec=line_spec)
    os.makedirs(PLOTDIR, exist_ok=True)
    cmap = mpl.cm.plasma
    norm = mpl.colors.Normalize(0, 80)

    # ---- 1. spectra at all inclinations ----
    fig, ax = plt.subplots(figsize=(12, 7))
    show = [0., 20., 40., 60., 80.]          # subset; full grid in the npz
    for k, inc in enumerate(INCLS):
        if inc not in show:
            continue
        ax.plot(wave, wave * spectra[k], '-', color=cmap(norm(inc)),
                lw=1.6, alpha=0.85, drawstyle='steps-mid')
    ax.axvline(3646, color='gray', ls=':', lw=1.2)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(1000, 9000)
    msel = (wave >= 1000) & (wave <= 9000)
    nuLnu = wave[None, :] * spectra
    ax.set_ylim(nuLnu[:, msel].min() * 0.5, nuLnu[:, msel].max() * 2)
    ax.set_xticks([1000, 2000, 3000, 5000, 9000])
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel(r'$\rm wavelength~[\AA]$', fontsize=18)
    ax.set_ylabel(r'$\lambda L_\lambda~[\rm erg~s^{-1}]$', fontsize=17)
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = add_colorbar(sm, ax, r'$\rm inclination~[deg]$')
    cb.set_ticks([0, 20, 40, 60, 80])
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=14)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()
    out = os.path.join(PLOTDIR, 'sim_incl_spectra.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")

    # ---- 1b. UV (top) + optical (bottom) zooms ----
    fig, (axU, axO) = plt.subplots(2, 1, figsize=(12, 12.5))
    for axz, (w1, w2), ticks in [
            (axU, (1000., 3200.), [1000, 1500, 2000, 3000]),
            (axO, (3800., 7200.), [4000, 5000, 6000, 7000])]:
        for k, inc in enumerate(INCLS):
            if inc not in show:
                continue
            axz.plot(wave, wave * spectra[k], '-', color=cmap(norm(inc)),
                     lw=3.0, alpha=0.65, drawstyle='steps-mid',
                     solid_joinstyle='round')
        axz.set_xscale('log')
        axz.set_yscale('log')
        axz.set_xlim(w1, w2)
        mz = (wave >= w1) & (wave <= w2)
        axz.set_ylim(nuLnu[:, mz].min() * 0.7, nuLnu[:, mz].max() * 1.5)
        axz.set_xticks(ticks)
        axz.xaxis.set_major_formatter(ScalarFormatter())
        axz.xaxis.set_minor_formatter(NullFormatter())
        axz.set_ylabel(r'$\lambda L_\lambda~[\rm erg~s^{-1}]$',
                       fontsize=17)
        smz = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        cbz = add_colorbar(smz, axz, r'$\rm inclination~[deg]$')
        cbz.set_ticks([0, 20, 40, 60, 80])
        axz.tick_params(which='major', direction='in', length=8, width=1.5,
                        top=True, right=True, labelsize=14)
        axz.tick_params(which='minor', direction='in', length=4, width=1.0,
                        top=True, right=True)
        axz.minorticks_on()
    axO.set_xlabel(r'$\rm wavelength~[\AA]$', fontsize=18)
    plt.tight_layout()
    out = os.path.join(PLOTDIR, 'sim_incl_spectra_zoom.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")

    # ---- 2. line fluxes vs inclination ----
    fig, ax = plt.subplots(figsize=(11, 7))
    for key, (tex, w1, w2, col) in LINES.items():
        m = (wave >= w1) & (wave <= w2)
        lf = np.trapz(line_spec[:, m], wave[m], axis=1)
        ax.plot(INCLS, lf / lf[0], '-o', color=col, lw=2.5, ms=6,
                alpha=0.85, label=tex)
        print(f"{key:7s}: L(i=0) = {lf[0]:.2e}, L(80)/L(0) = "
              f"{lf[-1]/lf[0]:.2f}")
    ax.set_xlabel(r'$\rm inclination~[deg]$', fontsize=18)
    ax.set_ylabel(r'$L_{\rm line}(i)/L_{\rm line}(0)$', fontsize=17)
    ax.axhline(1, color='gray', lw=1)
    ax.set_xlim(0, 80)
    ax.legend(fontsize=14, frameon=False)
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=14)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()
    out = os.path.join(PLOTDIR, 'sim_incl_linecurves.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")

    # ---- 3. H alpha profile vs inclination ----
    m = (wave >= 6350) & (wave <= 6800)
    vv = (wave[m] / 6563.0 - 1.0) * C_KMS
    fig, ax = plt.subplots(figsize=(10, 7))
    pc = ax.pcolormesh(vv, INCLS, spectra[:, m], cmap='inferno',
                       shading='auto')
    add_colorbar(pc, ax, r'$L_\lambda~[\rm erg~s^{-1}~\AA^{-1}]$')
    ax.set_xlabel(r'$\rm velocity~[km~s^{-1}]$', fontsize=17)
    ax.set_ylabel(r'$\rm inclination~[deg]$', fontsize=17)
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=14)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()
    out = os.path.join(PLOTDIR, 'sim_incl_trailed.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")


if __name__ == '__main__':
    main()
