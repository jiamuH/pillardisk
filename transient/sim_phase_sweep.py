#!/usr/bin/env python3
"""
sim_phase_sweep.py - phase-resolved spectral prediction: the observed
spectrum as the spiral pattern rotates past the observer.

Sweeping the observer azimuth phi_obs (equivalent to the pattern rotating
by -phi_obs, i.e. one full orbital phase cycle) at fixed inclination. The
Strommgren march and Cloudy slab lookup are phase-INDEPENDENT and done
once; per phase only the face selection (lit vs shadow side of each slab)
and the Doppler projection change.

Outputs:
  sim_phase_spectra.png    - full spectra at all phases (cyclic colors)
  sim_phase_linecurves.png - C IV / Mg II / H alpha flux vs orbital phase
  sim_phase_trailed.png    - trailed spectrogram of the H alpha region

Run:  python3 transient/sim_phase_sweep.py [dump.athdf]
"""

import os
import sys

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from transient.sim_spectrum import (  # noqa: E402
    load_sim, to_physical, vlos_cube, add_colorbar, CONFIG, C_KMS, PLOTDIR)
from transient.sim_cloudy_spectrum import load_grid, compute_rays  # noqa: E402

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

N_PHASE = 24
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
    nx, ny, nz = rays['normal']                       # (ph, th)
    matter = rays['matter']
    dA = rays['dA']                                   # (nray,)
    w = rays['w']
    wsum = rays['wsum']
    i_obs = np.radians(cfg['INCL_DEG'])
    si, ci = np.sin(i_obs), np.cos(i_obs)

    phases = np.linspace(0.0, 360.0, N_PHASE, endpoint=False)
    spectra = np.zeros((N_PHASE, wave.size))
    line_spec = np.zeros((N_PHASE, wave.size))
    for k, po in enumerate(phases):
        pr = np.radians(po)
        ox, oy = si * np.cos(pr), si * np.sin(pr)
        cosv = nx * ox + ny * oy + nz * ci
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
        _, (r_o, t_o, p_o) = vlos_cube(sim, cfg['INCL_DEG'], po)
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
        # total = reprocessed + thermal (phase-independent: flat-surface
        # Lambert factor has no azimuth dependence) + direct AGN
        spectra[k] = Lt + 4.0 * ci * L_th_base + L_agn
        line_spec[k] = Ll
        print(f"phase {po:5.1f} deg done")

    np.savez(os.path.join(HERE, 'data', 'sim_phase_spectra.npz'),
             phases=phases, wave=wave, spectra=spectra,
             line_spec=line_spec, incl=cfg['INCL_DEG'])
    os.makedirs(PLOTDIR, exist_ok=True)
    cmap = mpl.cm.twilight

    # ---- 1. spectra at all phases ----
    fig, ax = plt.subplots(figsize=(12, 7))
    for k, po in enumerate(phases):
        ax.plot(wave, wave * spectra[k], '-', color=cmap(po / 360.0),
                lw=1.2, alpha=0.75, drawstyle='steps-mid')
    ax.axvline(3646, color='gray', ls=':', lw=1.2)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(1000, 9000)
    msel = (wave >= 1000) & (wave <= 9000)
    nuLnu = wave[None, :] * spectra
    ax.set_ylim(nuLnu[:, msel].min() * 0.5, nuLnu[:, msel].max() * 2)
    from matplotlib.ticker import ScalarFormatter, NullFormatter
    ax.set_xticks([1000, 2000, 3000, 5000, 9000])
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel(r'$\rm wavelength~[\AA]$', fontsize=18)
    ax.set_ylabel(r'$\lambda L_\lambda~[\rm erg~s^{-1}]$',
                  fontsize=17)
    sm = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(0, 360),
                               cmap=cmap)
    cb = add_colorbar(sm, ax, r'$\rm orbital~phase~[deg]$')
    cb.set_ticks([0, 90, 180, 270, 360])
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=14)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()
    out = os.path.join(PLOTDIR, 'sim_phase_spectra.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")

    # ---- 2. line fluxes vs phase ----
    fig, ax = plt.subplots(figsize=(11, 7))
    for key, (tex, w1, w2, col) in LINES.items():
        m = (wave >= w1) & (wave <= w2)
        lf = np.trapz(line_spec[:, m], wave[m], axis=1)
        ax.plot(phases, lf / lf.mean(), '-o', color=col, lw=2.5, ms=6,
                alpha=0.85, label=tex)
        print(f"{key:7s}: mean L = {lf.mean():.2e} erg/s, "
              f"modulation (max-min)/mean = {(lf.max()-lf.min())/lf.mean():.2f}")
    ax.set_xlabel(r'$\rm orbital~phase~[deg]$', fontsize=18)
    ax.set_ylabel(r'$L_{\rm line}/\langle L_{\rm line}\rangle$', fontsize=17)
    ax.axhline(1, color='gray', lw=1)
    ax.set_xlim(0, 360)
    ax.legend(fontsize=14, frameon=False)
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=14)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()
    out = os.path.join(PLOTDIR, 'sim_phase_linecurves.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")

    # ---- 3. trailed spectrogram around H alpha ----
    m = (wave >= 6350) & (wave <= 6800)
    vv = (wave[m] / 6563.0 - 1.0) * C_KMS
    fig, ax = plt.subplots(figsize=(10, 7))
    pc = ax.pcolormesh(vv, phases, spectra[:, m], cmap='inferno',
                       shading='auto')
    add_colorbar(pc, ax, r'$L_\lambda~[\rm erg~s^{-1}~\AA^{-1}]$')
    ax.set_xlabel(r'$\rm velocity~[km~s^{-1}]$', fontsize=17)
    ax.set_ylabel(r'$\rm orbital~phase~[deg]$', fontsize=17)
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=14)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()
    out = os.path.join(PLOTDIR, 'sim_phase_trailed.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")


if __name__ == '__main__':
    main()
