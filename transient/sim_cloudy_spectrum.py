#!/usr/bin/env python3
"""
sim_cloudy_spectrum.py - Cloudy-based spectral prediction from the Athena++
disk simulation (stage 2 of the sim-to-spectrum pipeline).

Design: the central source (AGN, at the origin) ionizes the wedge from
inside. Each (theta, phi) direction is ONE radial ray = ONE Cloudy slab:

  1. Strommgren march along the ray: the per-steradian photon budget
     Q/4pi is eaten by recombinations n^2 alpha_B r^2 dr; the ionization
     front r_IF is where it is exhausted (radiation-bounded), else the ray
     is matter-bounded through the wedge.
  2. The ray's slab parameters: log phi at the illuminated face,
     emission-weighted log n_H of the ionized segment, and the ionized
     column log N_H (capped at the grid edge). The cached `arm_column`
     Cloudy grid (phi 17-21, hden 9-12, colden 21.5-23.5; reflected +
     outward, continuum + lines, nuF_nu per cm^2 of slab face) is
     interpolated there.
  3. The slab spectrum is weighted by the ray's illuminated-face area
     d^2 dOmega; the LINE spectrum is Doppler-shifted by the ray's
     emission-weighted line-of-sight velocity (observer at INCL_DEG).
  4. Sum over the 64 x 256 rays -> L_lambda with kinematic line profiles
     and photon conservation with the source built in.

Physical anchors are shared with sim_spectrum.py (r0, M_BH, n_mid).
The ionizing source is Q_ION [photons/s]; the default puts the face flux
at log phi = 21 (grid maximum).

Run:  python3 transient/sim_cloudy_spectrum.py [dump.athdf]
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from transient.sim_spectrum import (  # noqa: E402
    load_sim, to_physical, vlos_cube, CONFIG, LD_CM, C_KMS, PLOTDIR)

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

NPZ = os.path.join(HERE, 'data', 'cloudy_arm_spectra_arm_column.npz')
_NPZ2 = os.path.join(HERE, 'data', 'cloudy_arm_spectra_arm_column2.npz')
if os.path.exists(_NPZ2):
    NPZ = _NPZ2          # extended grid (phi 17-23, hden 6-12, N 21.5-24.5)
ALPHA_B = 2.59e-13          # case B recombination coefficient (1e4 K)
Q_ION = 8.5e52              # ionizing photon rate [1/s] -> log phi_face = 21


def load_grid():
    """separate interpolators for the reflected (illuminated-face) and
    outward (shielded back-face) components, continuum and lines."""
    d = np.load(NPZ)
    wave = d['wave']
    p, h, c = d['phi'], d['hden'], d['colden']
    pv, hv, cv = np.unique(p), np.unique(h), np.unique(c)
    i = np.searchsorted(pv, p)
    j = np.searchsorted(hv, h)
    m = np.searchsorted(cv, c)
    itp = {}
    for k in ('refl_cont', 'refl_line', 'out_cont', 'out_line'):
        arr = np.zeros((pv.size, hv.size, cv.size, wave.size))
        arr[i, j, m, :] = d[k]
        itp[k] = RegularGridInterpolator((pv, hv, cv),
                                         np.log10(arr + 1e-30))
    return wave, itp, (pv, hv, cv)


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else \
        os.path.join(HERE, 'data', 'sim', 'disk.out1.00012.athdf')
    cfg = CONFIG
    sim = load_sim(path)
    phys = to_physical(sim, cfg)
    grid = load_grid()
    rays = compute_rays(sim, phys, cfg, grid)
    wave = grid[0]
    F_cont, F_line = rays['F_cont'], rays['F_line']
    dA, v_ray, r_if = rays['dA'], rays['v_ray'], rays['r_if']
    th, ph = sim['th'], sim['ph']
    logphi = rays['logphi']
    make_figures(wave, F_cont, F_line, dA, v_ray, r_if, th, ph, logphi, cfg)


def compute_rays(sim, phys, cfg, grid):
    """Strommgren ray-march + Cloudy slab lookup for every (theta, phi)
    ray. Returns per-ray spectra, geometry, and the along-ray emission
    weights (for distributing a ray's emission back onto the disk).

    FACE SELECTION: each radiation-bounded slab shows the observer either
    its illuminated face (-> reflected spectrum) or its shielded back
    (-> outward spectrum), decided by the front-surface normal dotted with
    the observer direction (smooth tanh blend near edge-on). Non-resonant
    photons escape the back freely; resonant lines favor the lit face -
    Cloudy's refl/out split encodes that. Matter-bounded (transparent)
    rays emit both components regardless of viewing."""
    wave, itp, (pv, hv, cv) = grid

    r0 = cfg['R0_LD'] * LD_CM
    r, th, ph = sim['r'], sim['th'], sim['ph']
    nH = phys['nH']                                   # (ph, th, r)
    dr = np.diff(sim['rf']) * r0                      # cm
    r_cm = r * r0

    # ---- Strommgren march along each ray ----
    dN = nH * dr[None, None, :]                       # column per cell
    dS = nH ** 2 * ALPHA_B * r_cm ** 2 * dr[None, None, :]  # recomb/sr
    S = np.cumsum(dS, axis=2)
    budget = Q_ION / (4 * np.pi)
    ion = S <= budget                                  # ionized cells
    N_ion = np.where(ion, dN, 0.0).sum(axis=2)         # ionized column
    frac_matter = 1.0 - ion[:, :, -1].mean()
    # SUB-CELL front radius: interpolate where the cumulative
    # recombination S crosses the photon budget inside the front cell.
    # Snapping r_IF to whole cells makes the surface a staircase, and its
    # numerical gradients (the foreshortening mu, the face normal) then
    # flicker between adjacent columns - that quantization noise printed
    # radial 'spoke' streaks into the maps and weight jitter into the
    # spectra.
    idx_if = np.argmin(ion, axis=2)                # first shadowed cell
    S_hi = np.take_along_axis(S, idx_if[..., None], axis=2)[..., 0]
    S_lo = np.where(
        idx_if > 0,
        np.take_along_axis(S, np.maximum(idx_if - 1, 0)[..., None],
                           axis=2)[..., 0], 0.0)
    fcell = np.clip((budget - S_lo) / (S_hi - S_lo + 1e-300), 0.0, 1.0)
    rf_e = sim['rf']                               # cell faces, code units
    r_if_cont = rf_e[idx_if] + fcell * (rf_e[idx_if + 1] - rf_e[idx_if])
    r_if = np.where(ion[:, :, -1], rf_e[-1], r_if_cont)

    # ---- per-ray slab parameters ----
    # Lamp-post irradiation geometry: the slab's ionizing flux is evaluated
    # AT the ionization surface r_IF, including the foreshortening factor
    # mu = cos(incidence) from the surface orientation. A front that faces
    # the lamp (r_IF ~ constant with angle: the inner rim, arm walls) gets
    # mu ~ 1; where the front climbs steeply with angle (the flared disk
    # surface, seen at grazing incidence) mu << 1. Same physics as the
    # pillardisk T_irr foreshortening.
    lnf = np.log(r_if)                                 # (ph, th)
    dlnf_dph = np.gradient(lnf, ph, axis=0)
    dlnf_dth = np.gradient(lnf, th, axis=1)
    mu = 1.0 / np.sqrt(1.0 + dlnf_dth ** 2
                       + (dlnf_dph / np.sin(th)[None, :]) ** 2)
    logphi = np.log10(Q_ION / (4 * np.pi * (r_if * r0) ** 2) * mu)
    w = np.where(ion, dS, 0.0)                         # emission weight
    wsum = w.sum(axis=2)
    logn = np.where(wsum > 0,
                    (w * np.log10(nH)).sum(axis=2) / (wsum + 1e-300),
                    np.log10(nH[:, :, 0]))             # face cell fallback
    # matter-bounded rays: the true (finite) column; radiation-bounded
    # rays: hand Cloudy the grid-maximum column and let its own internal
    # ionization front terminate the slab.
    matter = ion[:, :, -1]
    logN = np.where(matter, np.log10(dN.sum(axis=2) + 1.0), cv[-1])
    print(f"ray log phi (at front, foreshortened): "
          f"{np.percentile(logphi, [5, 50, 95]).round(2)}   "
          f"matter-bounded rays: {100*frac_matter:.1f}%")
    print(f"ray log n_H: {np.percentile(logn, [5, 50, 95]).round(2)}   "
          f"ray log N_ion: {np.percentile(logN, [5, 50, 95]).round(2)}")
    # clip into grid coverage (report)
    pts = np.stack([logphi, logn, logN], axis=-1)
    nclip = ((pts[..., 0] < pv[0]) | (pts[..., 0] > pv[-1])
             | (pts[..., 1] < hv[0]) | (pts[..., 1] > hv[-1])
             | (pts[..., 2] < cv[0]) | (pts[..., 2] > cv[-1])).mean()
    print(f"rays clipped into grid box: {100*nclip:.1f}%")
    pts[..., 0] = np.clip(pts[..., 0], pv[0], pv[-1])
    pts[..., 1] = np.clip(pts[..., 1], hv[0], hv[-1])
    pts[..., 2] = np.clip(pts[..., 2], cv[0], cv[-1])

    flat = pts.reshape(-1, 3)
    F = {k: 10 ** itp[k](flat) for k in itp}           # (nray, nw) nuF_nu

    # ---- face selection: which side of each slab does the observer see?
    # back-face (away-from-lamp) normal of the front surface r = f(th,ph):
    #   n ~ r_hat - (dlnf/dth) th_hat - (dlnf/dph / sin th) ph_hat
    a = dlnf_dth
    b = dlnf_dph / np.sin(th)[None, :]
    st, ct = np.sin(th)[None, :], np.cos(th)[None, :]
    sp, cp = np.sin(ph)[:, None], np.cos(ph)[:, None]
    nx = st * cp - a * ct * cp + b * sp
    ny = st * sp - a * ct * sp - b * cp
    nz = ct + a * st
    nn = np.sqrt(nx ** 2 + ny ** 2 + nz ** 2)
    i_obs = np.radians(cfg['INCL_DEG'])
    cosv = (nx * np.sin(i_obs) + nz * np.cos(i_obs)) / nn
    s_out = 0.5 * (1.0 + np.tanh(cosv / 0.3))          # 1 = back visible
    # Lambert projection: an optically thick face viewed at angle
    # theta from its normal delivers flux ~ cos(theta) (projected
    # area); factor 4 = isotropic-equivalent of a Lambertian face.
    # Matter-bounded (transparent) slabs emit isotropically: factor 1.
    # OPACITY BLEND: the matter/radiation-bounded classification is a
    # hard threshold, but a slab that barely stays transparent is
    # physically almost opaque. Weighting matter-bounded slabs by a fixed
    # isotropic 1 makes each single-ray flip between neighboring columns
    # a factor ~2-3 jump (rectangular 'tab' artifacts in the maps).
    # Blend by f_used, the fraction of the photon budget the ray
    # consumed: f_used -> 1 matches the opaque Lambert limit exactly.
    f_used = np.clip(S[:, :, -1] / budget, 0.0, 1.0)
    proj_opq = 4.0 * np.abs(cosv)
    blend = np.where(matter, f_used, 1.0)
    proj = blend * proj_opq + (1.0 - blend)
    w_out = (blend * s_out * proj_opq + (1.0 - blend)).reshape(-1)
    w_ref = (blend * (1.0 - s_out) * proj_opq + (1.0 - blend)).reshape(-1)
    nrb = (~matter).sum()
    print(f"face selection (i={cfg['INCL_DEG']:.0f}): of "
          f"{nrb} radiation-bounded rays, "
          f"{100*np.mean(s_out[~matter] > 0.5):.0f}% show the shadow face, "
          f"{100*np.mean(s_out[~matter] <= 0.5):.0f}% the lit face")
    F_cont = w_ref[:, None] * F['refl_cont'] + w_out[:, None] * F['out_cont']
    F_line = w_ref[:, None] * F['refl_line'] + w_out[:, None] * F['out_line']

    # ---- ray geometry: face area and emission-weighted v_los ----
    dth = np.diff(sim['thf'])
    dph = np.diff(sim['phf'])
    dOm = dph[:, None] * (np.sin(th) * dth)[None, :]
    # slab surface area per ray: the patch the ray's solid angle subtends
    # at the front, enlarged by 1/mu for the tilt (so flux*area conserves
    # the per-steradian photon budget)
    dA = ((r_if * r0) ** 2 * dOm / mu).reshape(-1)     # cm^2 per ray
    # the observer sees ONE face of the opaque disk: keep only the
    # observer-side (upper, theta <= pi/2) hemisphere in the reprocessed
    # sum, matching the one-sided thermal treatment. The far-side strip
    # of the LOWER rim wall visible through the central funnel is
    # neglected (would need slab-on-slab occultation). Zeroing dA here
    # propagates to the spectra, both sweeps, and the line maps alike.
    upper = np.broadcast_to((th <= np.pi / 2)[None, :],
                            (len(ph), len(th))).reshape(-1)
    dA = np.where(upper, dA, 0.0)
    _, (r_o, t_o, p_o) = vlos_cube(sim, cfg['INCL_DEG'])
    vlos = phys['vr'] * r_o + phys['vt'] * t_o + phys['vp'] * p_o
    v_ray = np.where(wsum > 0,
                     (w * vlos).sum(axis=2) / (wsum + 1e-300),
                     vlos[:, :, 0]).reshape(-1)          # km/s
    return dict(F_cont=F_cont, F_line=F_line, dA=dA, v_ray=v_ray,
                f_used=f_used,
                r_if=r_if, w=w, wsum=wsum, logn=logn, logN=logN,
                logphi=logphi, ion=ion, matter=matter, F_raw=F,
                normal=(nx / nn, ny / nn, nz / nn))


def make_figures(wave, F_cont, F_line, dA, v_ray, r_if, th, ph, logphi,
                 cfg):
    # ---- assemble L_lambda ----
    # nuF_nu per cm^2 -> F_lambda = nuF_nu / lambda; each slab's ENTIRE
    # emitted spectrum (continuum, edges, and lines alike) is Doppler-
    # shifted by its emission-weighted line-of-sight velocity - the
    # rest-frame Cloudy spectrum observed from moving gas. (First-order
    # shift; beaming D^3 ~ 20% on the inner rim and gravitational
    # redshift ~1% are the neglected next-order terms.)
    L_cont = np.zeros_like(wave)
    L_line = np.zeros_like(wave)
    shift = 1.0 - v_ray / C_KMS
    for i in range(len(dA)):
        ws = wave * shift[i]
        L_cont += dA[i] * np.interp(wave, ws, F_cont[i] / ws,
                                    left=0, right=0)
        L_line += dA[i] * np.interp(wave, ws, F_line[i] / ws,
                                    left=0, right=0)
    L_tot = L_cont + L_line
    L_abs = Q_ION * 2.4e-11
    print(f"L(reprocessed, 1200-12000 A) = "
          f"{np.trapz(L_tot, wave):.2e} erg/s "
          f"(ionizing budget ~ {L_abs:.1e} erg/s)")

    # ---- figure: ionization-front map ----
    os.makedirs(PLOTDIR, exist_ok=True)
    from transient.sim_spectrum import add_colorbar
    fig, ax = plt.subplots(figsize=(11, 6))
    pc = ax.pcolormesh(np.degrees(ph), np.degrees(th), r_if.T,
                       cmap='viridis', shading='auto')
    add_colorbar(pc, ax, r'${\rm ionization~front}~r_{\rm IF}~[r_0]$')
    ax.set_xlabel(r'$\rm azimuth~\phi~[deg]$', fontsize=17)
    ax.set_ylabel(r'$\rm polar~\theta~[deg]$', fontsize=17)
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=13)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()
    out = os.path.join(PLOTDIR, 'sim_cloudy_ifront.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")

    # ---- figure: predicted spectrum ----
    fig, ax = plt.subplots(figsize=(12, 7))
    ax.plot(wave, wave * L_tot, '-', color='crimson', lw=2,
            drawstyle='steps-mid', label=r'$\rm total~(Cloudy~slabs)$')
    ax.plot(wave, wave * L_cont, '--', color='gray', lw=1.5,
            drawstyle='steps-mid', label=r'$\rm diffuse~continuum$')
    ax.axvline(3646, color='gray', ls=':', lw=1.2)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(1000, 9000)
    msel = (wave >= 1000) & (wave <= 9000)
    ymax = (wave * L_tot)[msel].max()
    ax.set_ylim(ymax / 1e3, ymax * 3)
    from matplotlib.ticker import ScalarFormatter, NullFormatter
    ax.set_xticks([1000, 2000, 3000, 5000, 9000])
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel(r'$\rm wavelength~[\AA]$', fontsize=18)
    ax.set_ylabel(r'$\lambda L_\lambda~[\rm erg~s^{-1}]$',
                  fontsize=17)
    ax.text(0.03, 0.95,
            rf'$\rm lamp\mbox{{-}}post~irradiation,~'
            rf'median~\log\phi={np.median(logphi):.1f},~i='
            rf'{cfg["INCL_DEG"]:.0f}^\circ$',
            transform=ax.transAxes, va='top', fontsize=14)
    ax.legend(fontsize=13, frameon=False)
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=14)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()
    out = os.path.join(PLOTDIR, 'sim_cloudy_spectrum.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")


if __name__ == '__main__':
    main()
