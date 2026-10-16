#!/usr/bin/env python3
"""
sim_full_spectrum.py - the complete predicted spectrum with all emission
components, decomposed:

  1. reprocessed  - Cloudy slab emission from AGN-photoionized gas
                    (face selection, Lambert projection, per-slab Doppler);
  2. TDE envelope - ENERGY-BUDGET effective temperature: the sim supplies
                    only the geometry (R_env from the hot-knot extent); the
                    radiated power is the physical input L_TDE, so
                    sigma T_eff^4 = L_TDE / (4 pi R_env^2) and the bump
                    luminosity equals L_TDE by construction. (The sim's
                    compression-heated gas temperature is a dynamical
                    closure feature, NOT a radiative flux - never use it
                    as T_eff.)
  3. disk         - energy-budget effective temperature per element:
                    sigma T_eff^4 = sigma T_visc^4 + F_irr, with T_visc =
                    the anchored (rescaled) midplane temperature and F_irr
                    the lamp's NON-IONIZING heating flux L_LAMP_HEAT,
                    deposited at each ray's ionization-front footprint
                    with the foreshortening factor (the ionizing part of
                    the lamp is already reradiated by component 1 - no
                    double counting). Lambert 4 cos(i) projection.
  4. AGN direct   - the lamp itself, with the SAME incident SED the
                    Cloudy decks use (their `interpolate` table),
                    normalized to the ionizing photon rate Q_ION.

Doppler shifts of the smooth thermal continua (< ~1%) are neglected.

Output: sim_full_spectrum.png (+ printed component luminosities)

Run:  python3 transient/sim_full_spectrum.py [dump.athdf]
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, NullFormatter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from transient.sim_spectrum import (  # noqa: E402
    load_sim, to_physical, star_position, CONFIG, LD_CM, C_KMS, PLOTDIR)
from transient.sim_cloudy_spectrum import (  # noqa: E402
    load_grid, compute_rays, Q_ION)

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

H_CGS, C_CGS, KB = 6.626e-27, 2.998e10, 1.381e-16
SB = 5.67e-5
NU0 = 3.29e15  # Lyman-limit frequency [Hz]

L_TDE = 1.0e42        # erg/s: radiated power of the TDE envelope (input;
                      # Eddington-embedded-star scale from the earlier work)
L_LAMP_HEAT = 4.0e42  # erg/s: the lamp's NON-ionizing luminosity that
                      # thermally heats the irradiated surfaces
                      # (alpha_nu = -1: ~2x the ionizing power)
THERMAL_MODE = 'sim_teff'  # 'sim_teff': ONE thermal component - the sim's
                      # temperature sampled at the tau~1 photosphere of each
                      # column (disk + TDE + arms together; assumes the
                      # closure's temperature contrast = emitted-flux
                      # contrast, i.e. the implied L_TDE is hidden in the
                      # anchor). 'energy_budget': explicit split -
                      # envelope from L_TDE, disk from T_visc + F_irr.
N_PHOT = 1.0e24       # cm^-2: column to the thermal photosphere


# incident AGN SED of the Cloudy grids (the decks' `interpolate` table:
# log E[Ryd] vs log f_nu, arbitrary normalization) - the direct lamp light
# MUST have the same shape the photoionized gas sees
SED_X = np.array([-8.0, -1.7, -1.3, -1.0, -0.19, 0.4, 1.87, 2.43,
                  2.98, 3.3, 3.7, 4.18, 4.99, 5.1, 7.0])
SED_Y = np.array([-13.6, -1.0, 1.48, 2.0, 1.65, 1.32, -1.87, -2.37,
                  -2.79, -2.96, -3.38, -4.05, -5.59, -6.67, -12.37])


def agn_direct(wave):
    """L_lambda of the lamp with the SAME SED Cloudy sees, normalized so
    the ionizing (E > 1 Ryd) photon rate equals Q_ION."""
    xg = np.linspace(0.0, 6.0, 2000)          # log E[Ryd], ionizing range
    integ = np.trapz(10.0 ** np.interp(xg, SED_X, SED_Y), xg) * np.log(10.0)
    A = Q_ION * H_CGS / integ                 # f_nu normalization
    x = np.log10(911.76 / wave)               # log E[Ryd] per wavelength
    fnu = A * 10.0 ** np.interp(x, SED_X, SED_Y)
    return fnu * C_CGS / (wave * 1e-8) ** 2 * 1e-8    # L_nu -> L_lambda/A


def b_lambda(wave_A, T):
    """Planck B_lambda [erg/s/cm^2/A/sr]; wave in Angstrom."""
    lam = wave_A * 1e-8
    x = np.clip(H_CGS * C_CGS / (lam * KB * T), 1e-6, 700.0)
    return 2 * H_CGS * C_CGS ** 2 / lam ** 5 / (np.exp(x) - 1.0) * 1e-8


def photosphere_T(sim, phys, cfg):
    """(phi, r) map of the sim temperature at the photosphere: the first
    cell from the top where the vertical column reaches N_PHOT (the
    midplane T where it never does)."""
    r, th = sim['r'], sim['th']
    r0_cm = cfg['R0_LD'] * LD_CM
    jm = np.argmin(np.abs(th - np.pi / 2))
    up = th <= np.pi / 2
    ds = (r[None, None, :] * r0_cm) * np.diff(sim['thf'])[None, :, None]
    Ncum = np.cumsum(phys['nH'][:, up, :] * ds[:, up, :], axis=1)
    above = Ncum >= N_PHOT
    j = np.argmax(above, axis=1)
    T_phot = np.take_along_axis(phys['T'][:, up, :], j[:, None, :],
                                axis=1)[:, 0, :]
    return np.where(above.any(axis=1), T_phot, phys['T'][:, jm, :])


def temperature_maps(sim, phys, cfg, path):
    """Face-on (phi, r) maps: photosphere T (the green thermal curve),
    midplane T (both the rescaled sim T) and the hotspot T_eff
    (hotspot_thermal), on one colour scale. The hotspot T_eff radiates
    only the heat in EXCESS of the phi-median at each (r, theta): it adds
    to the disk (T^4 sum); grey = no excess heat."""
    from matplotlib.colors import LogNorm
    from transient.sim_spectrum import add_colorbar, mark_star
    jm = np.argmin(np.abs(sim['th'] - np.pi / 2))
    T_ph = photosphere_T(sim, phys, cfg)
    T_hot = hotspot_thermal(sim, phys, cfg, np.array([5000.0, 5001.0]))[1]
    maps = [(T_ph, r'$\rm photosphere~T~(N_H=10^{24},~green~curve)$'),
            (phys['T'][:, jm, :], r'$\rm midplane~T$'),
            (T_hot, r'$\rm hotspot~T_{\rm eff}~(excess~heat/t_{\rm diff})$')]
    R, P = np.meshgrid(sim['rf'], sim['phf'])
    X, Y = R * np.cos(P), R * np.sin(P)
    norm = LogNorm(1e3, max(m.max() for m, _ in maps))
    fig, axs = plt.subplots(1, 3, figsize=(21, 6.6))
    for n, (ax, (m, title)) in enumerate(zip(axs, maps)):
        cmap = plt.get_cmap('inferno').copy()
        cmap.set_bad('0.8')
        val = np.where(m > 0, np.clip(m, 1e3, None), np.nan)
        pc = ax.pcolormesh(X, Y, val, norm=norm, cmap=cmap, shading='flat',
                           rasterized=True)
        ax.set_aspect('equal')
        ax.set_title(title, fontsize=17)
        ax.set_xlabel(r'$x~[r_0]$', fontsize=16)
        if n == 0:
            ax.set_ylabel(r'$y~[r_0]$', fontsize=16)
        mark_star(ax, sim)
        ax.tick_params(which='both', direction='in', top=True, right=True,
                       labelsize=13)
        ax.minorticks_on()
        add_colorbar(pc, ax, r'$T~[\rm K]$')
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"wrote {path}")


def thermal_and_agn(sim, phys, cfg, wave):
    """The ONE thermal component in 'sim_teff' mode - the sim temperature
    sampled at the tau~1 photosphere (vertical column N_H >= N_PHOT) of
    each (phi, r) column, disk + TDE + arms together - and the direct-AGN
    power law. Returns (L_th_base, L_agn); multiply L_th_base by
    4 cos(i) for the Lambert-projected observed value."""
    r = sim['r']
    r0_cm = cfg['R0_LD'] * LD_CM
    dA_disk = (r * r0_cm * np.diff(sim['rf']) * r0_cm)[None, :] * \
        np.diff(sim['phf'])[:, None]
    T_phot = photosphere_T(sim, phys, cfg)
    L_th_base = np.zeros_like(wave)
    tb = np.linspace(np.log10(T_phot.min()), np.log10(T_phot.max()) + 1e-6,
                     30)
    ibn = np.clip(np.digitize(np.log10(T_phot).ravel(), tb) - 1, 0, 28)
    for k in range(29):
        sel = ibn == k
        if not sel.any():
            continue
        Ak = dA_disk.ravel()[sel].sum()
        Tk = 10 ** (0.5 * (tb[k] + tb[k + 1]))
        L_th_base += Ak * np.pi * b_lambda(wave, Tk)
    L_agn = agn_direct(wave)
    print(f"thermal (sim T_eff photosphere, N={N_PHOT:.0e}): median "
          f"{np.median(T_phot):.0f} K, max {T_phot.max():.0f} K, "
          f"L_bol(2 faces) = {2*np.sum(dA_disk*SB*T_phot**4):.2e} erg/s")
    return L_th_base, L_agn


def hotspot_thermal(sim, phys, cfg, wave, gamma=5.0 / 3.0, hfrac=0.2):
    """Thermal emission of the TDE hotspot from its ENERGY, not the sim
    temperature (the sim is adiabatic, without cooling, so P/rho keeps
    every bit of dissipated heat). Per (phi, r) column: excess thermal
    energy E = sum (P - P_bg)/(gamma - 1) dV, P = rho (P/rho)_code v0^2
    (literal, set by the orbital speed), P_bg = phi-median at the same
    (r, theta); radiated over the photon diffusion time t = max(tau_es, 1)
    h / c (tau_es of the half column, h = hfrac r): L = E / t, two faces,
    T_eff = (L / 2 sigma A)^(1/4). E and tau both scale with n_mid, so L
    and T_eff do not depend on the density or temperature anchors. An
    upper bound: E is the heat stored over the whole run, which a cooling
    hotspot would have radiated. Returns (L_base as thermal_and_agn's,
    T_eff (phi, r), L_col)."""
    r, th = sim['r'], sim['th']
    r0 = cfg['R0_LD'] * LD_CM
    v0 = np.sqrt(6.674e-8 * cfg['MBH_MSUN'] * 1.989e33 / r0)
    rho = phys['nH'] * 1.4 * 1.6726e-24
    P = rho * sim['press'] / sim['rho'] * v0 ** 2
    dth = np.diff(sim['thf'])
    vol = (np.diff(sim['phf'])[:, None, None]
           * (np.sin(th) * dth)[None, :, None]
           * ((r * r0) ** 2 * np.diff(sim['rf']) * r0)[None, None, :])
    E = (np.clip(P - np.median(P, axis=0, keepdims=True), 0, None)
         / (gamma - 1) * vol).sum(axis=1)
    tau = 0.34 * (rho * (r * r0)[None, None, :]
                  * dth[None, :, None]).sum(axis=1) / 2
    L_col = E / (np.maximum(tau, 1.0) * hfrac * r * r0 / C_CGS)
    dA = (r * r0 * np.diff(sim['rf']) * r0)[None, :] * \
        np.diff(sim['phf'])[:, None]
    T = (L_col / (2 * SB * dA)) ** 0.25
    L_base = np.zeros_like(wave)
    hot = T > 1000.0
    tb = np.linspace(3.0, np.log10(T.max()) + 1e-6, 41)    # as thermal_and_agn
    ib = np.digitize(np.log10(T[hot]), tb) - 1
    for k in range(40):
        sel = ib == k
        if sel.any():
            L_base += dA[hot][sel].sum() * np.pi * b_lambda(
                wave, 10 ** (0.5 * (tb[k] + tb[k + 1])))
    Lt = L_col.sum()
    print(f"hotspot (excess heat / diffusion time, gamma = {gamma:.3g}): "
          f"L = {Lt:.2e} erg/s, T_eff L-weighted {np.sum(T * L_col) / Lt:.0f}"
          f" K, max {T.max():.0f} K")
    return L_base, T, L_col


def main():
    args = [x for x in sys.argv[1:] if x != '--tmaps']
    path = args[0] if args else \
        os.path.join(HERE, 'data', 'sim', 'disk.out1.00012.athdf')
    cfg = CONFIG
    sim = load_sim(path)
    phys = to_physical(sim, cfg)
    if '--tmaps' in sys.argv:          # temperature maps only
        temperature_maps(sim, phys, cfg,
                         os.path.join(PLOTDIR, 'sim_temperature_faceon.png'))
        return
    grid = load_grid()
    rays = compute_rays(sim, phys, cfg, grid)
    wave = grid[0]
    r0_cm = cfg['R0_LD'] * LD_CM
    i_obs = np.radians(cfg['INCL_DEG'])

    # ---- 1. reprocessed (Cloudy slabs), Doppler-assembled ----
    F_tot = rays['F_cont'] + rays['F_line']
    dA, v_ray = rays['dA'], rays['v_ray']
    L_rep = np.zeros_like(wave)
    shift = 1.0 - v_ray / C_KMS
    for i in range(len(dA)):
        ws = wave * shift[i]
        L_rep += dA[i] * np.interp(wave, ws, F_tot[i] / ws, left=0, right=0)

    # ---- thermal emission (mode switch) ----
    r, th, ph = sim['r'], sim['th'], sim['ph']
    jm = np.argmin(np.abs(th - np.pi / 2))
    dr = np.diff(sim['rf']) * r0_cm
    dphv = np.diff(sim['phf'])
    dth = np.diff(sim['thf'])
    dA_disk = (r * r0_cm * dr)[None, :] * dphv[:, None]   # cm^2
    proj = 4.0 * np.cos(i_obs)                        # Lambert, normal = z
    xs, ys = star_position(sim)
    rs, ps = np.hypot(xs, ys), np.arctan2(ys, xs) % (2 * np.pi)

    def bb_sum(Tmap, dAmap):
        """Lambert-projected blackbody sum over a (ph, r) T_eff map."""
        L = np.zeros_like(wave)
        tb = np.linspace(np.log10(Tmap.min()), np.log10(Tmap.max()) + 1e-6,
                         30)
        ibn = np.clip(np.digitize(np.log10(Tmap).ravel(), tb) - 1, 0, 28)
        for k in range(29):
            sel = ibn == k
            if not sel.any():
                continue
            Ak = dAmap.ravel()[sel].sum()
            Tk = 10 ** (0.5 * (tb[k] + tb[k + 1]))
            L += proj * Ak * np.pi * b_lambda(wave, Tk)
        return L

    if THERMAL_MODE == 'sim_teff':
        # ONE thermal component from the sim (disk + TDE + arms together)
        L_th = proj * thermal_and_agn(sim, phys, cfg, wave)[0]
    else:
        # ---- energy budget: TDE envelope from L_TDE ----
        R3, TH3, PH3 = np.meshgrid(r, th, ph, indexing='ij')
        T = np.transpose(phys['T'], (2, 1, 0))
        dphi = np.angle(np.exp(1j * (PH3 - ps)))
        dist2 = (R3 - rs) ** 2 + (R3 * np.cos(TH3)) ** 2 + (rs * dphi) ** 2
        near = dist2 < 0.15 ** 2
        ring = (np.abs(R3 - rs) < 0.1)
        T_bg = np.median(T[ring & ~near])
        env = near & (T > 4.0 * T_bg)
        R_env = 2.0 * np.sqrt(np.average(dist2[env],
                                         weights=T[env] ** 2)) * r0_cm
        T_env = (L_TDE / (4 * np.pi * R_env ** 2 * SB)) ** 0.25
        L_env = 4 * np.pi * R_env ** 2 * np.pi * b_lambda(wave, T_env)
        print(f"TDE envelope: R_env = {R_env/r0_cm:.3f} r0; "
              f"L_TDE = {L_TDE:.1e} erg/s -> T_eff = {T_env:.0f} K")
        # ---- disk: T_visc + lamp irradiation at the front footprints ----
        T_visc = phys['T'][:, jm, :]
        dOm = dphv[:, None] * (np.sin(th) * dth)[None, :]
        budget = Q_ION / (4 * np.pi)
        frac_abs = np.clip(rays['wsum'] / budget, 0.0, 1.0)
        P_ray = L_LAMP_HEAT / (4 * np.pi) * dOm * frac_abs
        R_fp = rays['r_if'] * np.sin(th)[None, :]
        ibin = np.clip(np.searchsorted(r, R_fp) - 1, 0, len(r) - 1)
        P_map = np.zeros_like(T_visc)
        for kph in range(len(ph)):
            np.add.at(P_map[kph], ibin[kph], P_ray[kph])
        Tmid = (T_visc ** 4 + (P_map / dA_disk) / SB) ** 0.25
        RM, PM = np.meshgrid(r, ph)
        dmid2 = (RM - rs) ** 2 \
            + (rs * np.angle(np.exp(1j * (PM - ps)))) ** 2
        dA_m = np.where(dmid2 < (R_env / r0_cm) ** 2, 0.0, dA_disk)
        L_disk = bb_sum(Tmid, dA_m)
        print(f"disk T_eff: visc median {np.median(T_visc):.0f} K, "
              f"with irradiation max {Tmid.max():.0f} K")

    # ---- 4. direct AGN: the Cloudy incident SED, normalized to Q_ION ----
    L_agn = agn_direct(wave)

    if THERMAL_MODE == 'sim_teff':
        L_total = L_rep + L_th + L_agn
    else:
        L_total = L_rep + L_env + L_disk + L_agn

    # ---- figure ----
    os.makedirs(PLOTDIR, exist_ok=True)
    fig, ax = plt.subplots(figsize=(12.5, 7.5))
    ax.plot(wave, wave * L_total, '-', color='black', lw=2.8, alpha=0.9,
            drawstyle='steps-mid', label=r'$\rm total$')
    ax.plot(wave, wave * L_rep, '-', color='crimson', lw=2, alpha=0.85,
            drawstyle='steps-mid', label=r'$\rm reprocessed~(Cloudy)$')
    if THERMAL_MODE == 'sim_teff':
        ax.plot(wave, wave * L_th, '-', color='seagreen', lw=2.4, alpha=0.9,
                label=r'${\rm thermal~(sim~}T_{\rm eff}{\rm ,~disk+TDE+arms)}$')
    else:
        ax.plot(wave, wave * L_env, '-', color='darkorange', lw=2.2,
                alpha=0.9, label=r'$\rm TDE~envelope~(thermal)$')
        ax.plot(wave, wave * L_disk, '-', color='seagreen', lw=2.2,
                alpha=0.9, label=r'$\rm disk~photosphere$')
    ax.plot(wave, wave * L_agn, '--', color='royalblue', lw=2.2, alpha=0.9,
            label=r'$\rm AGN~direct~(Cloudy~incident~SED)$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1000, 11000)
    msel = (wave >= 1000) & (wave <= 11000)
    ymax = (wave * L_total)[msel].max()
    ax.set_ylim(ymax / 3e3, ymax * 3)
    ax.set_xticks([1000, 2000, 3000, 5000, 10000])
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel(r'$\rm wavelength~[\AA]$', fontsize=18)
    ax.set_ylabel(r'$\lambda L_\lambda~[\rm erg~s^{-1}]$', fontsize=17)
    if THERMAL_MODE == 'sim_teff':
        note = (rf'$i={cfg["INCL_DEG"]:.0f}^\circ,~\rm sim~T_{{eff}}~at~'
                rf'N_H={N_PHOT:.0e}~cm^{{-2}}$'.replace('e+', r'\times10^{')
                .replace('~cm', r'}~cm'))
    else:
        note = (rf'$i={cfg["INCL_DEG"]:.0f}^\circ,~T_{{\rm env}}='
                rf'{T_env/1e3:.0f}~{{\rm kK}},~R_{{\rm env}}='
                rf'{R_env/r0_cm:.2f}~r_0$')
    ax.text(0.03, 0.96, note, transform=ax.transAxes, va='top', fontsize=14)
    ax.legend(fontsize=13, frameon=False, loc='upper right')
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=14)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()
    out = os.path.join(PLOTDIR, 'sim_full_spectrum.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")


if __name__ == '__main__':
    main()
