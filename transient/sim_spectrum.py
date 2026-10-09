#!/usr/bin/env python3
"""
sim_spectrum.py - spectral prediction from a real Athena++ disk simulation.

Pipeline (v1):
  1. Read the athdf dump (spherical polar; rho, press, vel1/2/3).
  2. Map code units -> physical units with explicit, configurable anchors:
       r = 1        -> R0_LD light days (embedded star's orbit)
       central mass -> MBH_MSUN (sets v0 = sqrt(GM/r0), the velocity unit)
       density      -> midplane phi-mean n_H at r = 1 equals NMID_CM3
       temperature  -> two modes:
         'rescaled' (default): keep the sim's RELATIVE P/rho structure but
             anchor the r = 1 midplane to TMID_K. The sim's H/R ~ 0.05 is
             numerically inflated (real AGN disk at 2 ld: H/R ~ 1e-3), so
             the literal EOS temperature is dynamically, not thermally, set.
         'eos': literal T = mu m_H (P/rho) v0^2 / k_B.
  3. Per-cell hydrogen emissivity in the ionized limit (x_ion configurable):
     free-bound recombination continuum (Balmer + Paschen edges with
     exp(-(E-E_edge)/kT) decline), Case B Balmer lines, free-free.
     This slot is pluggable - a Cloudy lookup grid replaces `cell_emissivity`
     later for full photoionization microphysics.
  4. Doppler-shift each cell's line emission by its line-of-sight velocity
     (observer at inclination INCL_DEG) and sum -> L_lambda with real
     velocity-resolved line profiles from the sim kinematics.

Outputs: transient/plots/sim_maps_nH.png, sim_maps_T.png,
         sim_spectrum.png, sim_line_profiles.png

Run:  python3 transient/sim_spectrum.py [dump.athdf]
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'data', 'sim'))
import athena_read  # noqa: E402

PLOTDIR = os.path.join(HERE, 'plots')

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

# ---------------- physical anchors (explicit, changeable) ----------------
CONFIG = dict(
    R0_LD=2.0,            # r = 1 in light days (star's orbital radius)
    MBH_MSUN=1.0e7,       # central black-hole mass (sets v0, t0)
    NMID_CM3=1.0e10,      # midplane phi-mean n_H at r = 1  [cm^-3]
    TMID_K=3000.0,        # midplane T anchor at r = 1 ('rescaled' mode) [K]
    TEMP_MODE='rescaled',  # 'rescaled' | 'eos'
    X_ION=1.0,            # ionized fraction (v1: uniform; Cloudy later)
    MU=0.6,               # mean molecular weight (ionized)
    INCL_DEG=45.0,        # observer inclination from the disk axis
)

# constants (cgs)
MSUN, G, KB, MH, C_KMS = 1.989e33, 6.674e-8, 1.381e-16, 1.673e-24, 2.998e5
LD_CM = 2.59e15
HK = 1.4388e8  # hc/k in Angstrom*K

CASE_B = {r'H\alpha': (6563., 2.86), r'H\beta': (4861., 1.0),
          r'H\gamma': (4340., 0.468), r'H\delta': (4102., 0.259)}
RATIO_BAC_HB = 8.0     # L(Balmer continuum)/L(Hbeta), Case B ~1e4 K
PASCHEN_FRAC = 0.35    # Paschen-continuum edge strength rel. to Balmer


def load_sim(path):
    ds = athena_read.athdf(path)
    f64 = lambda k: np.asarray(ds[k], dtype=np.float64)
    out = dict(r=f64('x1v'), th=f64('x2v'), ph=f64('x3v'),
               rf=f64('x1f'), thf=f64('x2f'), phf=f64('x3f'),
               rho=f64('rho'), press=f64('press'),
               v1=f64('vel1'), v2=f64('vel2'), v3=f64('vel3'),
               time=float(ds['Time']))
    return out


def to_physical(sim, cfg):
    """attach physical n_H, T, velocity [km/s] and cell volumes [cm^3]."""
    r0 = cfg['R0_LD'] * LD_CM
    v0 = np.sqrt(G * cfg['MBH_MSUN'] * MSUN / r0)      # cm/s
    t0 = r0 / v0
    r, th, ph = sim['r'], sim['th'], sim['ph']
    jm = np.argmin(np.abs(th - np.pi / 2))
    i1 = np.argmin(np.abs(r - 1.0))
    rho_ref = sim['rho'][:, jm, i1].mean()             # midplane, r=1
    nH = cfg['NMID_CM3'] * sim['rho'] / rho_ref
    cs2 = sim['press'] / sim['rho']                    # code units
    if cfg['TEMP_MODE'] == 'eos':
        T = cfg['MU'] * MH * cs2 * v0 ** 2 / KB
    else:
        cs2_ref = np.median(cs2[:, jm, i1])
        T = cfg['TMID_K'] * cs2 / cs2_ref
    # velocities in km/s
    vr = sim['v1'] * v0 / 1e5
    vt = sim['v2'] * v0 / 1e5
    vp = sim['v3'] * v0 / 1e5
    # cell volumes r^2 sin(th) dr dth dph  (cm^3); broadcast (ph,th,r)
    dr = np.diff(sim['rf']) * r0
    dth = np.diff(sim['thf'])
    dph = np.diff(sim['phf'])
    vol = (dph[:, None, None] * (np.sin(th) * dth)[None, :, None]
           * ((r * r0) ** 2 * dr)[None, None, :])
    print(f"units: v0 = {v0/1e5:.0f} km/s, t0 = {t0/86400:.1f} d, "
          f"P_orb(r=1) = {2*np.pi*t0/86400:.0f} d, "
          f"sim time = {sim['time']*t0/86400:.0f} d")
    print(f"n_H:  {nH.min():.2e} - {nH.max():.2e} cm^-3")
    print(f"T ({cfg['TEMP_MODE']}): {T.min():.0f} - {T.max():.0f} K "
          f"(midplane r=1: {np.median(T[:, jm, i1]):.0f} K)")
    return dict(nH=nH, T=T, vr=vr, vt=vt, vp=vp, vol=vol)


def vlos_cube(sim, incl_deg, phi_obs_deg=0.0):
    """line-of-sight velocity [km/s] toward the observer at inclination i
    and azimuth phi_obs (observer direction
    o = (sin i cos p, sin i sin p, cos i); positive = approaching)."""
    i = np.radians(incl_deg)
    po = np.radians(phi_obs_deg)
    si, ci = np.sin(i), np.cos(i)
    th, ph = sim['th'], sim['ph']
    st, ct = np.sin(th)[None, :, None], np.cos(th)[None, :, None]
    dph = (ph - po)[:, None, None]
    # unit vectors dotted with observer direction
    r_o = si * st * np.cos(dph) + ci * ct
    t_o = si * ct * np.cos(dph) - ci * st
    p_o = -si * np.sin(dph)
    return None, (r_o, t_o, p_o)


def predict_spectrum(sim, phys, cfg, wave):
    """L_lambda [erg/s/A] in the ionized limit, with Doppler line profiles."""
    nH, T, vol = phys['nH'], phys['T'], phys['vol']
    ne_np_V = (cfg['X_ION'] * nH) ** 2 * vol            # n_e n_p dV per cell
    # Case B Hbeta effective recombination (T-scaled)
    alpha_hb = 3.03e-14 * (np.clip(T, 500, 3e4) / 1e4) ** -0.87
    hnu_hb = 6.626e-27 * 2.998e10 / (4861e-8)
    L_hb_cell = ne_np_V * alpha_hb * hnu_hb             # erg/s per cell
    L_hb = L_hb_cell.sum()
    print(f"integrated L(Hbeta) = {L_hb:.2e} erg/s "
          f"(ionized limit, x={cfg['X_ION']})")

    _, (r_o, t_o, p_o) = vlos_cube(sim, cfg['INCL_DEG'])
    vlos = phys['vr'] * r_o + phys['vt'] * t_o + phys['vp'] * p_o  # km/s

    L_lam = np.zeros_like(wave)
    dw = np.gradient(wave)
    # ---- lines: histogram cell luminosities in Doppler-shifted wavelength
    w_l = L_hb_cell.ravel()
    v_l = vlos.ravel()
    edges = np.concatenate([wave - dw / 2, [wave[-1] + dw[-1] / 2]])
    for name, (lam0, ratio) in CASE_B.items():
        lam_obs = lam0 * (1.0 - v_l / C_KMS)
        h, _ = np.histogram(lam_obs, bins=edges, weights=w_l * ratio)
        L_lam += h / dw
    # ---- free-bound continuum: bin cells by T, one edge-shape per bin
    logT = np.log10(np.clip(T, 500, None)).ravel()
    tbins = np.linspace(logT.min(), logT.max(), 25)
    ib = np.clip(np.digitize(logT, tbins) - 1, 0, len(tbins) - 2)
    for k in range(len(tbins) - 1):
        sel = ib == k
        if not sel.any():
            continue
        Wk = (w_l[sel]).sum()                # Hbeta luminosity in this T bin
        Tk = 10 ** (0.5 * (tbins[k] + tbins[k + 1]))
        shape = np.zeros_like(wave)
        for ledge, s in ((3646.0, 1.0), (8204.0, PASCHEN_FRAC)):
            m = wave <= ledge
            dE_T = HK * (1.0 / wave[m] - 1.0 / ledge)   # (E-E_edge)/k [K]
            shape[m] += s * np.exp(-dE_T / Tk) / wave[m] ** 2
        norm = np.trapz(shape, wave)
        if norm > 0:
            L_lam += RATIO_BAC_HB * Wk * shape / norm
    # ---- free-free (smooth, Rayleigh-Jeans-to-Wien in the optical)
    # 4pi j_nu = 6.8e-38 T^-1/2 exp(-hnu/kT) n_e n_p ; convert to lambda
    for k in range(len(tbins) - 1):
        sel = ib == k
        if not sel.any():
            continue
        Tk = 10 ** (0.5 * (tbins[k] + tbins[k + 1]))
        env = (ne_np_V.ravel()[sel]).sum()
        jnu = 6.8e-38 / np.sqrt(Tk) * np.exp(-HK / (wave * Tk))
        L_lam += env * jnu * 2.998e18 / wave ** 2       # nu = c/lam
    return L_lam, vlos, L_hb_cell


def star_position(sim):
    """(x, y) of the embedded star: the compact HOT envelope, located as
    the midplane cs^2 = P/rho maximum (~10x the ring background) - the
    compression-heated gas in the star's potential well. Confirmed by the
    v_r sign flip (flow convergence) across the same point in both r and
    phi. NOTE: the global density maximum is NOT a reliable marker - the
    cold wake/arm crest can be denser than the envelope."""
    r, th, ph = sim['r'], sim['th'], sim['ph']
    jm = np.argmin(np.abs(th - np.pi / 2))
    cs2 = (sim['press'] / sim['rho'])[:, jm, :]
    mask = r > 0.6
    ip, ir = np.unravel_index(np.argmax(cs2[:, mask]), cs2[:, mask].shape)
    rs, ps = r[mask][ir], ph[ip]
    return rs * np.cos(ps), rs * np.sin(ps)


def mark_star(ax, sim):
    xs, ys = star_position(sim)
    ax.plot(xs, ys, marker='*', ms=22, mfc='white', mec='black', mew=1.5,
            ls='none', zorder=5)


def add_colorbar(pc, ax, label, fontsize=15):
    """Colorbar with exactly the height of the parent axes; major and
    minor ticks facing inward."""
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    cax = make_axes_locatable(ax).append_axes('right', size='4.5%',
                                              pad=0.12)
    cb = ax.figure.colorbar(pc, cax=cax)
    cb.set_label(label, fontsize=fontsize)
    cb.ax.minorticks_on()
    cb.ax.tick_params(which='major', direction='in', length=8, width=1.5,
                      labelsize=13)
    cb.ax.tick_params(which='minor', direction='in', length=4, width=1.0)
    return cb


def polar_map(ax, sim, field, jmid, label, log=True):
    r, ph = sim['r'], sim['ph']
    R, P = np.meshgrid(r, ph)
    X, Y = R * np.cos(P), R * np.sin(P)
    v = field[:, jmid, :]
    from matplotlib.colors import LogNorm
    pc = ax.pcolormesh(X, Y, v, norm=LogNorm() if log else None,
                       cmap='inferno', shading='auto')
    ax.set_aspect('equal')
    ax.set_xlabel(r'$x~[r_0]$', fontsize=16)
    ax.set_ylabel(r'$y~[r_0]$', fontsize=16)
    add_colorbar(pc, ax, label)


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else \
        os.path.join(HERE, 'data', 'sim', 'disk.out1.00012.athdf')
    cfg = CONFIG
    sim = load_sim(path)
    phys = to_physical(sim, cfg)
    jm = np.argmin(np.abs(sim['th'] - np.pi / 2))
    os.makedirs(PLOTDIR, exist_ok=True)

    # ---- midplane maps ----
    for field, lab, fname in [
            (phys['nH'], r'$n_{\rm H}~[\rm cm^{-3}]$', 'sim_maps_nH.png'),
            (phys['T'], r'$T~[\rm K]$', 'sim_maps_T.png')]:
        fig, ax = plt.subplots(figsize=(9, 8))
        polar_map(ax, sim, field, jm, lab)
        mark_star(ax, sim)
        ax.tick_params(which='major', direction='in', length=8, width=1.5,
                       top=True, right=True, labelsize=13)
        ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                       top=True, right=True)
        ax.minorticks_on()
        plt.savefig(os.path.join(PLOTDIR, fname), dpi=200,
                    bbox_inches='tight')
        plt.close()
        print(f"Saved {os.path.join(PLOTDIR, fname)}")

    # ---- predicted spectrum ----
    wave = np.linspace(2200., 7200., 2500)
    L_lam, vlos, L_hb_cell = predict_spectrum(sim, phys, cfg, wave)
    fig, ax = plt.subplots(figsize=(12, 7))
    ax.plot(wave, L_lam, '-', color='crimson', lw=2,
            drawstyle='steps-mid')
    ax.axvline(3646, color='gray', ls=':', lw=1.2)
    ax.set_yscale('log')
    ax.set_xlabel(r'$\rm wavelength~[\AA]$', fontsize=18)
    ax.set_ylabel(r'$L_\lambda~[\rm erg~s^{-1}~\AA^{-1}]$', fontsize=17)
    ax.text(0.03, 0.95,
            rf'$i={cfg["INCL_DEG"]:.0f}^\circ,~T_{{\rm mid}}='
            rf'{cfg["TMID_K"]:.0f}~{{\rm K}},~n_{{\rm mid}}=10^{{'
            rf'{np.log10(cfg["NMID_CM3"]):.0f}}}~{{\rm cm^{{-3}}}}$',
            transform=ax.transAxes, va='top', fontsize=14)
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=14)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()
    out = os.path.join(PLOTDIR, 'sim_spectrum.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")

    # ---- velocity-resolved Hbeta profile ----
    vgrid = np.linspace(-25000, 25000, 400)
    h, e = np.histogram(vlos.ravel(), bins=vgrid,
                        weights=L_hb_cell.ravel())
    vc = 0.5 * (e[1:] + e[:-1])
    fig, ax = plt.subplots(figsize=(11, 6.5))
    ax.plot(vc, h / np.diff(e), '-', color='royalblue', lw=2.5)
    ax.set_xlabel(r'$\rm line\mbox{-}of\mbox{-}sight~velocity~[km~s^{-1}]$',
                  fontsize=17)
    ax.set_ylabel(r'$dL({\rm H\beta})/dv~[\rm erg~s^{-1}/(km~s^{-1})]$',
                  fontsize=15)
    ax.axvline(0, color='gray', lw=1)
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=14)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()
    out = os.path.join(PLOTDIR, 'sim_line_profiles.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}")


if __name__ == '__main__':
    main()
