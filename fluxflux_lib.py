"""Shared machinery for the pillar-disc flux-flux experiments.

Holds the pieces that `test_fluxflux_embedded.py` (single geometry, embedded
mass scan) and `test_fluxflux_sigmaphi.py` (pillar width scan) both need: the
disc surface integration, the embedded-star hot spot temperature map, the
diffuse shadow glow, the lamp sweep, the PyROA-style driving light curve, the
observed PG 0844+349 decomposition, and the line-of-sight temperature render.

Nothing here executes on import.
"""
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D

from pillardisk.pillar_disk import C, DAY, FNU_TO_MJY, H, K, ANGSTROM

LD_CM = C * DAY                                          # 1 light day in cm
SIGMA_I = 2.0 * np.pi**4 * K**4 / (15.0 * H**3 * C**2)   # Planck intensity const
SIGMA_SB = np.pi * SIGMA_I                               # Stefan-Boltzmann, cgs
L_EDD_PER_MSUN = 1.26e38                                 # erg s^-1 per solar mass
AU_CM = 1.496e13


def b_nu(lam_A, T):
    """Planck function B_nu in erg cm^-2 s^-1 Hz^-1 ster^-1."""
    nu = C / (lam_A * ANGSTROM)
    x = np.clip(H * nu / (K * np.maximum(T, 1.0)), 1e-10, 700.0)
    return 2.0 * H * nu**3 / C**2 / (np.exp(x) - 1.0)


class DiscGrid:
    """Grid, projected-area weights and unit conversion for one disc, matching
    the integration used by PillarDisk.compute_sed."""

    def __init__(self, dk):
        self.dk = dk
        self.r2, self.p2 = np.meshgrid(dk.r, dk.phi, indexing='ij')
        self.h2 = dk.get_height(self.r2, self.p2)
        dr = np.diff(dk.r)[:, None]
        dh = self.h2[1:, :] - self.h2[:-1, :]
        ds = np.sqrt(dr**2 + dh**2)
        px = -(dh / (ds + 1e-10)) * np.cos(dk.phi)[None, :]
        pz = dr / (ds + 1e-10)
        dot = dk.ex * px + dk.ez * pz
        dot = np.where(dot > 0, dot, 0.0)
        self.da = ds * dk.r[:-1, None] * dk.dphi          # surface area, ld^2
        self.wgeo = self.da * dot                          # projected area, ld^2
        self.units = LD_CM**2 / (dk.d * LD_CM)**2 * FNU_TO_MJY
        self.tv2 = np.interp(self.r2.ravel(), dk.r,
                             dk.tv_base).reshape(self.r2.shape)

    def sed(self, T2d, wl):
        """F_nu (mJy) for an explicit temperature map on the full grid."""
        Tl = T2d[:-1, :] * self.dk.fcol
        f4 = self.dk.fcol**4
        return np.array([np.sum(b_nu(l, Tl) / f4 * self.wgeo)
                         for l in np.atleast_1d(wl)]) * self.units


def pillar_gauss(pillar, r2, p2):
    """The pillar's own Gaussian profile g_p(r, phi), periodic in phi, exactly
    as used by PillarDisk.get_height to raise the bump."""
    dr = r2 - pillar['r']
    dphi = np.mod(p2 - pillar['phi'] + np.pi, 2 * np.pi) - np.pi
    sp = pillar['sigma_phi']
    gphi = (np.exp(-0.5 * (dphi / sp)**2)
            + np.exp(-0.5 * ((dphi - 2 * np.pi) / sp)**2)
            + np.exp(-0.5 * ((dphi + 2 * np.pi) / sp)**2))
    return np.exp(-0.5 * (dr / pillar['sigma_r'])**2) * gphi


def spot_t4(grid, m_embedded):
    """T_spot^4 map (K^4, full grid) for embedded stars radiating at their
    Eddington luminosity and thermalized by their own pillar bump.

    m_embedded is the mass of the SINGLE embedded object sitting under one
    pillar: one pillar is raised by one object, so this is an individual mass,
    not a summed population. Reproducing the observed PG 0844 constant component
    needs about 1e4 to 1e5 Msun here, i.e. an intermediate-mass black hole
    rather than a star (see debug_one_star_per_pillar.py). Also returns
    the per-pillar footprint areas (cm^2) and peak spot temperatures (K)."""
    dk = grid.dk
    t4 = np.zeros_like(grid.r2)
    areas, tpeak = [], []
    l_star = m_embedded * L_EDD_PER_MSUN                 # erg/s per pillar
    da_cm2 = grid.da * LD_CM**2
    for pillar in dk.pillars:
        g = pillar_gauss(pillar, grid.r2, grid.p2)
        a_p = np.sum(g[:-1, :] * da_cm2)                 # cm^2
        flux = l_star * g / a_p                          # erg s^-1 cm^-2
        t4[:-1, :] += flux[:-1, :] / SIGMA_SB
        areas.append(a_p)
        tpeak.append((l_star * g.max() / a_p / SIGMA_SB) ** 0.25)
    t4[-1, :] = t4[-2, :]                                # unused by the integral
    return t4, np.array(areas), np.array(tpeak)


def spot_t4_rhd(grid, m_embedded, renormalize=True):
    """T_spot^4 map for a star buried at the midplane under each pillar, with
    the heated patch size following from the luminosity instead of from the
    pillar's own Gaussian.

    The star sits at (r_p, phi_p, z = 0) and the disc surface above it lies at
    height h(r, phi). For a surface whose normal is close to vertical, the
    normal flux from the buried point source is

        F = L h / (4 pi R^3),    R = sqrt(d^2 + h^2),

    with d the in-plane separation. Integrated over the whole upper surface this
    is exactly L/2, the half of the star's luminosity that escapes upward; the
    other half emerges from the underside and is never seen. The softening scale
    is the burial depth h itself, so the patch is broad and cool for a deeply
    buried star and compact and hot for a shallow one, with no free width.

    The radial grid only marginally resolves the core (about 29% of the power
    falls inside d < h), so with renormalize=True each pillar's discrete profile
    is rescaled to carry exactly L/2. Returns the T_spot^4 map, the per-pillar
    fraction of L/2 that the raw discrete sum captured, and the peak spot
    temperatures (K)."""
    dk = grid.dk
    t4 = np.zeros_like(grid.r2)
    da_cm2 = grid.da * LD_CM**2
    h_cm = grid.h2 * LD_CM
    l_half = 0.5 * m_embedded * L_EDD_PER_MSUN            # upward half, erg/s
    caught, tpeak = [], []
    for pillar in dk.pillars:
        d2 = ((grid.r2**2 + pillar['r']**2
               - 2.0 * grid.r2 * pillar['r']
               * np.cos(grid.p2 - pillar['phi'])) * LD_CM**2)
        R = np.sqrt(d2 + h_cm**2)
        F = (2.0 * l_half) * h_cm / (4.0 * np.pi * R**3)   # erg s^-1 cm^-2
        got = np.sum(F[:-1, :] * da_cm2)
        caught.append(got / l_half)
        if renormalize:
            F = F * (l_half / got)
        t4[:-1, :] += F[:-1, :] / SIGMA_SB
        tpeak.append((F.max() / SIGMA_SB) ** 0.25)
    t4[-1, :] = t4[-2, :]                                 # unused by the integral
    return t4, np.array(caught), np.array(tpeak)


def heating_radius(m_embedded, t_amb):
    """Radius (cm) at which a buried star's flux equals the ambient disc flux,
    L/(4 pi R^2) = sigma T_amb^4. The natural size of the warmed patch."""
    lum = m_embedded * L_EDD_PER_MSUN
    return np.sqrt(lum / (4.0 * np.pi * SIGMA_SB * t_amb**4))


def diffuse_amplitude(grid, f_diff, dnu):
    """Flat-F_nu amplitude (mJy) of the constant diffuse glow filling the pillar
    shadows: f_diff of the local direct-irradiation power spread over dnu,
    evaluated at the nominal lamp state so it is the same in every state. Zero
    for the bowl (no pillars, no shadows)."""
    dk = grid.dk
    Tns = dk.get_temperature(grid.r2, grid.p2, compute_shadows=False)
    Tdir4 = np.maximum(Tns**4 - grid.tv2**4, 0.0)
    if len(dk.pillars) > 0:
        smask = dk._compute_shadow_mask(grid.r2, grid.p2, grid.h2)
    else:
        smask = np.ones_like(grid.r2)
    level = f_diff * SIGMA_I * Tdir4[:-1, :] / dnu
    return np.sum(level * (1.0 - smask[:-1, :]) * grid.wgeo) * grid.units


def diffuse_flux(d0, wl):
    """The diffuse glow at the requested wavelengths (mJy): flat in F_nu with a
    Balmer bound-free rise blueward of 3646 A."""
    return d0 * np.where(np.atleast_1d(wl) < 3646.0, 1.5, 1.0)


def lamp_states(grid, s_vals, tag, verbose=True):
    """T^4 (K^4) of the viscous + irradiated disc for every lamp state. The map
    does not depend on the embedded stars, so it is computed once and reused for
    every embedded-mass value."""
    dk = grid.dk
    tx0 = dk.tx_base.copy()
    out = np.empty((len(s_vals),) + grid.r2.shape)
    for j, s in enumerate(s_vals):
        dk.tx_base = tx0 * s
        out[j] = dk.get_temperature(grid.r2, grid.p2,
                                    compute_shadows=True) ** 4
        if verbose:
            print(f"  [{tag}] lamp state {j+1}/{len(s_vals)} (s={s:.3f})",
                  flush=True)
    dk.tx_base = tx0
    return out


def eta_lp(grid):
    """Bolometric lamp / viscous power ratio at the nominal (s=1) state."""
    T = grid.dk.get_temperature(grid.r2, grid.p2, compute_shadows=True)
    tx4 = np.maximum(T**4 - grid.tv2**4, 0.0)
    return (np.sum(tx4[:-1, :] * grid.da)
            / np.sum(grid.tv2[:-1, :]**4 * grid.da))


def driving_lc(F, camp):
    """PyROA-style driving light curve: standardize each band over the campaign
    states, average across all bands, renormalize to zero mean / unit variance."""
    B = F[camp].mean(axis=0)
    A = F[camp].std(axis=0)
    X = ((F - B) / A).mean(axis=1)
    return (X - X[camp].mean()) / X[camp].std()


def nonvar_anchor(fit):
    """X_non-var: the least negative band F = 0 crossing, i.e. the first band to
    reach zero flux as the driver decreases."""
    return np.max(-fit[:, 1] / fit[:, 0])


# --- observed PG 0844+349 flux-flux decomposition --------------------------
# W2-anchored SED from jax_roa/uvot_fluxflux.py: rest-frame effective
# wavelength and the bright, faint and recovered-host states in mJy, already
# extinction corrected.
PG0844_SED = '/Users/jiamuh/python/pg0844/plots/UVOT/uvot_sed_W2gal.dat'
PG_KEYS = [Line2D([], [], color='darkred', marker='s', mec='k', mew=0.8, ms=8,
                  ls='none', label=r'$\rm PG\,0844:~host$'),
           Line2D([], [], color='royalblue', marker='o', mfc='white', mew=1.8,
                  ms=7, ls='none', label=r'$\rm PG\,0844:~bright-host$'),
           Line2D([], [], color='navy', marker='o', mfc='white', mew=1.8,
                  ms=7, ls='none', label=r'$\rm PG\,0844:~faint-host$'),
           Line2D([], [], color='gray', marker='o', mfc='white', mew=1.8,
                  ms=7, ls='none', alpha=0.35,
                  label=r'$\rm PG\,0844:~raw~bright,~faint$')]


def load_pg0844():
    """Observed PG 0844+349 bright / faint / host SED (mJy) versus rest-frame
    wavelength, sorted by wavelength. None if the file is not on this machine."""
    if not os.path.exists(PG0844_SED):
        print(f"NOTE: {PG0844_SED} not found, skipping the PG 0844 overlay")
        return None
    d = np.loadtxt(PG0844_SED, usecols=(1, 2, 3, 4, 5, 6, 7))
    d = d[np.argsort(d[:, 0])]
    out = {'lam': d[:, 0], 'nu': C / (d[:, 0] * ANGSTROM)}
    for i, key in enumerate(['bright', 'faint', 'host']):
        out[key] = d[:, 1 + 2 * i]
        out[key + '_e'] = d[:, 2 + 2 * i]
    return out


def overlay_pg0844(ax, o, nul, host_only=False):
    """Overplot the observed PG 0844 decomposition in nu L_nu (nul converts mJy
    to erg/s once multiplied by nu). The raw bright and faint states are drawn
    faded; the host-subtracted states and the flux-flux host are opaque."""
    if o is None:
        return
    nl = o['nu'] * nul
    if not host_only:
        for key, col in [('bright', 'royalblue'), ('faint', 'navy')]:
            ax.errorbar(o['lam'], nl*o[key], yerr=nl*o[key + '_e'], fmt='o',
                        mfc='white', mec=col, ms=7, mew=1.8, ecolor=col,
                        alpha=0.3, ls='none', zorder=8)
            sub = o[key] - o['host']
            err = np.hypot(o[key + '_e'], o['host_e'])
            ax.errorbar(o['lam'], nl*sub, yerr=nl*err, fmt='o', mfc='white',
                        mec=col, ms=7, mew=1.8, ecolor=col, ls='none', zorder=9)
    ax.errorbar(o['lam'], nl*o['host'], yerr=nl*o['host_e'], fmt='s',
                color='darkred', mec='k', mew=0.8, ms=8, ecolor='darkred',
                ls='none', zorder=11)


def mass_label(m):
    """Plain a x 10^b solar-mass label, avoiding log-exponent shorthand."""
    e = int(np.floor(np.log10(m)))
    return rf'${m/10**e:.1f} \times 10^{{{e}}}~M_\odot$'


def ticks(ax):
    """Inward major and minor ticks on all four sides."""
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=14)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()


def los_map(grid, T2d, fname, notes=(), lim=20.0, vmin=5.0e2, vmax=4.0e4,
            spots=None, spot_T=None):
    """Depth-sorted sky-plane render of a surface temperature map, as in
    transient/plot_los_view.py: the surface is projected for the observer
    direction e = (sin i, 0, cos i) and drawn far to near so raised pillars
    cover what they hide."""
    dk = grid.dk
    x = grid.r2 * np.cos(grid.p2)
    y = grid.r2 * np.sin(grid.p2)
    z = grid.h2
    nx, ny, nz = dk._surface_normal(grid.r2, grid.p2, z)
    dot = nx * dk.ex + ny * dk.ey + nz * dk.ez
    x_sky = y
    y_sky = -x * dk.cosi + z * dk.sini
    depth = x * dk.sini + z * dk.cosi          # larger = nearer the observer

    vis = dot > 0.0
    order = np.argsort(depth[vis].ravel())      # far first, near drawn last
    xs = x_sky[vis].ravel()[order]
    ys = y_sky[vis].ravel()[order]
    Ts = T2d[vis].ravel()[order]

    # marker area matched to the local cell size so the surface tiles with no
    # gaps: cells are the larger of the radial and azimuthal spacing
    cell = np.maximum(np.gradient(dk.r)[:, None] * np.ones_like(grid.p2),
                      grid.r2 * dk.dphi)[vis].ravel()[order]
    ax_w_in = 0.72 * 11.0                       # axes width after the colorbar
    pt_per_ld = ax_w_in * 72.0 / (2.0 * lim)
    sizes = np.clip((1.4 * cell * pt_per_ld)**2, 2.0, 400.0)

    norm = LogNorm(vmin=vmin, vmax=vmax)
    fig, ax = plt.subplots(figsize=(11, 3.0 + 8.0 * dk.cosi))
    ax.set_facecolor('black')
    ax.scatter(xs, ys, c=cm.inferno(norm(Ts)), s=sizes, marker='s',
               linewidths=0, rasterized=True)
    if spots is not None and spot_T is not None:
        # unresolved thermalizing spheres: far smaller than a grid cell, so they
        # are drawn as points at their own temperature rather than painted into
        # the surface map
        sr = np.array([p[0] for p in spots])
        sp_ = np.array([p[1] for p in spots])
        sh = np.interp(sr, dk.r, dk.h_base)
        sx, sy = sr * np.cos(sp_), sr * np.sin(sp_)
        vis_s = np.cos(sp_) * dk.sini >= -1.0        # all front-facing at i<90
        ax.scatter((sy)[vis_s], (-sx * dk.cosi + sh * dk.sini)[vis_s],
                   c=[cm.inferno(norm(spot_T))], s=14, marker='o',
                   edgecolors='white', linewidths=0.4, zorder=5)
    sm = cm.ScalarMappable(norm=norm, cmap='inferno')
    cb = plt.colorbar(sm, ax=ax, shrink=0.8, aspect=20, extend='neither')
    cb.set_label(r'$\rm surface~temperature~[K]$', fontsize=16)
    cb.ax.tick_params(direction='in', labelsize=13)

    inc_deg = np.degrees(np.arccos(dk.cosi))
    ax.text(0.02, 0.97, rf'$i = {inc_deg:.0f}^\circ$', transform=ax.transAxes,
            va='top', fontsize=20, color='white')
    for j, note in enumerate(notes):
        ax.text(0.02, 0.88 - 0.07 * j, note, transform=ax.transAxes, va='top',
                fontsize=13, color='white')
    ax.text(0.98, 0.03, r'$\rm near~side$', transform=ax.transAxes, ha='right',
            fontsize=13, color='lightgray')
    ax.text(0.98, 0.97, r'$\rm far~side$', transform=ax.transAxes, ha='right',
            va='top', fontsize=13, color='lightgray')

    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim * dk.cosi * 1.15, lim * dk.cosi * 1.15)
    ax.set_aspect('equal')
    ax.set_xlabel(r'$\rm sky~X~[light~days]$', fontsize=18)
    ax.set_ylabel(r'$\rm sky~Y~[light~days]$', fontsize=18)
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=14, color='gray')
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True, color='gray')
    ax.minorticks_on()
    plt.savefig(fname, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {fname}")
