#!/usr/bin/env python3
"""
sim_cloudy_linemaps.py - face-on emission-line maps of the simulated disk
for selected lines (C IV 1549, Mg II 2798, H alpha 6563).

Each (theta, phi) ray's line luminosity comes from integrating its Cloudy
slab LINE spectrum over a window around the line. The ray's luminosity is
then distributed back along its ionized cells proportionally to the local
recombination weight n^2 r^2 dr and summed over theta, giving a face-on
(r, phi) map of where each line forms. (The depth structure *within* a
slab is not resolved by the lookup - the along-ray distribution is the
recombination weight, an approximation good for the map's morphology.)

One figure per line: sim_linemap_<line>.png

Run:  python3 transient/sim_cloudy_linemaps.py [dump.athdf]
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from transient.sim_spectrum import (  # noqa: E402
    load_sim, to_physical, CONFIG, PLOTDIR, LD_CM)
from transient.sim_cloudy_spectrum import (  # noqa: E402
    load_grid, compute_rays, Q_ION, ALPHA_B)

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

# line label -> (tex, window w1, w2 in A)
LINES = {
    'civ': (r'$\rm C\,IV~\lambda1549$', 1515., 1585.),
    'mgii': (r'$\rm Mg\,II~\lambda2798$', 2755., 2845.),
    'halpha': (r'$\rm H\alpha~\lambda6563$', 6510., 6620.),
}

EMS_NPZ = os.path.join(HERE, 'data', 'cloudy_arm_ems_arm_column2.npz')


def geval(curve, xi, nxi):
    """Evaluate cumulative-emission curves (…, nxi) at xi in [0, 1]."""
    fx = xi * (nxi - 1)
    ix = np.clip(fx.astype(int), 0, nxi - 2)
    t = fx - ix
    return (np.take_along_axis(curve, ix, axis=-1) * (1 - t)
            + np.take_along_axis(curve, ix + 1, axis=-1) * t)


def ems_xi_curves(sim, phys, rays, cfg):
    """The smooth ingredients of the Cloudy-emissivity line placement:
    the per-cell cumulative absorbed-photon fractions (xi_lo, xi_hi)
    along every ray from the Strommgren march, and each ray's matched
    cumulative-emission curves G(xi) per line. Kept separate so callers
    can interpolate THESE (smooth in theta) rather than the final spiky
    deposition profiles - interpolating two sharply front-peaked
    profiles gives two half-spikes instead of one moving spike, which
    renders as nested 'rose petal' rings."""
    d = np.load(EMS_NPZ)
    xi_g = d['xi']
    pu = np.unique(d['phi'])
    hu = np.unique(d['hden'])
    cu = np.unique(d['colden'])
    M = np.zeros((pu.size, hu.size, cu.size), dtype=int)
    M[np.searchsorted(pu, d['phi']), np.searchsorted(hu, d['hden']),
      np.searchsorted(cu, d['colden'])] = np.arange(d['phi'].size)

    r, th, ph = sim['r'], sim['th'], sim['ph']
    nph, nth = len(ph), len(th)
    sh = (nph, nth)
    ipx = np.argmin(np.abs(rays['logphi'].reshape(sh)[..., None]
                           - pu[None, None, :]), axis=-1)
    ihx = np.argmin(np.abs(rays['logn'].reshape(sh)[..., None]
                           - hu[None, None, :]), axis=-1)
    icx = np.argmin(np.abs(rays['logN'].reshape(sh)[..., None]
                           - cu[None, None, :]), axis=-1)
    midx = M[ipx, ihx, icx]                            # (nph, nth)

    r0_cm = cfg['R0_LD'] * LD_CM
    drs = np.diff(sim['rf']) * r0_cm
    rcm = r * r0_cm
    dS = phys['nH'] ** 2 * ALPHA_B * rcm[None, None, :] ** 2 \
        * drs[None, None, :]
    S = np.cumsum(dS, axis=2)
    ion = rays['ion']
    Send = np.where(ion, S, 0.0).max(axis=2)
    Send = np.where(Send > 0, Send, 1.0)
    xi_hi = np.clip(S / Send[..., None], 0.0, 1.0)
    xi_lo = np.clip((S - dS) / Send[..., None], 0.0, 1.0)
    curves = {key: d[f'G_{key}'][midx] for key in LINES}
    return xi_lo, xi_hi, curves, ion, xi_g.size, S, dS


def ems_deposition(sim, phys, rays, cfg):
    """Per-line along-ray deposition fractions from Cloudy's OWN
    depth-resolved emissivities ('save lines emissivity', extracted by
    extract_cloudy_ems.py). Each ray is matched to its nearest grid
    model; the model's cumulative emission fraction G(xi), tabulated
    against the slab's cumulative RECOMBINATION fraction xi, is
    evaluated at the ray's cumulative absorbed-photon fraction from the
    Strommgren march. H alpha then reduces exactly to the validated
    n^2 r^2 dr recombination weight, and C IV / Mg II shift relative to
    it as Cloudy dictates. Returns {line: (nph, nth, nr) normalized}."""
    xi_lo, xi_hi, curves, ion, nxi, _, _ = ems_xi_curves(sim, phys, rays,
                                                         cfg)
    out = {}
    for key in LINES:
        frac = np.maximum(geval(curves[key], xi_hi, nxi)
                          - geval(curves[key], xi_lo, nxi), 0.0) * ion
        s = frac.sum(axis=2)
        zz = s <= 0
        frac[zz, :] = 0.0
        frac[zz, 0] = 1.0
        out[key] = frac / frac.sum(axis=2)[..., None]
    return out


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else \
        os.path.join(HERE, 'data', 'sim', 'disk.out1.00012.athdf')
    cfg = CONFIG
    sim = load_sim(path)
    phys = to_physical(sim, cfg)
    grid = load_grid()
    rays = compute_rays(sim, phys, cfg, grid)
    wave = grid[0]

    r, th, ph = sim['r'], sim['th'], sim['ph']
    nph, nth, nr = len(ph), len(th), len(r)
    F_line = rays['F_line'].reshape(nph, nth, -1)   # nuF_nu per cm^2
    dA = rays['dA'].reshape(nph, nth)

    # ---- along-ray deposition, theta-supersampled ----
    # Each ray's luminosity is spread over its ionized cells with the
    # recombination weight n^2 r^2 dr. Summing the 64 DISCRETE theta rays
    # directly renders as terraced 'onion layers' (each step = one polar
    # ray's ionization-front radius dropping out of the sum). The front
    # surface r_IF(theta) is physically continuous, so interpolate it and
    # the ray luminosities onto a NSUP-times finer theta grid, truncating
    # the radial cell containing the front fractionally. Pure rendering:
    # the total luminosity per line is unchanged.
    # NSUP fine rays between theta neighbors: with Cloudy's sharply
    # front-peaked emissivity profiles, consecutive fine rays drop
    # ~1-cell-wide spikes at positions (front gap)/NSUP apart, so NSUP
    # must be large enough for the spikes to overlap at the radial cell
    # scale or the map shows a picket fence of 'rose petal' rings
    # (verified in debug_rose_rings.py: peak spacing = front gap / NSUP)
    NSUP = 64
    rf_edges = sim['rf']
    drc = np.diff(rf_edges)
    w_full = phys['nH'] ** 2 * r[None, None, :] ** 2 * drc[None, None, :]
    r_if = rays['r_if'].reshape(nph, nth)
    matter = rays['matter'].reshape(nph, nth)
    rif_eff = np.where(matter, rf_edges[-1], r_if)
    th_fine = np.linspace(th[0], th[-1], NSUP * nth)
    # linear-in-theta interpolation weights for the along-ray profiles
    # (nearest-neighbor profiles snap the sharply peaked Mg II / C IV
    # deposition to the 64 discrete front radii and render as nested
    # 'rose petal' rings; the physical projection of the continuous bowl
    # is a smooth radial band)
    tpos = np.interp(th_fine, th, np.arange(nth))
    j0 = np.clip(np.floor(tpos).astype(int), 0, nth - 2)
    tw = (tpos - j0)[:, None]

    # along-ray line placement from Cloudy's OWN depth-resolved
    # emissivities (the arm_column2 'save lines emissivity' output):
    # each ray's line luminosity is distributed along its ionized cells
    # following its matched slab's cumulative emission-fraction curve
    # versus column density, so C IV sits where Cloudy forms it and
    # Mg II hugs the internal ionization front.
    L_lines = {}
    for key, (tex, w1, w2) in LINES.items():
        m = (wave >= w1) & (wave <= w2)
        Flam = F_line[:, :, m] / wave[None, None, m]
        L_lines[key] = dA * np.trapz(Flam, wave[m], axis=2)  # (ph, th)

    # fine-theta deposition: interpolate the SMOOTH ingredients (the
    # cumulative absorption fractions xi and the G curves) and evaluate
    # each fine ray's profile from them. Interpolating the final
    # profiles instead makes every sharply front-peaked line print
    # double images at both parents' front radii ('rose petals').
    # Each fine ray is built as a SELF-CONSISTENT pseudo-ray: blend the
    # raw cumulative absorption S(r) (smooth, unclipped) between the two
    # parent rays, and normalize by the fine ray's OWN budget, S at its
    # own theta-interpolated front. Blending the parents' clipped xi
    # profiles instead inherits each parent's hard saturation at ITS
    # front, so the per-theta front radii re-print as sharp layer edges
    # no matter how fine the supersampling (verified: the largest map
    # jumps sat exactly at individual coarse front radii).
    _, _, curves, ionc, nxi, S, dS = ems_xi_curves(sim, phys, rays, cfg)
    maps = {key: np.zeros((nph, nr)) for key in LINES}
    for i in range(nph):
        rif_f = np.interp(th_fine, th, rif_eff[i])
        trunc = np.clip((rif_f[:, None] - rf_edges[None, :-1])
                        / drc[None, :], 0.0, 1.0)
        Sf = (1.0 - tw) * S[i, j0, :] + tw * S[i, j0 + 1, :]
        dSf = (1.0 - tw) * dS[i, j0, :] + tw * dS[i, j0 + 1, :]
        idxr = np.clip(np.searchsorted(r, rif_f) - 1, 0, nr - 2)
        trr = np.clip((rif_f - r[idxr]) / (r[idxr + 1] - r[idxr]), 0, 1)
        Send_f = (np.take_along_axis(Sf, idxr[:, None], 1)[:, 0]
                  * (1 - trr)
                  + np.take_along_axis(Sf, (idxr + 1)[:, None], 1)[:, 0]
                  * trr)
        Send_f = np.maximum(Send_f, 1e-300)[:, None]
        xh = np.clip(Sf / Send_f, 0.0, 1.0)
        xl = np.clip((Sf - dSf) / Send_f, 0.0, 1.0)
        for key in LINES:
            cf = ((1.0 - tw) * curves[key][i, j0, :]
                  + tw * curves[key][i, j0 + 1, :])      # (nfine, nxi)
            prof = np.maximum(geval(cf, xh, nxi)
                              - geval(cf, xl, nxi), 0.0) * trunc
            s = prof.sum(axis=1)
            zz = s <= 0
            prof[zz, :] = 0.0
            prof[zz, 0] = 1.0
            prof /= prof.sum(axis=1)[:, None]
            Lf = np.interp(th_fine, th, L_lines[key][i]) / NSUP
            maps[key][i] = Lf @ prof

    os.makedirs(PLOTDIR, exist_ok=True)
    R, P = np.meshgrid(r, ph)
    X, Y = R * np.cos(P), R * np.sin(P)
    for key, (tex, w1, w2) in LINES.items():
        lmap = maps[key]
        Ltot = lmap.sum()
        print(f"{key:7s} L = {Ltot:.2e} erg/s")
        vmax = lmap.max()
        fig, ax = plt.subplots(figsize=(9.5, 8))
        pc = ax.pcolormesh(X, Y, np.clip(lmap, vmax / 3e3, None),
                           norm=LogNorm(vmin=vmax / 3e3, vmax=vmax),
                           cmap='magma', shading='auto')
        ax.set_aspect('equal')
        from transient.sim_spectrum import add_colorbar, mark_star
        add_colorbar(pc, ax, r'$L_{\rm line}~\rm per~cell~[erg~s^{-1}]$',
                     fontsize=14)
        mark_star(ax, sim)
        ax.set_xlabel(r'$x~[r_0]$', fontsize=16)
        ax.set_ylabel(r'$y~[r_0]$', fontsize=16)
        _m, _e = f'{Ltot:.0e}'.split('e')
        ax.text(0.02, 0.02, tex + '\n'
                + rf'$L = {_m}\times10^{{{int(_e)}}}~\rm erg~s^{{-1}}$',
                transform=ax.transAxes, va='bottom', fontsize=15,
                color='black')
        ax.tick_params(which='major', direction='in', length=8, width=1.5,
                       top=True, right=True, labelsize=13)
        ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                       top=True, right=True)
        ax.minorticks_on()
        out = os.path.join(PLOTDIR, f'sim_linemap_{key}.png')
        plt.savefig(out, dpi=200, bbox_inches='tight')
        plt.close()
        print(f"Saved {out}")

    # ---------------- RGB composite (R = Halpha, G = MgII, B = CIV) -------
    # each channel independently log-stretched over 3 dex and rendered with
    # its Reds/Greens/Blues colormap (high flux = deep saturated color);
    # channels combined multiplicatively on a white background.
    ngrid = 900
    xg = np.linspace(-r.max(), r.max(), ngrid)
    Xg, Yg = np.meshgrid(xg, xg)
    Rg = np.hypot(Xg, Yg)
    Pg = np.mod(np.arctan2(Yg, Xg), 2 * np.pi)
    inside = (Rg >= r.min()) & (Rg <= r.max())
    # BILINEAR sampling in (phi, r), periodic in phi: nearest-cell
    # sampling renders the 256 azimuthal columns as fine radial spokes
    # once the chroma gain amplifies them
    fr = np.clip(np.interp(Rg, r, np.arange(nr)), 0, nr - 1 - 1e-9)
    i0 = fr.astype(int)
    i1 = np.minimum(i0 + 1, nr - 1)
    fs = fr - i0
    phe = np.concatenate([[ph[-1] - 2 * np.pi], ph, [ph[0] + 2 * np.pi]])
    fp = np.interp(Pg, phe, np.arange(-1, nph + 1))
    j0 = np.floor(fp).astype(int) % nph
    j1 = (j0 + 1) % nph
    ft = fp - np.floor(fp)

    def sample(arr):
        return ((1 - ft) * (1 - fs) * arr[j0, i0]
                + ft * (1 - fs) * arr[j1, i0]
                + (1 - ft) * fs * arr[j0, i1]
                + ft * fs * arr[j1, i1])
    DR = 2.5                                   # dynamic range [dex]
    # additive astro-style composite: R = Halpha, G = MgII, B = CIV.
    # The three line maps are correlated at 0.89-0.98 (they all follow the
    # irradiation pattern), so a simple per-channel max normalization puts
    # R ~ G ~ B everywhere and renders as muddy olive. Instead each channel
    # is HISTOGRAM-EQUALIZED over the lit map pixels (rank transform of
    # log L): the hue then encodes where each line is strong relative to
    # its own distribution, and every channel gets a fair share of the
    # color balance (Mg II arms green, C IV rim blue, H alpha outskirts
    # red). Qualitative morphology only, like any per-channel-normalized
    # composite.
    # per-channel stretch: SMOOTH percentile-anchored log stretch (the
    # rank/histogram equalization used before has a jumpy slope wherever
    # values cluster, and with the chroma gain that manufactured blocky
    # steps that do not exist in the single-line maps)
    rgb = np.zeros((ngrid, ngrid, 3))
    for k, key in enumerate(['halpha', 'mgii', 'civ']):
        lm = maps[key]
        logl = np.log10(lm + 1e-30)
        lit_m = logl > logl.max() - DR
        # coarse-knot equalization: anchor each channel at 13 of its own
        # percentiles (channel balance, so no single hue dominates) but
        # interpolate smoothly between anchors (a full rank transform is
        # jumpy and manufactured blocky steps)
        qs = np.linspace(0.0, 100.0, 13)
        anchors = np.percentile(logl[lit_m], qs)
        anchors = np.maximum.accumulate(anchors)
        anchors += 1e-9 * np.arange(anchors.size)      # strictly increasing
        eq = np.interp(logl, anchors, qs / 100.0)
        eq[~lit_m] = 0.0
        samp = sample(eq)
        samp[~inside] = 0.0
        rgb[..., k] = samp
    # chroma boost on top of the equalization: amplify each pixel's color
    # difference from gray so the line-ratio structure reads. With the
    # extended (unclipped) grid the three lines' equalized maps track
    # each other closely, so the REAL ratio differences are small and
    # need a strong gain to be visible - read hues as qualitative.
    SAT = 5.0
    lum = rgb.mean(axis=2, keepdims=True)
    rgb = np.clip(lum + SAT * (rgb - lum), 0, 1)
    # luminance gamma: darken the diffuse mid-level disk so the arms and
    # the inner rim carry the picture
    scale = np.where(lum > 0, lum ** 1.4 / (lum + 1e-30), 0.0)
    rgb = np.clip(rgb * scale, 0, 1)
    fig, ax = plt.subplots(figsize=(9.5, 9))
    ax.imshow(rgb, extent=[-r.max(), r.max(), -r.max(), r.max()],
              origin='lower', interpolation='bilinear')
    ax.set_aspect('equal')
    from transient.sim_spectrum import mark_star
    mark_star(ax, sim)
    ax.set_xlabel(r'$x~[r_0]$', fontsize=16)
    ax.set_ylabel(r'$y~[r_0]$', fontsize=16)
    for i, (txt, col) in enumerate([
            (r'$\rm H\alpha$', 'crimson'),
            (r'$\rm Mg\,II$', 'seagreen'),
            (r'$\rm C\,IV$', 'royalblue')]):
        ax.text(0.02, 0.97 - 0.055 * i, txt, transform=ax.transAxes,
                va='top', fontsize=17, color=col)
    ax.tick_params(which='major', direction='in', length=8, width=1.5,
                   top=True, right=True, labelsize=13)
    ax.tick_params(which='minor', direction='in', length=4, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()
    out = os.path.join(PLOTDIR, 'sim_linemap_rgb.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved {out}  (smooth percentile log stretch over the top "
          f"{DR:.1f} dex, chroma boost {SAT:.1f})")


if __name__ == '__main__':
    main()
