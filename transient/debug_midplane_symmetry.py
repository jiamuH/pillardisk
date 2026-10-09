"""Midplane (theta -> pi - theta) symmetry check of the hotspot simulation.

Reports the theta grid coverage and the mirror asymmetry of density and
temperature (cs^2 = P/rho), globally and in the hotspot region, and saves
a meridional slice through the star's azimuth to transient/plots/.

Run:  python3 -m transient.debug_midplane_symmetry
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from transient.sim_spectrum import load_sim, star_position  # noqa: E402

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'


def main():
    sim = load_sim(os.path.join(HERE, 'data', 'sim', 'disk.out1.00012.athdf'))
    r, th, ph = sim['r'], sim['th'], sim['ph']
    rho, press = sim['rho'], sim['press']          # shape (ph, th, r)
    cs2 = press / rho

    print(f"grid: n_phi={ph.size} n_theta={th.size} n_r={r.size}")
    print(f"theta range: {np.degrees(th.min()):.2f} .. {np.degrees(th.max()):.2f} deg "
          f"(midplane = 90 deg); r range {r.min():.2f} .. {r.max():.2f}")
    print(f"phi range: {np.degrees(ph.min()):.1f} .. {np.degrees(ph.max()):.1f} deg")

    # is the theta grid itself mirror symmetric?
    th_m = np.pi - th[::-1]
    print(f"theta grid mirror symmetric: {np.allclose(th, th_m, atol=1e-4)} "
          f"(max |theta - mirror| = {np.degrees(np.abs(th - th_m).max()):.2e} deg)")
    if not np.allclose(th, th_m, atol=1e-4):
        print("  -> only one hemisphere (or asymmetric grid); mirror comparison skipped")
        return

    # mirror arrays in theta
    rho_m = rho[:, ::-1, :]
    cs2_m = cs2[:, ::-1, :]
    up = th < np.pi / 2                             # upper hemisphere cells

    def stats(a, a_m, label, sel=None):
        d = np.abs(a - a_m) / (0.5 * (a + a_m))
        d = d[:, up, :] if sel is None else d[sel][:, up, :] if sel.ndim == 1 else d[sel]
        print(f"{label}: median |dA/A| = {np.median(d):.2e}, 90% = {np.percentile(d, 90):.2e}, "
              f"99% = {np.percentile(d, 99):.2e}, max = {d.max():.2e}")

    print("\n-- whole domain (upper hemisphere vs mirrored lower) --")
    stats(rho, rho_m, "rho")
    stats(cs2, cs2_m, "cs2 (T)")

    # hotspot region: within dr = 0.15, dphi = 0.3 rad of the star
    xs, ys = star_position(sim)
    rs, ps = np.hypot(xs, ys), np.arctan2(ys, xs) % (2 * np.pi)
    dphi = np.angle(np.exp(1j * (ph - ps)))
    selp = np.abs(dphi) < 0.3
    selr = np.abs(r - rs) < 0.15
    print(f"\nstar at r = {rs:.3f}, phi = {np.degrees(ps):.1f} deg; "
          f"hotspot box: {selp.sum()} phi cells x {selr.sum()} r cells")
    print("-- hotspot region --")
    for a, a_m, label in [(rho, rho_m, "rho"), (cs2, cs2_m, "cs2 (T)")]:
        d = np.abs(a - a_m) / (0.5 * (a + a_m))
        d = d[np.ix_(selp, up, selr)]
        print(f"{label}: median |dA/A| = {np.median(d):.2e}, 90% = {np.percentile(d, 90):.2e}, "
              f"max = {d.max():.2e}")

    # vertical (theta) velocity: should be antisymmetric about the midplane
    v2, v2_m = sim['v2'], -sim['v2'][:, ::-1, :]
    d = np.abs(v2 - v2_m)[:, up, :]
    print(f"\nv_theta antisymmetry: median |v2 + v2_mirror| = {np.median(d):.2e} "
          f"(median |v2| = {np.median(np.abs(v2)):.2e}) [code units]")

    # meridional slice through the star's azimuth
    ip = np.argmin(np.abs(dphi))
    R = r[None, :] * np.sin(th)[:, None]
    Z = r[None, :] * np.cos(th)[:, None]
    fig, ax = plt.subplots(figsize=(9, 7))
    pc = ax.pcolormesh(R, Z, np.log10(rho[ip]), shading='auto', cmap='viridis')
    ax.axhline(0, color='w', lw=1.5, ls='--')
    ax.plot(rs, 0, marker='*', ms=18, mfc='white', mec='k', ls='none')
    ax.set_xlabel(r'$R~[r_0]$')
    ax.set_ylabel(r'$z~[r_0]$')
    ax.set_aspect('equal')
    ax.minorticks_on()
    ax.tick_params(which='both', direction='in', top=True, right=True)
    ax.tick_params(which='major', length=8, width=1.5)
    ax.tick_params(which='minor', length=4, width=1)
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    cax = make_axes_locatable(ax).append_axes('right', size='4.5%', pad=0.12)
    cb = fig.colorbar(pc, cax=cax)
    cb.set_label(r'$\log \rho~\rm [code]$')
    cb.ax.minorticks_on()
    cb.ax.tick_params(which='both', direction='in')
    os.makedirs(os.path.join(HERE, 'plots'), exist_ok=True)
    out = os.path.join(HERE, 'plots', 'debug_midplane_symmetry.png')
    fig.savefig(out, dpi=200, bbox_inches='tight')
    print(f"\nsaved {out}")


if __name__ == '__main__':
    main()
