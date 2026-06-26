#!/usr/bin/env python3
"""
Generate 3D geometry figures for the paper.
Calls the existing plot_3d_geometry(color_by='flux') for each case,
then combines into a two-panel figure.

Usage:
    python3 plot_geometry_figures.py [config_line.yaml]
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import sys, os

from pillardisk.pillar_disk import PillarDisk, load_config, resolve_rin


def make_disk(config):
    dp = config.get('disk', {}).copy()
    tp = config.get('temperature', {})
    lp = config.get('lamp', {})
    op = config.get('observation', {})
    resolve_rin(dp)
    ratio = float(tp.get('tirrad_tvisc_ratio', 10.0))
    tv1 = float(tp.get('tv1', 1000.0))
    return PillarDisk(
        rin=float(dp.get('rin', 0.1)), rout=float(dp.get('rout', 50.0)),
        nr=int(dp.get('nr', 300)), nphi=int(dp.get('nphi', 360)),
        h1=float(tp.get('h1', 0.05)), r0=float(tp.get('r0', 50.0)),
        beta=float(tp.get('beta', 10.0)), tv1=tv1, alpha=float(tp.get('alpha', 0.75)),
        hlamp=float(lp.get('hlamp', 0.05)), tx1=tv1*ratio, tirrad_tvisc_ratio=ratio,
        dmpc=float(op.get('dmpc', 100.0)),
        cosi=np.cos(np.radians(float(op.get('inclination', 45.0)))),
        fcol=float(op.get('fcol', 1.0)), redshift=float(op.get('redshift', 0.0)),
    )


def _save_horizontal_cbar(outpath, phi_min=16.0, phi_max=21.5):
    """Standalone horizontal cbar PNG with the same colormap/range used
    by the 3D geometry plots, so the .tex can place it once below the
    two stacked panels."""
    from matplotlib.colors import Normalize
    from matplotlib import cm
    fig, ax = plt.subplots(figsize=(8.0, 0.5))
    norm = Normalize(vmin=phi_min, vmax=phi_max)
    mappable = cm.ScalarMappable(cmap=cm.inferno, norm=norm)
    mappable.set_array([])
    cbar = plt.colorbar(mappable, cax=ax, orientation='horizontal')
    cbar.set_label(r'$\log\Phi~\rm [photons~cm^{-2}~s^{-1}]$', fontsize=18)
    cbar.ax.tick_params(labelsize=14, direction='in')
    plt.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  -> {outpath}")


def main(config_file='config_line.yaml', outdir=None):
    config = load_config(config_file)
    if outdir is None:
        outdir = 'pillar_disc_draft/figures'
    os.makedirs(outdir, exist_ok=True)

    # --- Generate individual plots using existing method ---
    # Single pillar - read from config
    d1 = make_disk(config)
    pc = config.get('pillars', {})
    # Parse pi expressions
    def parse_math_expr(value):
        if isinstance(value, (int, float)):
            return float(value)
        if isinstance(value, str):
            try:
                return float(eval(value.strip(), {"__builtins__": {}},
                                  {'pi': np.pi, 'e': np.e}))
            except:
                return float(value)
        return float(value)

    r_list = pc.get('r_pillar', [4.0])
    phi_list = pc.get('phi_pillar', ['pi/4'])
    h_list = pc.get('height', [0.05])
    sr_list = pc.get('sigma_r', [0.5])
    sp_list = pc.get('sigma_phi', [0.5])
    mh_list = pc.get('modify_height', [True])
    mt_list = pc.get('modify_temp', [True])
    tf_list = pc.get('temp_factor', [1])
    if not isinstance(r_list, list): r_list = [r_list]
    if not isinstance(phi_list, list): phi_list = [phi_list]
    if not isinstance(h_list, list): h_list = [h_list]
    if not isinstance(sr_list, list): sr_list = [sr_list]
    if not isinstance(sp_list, list): sp_list = [sp_list]
    if not isinstance(mh_list, list): mh_list = [mh_list]
    if not isinstance(mt_list, list): mt_list = [mt_list]
    if not isinstance(tf_list, list): tf_list = [tf_list]
    # Only the FIRST pillar for the left panel: the paper's "single test
    # pillar" configuration is the first manual list entry, the same one
    # picked up by _build_disk_with_one_pillar.
    i = 0
    d1.add_pillar(
        r_pillar=float(r_list[i]),
        phi_pillar=parse_math_expr(phi_list[i]),
        height=float(h_list[i]),
        sigma_r=float(sr_list[i]),
        sigma_phi=parse_math_expr(sp_list[i]),
        modify_height=bool(mh_list[i % len(mh_list)]),
        modify_temp=bool(mt_list[i % len(mt_list)]),
        temp_factor=float(tf_list[i % len(tf_list)])
    )
    f1 = os.path.join(outdir, 'fig_geometry_1pillar.png')
    d1.plot_3d_geometry(color_by='flux', filename=f1,
                        show_pillar_markers=False,
                        show_colorbar=False,
                        n_phi_plot=720, n_r_plot=1000)

    # 100 pillars
    d2 = make_disk(config)
    np.random.seed(42)
    for _ in range(int(pc.get('N_pillar', 100))):
        r_p = -1
        while r_p < float(pc.get('rmin', 2.0)) or r_p > d2.rout:
            r_p = np.random.normal(float(pc.get('r_mean', 10.0)), float(pc.get('sig_r', 5.0)))
        d2.add_pillar(r_pillar=r_p, phi_pillar=np.random.uniform(0, 2*np.pi),
                      height=float(pc.get('h_pillar', 0.04)) * 0.5,
                      sigma_r=float(pc.get('sigma_r_pillar', 0.2)),
                      sigma_phi=float(pc.get('sigma_phi_pillar', 0.1)),
                      modify_height=True, modify_temp=True,
                      temp_factor=5.0)
    f2 = os.path.join(outdir, 'fig_geometry_N100.png')
    d2.plot_3d_geometry(color_by='flux', show_light_rays=False,
                        show_pillar_markers=False, filename=f2,
                        show_colorbar=False,
                        n_phi_plot=720, n_r_plot=1000)

    # Two separate plots; the paper layout stacks them top-bottom via
    # LaTeX (see fig:geometry in agn_pillars.tex). The shared colorbar
    # is a separate horizontal PNG, placed once below both panels in
    # the .tex.
    print(f"\n  -> {f1}")
    print(f"  -> {f2}")
    _save_horizontal_cbar(os.path.join(outdir, 'fig_geometry_cbar.png'),
                          phi_min=16.0, phi_max=21.5)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'config_line.yaml')
