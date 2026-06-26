#!/usr/bin/env python3
"""
Regenerate all paper figures with updated NGC 5548 parameters.

Usage:
    python3 regenerate_figures.py [config_line.yaml]

Calls the existing main() functions from each module.
"""

import os
import shutil
import tempfile
import yaml


FIGDIR_BASE = 'pillar_disc_draft/figures'


def _run_cloudy_vdm(config_file, filename='velocity_delay_map.png'):
    """Compute Cloudy-weighted Hα VDM (with and without pillars) and plot."""
    from pillardisk.pillar_line import main as pillar_line_main, plot_velocity_delay_map
    from pillardisk.pillar_line_time_cloudy import compute_velocity_delay_map_cloudy, load_cloudy_models

    # Use pillar_line.main's disk setup by importing its internals
    from pillardisk.pillar_disk import PillarDisk, load_config, resolve_rin
    from pillardisk.pillar_line import parse_math_expr
    import numpy as np

    config = load_config(config_file)

    # --- Build disk (reuse pillar_line.py logic) ---
    disk_params = config.get('disk', {}).copy()
    temp_params = config.get('temperature', {}).copy()
    lamp_params = config.get('lamp', {}).copy()
    obs_params = config.get('observation', {}).copy()
    float_params = ['rin', 'rout', 'h1', 'r0', 'beta', 'tv1', 'alpha',
                    'tx1', 'tirrad_tvisc_ratio', 'fcol', 'hlamp', 'dmpc',
                    'cosi', 'inclination', 'redshift', 'M_BH', 'r_isco_rg', 'f_trans']
    int_params = ['nr', 'nphi']
    for d in [disk_params, temp_params, lamp_params, obs_params]:
        for k, v in d.items():
            if isinstance(v, str) and v.lower() == 'auto':
                continue
            if k in int_params:
                d[k] = int(float(v))
            elif k in float_params and v is not None:
                d[k] = float(v)
    resolve_rin(disk_params)
    if 'inclination' in obs_params and obs_params['inclination'] is not None:
        obs_params['cosi'] = np.cos(np.radians(obs_params['inclination']))
        del obs_params['inclination']
    disk = PillarDisk(**{**disk_params, **temp_params, **lamp_params, **obs_params})

    # --- Add pillars ---
    pillars_cfg = config.get('pillars', {})
    def to_list(x):
        return x if isinstance(x, list) else [x]
    r_list = to_list(pillars_cfg.get('r_pillar', []))
    phi_list = to_list(pillars_cfg.get('phi_pillar', []))
    h_list = to_list(pillars_cfg.get('height', [0.01]))
    sr_list = to_list(pillars_cfg.get('sigma_r', [1.0]))
    sp_list = to_list(pillars_cfg.get('sigma_phi', [0.1]))
    tf_list = to_list(pillars_cfg.get('temp_factor', [1.5]))
    mh_list = to_list(pillars_cfg.get('modify_height', [True]))
    mt_list = to_list(pillars_cfg.get('modify_temp', [False]))
    n_p = max(len(r_list), len(phi_list))
    def pad(lst, n, d):
        lst = lst if len(lst) >= n else lst * n
        return lst[:n]
    r_list = pad(r_list, n_p, 0); phi_list = pad(phi_list, n_p, 0)
    h_list = pad(h_list, n_p, 0.01); sr_list = pad(sr_list, n_p, 1.0)
    sp_list = pad(sp_list, n_p, 0.1); tf_list = pad(tf_list, n_p, 1.5)
    mh_list = pad(mh_list, n_p, True); mt_list = pad(mt_list, n_p, False)
    for i in range(n_p):
        disk.add_pillar(r_pillar=float(r_list[i]),
                        phi_pillar=parse_math_expr(phi_list[i]),
                        height=float(h_list[i]),
                        sigma_r=float(sr_list[i]),
                        sigma_phi=parse_math_expr(sp_list[i]),
                        temp_factor=float(tf_list[i]),
                        modify_height=bool(mh_list[i]),
                        modify_temp=bool(mt_list[i]))

    # --- Emission line params ---
    line_params = config.get('emission_line', {})
    lambda0 = 6562.8  # Hα
    from pillardisk.pillar_line import v_virial_from_mbh
    M_BH = float(line_params.get('M_BH', 7e7))
    v_virial = v_virial_from_mbh(M_BH, disk.r0)

    comp = config.get('computation', {})
    nlambda = int(comp.get('nlambda', 200))
    ntau = int(comp.get('ntau_line', comp.get('ntau', 100)))
    taumax = float(comp.get('taumax', 50.0))

    # --- Load Cloudy models (absolute flux + low-phi extension) ---
    cloudy_file = '/Users/jiamuh/c23.01/my_models/loc_metal_flux/strong_LOC_varym_N25_v100_lineflux_LineList_BLR_Fe2_flux.txt'
    cloudy_extension_file = '/Users/jiamuh/c23.01/my_models/loc_metal_flux/strong_LOC_varym_N25_v100_lineflux_extlow_LineList_BLR_Fe2_flux.txt'
    cloudy_interp_dict, phi_grid, Z_grid = load_cloudy_models(
        cloudy_file, Z_target=1.0, extension_file=cloudy_extension_file)

    # --- Compute Hα VDM with pillars ---
    print('Computing Cloudy Hα VDM with pillars...')
    lam, tau, psi = compute_velocity_delay_map_cloudy(
        disk, 'Halpha', lambda0, v_virial, cloudy_interp_dict,
        nlambda=nlambda, ntau=ntau, taumax=taumax,
        weighting='responsivity')

    # --- Compute per-pillar VDMs ---
    pillars_bak = disk.pillars.copy()
    psi_per_pillar = []
    for ip, p in enumerate(pillars_bak):
        print(f'Computing Cloudy Hα VDM for pillar {ip+1} alone...')
        disk.pillars = [p]
        _, _, psi_pp = compute_velocity_delay_map_cloudy(
            disk, 'Halpha', lambda0, v_virial, cloudy_interp_dict,
            nlambda=nlambda, ntau=ntau, taumax=taumax,
            weighting='responsivity')
        psi_per_pillar.append(psi_pp)

    # --- Compute Hα VDM without pillars ---
    print('Computing Cloudy Hα VDM without pillars...')
    disk.pillars = []
    _, _, psi_no = compute_velocity_delay_map_cloudy(
        disk, 'Halpha', lambda0, v_virial, cloudy_interp_dict,
        nlambda=nlambda, ntau=ntau, taumax=taumax,
        weighting='responsivity')
    disk.pillars = pillars_bak

    # --- Plot ---
    plot_velocity_delay_map(lam, tau, psi, lambda0,
                            psi_map_no_pillars=psi_no,
                            psi_per_pillar=psi_per_pillar,
                            filename=filename)
    print(f'Saved Cloudy Hα VDM to {filename}')


def make_config(config_file, make_many):
    """Load config and return a temp file with make_many toggled."""
    with open(config_file) as f:
        config = yaml.safe_load(f)
    config['pillars']['make_many'] = make_many
    # Also disable 3D geometry plot in pillar_line.py to save time
    config.setdefault('plotting', {})['plot_3d_geometry'] = False
    tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False)
    yaml.dump(config, tmp, default_flow_style=False)
    tmp.close()
    return tmp.name


def move(src, dst):
    """Move file if it exists."""
    if os.path.exists(src):
        shutil.move(src, dst)
        print(f'  -> {dst}')
    else:
        print(f'  [skip] {src} not found')


def cleanup_geometry_snapshots():
    """Remove orbital-phase geometry PNGs from plots/ (and cwd as fallback)."""
    for name in ['geometry_phi0_near.png', 'geometry_phi90_quad1.png',
                 'geometry_phi180_far.png', 'geometry_phi270_quad2.png']:
        for f in (os.path.join('plots', name), name):
            if os.path.exists(f):
                os.remove(f)


def _make_figdir(config_file):
    """Build a parameter-tagged output directory so runs don't overwrite."""
    with open(config_file) as f:
        cfg = yaml.safe_load(f)
    beta = cfg.get('temperature', {}).get('beta', '?')
    phi_in = cfg.get('lamp', {}).get('log_phi_inner', '?')
    nff = cfg.get('lamp', {}).get('no_fluxfloor', False)
    hlamp = cfg.get('lamp', {}).get('hlamp', '?')
    incl = cfg.get('observation', {}).get('inclination', '?')
    h_pillar = cfg.get('pillars', {}).get('h_pillar', '?')
    phi_min = 'noff' if nff else 'floor'
    transparent = cfg.get('lamp', {}).get('transparent', False)
    f_trans = cfg.get('lamp', {}).get('f_trans', 0.0)
    tag = (f'beta{beta}_phi_in{phi_in}_phi_min{phi_min}'
           f'_hlamp{hlamp}_hp{h_pillar}_inc{incl}')
    if transparent and f_trans != 0:
        tag += f'_f_trans{f_trans}'
    return os.path.join(FIGDIR_BASE, tag)


def main(config_file='config_line.yaml'):
    FIGDIR = _make_figdir(config_file)
    os.makedirs(FIGDIR, exist_ok=True)
    print(f'Output directory: {FIGDIR}')

    # ── Quick plots first (t=0 snapshots) ────────────────

    # ── 1. 3D geometry two-panel ──────────────────────────
    print('\n' + '='*50)
    print(' [1/9] 3D geometry two-panel figure')
    print('='*50)
    from plot_geometry_figures import main as geom_main
    geom_main(config_file=config_file, outdir=FIGDIR)

    # ── 2. Single-pillar ionflux (t=0 face-on) ───────────
    print('\n' + '='*50)
    print(' [2/9] Single-pillar ionflux map (plot-only)')
    print('='*50)
    tmp = make_config(config_file, make_many=False)
    try:
        from pillardisk.pillar_line_time_cloudy import main as time_cloudy_main
        time_cloudy_main(config_file=tmp, plot_only=True)
    finally:
        os.unlink(tmp)
    move('plots/geometry_t0_faceon.png', f'{FIGDIR}/fig_ionflux_1pillar.png')
    cleanup_geometry_snapshots()

    # ── 3. N100 ionflux (t=0 face-on) ────────────────────
    print('\n' + '='*50)
    print(' [3/9] 100-pillar ionflux map (plot-only)')
    print('='*50)
    tmp = make_config(config_file, make_many=True)
    try:
        time_cloudy_main(config_file=tmp, plot_only=True)
    finally:
        os.unlink(tmp)
    move('plots/geometry_t0_faceon.png', f'{FIGDIR}/fig_ionflux_N100.png')
    cleanup_geometry_snapshots()

    # ── 4. Single-pillar (two-pillar) RGB / line-emission map ─
    # Disabled to save time: the single/two-pillar RGB map is not used
    # in the current paper figures — only the 2-pillar VDM (step 6) and
    # the 100-pillar RGB (step 5) are kept. Re-enable if needed.
    # print('\n' + '='*50)
    # print(' [4/9] Single-pillar RGB map (pillar_line_cloudy.py)')
    # print('='*50)
    # tmp = make_config(config_file, make_many=False)
    # try:
    #     from pillardisk.pillar_line_cloudy import main as cloudy_main
    #     cloudy_main(config_file=tmp, use_absolute_flux=False)
    # finally:
    #     os.unlink(tmp)
    # for suffix in ['_xy', '_rphi', '_ionizing_flux', '_logF_logphi', '_EW_xy', '_EW_rphi']:
    #     move(f'plots/line_intensity_rgb{suffix}.png', f'{FIGDIR}/fig_rgb_1pillar{suffix}.png')

    # ── 5. 100-pillar RGB map ─────────────────────────────
    print('\n' + '='*50)
    print(' [5/9] 100-pillar RGB map (pillar_line_cloudy.py)')
    print('='*50)
    tmp = make_config(config_file, make_many=True)
    try:
        from pillardisk.pillar_line_cloudy import main as cloudy_main
        cloudy_main(config_file=tmp, use_absolute_flux=False)
    finally:
        os.unlink(tmp)
    # r-phi panels disabled to save iteration time; only x-y panels and
    # the ionizing-flux / logF-logphi maps are kept.
    for suffix in ['_xy', '_ionizing_flux', '_logF_logphi', '_EW_xy']:
        move(f'plots/line_intensity_rgb{suffix}.png', f'{FIGDIR}/fig_rgb_N100{suffix}.png')

    # ── 6. Single-pillar Cloudy-weighted Hα VDM (slower) ──
    # Pulled out into make_vdm_1pillar.py for faster iteration on this one
    # figure. Run `python3 -m pillardisk.make_vdm_1pillar` to regenerate it.

    # ── 7. Movie frame at t=0 (with geometry) ────────────
    # Disabled to save time: only the VDM-only movie frame (step 8)
    # is needed for the paper figure (fig_vdm_frame). Re-enable if the
    # geometry-included version is wanted.
    # print('\n' + '='*50)
    # print(' [7/10] Movie frame at t=0 (with geometry)')
    # print('='*50)
    # tmp = make_config(config_file, make_many=True)
    # try:
    #     from pillardisk.pillar_line_time_cloudy import generate_movie_frames
    #     generate_movie_frames(config_file=tmp, n_frames=1, output_dir='plots/_tmp_movie')
    # finally:
    #     os.unlink(tmp)
    # move('plots/_tmp_movie/frame_0000.png', f'{FIGDIR}/fig_movie_frame.png')
    # if os.path.exists('plots/_tmp_movie'):
    #     shutil.rmtree('plots/_tmp_movie')

    # ── 8. Movie frame at t=0 (VDM only, no geometry) ───
    print('\n' + '='*50)
    print(' [8/10] Movie frame at t=0 (VDM only)')
    print('='*50)
    tmp = make_config(config_file, make_many=True)
    try:
        from pillardisk.pillar_line_time_cloudy import generate_movie_frames
        generate_movie_frames(config_file=tmp, n_frames=1, output_dir='plots/_tmp_movie',
                              skip_geometry=True)
    finally:
        os.unlink(tmp)
    move('plots/_tmp_movie/frame_0000.png', f'{FIGDIR}/fig_vdm_frame.png')
    if os.path.exists('plots/_tmp_movie'):
        shutil.rmtree('plots/_tmp_movie')

    # ── Done ──────────────────────────────────────────────
    # Time-evolving (slow) maps live in make_time_evolving_figures.py.
    print('\n' + '='*50)
    print(f' All fast figures saved to {FIGDIR}/')
    print(' Time-evolving maps: run make_time_evolving_figures.py separately.')
    print('='*50)
    for f in sorted(os.listdir(FIGDIR)):
        if f.startswith('fig_') and f.endswith('.png'):
            print(f'  {f}')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Regenerate paper figures (fast plots only)')
    parser.add_argument('config', nargs='?', default='config_line.yaml',
                        help='Config YAML file (default: config_line.yaml)')
    args = parser.parse_args()
    main(config_file=args.config)
