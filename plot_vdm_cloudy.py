#!/usr/bin/env python3
"""
Plot Cloudy-weighted Hα velocity-delay map for paper figure.

Usage:
    python3 plot_vdm_cloudy.py [config_line.yaml]

Produces fig_vdm_1pillar.png with 4 panels:
  - Top left:  VDM with pillars
  - Top right: VDM without pillars
  - Bottom left: Mean delay (all pillars, per-pillar, no pillars)
  - Bottom right: Difference map
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


from pillardisk.pillar_disk import PillarDisk, load_config

VDM_CMAP = LinearSegmentedColormap.from_list('vdm_wrb',
                                              ['white', 'red', 'blue'])
from pillardisk.pillar_line import parse_math_expr, v_virial_from_mbh, C, KM_TO_CM
from pillardisk.pillar_line_time_cloudy import compute_velocity_delay_map_cloudy, load_cloudy_models

CLOUDY_FILE = '/Users/jiamuh/c23.01/my_models/loc_metal/strong_LOC_varym_N25_LineList_BLR_Fe2.txt'
FIGDIR = 'pillar_disc_draft/figures'


def build_disk(config):
    """Create PillarDisk from config dict and add pillars."""
    disk_params = config.get('disk', {}).copy()
    temp_params = config.get('temperature', {}).copy()
    lamp_params = config.get('lamp', {}).copy()
    obs_params = config.get('observation', {}).copy()

    float_keys = ['rin', 'rout', 'h1', 'r0', 'beta', 'tv1', 'alpha',
                  'tx1', 'tirrad_tvisc_ratio', 'fcol', 'hlamp', 'dmpc',
                  'cosi', 'inclination', 'redshift']
    int_keys = ['nr', 'nphi']
    for d in [disk_params, temp_params, lamp_params, obs_params]:
        for k, v in d.items():
            if v is None:
                continue
            if k in int_keys:
                d[k] = int(float(v))
            elif k in float_keys:
                d[k] = float(v)

    if 'inclination' in obs_params and obs_params['inclination'] is not None:
        obs_params['cosi'] = np.cos(np.radians(obs_params['inclination']))
        del obs_params['inclination']

    disk = PillarDisk(**{**disk_params, **temp_params, **lamp_params, **obs_params})

    # Add pillars
    pcfg = config.get('pillars', {})
    def to_list(x):
        return x if isinstance(x, list) else [x]
    r_list = to_list(pcfg.get('r_pillar', []))
    phi_list = to_list(pcfg.get('phi_pillar', []))
    h_list = to_list(pcfg.get('height', [0.01]))
    sr_list = to_list(pcfg.get('sigma_r', [1.0]))
    sp_list = to_list(pcfg.get('sigma_phi', [0.1]))
    tf_list = to_list(pcfg.get('temp_factor', [1.5]))
    mh_list = to_list(pcfg.get('modify_height', [True]))
    mt_list = to_list(pcfg.get('modify_temp', [False]))
    n_p = max(len(r_list), len(phi_list))
    def pad(lst, n, d):
        return (lst * n)[:n] if len(lst) < n else lst[:n]
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
    return disk


def compute_mean_delay(psi_map, tau_grid):
    """Compute mean delay <tau>(lambda) from a VDM."""
    nlam = psi_map.shape[1]
    tau_mean = np.zeros(nlam)
    for i in range(nlam):
        col = psi_map[:, i]
        if np.any(col > 0):
            tau_mean[i] = np.sum(col * tau_grid) / np.sum(col)
    return tau_mean


def plot_vdm(lambda_grid, tau_grid, psi, psi_no, psi_per_pillar,
             lambda0, filename):
    """Plot 4-panel VDM figure."""
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # Common color scale
    vmin = min(np.min(psi[psi > 0]), np.min(psi_no[psi_no > 0]))
    vmax = max(np.max(psi), np.max(psi_no))

    v_km_s = (C / KM_TO_CM) * (lambda_grid - lambda0) / lambda0

    def add_velocity_axis(ax):
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xlabel(r'$v~\rm [km\,s^{-1}]$', labelpad=10)
        v_ticks = np.linspace(v_km_s[0], v_km_s[-1], 5)
        lam_ticks = lambda0 * (1.0 + v_ticks * KM_TO_CM / C)
        ax2.set_xticks(lam_ticks)
        ax2.set_xticklabels([rf'${v:.0f}$' for v in v_ticks])
        ax2.minorticks_on()
        ax2.tick_params(top=True, right=True, which='major', length=8,
                        width=1.5, direction='in')
        ax2.tick_params(top=True, right=True, which='minor', length=4,
                        width=1, direction='in')

    def style_ax(ax):
        ax.minorticks_on()
        ax.tick_params(top=True, right=True, which='major', length=8,
                       width=1.5, direction='in')
        ax.tick_params(top=True, right=True, which='minor', length=4,
                       width=1, direction='in')

    # --- Top left: with pillars ---
    ax = axes[0, 0]
    im1 = ax.pcolormesh(lambda_grid, tau_grid, psi, shading='auto',
                        cmap=VDM_CMAP, vmin=vmin, vmax=vmax)
    ax.axvline(lambda0, color='red', linestyle='--', linewidth=2)
    ax.set_xlabel(r'$\lambda~\rm [\AA]$')
    ax.set_ylabel(r'$\tau~\rm [days]$')
    ax.set_title(r'$\rm With~Pillars~\Psi(\lambda,\tau)$', pad=10)
    ax.set_xlim(lambda_grid[0], lambda_grid[-1])
    ax.set_ylim(0, tau_grid[-1])
    plt.colorbar(im1, ax=ax, label=r'$\Psi(\lambda,\tau)$',
                 shrink=0.8, aspect=20)
    style_ax(ax)
    add_velocity_axis(ax)

    # --- Top right: without pillars ---
    ax = axes[0, 1]
    im2 = ax.pcolormesh(lambda_grid, tau_grid, psi_no, shading='auto',
                        cmap=VDM_CMAP, vmin=vmin, vmax=vmax)
    ax.axvline(lambda0, color='red', linestyle='--', linewidth=2)
    ax.set_xlabel(r'$\lambda~\rm [\AA]$')
    ax.set_ylabel(r'$\tau~\rm [days]$')
    ax.set_title(r'$\rm Without~Pillars~\Psi(\lambda,\tau)$', pad=10)
    ax.set_xlim(lambda_grid[0], lambda_grid[-1])
    ax.set_ylim(0, tau_grid[-1])
    plt.colorbar(im2, ax=ax, label=r'$\Psi(\lambda,\tau)$',
                 shrink=0.8, aspect=20)
    style_ax(ax)
    add_velocity_axis(ax)

    # --- Bottom left: mean delay ---
    ax = axes[1, 0]
    tau_mean = compute_mean_delay(psi, tau_grid)
    tau_mean_no = compute_mean_delay(psi_no, tau_grid)
    ax.plot(lambda_grid, tau_mean, 'b-', linewidth=2,
            label=r'$\rm all~pillars$')
    ax.plot(lambda_grid, tau_mean_no, 'k--', linewidth=2,
            label=r'$\rm no~pillars$')

    # Per-pillar contributions
    colors_pp = ['green', 'orange', 'purple', 'cyan', 'magenta']
    for ip, psi_pp in enumerate(psi_per_pillar):
        tau_pp = compute_mean_delay(psi_pp, tau_grid)
        c = colors_pp[ip % len(colors_pp)]
        ax.plot(lambda_grid, tau_pp, '-', color=c, linewidth=1.5,
                label=rf'$\rm pillar~{ip+1}$')

    ax.axvline(lambda0, color='red', linestyle=':', linewidth=1, alpha=0.5)
    ax.set_xlabel(r'$\lambda~\rm [\AA]$')
    ax.set_ylabel(r'$\langle\tau\rangle~\rm [days]$')
    ax.set_xlim(lambda_grid[0], lambda_grid[-1])
    all_tau = np.concatenate([tau_mean[tau_mean > 0],
                              tau_mean_no[tau_mean_no > 0]])
    if len(all_tau) > 0:
        ax.set_ylim(0, 1.3 * np.max(all_tau))
    ax.legend(fontsize=10)
    style_ax(ax)
    add_velocity_axis(ax)

    # --- Bottom right: difference ---
    ax = axes[1, 1]
    diff = psi - psi_no
    diff_vmax = float(np.nanmax(np.abs(diff))) if np.any(np.isfinite(diff)) else 1.0
    im3 = ax.pcolormesh(lambda_grid, tau_grid, diff, shading='auto',
                        cmap='RdBu', vmin=-diff_vmax, vmax=diff_vmax)
    ax.axvline(lambda0, color='black', linestyle='--', linewidth=2)
    ax.set_xlabel(r'$\lambda~\rm [\AA]$')
    ax.set_ylabel(r'$\tau~\rm [days]$')
    ax.set_title(r'$\rm Difference~\Psi_{\rm with} - \Psi_{\rm without}$')
    ax.set_xlim(lambda_grid[0], lambda_grid[-1])
    ax.set_ylim(0, tau_grid[-1])
    plt.colorbar(im3, ax=ax, label=r'$\Delta\Psi$',
                 shrink=0.8, aspect=20)
    style_ax(ax)
    add_velocity_axis(ax)

    plt.tight_layout()
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {filename}')


def main(config_file='config_line.yaml'):
    config = load_config(config_file)
    disk = build_disk(config)

    # Emission line params
    line_params = config.get('emission_line', {})
    lambda0 = 6562.8  # Hα
    M_BH = float(line_params.get('M_BH', 7e7))
    v_virial = v_virial_from_mbh(M_BH, disk.r0)

    comp = config.get('computation', {})
    nlambda = int(comp.get('nlambda', 200))
    ntau = int(comp.get('ntau_line', comp.get('ntau', 100)))
    taumax = float(comp.get('taumax', 50.0))

    # Load Cloudy
    print('Loading Cloudy models...')
    cloudy_interp_dict, phi_grid, Z_grid = load_cloudy_models(
        CLOUDY_FILE, Z_target=1.0)

    kw = dict(nlambda=nlambda, ntau=ntau, taumax=taumax)

    # All pillars
    print('Computing Hα VDM with all pillars...')
    lam, tau, psi = compute_velocity_delay_map_cloudy(
        disk, 'Halpha', lambda0, v_virial, cloudy_interp_dict, **kw)

    # Per-pillar
    pillars_bak = disk.pillars.copy()
    psi_per_pillar = []
    for ip, p in enumerate(pillars_bak):
        print(f'Computing Hα VDM for pillar {ip+1} alone...')
        disk.pillars = [p]
        _, _, psi_pp = compute_velocity_delay_map_cloudy(
            disk, 'Halpha', lambda0, v_virial, cloudy_interp_dict, **kw)
        psi_per_pillar.append(psi_pp)

    # No pillars
    print('Computing Hα VDM without pillars...')
    disk.pillars = []
    _, _, psi_no = compute_velocity_delay_map_cloudy(
        disk, 'Halpha', lambda0, v_virial, cloudy_interp_dict, **kw)
    disk.pillars = pillars_bak

    # Plot — filename encodes key parameters so runs don't overwrite
    os.makedirs(FIGDIR, exist_ok=True)
    beta = config.get('temperature', {}).get('beta', '?')
    phi_in = config.get('lamp', {}).get('log_phi_inner', '?')
    nff = config.get('lamp', {}).get('no_fluxfloor', False)
    hlamp = config.get('lamp', {}).get('hlamp', '?')
    incl = config.get('observation', {}).get('inclination', '?')
    phi_min = 'noff' if nff else 'floor'
    transparent = config.get('lamp', {}).get('transparent', False)
    f_trans = config.get('lamp', {}).get('f_trans', 0.0)
    tag = f'_beta{beta}_phi_in{phi_in}_phi_min{phi_min}_hlamp{hlamp}_inc{incl}'
    if transparent and f_trans != 0:
        tag += f'_f_trans{f_trans}'
    filename = os.path.join(FIGDIR, f'fig_vdm_1pillar{tag}.png')
    plot_vdm(lam, tau, psi, psi_no, psi_per_pillar, lambda0, filename)


if __name__ == '__main__':
    cfg = sys.argv[1] if len(sys.argv) > 1 else 'config_line.yaml'
    main(config_file=cfg)
