"""
Plot continuum lag spectrum and SED at different orbital epochs.
As pillars rotate, the shadow pattern changes with viewing angle,
modifying both the SED and lag spectrum.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from pillardisk.pillar_disk import PillarDisk, load_config, resolve_rin
from pillardisk.pillar_line import v_virial_from_mbh

config = load_config('config_line.yaml')

wavelengths = np.logspace(np.log10(1000), np.log10(10000), 50)
BETA = 50.0


def build_disk(hlamp=0.5):
    disk_params = config.get('disk', {}).copy()
    temp_params = config.get('temperature', {}).copy()
    lamp_params = config.get('lamp', {}).copy()
    obs_params = config.get('observation', {}).copy()
    temp_params['beta'] = BETA
    lamp_params['hlamp'] = hlamp
    float_params = ['rin', 'rout', 'h1', 'r0', 'beta', 'tv1', 'alpha',
                    'tx1', 'tirrad_tvisc_ratio', 'fcol', 'hlamp', 'dmpc',
                    'cosi', 'inclination', 'redshift', 'M_BH', 'r_isco_rg', 'f_trans',
                    'f_turb_line']
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
    if 'inclination' in obs_params:
        obs_params['cosi'] = np.cos(np.radians(obs_params['inclination']))
        del obs_params['inclination']
    return PillarDisk(**{**disk_params, **temp_params, **lamp_params, **obs_params})


def style_ax(ax):
    ax.tick_params(direction='in', which='major', length=8, width=1.8, labelsize=13,
                   top=True, right=True)
    ax.tick_params(direction='in', which='minor', length=5, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()


def scan_colors(n, cmap='plasma'):
    return [plt.colormaps[cmap](x) for x in np.linspace(0.15, 0.85, n)]


if __name__ == '__main__':
    DATAFILE = 'plots/lag_sed_epochs_data.npz'

    # Pillar config
    h_p, sr, sp = 0.5, 0.5, 0.3
    hlamp = 0.5
    r_p_mean = 15.0
    N_pillars = 100
    n_epochs = 6

    if not os.path.exists(DATAFILE):
        # Compute orbital period at mean pillar radius
        M_BH = float(config.get('disk', {}).get('M_BH', 7e7))
        r0 = float(config.get('temperature', {}).get('r0', 20.0))
        v_virial = v_virial_from_mbh(M_BH, r0)
        v_kep = v_virial * np.sqrt(r0 / r_p_mean) * 1e5  # cm/s
        r_cm = r_p_mean * 2.59e15
        T_orb = 2 * np.pi * r_cm / v_kep / 86400.0  # days

        # Phase angles
        dphi_epochs = np.linspace(0, 2*np.pi, n_epochs, endpoint=False)

        # Bowl baseline (no pillars)
        print("Computing bowl baseline...")
        disk_bowl = build_disk(hlamp=hlamp)
        tau_bowl, tau_grid, psi_bowl = disk_bowl.compute_lag_spectrum(
            wavelengths, ntau=500, taumax=50.0, parallel=True)
        sed_bowl = np.array([disk_bowl.compute_sed(wavelengths[iw:iw+1])[0]
                             for iw in range(len(wavelengths))])

        # Compute at each epoch
        all_tau = []
        all_psi = []
        all_sed = []

        # Build base pillar positions once
        rng = np.random.default_rng(42)
        r_vals = rng.normal(r_p_mean, 5.0, N_pillars)
        r_vals = np.clip(r_vals, 2.0, 20.0)
        phi_vals = rng.uniform(0, 2*np.pi, N_pillars)

        for ie, dphi in enumerate(dphi_epochs):
            phase_deg = np.degrees(dphi)
            print(f"Epoch {ie+1}/{n_epochs}: phase = {phase_deg:.0f}°")

            disk = build_disk(hlamp=hlamp)
            # Rotate all pillars by dphi (rigid rotation for simplicity)
            for i in range(N_pillars):
                disk.add_pillar(r_pillar=r_vals[i],
                                phi_pillar=np.mod(phi_vals[i] + dphi, 2*np.pi),
                                height=h_p, sigma_r=sr, sigma_phi=sp,
                                temp_factor=1.0, modify_height=True, modify_temp=False)

            tau_mean, _, psi = disk.compute_lag_spectrum(
                wavelengths, ntau=500, taumax=50.0, parallel=True)
            sed = np.array([disk.compute_sed(wavelengths[iw:iw+1])[0]
                            for iw in range(len(wavelengths))])

            all_tau.append(tau_mean)
            all_psi.append(psi)
            all_sed.append(sed)

        np.savez(DATAFILE,
                 wavelengths=wavelengths,
                 tau_grid=tau_grid,
                 tau_bowl=tau_bowl, psi_bowl=psi_bowl, sed_bowl=sed_bowl,
                 dphi_epochs=dphi_epochs, T_orb=T_orb,
                 all_tau=np.array(all_tau),
                 all_psi=np.array(all_psi),
                 all_sed=np.array(all_sed))
        print(f"Saved data to {DATAFILE}")
    else:
        print(f"Loading cached data from {DATAFILE}")

    # Load data
    data = np.load(DATAFILE)
    wavelengths = data['wavelengths']
    tau_grid = data['tau_grid']
    tau_bowl = data['tau_bowl']
    sed_bowl = data['sed_bowl']
    dphi_epochs = data['dphi_epochs']
    all_tau = data['all_tau']
    all_sed = data['all_sed']
    n_epochs = len(dphi_epochs)

    C_cgs = 2.998e10
    ANGSTROM = 1e-8
    nu_grid = C_cgs / (wavelengths * ANGSTROM)
    colors = scan_colors(n_epochs)

    # --- Plot 1: tau(lambda) at different epochs ---
    fig, axes = plt.subplots(2, 1, figsize=(10, 10), sharex=True)
    fig.subplots_adjust(hspace=0.08)

    ax = axes[0]
    for ie in range(n_epochs):
        phase_deg = np.degrees(dphi_epochs[ie])
        ax.plot(wavelengths, all_tau[ie], lw=3, alpha=0.7, color=colors[ie],
                label=rf'$\phi = {phase_deg:.0f}°$')
    ax.plot(wavelengths, tau_bowl, '--', color='dodgerblue', lw=4, alpha=1.0,
            label=r'$\rm bowl$')
    ax.set_ylabel(r'$\tau~\rm [days]$', fontsize=16)
    ax.set_title(r'$\rm Continuum~lag~at~different~orbital~phases$', fontsize=16, pad=10)
    ax.legend(fontsize=11, ncol=2)
    ax.set_xscale('log')
    style_ax(ax)

    # Relative to 1500 A
    ax = axes[1]
    i_ref = np.argmin(np.abs(wavelengths - 1500.0))
    for ie in range(n_epochs):
        phase_deg = np.degrees(dphi_epochs[ie])
        ax.plot(wavelengths, all_tau[ie] - all_tau[ie][i_ref], lw=3, alpha=0.7,
                color=colors[ie], label=rf'$\phi = {phase_deg:.0f}°$')
    ax.plot(wavelengths, tau_bowl - tau_bowl[i_ref], '--', color='dodgerblue', lw=4,
            alpha=1.0, label=r'$\rm bowl$')
    ax.set_xlabel(r'$\lambda~\rm [\AA]$', fontsize=16)
    ax.set_ylabel(r'$\tau(\lambda) - \tau(1500\,\rm \AA)~\rm [days]$', fontsize=16)
    ax.set_xscale('log')
    style_ax(ax)

    plt.savefig('plots/test_lag_sed_epochs_tau.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved plots/test_lag_sed_epochs_tau.png")

    # --- Plot 2: SED at different epochs ---
    fig, ax = plt.subplots(1, 1, figsize=(10, 7))
    for ie in range(n_epochs):
        phase_deg = np.degrees(dphi_epochs[ie])
        ax.plot(wavelengths, nu_grid * all_sed[ie], lw=3, alpha=0.7, color=colors[ie],
                label=rf'$\phi = {phase_deg:.0f}°$')
    ax.plot(wavelengths, nu_grid * sed_bowl, '--', color='dodgerblue', lw=4, alpha=1.0,
            label=r'$\rm bowl$')
    ax.set_xlabel(r'$\lambda~\rm [\AA]$', fontsize=16)
    ax.set_ylabel(r'$\nu F_\nu~\rm [arb.]$', fontsize=16)
    ax.set_title(r'$\rm SED~at~different~orbital~phases$', fontsize=16, pad=10)
    ax.legend(fontsize=11, ncol=2)
    ax.set_xscale('log')
    ax.set_yscale('log')
    style_ax(ax)

    plt.savefig('plots/test_lag_sed_epochs_sed.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved plots/test_lag_sed_epochs_sed.png")
