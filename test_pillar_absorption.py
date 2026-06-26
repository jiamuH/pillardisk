"""
Plot emission line profile at different orbital phases showing absorption
from a pillar crossing the LOS to the BLR.

Uses the existing time-evolving code which already handles pillar shadows,
observer occultation, and Cloudy emissivities. We just integrate Psi(lambda, tau)
over tau to get F(lambda, t) at each epoch.
"""
import numpy as np
import matplotlib.pyplot as plt
from pillardisk.pillar_disk import PillarDisk, load_config, resolve_rin
from pillardisk.pillar_line import parse_math_expr, v_virial_from_mbh
from pillardisk.pillar_line_time_cloudy import (compute_time_evolving_cloudy_map,
                                      load_cloudy_models, C, KM_TO_CM)

config = load_config('config_line.yaml')


def build_disk():
    disk_params = config.get('disk', {}).copy()
    temp_params = config.get('temperature', {}).copy()
    lamp_params = config.get('lamp', {}).copy()
    obs_params = config.get('observation', {}).copy()
    float_params = ['rin', 'rout', 'h1', 'r0', 'beta', 'tv1', 'alpha',
                    'tx1', 'tirrad_tvisc_ratio', 'fcol', 'hlamp', 'dmpc',
                    'cosi', 'inclination', 'redshift', 'M_BH', 'r_isco_rg',
                    'f_trans', 'f_turb_line']
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


# --- Setup ---
disk = build_disk()

# Add one pillar
pillars_cfg = config.get('pillars', {})
r_p = 10.0
phi_p = 0.0
tf = float(pillars_cfg.get('temp_factor', [1])[0]) if isinstance(pillars_cfg.get('temp_factor', [1]), list) else float(pillars_cfg.get('temp_factor', 1))
disk.add_pillar(r_pillar=r_p, phi_pillar=phi_p, height=0.15,
                sigma_r=0.1, sigma_phi=0.05, temp_factor=tf,
                modify_height=True, modify_temp=True)
print(f"  1 narrow pillar: r={r_p}, sigma_r=0.1, sigma_phi=0.05")

M_BH = float(config.get('disk', {}).get('M_BH', 7e7))
r0 = float(config.get('temperature', {}).get('r0', 20.0))
v_virial = v_virial_from_mbh(M_BH, r0)

cloudy_file = '/Users/jiamuh/c23.01/my_models/loc_metal/strong_LOC_varym_N25_LineList_BLR_Fe2.txt'
cloudy_interp_dict, phi_grid, Z_grid = load_cloudy_models(cloudy_file, Z_target=1.0)

# Use only Halpha for this test
lambda0_dict = {'Halpha': 6562.8}

comp = config.get('computation', {})
nlambda = int(comp.get('nlambda', 200))
ntime = 8  # 8 epochs over one orbital period

# Compute orbital period at pillar radius
v_kep = v_virial * np.sqrt(r0 / r_p)
r_p_cm = r_p * 2.59e15  # ld to cm
v_kep_cm = v_kep * 1e5   # km/s to cm/s
T_orb = 2 * np.pi * r_p_cm / v_kep_cm / 86400.0  # days
print(f"Pillar at r={r_p} ld, v_Kep={v_kep:.0f} km/s, T_orb={T_orb:.1f} days")

# Run time-evolving computation
# Force f_turb_line=0 for the computation (keep Keplerian),
# then convolve with a Gaussian afterwards to get smooth broad emission
# while preserving the narrow absorption structure.
disk.f_turb_line = 0.0

print(f"Computing {ntime} epochs (Keplerian, will convolve after)...")
lambda_grids, time_grid, flux_arrays = compute_time_evolving_cloudy_map(
    disk, cloudy_interp_dict, phi_grid, lambda0_dict,
    v_virial, nlambda=nlambda, ntime=ntime, tmax=T_orb,
    use_absolute_flux=False, Z_target=1.0, n_cores=1)

lambda_grid = lambda_grids['Halpha']
flux_2d_raw = flux_arrays['Halpha']  # shape: (ntime, nlambda)
v_grid = (C / KM_TO_CM) * (lambda_grid - 6562.8) / 6562.8  # km/s

# Also compute baseline (no pillars)
print("Computing baseline (no pillars)...")
disk_nopillar = build_disk()
disk_nopillar.f_turb_line = 0.0
lambda_grids_np, _, flux_arrays_np = compute_time_evolving_cloudy_map(
    disk_nopillar, cloudy_interp_dict, phi_grid, lambda0_dict,
    v_virial, nlambda=nlambda, ntime=1, tmax=1.0,
    use_absolute_flux=False, Z_target=1.0, n_cores=1)
F_baseline_raw = flux_arrays_np['Halpha'][0, :]

# Post-process: convolve baseline with a broad Gaussian to get smooth BLR profile,
# but keep the pillar absorption narrow by only convolving the baseline,
# then scale the raw profiles to match.
from scipy.ndimage import gaussian_filter1d

# Convolution kernel width in pixels: f_turb * v_Kep(r_mean) / dv_per_pixel
dv = v_grid[1] - v_grid[0]  # km/s per pixel
f_turb_conv = float(config.get('lamp', {}).get('f_turb_line', 0.2))
v_kep_mean = v_virial * np.sqrt(r0 / 5.0)  # typical BLR radius
sigma_pix = f_turb_conv * v_kep_mean / dv
print(f"  Convolving with sigma = {f_turb_conv} * v_Kep(5 ld) = {f_turb_conv*v_kep_mean:.0f} km/s ({sigma_pix:.1f} pixels)")

F_baseline = gaussian_filter1d(F_baseline_raw, sigma_pix)

# Add narrow NLR component: fixed Gaussian at line center
# Typical NLR: FWHM ~ 2000 km/s, ~50% of peak BLR flux
fwhm_nlr_kms = 2000.0  # FWHM in km/s
sigma_nlr_kms = fwhm_nlr_kms / 2.355  # convert to sigma
sigma_nlr_pix = sigma_nlr_kms / dv
F_nlr_peak = 0.50 * np.max(F_baseline)
F_nlr = F_nlr_peak * np.exp(-0.5 * (v_grid / sigma_nlr_kms)**2)
print(f"  NLR: sigma={sigma_nlr_kms:.0f} km/s, peak={F_nlr_peak:.2e}")

F_baseline = F_baseline + F_nlr

# For each epoch: F_absorbed = F_baseline + (F_raw_pillar - F_raw_nopillar)
# The difference is narrow (shadow structure), added on top of smooth baseline
# NLR is constant so it's already in F_baseline
flux_2d = np.zeros_like(flux_2d_raw)
for i in range(ntime):
    delta = flux_2d_raw[i, :] - F_baseline_raw  # narrow absorption/enhancement
    flux_2d[i, :] = F_baseline + delta

# --- Plot ---
fig, axes = plt.subplots(2, 1, figsize=(10, 10), sharex=True)
fig.subplots_adjust(hspace=0.08)

colors = plt.cm.coolwarm(np.linspace(0, 1, ntime))

# Top: line profiles at different epochs
ax = axes[0]
ax.plot(v_grid, F_baseline, '-', color='black', lw=3, alpha=0.8, label=r'$\rm no~pillar$')
for i in range(ntime):
    phase_deg = 360.0 * time_grid[i] / T_orb
    ax.plot(v_grid, flux_2d[i, :], '-', color=colors[i], lw=2, alpha=0.7,
            label=rf'$\phi = {phase_deg:.0f}°$')
ax.set_ylabel(r'$F_\lambda~\rm [arb.]$', fontsize=16)
ax.set_title(r'$\rm H\alpha~profile~with~orbiting~pillar$', fontsize=16, pad=10)
ax.legend(fontsize=10, ncol=3, loc='upper right')
style_ax(ax)

# Bottom: residuals
ax = axes[1]
for i in range(ntime):
    phase_deg = 360.0 * time_grid[i] / T_orb
    residual = np.where(F_baseline > 0.01 * np.max(F_baseline),
                        (flux_2d[i, :] - F_baseline) / F_baseline, 0)
    ax.plot(v_grid, residual, '-', color=colors[i], lw=2, alpha=0.8,
            label=rf'$\phi = {phase_deg:.0f}°$')
ax.axhline(0, color='black', ls='-', lw=0.5, alpha=0.5)
ax.set_xlabel(r'$v~\rm [km\,s^{-1}]$', fontsize=16)
ax.set_ylabel(r'$\Delta F / F_0$', fontsize=16)
ax.set_title(r'$\rm Fractional~change~(pillar~at~r=10~ld)$', fontsize=16, pad=10)
ax.legend(fontsize=10, ncol=3, loc='lower left')
style_ax(ax)

plt.savefig('plots/test_pillar_absorption.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved plots/test_pillar_absorption.png")

# --- Compute fixed axis ranges across all frames ---
ylim_top = (0, np.max(F_baseline) * 1.15)
all_residuals = []
for i in range(ntime):
    res = np.where(F_baseline > 0.01 * np.max(F_baseline),
                   (flux_2d[i, :] - F_baseline) / F_baseline, 0)
    all_residuals.append(res)
res_min = min(np.min(r) for r in all_residuals)
res_max = max(np.max(r) for r in all_residuals)
ylim_bot = (res_min * 1.2, res_max * 1.2)

# --- Individual phase plots ---
for i in range(ntime):
    phase_deg = 360.0 * time_grid[i] / T_orb

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    fig.subplots_adjust(hspace=0.08)

    ax = axes[0]
    ax.plot(v_grid, F_baseline, '-', color='black', lw=3, alpha=0.8, label=r'$\rm no~pillar$')
    ax.plot(v_grid, flux_2d[i, :], '-', color='orangered', lw=3, alpha=0.8,
            label=rf'$\phi = {phase_deg:.0f}°$')
    ax.set_ylabel(r'$F_\lambda~\rm [arb.]$', fontsize=16)
    ax.set_title(rf'$\rm H\alpha~profile~at~\phi = {phase_deg:.0f}°$', fontsize=16, pad=10)
    ax.legend(fontsize=13)
    ax.set_ylim(ylim_top)
    style_ax(ax)

    ax = axes[1]
    ax.plot(v_grid, all_residuals[i], '-', color='orangered', lw=3, alpha=0.8)
    ax.axhline(0, color='black', ls='-', lw=0.5, alpha=0.5)
    ax.set_xlabel(r'$v~\rm [km\,s^{-1}]$', fontsize=16)
    ax.set_ylabel(r'$\Delta F / F_0$', fontsize=16)
    ax.set_ylim(ylim_bot)
    style_ax(ax)

    fname = f'plots/test_pillar_absorption_phi{int(phase_deg):03d}.png'
    plt.savefig(fname, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved {fname}")
