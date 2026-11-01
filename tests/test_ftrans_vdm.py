"""Compare VDM for f_trans=0.1 vs f_trans=0.9 to check if transparent flag affects VDM."""
import tempfile, yaml
import numpy as np

from pillardisk.pillar_disk import PillarDisk, load_config, resolve_rin
from pillardisk.pillar_line import parse_math_expr, VDM_CMAP
from pillardisk.pillar_line_time_cloudy import compute_velocity_delay_map_cloudy, load_cloudy_models
import matplotlib.pyplot as plt

config = load_config('config_line.yaml')

def make_disk(f_trans_val):
    """Build disk with given f_trans."""
    disk_params = config.get('disk', {}).copy()
    temp_params = config.get('temperature', {}).copy()
    lamp_params = config.get('lamp', {}).copy()
    obs_params = config.get('observation', {}).copy()

    lamp_params['transparent'] = True
    lamp_params['f_trans'] = f_trans_val

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
    if 'inclination' in obs_params:
        obs_params['cosi'] = np.cos(np.radians(obs_params['inclination']))
        del obs_params['inclination']

    disk = PillarDisk(**{**disk_params, **temp_params, **lamp_params, **obs_params})
    print(f"  f_trans={disk.f_trans}, transparent={disk.transparent}")

    # Add pillars
    pillars_cfg = config.get('pillars', {})
    for phi_p in [0, 2*np.pi/3, 4*np.pi/3]:
        disk.add_pillar(r_pillar=4.0, phi_pillar=phi_p,
                        height=0.15, sigma_r=0.5, sigma_phi=0.2,
                        temp_factor=1, modify_height=True, modify_temp=True)
    return disk


# Load Cloudy models (absolute-flux grid + low-phi extension, as used by
# the main VDM figures)
cloudy_file = '/Users/jiamuh/c23.01/my_models/loc_metal_flux/strong_LOC_varym_N25_v100_lineflux_LineList_BLR_Fe2_flux.txt'
cloudy_extension_file = '/Users/jiamuh/c23.01/my_models/loc_metal_flux/strong_LOC_varym_N25_v100_lineflux_extlow_LineList_BLR_Fe2_flux.txt'
cloudy_interp_dict, phi_grid, Z_grid = load_cloudy_models(
    cloudy_file, Z_target=1.0, extension_file=cloudy_extension_file)

lambda0 = 6562.8  # Halpha
from pillardisk.pillar_line import v_virial_from_mbh
M_BH = 7e7
r0 = config.get('temperature', {}).get('r0', 20.0)
v_virial = v_virial_from_mbh(M_BH, float(r0))

comp = config.get('computation', {})
nlambda = int(comp.get('nlambda', 200))
ntau = int(comp.get('ntau_line', 100))
taumax = float(comp.get('taumax', 50.0))

# Compute VDM for two f_trans values
results = {}
disks = {}
for ft in [0.1, 0.6]:
    print(f"\nComputing VDM with f_trans={ft}...")
    disk = make_disk(ft)
    disks[ft] = disk
    lam, tau, psi = compute_velocity_delay_map_cloudy(
        disk, 'Halpha', lambda0, v_virial, cloudy_interp_dict,
        nlambda=nlambda, ntau=ntau, taumax=taumax)
    results[ft] = (lam, tau, psi)
    print(f"  psi sum = {np.sum(psi):.2f}, max = {np.max(psi):.2f}")

# Plot comparison
lam1, tau1, psi1 = results[0.1]
lam2, tau2, psi2 = results[0.6]

vmax = max(np.max(psi1), np.max(psi2))

fig = plt.figure(figsize=(13, 11))

# --- Top row: face-on ionizing flux ---
r_plot = np.linspace(disks[0.1].rin, disks[0.1].rout, 200)
phi_plot = np.linspace(0, 2*np.pi, 360)
r_2d_p, phi_2d_p = np.meshgrid(r_plot, phi_plot, indexing='ij')

for idx, ft in enumerate([0.1, 0.6]):
    ax = fig.add_subplot(2, 2, idx + 1, projection='polar')
    log_phi_2d = disks[ft].compute_log_ionizing_flux(r_2d_p, phi_2d_p)
    c = ax.pcolormesh(phi_2d_p, r_2d_p, log_phi_2d, shading='auto', cmap='inferno',
                      vmin=16.0, vmax=21.5)
    ax.set_title(rf'$f_{{\rm trans}} = {ft}$', fontsize=20, pad=15)
    cbar = fig.colorbar(c, ax=ax, pad=0.1, shrink=0.8, aspect=20)
    cbar.set_label(r'$\log\Phi~\rm [photons~cm^{-2}~s^{-1}]$', fontsize=15)
    ax.grid(True, alpha=0.3)
    ticks = ax.get_xticks()
    ax.set_xticks(ticks)
    ax.set_xticklabels([f'{int(np.degrees(t))}°' for t in ticks],
                       fontsize=13, fontweight='bold')

# --- Bottom row: VDMs ---
def style_ax(ax):
    ax.tick_params(direction='in', which='major', length=8, width=1.5, labelsize=13)
    ax.tick_params(direction='in', which='minor', length=5, width=1.0)
    ax.minorticks_on()

ax = fig.add_subplot(2, 2, 3)
im = ax.pcolormesh(lam1, tau1, psi1, shading='auto', cmap=VDM_CMAP, vmin=0, vmax=vmax)
ax.axvline(lambda0, color='red', ls='--', lw=2)
ax.set_xlabel(r'$\lambda~\rm [\AA]$', fontsize=18)
ax.set_ylabel(r'$\tau~\rm [days]$', fontsize=18)
ax.set_title(r'$f_{\rm trans} = 0.1$', fontsize=20)
plt.colorbar(im, ax=ax, label=r'$\Psi(\lambda,\tau)$', shrink=0.8, aspect=20)
style_ax(ax)

ax = fig.add_subplot(2, 2, 4)
im = ax.pcolormesh(lam2, tau2, psi2, shading='auto', cmap=VDM_CMAP, vmin=0, vmax=vmax)
ax.axvline(lambda0, color='red', ls='--', lw=2)
ax.set_xlabel(r'$\lambda~\rm [\AA]$', fontsize=18)
ax.set_ylabel(r'$\tau~\rm [days]$', fontsize=18)
ax.set_title(r'$f_{\rm trans} = 0.6$', fontsize=20)
plt.colorbar(im, ax=ax, label=r'$\Psi(\lambda,\tau)$', shrink=0.8, aspect=20)
style_ax(ax)

plt.tight_layout()
plt.savefig('plots/test_ftrans_vdm.png', dpi=150, bbox_inches='tight')
plt.close()

diff = psi2 - psi1
print(f"\nDifference stats: min={np.min(diff):.4f}, max={np.max(diff):.4f}, "
      f"mean={np.mean(np.abs(diff)):.6f}")
print("Saved plots/test_ftrans_vdm.png")
