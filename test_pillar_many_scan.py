"""
Scan: net T(r) enhancement with MANY pillars at different (r, phi).
Compare different pillar configs to find optimal parameters for boosting
temperature at large r without excessive shadow suppression.
"""
import numpy as np
import matplotlib.pyplot as plt
from pillardisk.pillar_disk import PillarDisk, load_config, resolve_rin

config = load_config('config_line.yaml')


def build_disk(hlamp=None):
    disk_params = config.get('disk', {}).copy()
    temp_params = config.get('temperature', {}).copy()
    lamp_params = config.get('lamp', {}).copy()
    obs_params = config.get('observation', {}).copy()
    if hlamp is not None:
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


def add_pillars(disk, N, h_pillar, sigma_r, sigma_phi, r_mean=10.0, sig_r_dist=5.0, seed=42):
    """Add N pillars with Gaussian radial distribution and uniform phi."""
    rng = np.random.default_rng(seed)
    r_vals = rng.normal(r_mean, sig_r_dist, N)
    r_vals = np.clip(r_vals, disk.rin + 0.5, disk.rout - 0.5)
    phi_vals = rng.uniform(0, 2*np.pi, N)
    for i in range(N):
        disk.add_pillar(r_pillar=r_vals[i], phi_pillar=phi_vals[i],
                        height=h_pillar, sigma_r=sigma_r, sigma_phi=sigma_phi,
                        modify_height=True, modify_temp=True)
    return r_vals


def get_T_ratio(disk):
    r, T_with, T_without, T_visc, _, _ = disk.get_temperature_profile()
    T_ratio = T_with / np.maximum(T_without, 1.0)
    return r, T_ratio, T_with, T_without


def style_ax(ax):
    ax.tick_params(direction='in', which='major', length=8, width=1.8, labelsize=13,
                   top=True, right=True)
    ax.tick_params(direction='in', which='minor', length=5, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()


N_pillars = 100
hlamp = 0.5

# Define configs to compare
configs = [
    {'label': r'$\rm small~steep~$' + r'$(h_p\!=\!0.2,\,\sigma_r\!=\!0.1,\,\sigma_\phi\!=\!0.1)$',
     'h_p': 0.2, 'sr': 0.1, 'sp': 0.1, 'color': 'royalblue'},
    {'label': r'$\rm small~wide~$' + r'$(h_p\!=\!0.2,\,\sigma_r\!=\!0.5,\,\sigma_\phi\!=\!0.2)$',
     'h_p': 0.2, 'sr': 0.5, 'sp': 0.2, 'color': 'orangered'},
    {'label': r'$\rm tall~steep~$' + r'$(h_p\!=\!0.5,\,\sigma_r\!=\!0.2,\,\sigma_\phi\!=\!0.1)$',
     'h_p': 0.5, 'sr': 0.2, 'sp': 0.1, 'color': 'forestgreen'},
    {'label': r'$\rm tall~wide~$' + r'$(h_p\!=\!0.5,\,\sigma_r\!=\!0.5,\,\sigma_\phi\!=\!0.3)$',
     'h_p': 0.5, 'sr': 0.5, 'sp': 0.3, 'color': 'purple'},
    {'label': r'$\rm very~tall~$' + r'$(h_p\!=\!1.0,\,\sigma_r\!=\!0.2,\,\sigma_\phi\!=\!0.1)$',
     'h_p': 1.0, 'sr': 0.2, 'sp': 0.1, 'color': 'darkorange'},
]

fig, axes = plt.subplots(2, 1, figsize=(10, 10), sharex=True)
fig.subplots_adjust(hspace=0.08)

# Top: T(r) / T_flat ratio
ax = axes[0]
for cfg in configs:
    print(f"  {cfg['label']}: h_p={cfg['h_p']}, sr={cfg['sr']}, sp={cfg['sp']}")
    disk = build_disk(hlamp=hlamp)
    add_pillars(disk, N_pillars, cfg['h_p'], cfg['sr'], cfg['sp'])
    r, T_ratio, T_with, T_without = get_T_ratio(disk)
    ax.semilogx(r, T_ratio, '-', color=cfg['color'], lw=3, alpha=0.8, label=cfg['label'])

ax.axhline(1.0, color='black', ls=':', lw=1, alpha=0.5)
ax.set_ylabel(r'$\langle T \rangle_\phi / T_{\rm bowl}$', fontsize=16)
ax.set_title(rf'$\rm {N_pillars}~pillars,~h_{{LP}} = {hlamp}~ld$', fontsize=16, pad=10)
ax.legend(fontsize=11, loc='upper left')
ax.set_xlim(1, 20)
style_ax(ax)

# Bottom: T(r) absolute
ax = axes[1]
# Baseline
disk0 = build_disk(hlamp=hlamp)
r0, _, _, T_visc0, _, _ = disk0.get_temperature_profile()
_, T0, _, _ = get_T_ratio(disk0)
ax.loglog(r0, disk0.get_temperature_profile()[1], '--', color='black', lw=2, alpha=0.5,
          label=r'$\rm bowl~disk$')

for cfg in configs:
    disk = build_disk(hlamp=hlamp)
    add_pillars(disk, N_pillars, cfg['h_p'], cfg['sr'], cfg['sp'])
    r, _, T_with, _ = get_T_ratio(disk)
    ax.loglog(r, T_with, '-', color=cfg['color'], lw=3, alpha=0.8, label=cfg['label'])

ax.set_xlabel(r'$r~\rm [ld]$', fontsize=16)
ax.set_ylabel(r'$\langle T \rangle_\phi~\rm [K]$', fontsize=16)
ax.set_xlim(1, 20)
style_ax(ax)

plt.savefig('plots/test_pillar_many_scan.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved plots/test_pillar_many_scan.png")

# --- Summary ---
print(f"\n--- Net effect with {N_pillars} pillars, h_LP={hlamp} ---")
print(f"{'config':>20s} {'mean_boost(5-15ld)':>18s} {'max_boost':>10s} {'min_ratio':>10s}")
for cfg in configs:
    disk = build_disk(hlamp=hlamp)
    add_pillars(disk, N_pillars, cfg['h_p'], cfg['sr'], cfg['sp'])
    r, T_ratio, _, _ = get_T_ratio(disk)
    mask = (r > 5) & (r < 15)
    mean_boost = np.mean(T_ratio[mask])
    print(f"{cfg['label']:>20s} {mean_boost:18.4f} {np.max(T_ratio):10.4f} {np.min(T_ratio[r>2]):10.4f}")
