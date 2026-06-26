"""Plot azimuthally-averaged T(r) with pillar irradiation boost included."""
import numpy as np
import matplotlib.pyplot as plt
from pillardisk.pillar_disk import PillarDisk, load_config, resolve_rin
from pillardisk.pillar_line import parse_math_expr

config = load_config('config_line.yaml')


def build_disk():
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
    if 'inclination' in obs_params:
        obs_params['cosi'] = np.cos(np.radians(obs_params['inclination']))
        del obs_params['inclination']
    return PillarDisk(**{**disk_params, **temp_params, **lamp_params, **obs_params})


def add_random_pillars(disk, n, seed=42):
    rng = np.random.default_rng(seed)
    r_vals = rng.uniform(disk.rin + 1, disk.rout - 1, n)
    phi_vals = rng.uniform(0, 2 * np.pi, n)
    tf_rand = float(config.get('pillars', {}).get('temp_factor_pillar', 1))
    print(f"  Random pillars: temp_factor={tf_rand}")
    for i in range(n):
        disk.add_pillar(r_pillar=r_vals[i], phi_pillar=phi_vals[i],
                        height=0.15, sigma_r=0.5, sigma_phi=0.2,
                        temp_factor=tf_rand,
                        modify_height=True, modify_temp=True)


def derive_T_irrad(T_total, T_visc):
    return np.maximum(T_total**4 - T_visc**4, 0)**0.25


def style_ax(ax):
    ax.tick_params(direction='in', which='major', length=8, width=1.8, labelsize=13,
                   top=True, right=True)
    ax.tick_params(direction='in', which='minor', length=5, width=1.0,
                   top=True, right=True)


# --- No pillars ---
disk0 = build_disk()
r0, T_tot0, _, T_visc0, _, _ = disk0.get_temperature_profile()
T_irrad0 = derive_T_irrad(T_tot0, T_visc0)

# --- 2 pillars (read params from config) ---
pillars_cfg = config.get('pillars', {})
tf_list = pillars_cfg.get('temp_factor', [1])
if not isinstance(tf_list, list):
    tf_list = [tf_list]

disk2 = build_disk()
for i, (r_p, phi_p) in enumerate(zip([4.0, 4.0], [parse_math_expr('3*pi/4'), parse_math_expr('0')])):
    tf = float(tf_list[i] if i < len(tf_list) else tf_list[-1])
    disk2.add_pillar(r_pillar=r_p, phi_pillar=phi_p, height=0.15,
                     sigma_r=0.5, sigma_phi=0.2,
                     temp_factor=tf, modify_height=True, modify_temp=True)
    print(f"  Pillar {i+1}: r={r_p}, phi={phi_p:.2f}, temp_factor={tf}")
r2, T_tot2, _, T_visc2, _, _ = disk2.get_temperature_profile()
T_irrad2 = derive_T_irrad(T_tot2, T_visc2)

# --- 100 random pillars ---
disk100 = build_disk()
add_random_pillars(disk100, 100)
r100, T_tot100, _, T_visc100, _, _ = disk100.get_temperature_profile()
T_irrad100 = derive_T_irrad(T_tot100, T_visc100)

r_isco = disk0.rin
xlim = (r_isco * 0.3, r0[-1])


def add_isco_shade(ax):
    ylim = ax.get_ylim()
    ax.axvspan(xlim[0], r_isco, color='lightgray', alpha=0.5, zorder=0)
    ax.axvline(r_isco, color='gray', ls='-', lw=1.5, alpha=0.7)
    ax.text(r_isco * 0.65, 10**(0.5 * (np.log10(ylim[0]) + np.log10(ylim[1]))),
            r'$r < r_{\rm ISCO}$', fontsize=12, color='gray', rotation=90,
            ha='center', va='center')
    ax.set_ylim(ylim)


# --- Plot ---
fig, ax = plt.subplots(1, 1, figsize=(10, 7))

ax.loglog(r0, T_visc0, '--', color='gray', lw=3, alpha=0.7, label=r'$T_{\rm visc}$')
ax.loglog(r0, T_tot0, '-', color='black', lw=3, alpha=0.8, label=r'$\rm no~pillars$')
ax.loglog(r2, T_tot2, '-', color='orangered', lw=3, alpha=0.8, label=r'$\rm 2~pillars$')
ax.loglog(r100, T_tot100, '-', color='royalblue', lw=3, alpha=0.8, label=r'$\rm 100~pillars$')
ax.set_xlabel(r'$r~\rm [ld]$', fontsize=16)
ax.set_ylabel(r'$\langle T \rangle_\phi~\rm [K]$', fontsize=16)
ax.legend(fontsize=13)
ax.set_xlim(xlim)
style_ax(ax)
add_isco_shade(ax)

plt.savefig('plots/test_pillar_boost.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved plots/test_pillar_boost.png")
