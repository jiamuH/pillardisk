"""
Continuum lag spectrum tau(lambda) and T(r) for parameter scan around best config.
Saves data to .npz so plots can be tweaked without rerunning.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from pillardisk.pillar_disk import PillarDisk, load_config, resolve_rin

config = load_config('config_line.yaml')

wavelengths = np.logspace(np.log10(1000), np.log10(10000), 60)


BETA = 50.0  # bowl steepness — high beta keeps outer edge out of shadow

# Ripple-disc reference (axisymmetric Starkey+23-style structure).
# Because our bowl uses a steep beta (h_bowl ~ 0 except near r_out), we
# add the ripple ADDITIVELY on top of the bowl, like the pillars, with
# amplitude in light-days. Set h_ripple = 0 to disable.
# Parameters otherwise follow bowl.for / Starkey+23 conventions.
RIPPLE_PARAMS = dict(
    h_ripple=0.5,  # ripple amplitude in light days (analogous to h_p)
    wpow=0.0,      # amplitude scaling: amp(r) = h_ripple * (r/r0)^wpow
    wnum=2.0,      # ripple wave number k1 at r=r0 (rad / ld); wavelength = 2pi/wnum
    wnpow=0.0,     # k(r) = k1 * (r/r0)^wnpow
    crest=2.0,     # crest peakiness; ripple_norm = ((sin+1)/2)^crest
    r_in_ripple=2.0,  # ripples vanish for r < r_in_ripple (matches the
                      # pillar inner clip in add_pillars())
)
RIPPLE_CACHE = 'plots/lag_spectrum_ripple.npz'
# Bump this when the ripple formula changes so old caches are invalidated
# even if RIPPLE_PARAMS values are unchanged.
RIPPLE_FORMULA_VERSION = 2


def make_ripple_height(r, h1, r0, beta,
                       h_ripple=0.5, wpow=0.0, wnum=2.0, wnpow=0.0,
                       crest=2.0, r_in_ripple=2.0):
    """Bowl + additive axisymmetric ripple.

    h(r) = h_bowl(r) + h_ripple * (r/r0)^wpow * ripple_norm(r)  for r >= r_in_ripple
    h(r) = h_bowl(r)                                            for r < r_in_ripple

    where h_bowl(r) = h1 (r/r0)^beta and the dimensionless ripple
    profile is normalised to lie in [0, 1] so h_ripple is the peak
    bump height in light days, matching the convention used for
    pillars (where h_p is also the peak above the bowl).

    To avoid a discontinuity at r = r_in_ripple, the phase is chosen so
    that the ripple is at its trough (sin = -1, ripple_norm = 0) exactly
    at r_in_ripple:
        phase(r) = wnum * r0 * [(r/r0)^(wnpow+1) - (r_in/r0)^(wnpow+1)]
                   / (wnpow+1)  -  pi/2 .
    Then ripple_norm(r) = ((sin(phase)+1)/2)^crest, which is the
    bowl.for transformation but normalised to peak = 1. crest = 1 gives
    a smooth (1+sin)/2 oscillation; crest > 1 sharpens the crests and
    suppresses the troughs.

    The r_in_ripple cut matches the inner radius at which pillars start
    in add_pillars() so the rippled and pillared discs share the same
    structured belt.
    """
    h_bowl = h1 * (r / r0) ** beta
    if h_ripple == 0.0:
        return h_bowl
    x = r / r0
    x_in = r_in_ripple / r0
    amp = h_ripple * x ** wpow
    pow_ = wnpow + 1.0
    # Phase 0 anchored at r_in_ripple, then offset by -pi/2 so the
    # ripple is at trough (smooth zero) at that radius.
    phase = wnum * r0 * (x ** pow_ - x_in ** pow_) / pow_ - 0.5 * np.pi
    ripple = ((np.sin(phase) + 1.0) / 2.0) ** crest
    delta_h = amp * ripple
    delta_h = np.where(r < r_in_ripple, 0.0, delta_h)
    return h_bowl + delta_h


def build_ripple_disk(hlamp=0.5, ripple_params=None):
    """Build a bowl disc and overwrite h_base with the rippled profile."""
    rp = ripple_params or RIPPLE_PARAMS
    disk = build_disk(hlamp=hlamp)
    disk.h_base = make_ripple_height(
        disk.r, disk.h1, disk.r0, disk.beta,
        h_ripple=rp['h_ripple'], wpow=rp['wpow'],
        r_in_ripple=rp['r_in_ripple'],
        wnum=rp['wnum'], wnpow=rp['wnpow'], crest=rp['crest'],
    )
    return disk


def compute_ripple_tau_and_T(hlamp, ripple_params=None):
    """Compute tau(lambda) and T(r) for the ripple-disc reference.

    Cached separately from the bowl/pillar data so existing caches don't
    need to be invalidated when ripple parameters change.
    """
    disk = build_ripple_disk(hlamp=hlamp, ripple_params=ripple_params)
    tau_mean, tau_grid, _ = disk.compute_lag_spectrum(
        wavelengths, ntau=500, taumax=50.0, parallel=True)
    r, T_with, _, _, _, _ = disk.get_temperature_profile()
    return tau_mean, tau_grid, r, T_with

def build_disk(hlamp=None):
    disk_params = config.get('disk', {}).copy()
    temp_params = config.get('temperature', {}).copy()
    lamp_params = config.get('lamp', {}).copy()
    obs_params = config.get('observation', {}).copy()
    temp_params['beta'] = BETA
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


def add_pillars(disk, N, h_pillar, sigma_r, sigma_phi, r_mean=15.0, sig_r_dist=5.0,
                seed=42, r_min=2.0):
    rng = np.random.default_rng(seed)
    r_vals = rng.normal(r_mean, sig_r_dist, N)
    r_vals = np.clip(r_vals, r_min, disk.rout)
    phi_vals = rng.uniform(0, 2*np.pi, N)
    for i in range(N):
        disk.add_pillar(r_pillar=r_vals[i], phi_pillar=phi_vals[i],
                        height=h_pillar, sigma_r=sigma_r, sigma_phi=sigma_phi,
                        temp_factor=1.0, modify_height=True, modify_temp=False)


def compute_tau_and_T(hlamp, h_p, sr, sp):
    """Compute tau(lambda), psi(lambda, tau), and T(r) for a given config."""
    disk = build_disk(hlamp=hlamp)
    if h_p > 0:
        add_pillars(disk, 100, h_p, sr, sp)
    tau_mean, tau_grid, psi = disk.compute_lag_spectrum(wavelengths, ntau=500, taumax=50.0,
                                                        parallel=True)
    r, T_with, T_without, T_visc, _, _ = disk.get_temperature_profile()
    return tau_mean, tau_grid, psi, r, T_with, T_without


def style_ax(ax):
    ax.tick_params(direction='in', which='major', length=8, width=1.8, labelsize=13,
                   top=True, right=True)
    ax.tick_params(direction='in', which='minor', length=5, width=1.0,
                   top=True, right=True)
    ax.minorticks_on()


def scan_colors(n, cmap='plasma'):
    """Generate n colors from a sequential colormap."""
    return [plt.colormaps[cmap](x) for x in np.linspace(0.15, 0.85, n)]


# Scan definitions (same as scan2)
scans = {
    'vary_hp': {
        'title': r'$\mathrm{Vary}~h_p~(\sigma_r=0.5,~\sigma_\phi=0.3)$',
        'param_name': 'h_p', 'values': [0.1, 0.3, 0.5, 1.0, 2.0],
        'fixed': {'sr': 0.5, 'sp': 0.3, 'hlamp': 0.5},
        'labels': [rf'$h_p = {v}$' for v in [0.1, 0.3, 0.5, 1.0, 2.0]],
    },
    'vary_sr': {
        'title': r'$\mathrm{Vary}~\sigma_r~(h_p=0.5,~\sigma_\phi=0.3)$',
        'param_name': 'sr', 'values': [0.1, 0.3, 0.5, 1.0, 2.0],
        'fixed': {'h_p': 0.5, 'sp': 0.3, 'hlamp': 0.5},
        'labels': [rf'$\sigma_r = {v}$' for v in [0.1, 0.3, 0.5, 1.0, 2.0]],
    },
    'vary_sp': {
        'title': r'$\mathrm{Vary}~\sigma_\phi~(h_p=0.5,~\sigma_r=0.5)$',
        'param_name': 'sp', 'values': [0.05, 0.1, 0.3, 0.5, 1.0],
        'fixed': {'h_p': 0.5, 'sr': 0.5, 'hlamp': 0.5},
        'labels': [rf'$\sigma_\phi = {v}$' for v in [0.05, 0.1, 0.3, 0.5, 1.0]],
    },
    'vary_hlamp': {
        'title': r'$\mathrm{Vary}~h_{\rm LP}~(h_p=0.5,~\sigma_r=0.5,~\sigma_\phi=0.3)$',
        'param_name': 'hlamp', 'values': [0.1, 0.3, 0.5, 1.0, 2.0],
        'fixed': {'h_p': 0.5, 'sr': 0.5, 'sp': 0.3},
        'labels': [rf'$h_{{\rm LP}} = {v}$' for v in [0.1, 0.3, 0.5, 1.0, 2.0]],
    },
}

DATAFILE = 'plots/lag_spectrum_scan_data.npz'


def run_scan():
    """Compute all tau(lambda) and T(r) and save to npz."""
    all_data = {'wavelengths': wavelengths}

    # Bowl baseline (hlamp=0.5)
    print("Computing bowl baseline (hlamp=0.5)...")
    tau_bowl, tau_grid, psi_bowl, r_grid, T_bowl, _ = compute_tau_and_T(0.5, 0, 0, 0)
    all_data['tau_bowl_0.5'] = tau_bowl
    all_data['tau_grid'] = tau_grid
    all_data['psi_bowl_0.5'] = psi_bowl
    all_data['r_grid'] = r_grid
    all_data['T_bowl_0.5'] = T_bowl

    # Bowl baselines for other hlamp values
    for hl in [0.1, 0.3, 1.0, 2.0]:
        print(f"Computing bowl baseline (hlamp={hl})...")
        tau_b, _, psi_b, _, T_b, _ = compute_tau_and_T(hl, 0, 0, 0)
        all_data[f'tau_bowl_{hl}'] = tau_b
        all_data[f'psi_bowl_{hl}'] = psi_b
        all_data[f'T_bowl_{hl}'] = T_b

    for scan_name, scan in scans.items():
        for i, val in enumerate(scan['values']):
            params = dict(scan['fixed'])
            params[scan['param_name']] = val
            key = f"{scan_name}_{i}"
            print(f"  {scan_name}: {scan['param_name']}={val}")
            tau, _, psi, _, T_with, _ = compute_tau_and_T(
                params['hlamp'], params['h_p'], params['sr'], params['sp'])
            all_data[f'tau_{key}'] = tau
            all_data[f'psi_{key}'] = psi
            all_data[f'T_{key}'] = T_with

    np.savez(DATAFILE, **all_data)
    print(f"Saved data to {DATAFILE}")
    return all_data


def load_data():
    data = dict(np.load(DATAFILE, allow_pickle=True))
    return data


def load_or_compute_ripple(hlamp=0.5, ripple_params=None):
    """Load cached ripple-disc tau/T or compute and cache them.

    Cached separately from the bowl/pillar scan so the heavy scan cache
    doesn't have to be invalidated when ripple parameters change.
    """
    rp = ripple_params or RIPPLE_PARAMS
    if os.path.exists(RIPPLE_CACHE):
        d = dict(np.load(RIPPLE_CACHE, allow_pickle=True))
        # invalidate cache if ripple parameters or the formula version changed
        same_params = all(float(d.get(k, np.nan)) == float(rp[k]) for k in rp)
        same_hlamp = float(d.get('hlamp', np.nan)) == float(hlamp)
        same_version = int(d.get('formula_version', -1)) == RIPPLE_FORMULA_VERSION
        if same_params and same_hlamp and same_version:
            return d
    print(f"Computing ripple-disc reference (hlamp={hlamp}, params={rp})...")
    tau, tau_grid, r, T = compute_ripple_tau_and_T(hlamp=hlamp, ripple_params=rp)
    out = {'tau_ripple': tau, 'tau_grid': tau_grid, 'r_grid': r,
           'T_ripple': T, 'hlamp': hlamp,
           'formula_version': RIPPLE_FORMULA_VERSION,
           **{k: float(v) for k, v in rp.items()}}
    np.savez(RIPPLE_CACHE, **out)
    print(f"Saved ripple cache to {RIPPLE_CACHE}")
    return out


def plot_all(data):
    wavelengths = data['wavelengths']
    r_grid = data['r_grid']
    ripple = load_or_compute_ripple(hlamp=0.5)

    # --- Plot 1a: absolute tau(lambda) ---
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    for idx, (scan_name, scan) in enumerate(scans.items()):
        ax = axes[idx // 2, idx % 2]
        colors = scan_colors(len(scan['values']))
        for i, val in enumerate(scan['values']):
            key = f"{scan_name}_{i}"
            tau = data[f'tau_{key}']
            ax.plot(wavelengths, tau, lw=3.5, alpha=0.7, color=colors[i], label=scan['labels'][i])
        hlamp_val = scan['fixed'].get('hlamp', 0.5)
        tau_bowl = data[f'tau_bowl_{hlamp_val}']
        ax.plot(wavelengths, tau_bowl, '--', color='dodgerblue', lw=4, alpha=1.0,
                label=r'$\rm bowl$')
        ax.set_title(scan['title'], fontsize=14, pad=10)
        ax.set_ylabel(r'$\tau~\rm [days]$', fontsize=16)
        if idx >= 2:
            ax.set_xlabel(r'$\lambda~\rm [\AA]$', fontsize=16)
        ax.legend(fontsize=11)
        ax.set_xlim(wavelengths[0], wavelengths[-1])
        ax.set_xscale('log')
        style_ax(ax)
    plt.tight_layout()
    plt.savefig('plots/test_lag_spectrum_tau.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved plots/test_lag_spectrum_tau.png")

    # --- Plot 1b: tau(lambda) - tau(lambda_ref) ---
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    for idx, (scan_name, scan) in enumerate(scans.items()):
        ax = axes[idx // 2, idx % 2]
        colors = scan_colors(len(scan['values']))
        for i, val in enumerate(scan['values']):
            key = f"{scan_name}_{i}"
            tau = data[f'tau_{key}']
            i_ref = np.argmin(np.abs(wavelengths - 1500.0))
            ax.plot(wavelengths, tau - tau[i_ref], lw=3.5, alpha=0.7, color=colors[i], label=scan['labels'][i])
        hlamp_val = scan['fixed'].get('hlamp', 0.5)
        tau_bowl = data[f'tau_bowl_{hlamp_val}']
        i_ref = np.argmin(np.abs(wavelengths - 1500.0))
        ax.plot(wavelengths, tau_bowl - tau_bowl[i_ref], '--', color='dodgerblue', lw=4, alpha=1.0,
                label=(r'$\rm bowl$' if idx == 0 else None))
        # Rippled-disc reference (axisymmetric Starkey+23-style)
        tau_rip = ripple['tau_ripple']
        ax.plot(wavelengths, tau_rip - tau_rip[i_ref], '-', color='forestgreen', lw=3.5,
                alpha=0.5, zorder=-10,
                label=(r'$\rm ripple$' if idx == 0 else None))
        ax.set_title(scan['title'], fontsize=14, pad=10)
        ax.text(0.03, 0.95, rf'$\rm ({chr(ord("a") + idx)})$',
                transform=ax.transAxes, fontsize=26, ha='left', va='top',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                          alpha=0.85, edgecolor='black', linewidth=1.0))
        if idx % 2 == 0:
            ax.set_ylabel(r'$\tau(\lambda) - \tau(1500\,\rm \AA)~\rm [days]$', fontsize=16)
        else:
            ax.tick_params(labelleft=False)
        if idx >= 2:
            ax.set_xlabel(r'$\lambda~\rm [\AA]$', fontsize=16)
        ax.legend(fontsize=15)
        ax.set_xlim(wavelengths[0], wavelengths[-1])
        ax.set_xscale('log')
        style_ax(ax)
    plt.tight_layout()
    plt.savefig('plots/test_lag_spectrum_tau_rel.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved plots/test_lag_spectrum_tau_rel.png")

    # --- Plot 2: T(r) ---
    # Get T_visc from a bare disk (no irradiation dependence on pillars)
    disk_tmp = build_disk(hlamp=0.5)
    T_visc = np.interp(r_grid, disk_tmp.r, disk_tmp.tv_base)
    r_isco_ld = disk_tmp.rin

    def _plot_T(xscale, savepath, xlim):
        fig, axes = plt.subplots(2, 2, figsize=(14, 12))
        for idx, (scan_name, scan) in enumerate(scans.items()):
            ax = axes[idx // 2, idx % 2]
            colors = scan_colors(len(scan['values']))
            for i, val in enumerate(scan['values']):
                key = f"{scan_name}_{i}"
                T = data[f'T_{key}']
                ax.plot(r_grid, T, lw=3.5, alpha=0.7,
                        color=colors[i], label=scan['labels'][i])
            hlamp_val = scan['fixed'].get('hlamp', 0.5)
            T_bowl = data[f'T_bowl_{hlamp_val}']
            ax.plot(r_grid, T_bowl, '--', color='dodgerblue', lw=4, alpha=1.0,
                    label=(r'$\rm bowl$' if idx == 0 else None))
            # Rippled-disc reference (axisymmetric Starkey+23-style)
            r_rip = ripple['r_grid']
            T_rip = ripple['T_ripple']
            ax.plot(r_rip, T_rip, '-', color='forestgreen', lw=3.5, alpha=0.5,
                    zorder=-10,
                    label=(r'$\rm ripple$' if idx == 0 else None))
            ax.plot(r_grid, T_visc, '-', color='gray', lw=3.5, alpha=0.5,
                    label=(r'$T_{\rm visc}$' if idx == 0 else None))
            # Reference power laws normalized to total bowl T at r = 1 ld
            i_ref = np.argmin(np.abs(r_grid - 1.0))
            T_ref = T_bowl[i_ref]
            r_ref = r_grid[i_ref]
            T_p34 = T_ref * (r_grid / r_ref) ** (-0.75)
            T_p78 = T_ref * (r_grid / r_ref) ** (-0.875)
            ax.plot(r_grid, T_p34, ':', color='black', lw=4, alpha=0.7,
                    label=(r'$T \propto r^{-3/4}$' if idx == 0 else None))
            ax.plot(r_grid, T_p78, ':', color='saddlebrown', lw=4, alpha=0.7,
                    label=(r'$T \propto r^{-7/8}$' if idx == 0 else None))
            # Shade the unphysical region inside r_ISCO
            ax.axvspan(0, r_isco_ld, color='gray', alpha=0.25, zorder=0)
            ax.set_xscale(xscale)
            ax.set_yscale('log')
            ax.set_title(scan['title'], fontsize=14, pad=10)
            ax.text(0.03, 0.95, rf'$\rm ({chr(ord("a") + idx)})$',
                    transform=ax.transAxes, fontsize=26, ha='left', va='top',
                    bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                              alpha=0.85, edgecolor='black', linewidth=1.0))
            if idx % 2 == 0:
                ax.set_ylabel(r'$\langle T \rangle_\phi~\rm [K]$', fontsize=16)
            else:
                ax.tick_params(labelleft=False)
            if idx >= 2:
                ax.set_xlabel(r'$r~\rm [ld]$', fontsize=16)
            ax.legend(fontsize=14, loc='lower left')
            ax.set_xlim(*xlim)
            ax.set_ylim(1e2, 1e5)
            style_ax(ax)
        plt.tight_layout()
        plt.savefig(savepath, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Saved {savepath}")

    # Pillar belt is clipped to r_min=2 ld in add_pillars; zoom to that region
    _plot_T('log', 'plots/test_lag_spectrum_T.png', xlim=(1, 22))

    # --- Plot 3: 2D response functions per scan ---
    tau_grid = data['tau_grid']
    for scan_name, scan in scans.items():
        n_vals = len(scan['values'])
        from matplotlib.colors import PowerNorm
        from matplotlib.gridspec import GridSpec

        ncols = n_vals + 1
        fig = plt.figure(figsize=(4 * ncols, 14))
        gs = GridSpec(3, ncols, figure=fig, height_ratios=[1.2, 1, 1],
                      hspace=0.1, wspace=0.3)

        hlamp_val = scan['fixed'].get('hlamp', 0.5)

        # --- Row 0: face-on polar T(r,phi) maps ---
        r_plot = np.linspace(0.5, 20, 150)
        phi_plot = np.linspace(0, 2*np.pi, 360)
        r_2d_p, phi_2d_p = np.meshgrid(r_plot, phi_plot, indexing='ij')

        # Compute all T maps and find shared color scale
        T_maps = []
        # Bowl
        disk_bowl = build_disk(hlamp=hlamp_val)
        T_bowl_2d = disk_bowl.get_temperature(r_2d_p, phi_2d_p)
        T_maps.append(T_bowl_2d)
        for i, val in enumerate(scan['values']):
            params = dict(scan['fixed'])
            params[scan['param_name']] = val
            disk_i = build_disk(hlamp=params['hlamp'])
            add_pillars(disk_i, 100, params['h_p'], params['sr'], params['sp'])
            T_i = disk_i.get_temperature(r_2d_p, phi_2d_p)
            T_maps.append(T_i)

        from matplotlib.colors import LogNorm, PowerNorm as PN2
        T_all = np.concatenate([t.ravel() for t in T_maps])
        T_vmin = max(np.percentile(T_all[T_all > 100], 1), 500)
        T_vmax = np.percentile(T_all, 99)
        T_norm = PN2(gamma=0.5, vmin=T_vmin, vmax=T_vmax)

        # Plot bowl
        ax = fig.add_subplot(gs[0, 0], projection='polar')
        c = ax.pcolormesh(phi_2d_p, r_2d_p, np.clip(T_maps[0], T_vmin, None),
                          shading='auto', cmap='hot', norm=T_norm)
        ax.set_title(r'$\rm bowl$', fontsize=30, pad=10)
        ax.grid(True, alpha=0.7, linewidth=1.2, color='white')
        ticks = ax.get_xticks()
        ax.set_xticks(ticks)
        ax.set_xticklabels(
            [rf'${int(np.degrees(t))}^{{\circ}}$' for t in ticks],
            fontsize=16)
        ax.set_yticks([10, 20])
        ytl = ax.set_yticklabels([r'$10~\rm ld$', r'$20~\rm ld$'], fontsize=18)
        ytl[0].set_color('white')
        ytl[1].set_color('black')

        # Plot each config — suppress angle/radius labels (shown only on bowl)
        for i in range(n_vals):
            ax = fig.add_subplot(gs[0, i + 1], projection='polar')
            c = ax.pcolormesh(phi_2d_p, r_2d_p, np.clip(T_maps[i + 1], T_vmin, None),
                              shading='auto', cmap='hot', norm=T_norm)
            ax.set_title(scan['labels'][i], fontsize=28, pad=10)
            ax.grid(True, alpha=0.7, linewidth=1.2, color='white')
            ax.set_xticklabels([])
            ax.set_yticks([10, 20])
            ax.set_yticklabels([])

        # Shared colorbar for T maps (units of 10^4 K)
        cbar_ax = fig.add_axes([0.92, 0.68, 0.018, 0.2])
        cbar = fig.colorbar(c, cax=cbar_ax)
        cbar.set_label(r'$T~\rm [10^{4}~K]$', fontsize=30)
        cbar.set_ticks([1e3, 5e3, 1e4, 2e4, 4e4])
        cbar.set_ticklabels([r'$0.1$', r'$0.5$', r'$1$', r'$2$', r'$4$'])
        cbar.ax.tick_params(labelsize=26)

        # --- Rows 1-2: psi and difference (use gridspec axes, not subplot) ---
        # Build regular axes for rows 1 and 2
        axes_psi = [fig.add_subplot(gs[1, j]) for j in range(ncols)]
        axes_diff = [fig.add_subplot(gs[2, j]) for j in range(ncols)]

        # Bowl baseline
        psi_bowl = data[f'psi_bowl_{hlamp_val}']
        tau_bowl_scan = data[f'tau_bowl_{hlamp_val}']

        # Collect all psi for consistent color scale
        all_psi = [psi_bowl]
        all_diff = []
        for i in range(len(scan['values'])):
            key = f"{scan_name}_{i}"
            all_psi.append(data[f'psi_{key}'])
            all_diff.append(data[f'psi_{key}'] - psi_bowl)

        # Log normalization for psi panels (much better contrast than PowerNorm)
        tau_mask = tau_grid < 40
        psi_pos_all = np.concatenate([p[tau_mask, :][p[tau_mask, :] > 0].ravel() for p in all_psi])
        from matplotlib.colors import LogNorm as LN
        psi_vmin = np.percentile(psi_pos_all, 5)
        psi_vmax = np.percentile(psi_pos_all, 99.5)
        norm_pow = LN(vmin=psi_vmin, vmax=psi_vmax)

        # Absolute difference with power-law normalization for visibility
        all_diff_flat = np.concatenate([d.ravel() for d in all_diff])
        dmax = np.percentile(np.abs(all_diff_flat[all_diff_flat != 0]), 95) if np.any(all_diff_flat != 0) else 1

        TAU_YLIM = 40

        # Row 1: psi(lambda, tau)
        ax = axes_psi[0]
        im = ax.pcolormesh(wavelengths, tau_grid, np.maximum(psi_bowl, psi_vmin),
                           shading='auto', cmap='inferno_r', norm=norm_pow)
        ax.plot(wavelengths, tau_bowl_scan, 'k--', lw=1.5, alpha=0.8)
        ax.set_ylabel(r'$\tau~\rm [days]$', fontsize=30)
        ax.set_xscale('log')
        ax.set_xticklabels([])
        ax.set_ylim(0, TAU_YLIM)
        style_ax(ax)
        ax.tick_params(which='major', labelsize=27)

        # Row 2: blank for bowl column
        axes_diff[0].axis('off')

        for i, val in enumerate(scan['values']):
            key = f"{scan_name}_{i}"
            psi = data[f'psi_{key}']
            tau_mean = data[f'tau_{key}']

            ax = axes_psi[i + 1]
            im = ax.pcolormesh(wavelengths, tau_grid, np.maximum(psi, psi_vmin),
                               shading='auto', cmap='inferno_r', norm=norm_pow)
            ax.plot(wavelengths, tau_mean, 'k--', lw=1.5, alpha=0.8)
            ax.set_xscale('log')
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_ylim(0, TAU_YLIM)
            style_ax(ax)
            ax.tick_params(which='major', labelsize=27)

            ax = axes_diff[i + 1]
            im_diff_last = ax.pcolormesh(
                wavelengths, tau_grid, all_diff[i], shading='auto',
                cmap='RdBu_r', vmin=-dmax, vmax=dmax)
            ax.set_xscale('log')
            ax.set_xlabel(r'$\lambda~\rm [\AA]$', fontsize=30)
            ax.set_yticklabels([])
            ax.set_ylim(0, TAU_YLIM)
            style_ax(ax)
            ax.tick_params(which='major', labelsize=27)

        # Shared colorbars for the psi (row 1) and diff (row 2) panels.
        # Position dynamically against the rightmost axis of each row so the
        # cbar height/y-range matches the panels exactly.
        fig.canvas.draw()  # ensure axis positions are finalised
        pos_psi = axes_psi[-1].get_position()
        cax_psi = fig.add_axes([pos_psi.x1 + 0.012, pos_psi.y0,
                                0.018, pos_psi.height])
        cb_psi = fig.colorbar(im, cax=cax_psi)
        cb_psi.set_label(r'$\Psi(\lambda,\tau)$', fontsize=30)
        cb_psi.ax.tick_params(labelsize=26)
        pos_diff = axes_diff[-1].get_position()
        cax_diff = fig.add_axes([pos_diff.x1 + 0.012, pos_diff.y0,
                                 0.018, pos_diff.height])
        cb_diff = fig.colorbar(im_diff_last, cax=cax_diff)
        cb_diff.set_label(r'$\Delta\Psi$', fontsize=30)
        cb_diff.ax.tick_params(labelsize=26)
        # Don't call tight_layout since we placed colorbars with add_axes.
        # Save inverted version (shadows = white)
        fname_inv = f'plots/test_lag_spectrum_psi_{scan_name}_inv.png'
        plt.savefig(fname_inv, dpi=200, bbox_inches='tight')
        plt.close()
        print(f"Saved {fname_inv}")

        # --- Re-plot with original (non-inverted) colormap ---
        fig2 = plt.figure(figsize=(4 * ncols, 14))
        gs2 = GridSpec(3, ncols, figure=fig2, height_ratios=[1.2, 1, 1],
                       hspace=0.25, wspace=0.15)

        # Row 0: T maps (hot_r, non-inverted version)
        T_norm2 = PN2(gamma=0.5, vmin=T_vmin, vmax=T_vmax)
        ax = fig2.add_subplot(gs2[0, 0], projection='polar')
        c = ax.pcolormesh(phi_2d_p, r_2d_p, np.clip(T_maps[0], T_vmin, None),
                          shading='auto', cmap='hot_r', norm=T_norm2)
        ax.set_title(r'$\rm bowl$', fontsize=13, pad=10)
        ax.grid(True, alpha=0.7, linewidth=1.2, color='white')
        ticks = ax.get_xticks()
        ax.set_xticks(ticks)
        ax.set_xticklabels([f'{int(np.degrees(t))}°' for t in ticks],
                           fontsize=10, fontweight='bold')
        ax.set_yticks([10, 20])
        ax.set_yticklabels([r'$10~\rm ld$', r'$20~\rm ld$'], fontsize=9)
        for ii in range(n_vals):
            ax = fig2.add_subplot(gs2[0, ii + 1], projection='polar')
            c = ax.pcolormesh(phi_2d_p, r_2d_p, np.clip(T_maps[ii + 1], T_vmin, None),
                              shading='auto', cmap='hot_r', norm=T_norm2)
            ax.set_title(scan['labels'][ii], fontsize=11, pad=10)
            ax.grid(True, alpha=0.7, linewidth=1.2, color='white')
            ticks = ax.get_xticks()
            ax.set_xticks(ticks)
            ax.set_xticklabels([f'{int(np.degrees(t))}°' for t in ticks],
                               fontsize=10, fontweight='bold')
            ax.set_yticks([10, 20])
            ax.set_yticklabels([r'$10~\rm ld$', r'$20~\rm ld$'], fontsize=9)
        cbar_ax2 = fig2.add_axes([0.92, 0.68, 0.01, 0.2])
        fig2.colorbar(c, cax=cbar_ax2).set_label(r'$T~\rm [K]$', fontsize=13)

        # Row 1: psi (inferno, not inverted)
        axes_psi2 = [fig2.add_subplot(gs2[1, j]) for j in range(ncols)]
        axes_diff2 = [fig2.add_subplot(gs2[2, j]) for j in range(ncols)]

        ax = axes_psi2[0]
        im_psi2 = ax.pcolormesh(wavelengths, tau_grid,
                                 np.maximum(psi_bowl, psi_vmin),
                                 shading='auto', cmap='inferno',
                                 norm=norm_pow)
        ax.plot(wavelengths, tau_bowl_scan, 'w--', lw=1.5, alpha=0.8)
        ax.set_ylabel(r'$\tau~\rm [days]$', fontsize=14)
        ax.set_xscale('log')
        ax.set_xticklabels([])
        ax.set_ylim(0, TAU_YLIM)
        style_ax(ax)
        axes_diff2[0].axis('off')

        im_diff2_last = None
        for ii in range(n_vals):
            key = f"{scan_name}_{ii}"
            psi = data[f'psi_{key}']
            tau_m = data[f'tau_{key}']
            ax = axes_psi2[ii + 1]
            ax.pcolormesh(wavelengths, tau_grid, np.maximum(psi, psi_vmin),
                          shading='auto', cmap='inferno', norm=norm_pow)
            ax.plot(wavelengths, tau_m, 'w--', lw=1.5, alpha=0.8)
            ax.set_xscale('log')
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_ylim(0, TAU_YLIM)
            style_ax(ax)

            ax = axes_diff2[ii + 1]
            im_diff2_last = ax.pcolormesh(
                wavelengths, tau_grid, all_diff[ii], shading='auto',
                cmap='RdBu_r', vmin=-dmax, vmax=dmax)
            ax.set_xscale('log')
            ax.set_xlabel(r'$\lambda~\rm [\AA]$', fontsize=14)
            ax.set_yticklabels([])
            ax.set_ylim(0, TAU_YLIM)
            style_ax(ax)

        # Shared colorbars for the psi (row 1) and diff (row 2) panels.
        # Position dynamically against the rightmost axis of each row.
        fig2.canvas.draw()
        pos_psi2 = axes_psi2[-1].get_position()
        cax_psi2 = fig2.add_axes([pos_psi2.x1 + 0.012, pos_psi2.y0,
                                  0.012, pos_psi2.height])
        fig2.colorbar(im_psi2, cax=cax_psi2).set_label(
            r'$\Psi(\lambda,\tau)$', fontsize=13)
        if im_diff2_last is not None:
            pos_diff2 = axes_diff2[-1].get_position()
            cax_diff2 = fig2.add_axes([pos_diff2.x1 + 0.012, pos_diff2.y0,
                                       0.012, pos_diff2.height])
            fig2.colorbar(im_diff2_last, cax=cax_diff2).set_label(
                r'$\Delta\Psi$', fontsize=13)

        fig2.suptitle(scan['title'], fontsize=16, y=1.02)
        fname = f'plots/test_lag_spectrum_psi_{scan_name}.png'
        plt.savefig(fname, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Saved {fname}")

    # --- Plot 4: tau ratio (kept from before) ---
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    for idx, (scan_name, scan) in enumerate(scans.items()):
        ax = axes[idx // 2, idx % 2]
        colors = scan_colors(len(scan['values']))
        hlamp_val = scan['fixed'].get('hlamp', 0.5)
        for i, val in enumerate(scan['values']):
            key = f"{scan_name}_{i}"
            tau = data[f'tau_{key}']
            if scan_name == 'vary_hlamp':
                tau_bowl = data[f'tau_bowl_{val}']
            else:
                tau_bowl = data[f'tau_bowl_{hlamp_val}']
            ratio = np.where(tau_bowl > 0, tau / tau_bowl, 1.0)
            ax.plot(wavelengths, ratio, lw=3.5, alpha=0.7, color=colors[i], label=scan['labels'][i])
        ax.axhline(1.0, color='black', ls=':', lw=1, alpha=0.5)
        ax.set_title(scan['title'], fontsize=14, pad=10)
        ax.set_ylabel(r'$\tau_{\rm pillars} / \tau_{\rm bowl}$', fontsize=16)
        if idx >= 2:
            ax.set_xlabel(r'$\lambda~\rm [\AA]$', fontsize=16)
        ax.legend(fontsize=11)
        ax.set_xlim(wavelengths[0], wavelengths[-1])
        ax.set_xscale('log')
        style_ax(ax)
    plt.tight_layout()
    plt.savefig('plots/test_lag_spectrum_ratio.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved plots/test_lag_spectrum_ratio.png")

    # --- Plot 5: Δτ ratio (relative lags) ---
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    for idx, (scan_name, scan) in enumerate(scans.items()):
        ax = axes[idx // 2, idx % 2]
        hlamp_val = scan['fixed'].get('hlamp', 0.5)
        for i, val in enumerate(scan['values']):
            key = f"{scan_name}_{i}"
            tau = data[f'tau_{key}']
            if scan_name == 'vary_hlamp':
                tau_bowl = data[f'tau_bowl_{val}']
            else:
                tau_bowl = data[f'tau_bowl_{hlamp_val}']
            # Relative lags: Δτ = τ(λ) - τ(λ_ref)
            i_ref = np.argmin(np.abs(wavelengths - 1500.0))
            dtau_pillars = tau - tau[i_ref]
            dtau_bowl = tau_bowl - tau_bowl[i_ref]
            ratio = np.where(dtau_bowl > 0.01, dtau_pillars / dtau_bowl, np.nan)
            ax.plot(wavelengths, ratio, lw=3.5, alpha=0.7, label=scan['labels'][i])
        ax.axhline(1.0, color='black', ls=':', lw=1, alpha=0.5)
        ax.set_title(scan['title'], fontsize=14, pad=10)
        ax.set_ylabel(r'$\Delta\tau_{\rm pillars} / \Delta\tau_{\rm bowl}$', fontsize=16)
        if idx >= 2:
            ax.set_xlabel(r'$\lambda~\rm [\AA]$', fontsize=16)
        ax.legend(fontsize=11)
        ax.set_xlim(wavelengths[0], wavelengths[-1])
        ax.set_xscale('log')
        style_ax(ax)
    plt.tight_layout()
    plt.savefig('plots/test_lag_spectrum_delta_ratio.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved plots/test_lag_spectrum_delta_ratio.png")


def compute_sed_from_T(r_grid, T_profile, wavelengths, cosi=0.707):
    """Compute approximate SED from azimuthally-averaged T(r) profile.
    F_nu(lambda) = 2*pi*cosi * integral[ B_nu(lambda, T(r)) * r * dr ]
    """
    H_cgs = 6.626e-27    # erg s
    K_cgs = 1.381e-16    # erg/K
    C_cgs = 2.998e10     # cm/s
    LD_CM = 2.59e15      # light-day in cm
    ANGSTROM = 1e-8      # cm

    F_nu = np.zeros(len(wavelengths))
    for iw, lam in enumerate(wavelengths):
        nu = C_cgs / (lam * ANGSTROM)
        for ir in range(len(r_grid) - 1):
            T = T_profile[ir]
            if T < 10:
                continue
            r_cm = r_grid[ir] * LD_CM
            dr_cm = (r_grid[ir+1] - r_grid[ir]) * LD_CM
            x = H_cgs * nu / (K_cgs * T)
            if x > 500:
                continue
            B_nu = 2 * H_cgs * nu**3 / C_cgs**2 / (np.exp(x) - 1)
            F_nu[iw] += B_nu * r_cm * dr_cm
    F_nu *= 2 * np.pi * cosi
    return F_nu


def plot_sed(data):
    """Plot SED from cached T(r) profiles."""
    wavelengths = np.logspace(np.log10(500), np.log10(50000), 100)  # wider range for SED
    r_grid = data['r_grid']

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    for idx, (scan_name, scan) in enumerate(scans.items()):
        ax = axes[idx // 2, idx % 2]
        hlamp_val = scan['fixed'].get('hlamp', 0.5)

        T_bowl = data[f'T_bowl_{hlamp_val}']
        F_bowl = compute_sed_from_T(r_grid, T_bowl, wavelengths)
        C_cgs = 2.998e10
        ANGSTROM = 1e-8
        nu_grid = C_cgs / (wavelengths * ANGSTROM)

        colors = scan_colors(len(scan['values']))
        for i, val in enumerate(scan['values']):
            key = f"{scan_name}_{i}"
            T = data[f'T_{key}']
            F = compute_sed_from_T(r_grid, T, wavelengths)
            ax.plot(wavelengths, nu_grid * F, lw=3.5, alpha=0.7, color=colors[i], label=scan['labels'][i])

        # Bowl on top — plot last with dashed so it's visible through overlaps
        ax.plot(wavelengths, nu_grid * F_bowl, '--', color='dodgerblue', lw=4, alpha=1.0,
                label=r'$\rm bowl$')

        ax.set_title(scan['title'], fontsize=14, pad=10)
        ax.set_ylabel(r'$\nu F_\nu~\rm [arb.]$', fontsize=16)
        if idx >= 2:
            ax.set_xlabel(r'$\lambda~\rm [\AA]$', fontsize=16)
        ax.legend(fontsize=11)
        ax.set_xlim(wavelengths[0], wavelengths[-1])
        ax.set_xscale('log')
        ax.set_yscale('log')
        style_ax(ax)
    plt.tight_layout()
    plt.savefig('plots/test_lag_spectrum_sed.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved plots/test_lag_spectrum_sed.png")


if __name__ == '__main__':
    if os.path.exists(DATAFILE):
        print(f"Loading cached data from {DATAFILE}")
        data = load_data()
    else:
        data = run_scan()

    plot_all(data)
    plot_sed(data)
