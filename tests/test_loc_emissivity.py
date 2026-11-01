"""Plot LOC-averaged line emissivity, equivalent width, and responsivity
as a function of ionizing flux over the extended Cloudy grid (log Phi =
15 to 21).

Each line is interpolated in log-F space with a natural cubic spline
(C^2 continuous, so the responsivity eta = d log F / d log Phi has a
smooth derivative). This avoids the wild overshoot a linear-space
RectBivariateSpline shows for CIV across its ~7 dex dynamic range.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

from pillardisk.pillar_line_cloudy import load_cloudy_models

cloudy_file = '/Users/jiamuh/c23.01/my_models/loc_metal_flux/strong_LOC_varym_N25_v100_lineflux_LineList_BLR_Fe2_flux.txt'
cloudy_extension_file = '/Users/jiamuh/c23.01/my_models/loc_metal_flux/strong_LOC_varym_N25_v100_lineflux_extlow_LineList_BLR_Fe2_flux.txt'
interp_dict, phi_grid, Z_grid = load_cloudy_models(
    cloudy_file, Z_target=1.0, extension_file=cloudy_extension_file)

# Lines to plot
lines = {
    'Halpha': {'color': 'orangered', 'ls': '-', 'label': r'$\rm H\alpha$'},
    'Mg2':    {'color': 'forestgreen', 'ls': '-', 'label': r'$\rm MgII$'},
    'C4':     {'color': 'royalblue', 'ls': '-', 'label': r'$\rm CIV$'},
    'HI':     {'color': 'darkorange', 'ls': '--', 'label': r'$\rm H\beta$'},
    'HeII':   {'color': 'purple', 'ls': '--', 'label': r'$\rm HeII$'},
    'LyA':    {'color': 'crimson', 'ls': ':', 'label': r'$\rm Ly\alpha$'},
}


def load_vdb_ews(path='data/vandenberk2001_table2.txt'):
    """Parse the Vanden Berk et al. (2001) Table 2 file and return the
    observed composite EWs (Angstroms) for the six lines plotted here.

    Rows are matched by the ID field and, where the same ID appears
    twice (HeII 1640 vs 4687), by the rest wavelength in f_ID.
    """
    # (ID, f_ID rest wavelength or None) -> our line key
    want = {('Ly{alpha}', None): 'LyA',
            ('CIV', None): 'C4',
            ('HeII', '1640.42'): 'HeII',
            ('MgII', None): 'Mg2',
            ('H{beta}', None): 'HI',
            ('H{alpha}', None): 'Halpha'}
    ews = {}
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 12:
                continue
            try:
                float(parts[0])
            except ValueError:
                continue
            line_id, f_id = parts[-2], parts[-1]
            for (wid, wfid), key in want.items():
                if line_id == wid and (wfid is None or f_id == wfid):
                    ews[key] = (float(parts[6]), float(parts[7]))
    return ews


VDB_EW = load_vdb_ews()
print('VdB01 composite EWs [A]:', {k: v[0] for k, v in VDB_EW.items()})

Z_target = 1.0
phi_lo, phi_hi = float(phi_grid[0]), float(phi_grid[-1])
print(f"phi range from Cloudy load: {phi_lo} to {phi_hi}")

# Evaluate emissivity and continuum AT the grid knots, then build a
# 1-D natural cubic spline on log_10 F per line. Sampling at the actual
# grid values just gives back the underlying data, which is what we
# want to interpolate in log space.
def _grid_values(key):
    return np.array([float(interp_dict[key](p, Z_target, grid=False))
                     for p in phi_grid])

line_grid = {key: _grid_values(key) for key in lines if key in interp_dict}
# Incident AGN continuum nu*F_nu at 1215 A [erg cm^-2 s^-1]. Korista-style
# EW convention: EW = F_line * 1215 A / Inci(1215), i.e. the line flux
# relative to the incident continuum flux density at 1215 A, for full
# cloud coverage.
inci_grid = _grid_values('inci1215')
LAM_REF = 1215.0  # A

log_F_interps = {}
log_EW_interps = {}
for key in lines:
    if key not in line_grid:
        continue
    log_F = np.log10(np.maximum(line_grid[key], 1e-300))
    log_F_interps[key] = CubicSpline(phi_grid, log_F, bc_type='natural',
                                      extrapolate=False)
    ew = np.where(inci_grid > 0, line_grid[key] * LAM_REF / inci_grid, np.nan)
    log_ew = np.log10(np.maximum(ew, 1e-300))
    log_EW_interps[key] = CubicSpline(phi_grid, log_ew, bc_type='natural',
                                       extrapolate=False)

# Fine grid for plotting
log_phi = np.linspace(phi_lo, phi_hi, 401)


def make_plot(use_log, filename):
    fig, axes = plt.subplots(3, 1, figsize=(8, 14), sharex=True)

    # ---- Panel 1: surface flux ----
    ax = axes[0]
    for key, props in lines.items():
        if key not in log_F_interps:
            continue
        log_F_fine = log_F_interps[key](log_phi)
        if use_log:
            ax.plot(log_phi, log_F_fine,
                    color=props['color'], ls=props['ls'], lw=3, alpha=0.8,
                    label=props['label'])
        else:
            ax.plot(log_phi, 10.0**log_F_fine,
                    color=props['color'], ls=props['ls'], lw=3, alpha=0.8,
                    label=props['label'])
    if use_log:
        ax.set_ylabel(r'$\log\,F_{\rm line}~\rm [erg~cm^{-2}~s^{-1}]$',
                      fontsize=16)
    else:
        ax.set_ylabel(r'$F_{\rm line}~\rm [erg~cm^{-2}~s^{-1}]$', fontsize=16)
    ax.text(0.03, 0.94, r'$\rm LOC~surface~flux$', transform=ax.transAxes,
            fontsize=15, va='top', ha='left',
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3', alpha=0.9))
    ax.set_xlim(phi_lo, phi_hi)
    ax.tick_params(top=True, right=True, direction='in', which='major',
                   length=8, width=1.5, labelsize=13)
    ax.tick_params(top=True, right=True, direction='in', which='minor',
                   length=5, width=1.0)
    ax.minorticks_on()

    # ---- Panel 2: equivalent width ----
    ax = axes[1]
    for key, props in lines.items():
        if key not in log_EW_interps:
            continue
        log_ew_fine = log_EW_interps[key](log_phi)
        if use_log:
            ax.plot(log_phi, log_ew_fine,
                    color=props['color'], ls=props['ls'], lw=3, alpha=0.8,
                    label=props['label'])
        else:
            ax.plot(log_phi, 10.0**log_ew_fine,
                    color=props['color'], ls=props['ls'], lw=3, alpha=0.8,
                    label=props['label'])
    # Observed composite EWs (Vanden Berk et al. 2001, Table 2),
    # drawn as horizontal dotted lines in the matching line colors,
    # for CIV, MgII and Halpha only.
    # for key in ['C4', 'Mg2', 'Halpha']:
    #     if key not in VDB_EW:
    #         continue
    #     props = lines[key]
    #     ew_obs = VDB_EW[key][0]
    #     y = np.log10(ew_obs) if use_log else ew_obs
    #     ax.axhline(y, color=props['color'], ls=':', lw=1.8, alpha=0.7)
    if use_log:
        ax.set_ylabel(r'$\log\,\rm EW~[\AA]$', fontsize=16)
    else:
        ax.set_ylabel(r'$\rm EW~[\AA]$', fontsize=16)
    ax.text(0.03, 0.94, r'$\rm Equivalent~width$', transform=ax.transAxes,
            fontsize=15, va='top', ha='left',
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3', alpha=0.9))
    ax.set_xlim(phi_lo, phi_hi)
    if use_log:
        ax.set_ylim(bottom=-1)
    ax.tick_params(top=True, right=True, direction='in', which='major',
                   length=8, width=1.5, labelsize=13)
    ax.tick_params(top=True, right=True, direction='in', which='minor',
                   length=5, width=1.0)
    ax.minorticks_on()

    # ---- Panel 3: responsivity ----
    ax = axes[2]
    for key, props in lines.items():
        if key not in log_F_interps:
            continue
        eta_fine = log_F_interps[key].derivative()(log_phi)
        ax.plot(log_phi, eta_fine,
                color=props['color'], ls=props['ls'], lw=3, alpha=0.8,
                label=props['label'])
    ax.set_xlabel(r'$\log\Phi~\rm [photons~cm^{-2}~s^{-1}]$', fontsize=16)
    ax.set_ylabel(r'$\eta = d\log F / d\log\Phi$', fontsize=16)
    ax.text(0.03, 0.94, r'$\rm Responsivity$', transform=ax.transAxes,
            fontsize=15, va='top', ha='left',
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3', alpha=0.9))
    ax.set_xlim(phi_lo, phi_hi)
    ax.tick_params(top=True, right=True, direction='in', which='major',
                   length=8, width=1.5, labelsize=13)
    ax.tick_params(top=True, right=True, direction='in', which='minor',
                   length=5, width=1.0)
    ax.minorticks_on()

    handles, labels = axes[0].get_legend_handles_labels()
    ncol = int(np.ceil(len(labels) / 2))
    # Reorder so the solid-line entries occupy the top row of a column-major
    # legend (matplotlib fills column-by-column, so we interleave).
    n = len(labels)
    top = list(range(ncol))
    bot = list(range(ncol, n))
    while len(bot) < ncol:
        bot.append(None)
    order = []
    for t, b in zip(top, bot):
        order.append(t)
        if b is not None:
            order.append(b)
    handles = [handles[i] for i in order]
    labels = [labels[i] for i in order]
    fig.legend(handles, labels, fontsize=16, ncol=ncol,
               loc='upper center', bbox_to_anchor=(0.5, 1.04),
               frameon=False)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.05)
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved {filename}")


make_plot(use_log=False, filename='plots/test_loc_emissivity_linear.png')
make_plot(use_log=True, filename='plots/test_loc_emissivity_log.png')
