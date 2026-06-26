#!/usr/bin/env python3
"""
pillar_line_cloudy.py - Integrate Cloudy models with pillar disk geometry

This script loads Cloudy model data and integrates it with the pillar disk model.
The ionizing flux is mapped to temperature (highest at illuminated pillars, lowest in shadows).
Line intensities are interpolated from Cloudy models as a function of ionizing flux and density.

Output: RGB 2D map where R=Halpha, G=MgII, B=CIV
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import pandas as pd
from scipy.interpolate import RectBivariateSpline, griddata

# Try to import YAML for config file
try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False
    print("Warning: PyYAML not installed. Install with 'pip install pyyaml' to use config files.")

# Try to import tqdm for progress bars
try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False
    def tqdm(iterable, *args, **kwargs):
        return iterable

# Import the PillarDisk class and utilities

try:
    from pillardisk.pillar_disk import PillarDisk, load_config, resolve_rin
except ImportError as e:
    print(f"Error: pillar_disk.py not found. Make sure it's in the same directory.")
    print(f"Import error: {e}")
    sys.exit(1)

# Define parse_math_expr
def parse_math_expr(value):
    """Parse mathematical expressions like 'pi', 'pi/2', etc."""
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        expr = value.lower().strip()
        expr = expr.replace('pi', str(np.pi))
        expr = expr.replace('π', str(np.pi))
        try:
            safe_dict = {
                'pi': np.pi, 'π': np.pi, 'e': np.e,
                'sin': np.sin, 'cos': np.cos, 'tan': np.tan,
                'sqrt': np.sqrt, 'exp': np.exp, 'log': np.log,
            }
            result = eval(expr, {"__builtins__": {}}, safe_dict)
            return float(result)
        except:
            return float(value)
    return float(value)

# Physical constants (cgs units)
C = 2.997925e10  # speed of light (cm/s)
PC = 3.0857e18   # parsec (cm)
DAY = 24. * 3600.  # day (s)
KM_TO_CM = 1e5    # km to cm conversion

# Conversion factors
PC_TO_LD = 1e6 * PC / (C * DAY)  # parsec to light days
LD_TO_CM = C * DAY  # light days to cm


def load_cloudy_models(model_file, Z_target=1.0, extension_file=None):
    """
    Load Cloudy model data and create interpolation functions for solar metallicity (Z=1).

    Parameters:
    -----------
    model_file : str
        Path to the main Cloudy model data file (covers log phi in [17, 21]).
    Z_target : float
        Target metallicity (default: 1.0 for solar).
    extension_file : str or None
        Optional path to a low-phi extension Cloudy file covering log phi
        in [15, 16.5] with Z = Z_sun only. When supplied, its rows are
        spliced onto the front of the main grid so the resulting grid
        spans log phi in [15, 21]. If None, the loader falls back to the
        previous log-linear extrapolation between 16 and the main-file
        floor (default 17). Because the extension is Z = Z_sun only,
        non-solar Z columns at the extended phi values are populated
        with the same Z_sun value (this is fine for analyses that only
        use Z = Z_sun).

    Returns:
    --------
    interp_dict : dict
        Dictionary of interpolation functions: {'Mg2': interp, 'C4': interp, ...}
    phi_grid : array
        Ionizing flux grid (log phi)
    Z_grid : array
        Metallicity grid
    """
    print(f"Loading Cloudy models from {model_file}...")

    if not os.path.exists(model_file):
        raise FileNotFoundError(f"Cloudy model file not found: {model_file}")

    # Grid definitions — phi and n are fixed; Z varies by model file
    phi_grid = np.arange(17, 21.5, 0.5)       # 9 φ values
    n_grid = np.arange(9, 12.5, 0.5)          # 7 n values
    
    # Line names mapping
    line_names = {
        "LyA": 'H  1 1215.67A',
        'Halpha': 'h  1 6562.80A',  # Halpha for RGB red channel
        'HI':  'H  1 4861.32A',  # Hbeta
        'HeII': 'he 2 1640.41A',
        'Mg2': 'blnd 2798.00A',  # MgII for RGB green channel
        'C4':  'blnd 1549.00A',  # CIV for RGB blue channel
        'Si4': 'BLND 1397.00A',
        'O4':  'blnd 1402.00A',
        'C3':  'blnd 1909.00A',
        'N5':  'blnd 1240.00A',
        'N4':  'blnd 1486.00A',
        'N3':  'blnd 1750.00A',
        'nufnu4885': 'nFnu 4885.36A',
        'nufnu1793': 'nFnu 1793.44A',
        'nufnu1458': 'nFnu 1458.33A',
        'nufnu1262': 'nFnu 1262.79A',
        'inci1215': 'Inci 1215.00A',
    }
    
    # Load data — handle both '#lineslist' and 'lineslist' headers
    with open(model_file) as fh:
        header_line = fh.readline().lstrip('#').strip()
    col_names = header_line.split('\t')
    df = pd.read_csv(model_file, sep='\t', names=col_names, comment='#', skiprows=0)
    # Keep only 'iteration 1' rows (drop GRID_DELIMIT separators)
    df = df[df.iloc[:, 0] == 'iteration 1'].reset_index(drop=True)
    # Auto-detect Z grid size from row count
    n_Z = len(df) // (len(phi_grid) * len(n_grid))
    Z_grid = np.arange(1, 1 + n_Z * 0.5, 0.5)
    expected = len(Z_grid) * len(phi_grid) * len(n_grid)
    assert len(df) == expected, f"Row count mismatch: got {len(df)}, expected {expected} (Z:{n_Z}, phi:{len(phi_grid)}, n:{len(n_grid)})"
    
    # Initialize 3D arrays: shape = (len(n), len(phi), len(Z))
    shape = (len(n_grid), len(phi_grid), len(Z_grid))
    intensity_data = {key: np.zeros(shape) for key in line_names}
    
    # Fill 3D arrays
    for j_phi, phi_val in enumerate(phi_grid):
        for i_n, n_val in enumerate(n_grid):
            index = j_phi * len(n_grid) + i_n
            start = index * len(Z_grid)
            end = start + len(Z_grid)
            df_slice = df.iloc[start:end]
            for key, colname in line_names.items():
                intensity_data[key][i_n, j_phi, :] = df_slice[colname].values
    
    print("Collapsing over density axis...")
    # Collapse over density (marginalize): average along axis=0
    collapsed_data = {}
    for key in line_names:
        collapsed_data[key] = np.mean(intensity_data[key], axis=0)  # Shape: (len(phi), len(Z))
    
    # For Z=1 (solar), extract the column at Z=1.0
    # Find closest Z value to Z_target
    Z_idx = np.argmin(np.abs(Z_grid - Z_target))
    Z_actual = Z_grid[Z_idx]
    print(f"Using Z = {Z_actual:.1f} (closest to target Z = {Z_target:.1f})")
    
    # Extract line intensities at Z=1 for all phi values
    line_intensities = {}
    for key in line_names:
        # Extract column at Z=Z_idx: shape (len(phi),)
        line_intensities[key] = collapsed_data[key][:, Z_idx]
    
    # ----- Extend the grid down to log phi = 15 -----
    # Preferred: splice in a real low-phi Cloudy extension file (Z = Z_sun only).
    # Fallback: log-linear extrapolation between 16.0 and the first main-grid
    # point (used when no extension file is given).
    if extension_file is not None:
        if not os.path.exists(extension_file):
            raise FileNotFoundError(f"Cloudy extension file not found: {extension_file}")
        print(f"Loading low-phi extension from {extension_file}...")

        # Extension grid: phi in [15, 16.5] step 0.5 (4 values); Z = Z_sun only
        phi_grid_extlow = np.arange(15.0, 17.0, 0.5)
        n_Z_ext = 1

        with open(extension_file) as fh_ext:
            header_line_ext = fh_ext.readline().lstrip('#').strip()
        col_names_ext = header_line_ext.split('\t')
        df_ext = pd.read_csv(extension_file, sep='\t', names=col_names_ext, comment='#', skiprows=0)
        df_ext = df_ext[df_ext.iloc[:, 0] == 'iteration 1'].reset_index(drop=True)
        expected_ext = n_Z_ext * len(phi_grid_extlow) * len(n_grid)
        assert len(df_ext) == expected_ext, (
            f"Extension row count mismatch: got {len(df_ext)}, expected {expected_ext} "
            f"(phi:{len(phi_grid_extlow)}, n:{len(n_grid)}, Z:{n_Z_ext})"
        )

        # Pack into 3D array (n_density, n_phi_ext, 1) following the same
        # iteration order as the main file: outermost = phi, then n, then Z.
        ext_shape = (len(n_grid), len(phi_grid_extlow), n_Z_ext)
        ext_intensity = {key: np.zeros(ext_shape) for key in line_names}
        for j_phi, _ in enumerate(phi_grid_extlow):
            for i_n, _ in enumerate(n_grid):
                index = j_phi * len(n_grid) + i_n
                start = index * n_Z_ext
                df_slice = df_ext.iloc[start:start + n_Z_ext]
                for key, colname in line_names.items():
                    ext_intensity[key][i_n, j_phi, :] = df_slice[colname].values

        # Average over density (same as main file)
        ext_collapsed = {key: np.mean(ext_intensity[key], axis=0)  # shape (n_phi_ext, 1)
                          for key in line_names}

        # Splice: prepend extension to main grid. The extension has only
        # the Z = Z_sun column, so for the existing Z grid (which has
        # multiple Z values) we replicate the Z_sun value across all
        # columns. This is fine because the paper only ever queries
        # Z = 1; non-solar Z below log phi = 17 is undefined.
        extended_data = {}
        for key in line_names:
            ext_block = ext_collapsed[key]  # (n_phi_ext, 1)
            ext_replicated = np.tile(ext_block, (1, len(Z_grid)))  # (n_phi_ext, n_Z)
            extended_data[key] = np.vstack([ext_replicated, collapsed_data[key]])
        phi_grid_ext = np.concatenate([phi_grid_extlow, phi_grid])
        print(f"Spliced extension: phi grid now {phi_grid_ext[0]:.1f} -> {phi_grid_ext[-1]:.1f} "
              f"({len(phi_grid_ext)} values; ext Z = Z_sun replicated to all Z columns)")
    else:
        # Fallback: log-linear extrapolation between 16.0 and phi_grid[0]
        phi_extrap = np.arange(16.0, phi_grid[0], 0.5)  # [16.0, 16.5]
        if len(phi_extrap) > 0:
            extended_data = {}
            for key in line_names:
                vals = collapsed_data[key]  # shape: (len(phi), len(Z))
                v0 = vals[0, :]  # at phi_grid[0]
                v1 = vals[1, :]  # at phi_grid[1]
                safe_v0 = np.maximum(v0, 1e-30)
                safe_v1 = np.maximum(v1, 1e-30)
                log_slope = (np.log10(safe_v1) - np.log10(safe_v0)) / (phi_grid[1] - phi_grid[0])
                extrap_rows = []
                for phi_e in phi_extrap:
                    dphi = phi_e - phi_grid[0]
                    log_extrap = np.log10(safe_v0) + log_slope * dphi
                    extrap_rows.append(np.maximum(10.0**log_extrap, 0.0))
                extended_data[key] = np.vstack(extrap_rows + [vals])  # prepend
            phi_grid_ext = np.concatenate([phi_extrap, phi_grid])
            print(f"Extended Cloudy grid: phi {phi_grid[0]:.1f} -> {phi_grid_ext[0]:.1f} (log-linear extrapolation)")
        else:
            extended_data = collapsed_data
            phi_grid_ext = phi_grid

    # Create interpolation functions over extended phi grid
    interp_dict = {}
    for key in line_names:
        interp_dict[key] = RectBivariateSpline(phi_grid_ext, Z_grid, extended_data[key], kx=2, ky=2)

    print(f"Loaded Cloudy models: {list(line_names.keys())}")
    print(f"  phi range: {phi_grid_ext[0]:.1f} to {phi_grid_ext[-1]:.1f} (log)")
    print(f"  Z range: {Z_grid[0]:.1f} to {Z_grid[-1]:.1f}")

    # Debug table (verbose only)
    if os.environ.get('PILLAR_DEBUG'):
        print(f"\n  [DEBUG] Line emissivity vs log_phi at Z={Z_actual:.1f}:")
        print(f"  log_phi  |   CIV    |   MgII   |  Halpha")
        print(f"  ---------|----------|----------|----------")
        for i, phi_val in enumerate(phi_grid):
            c4_val = collapsed_data['C4'][i, Z_idx]
            mg2_val = collapsed_data['Mg2'][i, Z_idx]
            ha_val = collapsed_data['Halpha'][i, Z_idx]
            print(f"    {phi_val:.1f}   | {c4_val:8.2e} | {mg2_val:8.2e} | {ha_val:8.2e}")
        print()

    return interp_dict, phi_grid_ext, Z_grid


def compute_ionizing_flux(disk, r, phi):
    """Thin wrapper around PillarDisk.compute_log_ionizing_flux for back-compat."""
    return disk.compute_log_ionizing_flux(r, phi)


def compute_line_intensity_map(disk, cloudy_interp_dict, phi_grid, Z_target=1.0, 
                               n_r_plot=200, n_phi_plot=360):
    """
    Compute line intensity map and equivalent widths using Cloudy models.
    
    Parameters:
    -----------
    disk : PillarDisk
        Disk model
    cloudy_interp_dict : dict
        Dictionary of Cloudy interpolation functions
    phi_grid : array
        Ionizing flux grid (log phi)
    Z_target : float
        Target metallicity (default: 1.0 for solar)
    n_r_plot : int
        Number of radial points for plotting
    n_phi_plot : int
        Number of azimuthal points for plotting
    
    Returns:
    --------
    r_plot : array
        Radial grid for plotting (1D)
    phi_plot : array
        Azimuthal grid for plotting (1D)
    intensity_map : dict
        Dictionary of line intensity maps: {'Mg2': map, 'C4': map, ...}
    ew_map : dict
        Dictionary of equivalent width maps: EW = F_line / F_continuum
    log_phi_2d : array
        Ionizing flux map (2D)
    """
    print("Computing line intensity map...")
    
    # Create fine grid for plotting
    r_plot = np.logspace(np.log10(disk.rin), np.log10(disk.rout), n_r_plot)
    phi_plot = np.linspace(0, 2*np.pi, n_phi_plot, endpoint=False)
    r_2d, phi_2d = np.meshgrid(r_plot, phi_plot, indexing='ij')
    
    # Compute ionizing flux at each position
    log_phi_2d = compute_ionizing_flux(disk, r_2d, phi_2d)
    
    # Interpolate line intensities from Cloudy models
    intensity_map = {}
    for key, interp_func in cloudy_interp_dict.items():
        print(f"  Interpolating {key}...")
        # Interpolate at (log_phi, Z=1.0)
        # grid=False returns scalar for each point
        intensity_2d = np.zeros_like(log_phi_2d)
        no_fluxfloor = getattr(disk, 'no_fluxfloor', False)
        for ir in range(n_r_plot):
            for iphi in range(n_phi_plot):
                log_phi_val = log_phi_2d[ir, iphi]
                if no_fluxfloor and log_phi_val < 15.0:
                    intensity_2d[ir, iphi] = 0.0
                    continue
                log_phi_val = np.clip(log_phi_val, phi_grid[0], phi_grid[-1])
                intensity_2d[ir, iphi] = float(interp_func(log_phi_val, Z_target, grid=False))
        intensity_map[key] = intensity_2d
    
    # Compute equivalent widths: EW = F_line / F_continuum
    # Use nufnu4885 as continuum normalization
    print("Computing equivalent widths...")
    ew_map = {}
    if 'nufnu4885' in intensity_map:
        continuum = intensity_map['nufnu4885']
        # Define emission lines (exclude continuum measurements)
        emission_lines = ['Halpha', 'LyA', 'HI', 'HeII', 'Mg2', 'C4', 'Si4', 'O4', 'C3', 'N5', 'N4', 'N3']
        for line_key in emission_lines:
            if line_key in intensity_map:
                # EW = F_line / F_continuum (in units of continuum)
                # Avoid division by zero
                ew_2d = np.zeros_like(continuum)
                mask = continuum > 0
                ew_2d[mask] = intensity_map[line_key][mask] / continuum[mask]
                ew_2d[~mask] = 0.0
                ew_map[line_key] = ew_2d
                print(f"  Computed EW for {line_key}")
    else:
        print("Warning: nufnu4885 continuum not found. Cannot compute equivalent widths.")
    
    return r_plot, phi_plot, intensity_map, ew_map, log_phi_2d


def plot_rgb_line_map(r_plot, phi_plot, intensity_map, ew_map, log_phi_2d, phi_grid, use_absolute_flux=True, filename='line_intensity_rgb.png', rgb_sat_mode='shared', composite_intensity_map=None, composite_brightness='none', composite_sat_boost=1.0, composite_brightness_gamma=1.0, composite_brightness_floor=0.0, composite_brightness_boost=1.0, composite_brightness_phi_min=16.0, composite_brightness_phi_max=21.0, composite_brightness_add=0.0, pillar_rmin=None):
    """
    Plot RGB map of line intensities: R=Halpha, G=MgII, B=CIV.
    
    Note: If Halpha is not available, we'll use a proxy or note it in the plot.
    
    Parameters:
    -----------
    r_plot : array
        Radial grid (1D)
    phi_plot : array
        Azimuthal grid (1D)
    intensity_map : dict
        Dictionary of line intensity maps
    filename : str
        Output filename
    """
    print("Creating RGB line intensity map...")
    
    # Check available lines for RGB: Halpha (R), MgII (G), CIV (B)
    has_halpha = 'Halpha' in intensity_map
    has_mg2 = 'Mg2' in intensity_map
    has_c4 = 'C4' in intensity_map
    
    if not has_halpha:
        print("Error: Halpha not found in Cloudy models.")
        return
    
    if not has_mg2:
        print("Error: MgII not found in Cloudy models.")
        return
    
    if not has_c4:
        print("Error: CIV not found in Cloudy models.")
        return
    
    # Use correct lines for RGB.
    # `intensity_map` drives the per-line panels. `composite_intensity_map`
    # (if given) drives the RGB composite — used to mix Hβ-normalized
    # per-line panels with an absolute-flux composite in a single figure.
    R = intensity_map['Halpha']  # Red channel: Halpha
    G = intensity_map['Mg2']     # Green channel: MgII
    B = intensity_map['C4']      # Blue channel: CIV
    if composite_intensity_map is not None:
        R_comp = composite_intensity_map['Halpha']
        G_comp = composite_intensity_map['Mg2']
        B_comp = composite_intensity_map['C4']
    else:
        R_comp, G_comp, B_comp = R, G, B
    
    # Compute log scale values for plotting (absolute values, not normalized)
    # For RGB display, we still need normalized values [0,1]
    def compute_log_values(channel):
        """Compute log10 of channel values for plotting."""
        channel_pos = channel[channel > 0]
        if len(channel_pos) == 0:
            return np.full_like(channel, np.nan), np.nan, np.nan
        channel_min = np.nanmin(channel_pos)
        channel_max = np.nanmax(channel)
        # Log scale values
        log_channel = np.log10(channel + 1e-10 * channel_min)  # Add small offset to avoid log(0)
        log_min = np.log10(channel_min + 1e-10 * channel_min)
        log_max = np.log10(channel_max + 1e-10 * channel_min)
        return log_channel, log_min, log_max
    
    # Shared color range used by the per-channel panels.
    if use_absolute_flux and composite_intensity_map is None:
        all_log = np.concatenate([
            compute_log_values(R)[0].ravel(),
            compute_log_values(G)[0].ravel(),
            compute_log_values(B)[0].ravel(),
        ])
        all_log = all_log[np.isfinite(all_log)]
        VMIN_SHARED = float(np.floor(np.nanpercentile(all_log, 2)))
        VMAX_SHARED = float(np.ceil(np.nanpercentile(all_log, 98)))
    else:
        # Per-line panels are F/F_Hβ (either because use_absolute_flux=False,
        # or because we're in hybrid mode where intensity_map is normalized
        # and only composite_intensity_map is absolute).
        VMIN_SHARED = -0.2
        VMAX_SHARED = 1.0

    def _channel_range(channel, mode, fallback_vmin, fallback_vmax):
        log_c = compute_log_values(channel)[0]
        finite = log_c[np.isfinite(log_c)]
        if finite.size == 0:
            return fallback_vmin, fallback_vmax
        if mode == 'percentile':
            return float(np.nanpercentile(finite, 2)), float(np.nanpercentile(finite, 98))
        elif mode == 'minmax':
            return float(np.min(finite)), float(np.max(finite))
        else:
            return fallback_vmin, fallback_vmax

    def get_clipped_range(log_channel):
        """Return the shared (vmin, vmax) used for the per-channel panels."""
        return VMIN_SHARED, VMAX_SHARED

    def _normalize_with_range(channel, vmin, vmax):
        log_channel, log_min, log_max = compute_log_values(channel)
        if np.isnan(log_min) or vmax <= vmin:
            return np.zeros_like(channel)
        return np.clip((log_channel - vmin) / (vmax - vmin), 0, 1)

    def normalize_channel_log(channel):
        """Normalize channel to [0,1] for RGB display using the shared
        log10 range so the three RGB planes are directly comparable."""
        return _normalize_with_range(channel, VMIN_SHARED, VMAX_SHARED)

    # --- RGB composite normalization ---
    # The composite uses `R_comp/G_comp/B_comp` (absolute flux in hybrid
    # mode, otherwise the same as R/G/B). When rgb_sat_mode != 'shared',
    # each channel is normalized over its own dynamic range so each color
    # saturates independently.
    if composite_intensity_map is not None or use_absolute_flux:
        comp_all_log = np.concatenate([
            compute_log_values(R_comp)[0].ravel(),
            compute_log_values(G_comp)[0].ravel(),
            compute_log_values(B_comp)[0].ravel(),
        ])
        comp_all_log = comp_all_log[np.isfinite(comp_all_log)]
        VMIN_COMP_SHARED = float(np.floor(np.nanpercentile(comp_all_log, 2)))
        VMAX_COMP_SHARED = float(np.ceil(np.nanpercentile(comp_all_log, 98)))
    else:
        VMIN_COMP_SHARED, VMAX_COMP_SHARED = VMIN_SHARED, VMAX_SHARED

    if rgb_sat_mode in ('percentile', 'minmax'):
        Rv = _channel_range(R_comp, rgb_sat_mode, VMIN_COMP_SHARED, VMAX_COMP_SHARED)
        Gv = _channel_range(G_comp, rgb_sat_mode, VMIN_COMP_SHARED, VMAX_COMP_SHARED)
        Bv = _channel_range(B_comp, rgb_sat_mode, VMIN_COMP_SHARED, VMAX_COMP_SHARED)
        R_norm = _normalize_with_range(R_comp, Rv[0], Rv[1])
        G_norm = _normalize_with_range(G_comp, Gv[0], Gv[1])
        B_norm = _normalize_with_range(B_comp, Bv[0], Bv[1])
    else:
        R_norm = _normalize_with_range(R_comp, VMIN_COMP_SHARED, VMAX_COMP_SHARED)
        G_norm = _normalize_with_range(G_comp, VMIN_COMP_SHARED, VMAX_COMP_SHARED)
        B_norm = _normalize_with_range(B_comp, VMIN_COMP_SHARED, VMAX_COMP_SHARED)

    # Optional brightness weighting via HSV: keep the hue and saturation
    # of the F/F_Hβ-derived color (the differential-response info), and
    # set the V (value/brightness) channel from a separate mask such as
    # log Φ. Simply multiplying R/G/B by brightness preserves the ratios
    # mathematically but compresses chroma toward black in dim regions,
    # which hides the color contrast visually.
    brightness_mask = None
    if composite_brightness == 'phi':
        # V = (logΦ − phi_min) / (phi_max − phi_min), clipped to [0, 1].
        # Then optional gamma, then floor lift + boost.
        phi_min_v = float(composite_brightness_phi_min)
        phi_max_v = float(composite_brightness_phi_max)
        brightness_mask = np.clip((log_phi_2d - phi_min_v) / (phi_max_v - phi_min_v), 0, 1)
        if composite_brightness_gamma != 1.0:
            brightness_mask = brightness_mask ** float(composite_brightness_gamma)
        floor = float(composite_brightness_floor)
        if floor > 0:
            brightness_mask = floor + (1.0 - floor) * brightness_mask
        if composite_brightness_boost != 1.0:
            brightness_mask = np.clip(brightness_mask * float(composite_brightness_boost), 0, 1)
        if composite_brightness_add != 0.0:
            brightness_mask = np.clip(brightness_mask + float(composite_brightness_add), 0, 1)
    
    # For plotting absolute values (log scale)
    R_log, R_log_min, R_log_max = compute_log_values(R)
    G_log, G_log_min, G_log_max = compute_log_values(G)
    B_log, B_log_min, B_log_max = compute_log_values(B)
    
    # Create RGB image
    rgb_image = np.zeros((len(r_plot), len(phi_plot), 3))
    rgb_image[:, :, 0] = R_norm  # Red
    rgb_image[:, :, 1] = G_norm  # Green
    rgb_image[:, :, 2] = B_norm  # Blue

    # If a brightness mask was supplied, replace the V (value) channel of
    # the RGB-as-HSV so hue + saturation (the color contrast) survive
    # while brightness encodes the supplied scalar (e.g. log Φ). Optional
    # composite_sat_boost (>1) inflates the S channel for more vivid hues.
    if brightness_mask is not None or composite_sat_boost != 1.0:
        from matplotlib.colors import rgb_to_hsv, hsv_to_rgb
        hsv = rgb_to_hsv(np.clip(rgb_image, 0, 1))
        if composite_sat_boost != 1.0:
            hsv[:, :, 1] = np.clip(hsv[:, :, 1] * composite_sat_boost, 0, 1)
        if brightness_mask is not None:
            hsv[:, :, 2] = brightness_mask
        rgb_image = hsv_to_rgb(hsv)
    
    # Convert to Cartesian coordinates (x-y plane) for plotting
    # Create meshgrid in polar coordinates, then convert to Cartesian
    R_mesh, Phi_mesh = np.meshgrid(r_plot, phi_plot, indexing='ij')
    X = R_mesh * np.cos(Phi_mesh)  # x coordinates (light days)
    Y = R_mesh * np.sin(Phi_mesh)  # y coordinates (light days)
    
    # Determine plot extent
    x_min, x_max = X.min(), X.max()
    y_min, y_max = Y.min(), Y.max()
    
    # For Cloudy interpolation, clip to grid range; for plotting, use actual data range
    phi_min_grid = phi_grid[0]  # Cloudy grid minimum (17.0)
    log_phi_2d_clipped = np.clip(log_phi_2d, phi_min_grid, None)
    # Dynamic colorbar range: exclude unphysical pillar-wall values (log_phi < 0)
    finite_phi = log_phi_2d[np.isfinite(log_phi_2d) & (log_phi_2d > 0)]
    phi_plot_min = np.floor(np.min(finite_phi) * 2) / 2
    phi_plot_max = np.ceil(np.max(finite_phi) * 2) / 2
    
    # Prepare interpolation grid (used for both plots)
    x_flat = X.flatten()
    y_flat = Y.flatten()
    x_reg = np.linspace(x_min, x_max, len(r_plot))
    y_reg = np.linspace(y_min, y_max, len(r_plot))
    X_reg, Y_reg = np.meshgrid(x_reg, y_reg)
    
    # Colorbar labels based on what the per-line panels actually show.
    # In hybrid mode the per-line panels use F/F_Hβ even when the composite
    # uses absolute flux.
    per_line_is_absolute = use_absolute_flux and (composite_intensity_map is None)
    if per_line_is_absolute:
        # Cloudy "save line list absolute" → erg s^-1 cm^-2 at the cloud surface
        label_R = r'$\rm \log(F_{H\alpha})$'
        label_G = r'$\rm \log(F_{MgII})$'
        label_B = r'$\rm \log(F_{CIV})$'
    else:
        # Cloudy normalized output: F_Hbeta = 1, so each line flux is F/F_Hbeta
        label_R = r'$\rm \log(F_{H\alpha}/F_{H\beta})$'
        label_G = r'$\rm \log(F_{MgII}/F_{H\beta})$'
        label_B = r'$\rm \log(F_{CIV}/F_{H\beta})$'
    
    # Separate into three plots:
    # 1. Ionizing flux (separate plot)
    # 2. RGB + 3 channels in x-y plane
    # 3. RGB + 3 channels in r-phi plane
    
    # Per-block enable flag: set _PLOT_EXTRAS=True to re-enable the
    # ionizing-flux, r-phi, logF-logphi, and EW figures. While iterating
    # on the xy RGB composite we keep them off to save plotting time.
    _PLOT_EXTRAS = False

    # ====================================================
    # Plot 1: Ionizing flux (separate plot)
    # ====================================================
    if _PLOT_EXTRAS:
        fig_phi, axes_phi = plt.subplots(1, 2, figsize=(14, 6))

        # Left: x-y plane
        ax = axes_phi[0]
        log_phi_flat = log_phi_2d.flatten()
        log_phi_reg = griddata((x_flat, y_flat), log_phi_flat, (X_reg, Y_reg), method='linear', fill_value=phi_plot_min)
        im_phi_xy = ax.contourf(X_reg, Y_reg, log_phi_reg, levels=50, cmap='viridis', extend='neither',
                                vmin=phi_plot_min, vmax=phi_plot_max)
        ax.set_xlabel(r'$\rm X~[light~days]$', fontsize=14)
        ax.set_ylabel(r'$\rm Y~[light~days]$', fontsize=14)
        ax.set_title(r'$\rm Log~Ionizing~Flux~\log\phi~(x-y~plane)$', fontsize=14)
        plt.colorbar(im_phi_xy, ax=ax, label=r'$\rm \log\phi$', shrink=0.8, aspect=20)
        ax.set_aspect('equal')

        # Right: r-phi plane
        ax = axes_phi[1]
        im_phi_rphi = ax.imshow(log_phi_2d, extent=[phi_plot[0], phi_plot[-1], r_plot[0], r_plot[-1]],
                               aspect='auto', origin='lower', interpolation='bilinear', cmap='viridis',
                               vmin=phi_plot_min, vmax=phi_plot_max)
        ax.set_xlabel(r'$\rm Azimuthal~Angle~\phi~[rad]$', fontsize=14)
        ax.set_ylabel(r'$\rm Radius~[light~days]$', fontsize=14)
        ax.set_title(r'$\rm Log~Ionizing~Flux~\log\phi~(r-\phi~plane)$', fontsize=14)
        plt.colorbar(im_phi_rphi, ax=ax, label=r'$\rm \log\phi$', shrink=0.8, aspect=20)

        plt.tight_layout()
        filename_phi = filename.replace('.png', '_ionizing_flux.png')
        plt.savefig(filename_phi, dpi=150, bbox_inches='tight')
        print(f"Ionizing flux map saved to {filename_phi}")
        plt.close()
    
    # ====================================================
    # Plot 2: RGB + 3 channels in x-y plane
    # ====================================================
    fig1, axes1 = plt.subplots(2, 2, figsize=(14, 13),
                                gridspec_kw={'hspace': -0.05, 'wspace': 0.35})
    
    # Interpolate RGB image to regular grid
    rgb_r_flat = rgb_image[:, :, 0].flatten()
    rgb_g_flat = rgb_image[:, :, 1].flatten()
    rgb_b_flat = rgb_image[:, :, 2].flatten()
    
    rgb_r_reg = griddata((x_flat, y_flat), rgb_r_flat, (X_reg, Y_reg), method='linear', fill_value=0)
    rgb_g_reg = griddata((x_flat, y_flat), rgb_g_flat, (X_reg, Y_reg), method='linear', fill_value=0)
    rgb_b_reg = griddata((x_flat, y_flat), rgb_b_flat, (X_reg, Y_reg), method='linear', fill_value=0)
    
    rgb_image_reg = np.zeros((rgb_r_reg.shape[0], rgb_r_reg.shape[1], 3))
    rgb_image_reg[:, :, 0] = np.clip(rgb_r_reg, 0, 1)
    rgb_image_reg[:, :, 1] = np.clip(rgb_g_reg, 0, 1)
    rgb_image_reg[:, :, 2] = np.clip(rgb_b_reg, 0, 1)
    
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from matplotlib.ticker import AutoLocator

    TITLE_FS = 22
    TICK_FS = 16
    CBAR_LABEL_FS = 20
    SUPLABEL_FS = 22

    def _attach_cax(ax, im, label):
        """Attach a colorbar matched to ax height; return the cax."""
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.1)
        if im is None:
            cax.set_visible(False)
            return cax
        cbar = plt.colorbar(im, cax=cax)
        cbar.set_label(label, fontsize=CBAR_LABEL_FS)
        cbar.ax.tick_params(labelsize=TICK_FS)
        # Override the default contour-level ticks (which give ugly fractional
        # numbers) with matplotlib's standard auto-locator → clean values
        # like 0, 0.2, 0.4, 0.6, 0.8, 1.0.
        cbar.locator = AutoLocator()
        cbar.update_ticks()
        return cax

    def _style_panel(ax, title, hide_xticks=False, hide_yticks=False):
        ax.set_title(title, fontsize=TITLE_FS, pad=8)
        ax.tick_params(direction='in', which='major', length=8, width=1.5,
                       labelsize=TICK_FS, top=True, right=True)
        ax.tick_params(direction='in', which='minor', length=5, width=1.0,
                       top=True, right=True)
        ax.minorticks_on()
        if hide_xticks:
            ax.tick_params(labelbottom=False)
        if hide_yticks:
            ax.tick_params(labelleft=False)
        ax.set_aspect('equal')

    # Per-channel color ranges (shared, set by VMIN_SHARED/VMAX_SHARED)
    R_vmin, R_vmax = get_clipped_range(R_log)
    G_vmin, G_vmax = get_clipped_range(G_log)
    B_vmin, B_vmax = get_clipped_range(B_log)

    # Top left: RGB composite (invisible cax keeps the panel size matched)
    ax = axes1[0, 0]
    ax.imshow(rgb_image_reg, extent=[x_min, x_max, y_min, y_max],
              origin='lower', aspect='equal', interpolation='bilinear')
    _style_panel(ax, r'$\rm RGB~(R=H\alpha,~G=MgII,~B=CIV)$', hide_xticks=True)
    # White ticks for the RGB composite so they remain visible over the
    # dark composite background.
    ax.tick_params(which='both', color='white')
    for spine in ax.spines.values():
        spine.set_edgecolor('white')
    # Dashed white circle marking the minimum allowed pillar radius.
    if pillar_rmin is not None:
        from matplotlib.patches import Circle
        ax.add_patch(Circle((0, 0), float(pillar_rmin),
                            fill=False, edgecolor='white',
                            linestyle='--', linewidth=1.5))
    _attach_cax(ax, None, None)

    # Top right: Halpha (shares y axis with RGB → hide y-tick labels)
    ax = axes1[0, 1]
    R_log_flat = R_log.flatten()
    R_log_reg = griddata((x_flat, y_flat), R_log_flat, (X_reg, Y_reg), method='linear', fill_value=np.nan)
    im1 = ax.contourf(X_reg, Y_reg, np.clip(R_log_reg, R_vmin, R_vmax),
                      levels=np.linspace(R_vmin, R_vmax, 50),
                      cmap='Reds', extend='neither',
                      vmin=R_vmin, vmax=R_vmax)
    _style_panel(ax, r'$\rm H\alpha$', hide_xticks=True, hide_yticks=True)
    _attach_cax(ax, im1, label_R)

    # Bottom left: MgII
    ax = axes1[1, 0]
    G_log_flat = G_log.flatten()
    G_log_reg = griddata((x_flat, y_flat), G_log_flat, (X_reg, Y_reg), method='linear', fill_value=np.nan)
    im2 = ax.contourf(X_reg, Y_reg, np.clip(G_log_reg, G_vmin, G_vmax),
                      levels=np.linspace(G_vmin, G_vmax, 50),
                      cmap='Greens', extend='neither',
                      vmin=G_vmin, vmax=G_vmax)
    _style_panel(ax, r'$\rm MgII$')
    _attach_cax(ax, im2, label_G)

    # Bottom right: CIV (shares y axis with MgII → hide y-tick labels)
    ax = axes1[1, 1]
    B_log_flat = B_log.flatten()
    B_log_reg = griddata((x_flat, y_flat), B_log_flat, (X_reg, Y_reg), method='linear', fill_value=np.nan)
    im3 = ax.contourf(X_reg, Y_reg, np.clip(B_log_reg, B_vmin, B_vmax),
                      levels=np.linspace(B_vmin, B_vmax, 50),
                      cmap='Blues', extend='neither',
                      vmin=B_vmin, vmax=B_vmax)
    _style_panel(ax, r'$\rm CIV$', hide_yticks=True)
    _attach_cax(ax, im3, label_B)

    # Shared X / Y labels at the figure edges; bring the X label up close
    # to the bottom panels to reduce vertical whitespace. Skip tight_layout
    # — it overrides the negative gridspec hspace.
    fig1.supxlabel(r'$\rm X~[light~days]$', fontsize=SUPLABEL_FS, y=0.10)
    fig1.supylabel(r'$\rm Y~[light~days]$', fontsize=SUPLABEL_FS, x=0.06)
    filename_xy = filename.replace('.png', '_xy.png')
    plt.savefig(filename_xy, dpi=200, bbox_inches='tight')
    print(f"RGB line intensity map (x-y plane) saved to {filename_xy}")
    plt.close()

    # Skip the remaining plots (r-phi RGB, logF-logphi, EW xy, EW r-phi)
    # while iterating on the xy composite. Flip _PLOT_EXTRAS=True above to
    # re-enable them all.
    if not _PLOT_EXTRAS:
        return

    # ====================================================
    # Plot 3: RGB + 3 channels in r-phi plane
    # ====================================================
    fig2, axes2 = plt.subplots(2, 2, figsize=(14, 14))
    
    # Top left: RGB composite in r-phi
    ax = axes2[0, 0]
    ax.imshow(rgb_image, extent=[phi_plot[0], phi_plot[-1], r_plot[0], r_plot[-1]], 
              aspect='auto', origin='lower', interpolation='bilinear')
    ax.set_xlabel(r'$\rm Azimuthal~Angle~\phi~[rad]$', fontsize=14)
    ax.set_ylabel(r'$\rm Radius~[light~days]$', fontsize=14)
    ax.set_title(r'$\rm RGB~Line~Intensity~Map~(R=H\alpha,~G=MgII,~B=CIV)$', fontsize=14)
    
    # Top right: Red channel r-phi - log scale absolute values (inverted colormap)
    ax = axes2[0, 1]
    im1b = ax.imshow(R_log, extent=[phi_plot[0], phi_plot[-1], r_plot[0], r_plot[-1]], 
                     aspect='auto', origin='lower', cmap='Reds', interpolation='bilinear')
    ax.set_xlabel(r'$\rm Azimuthal~Angle~\phi~[rad]$', fontsize=14)
    ax.set_ylabel(r'$\rm Radius~[light~days]$', fontsize=14)
    ax.set_title(r'$\rm Red~Channel~(H\alpha)$', fontsize=14)
    plt.colorbar(im1b, ax=ax, label=label_R, shrink=0.8, aspect=20)
    
    # Bottom left: Green channel r-phi - log scale absolute values (inverted colormap)
    ax = axes2[1, 0]
    im2b = ax.imshow(G_log, extent=[phi_plot[0], phi_plot[-1], r_plot[0], r_plot[-1]], 
                     aspect='auto', origin='lower', cmap='Greens', interpolation='bilinear')
    ax.set_xlabel(r'$\rm Azimuthal~Angle~\phi~[rad]$', fontsize=14)
    ax.set_ylabel(r'$\rm Radius~[light~days]$', fontsize=14)
    ax.set_title(r'$\rm Green~Channel~(MgII)$', fontsize=14)
    plt.colorbar(im2b, ax=ax, label=label_G, shrink=0.8, aspect=20)
    
    # Bottom right: Blue channel r-phi - log scale absolute values (inverted colormap)
    ax = axes2[1, 1]
    im3b = ax.imshow(B_log, extent=[phi_plot[0], phi_plot[-1], r_plot[0], r_plot[-1]], 
                     aspect='auto', origin='lower', cmap='Blues', interpolation='bilinear')
    ax.set_xlabel(r'$\rm Azimuthal~Angle~\phi~[rad]$', fontsize=14)
    ax.set_ylabel(r'$\rm Radius~[light~days]$', fontsize=14)
    ax.set_title(r'$\rm Blue~Channel~(CIV)$', fontsize=14)
    plt.colorbar(im3b, ax=ax, label=label_B, shrink=0.8, aspect=20)
    
    plt.tight_layout()
    filename_rphi = filename.replace('.png', '_rphi.png')
    plt.savefig(filename_rphi, dpi=150, bbox_inches='tight')
    print(f"RGB line intensity map (r-phi plane) saved to {filename_rphi}")
    plt.close()
    
    # ====================================================
    # Plot 4: logF - log phi maps in x-y plane + RGB composite (no r-phi panel)
    # ====================================================
    # Direct normalization check: (logF - logφ) for each of the RGB lines
    # We only plot x-y maps here (r-phi view is intentionally omitted).
    fig_logf, axes_logf = plt.subplots(2, 2, figsize=(14, 14))
    
    R_logF_minus_logphi = R_log - log_phi_2d_clipped
    G_logF_minus_logphi = G_log - log_phi_2d_clipped
    B_logF_minus_logphi = B_log - log_phi_2d_clipped
    
    # Interpolate to a regular x-y grid
    R_diff_reg = griddata((x_flat, y_flat), R_logF_minus_logphi.flatten(), (X_reg, Y_reg), method='linear', fill_value=np.nan)
    G_diff_reg = griddata((x_flat, y_flat), G_logF_minus_logphi.flatten(), (X_reg, Y_reg), method='linear', fill_value=np.nan)
    B_diff_reg = griddata((x_flat, y_flat), B_logF_minus_logphi.flatten(), (X_reg, Y_reg), method='linear', fill_value=np.nan)
    
    def normalize_linear(arr):
        """Normalize a 2D array to [0,1] ignoring NaNs."""
        finite = np.isfinite(arr)
        if not np.any(finite):
            return np.zeros_like(arr)
        vmin = np.nanmin(arr[finite])
        vmax = np.nanmax(arr[finite])
        if vmax <= vmin:
            return np.zeros_like(arr)
        out = (arr - vmin) / (vmax - vmin)
        out = np.clip(out, 0.0, 1.0)
        out[~finite] = 0.0
        return out
    
    # RGB composite of (logF - logφ)
    rgb_diff = np.zeros((X_reg.shape[0], X_reg.shape[1], 3))
    rgb_diff[:, :, 0] = normalize_linear(R_diff_reg)
    rgb_diff[:, :, 1] = normalize_linear(G_diff_reg)
    rgb_diff[:, :, 2] = normalize_linear(B_diff_reg)
    
    # Top left: RGB composite (logF - logφ) in x-y plane
    ax = axes_logf[0, 0]
    ax.imshow(rgb_diff, extent=[x_min, x_max, y_min, y_max], origin='lower', aspect='equal', interpolation='bilinear')
    ax.set_xlabel(r'$\rm X~[light~days]$', fontsize=14)
    ax.set_ylabel(r'$\rm Y~[light~days]$', fontsize=14)
    ax.set_title(r'$\rm RGB~Map~of~(\log F - \log\phi)~(R=H\alpha,~G=MgII,~B=CIV)$', fontsize=14)
    ax.set_aspect('equal')
    
    # Top right: Halpha (logF - logφ)
    ax = axes_logf[0, 1]
    im_R_diff = ax.contourf(X_reg, Y_reg, R_diff_reg, levels=50, cmap='Reds', extend='neither')
    ax.set_xlabel(r'$\rm X~[light~days]$', fontsize=14)
    ax.set_ylabel(r'$\rm Y~[light~days]$', fontsize=14)
    ax.set_title(r'$\rm \log(F_{H\alpha}) - \log\phi$', fontsize=14)
    plt.colorbar(im_R_diff, ax=ax, label=r'$\rm \log(F_{H\alpha}) - \log\phi$', shrink=0.8, aspect=20)
    ax.set_aspect('equal')
    
    # Bottom left: MgII (logF - logφ)
    ax = axes_logf[1, 0]
    im_G_diff = ax.contourf(X_reg, Y_reg, G_diff_reg, levels=50, cmap='Greens', extend='neither')
    ax.set_xlabel(r'$\rm X~[light~days]$', fontsize=14)
    ax.set_ylabel(r'$\rm Y~[light~days]$', fontsize=14)
    ax.set_title(r'$\rm \log(F_{MgII}) - \log\phi$', fontsize=14)
    plt.colorbar(im_G_diff, ax=ax, label=r'$\rm \log(F_{MgII}) - \log\phi$', shrink=0.8, aspect=20)
    ax.set_aspect('equal')
    
    # Bottom right: CIV (logF - logφ)
    ax = axes_logf[1, 1]
    im_B_diff = ax.contourf(X_reg, Y_reg, B_diff_reg, levels=50, cmap='Blues', extend='neither')
    ax.set_xlabel(r'$\rm X~[light~days]$', fontsize=14)
    ax.set_ylabel(r'$\rm Y~[light~days]$', fontsize=14)
    ax.set_title(r'$\rm \log(F_{CIV}) - \log\phi$', fontsize=14)
    plt.colorbar(im_B_diff, ax=ax, label=r'$\rm \log(F_{CIV}) - \log\phi$', shrink=0.8, aspect=20)
    ax.set_aspect('equal')
    
    plt.tight_layout()
    filename_logf = filename.replace('.png', '_logF_logphi.png')
    plt.savefig(filename_logf, dpi=150, bbox_inches='tight')
    print(f"logF - log phi (x-y + RGB) map saved to {filename_logf}")
    plt.close()
    
    # ====================================================
    # Plot 5: Equivalent Width maps (4-panel: RGB composite + 3 individual lines)
    # ====================================================
    if len(ew_map) > 0:
        # Only plot RGB lines: Halpha, MgII, CIV
        ew_lines_rgb = ['Halpha', 'Mg2', 'C4']
        available_ew_lines = [line for line in ew_lines_rgb if line in ew_map]
        
        if len(available_ew_lines) == 3:
            # Compute log scale EW values for RGB display
            ew_R_log, ew_R_log_min, ew_R_log_max = compute_log_values(ew_map['Halpha'])
            ew_G_log, ew_G_log_min, ew_G_log_max = compute_log_values(ew_map['Mg2'])
            ew_B_log, ew_B_log_min, ew_B_log_max = compute_log_values(ew_map['C4'])
            
            # Normalize for RGB display: high EW = deep/saturated color
            ew_R_norm = normalize_channel_log(ew_map['Halpha'])
            ew_G_norm = normalize_channel_log(ew_map['Mg2'])
            ew_B_norm = normalize_channel_log(ew_map['C4'])
            
            # Create RGB EW image
            ew_rgb_image = np.zeros((len(r_plot), len(phi_plot), 3))
            ew_rgb_image[:, :, 0] = ew_R_norm  # Red
            ew_rgb_image[:, :, 1] = ew_G_norm  # Green
            ew_rgb_image[:, :, 2] = ew_B_norm  # Blue
            
            # ====================================================
            # EW Plot 1: x-y plane (4-panel)
            # ====================================================
            fig_ew_xy, axes_ew_xy = plt.subplots(2, 2, figsize=(14, 14))
            
            # Top left: RGB composite EW in x-y plane
            ax = axes_ew_xy[0, 0]
            ew_rgb_r_flat = ew_rgb_image[:, :, 0].flatten()
            ew_rgb_g_flat = ew_rgb_image[:, :, 1].flatten()
            ew_rgb_b_flat = ew_rgb_image[:, :, 2].flatten()
            
            ew_rgb_r_reg = griddata((x_flat, y_flat), ew_rgb_r_flat, (X_reg, Y_reg), method='linear', fill_value=0)
            ew_rgb_g_reg = griddata((x_flat, y_flat), ew_rgb_g_flat, (X_reg, Y_reg), method='linear', fill_value=0)
            ew_rgb_b_reg = griddata((x_flat, y_flat), ew_rgb_b_flat, (X_reg, Y_reg), method='linear', fill_value=0)
            
            ew_rgb_image_reg = np.zeros((ew_rgb_r_reg.shape[0], ew_rgb_r_reg.shape[1], 3))
            ew_rgb_image_reg[:, :, 0] = np.clip(ew_rgb_r_reg, 0, 1)
            ew_rgb_image_reg[:, :, 1] = np.clip(ew_rgb_g_reg, 0, 1)
            ew_rgb_image_reg[:, :, 2] = np.clip(ew_rgb_b_reg, 0, 1)
            
            ax.imshow(ew_rgb_image_reg, extent=[x_min, x_max, y_min, y_max], 
                     origin='lower', aspect='equal', interpolation='bilinear')
            ax.set_xlabel(r'$\rm X~[light~days]$', fontsize=14)
            ax.set_ylabel(r'$\rm Y~[light~days]$', fontsize=14)
            ax.set_title(r'$\rm RGB~EW~Map~(R=H\alpha,~G=MgII,~B=CIV)$', fontsize=14)
            ax.set_aspect('equal')
            
            # Top right: Halpha EW in x-y plane (inverted colormap)
            ax = axes_ew_xy[0, 1]
            ew_R_log_flat = ew_R_log.flatten()
            ew_R_log_reg = griddata((x_flat, y_flat), ew_R_log_flat, (X_reg, Y_reg), method='linear', fill_value=np.nan)
            im_ew_R = ax.contourf(X_reg, Y_reg, ew_R_log_reg, levels=50, cmap='Reds', extend='neither')
            ax.set_xlabel(r'$\rm X~[light~days]$', fontsize=14)
            ax.set_ylabel(r'$\rm Y~[light~days]$', fontsize=14)
            ax.set_title(r'$\rm EW(H\alpha)$', fontsize=14)
            plt.colorbar(im_ew_R, ax=ax, label=r'$\rm \log(EW)$', shrink=0.8, aspect=20)
            ax.set_aspect('equal')
            
            # Bottom left: MgII EW in x-y plane (inverted colormap)
            ax = axes_ew_xy[1, 0]
            ew_G_log_flat = ew_G_log.flatten()
            ew_G_log_reg = griddata((x_flat, y_flat), ew_G_log_flat, (X_reg, Y_reg), method='linear', fill_value=np.nan)
            im_ew_G = ax.contourf(X_reg, Y_reg, ew_G_log_reg, levels=50, cmap='Greens', extend='neither')
            ax.set_xlabel(r'$\rm X~[light~days]$', fontsize=14)
            ax.set_ylabel(r'$\rm Y~[light~days]$', fontsize=14)
            ax.set_title(r'$\rm EW(MgII)$', fontsize=14)
            plt.colorbar(im_ew_G, ax=ax, label=r'$\rm \log(EW)$', shrink=0.8, aspect=20)
            ax.set_aspect('equal')
            
            # Bottom right: CIV EW in x-y plane (inverted colormap)
            ax = axes_ew_xy[1, 1]
            ew_B_log_flat = ew_B_log.flatten()
            ew_B_log_reg = griddata((x_flat, y_flat), ew_B_log_flat, (X_reg, Y_reg), method='linear', fill_value=np.nan)
            im_ew_B = ax.contourf(X_reg, Y_reg, ew_B_log_reg, levels=50, cmap='Blues', extend='neither')
            ax.set_xlabel(r'$\rm X~[light~days]$', fontsize=14)
            ax.set_ylabel(r'$\rm Y~[light~days]$', fontsize=14)
            ax.set_title(r'$\rm EW(CIV)$', fontsize=14)
            plt.colorbar(im_ew_B, ax=ax, label=r'$\rm \log(EW)$', shrink=0.8, aspect=20)
            ax.set_aspect('equal')
            
            plt.tight_layout()
            filename_ew_xy = filename.replace('.png', '_EW_xy.png')
            plt.savefig(filename_ew_xy, dpi=150, bbox_inches='tight')
            print(f"Equivalent width map (x-y plane) saved to {filename_ew_xy}")
            plt.close()
            
            # ====================================================
            # EW Plot 2: r-phi plane (4-panel)
            # ====================================================
            fig_ew_rphi, axes_ew_rphi = plt.subplots(2, 2, figsize=(14, 14))
            
            # Top left: RGB composite EW in r-phi plane
            ax = axes_ew_rphi[0, 0]
            ax.imshow(ew_rgb_image, extent=[phi_plot[0], phi_plot[-1], r_plot[0], r_plot[-1]], 
                     aspect='auto', origin='lower', interpolation='bilinear')
            ax.set_xlabel(r'$\rm Azimuthal~Angle~\phi~[rad]$', fontsize=14)
            ax.set_ylabel(r'$\rm Radius~[light~days]$', fontsize=14)
            ax.set_title(r'$\rm RGB~EW~Map~(R=H\alpha,~G=MgII,~B=CIV)$', fontsize=14)
            
            # Top right: Halpha EW in r-phi plane (inverted colormap)
            ax = axes_ew_rphi[0, 1]
            im_ew_R_rphi = ax.imshow(ew_R_log, extent=[phi_plot[0], phi_plot[-1], r_plot[0], r_plot[-1]], 
                                    aspect='auto', origin='lower', interpolation='bilinear', cmap='Reds')
            ax.set_xlabel(r'$\rm Azimuthal~Angle~\phi~[rad]$', fontsize=14)
            ax.set_ylabel(r'$\rm Radius~[light~days]$', fontsize=14)
            ax.set_title(r'$\rm EW(H\alpha)$', fontsize=14)
            plt.colorbar(im_ew_R_rphi, ax=ax, label=r'$\rm \log(EW)$', shrink=0.8, aspect=20)
            
            # Bottom left: MgII EW in r-phi plane (inverted colormap)
            ax = axes_ew_rphi[1, 0]
            im_ew_G_rphi = ax.imshow(ew_G_log, extent=[phi_plot[0], phi_plot[-1], r_plot[0], r_plot[-1]], 
                                    aspect='auto', origin='lower', interpolation='bilinear', cmap='Greens')
            ax.set_xlabel(r'$\rm Azimuthal~Angle~\phi~[rad]$', fontsize=14)
            ax.set_ylabel(r'$\rm Radius~[light~days]$', fontsize=14)
            ax.set_title(r'$\rm EW(MgII)$', fontsize=14)
            plt.colorbar(im_ew_G_rphi, ax=ax, label=r'$\rm \log(EW)$', shrink=0.8, aspect=20)
            
            # Bottom right: CIV EW in r-phi plane (inverted colormap)
            ax = axes_ew_rphi[1, 1]
            im_ew_B_rphi = ax.imshow(ew_B_log, extent=[phi_plot[0], phi_plot[-1], r_plot[0], r_plot[-1]], 
                                    aspect='auto', origin='lower', interpolation='bilinear', cmap='Blues')
            ax.set_xlabel(r'$\rm Azimuthal~Angle~\phi~[rad]$', fontsize=14)
            ax.set_ylabel(r'$\rm Radius~[light~days]$', fontsize=14)
            ax.set_title(r'$\rm EW(CIV)$', fontsize=14)
            plt.colorbar(im_ew_B_rphi, ax=ax, label=r'$\rm \log(EW)$', shrink=0.8, aspect=20)
            
            plt.tight_layout()
            filename_ew_rphi = filename.replace('.png', '_EW_rphi.png')
            plt.savefig(filename_ew_rphi, dpi=150, bbox_inches='tight')
            print(f"Equivalent width map (r-phi plane) saved to {filename_ew_rphi}")
            plt.close()
        else:
            print(f"Warning: Not all RGB EW lines available. Found: {available_ew_lines}")


def main(config_file='config_line.yaml', cloudy_file=None,
         cloudy_extension_file=None, use_absolute_flux=True,
         rgb_sat_mode='shared', hybrid_composite=False,
         composite_brightness='none', composite_sat_boost=1.0,
         composite_brightness_gamma=1.0,
         composite_brightness_floor=0.0, composite_brightness_boost=1.0,
         composite_brightness_phi_min=16.0, composite_brightness_phi_max=21.0,
         composite_brightness_add=0.0,
         draw_pillar_rmin=True):
    """
    Main function to compute and plot RGB line intensity map.

    Parameters:
    -----------
    config_file : str
        Path to configuration YAML file
    cloudy_file : str
        Path to Cloudy model data file. If None, use default path based on use_absolute_flux
    cloudy_extension_file : str or None
        Optional path to a low-phi (log Phi in [15, 16.5], Z = Z_sun)
        Cloudy extension file that gets spliced onto the main grid.
        Only meaningful when use_absolute_flux=True (the absolute-flux
        grid is the only one with a matching extension file). Default
        (None at the call) auto-fills with the standard extension file
        living next to the main absolute-flux file.
    use_absolute_flux : bool
        If True, use absolute line flux models.
        If False, use flux normalized by Hbeta.
    """
    # Load configuration
    config = load_config(config_file)

    if config is None:
        print("Error: Could not load config file. Using defaults.")
        return

    # Default Cloudy file path based on use_absolute_flux
    if cloudy_file is None:
        if use_absolute_flux:
            basepath = '/Users/jiamuh/c23.01/my_models/loc_metal_flux'
            cloudy_file = os.path.join(basepath, 'strong_LOC_varym_N25_v100_lineflux_LineList_BLR_Fe2_flux.txt')
            if cloudy_extension_file is None:
                cloudy_extension_file = os.path.join(
                    basepath,
                    'strong_LOC_varym_N25_v100_lineflux_extlow_LineList_BLR_Fe2_flux.txt')
        else:
            basepath = '/Users/jiamuh/c23.01/my_models/loc_metal'
            cloudy_file = os.path.join(basepath, 'strong_LOC_varym_N25_LineList_BLR_Fe2.txt')

    if not os.path.exists(cloudy_file):
        print(f"Error: Cloudy model file not found: {cloudy_file}")
        print("Please specify the path to the Cloudy model file.")
        return

    # Load Cloudy models (Z=1, solar metallicity), optionally with low-phi extension
    cloudy_interp_dict, phi_grid, Z_grid = load_cloudy_models(
        cloudy_file, Z_target=1.0, extension_file=cloudy_extension_file)
    
    # Extract disk parameters (same as pillar_line.py)
    disk_params = config.get('disk', {}).copy()
    temp_params = config.get('temperature', {}).copy()
    lamp_params = config.get('lamp', {}).copy()
    obs_params = config.get('observation', {}).copy()
    
    float_params = ['rin', 'rout', 'h1', 'r0', 'beta', 'tv1', 'alpha',
                   'tx1', 'tirrad_tvisc_ratio', 'fcol', 'hlamp', 'dmpc', 'cosi', 'inclination', 'redshift',
                   'M_BH', 'r_isco_rg', 'f_trans']
    int_params = ['nr', 'nphi']

    def convert_value(key, value):
        if value is None:
            return None
        if isinstance(value, str):
            if value.lower() == 'auto':
                return 'auto'
            try:
                if key in int_params:
                    return int(float(value))
                elif key in float_params:
                    return float(value)
            except (ValueError, TypeError):
                if value.lower() in ('true', '1', 'yes'):
                    return True
                elif value.lower() in ('false', '0', 'no'):
                    return False
                elif value.lower() in ('null', 'none'):
                    return None
        if key in int_params:
            return int(float(value))
        elif key in float_params:
            return float(value)
        return value

    for params_dict in [disk_params, temp_params, lamp_params, obs_params]:
        for key, value in params_dict.items():
            params_dict[key] = convert_value(key, value)

    # Convert inclination (degrees) to cosi if provided
    if 'inclination' in obs_params and obs_params['inclination'] is not None:
        inclination_rad = np.radians(obs_params['inclination'])
        obs_params['cosi'] = np.cos(inclination_rad)
        del obs_params['inclination']

    disk_params = {**disk_params, **temp_params, **lamp_params, **obs_params}

    # Resolve rin='auto' to ISCO if needed
    resolve_rin(disk_params)

    # Create disk model
    disk = PillarDisk(**disk_params)
    
    # Process pillars (same logic as pillar_line.py)
    pillars_config = config.get('pillars', {})
    make_many = pillars_config.get('make_many', False)
    if isinstance(make_many, str):
        make_many = make_many.lower() in ('true', '1', 'yes')
    make_many = bool(make_many)
    
    if make_many:
        print("Generating random pillars...")
        N_pillar = int(pillars_config.get('N_pillar', 100))
        r_mean = float(pillars_config.get('r_mean', 25.0))
        sig_r = float(pillars_config.get('sig_r', 10.0))
        rmin = pillars_config.get('rmin', None)
        rmin = float(rmin) if rmin is not None else None
        
        def get_pillar_param(key, default):
            val = pillars_config.get(key, default)
            if isinstance(val, list):
                return val
            return [val]
        
        h_pillar_list = get_pillar_param('h_pillar', 0.12)
        sigma_r_pillar_list = get_pillar_param('sigma_r_pillar', 2.0)
        sigma_phi_pillar_list = get_pillar_param('sigma_phi_pillar', 0.2)
        modify_height_pillar_list = get_pillar_param('modify_height_pillar', True)
        modify_temp_pillar_list = get_pillar_param('modify_temp_pillar', True)
        temp_factor_pillar_list = get_pillar_param('temp_factor_pillar', 1.5)
        
        def expand_list(lst, n, default):
            if len(lst) == 0:
                return [default] * n
            if len(lst) == 1:
                return lst * n
            while len(lst) < n:
                lst.append(lst[-1])
            return lst[:n]
        
        h_pillar_list = expand_list(h_pillar_list, N_pillar, 0.12)
        sigma_r_pillar_list = expand_list(sigma_r_pillar_list, N_pillar, 2.0)
        sigma_phi_pillar_list = expand_list(sigma_phi_pillar_list, N_pillar, 0.2)
        modify_height_pillar_list = expand_list(modify_height_pillar_list, N_pillar, True)
        modify_temp_pillar_list = expand_list(modify_temp_pillar_list, N_pillar, True)
        temp_factor_pillar_list = expand_list(temp_factor_pillar_list, N_pillar, 1.5)
        
        np.random.seed(42)
        r_pillar_random = np.random.normal(r_mean, sig_r, N_pillar)
        rmin_eff = max(disk.rin, rmin) if rmin is not None else disk.rin
        r_pillar_random = np.clip(r_pillar_random, rmin_eff, disk.rout)
        phi_pillar_random = np.random.uniform(0, 2*np.pi, N_pillar)
        
        def convert_bool(val):
            if isinstance(val, str):
                return val.lower() in ('true', '1', 'yes')
            return bool(val)
        
        expanded_pillars = []
        for i in range(N_pillar):
            pillar = {
                'r_pillar': float(r_pillar_random[i]),
                'phi_pillar': float(phi_pillar_random[i]),
                'height': float(h_pillar_list[i]),
                'sigma_r': float(sigma_r_pillar_list[i]),
                'sigma_phi': float(sigma_phi_pillar_list[i]),
                'temp_factor': float(temp_factor_pillar_list[i]),
                'modify_height': convert_bool(modify_height_pillar_list[i]),
                'modify_temp': convert_bool(modify_temp_pillar_list[i])
            }
            expanded_pillars.append(pillar)
        
        pillars_config = expanded_pillars
        print(f"Generated {N_pillar} random pillars")
    else:
        # Manual pillar placement
        if isinstance(pillars_config, dict):
            expanded_pillars = []
            
            def to_list(x):
                if isinstance(x, list):
                    return x
                return [x]
            
            r_pillar_list = to_list(pillars_config.get('r_pillar', []))
            phi_pillar_list = to_list(pillars_config.get('phi_pillar', []))
            height_list = to_list(pillars_config.get('height', [0.01]))
            sigma_r_list = to_list(pillars_config.get('sigma_r', [1.0]))
            sigma_phi_list = to_list(pillars_config.get('sigma_phi', [0.1]))
            temp_factor_list = to_list(pillars_config.get('temp_factor', [1.5]))
            modify_height_list = to_list(pillars_config.get('modify_height', [True]))
            modify_temp_list = to_list(pillars_config.get('modify_temp', [False]))
            
            n_pillars = max(len(r_pillar_list), len(phi_pillar_list), len(height_list),
                          len(sigma_r_list), len(sigma_phi_list), len(temp_factor_list),
                          len(modify_height_list), len(modify_temp_list))
            
            def expand_list(lst, n, default):
                if len(lst) == 0:
                    return [default] * n
                if len(lst) == 1:
                    return lst * n
                while len(lst) < n:
                    lst.append(lst[-1])
                return lst[:n]
            
            r_pillar_list = expand_list(r_pillar_list, n_pillars, 20.0)
            phi_pillar_list = expand_list(phi_pillar_list, n_pillars, 0.0)
            height_list = expand_list(height_list, n_pillars, 0.01)
            sigma_r_list = expand_list(sigma_r_list, n_pillars, 1.0)
            sigma_phi_list = expand_list(sigma_phi_list, n_pillars, 0.1)
            temp_factor_list = expand_list(temp_factor_list, n_pillars, 1.5)
            modify_height_list = expand_list(modify_height_list, n_pillars, True)
            modify_temp_list = expand_list(modify_temp_list, n_pillars, False)
            
            for i in range(n_pillars):
                pillar = {
                    'r_pillar': float(r_pillar_list[i]),
                    'phi_pillar': parse_math_expr(phi_pillar_list[i]),
                    'height': float(height_list[i]),
                    'sigma_r': float(sigma_r_list[i]),
                    'sigma_phi': parse_math_expr(sigma_phi_list[i]),
                    'temp_factor': float(temp_factor_list[i]),
                    'modify_height': bool(modify_height_list[i]) if not isinstance(modify_height_list[i], str) 
                                   else modify_height_list[i].lower() in ('true', '1', 'yes'),
                    'modify_temp': bool(modify_temp_list[i]) if not isinstance(modify_temp_list[i], str)
                                 else modify_temp_list[i].lower() in ('true', '1', 'yes')
                }
                expanded_pillars.append(pillar)
            
            pillars_config = expanded_pillars
    
    # Add pillars
    for pillar_cfg in pillars_config:
        disk.add_pillar(**pillar_cfg)
    
    print(f"\nNumber of pillars: {len(disk.pillars)}")
    
    # Compute line intensity map
    plot_params = config.get('plotting', {})
    n_r_plot = int(plot_params.get('n_r_plot', 200))
    n_phi_plot = int(plot_params.get('n_phi_plot', 360))
    
    r_plot, phi_plot, intensity_map, ew_map, log_phi_2d = compute_line_intensity_map(
        disk, cloudy_interp_dict, phi_grid, Z_target=1.0,
        n_r_plot=n_r_plot, n_phi_plot=n_phi_plot)

    # In hybrid_composite mode, load the OTHER cloudy grid so the RGB
    # composite uses absolute flux while the per-line panels stay on
    # F/F_Hβ (or vice versa). Detect the partner grid from the absolute
    # vs normalized basepath.
    composite_intensity_map = None
    if hybrid_composite:
        if use_absolute_flux:
            partner_basepath = '/Users/jiamuh/c23.01/my_models/loc_metal'
            partner_file = os.path.join(partner_basepath, 'strong_LOC_varym_N25_LineList_BLR_Fe2.txt')
            partner_ext = None
        else:
            partner_basepath = '/Users/jiamuh/c23.01/my_models/loc_metal_flux'
            partner_file = os.path.join(partner_basepath, 'strong_LOC_varym_N25_v100_lineflux_LineList_BLR_Fe2_flux.txt')
            partner_ext = os.path.join(partner_basepath,
                                       'strong_LOC_varym_N25_v100_lineflux_extlow_LineList_BLR_Fe2_flux.txt')
        if os.path.exists(partner_file):
            partner_interp, partner_phi_grid, _ = load_cloudy_models(
                partner_file, Z_target=1.0, extension_file=partner_ext)
            _, _, composite_intensity_map, _, _ = compute_line_intensity_map(
                disk, partner_interp, partner_phi_grid, Z_target=1.0,
                n_r_plot=n_r_plot, n_phi_plot=n_phi_plot)
        else:
            print(f"Warning: hybrid_composite requested but partner grid not found: {partner_file}")

    # Plot RGB map; save under plots/ rather than the project root.
    filename = plot_params.get('line_intensity_rgb_filename',
                               'plots/line_intensity_rgb.png')
    out_dir = os.path.dirname(filename)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    pillar_rmin_arg = None
    if draw_pillar_rmin:
        _rmin_cfg = config.get('pillars', {}).get('rmin', None)
        if _rmin_cfg is not None:
            pillar_rmin_arg = float(_rmin_cfg)
    plot_rgb_line_map(r_plot, phi_plot, intensity_map, ew_map, log_phi_2d, phi_grid, use_absolute_flux=use_absolute_flux, filename=filename, rgb_sat_mode=rgb_sat_mode, composite_intensity_map=composite_intensity_map, composite_brightness=composite_brightness, composite_sat_boost=composite_sat_boost, composite_brightness_gamma=composite_brightness_gamma, composite_brightness_floor=composite_brightness_floor, composite_brightness_boost=composite_brightness_boost, composite_brightness_phi_min=composite_brightness_phi_min, composite_brightness_phi_max=composite_brightness_phi_max, composite_brightness_add=composite_brightness_add, pillar_rmin=pillar_rmin_arg)
    
    print("\nSummary:")
    print(f"Radial range: {r_plot[0]:.2f} to {r_plot[-1]:.2f} light days")
    print(f"Azimuthal range: {phi_plot[0]:.2f} to {phi_plot[-1]:.2f} rad")
    print(f"Available lines: {list(intensity_map.keys())}")


if __name__ == '__main__':
    import sys
    
    # Simple argument parsing
    config_file = sys.argv[1] if len(sys.argv) > 1 else 'config_line.yaml'
    cloudy_file = sys.argv[2] if len(sys.argv) > 2 else None
    
    # Check for use_absolute_flux flag
    use_absolute_flux = False  # default
    if len(sys.argv) > 3:
        flag = sys.argv[3].lower()
        if flag in ('false', '0', 'no', '--use_normalized_flux'):
            use_absolute_flux = False
        elif flag in ('true', '1', 'yes', '--use_absolute_flux'):
            use_absolute_flux = True
    
    main(config_file=config_file, cloudy_file=cloudy_file, use_absolute_flux=use_absolute_flux)
