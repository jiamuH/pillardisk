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
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    from pillar_disk import PillarDisk, load_config
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


def load_cloudy_models(model_file, Z_target=1.0):
    """
    Load Cloudy model data and create interpolation functions for solar metallicity (Z=1).
    
    Parameters:
    -----------
    model_file : str
        Path to Cloudy model data file
    Z_target : float
        Target metallicity (default: 1.0 for solar)
    
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
    
    # Grid definitions (from line_ratio_breathing_effect.py)
    Z_grid = np.arange(1, 15.5, 0.5)          # 29 metallicity values
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
    
    # Load data
    df = pd.read_csv(model_file, sep='\t', header=0, comment='#')
    assert len(df) == len(Z_grid) * len(phi_grid) * len(n_grid), "Row count mismatch!"
    
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
    
    # Create 1D interpolation functions over phi (for Z=1)
    interp_dict = {}
    for key in line_names:
        interp_dict[key] = RectBivariateSpline(phi_grid, Z_grid, collapsed_data[key], kx=2, ky=2)
    
    print(f"Loaded Cloudy models: {list(line_names.keys())}")
    print(f"  phi range: {phi_grid[0]:.1f} to {phi_grid[-1]:.1f} (log)")
    print(f"  Z range: {Z_grid[0]:.1f} to {Z_grid[-1]:.1f}")
    
    return interp_dict, phi_grid, Z_grid


def compute_ionizing_flux(disk, r, phi):
    """
    Compute ionizing flux at disk positions.
    
    The ionizing flux is proportional to the irradiation temperature.
    Highest flux is at illuminated pillars, lowest in shadows.
    
    Parameters:
    -----------
    disk : PillarDisk
        Disk model
    r : array
        Radial positions (light days), 2D
    phi : array
        Azimuthal positions (radians), 2D
    
    Returns:
    --------
    log_phi : array
        Log ionizing flux (2D, same shape as r)
    """
    # Get temperature (with shadows)
    T_2d = disk.get_temperature(r, phi, compute_shadows=True)
    
    # Get base irradiation temperature (without pillars, for normalization)
    # This gives us the baseline flux
    T_base_2d = disk.get_temperature(r, phi, compute_shadows=False)
    
    # Ionizing flux is proportional to T^4 (Stefan-Boltzmann)
    # But we want to map to Cloudy's phi grid (log phi from 17 to 21.5)
    # Use the ratio of T to reference temperature to scale flux
    
    # Reference: at r0, T_irrad = tv1 * tirrad_tvisc_ratio
    T_ref = disk.tv1 * disk.tirrad_tvisc_ratio
    
    # Normalize by reference temperature
    T_ratio = T_2d / T_ref
    
    # Map to log phi range [17, 21.5]
    # Use T^4 scaling for flux, then map to log space
    # phi ∝ T^4, so log phi = log phi_min + 4 * log(T/T_min)
    # But we want to map T_ratio to the phi grid range
    
    # Simple mapping: log_phi = phi_min + (phi_max - phi_min) * (T_ratio^4 - T_min^4) / (T_max^4 - T_min^4)
    # But T_ratio can be < 1 in shadows, so we need to handle that
    
    # Better: map T_ratio to [0, 1], then scale to [phi_min, phi_max]
    T_ratio_clipped = np.clip(T_ratio, 0.1, 10.0)  # Clip to reasonable range
    T_ratio_normalized = (T_ratio_clipped - 0.1) / (10.0 - 0.1)  # Normalize to [0, 1]
    
    phi_min = 17.0
    phi_max = 21.5
    log_phi = phi_min + T_ratio_normalized * (phi_max - phi_min)
    
    # In shadows, T is very low, so set log_phi to minimum
    shadow_mask = T_2d < (T_ref * 0.1)
    log_phi[shadow_mask] = phi_min
    
    return log_phi


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
        for ir in range(n_r_plot):
            for iphi in range(n_phi_plot):
                log_phi_val = log_phi_2d[ir, iphi]
                # Clip to grid range
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


def plot_rgb_line_map(r_plot, phi_plot, intensity_map, ew_map, log_phi_2d, phi_grid, use_absolute_flux=True, filename='line_intensity_rgb.png'):
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
    
    # Use correct lines for RGB
    R = intensity_map['Halpha']  # Red channel: Halpha
    G = intensity_map['Mg2']     # Green channel: MgII
    B = intensity_map['C4']      # Blue channel: CIV
    
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
    
    def normalize_channel_log(channel):
        """Normalize channel to [0,1] for RGB display using log scale."""
        log_channel, log_min, log_max = compute_log_values(channel)
        if np.isnan(log_min):
            return np.zeros_like(channel)
        if log_max > log_min:
            normalized = (log_channel - log_min) / (log_max - log_min)
            normalized = np.clip(normalized, 0, 1)
        else:
            normalized = np.zeros_like(channel)
        return normalized
    
    # For RGB display (normalized, inverted so high values = bright colors)
    R_norm = 1.0 - normalize_channel_log(R)  # Invert: high flux = bright
    G_norm = 1.0 - normalize_channel_log(G)  # Invert: high flux = bright
    B_norm = 1.0 - normalize_channel_log(B)  # Invert: high flux = bright
    
    # For plotting absolute values (log scale)
    R_log, R_log_min, R_log_max = compute_log_values(R)
    G_log, G_log_min, G_log_max = compute_log_values(G)
    B_log, B_log_min, B_log_max = compute_log_values(B)
    
    # Create RGB image
    rgb_image = np.zeros((len(r_plot), len(phi_plot), 3))
    rgb_image[:, :, 0] = R_norm  # Red
    rgb_image[:, :, 1] = G_norm  # Green
    rgb_image[:, :, 2] = B_norm  # Blue
    
    # Convert to Cartesian coordinates (x-y plane) for plotting
    # Create meshgrid in polar coordinates, then convert to Cartesian
    R_mesh, Phi_mesh = np.meshgrid(r_plot, phi_plot, indexing='ij')
    X = R_mesh * np.cos(Phi_mesh)  # x coordinates (light days)
    Y = R_mesh * np.sin(Phi_mesh)  # y coordinates (light days)
    
    # Determine plot extent
    x_min, x_max = X.min(), X.max()
    y_min, y_max = Y.min(), Y.max()
    
    # Clip log_phi to minimum of model grid (17.0)
    phi_min_grid = phi_grid[0]  # Should be 17.0
    log_phi_2d_clipped = np.clip(log_phi_2d, phi_min_grid, None)
    
    # Prepare interpolation grid (used for both plots)
    x_flat = X.flatten()
    y_flat = Y.flatten()
    x_reg = np.linspace(x_min, x_max, len(r_plot))
    y_reg = np.linspace(y_min, y_max, len(r_plot))
    X_reg, Y_reg = np.meshgrid(x_reg, y_reg)
    
    # Colorbar labels based on use_absolute_flux
    # When use_absolute_flux=False, Cloudy output already has F_Hbeta=1, so F_Halpha/F_Hbeta = F_Halpha
    if use_absolute_flux:
        label_R = r'$\rm \log(F_{H\alpha})$'
        label_G = r'$\rm \log(F_{MgII})$'
        label_B = r'$\rm \log(F_{CIV})$'
    else:
        # When normalized, F_Hbeta=1, so F_Halpha/F_Hbeta = F_Halpha (actual value)
        label_R = r'$\rm \log(F_{H\alpha})$'
        label_G = r'$\rm \log(F_{MgII})$'
        label_B = r'$\rm \log(F_{CIV})$'
    
    # Separate into three plots:
    # 1. Ionizing flux (separate plot)
    # 2. RGB + 3 channels in x-y plane
    # 3. RGB + 3 channels in r-phi plane
    
    # ====================================================
    # Plot 1: Ionizing flux (separate plot)
    # ====================================================
    fig_phi, axes_phi = plt.subplots(1, 2, figsize=(14, 6))
    
    # Left: x-y plane
    ax = axes_phi[0]
    log_phi_flat_clipped = log_phi_2d_clipped.flatten()
    log_phi_reg = griddata((x_flat, y_flat), log_phi_flat_clipped, (X_reg, Y_reg), method='linear', fill_value=phi_min_grid)
    im_phi_xy = ax.contourf(X_reg, Y_reg, log_phi_reg, levels=50, cmap='viridis', extend='both', vmin=phi_min_grid)
    ax.set_xlabel(r'$\rm X~[light~days]$', fontsize=14)
    ax.set_ylabel(r'$\rm Y~[light~days]$', fontsize=14)
    ax.set_title(r'$\rm Log~Ionizing~Flux~\log\phi~(x-y~plane)$', fontsize=14)
    plt.colorbar(im_phi_xy, ax=ax, label=r'$\rm \log\phi$')
    ax.set_aspect('equal')
    
    # Right: r-phi plane
    ax = axes_phi[1]
    im_phi_rphi = ax.imshow(log_phi_2d_clipped, extent=[phi_plot[0], phi_plot[-1], r_plot[0], r_plot[-1]], 
                           aspect='auto', origin='lower', interpolation='bilinear', cmap='viridis', vmin=phi_min_grid)
    ax.set_xlabel(r'$\rm Azimuthal~Angle~\phi~[rad]$', fontsize=14)
    ax.set_ylabel(r'$\rm Radius~[light~days]$', fontsize=14)
    ax.set_title(r'$\rm Log~Ionizing~Flux~\log\phi~(r-\phi~plane)$', fontsize=14)
    plt.colorbar(im_phi_rphi, ax=ax, label=r'$\rm \log\phi$')
    
    plt.tight_layout()
    filename_phi = filename.replace('.png', '_ionizing_flux.png')
    plt.savefig(filename_phi, dpi=150, bbox_inches='tight')
    print(f"Ionizing flux map saved to {filename_phi}")
    plt.close()
    
    # ====================================================
    # Plot 2: RGB + 3 channels in x-y plane
    # ====================================================
    fig1, axes1 = plt.subplots(2, 2, figsize=(14, 14))
    
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
    
    # Top left: RGB composite in x-y plane
    ax = axes1[0, 0]
    ax.imshow(rgb_image_reg, extent=[x_min, x_max, y_min, y_max], 
              origin='lower', aspect='equal', interpolation='bilinear')
    ax.set_xlabel(r'$\rm X~[light~days]$', fontsize=14)
    ax.set_ylabel(r'$\rm Y~[light~days]$', fontsize=14)
    ax.set_title(r'$\rm RGB~Line~Intensity~Map~(R=H\alpha,~G=MgII,~B=CIV)$', fontsize=14)
    ax.set_aspect('equal')
    
    # Top right: Red channel (Halpha) in x-y plane - log scale absolute values (inverted colormap)
    ax = axes1[0, 1]
    R_log_flat = R_log.flatten()
    R_log_reg = griddata((x_flat, y_flat), R_log_flat, (X_reg, Y_reg), method='linear', fill_value=np.nan)
    im1 = ax.contourf(X_reg, Y_reg, R_log_reg, levels=50, cmap='Reds_r', extend='both')
    ax.set_xlabel(r'$\rm X~[light~days]$', fontsize=14)
    ax.set_ylabel(r'$\rm Y~[light~days]$', fontsize=14)
    ax.set_title(r'$\rm Red~Channel~(H\alpha)$', fontsize=14)
    plt.colorbar(im1, ax=ax, label=label_R)
    ax.set_aspect('equal')
    
    # Bottom left: Green channel (MgII) in x-y plane - log scale absolute values (inverted colormap)
    ax = axes1[1, 0]
    G_log_flat = G_log.flatten()
    G_log_reg = griddata((x_flat, y_flat), G_log_flat, (X_reg, Y_reg), method='linear', fill_value=np.nan)
    im2 = ax.contourf(X_reg, Y_reg, G_log_reg, levels=50, cmap='Greens_r', extend='both')
    ax.set_xlabel(r'$\rm X~[light~days]$', fontsize=14)
    ax.set_ylabel(r'$\rm Y~[light~days]$', fontsize=14)
    ax.set_title(r'$\rm Green~Channel~(MgII)$', fontsize=14)
    plt.colorbar(im2, ax=ax, label=label_G)
    ax.set_aspect('equal')
    
    # Bottom right: Blue channel (CIV) in x-y plane - log scale absolute values (inverted colormap)
    ax = axes1[1, 1]
    B_log_flat = B_log.flatten()
    B_log_reg = griddata((x_flat, y_flat), B_log_flat, (X_reg, Y_reg), method='linear', fill_value=np.nan)
    im3 = ax.contourf(X_reg, Y_reg, B_log_reg, levels=50, cmap='Blues_r', extend='both')
    ax.set_xlabel(r'$\rm X~[light~days]$', fontsize=14)
    ax.set_ylabel(r'$\rm Y~[light~days]$', fontsize=14)
    ax.set_title(r'$\rm Blue~Channel~(CIV)$', fontsize=14)
    plt.colorbar(im3, ax=ax, label=label_B)
    ax.set_aspect('equal')
    
    plt.tight_layout()
    filename_xy = filename.replace('.png', '_xy.png')
    plt.savefig(filename_xy, dpi=150, bbox_inches='tight')
    print(f"RGB line intensity map (x-y plane) saved to {filename_xy}")
    plt.close()
    
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
                     aspect='auto', origin='lower', cmap='Reds_r', interpolation='bilinear')
    ax.set_xlabel(r'$\rm Azimuthal~Angle~\phi~[rad]$', fontsize=14)
    ax.set_ylabel(r'$\rm Radius~[light~days]$', fontsize=14)
    ax.set_title(r'$\rm Red~Channel~(H\alpha)$', fontsize=14)
    plt.colorbar(im1b, ax=ax, label=label_R)
    
    # Bottom left: Green channel r-phi - log scale absolute values (inverted colormap)
    ax = axes2[1, 0]
    im2b = ax.imshow(G_log, extent=[phi_plot[0], phi_plot[-1], r_plot[0], r_plot[-1]], 
                     aspect='auto', origin='lower', cmap='Greens_r', interpolation='bilinear')
    ax.set_xlabel(r'$\rm Azimuthal~Angle~\phi~[rad]$', fontsize=14)
    ax.set_ylabel(r'$\rm Radius~[light~days]$', fontsize=14)
    ax.set_title(r'$\rm Green~Channel~(MgII)$', fontsize=14)
    plt.colorbar(im2b, ax=ax, label=label_G)
    
    # Bottom right: Blue channel r-phi - log scale absolute values (inverted colormap)
    ax = axes2[1, 1]
    im3b = ax.imshow(B_log, extent=[phi_plot[0], phi_plot[-1], r_plot[0], r_plot[-1]], 
                     aspect='auto', origin='lower', cmap='Blues_r', interpolation='bilinear')
    ax.set_xlabel(r'$\rm Azimuthal~Angle~\phi~[rad]$', fontsize=14)
    ax.set_ylabel(r'$\rm Radius~[light~days]$', fontsize=14)
    ax.set_title(r'$\rm Blue~Channel~(CIV)$', fontsize=14)
    plt.colorbar(im3b, ax=ax, label=label_B)
    
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
    im_R_diff = ax.contourf(X_reg, Y_reg, R_diff_reg, levels=50, cmap='Reds_r', extend='both')
    ax.set_xlabel(r'$\rm X~[light~days]$', fontsize=14)
    ax.set_ylabel(r'$\rm Y~[light~days]$', fontsize=14)
    ax.set_title(r'$\rm \log(F_{H\alpha}) - \log\phi$', fontsize=14)
    plt.colorbar(im_R_diff, ax=ax, label=r'$\rm \log(F_{H\alpha}) - \log\phi$')
    ax.set_aspect('equal')
    
    # Bottom left: MgII (logF - logφ)
    ax = axes_logf[1, 0]
    im_G_diff = ax.contourf(X_reg, Y_reg, G_diff_reg, levels=50, cmap='Greens_r', extend='both')
    ax.set_xlabel(r'$\rm X~[light~days]$', fontsize=14)
    ax.set_ylabel(r'$\rm Y~[light~days]$', fontsize=14)
    ax.set_title(r'$\rm \log(F_{MgII}) - \log\phi$', fontsize=14)
    plt.colorbar(im_G_diff, ax=ax, label=r'$\rm \log(F_{MgII}) - \log\phi$')
    ax.set_aspect('equal')
    
    # Bottom right: CIV (logF - logφ)
    ax = axes_logf[1, 1]
    im_B_diff = ax.contourf(X_reg, Y_reg, B_diff_reg, levels=50, cmap='Blues_r', extend='both')
    ax.set_xlabel(r'$\rm X~[light~days]$', fontsize=14)
    ax.set_ylabel(r'$\rm Y~[light~days]$', fontsize=14)
    ax.set_title(r'$\rm \log(F_{CIV}) - \log\phi$', fontsize=14)
    plt.colorbar(im_B_diff, ax=ax, label=r'$\rm \log(F_{CIV}) - \log\phi$')
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
            
            # Normalize for RGB display (invert so high EW = bright colors)
            ew_R_norm = 1.0 - normalize_channel_log(ew_map['Halpha'])  # Invert normalization
            ew_G_norm = 1.0 - normalize_channel_log(ew_map['Mg2'])  # Invert normalization
            ew_B_norm = 1.0 - normalize_channel_log(ew_map['C4'])  # Invert normalization
            
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
            im_ew_R = ax.contourf(X_reg, Y_reg, ew_R_log_reg, levels=50, cmap='Reds_r', extend='both')
            ax.set_xlabel(r'$\rm X~[light~days]$', fontsize=14)
            ax.set_ylabel(r'$\rm Y~[light~days]$', fontsize=14)
            ax.set_title(r'$\rm EW(H\alpha)$', fontsize=14)
            plt.colorbar(im_ew_R, ax=ax, label=r'$\rm \log(EW)$')
            ax.set_aspect('equal')
            
            # Bottom left: MgII EW in x-y plane (inverted colormap)
            ax = axes_ew_xy[1, 0]
            ew_G_log_flat = ew_G_log.flatten()
            ew_G_log_reg = griddata((x_flat, y_flat), ew_G_log_flat, (X_reg, Y_reg), method='linear', fill_value=np.nan)
            im_ew_G = ax.contourf(X_reg, Y_reg, ew_G_log_reg, levels=50, cmap='Greens_r', extend='both')
            ax.set_xlabel(r'$\rm X~[light~days]$', fontsize=14)
            ax.set_ylabel(r'$\rm Y~[light~days]$', fontsize=14)
            ax.set_title(r'$\rm EW(MgII)$', fontsize=14)
            plt.colorbar(im_ew_G, ax=ax, label=r'$\rm \log(EW)$')
            ax.set_aspect('equal')
            
            # Bottom right: CIV EW in x-y plane (inverted colormap)
            ax = axes_ew_xy[1, 1]
            ew_B_log_flat = ew_B_log.flatten()
            ew_B_log_reg = griddata((x_flat, y_flat), ew_B_log_flat, (X_reg, Y_reg), method='linear', fill_value=np.nan)
            im_ew_B = ax.contourf(X_reg, Y_reg, ew_B_log_reg, levels=50, cmap='Blues_r', extend='both')
            ax.set_xlabel(r'$\rm X~[light~days]$', fontsize=14)
            ax.set_ylabel(r'$\rm Y~[light~days]$', fontsize=14)
            ax.set_title(r'$\rm EW(CIV)$', fontsize=14)
            plt.colorbar(im_ew_B, ax=ax, label=r'$\rm \log(EW)$')
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
                                    aspect='auto', origin='lower', interpolation='bilinear', cmap='Reds_r')
            ax.set_xlabel(r'$\rm Azimuthal~Angle~\phi~[rad]$', fontsize=14)
            ax.set_ylabel(r'$\rm Radius~[light~days]$', fontsize=14)
            ax.set_title(r'$\rm EW(H\alpha)$', fontsize=14)
            plt.colorbar(im_ew_R_rphi, ax=ax, label=r'$\rm \log(EW)$')
            
            # Bottom left: MgII EW in r-phi plane (inverted colormap)
            ax = axes_ew_rphi[1, 0]
            im_ew_G_rphi = ax.imshow(ew_G_log, extent=[phi_plot[0], phi_plot[-1], r_plot[0], r_plot[-1]], 
                                    aspect='auto', origin='lower', interpolation='bilinear', cmap='Greens_r')
            ax.set_xlabel(r'$\rm Azimuthal~Angle~\phi~[rad]$', fontsize=14)
            ax.set_ylabel(r'$\rm Radius~[light~days]$', fontsize=14)
            ax.set_title(r'$\rm EW(MgII)$', fontsize=14)
            plt.colorbar(im_ew_G_rphi, ax=ax, label=r'$\rm \log(EW)$')
            
            # Bottom right: CIV EW in r-phi plane (inverted colormap)
            ax = axes_ew_rphi[1, 1]
            im_ew_B_rphi = ax.imshow(ew_B_log, extent=[phi_plot[0], phi_plot[-1], r_plot[0], r_plot[-1]], 
                                    aspect='auto', origin='lower', interpolation='bilinear', cmap='Blues_r')
            ax.set_xlabel(r'$\rm Azimuthal~Angle~\phi~[rad]$', fontsize=14)
            ax.set_ylabel(r'$\rm Radius~[light~days]$', fontsize=14)
            ax.set_title(r'$\rm EW(CIV)$', fontsize=14)
            plt.colorbar(im_ew_B_rphi, ax=ax, label=r'$\rm \log(EW)$')
            
            plt.tight_layout()
            filename_ew_rphi = filename.replace('.png', '_EW_rphi.png')
            plt.savefig(filename_ew_rphi, dpi=150, bbox_inches='tight')
            print(f"Equivalent width map (r-phi plane) saved to {filename_ew_rphi}")
            plt.close()
        else:
            print(f"Warning: Not all RGB EW lines available. Found: {available_ew_lines}")


def main(config_file='config_line.yaml', cloudy_file=None, use_absolute_flux=True):
    """
    Main function to compute and plot RGB line intensity map.
    
    Parameters:
    -----------
    config_file : str
        Path to configuration YAML file
    cloudy_file : str
        Path to Cloudy model data file. If None, use default path based on use_absolute_flux
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
        else:
            basepath = '/Users/jiamuh/c23.01/my_models/loc_metal'
            cloudy_file = os.path.join(basepath, 'strong_LOC_varym_N25_LineList_BLR_Fe2.txt')
    
    if not os.path.exists(cloudy_file):
        print(f"Error: Cloudy model file not found: {cloudy_file}")
        print("Please specify the path to the Cloudy model file.")
        return
    
    # Load Cloudy models (Z=1, solar metallicity)
    cloudy_interp_dict, phi_grid, Z_grid = load_cloudy_models(cloudy_file, Z_target=1.0)
    
    # Extract disk parameters (same as pillar_line.py)
    disk_params = config.get('disk', {}).copy()
    temp_params = config.get('temperature', {}).copy()
    lamp_params = config.get('lamp', {}).copy()
    obs_params = config.get('observation', {}).copy()
    
    float_params = ['rin', 'rout', 'h1', 'r0', 'beta', 'tv1', 'alpha', 
                   'tx1', 'tirrad_tvisc_ratio', 'fcol', 'hlamp', 'dmpc', 'cosi', 'redshift']
    int_params = ['nr', 'nphi']
    
    def convert_value(key, value):
        if value is None:
            return None
        if isinstance(value, str):
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
    
    disk_params = {**disk_params, **temp_params, **lamp_params, **obs_params}
    
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
    
    # Plot RGB map
    filename = plot_params.get('line_intensity_rgb_filename', 'line_intensity_rgb.png')
    plot_rgb_line_map(r_plot, phi_plot, intensity_map, ew_map, log_phi_2d, phi_grid, use_absolute_flux=use_absolute_flux, filename=filename)
    
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
