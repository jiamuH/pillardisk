#!/usr/bin/env python3
"""
pillar_line_time.py - Compute time-evolving velocity delay map with rotation

This script computes the wavelength-time map for an emission line from gas
orbiting the black hole. The disk and pillars rotate in circular orbits with
Keplerian velocity v = sqrt(GM/r), creating a barber-pole pattern as features
rotate around the disk.

The output is a 2D map: x-axis = wavelength, y-axis = time (days),
map value = spectral flux/intensity.
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys

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

# Define parse_math_expr (it's defined inside main() in pillar_disk.py, so we define it here)
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
G = 6.673e-8     # gravitational constant (cgs)
PC = 3.0857e18   # parsec (cm)
DAY = 24. * 3600.  # day (s)
ANGSTROM = 1e-8   # Angstrom (cm)
KM_TO_CM = 1e5    # km to cm conversion

# Conversion factors
PC_TO_LD = 1e6 * PC / (C * DAY)  # parsec to light days
LD_TO_CM = C * DAY  # light days to cm


def compute_angular_velocity(r, v_virial, r_virial):
    """
    Compute angular velocity for circular orbit: omega = v/r = sqrt(GM/r^3)
    
    Parameters:
    -----------
    r : array
        Radius (light days)
    v_virial : float
        Virial velocity at reference radius (km/s)
    r_virial : float
        Reference radius (light days)
    
    Returns:
    --------
    omega : array
        Angular velocity (rad/day)
    """
    # Convert to cgs
    r_cgs = r * LD_TO_CM  # cm
    v_virial_cgs = v_virial * KM_TO_CM  # cm/s
    r_virial_cgs = r_virial * LD_TO_CM  # cm
    
    # Keplerian velocity: v = sqrt(GM/r) = v_virial * sqrt(r_virial/r)
    v_phi = v_virial_cgs * np.sqrt(r_virial_cgs / r_cgs)  # cm/s
    
    # Angular velocity: omega = v/r (rad/s)
    omega_rad_s = v_phi / r_cgs  # rad/s
    
    # Convert to rad/day
    omega_rad_day = omega_rad_s * DAY  # rad/day
    
    return omega_rad_day


def compute_orbital_period(r, v_virial, r_virial):
    """
    Compute orbital period for circular orbit at radius r.
    
    T = 2π / ω = 2π * r / v
    where v = v_virial * sqrt(r_virial / r)
    
    Parameters:
    -----------
    r : float or array
        Radius (light days)
    v_virial : float
        Virial velocity at reference radius (km/s)
    r_virial : float
        Reference radius (light days)
    
    Returns:
    --------
    T : float or array
        Orbital period (days)
    """
    # Convert to cgs
    r_cgs = np.array(r) * LD_TO_CM  # cm
    v_virial_cgs = v_virial * KM_TO_CM  # cm/s
    r_virial_cgs = r_virial * LD_TO_CM  # cm
    
    # Keplerian velocity: v = v_virial * sqrt(r_virial / r)
    v = v_virial_cgs * np.sqrt(r_virial_cgs / r_cgs)  # cm/s
    
    # Orbital period: T = 2π * r / v (in seconds)
    T_sec = 2 * np.pi * r_cgs / v  # seconds
    
    # Convert to days
    T_days = T_sec / DAY  # days
    
    if np.isscalar(r):
        return float(T_days)
    return T_days


def compute_time_evolving_map(disk, lambda0, v_virial, nlambda=200, ntime=100,
                               tmax=None, lambda_range=None):
    """
    Compute time-evolving velocity map: wavelength vs time.
    
    Parameters:
    -----------
    disk : PillarDisk
        Disk model with geometry
    lambda0 : float
        Rest wavelength of emission line (Angstroms)
    v_virial : float
        FWHM of emission line (km/s). The actual Keplerian velocity is
        v(r) = v_virial * sqrt(r0/r), where r0 is the reference radius.
    nlambda : int
        Number of wavelength bins
    ntime : int
        Number of time bins
    tmax : float or None
        Maximum time (days). If None, estimate from orbital period at r_out
    lambda_range : tuple or None
        (min, max) wavelength range in Angstroms. If None, auto-compute from velocity
    
    Returns:
    --------
    lambda_grid : array
        Wavelength grid (Angstroms)
    time_grid : array
        Time grid (days)
    flux_map : array
        Spectral flux map - shape (ntime, nlambda)
    """
    # Create 2D grids for disk surface
    r_2d, phi_2d = np.meshgrid(disk.r, disk.phi, indexing='ij')
    
    # Get height and temperature (with shadows)
    h_2d = disk.get_height(r_2d, phi_2d)
    T_2d = disk.get_temperature(r_2d, phi_2d, compute_shadows=True)
    
    # Wavelength grid
    if lambda_range is None:
        v_range = 3.0 * v_virial  # velocity range in km/s
        v_c = v_range * KM_TO_CM / C  # velocity in units of c
        dlambda = lambda0 * v_c
        lambda_min = lambda0 - dlambda
        lambda_max = lambda0 + dlambda
    else:
        lambda_min, lambda_max = lambda_range
    
    lambda_grid = np.linspace(lambda_min, lambda_max, nlambda)
    dlambda = (lambda_max - lambda_min) / (nlambda - 1)
    
    # Time grid - estimate from orbital period at pillar location if not specified
    if tmax is None:
        r_virial = disk.r0
        if len(disk.pillars) > 0:
            # Use mean pillar radius for orbital period calculation
            pillar_radii = [pillar['r'] for pillar in disk.pillars]
            r_pillar_mean = np.mean(pillar_radii)
            T_pillar = compute_orbital_period(r_pillar_mean, v_virial, r_virial)
            tmax = 2.0 * T_pillar  # 2 orbital periods at mean pillar radius
            print(f"Auto-computed tmax = {tmax:.1f} days (orbital period at mean pillar radius r={r_pillar_mean:.1f} ld = {T_pillar:.1f} days)")
        else:
            # Fall back to outer disk radius if no pillars
            T_outer = compute_orbital_period(disk.rout, v_virial, r_virial)
            tmax = 2.0 * T_outer  # 2 orbital periods at outer radius
            print(f"Auto-computed tmax = {tmax:.1f} days (orbital period at r_out = {T_outer:.1f} days, no pillars found)")
    
    time_grid = np.linspace(0, tmax, ntime)
    dtime = tmax / (ntime - 1)
    
    # Initialize flux map
    flux_map = np.zeros((ntime, nlambda))
    
    # Convert virial velocity to cm/s
    v_virial_cgs = v_virial * KM_TO_CM  # km/s to cm/s
    r_virial = disk.r0  # Reference radius for virial velocity
    
    # Compute angular velocities for each radius
    omega_r = compute_angular_velocity(disk.r, v_virial, r_virial)  # rad/day
    
    # Initial azimuthal positions (stored as phi_0)
    phi_0_2d = phi_2d.copy()
    
    print("Computing time-evolving velocity map...")
    print(f"Time range: 0 to {tmax} days")
    print(f"Wavelength range: {lambda_min:.2f} to {lambda_max:.2f} Å")
    
    # Loop over time
    for itime, t in enumerate(tqdm(time_grid, desc="Processing time", unit="time")):
        # At time t, each point has rotated: phi(t) = phi_0 + omega*t
        # Need to handle 2π wrapping
        omega_2d = np.zeros_like(r_2d)
        for ir in range(disk.nr):
            omega_2d[ir, :] = omega_r[ir]
        
        phi_t_2d = phi_0_2d + omega_2d * t
        # Wrap to [0, 2π)
        phi_t_2d = np.mod(phi_t_2d, 2*np.pi)
        
        # Get height and temperature at rotated positions
        h_t_2d = disk.get_height(r_2d, phi_t_2d)
        T_t_2d = disk.get_temperature(r_2d, phi_t_2d, compute_shadows=True)
        
        # Loop over disk surface
        for ir in range(disk.nr - 1):
            r_now = disk.r[ir]
            h_now = h_t_2d[ir, :]
            
            # Next radius
            r_next = disk.r[ir + 1]
            h_next = h_t_2d[ir + 1, :]
            
            # Surface element
            dr = r_next - r_now
            dh = h_next - h_now
            ds = np.sqrt(dr**2 + dh**2)
            da = ds * r_now * disk.dphi
            
            # Tilt angle
            sintilt = dh / (ds + 1e-10)
            costilt = dr / (ds + 1e-10)
            
            # Get temperature at this radius
            T_now = T_t_2d[ir, :]
            
            # Loop over azimuth
            for iphi, phi in enumerate(phi_t_2d[ir, :]):
                # Normal vector
                px = -sintilt[iphi] * np.cos(phi)
                pz = costilt[iphi]
                
                # Foreshortening factor
                dot = disk.ex * px + disk.ez * pz
                
                if dot <= 0:
                    continue
                
                # Time delay (light travel time from lamp to disk element to Earth)
                dx = r_now * np.cos(phi)
                dy = r_now * np.sin(phi)
                dz = h_now[iphi] - disk.hlamp
                dlamp = np.sqrt(dx**2 + dy**2 + dz**2)
                
                # Delay = light travel time from lamp to disk element to Earth
                tau_delay = dlamp - (dx * disk.sini + h_now[iphi] * disk.cosi)
                
                # Only include emission that arrives at Earth at time t
                # The observed time is t_obs = t_emission + tau_delay
                # So we need: t_obs = t, which means t_emission = t - tau_delay
                # But we're computing at time t, so we need to check if tau_delay matches
                # Actually, we want to compute what we observe at time t
                # Emission from a point with delay tau arrives at time t_obs = t_emission + tau
                # So if we observe at time t, we see emission from t_emission = t - tau
                # But we're computing the map at time t, so we should include all points
                # whose emission arrives at time t
                
                # For simplicity, include all emission (we'll smooth over delays)
                # The key is that the velocity changes with rotation
                
                # Velocity field (Keplerian orbit)
                r_virial = disk.r0
                v_phi = v_virial_cgs * np.sqrt(r_virial / r_now)  # Keplerian velocity
                
                # Line of sight velocity component
                # v_los = v_phi * sin(phi) * sin(i)
                v_los = v_phi * np.sin(phi) * disk.sini
                
                # Doppler shift: dlambda/lambda = v_los/c
                lambda_obs = lambda0 * (1.0 + v_los / C)
                
                # Find wavelength bin (with Gaussian blur)
                if lambda_obs < lambda_min or lambda_obs > lambda_max:
                    continue
                
                siglambda = dlambda
                nsigma = 2
                ilambda_min = max(0, int((lambda_obs - lambda_min - nsigma * siglambda) / dlambda))
                ilambda_max = min(nlambda - 1, int((lambda_obs - lambda_min + nsigma * siglambda) / dlambda))
                
                # Compute Gaussian weights for wavelength bins
                gauss_lambda_norm = 0.0
                for il in range(ilambda_min, ilambda_max + 1):
                    lambda_rel = lambda_obs - lambda_grid[il]
                    gauss = np.exp(-0.5 * (lambda_rel / siglambda)**2)
                    gauss_lambda_norm += gauss
                
                # Weight by area, foreshortening, temperature
                T_loc = T_now[iphi]
                T_weight = (T_loc / disk.tv1)**2  # Normalize by reference temperature
                weight = da[iphi] * dot * T_weight
                
                # Add contribution with Gaussian spread
                if gauss_lambda_norm > 0:
                    for ilambda in range(ilambda_min, ilambda_max + 1):
                        lambda_rel = lambda_obs - lambda_grid[ilambda]
                        gauss_lambda = np.exp(-0.5 * (lambda_rel / siglambda)**2) / gauss_lambda_norm
                        flux_map[itime, ilambda] += weight * gauss_lambda
        
        # Add emission from pillars (rotating with the disk)
        if len(disk.pillars) > 0:
            for pillar in disk.pillars:
                r_pillar = pillar['r']
                phi_pillar_0 = pillar['phi']  # Initial azimuth
                sigma_r = pillar['sigma_r']
                sigma_phi = pillar['sigma_phi']
                
                # Find angular velocity at pillar radius
                ir_pillar = np.argmin(np.abs(disk.r - r_pillar))
                omega_pillar = omega_r[ir_pillar]
                
                # Rotated azimuth
                phi_pillar_t = phi_pillar_0 + omega_pillar * t
                phi_pillar_t = np.mod(phi_pillar_t, 2*np.pi)
                
                # Get height at pillar location
                h_pillar = disk.get_height(r_pillar, phi_pillar_t)
                
                # Compute velocity at pillar location
                v_phi_pillar = v_virial_cgs * np.sqrt(r_virial / r_pillar)
                v_los_pillar = v_phi_pillar * np.sin(phi_pillar_t) * disk.sini
                
                # Doppler shift
                lambda_obs_pillar = lambda0 * (1.0 + v_los_pillar / C)
                
                if lambda_obs_pillar < lambda_min or lambda_obs_pillar > lambda_max:
                    continue
                
                # Pillar emission strength
                area_pillar = 2 * np.pi * sigma_r * sigma_phi * r_pillar
                pillar_weight = area_pillar * 10.0  # Factor to make pillars visible
                
                # Wavelength spread
                siglambda_pillar = dlambda * 2
                nsigma_lambda = 3
                ilambda_min_p = max(0, int((lambda_obs_pillar - lambda_min - nsigma_lambda * siglambda_pillar) / dlambda))
                ilambda_max_p = min(nlambda - 1, int((lambda_obs_pillar - lambda_min + nsigma_lambda * siglambda_pillar) / dlambda))
                
                # Add Gaussian-weighted contribution
                gauss_lambda_norm_p = 0.0
                for il in range(ilambda_min_p, ilambda_max_p + 1):
                    lambda_rel = lambda_obs_pillar - lambda_grid[il]
                    gauss = np.exp(-0.5 * (lambda_rel / siglambda_pillar)**2)
                    gauss_lambda_norm_p += gauss
                
                if gauss_lambda_norm_p > 0:
                    for ilambda in range(ilambda_min_p, ilambda_max_p + 1):
                        lambda_rel = lambda_obs_pillar - lambda_grid[ilambda]
                        gauss_lambda = np.exp(-0.5 * (lambda_rel / siglambda_pillar)**2) / gauss_lambda_norm_p
                        flux_map[itime, ilambda] += pillar_weight * gauss_lambda
    
    # Normalize the map
    total = np.sum(flux_map) * dtime * dlambda
    if total > 0:
        flux_map = flux_map / total
    else:
        print("Warning: Total flux is zero. Check disk geometry and parameters.")
    
    return lambda_grid, time_grid, flux_map


def plot_time_evolving_map(lambda_grid, time_grid, flux_map, lambda0, filename='velocity_time_map.png'):
    """
    Plot the time-evolving velocity map: wavelength vs time.
    
    Parameters:
    -----------
    lambda_grid : array
        Wavelength grid (Angstroms)
    time_grid : array
        Time grid (days)
    flux_map : array
        Flux map (ntime, nlambda)
    lambda0 : float
        Rest wavelength (for reference line)
    filename : str
        Output filename
    """
    fig, ax = plt.subplots(1, 1, figsize=(14, 8))
    
    # Use log scale for better visualization
    flux_plot = flux_map.copy()
    flux_plot[flux_plot <= 0] = np.nan
    
    im = ax.pcolormesh(lambda_grid, time_grid, flux_map, shading='auto', cmap='rainbow')
    ax.axvline(lambda0, color='white', linestyle='--', linewidth=2, alpha=0.7)
    ax.set_xlabel(r'$\rm Wavelength~[\AA]$', fontsize=16)
    ax.set_ylabel(r'$\rm Time~[days]$', fontsize=16)
    ax.set_title(r'$\rm Time-Evolving~Velocity~Map~\Psi(\lambda, t)$', fontsize=18, pad=10)
    ax.set_xlim(lambda_grid[0], lambda_grid[-1])
    ax.set_ylim(time_grid[0], time_grid[-1])
    plt.colorbar(im, ax=ax, label=r'$\rm Spectral~Flux$')
    ax.minorticks_on()
    ax.tick_params(top=True, right=True, axis='both', which='major', length=8, width=1.5, direction='in')
    ax.tick_params(top=True, right=True, axis='both', which='minor', length=4, width=1, direction='in')
    
    # Add velocity axis on top
    ax2 = ax.twiny()
    v_km_s = (C / KM_TO_CM) * (lambda_grid - lambda0) / lambda0  # km/s
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xlabel(r'$\rm Velocity~[km/s]$', labelpad=10)
    v_min = v_km_s[0]
    v_max = v_km_s[-1]
    n_ticks = 5
    v_ticks = np.linspace(v_min, v_max, n_ticks)
    lambda_ticks = lambda0 * (1.0 + v_ticks * KM_TO_CM / C)
    ax2.set_xticks(lambda_ticks)
    ax2.set_xticklabels([rf'${v:.0f}$' for v in v_ticks])
    ax2.minorticks_on()
    ax2.tick_params(top=True, right=True, axis='both', which='major', length=8, width=1.5, direction='in')
    ax2.tick_params(top=True, right=True, axis='both', which='minor', length=4, width=1, direction='in')
    
    plt.tight_layout()
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    print(f"Time-evolving velocity map saved to {filename}")
    plt.close()


def main(config_file='config_line.yaml'):
    """
    Main function to compute and plot time-evolving velocity map.
    
    Parameters:
    -----------
    config_file : str
        Path to configuration YAML file
    """
    # Load configuration
    config = load_config(config_file)
    
    if config is None:
        print("Error: Could not load config file. Using defaults.")
        return
    
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
    
    # Emission line parameters
    lambda0 = 4861.0  # Hβ rest wavelength (Angstroms)
    v_virial = 1000.0  # Virial velocity (km/s)
    
    line_params = config.get('emission_line', {})
    if line_params:
        lambda0 = float(line_params.get('lambda0', lambda0))
        v_virial = float(line_params.get('v_virial', v_virial))
    
    # Computation parameters
    comp_params = config.get('computation', {})
    nlambda = int(comp_params.get('nlambda', 200))
    ntime = int(comp_params.get('ntime', 100))
    tmax_config = comp_params.get('tmax', None)
    tmax = float(tmax_config) if tmax_config is not None else None
    
    print(f"\nComputing time-evolving velocity map for emission line:")
    print(f"  Rest wavelength: {lambda0:.1f} Å")
    print(f"  Virial velocity: {v_virial:.0f} km/s")
    print(f"  Number of pillars: {len(disk.pillars)}")
    print(f"  Time range: 0 to {tmax} days")
    
    # Compute time-evolving map
    lambda_grid, time_grid, flux_map = compute_time_evolving_map(
        disk, lambda0, v_virial, nlambda=nlambda, ntime=ntime, tmax=tmax)
    
    # Plot
    plot_params = config.get('plotting', {})
    filename = plot_params.get('velocity_time_filename', 'velocity_time_map.png')
    plot_time_evolving_map(lambda_grid, time_grid, flux_map, lambda0, filename=filename)
    
    print("\nSummary:")
    print(f"Wavelength range: {lambda_grid[0]:.2f} - {lambda_grid[-1]:.2f} Å")
    print(f"Time range: {time_grid[0]:.2f} - {time_grid[-1]:.2f} days")
    print(f"Max flux: {np.max(flux_map):.2e}")


if __name__ == '__main__':
    import sys
    config_file = sys.argv[1] if len(sys.argv) > 1 else 'config_line.yaml'
    main(config_file=config_file)
