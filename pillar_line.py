#!/usr/bin/env python3
"""
pillar_line.py - Compute velocity delay map for emission line from AGN accretion disk with pillars

This script computes the wavelength-delay map Psi(lambda, tau) for an emission line from gas
orbiting the black hole. The gas follows Keplerian orbits with virial velocity.

The disk geometry (including pillars) is loaded from the same config file as pillar_disk.py.
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

# Try to import 3D plotting
try:
    from mpl_toolkits.mplot3d import Axes3D
    HAS_3D = True
except ImportError:
    HAS_3D = False

# Import the PillarDisk class and utilities
# Add current directory to path to ensure we can import pillar_disk
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    from pillar_disk import PillarDisk, load_config
except ImportError as e:
    # Fallback if pillar_disk is not available
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


def compute_velocity_delay_map(disk, lambda0, v_virial, nlambda=200, ntau=200, 
                                taumax=150.0, lambda_range=None):
    """
    Compute velocity delay map Psi(lambda, tau) for an emission line.
    
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
    ntau : int
        Number of delay bins
    taumax : float
        Maximum delay (days)
    lambda_range : tuple or None
        (min, max) wavelength range in Angstroms. If None, auto-compute from velocity
    
    Returns:
    --------
    lambda_grid : array
        Wavelength grid (Angstroms)
    tau_grid : array
        Delay grid (days)
    psi_map : array
        Response function Psi(lambda, tau) - shape (ntau, nlambda)
    """
    # Create 2D grids for disk surface
    r_2d, phi_2d = np.meshgrid(disk.r, disk.phi, indexing='ij')
    
    # Get height and temperature (with shadows)
    h_2d = disk.get_height(r_2d, phi_2d)
    T_2d = disk.get_temperature(r_2d, phi_2d, compute_shadows=True)
    
    # Wavelength grid
    # v_virial is the FWHM of the line, not the edge
    # For a Gaussian line profile, FWHM ≈ 2.35 * sigma
    # For velocity range, we want to cover the full line profile
    # Use ±3*FWHM for wavelength range to show larger range
    if lambda_range is None:
        # Auto-compute from velocity: dlambda/lambda = v/c
        # v_virial is FWHM, so use ±3*FWHM for wavelength range
        v_range = 3.0 * v_virial  # velocity range in km/s (extended beyond FWHM)
        v_c = v_range * KM_TO_CM / C  # velocity in units of c
        dlambda = lambda0 * v_c  # range for ±v_range
        lambda_min = lambda0 - dlambda
        lambda_max = lambda0 + dlambda
    else:
        lambda_min, lambda_max = lambda_range
    
    lambda_grid = np.linspace(lambda_min, lambda_max, nlambda)
    dlambda = (lambda_max - lambda_min) / (nlambda - 1)
    
    # Delay grid
    # If taumax not specified, estimate from disk outer radius
    # Maximum delay occurs at outer edge: tau_max ≈ sqrt(rout^2 + h^2) - (rout*sini + h*cosi)
    if taumax is None:
        # More accurate: tau_max = sqrt(rout^2 + (h_outer - hlamp)^2) - (rout*sini + h_outer*cosi)
        # For safety, use a factor larger than rout
        h_outer = disk.get_height(disk.rout, 0.0)  # Height at outer radius
        dlamp_max = np.sqrt(disk.rout**2 + (h_outer - disk.hlamp)**2)
        tau_max_est = dlamp_max - (disk.rout * disk.sini + h_outer * disk.cosi)
        taumax = max(1.5 * disk.rout, tau_max_est * 1.2)  # Add 20% safety margin
        print(f"Auto-computed taumax = {taumax:.1f} days (rout = {disk.rout:.1f} ld)")
    
    tau_grid = np.linspace(0, taumax, ntau)
    dtau = taumax / (ntau - 1)
    
    # Maximum reasonable delay based on disk size
    # Cap delays to prevent unrealistic values beyond disk edge
    tau_max_reasonable = taumax * 1.1  # Allow 10% beyond taumax for Gaussian blur
    
    # Initialize response map
    psi_map = np.zeros((ntau, nlambda))
    
    # Convert virial velocity to cm/s
    # v_virial is the FWHM of the line profile (km/s)
    v_virial_cgs = v_virial * KM_TO_CM  # km/s to cm/s
    
    # Loop over disk surface
    print("Computing velocity delay map...")
    for ir in tqdm(range(disk.nr - 1), desc="Processing radii", unit="radius"):
        r_now = disk.r[ir]
        h_now = h_2d[ir, :]
        
        # Next radius
        r_next = disk.r[ir + 1]
        h_next = h_2d[ir + 1, :]
        
        # Surface element
        dr = r_next - r_now
        dh = h_next - h_now
        ds = np.sqrt(dr**2 + dh**2)
        da = ds * r_now * disk.dphi
        
        # Tilt angle
        sintilt = dh / (ds + 1e-10)
        costilt = dr / (ds + 1e-10)
        
        # Get temperature at this radius
        T_now = T_2d[ir, :]
        
        # Loop over azimuth
        for iphi, phi in enumerate(disk.phi):
            # Normal vector
            px = -sintilt[iphi] * np.cos(phi)
            pz = costilt[iphi]
            
            # Foreshortening factor
            dot = disk.ex * px + disk.ez * pz
            
            if dot <= 0:
                continue
            
            # Time delay
            # Distance from lamp
            dx = r_now * np.cos(phi)
            dy = r_now * np.sin(phi)
            dz = h_now[iphi] - disk.hlamp
            dlamp = np.sqrt(dx**2 + dy**2 + dz**2)
            
            # Delay = light travel time from lamp to disk element to Earth
            # Earth direction: (sini, 0, cosi)
            tau_delay = dlamp - (dx * disk.sini + h_now[iphi] * disk.cosi)
            
            if tau_delay < 0 or tau_delay > tau_max_reasonable:
                continue
            
            # Get temperature at this location (includes shadow effects)
            T_loc = T_now[iphi]
            
            # Velocity field (Keplerian orbit)
            # For circular orbit in disk plane: v = sqrt(GM/r)
            # The velocity is tangential (in azimuthal direction)
            # For a circular orbit at radius r:
            # - Tangential velocity: v_phi = sqrt(GM/r) (in disk plane)
            # - Line of sight velocity: v_los = v_phi * sin(phi) * sin(i)
            #   where phi is the azimuthal angle (0 = receding, pi = approaching)
            #   and i is the inclination angle
            
            # Compute Keplerian velocity at this radius
            # v(r) = sqrt(GM/r)
            # We use v_virial as the velocity at a reference radius r_virial
            # v_virial = sqrt(GM/r_virial), so GM = v_virial^2 * r_virial
            # Use r0 (reference radius) as r_virial for scaling
            r_virial = disk.r0  # Reference radius for virial velocity
            v_phi = v_virial_cgs * np.sqrt(r_virial / r_now)  # Keplerian velocity at radius r
            
            # Line of sight velocity component
            # v_los = v_phi * sin(phi) * sin(i)
            # For phi=0: v_los = 0 (sideways motion)
            # For phi=pi/2: v_los = v_phi * sin(i) (maximum receding)
            # For phi=-pi/2: v_los = -v_phi * sin(i) (maximum approaching)
            v_los = v_phi * np.sin(phi) * disk.sini
            
            # Doppler shift: dlambda/lambda = v_los/c
            # lambda_obs = lambda0 * (1 + v_los/c)
            lambda_obs = lambda0 * (1.0 + v_los / C)
            
            # Find wavelength bin (with Gaussian blur for smoothness)
            if lambda_obs < lambda_min or lambda_obs > lambda_max:
                continue
            
            # Wavelength bin with Gaussian spread
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
            
            # Find delay bin (with Gaussian blur)
            sigtau = dtau
            nsigma_tau = 3
            ibin_min = max(0, int((tau_delay - nsigma_tau * sigtau) / dtau))
            ibin_max = min(ntau - 1, int((tau_delay + nsigma_tau * sigtau) / dtau))
            
            # Add contribution with 2D Gaussian spread (both wavelength and delay)
            if gauss_lambda_norm > 0:
                for ibin in range(ibin_min, ibin_max + 1):
                    tau_rel = tau_delay - tau_grid[ibin]
                    gauss_tau = np.exp(-0.5 * (tau_rel / sigtau)**2)
                    
                    for ilambda in range(ilambda_min, ilambda_max + 1):
                        lambda_rel = lambda_obs - lambda_grid[ilambda]
                        gauss_lambda = np.exp(-0.5 * (lambda_rel / siglambda)**2) / gauss_lambda_norm
                        
                        # Weight by area, foreshortening, temperature, and both Gaussians
                        # Temperature affects emission line strength (shadows reduce T, so reduce emission)
                        # Use T^2 scaling for line emission (proportional to density * temperature)
                        T_loc = T_now[iphi]
                        T_weight = (T_loc / disk.tv1)**2  # Normalize by reference temperature
                        weight = da[iphi] * dot * T_weight * gauss_tau * gauss_lambda
                        psi_map[ibin, ilambda] += weight
    
    # Add emission from pillars (if any)
    # Pillars emit the line directly, adding to the response map
    if len(disk.pillars) > 0:
        print("Adding pillar emission to velocity delay map...")
        for pillar in disk.pillars:
            # Pillar parameters - pillars are stored with keys 'r' and 'phi'
            r_pillar = pillar['r']
            phi_pillar = pillar['phi']
            sigma_r = pillar['sigma_r']
            sigma_phi = pillar['sigma_phi']
            
            # Find the radial bin closest to pillar position
            ir_pillar = np.argmin(np.abs(disk.r - r_pillar))
            r_pillar_actual = disk.r[ir_pillar]
            
            # Find the azimuthal bin closest to pillar position
            # Handle phi wrapping
            phi_diff = np.abs(disk.phi - phi_pillar)
            phi_diff = np.minimum(phi_diff, 2*np.pi - phi_diff)
            iphi_pillar = np.argmin(phi_diff)
            phi_pillar_actual = disk.phi[iphi_pillar]
            
            # Get height at pillar location
            h_pillar = disk.get_height(r_pillar_actual, phi_pillar_actual)
            
            # Compute velocity at pillar location (Keplerian)
            # v(r) = sqrt(GM/r) = v_virial * sqrt(r_virial/r)
            r_virial = disk.r0  # Reference radius for virial velocity
            v_phi_pillar = v_virial_cgs * np.sqrt(r_virial / r_pillar_actual)  # Keplerian velocity
            # Line of sight velocity: v_los = v_phi * sin(phi) * sin(i)
            v_los_pillar = v_phi_pillar * np.sin(phi_pillar_actual) * disk.sini
            
            # Doppler shift
            lambda_obs_pillar = lambda0 * (1.0 + v_los_pillar / C)
            
            # Skip if outside wavelength range
            if lambda_obs_pillar < lambda_min or lambda_obs_pillar > lambda_max:
                continue
            
            # Time delay at pillar location
            dx_pillar = r_pillar_actual * np.cos(phi_pillar_actual)
            dy_pillar = r_pillar_actual * np.sin(phi_pillar_actual)
            dz_pillar = h_pillar - disk.hlamp
            dlamp_pillar = np.sqrt(dx_pillar**2 + dy_pillar**2 + dz_pillar**2)
            tau_delay_pillar = dlamp_pillar - (dx_pillar * disk.sini + h_pillar * disk.cosi)
            
            if tau_delay_pillar < 0 or tau_delay_pillar > tau_max_reasonable:
                continue
            
            # Pillar emission strength (weighted by Gaussian profile and area)
            # The pillar contributes emission proportional to its area
            # Use the Gaussian widths to estimate the effective area
            area_pillar = 2 * np.pi * sigma_r * sigma_phi * r_pillar_actual
            
            # Add emission with Gaussian spread in both wavelength and delay
            # Wavelength spread
            siglambda_pillar = dlambda * 2  # Wider spread for pillar emission
            nsigma_lambda = 3
            ilambda_min_p = max(0, int((lambda_obs_pillar - lambda_min - nsigma_lambda * siglambda_pillar) / dlambda))
            ilambda_max_p = min(nlambda - 1, int((lambda_obs_pillar - lambda_min + nsigma_lambda * siglambda_pillar) / dlambda))
            
            # Delay spread
            sigtau_pillar = dtau * 2  # Wider spread for pillar emission
            nsigma_tau = 3
            ibin_min_p = max(0, int((tau_delay_pillar - nsigma_tau * sigtau_pillar) / dtau))
            ibin_max_p = min(ntau - 1, int((tau_delay_pillar + nsigma_tau * sigtau_pillar) / dtau))
            
            # Add Gaussian-weighted contribution
            gauss_lambda_norm_p = 0.0
            for il in range(ilambda_min_p, ilambda_max_p + 1):
                lambda_rel = lambda_obs_pillar - lambda_grid[il]
                gauss = np.exp(-0.5 * (lambda_rel / siglambda_pillar)**2)
                gauss_lambda_norm_p += gauss
            
            if gauss_lambda_norm_p > 0:
                for ibin in range(ibin_min_p, ibin_max_p + 1):
                    tau_rel = tau_delay_pillar - tau_grid[ibin]
                    gauss_tau = np.exp(-0.5 * (tau_rel / sigtau_pillar)**2)
                    
                    for ilambda in range(ilambda_min_p, ilambda_max_p + 1):
                        lambda_rel = lambda_obs_pillar - lambda_grid[ilambda]
                        gauss_lambda = np.exp(-0.5 * (lambda_rel / siglambda_pillar)**2) / gauss_lambda_norm_p
                        
                        # Weight by pillar area (pillars are bright emitters)
                        # Use a factor to make pillar emission significant
                        pillar_weight = area_pillar * 10.0  # Factor to make pillars visible
                        weight = pillar_weight * gauss_tau * gauss_lambda
                        psi_map[ibin, ilambda] += weight
    
    # Normalize the map
    # Normalize so that sum over all wavelengths and delays gives total response
    total = np.sum(psi_map) * dtau * dlambda
    if total > 0:
        psi_map = psi_map / total
    else:
        print("Warning: Total response is zero. Check disk geometry and parameters.")
    
    return lambda_grid, tau_grid, psi_map


def plot_temperature_profile(disk, filename='temperature_profile.png'):
    """
    Plot azimuthally-averaged temperature profile T(r) with and without pillars,
    showing T_visc, T_tot, and T_irrad on the same plot.
    
    Parameters:
    -----------
    disk : PillarDisk
        Disk model
    filename : str
        Output filename
    """
    print("Computing temperature profiles...")
    r_profile, T_with, T_without, T_visc, T_irrad_with, T_irrad_without = disk.get_temperature_profile(compute_shadows=True)
    
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    
    # Plot all components on the same plot
    ax.plot(r_profile, T_visc, 'g-', linewidth=2, label=r'$\rm T_{visc}$')
    ax.plot(r_profile, T_without, 'r--', linewidth=2, label=r'$\rm T_{total}~without~pillars$')
    ax.plot(r_profile, T_with, 'b-', linewidth=2, label=r'$\rm T_{total}~with~pillars$')
    ax.plot(r_profile, T_irrad_without, 'r:', linewidth=2, alpha=0.7, label=r'$\rm T_{irrad}~without~pillars$')
    ax.plot(r_profile, T_irrad_with, 'b:', linewidth=2, alpha=0.7, label=r'$\rm T_{irrad}~with~pillars$')
    
    ax.set_xlabel(r'$\rm Radius~[light~days]$')
    ax.set_ylabel(r'$\rm Temperature~[K]$')
    ax.set_title(r'$\rm Azimuthally-Averaged~Temperature~Profile~T(r)$')
    ax.legend()
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3)
    ax.minorticks_on()
    ax.tick_params(top=True, right=True, axis='both', which='major', length=8, width=1.5, direction='in')
    ax.tick_params(top=True, right=True, axis='both', which='minor', length=4, width=1, direction='in')
    
    plt.tight_layout()
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    print(f"Saved temperature profile plot to {filename}")
    plt.close()


def plot_velocity_delay_map(lambda_grid, tau_grid, psi_map, lambda0, 
                            psi_map_no_pillars=None, filename='velocity_delay_map.png'):
    """
    Plot the velocity delay map Psi(lambda, tau).
    
    Parameters:
    -----------
    lambda_grid : array
        Wavelength grid (Angstroms)
    tau_grid : array
        Delay grid (days)
    psi_map : array
        Response function with pillars (ntau, nlambda)
    lambda0 : float
        Rest wavelength (for reference line)
    psi_map_no_pillars : array, optional
        Response function without pillars (ntau, nlambda)
    filename : str
        Output filename
    """
    if psi_map_no_pillars is not None:
        # Comparison plot: with and without pillars
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        # Determine common color scale for both maps
        vmin = min(np.min(psi_map[psi_map > 0]), np.min(psi_map_no_pillars[psi_map_no_pillars > 0]))
        vmax = 2*np.max(psi_map_no_pillars)
        
        # Top left: With pillars
        ax = axes[0, 0]
        psi_plot = psi_map.copy()
        psi_plot[psi_plot <= 0] = np.nan
        im1 = ax.pcolormesh(lambda_grid, tau_grid, psi_map, shading='auto', cmap='rainbow', vmin=vmin, vmax=vmax)
        ax.axvline(lambda0, color='red', linestyle='--', linewidth=2)
        ax.set_xlabel(r'$\rm Wavelength~[\AA]$')
        ax.set_ylabel(r'$\rm Delay~[days]$')
        ax.set_title(r'$\rm With~Pillars~\Psi(\lambda, \tau)$', pad=10)
        ax.set_xlim(lambda_grid[0], lambda_grid[-1])
        ax.set_ylim(0, 125)
        plt.colorbar(im1, ax=ax, label=r'$\rm Response~Function$')
        ax.minorticks_on()
        ax.tick_params(top=True, right=True, axis='both', which='major', length=8, width=1.5, direction='in')
        ax.tick_params(top=True, right=True, axis='both', which='minor', length=4, width=1, direction='in')
        
        # Add velocity axis on top
        # v = c * (λ - λ0) / λ0, convert to km/s
        ax2 = ax.twiny()
        v_km_s = (C / KM_TO_CM) * (lambda_grid - lambda0) / lambda0  # km/s
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xlabel(r'$\rm Velocity~[km/s]$', labelpad=10)
        # Set velocity ticks at reasonable intervals
        v_min = v_km_s[0]
        v_max = v_km_s[-1]
        n_ticks = 5
        v_ticks = np.linspace(v_min, v_max, n_ticks)
        # Convert velocity ticks back to wavelength for positioning
        lambda_ticks = lambda0 * (1.0 + v_ticks * KM_TO_CM / C)
        ax2.set_xticks(lambda_ticks)
        ax2.set_xticklabels([rf'${v:.0f}$' for v in v_ticks])
        ax2.minorticks_on()
        ax2.tick_params(top=True, right=True, axis='both', which='major', length=8, width=1.5, direction='in')
        ax2.tick_params(top=True, right=True, axis='both', which='minor', length=4, width=1, direction='in')
        
        # Top right: Without pillars
        ax = axes[0, 1]
        psi_plot_no = psi_map_no_pillars.copy()
        psi_plot_no[psi_plot_no <= 0] = np.nan
        im2 = ax.pcolormesh(lambda_grid, tau_grid, psi_map_no_pillars, shading='auto', cmap='rainbow', vmin=vmin, vmax=vmax)
        ax.axvline(lambda0, color='red', linestyle='--', linewidth=2)
        ax.set_xlabel(r'$\rm Wavelength~[\AA]$')
        ax.set_ylabel(r'$\rm Delay~[days]$')
        ax.set_title(r'$\rm Without~Pillars~\Psi(\lambda, \tau)$', pad=10)
        ax.set_xlim(lambda_grid[0], lambda_grid[-1])
        ax.set_ylim(0, 125)
        plt.colorbar(im2, ax=ax, label=r'$\rm Response~Function$')
        ax.minorticks_on()
        ax.tick_params(top=True, right=True, axis='both', which='major', length=8, width=1.5, direction='in')
        ax.tick_params(top=True, right=True, axis='both', which='minor', length=4, width=1, direction='in')
        
        # Add velocity axis on top
        ax2 = ax.twiny()
        v_km_s = (C / KM_TO_CM) * (lambda_grid - lambda0) / lambda0  # km/s
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xlabel(r'$\rm Velocity~[km/s]$', labelpad=10)
        # Set velocity ticks at reasonable intervals
        v_min = v_km_s[0]
        v_max = v_km_s[-1]
        n_ticks = 5
        v_ticks = np.linspace(v_min, v_max, n_ticks)
        # Convert velocity ticks back to wavelength for positioning
        lambda_ticks = lambda0 * (1.0 + v_ticks * KM_TO_CM / C)
        ax2.set_xticks(lambda_ticks)
        ax2.set_xticklabels([rf'${v:.0f}$' for v in v_ticks])
        ax2.minorticks_on()
        ax2.tick_params(top=True, right=True, axis='both', which='major', length=8, width=1.5, direction='in')
        ax2.tick_params(top=True, right=True, axis='both', which='minor', length=4, width=1, direction='in')
        
        # Bottom left: Mean delay with pillars
        ax = axes[1, 0]
        tau_mean = np.zeros(len(lambda_grid))
        for i in range(len(lambda_grid)):
            mask = psi_map[:, i] > 0
            if np.any(mask):
                tau_mean[i] = np.sum(psi_map[:, i] * tau_grid) / np.sum(psi_map[:, i])
        ax.plot(lambda_grid, tau_mean, 'b-', linewidth=2, label=r'$\rm with~pillars$')
        tau_mean_no = np.zeros(len(lambda_grid))
        for i in range(len(lambda_grid)):
            mask = psi_map_no_pillars[:, i] > 0
            if np.any(mask):
                tau_mean_no[i] = np.sum(psi_map_no_pillars[:, i] * tau_grid) / np.sum(psi_map_no_pillars[:, i])
        ax.plot(lambda_grid, tau_mean_no, 'r--', linewidth=2, label=r'$\rm without~pillars$')
        ax.axvline(lambda0, color='red', linestyle=':', linewidth=1, alpha=0.5)
        ax.set_xlabel(r'$\rm Wavelength~[\AA]$')
        ax.set_ylabel(r'$\rm Mean~Delay~[days]$')
        ax.set_xlim(lambda_grid[0], lambda_grid[-1])
        ax.set_ylim(0, 125)
        ax.legend()
        ax.minorticks_on()
        ax.tick_params(top=True, right=True, axis='both', which='major', length=8, width=1.5, direction='in')
        ax.tick_params(top=True, right=True, axis='both', which='minor', length=4, width=1, direction='in')
        
        # Add velocity axis on top
        ax2 = ax.twiny()
        v_km_s = (C / KM_TO_CM) * (lambda_grid - lambda0) / lambda0  # km/s
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xlabel(r'$\rm Velocity~[km/s]$', labelpad=10)
        # Set velocity ticks at reasonable intervals
        v_min = v_km_s[0]
        v_max = v_km_s[-1]
        n_ticks = 5
        v_ticks = np.linspace(v_min, v_max, n_ticks)
        # Convert velocity ticks back to wavelength for positioning
        lambda_ticks = lambda0 * (1.0 + v_ticks * KM_TO_CM / C)
        ax2.set_xticks(lambda_ticks)
        ax2.set_xticklabels([rf'${v:.0f}$' for v in v_ticks])
        ax2.minorticks_on()
        ax2.tick_params(top=True, right=True, axis='both', which='major', length=8, width=1.5, direction='in')
        ax2.tick_params(top=True, right=True, axis='both', which='minor', length=4, width=1, direction='in')
        
        # Bottom right: Difference map
        ax = axes[1, 1]
        diff_map = psi_map - psi_map_no_pillars
        diff_plot = diff_map.copy()
        diff_plot[diff_plot == 0] = np.nan
        im3 = ax.pcolormesh(lambda_grid, tau_grid, diff_map, shading='auto', cmap='RdBu_r')
        ax.axvline(lambda0, color='black', linestyle='--', linewidth=2)
        ax.set_xlabel(r'$\rm Wavelength~[\AA]$')
        ax.set_ylabel(r'$\rm Delay~[days]$')
        ax.set_title(r'$\rm Difference~\Psi_{\rm with} - \Psi_{\rm without}$')
        ax.set_xlim(lambda_grid[0], lambda_grid[-1])
        ax.set_ylim(0, 125)
        plt.colorbar(im3, ax=ax, label=r'$\rm \Delta Response~Function$')
        ax.minorticks_on()
        ax.tick_params(top=True, right=True, axis='both', which='major', length=8, width=1.5, direction='in')
        ax.tick_params(top=True, right=True, axis='both', which='minor', length=4, width=1, direction='in')
        
        # Add velocity axis on top
        ax2 = ax.twiny()
        v_km_s = (C / KM_TO_CM) * (lambda_grid - lambda0) / lambda0  # km/s
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xlabel(r'$\rm Velocity~[km/s]$', labelpad=10)
        # Set velocity ticks at reasonable intervals
        v_min = v_km_s[0]
        v_max = v_km_s[-1]
        n_ticks = 5
        v_ticks = np.linspace(v_min, v_max, n_ticks)
        # Convert velocity ticks back to wavelength for positioning
        lambda_ticks = lambda0 * (1.0 + v_ticks * KM_TO_CM / C)
        ax2.set_xticks(lambda_ticks)
        ax2.set_xticklabels([rf'${v:.0f}$' for v in v_ticks])
        ax2.minorticks_on()
        ax2.tick_params(top=True, right=True, axis='both', which='major', length=8, width=1.5, direction='in')
        ax2.tick_params(top=True, right=True, axis='both', which='minor', length=4, width=1, direction='in')
    else:
        # Single plot: only with pillars
        fig, axes = plt.subplots(2, 1, figsize=(12, 10))
        
        # 2D map
        ax = axes[0]
        # Use log scale for better visualization
        psi_plot = psi_map.copy()
        psi_plot[psi_plot <= 0] = np.nan
        im = ax.pcolormesh(lambda_grid, tau_grid, psi_map, shading='auto', cmap='rainbow')
        ax.axvline(lambda0, color='red', linestyle='--', linewidth=2)
        ax.set_xlabel(r'$\rm Wavelength~[\AA]$')
        ax.set_ylabel(r'$\rm Delay~[days]$')
        ax.set_title(r'$\rm Velocity~Delay~Map~\Psi(\lambda, \tau)$')
        ax.set_xlim(lambda_grid[0], lambda_grid[-1])
        ax.set_ylim(0, 125)
        plt.colorbar(im, ax=ax, label=r'$\rm Response~Function$')
        ax.minorticks_on()
        ax.tick_params(top=True, right=True, axis='both', which='major', length=8, width=1.5, direction='in')
        ax.tick_params(top=True, right=True, axis='both', which='minor', length=4, width=1, direction='in')
        
        # Add velocity axis on top
        ax2 = ax.twiny()
        v_km_s = (C / KM_TO_CM) * (lambda_grid - lambda0) / lambda0  # km/s
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xlabel(r'$\rm Velocity~[km/s]$', labelpad=10)
        # Set velocity ticks at reasonable intervals
        v_min = v_km_s[0]
        v_max = v_km_s[-1]
        n_ticks = 5
        v_ticks = np.linspace(v_min, v_max, n_ticks)
        # Convert velocity ticks back to wavelength for positioning
        lambda_ticks = lambda0 * (1.0 + v_ticks * KM_TO_CM / C)
        ax2.set_xticks(lambda_ticks)
        ax2.set_xticklabels([rf'${v:.0f}$' for v in v_ticks])
        ax2.minorticks_on()
        ax2.tick_params(top=True, right=True, axis='both', which='major', length=8, width=1.5, direction='in')
        ax2.tick_params(top=True, right=True, axis='both', which='minor', length=4, width=1, direction='in')
        
        # Integrated delay spectrum (sum over wavelengths)
        ax = axes[1]
        tau_mean = np.zeros(len(lambda_grid))
        for i in range(len(lambda_grid)):
            mask = psi_map[:, i] > 0
            if np.any(mask):
                tau_mean[i] = np.sum(psi_map[:, i] * tau_grid) / np.sum(psi_map[:, i])
        
        ax.plot(lambda_grid, tau_mean, 'b-', linewidth=2)
        ax.axvline(lambda0, color='red', linestyle='--', linewidth=2)
        ax.set_xlabel(r'$\rm Wavelength~[\AA]$')
        ax.set_ylabel(r'$\rm Mean~Delay~[days]$')
        ax.set_xlim(lambda_grid[0], lambda_grid[-1])
        ax.set_ylim(0, 125)
        ax.minorticks_on()
        ax.tick_params(top=True, right=True, axis='both', which='major', length=8, width=1.5, direction='in')
        ax.tick_params(top=True, right=True, axis='both', which='minor', length=4, width=1, direction='in')
    
    plt.tight_layout()
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    print(f"Velocity delay map saved to {filename}")
    plt.close()


def main(config_file='config_line.yaml'):
    """
    Main function to compute and plot velocity delay map.
    
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
    
    # Extract disk parameters
    disk_params = config.get('disk', {}).copy()
    temp_params = config.get('temperature', {}).copy()
    lamp_params = config.get('lamp', {}).copy()
    obs_params = config.get('observation', {}).copy()
    
    # Convert parameters (same as in pillar_disk.py)
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
    
    # Process pillars
    pillars_config = config.get('pillars', {})
    
    # Check if we should generate random pillars
    make_many = pillars_config.get('make_many', False)
    if isinstance(make_many, str):
        make_many = make_many.lower() in ('true', '1', 'yes')
    make_many = bool(make_many)
    
    if make_many:
        # Generate random pillars
        print("Generating random pillars...")
        N_pillar = int(pillars_config.get('N_pillar', 100))
        r_mean = float(pillars_config.get('r_mean', 25.0))
        sig_r = float(pillars_config.get('sig_r', 10.0))
        rmin = pillars_config.get('rmin', None)
        rmin = float(rmin) if rmin is not None else None
        
        # Get pillar properties (can be single values or lists)
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
        
        # Expand lists to match N_pillar
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
        
        # Generate random pillar positions
        # r: Gaussian distribution, clipped to disk range
        r_pillar_random = np.random.normal(r_mean, sig_r, N_pillar)
        rmin_eff = max(disk.rin, rmin) if rmin is not None else disk.rin
        r_pillar_random = np.clip(r_pillar_random, rmin_eff, disk.rout)
        
        # phi: Uniform distribution from 0 to 2*pi
        phi_pillar_random = np.random.uniform(0, 2*np.pi, N_pillar)
        
        # Convert boolean strings if needed
        def convert_bool(val):
            if isinstance(val, str):
                return val.lower() in ('true', '1', 'yes')
            return bool(val)
        
        # Create pillar configurations
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
        print(f"  r: mean={r_mean:.1f}, sigma={sig_r:.1f} (clipped to [{rmin_eff:.1f}, {disk.rout:.1f}])")
        print(f"  phi: uniform [0, 2π]")
    else:
        # Manual pillar placement (same as before)
        if isinstance(pillars_config, dict):
            # New format: dict with lists
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
        else:
            # Old format: list of dicts
            for pillar in pillars_config:
                for key in ['r_pillar', 'height', 'sigma_r', 'temp_factor']:
                    if key in pillar:
                        pillar[key] = float(pillar[key])
                for key in ['phi_pillar', 'sigma_phi']:
                    if key in pillar:
                        pillar[key] = parse_math_expr(pillar[key])
                for key in ['modify_height', 'modify_temp']:
                    if key in pillar:
                        if isinstance(pillar[key], str):
                            pillar[key] = pillar[key].lower() in ('true', '1', 'yes')
                        else:
                            pillar[key] = bool(pillar[key])
    
    # Add pillars
    for pillar_cfg in pillars_config:
        disk.add_pillar(**pillar_cfg)
    
    # Emission line parameters (can be added to config later)
    lambda0 = 4861.0  # Hβ rest wavelength (Angstroms) - can be changed
    v_virial = 1000.0  # Virial velocity (km/s)
    
    # Check if line parameters are in config
    line_params = config.get('emission_line', {})
    if line_params:
        lambda0 = float(line_params.get('lambda0', lambda0))
        v_virial = float(line_params.get('v_virial', v_virial))
    
    # Computation parameters
    comp_params = config.get('computation', {})
    nlambda = int(comp_params.get('nlambda', 200))
    ntau = int(comp_params.get('ntau_line', comp_params.get('ntau', 200)))
    # If taumax is specified, use it; otherwise auto-compute from disk size
    taumax_config = comp_params.get('taumax', None)
    taumax = float(taumax_config) if taumax_config is not None else None
    
    print(f"\nComputing velocity delay map for emission line:")
    print(f"  Rest wavelength: {lambda0:.1f} Å")
    print(f"  Virial velocity: {v_virial:.0f} km/s")
    print(f"  Number of pillars: {len(disk.pillars)}")
    
    # Compute velocity delay map with pillars
    print("Computing velocity delay map with pillars...")
    lambda_grid, tau_grid, psi_map = compute_velocity_delay_map(
        disk, lambda0, v_virial, nlambda=nlambda, ntau=ntau, taumax=taumax)
    
    # Compute velocity delay map without pillars
    print("Computing velocity delay map without pillars...")
    # Temporarily disable pillars
    pillars_backup = disk.pillars.copy()
    disk.pillars = []
    lambda_grid_no, tau_grid_no, psi_map_no_pillars = compute_velocity_delay_map(
        disk, lambda0, v_virial, nlambda=nlambda, ntau=ntau, taumax=taumax)
    # Restore pillars
    disk.pillars = pillars_backup
    
    # Plot disk geometry
    plot_params = config.get('plotting', {})
    plot_3d = plot_params.get('plot_3d_geometry', True)
    if plot_3d:
        print("\nPlotting 3D disk geometry...")
        n_phi_plot = plot_params.get('n_phi_plot', 360)
        n_r_plot = plot_params.get('n_r_plot', 500)
        show_light_rays = plot_params.get('show_light_rays', True)
        plot3d_filename = plot_params.get('plot3d_filename', 'pillar_line_3d.png')
        disk.plot_3d_geometry(n_phi_plot=n_phi_plot, n_r_plot=n_r_plot, 
                              show_light_rays=show_light_rays, filename=plot3d_filename)
    
    # Plot velocity delay map
    filename = plot_params.get('velocity_delay_filename', 'velocity_delay_map.png')
    plot_velocity_delay_map(lambda_grid, tau_grid, psi_map, lambda0, 
                           psi_map_no_pillars=psi_map_no_pillars, filename=filename)
    
    # Plot temperature profile
    temp_filename = plot_params.get('temperature_profile_filename', 'temperature_profile.png')
    plot_temperature_profile(disk, filename=temp_filename)
    
    print("\nSummary:")
    print(f"Wavelength range: {lambda_grid[0]:.2f} - {lambda_grid[-1]:.2f} Å")
    print(f"Delay range: {tau_grid[0]:.2f} - {tau_grid[-1]:.2f} days")
    print(f"Max response: {np.max(psi_map):.2e}")
    print(f"Disk outer radius: {disk.rout:.1f} light days")


if __name__ == '__main__':
    import sys
    # Allow config file to be specified as command line argument
    # Default to config_line.yaml for emission line computations
    config_file = sys.argv[1] if len(sys.argv) > 1 else 'config_line.yaml'
    main(config_file=config_file)

