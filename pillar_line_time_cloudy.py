#!/usr/bin/env python3
"""
pillar_line_time_cloudy.py - Time-evolving velocity delay maps with Cloudy models

This script combines the time evolution from pillar_line_time.py with Cloudy model
integration from pillar_line_cloudy.py to produce barber-pole patterns for multiple
emission lines (Halpha, MgII, CIV).

The disk and pillars rotate in circular orbits, and at each time step:
1. Compute ionizing flux map (based on temperature/shadowing)
2. Interpolate Cloudy models to get line intensities
3. Compute velocity-delay map for each line
4. Plot barber-pole patterns showing how line profiles evolve with time
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from multiprocessing import Pool, cpu_count

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

# Import Cloudy functions from pillar_line_cloudy.py
try:
    from pillar_line_cloudy import load_cloudy_models, compute_ionizing_flux
except ImportError as e:
    print(f"Error: pillar_line_cloudy.py not found. Make sure it's in the same directory.")
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


def _process_time_step(args):
    """
    Worker function to process a single time step (for multiprocessing).
    
    Parameters:
    -----------
    args : tuple
        (itime, t, disk_params, cloudy_interp_dict, phi_grid, lambda0_dict, 
         v_virial_cgs, r_virial, omega_r, r_2d, phi_0_2d, lambda_grid_dict, 
         nlambda, Z_target)
    
    Returns:
    --------
    itime : int
        Time step index
    flux_slice_dict : dict
        Dictionary of flux slices for each line (1D arrays of length nlambda)
    """
    (itime, t, disk_params, pillars_list, cloudy_interp_dict, phi_grid, lambda0_dict, 
     v_virial_cgs, r_virial, omega_r, r_2d, phi_0_2d, lambda_grid_dict, 
     nlambda, Z_target) = args
    
    # Recreate disk object (needed for multiprocessing)
    import sys
    import os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from pillar_disk import PillarDisk
    from pillar_line_cloudy import compute_ionizing_flux
    
    # Create disk object
    disk = PillarDisk(**disk_params)
    
    # Add pillars to the disk
    # Note: pillar dictionaries use 'r' and 'phi', but add_pillar expects 'r_pillar' and 'phi_pillar'
    for pillar in pillars_list:
        pillar_kwargs = {
            'r_pillar': pillar['r'],
            'phi_pillar': pillar['phi'],
            'height': pillar.get('height', 0.01),
            'sigma_r': pillar.get('sigma_r', 1.0),
            'sigma_phi': pillar.get('sigma_phi', 0.1),
            'modify_height': pillar.get('modify_height', True),
            'modify_temp': pillar.get('modify_temp', False),
            'temp_factor': pillar.get('temp_factor', 1.5)
        }
        disk.add_pillar(**pillar_kwargs)
    
    # At time t, each point has rotated: phi(t) = phi_0 + omega*t
    omega_2d = np.zeros_like(r_2d)
    for ir in range(disk.nr):
        omega_2d[ir, :] = omega_r[ir]

    phi_t_2d = phi_0_2d + omega_2d * t
    # Wrap to [0, 2π)
    phi_t_2d = np.mod(phi_t_2d, 2*np.pi)

    # Rotate pillar positions to their current azimuth so that
    # get_height / get_temperature / _compute_shadow_mask see the
    # pillar at the correct lab-frame position at this time step.
    original_pillar_phis = [p['phi'] for p in disk.pillars]
    for pillar in disk.pillars:
        ir_p = np.argmin(np.abs(disk.r - pillar['r']))
        pillar['phi'] = np.mod(pillar['phi'] + omega_r[ir_p] * t, 2*np.pi)

    # Get height at rotated positions
    h_t_2d = disk.get_height(r_2d, phi_t_2d)

    # Get temperature at rotated positions (needed for fallback)
    T_t_2d = disk.get_temperature(r_2d, phi_t_2d, compute_shadows=True)

    # Compute ionizing flux map at current time (with shadows)
    log_phi_2d = compute_ionizing_flux(disk, r_2d, phi_t_2d)

    # Observer occultation: pillar wall blocks the observer's line of sight
    # to shadow cells directly behind it.  Only cells that are ACTUALLY IN
    # SHADOW should be occluded — occluding non-shadow (baseline) cells
    # removes their emission without compensation and creates a spurious
    # deficit.
    #
    # Geometric model: the face wall (height h_p) at inclination i blocks
    # the view of disk cells at r > r_p out to:
    #   dr_occ = h_p × tan(i) × |cos(phi_p)|
    # This is a hard geometric cutoff, not a Gaussian.
    shadow_mask_2d = disk._compute_shadow_mask(r_2d, phi_t_2d, h_t_2d)
    is_shadowed = shadow_mask_2d < 0.99   # cells at least partially shadowed

    observer_occult = np.ones_like(log_phi_2d)

    for pillar in disk.pillars:
        r_p = pillar['r']
        phi_p = pillar['phi']          # already rotated to lab frame
        h_p = pillar.get('height', 0.01)
        sigma_phi_p = pillar.get('sigma_phi', 0.1)

        # Geometric occultation depth: how far behind the pillar the
        # observer's line of sight is blocked by the wall of height h_p.
        cos_phi_p = np.cos(phi_p)
        face_foreshort = max(0.0, -cos_phi_p) * disk.sini
        if face_foreshort < 1e-4:
            continue                   # face edge-on or on near side

        # Maximum radial extent blocked by the wall
        tani = disk.sini / (disk.cosi + 1e-10)
        dr_occ_max = h_p * tani * max(0.0, -cos_phi_p)

        # Vectorised over grid
        dr_2d_occ = r_2d - r_p
        dphi_2d_occ = (phi_t_2d - phi_p + np.pi) % (2 * np.pi) - np.pi

        # Only occult cells that are:
        #  1. Behind the pillar (r > r_p)
        #  2. Within the geometric occultation depth
        #  3. Within the pillar's azimuthal extent
        #  4. Actually in shadow (not baseline cells)
        gauss_phi = np.exp(-0.5 * (dphi_2d_occ / sigma_phi_p) ** 2)

        eligible = (is_shadowed
                    & (dr_2d_occ > 0)
                    & (dr_2d_occ < dr_occ_max)
                    & (np.abs(dphi_2d_occ) < 3 * sigma_phi_p))

        # Smooth radial ramp within the geometric extent
        radial_weight = np.where(dr_occ_max > 0,
                                 np.clip(1.0 - dr_2d_occ / (dr_occ_max + 1e-10), 0, 1),
                                 0.0)

        occ_frac = np.where(eligible,
                            face_foreshort * radial_weight * gauss_phi,
                            0.0)
        occ_frac = np.clip(occ_frac, 0.0, 1.0)
        observer_occult *= (1.0 - occ_frac)

        if itime % 10 == 0:
            n_occ = int(np.sum(occ_frac > 0.01))
            max_occ = float(np.max(occ_frac))
            print(f"  Occultation: phi_p={np.degrees(phi_p):.1f}deg, "
                  f"foreshort={face_foreshort:.3f}, dr_occ={dr_occ_max:.3f}, "
                  f"max_occ={max_occ:.3f}, n_cells(>1%)={n_occ}")

    # Restore original pillar positions for the explicit-pillar loop
    # below, which tracks rotation independently.
    for i, pillar in enumerate(disk.pillars):
        pillar['phi'] = original_pillar_phis[i]
    
    # Interpolate Cloudy models to get line intensities
    line_intensity_2d = {}
    for line_key in lambda0_dict.keys():
        # Map line_key to Cloudy key
        cloudy_key = line_key
        if line_key == 'Halpha':
            cloudy_key = 'Halpha'
        elif line_key == 'Mg2':
            cloudy_key = 'Mg2'
        elif line_key == 'C4':
            cloudy_key = 'C4'
        
        if cloudy_key not in cloudy_interp_dict:
            continue
        
        # Interpolate at each point
        intensity_2d = np.zeros_like(log_phi_2d)
        for ir in range(disk.nr):
            for iphi in range(disk.nphi):
                log_phi_val = log_phi_2d[ir, iphi]
                log_phi_val = np.clip(log_phi_val, phi_grid[0], phi_grid[-1])
                intensity_2d[ir, iphi] = float(cloudy_interp_dict[cloudy_key](log_phi_val, Z_target, grid=False))
        
        line_intensity_2d[line_key] = intensity_2d
    
    # Initialize flux slices for each line
    flux_slice_dict = {}
    for line_key in lambda0_dict.keys():
        flux_slice_dict[line_key] = np.zeros(nlambda)
    
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
        
        # Loop over azimuth
        for iphi, phi in enumerate(phi_t_2d[ir, :]):
            # Normal vector
            px = -sintilt[iphi] * np.cos(phi)
            pz = costilt[iphi]
            
            # Foreshortening factor
            dot = disk.ex * px + disk.ez * pz
            
            if dot <= 0:
                continue
            
            # Velocity field (Keplerian orbit)
            v_phi = v_virial_cgs * np.sqrt(r_virial / r_now)  # Keplerian velocity
            
            # Line of sight velocity component
            v_los = v_phi * np.sin(phi) * disk.sini
            
            # Process each emission line
            for line_key, lambda0 in lambda0_dict.items():
                # Doppler shift: dlambda/lambda = v_los/c
                lambda_obs = lambda0 * (1.0 + v_los / C)
                
                lambda_grid = lambda_grid_dict[line_key]
                lambda_min = lambda_grid[0]
                lambda_max = lambda_grid[-1]
                dlambda = (lambda_max - lambda_min) / (nlambda - 1)
                
                if lambda_obs < lambda_min or lambda_obs > lambda_max:
                    continue
                
                # Get line intensity from Cloudy
                if line_key in line_intensity_2d:
                    line_intensity = line_intensity_2d[line_key][ir, iphi]
                else:
                    # Fallback: use temperature weighting
                    T_loc = T_t_2d[ir, iphi]
                    line_intensity = (T_loc / disk.tv1)**2
                
                # Wavelength spread
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
                
                # Weight by area, foreshortening, line intensity, and observer occultation
                weight = da[iphi] * dot * line_intensity * observer_occult[ir, iphi]
                
                # Add contribution with Gaussian spread
                if gauss_lambda_norm > 0:
                    for ilambda in range(ilambda_min, ilambda_max + 1):
                        lambda_rel = lambda_obs - lambda_grid[ilambda]
                        gauss_lambda = np.exp(-0.5 * (lambda_rel / siglambda)**2) / gauss_lambda_norm
                        flux_slice_dict[line_key][ilambda] += weight * gauss_lambda
    
    # Add emission from pillars (distributed across azimuthal extent).
    #
    # Physical model: the pillar is an opaque wall on the disk.  It
    # REPLACES a patch of normal disk surface with its visible face.
    #
    #   ΔF = (face_emission − baseline_disk_emission) × foreshortening × area
    #
    # Foreshortening of the radial face:  |cos(φ) · sin(i)|
    #   • cos(φ)>0  (near side, v≈0) → observer sees outward/shadow face
    #   • cos(φ)<0  (far  side, v≈0) → observer sees inward/illuminated face
    #   • cos(φ)≈0  (quadrature, |v|=max) → face is edge-on, ΔF→0
    #
    # Shadow face has LOW ionizing flux → line emission BELOW baseline → ΔF<0
    # Illuminated face has HIGH flux   → line emission ABOVE baseline → ΔF>0
    if len(disk.pillars) > 0:
        # Temporarily rotate pillars to current lab-frame positions
        for j, p in enumerate(disk.pillars):
            ir_j = np.argmin(np.abs(disk.r - p['r']))
            p['phi'] = np.mod(original_pillar_phis[j] + omega_r[ir_j] * t, 2*np.pi)

        for pillar in disk.pillars:
            r_pillar = pillar['r']
            phi_pillar_t = pillar['phi']  # Already rotated
            sigma_r_p = pillar.get('sigma_r', 2.0)
            sigma_phi_p = pillar.get('sigma_phi', 0.2)
            h_pillar_p = pillar.get('height', 0.01)

            # Baseline ionizing flux at r_pillar (smooth disk, no pillar)
            T_visc_rp = np.interp(r_pillar, disk.r, disk.tv_base)
            T_irrad_rp = np.interp(r_pillar, disk.r, disk.tx_base)
            T_base_rp = (T_visc_rp**4 + T_irrad_rp**4)**0.25
            T_ref = disk.tv1 * disk.tirrad_tvisc_ratio
            T_ratio_base = np.clip(T_base_rp / T_ref, 0.1, 10.0)
            T_ratio_norm = (T_ratio_base - 0.1) / (10.0 - 0.1)
            log_phi_baseline = 17.0 + T_ratio_norm * (21.5 - 17.0)
            log_phi_baseline = np.clip(log_phi_baseline, phi_grid[0], phi_grid[-1])

            # Geometric irradiation enhancement for illuminated face
            # The pillar face is roughly perpendicular to the lamp direction,
            # so it intercepts much more flux per unit area than the nearly
            # face-on flat disk at the same radius.
            # cos(theta_face) ~ r_p / d  (face normal points radially toward lamp)
            # cos(theta_disk) ~ (hlamp - h_disk) / d  (disk normal points up)
            h_pillar_val = float(np.asarray(disk.get_height(
                np.atleast_1d(r_pillar), np.atleast_1d(phi_pillar_t))).flat[0])
            d_lamp = np.sqrt(r_pillar**2 + (disk.hlamp - h_pillar_val)**2)
            cos_face = r_pillar / d_lamp  # face perpendicular to radial direction
            h_disk_base = np.interp(r_pillar, disk.r, disk.h_base)
            d_disk_base = np.sqrt(r_pillar**2 + (disk.hlamp - h_disk_base)**2)
            cos_disk = (disk.hlamp - h_disk_base) / d_disk_base  # flat disk
            irrad_enhancement = cos_face / (cos_disk + 1e-10)
            log_phi_illuminated = np.clip(
                log_phi_baseline + np.log10(max(irrad_enhancement, 1.0)),
                phi_grid[0], phi_grid[-1])
            # Shadow face: only viscous temperature (no irradiation)
            T_ratio_shadow = np.clip(T_visc_rp / T_ref, 0.1, 10.0)
            T_ratio_norm_shadow = (T_ratio_shadow - 0.1) / (10.0 - 0.1)
            log_phi_shadow = np.clip(
                17.0 + T_ratio_norm_shadow * (21.5 - 17.0),
                phi_grid[0], phi_grid[-1])

            # Keplerian velocity at pillar radius
            v_phi_pillar = v_virial_cgs * np.sqrt(r_virial / r_pillar)

            # Sample the pillar across ±3 sigma_phi
            n_phi_samples = max(20, int(6 * sigma_phi_p / 0.05))
            n_phi_samples = min(n_phi_samples, 200)
            dphi_arr = np.linspace(-3 * sigma_phi_p, 3 * sigma_phi_p, n_phi_samples)
            delta_dphi = dphi_arr[1] - dphi_arr[0] if n_phi_samples > 1 else 6 * sigma_phi_p

            # Debug output for first pillar at selected time steps
            if itime % 10 == 0 and pillar == disk.pillars[0]:
                cos_phi_c = np.cos(phi_pillar_t)
                foreshort_c = np.abs(cos_phi_c * disk.sini)
                seeing = "shadow" if cos_phi_c * disk.sini > 0 else "illuminated"
                print(f"  t={t:.1f}d, phi_t={np.degrees(phi_pillar_t):.1f}deg, "
                      f"foreshort={foreshort_c:.3f}, seeing={seeing}, "
                      f"log_phi_base={log_phi_baseline:.2f}, "
                      f"log_phi_illum={log_phi_illuminated:.2f} (enhance={irrad_enhancement:.1f}x), "
                      f"n_samples={n_phi_samples}")

            for dphi_s in dphi_arr:
                phi_s = phi_pillar_t + dphi_s
                gauss_weight = np.exp(-0.5 * (dphi_s / sigma_phi_p)**2)

                # Foreshortening of the radial face: |cos(φ) · sin(i)|
                # Sign of cos(φ)·sin(i) determines which face is visible:
                #   > 0 → observer sees outward (shadow) face
                #   < 0 → observer sees inward (illuminated) face
                # Use a smooth blend so the transition through quadrature
                # (edge-on, foreshortening=0) is gradual, not a step.
                face_dot_obs = np.cos(phi_s) * disk.sini
                foreshortening = np.abs(face_dot_obs)

                # Smooth blend: 0 = illuminated face, 1 = shadow face
                # tanh gives a smooth sigmoid; width ~0.05 in face_dot_obs
                blend = 0.5 * (1.0 + np.tanh(face_dot_obs / 0.05))
                log_phi_face = (1.0 - blend) * log_phi_illuminated + blend * log_phi_shadow

                # Area element: wall height × azimuthal arc × Gaussian weight
                # The face is a vertical wall of height h_pillar (not sigma_r).
                area_s = h_pillar_p * r_pillar * delta_dphi * gauss_weight

                # LOS velocity at this azimuth
                v_los_s = v_phi_pillar * np.sin(phi_s) * disk.sini

                for line_key, lambda0 in lambda0_dict.items():
                    lambda_obs_s = lambda0 * (1.0 + v_los_s / C)

                    lambda_grid_line = lambda_grid_dict[line_key]
                    lambda_min = lambda_grid_line[0]
                    lambda_max = lambda_grid_line[-1]
                    dlambda = (lambda_max - lambda_min) / (nlambda - 1)

                    if lambda_obs_s < lambda_min or lambda_obs_s > lambda_max:
                        continue

                    cloudy_key = line_key
                    if cloudy_key in cloudy_interp_dict:
                        face_intensity = float(cloudy_interp_dict[cloudy_key](log_phi_face, Z_target, grid=False))
                        shadow_intensity = float(cloudy_interp_dict[cloudy_key](log_phi_shadow, Z_target, grid=False))
                        # The face replaces the view of disk cells behind
                        # the pillar.  Those cells are always in shadow
                        # (the shadow extends radially outward from the
                        # pillar regardless of observer direction).
                        replaced_intensity = shadow_intensity
                    else:
                        face_intensity = 10.0
                        replaced_intensity = 0.0

                    # Net contribution: face emission minus what it covers.
                    # Far side: illuminated face − shadow ≈ large positive
                    # Near side: shadow face − shadow ≈ 0 (no net change)
                    delta_intensity = face_intensity - replaced_intensity
                    pillar_weight = area_s * delta_intensity * foreshortening

                    siglambda_s = dlambda
                    nsigma_lambda = 2
                    il_min = max(0, int((lambda_obs_s - lambda_min - nsigma_lambda * siglambda_s) / dlambda))
                    il_max = min(nlambda - 1, int((lambda_obs_s - lambda_min + nsigma_lambda * siglambda_s) / dlambda))

                    gauss_norm = 0.0
                    for il in range(il_min, il_max + 1):
                        gauss_norm += np.exp(-0.5 * ((lambda_obs_s - lambda_grid_line[il]) / siglambda_s)**2)

                    if gauss_norm > 0:
                        for ilambda in range(il_min, il_max + 1):
                            gl = np.exp(-0.5 * ((lambda_obs_s - lambda_grid_line[ilambda]) / siglambda_s)**2) / gauss_norm
                            flux_slice_dict[line_key][ilambda] += pillar_weight * gl

        # Restore original pillar positions
        for j, p in enumerate(disk.pillars):
            p['phi'] = original_pillar_phis[j]

    return itime, flux_slice_dict


def compute_time_evolving_cloudy_map(disk, cloudy_interp_dict, phi_grid, lambda0_dict, 
                                     v_virial, nlambda=200, ntime=50, tmax=3000.0, 
                                     lambda_range_dict=None, use_absolute_flux=True, Z_target=1.0,
                                     n_cores=None, normalize_output=False):
    """
    Compute time-evolving velocity map with Cloudy line intensities: wavelength vs time.
    
    Parameters:
    -----------
    disk : PillarDisk
        Disk model with geometry
    cloudy_interp_dict : dict
        Dictionary of Cloudy interpolation functions
    phi_grid : array
        Ionizing flux grid (log phi)
    lambda0_dict : dict
        Dictionary of rest wavelengths: {'Halpha': 6562.8, 'Mg2': 2798.0, 'C4': 1549.0}
    v_virial : float
        FWHM of emission line (km/s). The actual Keplerian velocity is
        v(r) = v_virial * sqrt(r0/r), where r0 is the reference radius.
    nlambda : int
        Number of wavelength bins per line
    ntime : int
        Number of time bins
    tmax : float or None
        Maximum time (days). If None, estimate from orbital period at r_out
    lambda_range_dict : dict or None
        Dictionary of (min, max) wavelength ranges per line. If None, auto-compute
    use_absolute_flux : bool
        Whether to use absolute flux (True) or Hbeta-normalized flux (False)
    Z_target : float
        Target metallicity for Cloudy interpolation (default: 1.0)
    
    Returns:
    --------
    lambda_grid_dict : dict
        Dictionary of wavelength grids per line
    time_grid : array
        Time grid (days)
    flux_map_dict : dict
        Dictionary of flux maps per line: {'Halpha': map, 'Mg2': map, 'C4': map}
        Each map has shape (ntime, nlambda)
    """
    # Create 2D grids for disk surface
    r_2d, phi_2d = np.meshgrid(disk.r, disk.phi, indexing='ij')
    
    # Get height (with shadows)
    h_2d = disk.get_height(r_2d, phi_2d)
    
    # Wavelength grids for each line
    lambda_grid_dict = {}
    if lambda_range_dict is None:
        lambda_range_dict = {}
    
    for line_key, lambda0 in lambda0_dict.items():
        if line_key in lambda_range_dict:
            lambda_min, lambda_max = lambda_range_dict[line_key]
        else:
            v_range = 3.0 * v_virial  # velocity range in km/s
            v_c = v_range * KM_TO_CM / C  # velocity in units of c
            dlambda = lambda0 * v_c
            lambda_min = lambda0 - dlambda
            lambda_max = lambda0 + dlambda
        
        lambda_grid = np.linspace(lambda_min, lambda_max, nlambda)
        lambda_grid_dict[line_key] = lambda_grid
    
    # Time grid
    if tmax is None:
        r_virial_loc = disk.r0
        if len(disk.pillars) > 0:
            pillar_radii = [pillar['r'] for pillar in disk.pillars]
            r_pillar_mean = np.mean(pillar_radii)
            T_orbital = compute_orbital_period(r_pillar_mean, v_virial, r_virial_loc)
        else:
            T_orbital = compute_orbital_period(disk.rout, v_virial, r_virial_loc)
        tmax = 2.0 * T_orbital
        print(f"Auto-computed tmax = {tmax:.1f} days (2 × T_orbital = {T_orbital:.1f} days)")

    time_grid = np.linspace(0, tmax, ntime, endpoint=False)
    dtime = tmax / ntime
    
    # Initialize flux maps for each line
    flux_map_dict = {}
    for line_key in lambda0_dict.keys():
        flux_map_dict[line_key] = np.zeros((ntime, nlambda))
    
    # Convert virial velocity to cm/s
    v_virial_cgs = v_virial * KM_TO_CM  # km/s to cm/s
    r_virial = disk.r0  # Reference radius for virial velocity
    
    # Compute angular velocities for each radius
    omega_r = compute_angular_velocity(disk.r, v_virial, r_virial)  # rad/day
    
    # Initial azimuthal positions (stored as phi_0)
    phi_0_2d = phi_2d.copy()
    
    print("Computing time-evolving velocity map with Cloudy models...")
    print(f"Time range: 0 to {tmax} days")
    print(f"Lines: {list(lambda0_dict.keys())}")
    
    # Prepare disk parameters for multiprocessing (need to serialize disk state)
    # Note: pillars are stored separately and added after disk creation
    disk_params = {
        'rin': disk.rin, 'rout': disk.rout, 'nr': disk.nr, 'nphi': disk.nphi,
        'h1': disk.h1, 'r0': disk.r0, 'beta': disk.beta,
        'tv1': disk.tv1, 'alpha': disk.alpha, 'hlamp': disk.hlamp,
        'tx1': disk.tx1, 'tirrad_tvisc_ratio': disk.tirrad_tvisc_ratio,
        'dmpc': disk.dmpc, 'cosi': disk.cosi, 'fcol': disk.fcol, 'redshift': disk.redshift
    }
    # Store pillars separately (they'll be added after disk creation in worker)
    pillars_list = disk.pillars.copy()
    
    # Determine number of cores to use
    if n_cores is None:
        n_cores = max(1, cpu_count() - 1)  # Use all cores except one
    n_cores = min(n_cores, ntime)  # Don't use more cores than time steps
    
    print(f"Using {n_cores} CPU cores for parallel processing")
    
    # Prepare arguments for multiprocessing
    process_args = []
    for itime, t in enumerate(time_grid):
        args = (itime, t, disk_params, pillars_list, cloudy_interp_dict, phi_grid, lambda0_dict,
                v_virial_cgs, r_virial, omega_r, r_2d, phi_0_2d, lambda_grid_dict,
                nlambda, Z_target)
        process_args.append(args)
    
    # Process time steps in parallel
    if n_cores > 1:
        with Pool(processes=n_cores) as pool:
            results = list(tqdm(pool.imap(_process_time_step, process_args), 
                              total=ntime, desc="Processing time", unit="time"))
    else:
        # Serial processing (for debugging or single core)
        results = []
        for args in tqdm(process_args, desc="Processing time", unit="time"):
            results.append(_process_time_step(args))
    
    # Collect results
    print("Collecting results from parallel workers...")
    for itime, flux_slice_dict in results:
        for line_key in lambda0_dict.keys():
            flux_map_dict[line_key][itime, :] = flux_slice_dict[line_key]
    print("Results collected. Normalizing maps...")
    
    # Optional: normalize each map (dimensionless), so integral over (t, lambda) ~ 1
    # Default is False so the map reflects the (disk-integrated) Cloudy-weighted flux scale.
    if normalize_output:
        print("Normalizing maps...")
        for line_key in flux_map_dict.keys():
            dlambda = lambda_grid_dict[line_key][1] - lambda_grid_dict[line_key][0]
            total = np.sum(flux_map_dict[line_key]) * dtime * dlambda
            if total > 0:
                flux_map_dict[line_key] = flux_map_dict[line_key] / total
                print(f"  Normalized {line_key}: total = {total:.2e}")
            else:
                print(f"Warning: Total flux is zero for {line_key}. Check disk geometry and parameters.")
        print("Normalization complete.")
    return lambda_grid_dict, time_grid, flux_map_dict


def plot_time_evolving_cloudy_map(lambda_grid_dict, time_grid, flux_map_dict, lambda0_dict, 
                                  filename='velocity_time_cloudy.png', use_absolute_flux=False,
                                  normalize_output=False):
    """
    Plot the time-evolving velocity maps for multiple lines: wavelength vs time.
    
    Parameters:
    -----------
    lambda_grid_dict : dict
        Dictionary of wavelength grids per line
    time_grid : array
        Time grid (days)
    flux_map_dict : dict
        Dictionary of flux maps per line
    lambda0_dict : dict
        Dictionary of rest wavelengths per line
    filename : str
        Output filename
    """
    print(f"  Creating figure with {len(lambda_grid_dict)} panels...")
    n_lines = len(lambda_grid_dict)
    fig, axes = plt.subplots(n_lines, 1, figsize=(14, 5*n_lines))
    
    if n_lines == 1:
        axes = [axes]
    
    line_colors = {'Halpha': 'red', 'Mg2': 'green', 'C4': 'blue'}
    line_labels = {'Halpha': r'H\alpha', 'Mg2': r'MgII', 'C4': r'CIV'}

    # Colorbar label: we just call it "Flux" regardless of absolute vs Hβ-normalized
    cbar_label = r'$\rm Flux$'
    
    for idx, (line_key, lambda_grid) in enumerate(lambda_grid_dict.items()):
        ax = axes[idx]
        flux_map = flux_map_dict[line_key]
        lambda0 = lambda0_dict[line_key]
        
        # Use log scale for better visualization
        flux_plot = flux_map.copy()
        flux_plot[flux_plot <= 0] = np.nan
        
        im = ax.pcolormesh(lambda_grid, time_grid, flux_map, shading='auto', cmap='rainbow')
        ax.axvline(lambda0, color='white', linestyle='--', linewidth=2, alpha=0.7)
        ax.set_xlabel(r'$\rm Wavelength~[\AA]$', fontsize=14)
        ax.set_ylabel(r'$\rm Time~[days]$', fontsize=14)
        # Panel title: just the line name (no repeated "Time-Evolving Velocity Map")
        line_label = line_labels.get(line_key, line_key)
        ax.set_title(rf'$\rm {line_label}$', fontsize=16, pad=10)
        ax.set_xlim(lambda_grid[0], lambda_grid[-1])
        ax.set_ylim(time_grid[0], time_grid[-1])
        plt.colorbar(im, ax=ax, label=cbar_label)
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
    print(f"Time-evolving velocity maps saved to {filename}")
    plt.close()


def plot_time_evolving_diff_map(lambda_grid_dict, time_grid,
                                flux_with_pillars, flux_without_pillars,
                                lambda0_dict,
                                filename='velocity_time_cloudy_diff.png'):
    """
    Plot the difference time-evolving velocity maps: (with pillars - without pillars).

    This shows the contribution of pillars to the emission.

    Parameters:
    -----------
    lambda_grid_dict : dict
        Dictionary of wavelength grids per line
    time_grid : array
        Time grid (days)
    flux_with_pillars : dict
        Flux maps including pillars
    flux_without_pillars : dict
        Flux maps without pillars
    lambda0_dict : dict
        Dictionary of rest wavelengths per line
    filename : str
        Output filename
    """
    print(f"  Creating difference figure with {len(lambda_grid_dict)} panels...")
    n_lines = len(lambda_grid_dict)
    fig, axes = plt.subplots(n_lines, 1, figsize=(14, 5*n_lines))

    if n_lines == 1:
        axes = [axes]

    line_labels = {'Halpha': r'H\alpha', 'Mg2': r'MgII', 'C4': r'CIV'}

    for idx, (line_key, lambda_grid) in enumerate(lambda_grid_dict.items()):
        ax = axes[idx]
        flux_with = flux_with_pillars[line_key]
        flux_without = flux_without_pillars[line_key]

        diff = flux_with - flux_without

        # Choose symmetric color scale around zero
        vmax = np.nanmax(np.abs(diff))
        if not np.isfinite(vmax) or vmax == 0:
            vmax = 1.0
        vmin = -vmax

        im = ax.pcolormesh(lambda_grid, time_grid, diff, shading='auto',
                           cmap='seismic', vmin=vmin, vmax=vmax)
        lambda0 = lambda0_dict[line_key]
        ax.axvline(lambda0, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
        ax.set_xlabel(r'$\rm Wavelength~[\AA]$', fontsize=14)
        ax.set_ylabel(r'$\rm Time~[days]$', fontsize=14)
        line_label = line_labels.get(line_key, line_key)
        ax.set_title(rf'$\rm {line_label}~(with~pillars - without~pillars)$', fontsize=16, pad=10)
        ax.set_xlim(lambda_grid[0], lambda_grid[-1])
        ax.set_ylim(time_grid[0], time_grid[-1])
        plt.colorbar(im, ax=ax, label=r'$\rm \Delta Flux$')
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
    print(f"Difference time-evolving velocity maps saved to {filename}")
    plt.close()


def plot_time_evolving_residual_map(lambda_grid_dict, time_grid, flux_map_dict, lambda0_dict,
                                    filename='velocity_time_cloudy_residual.png'):
    """
    Plot the true residual time-evolving velocity maps: flux(t) - <flux>_t.

    The residual is computed by subtracting the time-averaged spectrum at each wavelength.
    This reveals time-variable features like the barber-pole pattern from rotating pillars.

    Parameters:
    -----------
    lambda_grid_dict : dict
        Dictionary of wavelength grids per line
    time_grid : array
        Time grid (days)
    flux_map_dict : dict
        Dictionary of flux maps per line (ntime x nlambda)
    lambda0_dict : dict
        Dictionary of rest wavelengths per line
    filename : str
        Output filename
    """
    print(f"  Creating residual (time-averaged subtracted) figure with {len(lambda_grid_dict)} panels...")
    n_lines = len(lambda_grid_dict)
    fig, axes = plt.subplots(n_lines, 1, figsize=(14, 5*n_lines))

    if n_lines == 1:
        axes = [axes]

    line_labels = {'Halpha': r'H\alpha', 'Mg2': r'MgII', 'C4': r'CIV'}

    for idx, (line_key, lambda_grid) in enumerate(lambda_grid_dict.items()):
        ax = axes[idx]
        flux_map = flux_map_dict[line_key]  # Shape: (ntime, nlambda)

        # Compute time-averaged spectrum at each wavelength
        time_avg = np.mean(flux_map, axis=0)  # Shape: (nlambda,)

        # Compute residual: flux(t, lambda) - <flux(lambda)>_t
        residual = flux_map - time_avg[np.newaxis, :]  # Broadcasting: (ntime, nlambda)

        # Choose symmetric color scale around zero
        vmax = np.nanmax(np.abs(residual))
        if not np.isfinite(vmax) or vmax == 0:
            vmax = 1.0
        vmin = -vmax

        im = ax.pcolormesh(lambda_grid, time_grid, residual, shading='auto',
                           cmap='seismic', vmin=vmin, vmax=vmax)
        lambda0 = lambda0_dict[line_key]
        ax.axvline(lambda0, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
        ax.set_xlabel(r'$\rm Wavelength~[\AA]$', fontsize=14)
        ax.set_ylabel(r'$\rm Time~[days]$', fontsize=14)
        line_label = line_labels.get(line_key, line_key)
        ax.set_title(rf'$\rm {line_label}~Residual~(F(t,\lambda) - \langle F(\lambda) \rangle_t)$', fontsize=16, pad=10)
        ax.set_xlim(lambda_grid[0], lambda_grid[-1])
        ax.set_ylim(time_grid[0], time_grid[-1])
        plt.colorbar(im, ax=ax, label=r'$\rm Residual~Flux$')
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
    print(f"Residual (time-averaged subtracted) maps saved to {filename}")
    plt.close()


def main(config_file='config_line.yaml', use_absolute_flux=False):
    """
    Main function to compute and plot time-evolving velocity maps with Cloudy models.
    
    Parameters:
    -----------
    config_file : str
        Path to configuration YAML file
    use_absolute_flux : bool
        Whether to use absolute flux (True) or Hbeta-normalized flux (False)
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
    
    # Process pillars (same logic as pillar_line_time.py)
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

    # ---- Geometry visualization at t=0 ----
    print("\n--- Geometry check at t=0 ---")

    # 3D geometry plot
    disk.plot_3d_geometry(show_light_rays=True, filename='geometry_t0_3d.png')
    print("  Saved 3D geometry: geometry_t0_3d.png")

    # Face-on ionizing flux + height map (r-phi plane)
    # Build radial grid with refinement around each pillar so that narrow
    # pillars (sigma_r << global grid spacing) are properly resolved.
    nphi_plot = 360
    r_coarse = np.linspace(disk.rin, disk.rout, 200)
    r_refined_parts = [r_coarse]
    for p in disk.pillars:
        sr = max(p['sigma_r'], 0.01)
        r_lo = max(disk.rin, p['r'] - 5*sr)
        r_hi = min(disk.rout, p['r'] + 5*sr)
        r_refined_parts.append(np.linspace(r_lo, r_hi, 60))
    r_plot = np.unique(np.concatenate(r_refined_parts))
    phi_plot = np.linspace(0, 2*np.pi, nphi_plot)
    r_2d_plot, phi_2d_plot = np.meshgrid(r_plot, phi_plot, indexing='ij')

    h_2d_plot = disk.get_height(r_2d_plot, phi_2d_plot)
    log_phi_2d_plot = compute_ionizing_flux(disk, r_2d_plot, phi_2d_plot)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5), subplot_kw={'projection': 'polar'})

    ax0 = axes[0]
    c0 = ax0.pcolormesh(phi_2d_plot, r_2d_plot, h_2d_plot, shading='auto', cmap='terrain')
    ax0.set_title(r'$\rm Disk~Height$', pad=15)
    fig.colorbar(c0, ax=ax0, pad=0.1)
    for p in disk.pillars:
        ax0.plot(p['phi'], p['r'], 'r*', markersize=5)

    ax1 = axes[1]
    c1 = ax1.pcolormesh(phi_2d_plot, r_2d_plot, log_phi_2d_plot, shading='auto', cmap='inferno')
    ax1.set_title(r'$\rm Ionizing~Flux$', pad=15)
    fig.colorbar(c1, ax=ax1, pad=0.1)
    for p in disk.pillars:
        ax1.plot(p['phi'], p['r'], 'c*', markersize=5)

    plt.tight_layout()
    plt.savefig('geometry_t0_faceon.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved face-on map: geometry_t0_faceon.png")
    print("--- End geometry check ---\n")

    # Load Cloudy models
    if use_absolute_flux:
        cloudy_file = '/Users/jiamuh/c23.01/my_models/loc_metal_flux/strong_LOC_varym_N25_v100_lineflux_LineList_BLR_Fe2_flux.txt'
    else:
        cloudy_file = '/Users/jiamuh/c23.01/my_models/loc_metal/strong_LOC_varym_N25_LineList_BLR_Fe2.txt'
    
    print(f"\nLoading Cloudy models from: {cloudy_file}")
    cloudy_interp_dict, phi_grid, Z_grid = load_cloudy_models(cloudy_file, Z_target=1.0)
    
    # Emission line parameters
    lambda0_dict = {
        'Halpha': 6562.8,   # Halpha rest wavelength
        'Mg2': 2798.0,      # MgII rest wavelength
        'C4': 1549.0        # CIV rest wavelength
    }
    
    v_virial = 1000.0  # Virial velocity (km/s)
    
    line_params = config.get('emission_line', {})
    if line_params:
        v_virial = float(line_params.get('v_virial', v_virial))
        # Allow custom lambda0 per line
        if 'lambda0_Halpha' in line_params:
            lambda0_dict['Halpha'] = float(line_params['lambda0_Halpha'])
        if 'lambda0_Mg2' in line_params:
            lambda0_dict['Mg2'] = float(line_params['lambda0_Mg2'])
        if 'lambda0_C4' in line_params:
            lambda0_dict['C4'] = float(line_params['lambda0_C4'])
    
    # Computation parameters
    comp_params = config.get('computation', {})
    nlambda = int(comp_params.get('nlambda', 200))
    ntime = int(comp_params.get('ntime', 50))  # Default: 50 steps
    tmax_config = comp_params.get('tmax', None)
    tmax = float(tmax_config) if tmax_config is not None else None
    n_cores = comp_params.get('n_cores', None)  # None = auto-detect
    if n_cores is not None:
        n_cores = int(n_cores)

    # Round tmax to integer orbital periods for exact periodicity
    if len(disk.pillars) > 0:
        pillar_radii = [p['r'] for p in disk.pillars]
        r_pillar_mean = np.mean(pillar_radii)
        T_orbital = compute_orbital_period(r_pillar_mean, v_virial, disk.r0)
        if tmax is not None:
            n_orbits = max(1, round(tmax / T_orbital))
            tmax_orig = tmax
            tmax = n_orbits * T_orbital
            print(f"Adjusted tmax: {tmax_orig:.1f} → {tmax:.1f} days "
                  f"({n_orbits} × T_orbital={T_orbital:.1f}d at r={r_pillar_mean:.1f})")
        else:
            tmax = 2.0 * T_orbital
            print(f"Auto tmax = {tmax:.1f} days (2 × T_orbital={T_orbital:.1f}d)")
    elif tmax is None:
        tmax = 3000.0

    print(f"\nComputing time-evolving velocity maps with Cloudy models:")
    print(f"  Lines: {list(lambda0_dict.keys())}")
    print(f"  Rest wavelengths: {lambda0_dict}")
    print(f"  Virial velocity: {v_virial:.0f} km/s")
    print(f"  Number of pillars: {len(disk.pillars)}")
    
    # Compute time-evolving maps: with and without pillars
    print("\nStarting computation with pillars...")
    lambda_grid_dict, time_grid, flux_with_pillars = compute_time_evolving_cloudy_map(
        disk, cloudy_interp_dict, phi_grid, lambda0_dict, v_virial, 
        nlambda=nlambda, ntime=ntime, tmax=tmax, use_absolute_flux=use_absolute_flux,
        n_cores=n_cores)
    
    print("Starting computation without pillars...")
    saved_pillars = disk.pillars.copy()
    disk.pillars = []  # temporarily remove pillars
    lambda_grid_dict_no, time_grid_no, flux_without_pillars = compute_time_evolving_cloudy_map(
        disk, cloudy_interp_dict, phi_grid, lambda0_dict, v_virial, 
        nlambda=nlambda, ntime=ntime, tmax=tmax, use_absolute_flux=use_absolute_flux,
        n_cores=n_cores)
    disk.pillars = saved_pillars  # restore pillars
    
    # Sanity check: grids should match
    # (we reuse lambda_grid_dict and time_grid for plotting)
    
    print("Computation complete. Creating plots...")
    # Plot main maps (with pillars)
    plot_params = config.get('plotting', {})
    filename = plot_params.get('velocity_time_cloudy_filename', 'velocity_time_cloudy.png')
    plot_time_evolving_cloudy_map(
        lambda_grid_dict, time_grid, flux_with_pillars, lambda0_dict,
        filename=filename, use_absolute_flux=use_absolute_flux,
        normalize_output=False
    )
    
    # Plot difference maps (with pillars - without pillars)
    diff_filename = plot_params.get('velocity_time_cloudy_diff_filename',
                                    'velocity_time_cloudy_diff.png')
    plot_time_evolving_diff_map(
        lambda_grid_dict, time_grid,
        flux_with_pillars, flux_without_pillars,
        lambda0_dict,
        filename=diff_filename
    )

    # Plot true residual maps (flux - time-averaged flux)
    residual_filename = plot_params.get('velocity_time_cloudy_residual_filename',
                                        'velocity_time_cloudy_residual.png')
    plot_time_evolving_residual_map(
        lambda_grid_dict, time_grid,
        flux_with_pillars, lambda0_dict,
        filename=residual_filename
    )
    print("Plotting complete.")
    
    print("\nSummary:")
    for line_key, lambda_grid in lambda_grid_dict.items():
        print(f"{line_key}:")
        print(f"  Wavelength range: {lambda_grid[0]:.2f} - {lambda_grid[-1]:.2f} Å")
        print(f"  Max flux (with pillars): {np.max(flux_with_pillars[line_key]):.2e}")
        print(f"  Max flux (without pillars): {np.max(flux_without_pillars[line_key]):.2e}")
    print(f"Time range: {time_grid[0]:.2f} - {time_grid[-1]:.2f} days")


if __name__ == '__main__':
    import sys
    config_file = sys.argv[1] if len(sys.argv) > 1 else 'config_line.yaml'
    use_absolute_flux = len(sys.argv) > 2 and sys.argv[2].lower() in ('true', '1', 'yes')
    main(config_file=config_file, use_absolute_flux=use_absolute_flux)
