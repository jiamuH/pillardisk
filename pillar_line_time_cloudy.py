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
from matplotlib.colors import LinearSegmentedColormap
import os
import sys
from multiprocessing import Pool, cpu_count

VDM_CMAP = LinearSegmentedColormap.from_list('vdm_wrb',
                                              ['white', 'red', 'blue'])

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

# Import Cloudy functions from pillar_line_cloudy.py
try:
    from pillardisk.pillar_line_cloudy import load_cloudy_models, compute_ionizing_flux
except ImportError as e:
    print(f"Error: pillar_line_cloudy.py not found. Make sure it's in the same directory.")
    print(f"Import error: {e}")
    sys.exit(1)

# Import velocity-delay map computation from pillar_line.py
try:
    from pillardisk.pillar_line import compute_velocity_delay_map
except ImportError as e:
    print(f"Error: pillar_line.py not found. Make sure it's in the same directory.")
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
M_SUN = 1.989e33 # solar mass (g)
PC = 3.0857e18   # parsec (cm)
DAY = 24. * 3600.  # day (s)
ANGSTROM = 1e-8   # Angstrom (cm)
KM_TO_CM = 1e5    # km to cm conversion

# Conversion factors
PC_TO_LD = 1e6 * PC / (C * DAY)  # parsec to light days
LD_TO_CM = C * DAY  # light days to cm


def v_virial_from_mbh(M_BH_solar, r0_ld):
    """Compute virial velocity at reference radius r0 from SMBH mass.

    v_vir = sqrt(G * M_BH / r0)

    Parameters
    ----------
    M_BH_solar : float
        Black hole mass in solar masses
    r0_ld : float
        Reference radius in light-days

    Returns
    -------
    v_virial : float
        Virial velocity in km/s
    """
    M_cgs = M_BH_solar * M_SUN
    r0_cgs = r0_ld * LD_TO_CM
    v_cgs = np.sqrt(G * M_cgs / r0_cgs)
    return v_cgs / KM_TO_CM


def compute_velocity_delay_map_cloudy(disk, line_key, lambda0, v_virial,
                                       cloudy_interp_dict, Z_target=1.0,
                                       nlambda=200, ntau=100, taumax=150.0,
                                       phi_grid_range=None,
                                       weighting='emissivity'):
    """
    Compute Cloudy-weighted velocity-delay map Ψ(λ, τ) for an emission line.

    ``weighting`` selects how each cell is weighted:

    * ``'emissivity'`` (default): per-cell weight is
      F_line(log Φ_H), i.e. the steady-state line surface flux. The
      resulting map is the steady-state line emission distribution in
      (λ, τ) space.
    * ``'responsivity'``: per-cell weight is
      F_line · η, where η ≡ ∂log F_line / ∂log Φ_H is the local
      logarithmic slope of the F–L curve, computed by finite differencing
      the Cloudy spline. This is the linear response of the line flux to
      a small fractional perturbation in the lamp luminosity (the
      Blandford-McKee response function) for cells whose F–L curve is
      smooth at the local log Φ_H.

    Parameters
    ----------
    disk : PillarDisk
        Disk model with geometry.
    line_key : str
        Line name for Cloudy lookup ('Halpha', 'Mg2', 'C4', etc.).
    lambda0 : float
        Rest wavelength of emission line (Angstroms).
    v_virial : float
        FWHM of emission line (km/s).
    cloudy_interp_dict : dict
        Dictionary of Cloudy interpolation functions from
        ``load_cloudy_models``.
    Z_target : float
        Target metallicity (default: 1.0 for solar).
    nlambda, ntau : int
        Grid sizes.
    taumax : float
        Maximum delay (days).
    phi_grid_range : tuple or None
        (phi_min, phi_max) bounds for the Cloudy ionizing-flux axis.
    weighting : {'emissivity', 'responsivity'}
        Per-cell weighting (see docstring above).

    Returns
    -------
    lambda_grid, tau_grid, psi_map
    """
    if weighting not in ('emissivity', 'responsivity'):
        raise ValueError(
            f"weighting must be 'emissivity' or 'responsivity', "
            f"got {weighting!r}")
    # Create 2D grids for disk surface
    r_2d, phi_2d = np.meshgrid(disk.r, disk.phi, indexing='ij')

    # Get height and temperature (with shadows)
    h_2d = disk.get_height(r_2d, phi_2d)
    T_2d = disk.get_temperature(r_2d, phi_2d, compute_shadows=True)

    # Compute ionizing flux map (log phi)
    log_phi_2d = compute_ionizing_flux(disk, r_2d, phi_2d)

    # Get Cloudy interpolation function for this line (need early for debug)
    if line_key not in cloudy_interp_dict:
        raise ValueError(f"Line '{line_key}' not found in Cloudy models. "
                        f"Available: {list(cloudy_interp_dict.keys())}")
    cloudy_interp = cloudy_interp_dict[line_key]
    # Use the actual phi-grid range from the spline rather than the legacy
    # hard-coded [17, 21]. The loader extends the lower edge to 15 when an
    # extension file is provided; honouring that range lets the VDM see
    # shadow-zone Cloudy emissivities instead of clipping them away.
    if phi_grid_range is None:
        try:
            tx = cloudy_interp.get_knots()[0]
            phi_grid_min = float(tx.min())
            phi_grid_max = float(tx.max())
        except Exception:
            phi_grid_min, phi_grid_max = 17.0, 21.0
    else:
        phi_grid_min, phi_grid_max = float(phi_grid_range[0]), float(phi_grid_range[1])
    no_fluxfloor = getattr(disk, 'no_fluxfloor', False)

    # Debug: check if shadow is reducing ionizing flux
    if len(disk.pillars) > 0:
        if os.environ.get('PILLAR_DEBUG'):
            shadow_mask = disk._compute_shadow_mask(r_2d, phi_2d, h_2d)
            n_shadowed = np.sum(shadow_mask < 0.5)
            log_phi_in_shadow = log_phi_2d[shadow_mask < 0.5]
            log_phi_not_shadow = log_phi_2d[shadow_mask >= 0.5]
            print(f"  [DEBUG Shadow] {n_shadowed} cells in shadow (line={line_key})")
            if n_shadowed > 0:
                print(f"  [DEBUG Shadow] log_phi in shadow: min={np.min(log_phi_in_shadow):.2f}, mean={np.mean(log_phi_in_shadow):.2f}")
                print(f"  [DEBUG Shadow] log_phi outside shadow: min={np.min(log_phi_not_shadow):.2f}, mean={np.mean(log_phi_not_shadow):.2f}")

    # Wavelength grid
    v_range = 3.0 * v_virial  # velocity range in km/s
    v_c = v_range * KM_TO_CM / C
    dlambda_range = lambda0 * v_c
    lambda_min = lambda0 - dlambda_range
    lambda_max = lambda0 + dlambda_range
    lambda_grid = np.linspace(lambda_min, lambda_max, nlambda)
    dlambda = (lambda_max - lambda_min) / (nlambda - 1)

    # Delay grid
    tau_grid = np.linspace(0, taumax, ntau)
    dtau = taumax / (ntau - 1)
    tau_max_reasonable = taumax * 1.1

    # Initialize response map
    psi_map = np.zeros((ntau, nlambda))

    # Convert virial velocity to cm/s
    v_virial_cgs = v_virial * KM_TO_CM
    r_virial = disk.r0

    # Loop over disk surface
    for ir in range(disk.nr - 1):
        r_now = disk.r[ir]
        h_now = h_2d[ir, :]
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

        # Get temperature and ionizing flux at this radius
        T_now = T_2d[ir, :]
        log_phi_now = log_phi_2d[ir, :]

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
            dx = r_now * np.cos(phi)
            dy = r_now * np.sin(phi)
            dz = h_now[iphi] - disk.hlamp
            dlamp = np.sqrt(dx**2 + dy**2 + dz**2)
            tau_delay = dlamp - (dx * disk.sini + h_now[iphi] * disk.cosi)

            if tau_delay < 0 or tau_delay > tau_max_reasonable:
                continue

            # Keplerian velocity
            v_phi = v_virial_cgs * np.sqrt(r_virial / r_now)
            v_los = v_phi * np.sin(phi) * disk.sini
            lambda_obs = lambda0 * (1.0 + v_los / C)

            if lambda_obs < lambda_min or lambda_obs > lambda_max:
                continue

            # Wavelength spread: turbulent broadening + numerical floor
            f_turb = getattr(disk, 'f_turb_line', 0.0)
            v_phi = v_virial_cgs * np.sqrt(r_virial / r_now)
            if f_turb > 0:
                siglambda = max(lambda0 * f_turb * v_phi / C, dlambda)
            else:
                siglambda = dlambda
            nsigma = 3 if f_turb > 0 else 2
            ilambda_min = max(0, int((lambda_obs - lambda_min - nsigma * siglambda) / dlambda))
            ilambda_max = min(nlambda - 1, int((lambda_obs - lambda_min + nsigma * siglambda) / dlambda))

            gauss_lambda_norm = 0.0
            for il in range(ilambda_min, ilambda_max + 1):
                lambda_rel = lambda_obs - lambda_grid[il]
                gauss = np.exp(-0.5 * (lambda_rel / siglambda)**2)
                gauss_lambda_norm += gauss

            # Delay bin with Gaussian spread
            sigtau = dtau
            nsigma_tau = 3
            ibin_min = max(0, int((tau_delay - nsigma_tau * sigtau) / dtau))
            ibin_max = min(ntau - 1, int((tau_delay + nsigma_tau * sigtau) / dtau))

            # Get Cloudy line emissivity at this location
            log_phi_loc = log_phi_now[iphi]
            if no_fluxfloor and log_phi_loc < phi_grid_min - 1.0:
                continue  # zero emissivity well below Cloudy grid
            log_phi_loc = np.clip(log_phi_loc, phi_grid_min, phi_grid_max)
            # Cloudy interp evaluation. Use grid=False and wrap in
            # try/except so that any scipy/dfitpack errors degrade
            # gracefully to "skip this cell" rather than propagating an
            # unpicklable Fortran exception across the multiprocessing
            # boundary.
            try:
                em_val = cloudy_interp(log_phi_loc, Z_target, grid=False)
                line_emissivity = float(np.asarray(em_val).ravel()[0])
            except Exception:
                continue

            # Ensure positive, finite emissivity
            if not np.isfinite(line_emissivity) or line_emissivity <= 0:
                continue

            if weighting == 'responsivity':
                # Finite-difference responsivity:
                # eta = d log F_line / d log Phi_H, computed via a
                # symmetric difference on the spline around log_phi_loc.
                # Cells touching the grid edge are clamped; cells where
                # F_line crosses zero get eta=0.
                dphi_fd = 0.05
                lp_lo = max(phi_grid_min, log_phi_loc - dphi_fd)
                lp_hi = min(phi_grid_max, log_phi_loc + dphi_fd)
                try:
                    f_lo_val = cloudy_interp(lp_lo, Z_target, grid=False)
                    f_hi_val = cloudy_interp(lp_hi, Z_target, grid=False)
                    f_lo = float(np.asarray(f_lo_val).ravel()[0])
                    f_hi = float(np.asarray(f_hi_val).ravel()[0])
                except Exception:
                    f_lo, f_hi = 0.0, 0.0
                if (np.isfinite(f_lo) and np.isfinite(f_hi)
                        and f_lo > 0 and f_hi > 0
                        and (lp_hi - lp_lo) > 1e-9):
                    eta = (np.log10(f_hi) - np.log10(f_lo)) / (lp_hi - lp_lo)
                else:
                    eta = 0.0
                # Linear response weight: F_line * eta is the cell's
                # contribution to dF_line/dlog(Phi_H). Drop the global
                # 1/L_ion factor (constant across cells).
                line_emissivity = line_emissivity * eta
                if line_emissivity == 0.0 or not np.isfinite(line_emissivity):
                    continue

            # Add contribution with 2D Gaussian spread
            if gauss_lambda_norm > 0:
                for ibin in range(ibin_min, ibin_max + 1):
                    tau_rel = tau_delay - tau_grid[ibin]
                    gauss_tau = np.exp(-0.5 * (tau_rel / sigtau)**2)

                    for ilambda in range(ilambda_min, ilambda_max + 1):
                        lambda_rel = lambda_obs - lambda_grid[ilambda]
                        gauss_lambda = np.exp(-0.5 * (lambda_rel / siglambda)**2) / gauss_lambda_norm

                        # Weight by area, foreshortening, Cloudy emissivity, and Gaussians
                        weight = da[iphi] * dot * line_emissivity * gauss_tau * gauss_lambda
                        psi_map[ibin, ilambda] += weight

    # Pillar face emission is no longer treated separately — the ionizing flux
    # computation now uses the full surface normal (including azimuthal gradient
    # from pillar walls), so pillar geometry is handled uniformly with the disk.

    # Do NOT normalize the map - we want to see absolute variations in response
    # as the pillar moves through different orbital phases. Normalization would
    # hide the foreshortening effect.
    # (The disk contribution provides a stable reference; pillar adds on top.)

    return lambda_grid, tau_grid, psi_map


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
    from pillardisk.pillar_disk import PillarDisk, resolve_rin
    from pillardisk.pillar_line_cloudy import compute_ionizing_flux
    
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
    
    # Static lab-frame grid. Cells stay put at the fixed (r, phi) grid; only
    # the pillars rotate (next block) at their own omega(r_p). The previous
    # version rotated each cell by omega(r) so that finite-difference dh/dr
    # at adjacent ir / same iphi spanned DIFFERENT lab azimuths (because
    # omega(r_ir) != omega(r_ir+1)), conflating the radial gradient with an
    # azimuthal one whose magnitude grew linearly in t. That contaminated
    # the surface-tilt and produced a secular drift in the lightcurve.
    # The lab-frame view is physically equivalent (sum over the full ring is
    # the same) and gives a clean radial gradient at fixed lab phi.
    phi_t_2d = phi_0_2d.copy()

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

    # TODO DEBUG: The foreshortening and observer-occultation effects below
    # need careful validation.  The direct spatial sum F(lambda,t) should
    # agree with the VDM collapsed over delay (int Psi(lambda,tau) dtau)
    # for a steady continuum, but currently the two give different results.
    # Possible issues:
    #   1. Observer occultation geometry (dr_occ, face_foreshort) may be
    #      too aggressive or applied to wrong cells.
    #   2. Pillar face emission (lines 760+) double-counts or mis-weights
    #      the illuminated vs shadow face contribution.
    #   3. The foreshortening dot product (line 705) should match the VDM's
    #      cos(theta) weighting but may differ in sign convention.
    # To debug: compare F(lambda) from this code vs sum_tau Psi(lambda,tau)
    # at t=0 with identical pillar configuration.

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

        if itime % 10 == 0 and os.environ.get('PILLAR_DEBUG'):
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
        no_fluxfloor = getattr(disk, 'no_fluxfloor', False)
        floor_cutoff = float(phi_grid[0]) - 1.0
        for ir in range(disk.nr):
            for iphi in range(disk.nphi):
                log_phi_val = log_phi_2d[ir, iphi]
                if no_fluxfloor and log_phi_val < floor_cutoff:
                    intensity_2d[ir, iphi] = 0.0
                    continue
                log_phi_val = np.clip(log_phi_val, phi_grid[0], phi_grid[-1])
                intensity_2d[ir, iphi] = float(cloudy_interp_dict[cloudy_key](log_phi_val, Z_target, grid=False))
        
        line_intensity_2d[line_key] = intensity_2d
    
    # Initialize flux slices for each line
    flux_slice_dict = {}
    for line_key in lambda0_dict.keys():
        flux_slice_dict[line_key] = np.zeros(nlambda)

    # Full 3D outward surface normal (uses both dh/dr AND dh/dphi). Same
    # function as the lamp-illumination path so the observer foreshortening
    # is computed self-consistently with the lamp's incidence angle.
    nx_2d, ny_2d, nz_2d = disk._surface_normal(r_2d, phi_t_2d, h_t_2d)
    dot_2d = disk.ex * nx_2d + disk.ey * ny_2d + disk.ez * nz_2d
    # DEBUG: amplify the *deviation* of the dot product from the face-on
    # baseline cos(i). A uniform multiplier on dot_2d cancels in any
    # fractional-flux quantity; only amplifying the deviation increases the
    # modulation depth. gain = 1 reproduces the physical foreshortening.
    foreshortening_gain = 1.0
    dot_baseline = disk.cosi   # dot for a flat surface (face-on contribution)
    dot_2d = dot_baseline + foreshortening_gain * (dot_2d - dot_baseline)

    # Loop over disk surface
    for ir in range(disk.nr - 1):
        r_now = disk.r[ir]
        h_now = h_t_2d[ir, :]

        # Next radius
        r_next = disk.r[ir + 1]
        h_next = h_t_2d[ir + 1, :]

        # Radial surface element (used for area only; the foreshortening dot
        # product uses the full 3D normal computed above).
        dr = r_next - r_now
        dh = h_next - h_now
        ds = np.sqrt(dr**2 + dh**2)
        da = ds * r_now * disk.dphi

        # Loop over azimuth
        for iphi, phi in enumerate(phi_t_2d[ir, :]):
            # Foreshortening factor from the full 3D normal
            dot = dot_2d[ir, iphi]

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
                
                # Wavelength spread: turbulent broadening + numerical floor
                # sigma_turb = f_turb * v_Kep(r) * lambda0 / c
                f_turb = disk_params.get('f_turb_line', 0.0)
                if f_turb > 0:
                    siglambda = max(lambda0 * f_turb * v_phi / C, dlambda)
                else:
                    siglambda = dlambda
                nsigma = 3 if f_turb > 0 else 2
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
    
    # Pillar face emission is no longer treated separately — the ionizing flux
    # computation now uses the full surface normal (including azimuthal gradient
    # from pillar walls), so pillar geometry is handled uniformly with the disk.

    # Restore original pillar positions (rotated earlier for this time step)
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
        'dmpc': disk.dmpc, 'cosi': disk.cosi, 'fcol': disk.fcol, 'redshift': disk.redshift,
        'flux_method': disk.flux_method, 'log_phi_inner': disk.log_phi_inner,
        'no_fluxfloor': disk.no_fluxfloor,
        'transparent': disk.transparent, 'f_trans': disk.f_trans,
        'f_turb_line': disk.f_turb_line,
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
    line_labels = {'Halpha': r'$\rm H\alpha$', 'Mg2': r'$\rm MgII$', 'C4': r'$\rm CIV$'}

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
        ax.set_title(line_label, fontsize=16, pad=10)
        ax.set_xlim(lambda_grid[0], lambda_grid[-1])
        ax.set_ylim(time_grid[0], time_grid[-1])
        plt.colorbar(im, ax=ax, label=cbar_label)
        ax.minorticks_on()
        ax.tick_params(top=True, right=True, axis='both', which='major', length=8, width=1.5, direction='in')
        ax.tick_params(top=True, right=True, axis='both', which='minor', length=4, width=1, direction='in')
        
        # Add velocity axis on top using secondary_xaxis for auto tick placement
        lam2v = lambda lam: (C / KM_TO_CM) * (lam - lambda0) / lambda0
        v2lam = lambda v: lambda0 * (1.0 + v * KM_TO_CM / C)
        ax2 = ax.secondary_xaxis('top', functions=(lam2v, v2lam))
        ax2.set_xlabel(r'$v~\rm [km\,s^{-1}]$', labelpad=10)
        ax2.xaxis.set_major_locator(plt.MaxNLocator(nbins=5, integer=True, steps=[1, 2, 3, 5, 6, 10]))
        ax2.minorticks_on()
        ax2.tick_params(axis='both', which='major', length=8, width=1.5, direction='in')
        ax2.tick_params(axis='both', which='minor', length=4, width=1, direction='in')
    
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

    line_labels = {'Halpha': r'$\rm H\alpha$', 'Mg2': r'$\rm MgII$', 'C4': r'$\rm CIV$'}

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
        ax.set_title(line_label + r'$\rm ~(with~pillars - without~pillars)$', fontsize=16, pad=10)
        ax.set_xlim(lambda_grid[0], lambda_grid[-1])
        ax.set_ylim(time_grid[0], time_grid[-1])
        plt.colorbar(im, ax=ax, label=r'$\rm \Delta Flux$')
        ax.minorticks_on()
        ax.tick_params(top=True, right=True, axis='both', which='major', length=8, width=1.5, direction='in')
        ax.tick_params(top=True, right=True, axis='both', which='minor', length=4, width=1, direction='in')

        # Add velocity axis on top using secondary_xaxis for auto tick placement
        lam2v = lambda lam: (C / KM_TO_CM) * (lam - lambda0) / lambda0
        v2lam = lambda v: lambda0 * (1.0 + v * KM_TO_CM / C)
        ax2 = ax.secondary_xaxis('top', functions=(lam2v, v2lam))
        ax2.set_xlabel(r'$v~\rm [km\,s^{-1}]$', labelpad=10)
        ax2.xaxis.set_major_locator(plt.MaxNLocator(nbins=5, integer=True, steps=[1, 2, 3, 5, 6, 10]))
        ax2.minorticks_on()
        ax2.tick_params(axis='both', which='major', length=8, width=1.5, direction='in')
        ax2.tick_params(axis='both', which='minor', length=4, width=1, direction='in')

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

    line_labels = {'Halpha': r'$\rm H\alpha$', 'Mg2': r'$\rm MgII$', 'C4': r'$\rm CIV$'}

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
        ax.set_title(line_label + r'$\rm ~Residual~(F(t,\lambda) - \langle F(\lambda) \rangle_t)$', fontsize=16, pad=10)
        ax.set_xlim(lambda_grid[0], lambda_grid[-1])
        ax.set_ylim(time_grid[0], time_grid[-1])
        plt.colorbar(im, ax=ax, label=r'$\rm Residual~Flux$')
        ax.minorticks_on()
        ax.tick_params(top=True, right=True, axis='both', which='major', length=8, width=1.5, direction='in')
        ax.tick_params(top=True, right=True, axis='both', which='minor', length=4, width=1, direction='in')

        # Add velocity axis on top using secondary_xaxis for auto tick placement
        lam2v = lambda lam: (C / KM_TO_CM) * (lam - lambda0) / lambda0
        v2lam = lambda v: lambda0 * (1.0 + v * KM_TO_CM / C)
        ax2 = ax.secondary_xaxis('top', functions=(lam2v, v2lam))
        ax2.set_xlabel(r'$v~\rm [km\,s^{-1}]$', labelpad=10)
        ax2.xaxis.set_major_locator(plt.MaxNLocator(nbins=5, integer=True, steps=[1, 2, 3, 5, 6, 10]))
        ax2.minorticks_on()
        ax2.tick_params(axis='both', which='major', length=8, width=1.5, direction='in')
        ax2.tick_params(axis='both', which='minor', length=4, width=1, direction='in')

    plt.tight_layout()
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    print(f"Residual (time-averaged subtracted) maps saved to {filename}")
    plt.close()


def plot_time_evolving_raw_residual_combined(lambda_grid_dict, time_grid, flux_map_dict,
                                              lambda0_dict,
                                              filename='velocity_time_cloudy_combined.png'):
    """
    Combined 3x2 figure: rows = emission lines (CIV, MgII, Halpha by increasing
    rest wavelength), columns = (left: raw F(lambda, t); right: residual
    F(lambda, t) - <F(lambda)>_t).

    Parameters
    ----------
    lambda_grid_dict : dict[str, ndarray]
        Wavelength grid per line.
    time_grid : ndarray
        Time grid (days).
    flux_map_dict : dict[str, ndarray]
        Flux map per line, shape (ntime, nlambda).
    lambda0_dict : dict[str, float]
        Rest wavelength per line.
    filename : str
        Output filename.
    """
    line_order = [k for k in ['C4', 'Mg2', 'Halpha'] if k in flux_map_dict]
    n_rows = len(line_order)
    if n_rows == 0:
        return

    print(f"  Creating combined raw + residual figure with {n_rows*2} panels...")
    fig, axes = plt.subplots(n_rows, 2, figsize=(15, 3.6*n_rows + 0.8))
    if n_rows == 1:
        axes = axes[np.newaxis, :]

    line_labels = {'Halpha': r'$\rm H\alpha$',
                   'Mg2':    r'$\rm Mg\,II$',
                   'C4':     r'$\rm C\,IV$'}

    from matplotlib.ticker import FixedLocator, NullFormatter
    # Round velocity ticks (km/s) on the top secondary axis.
    v_ticks_major = np.array([-9000, -6000, -3000, 0, 3000, 6000, 9000])
    v_ticks_minor = np.array([-7500, -4500, -1500, 1500, 4500, 7500])

    for i, line_key in enumerate(line_order):
        lambda_grid = lambda_grid_dict[line_key]
        flux_map = flux_map_dict[line_key]
        lambda0 = lambda0_dict[line_key]
        line_label = line_labels.get(line_key, line_key)

        time_avg = np.mean(flux_map, axis=0)
        residual = flux_map - time_avg[np.newaxis, :]

        # Per-line normalisation by the global mean <F> over both lambda and
        # t (a single scalar per line). Same denominator on the residual.
        norm = float(np.mean(flux_map))
        if not np.isfinite(norm) or norm == 0:
            norm = 1.0
        flux_norm = flux_map / norm
        residual_norm = residual / norm

        ax_raw = axes[i, 0]
        vmax_raw = float(np.nanmax(flux_norm))
        if not np.isfinite(vmax_raw) or vmax_raw == 0:
            vmax_raw = 1.0
        im_raw = ax_raw.pcolormesh(lambda_grid, time_grid, flux_norm,
                                    shading='auto', cmap='rainbow',
                                    vmin=0.0, vmax=vmax_raw)
        ax_raw.axvline(lambda0, color='white', linestyle='--',
                       linewidth=1.5, alpha=0.7)
        cb_raw = plt.colorbar(im_raw, ax=ax_raw, shrink=0.85, aspect=20,
                              extend='neither')
        cb_raw.set_label(r'$F\,/\,\langle F\rangle$', fontsize=15)
        cb_raw.ax.tick_params(labelsize=12)
        # Line name as in-panel text instead of title.
        ax_raw.text(0.03, 0.93, line_label,
                    transform=ax_raw.transAxes, fontsize=18,
                    color='white', ha='left', va='top',
                    bbox=dict(boxstyle='round,pad=0.25',
                              facecolor='black', alpha=0.45,
                              edgecolor='none'))

        ax_res = axes[i, 1]
        vmax_res = float(np.nanmax(np.abs(residual_norm)))
        if not np.isfinite(vmax_res) or vmax_res == 0:
            vmax_res = 1.0
        im_res = ax_res.pcolormesh(lambda_grid, time_grid, residual_norm,
                                    shading='auto', cmap='seismic',
                                    vmin=-vmax_res, vmax=vmax_res)
        ax_res.axvline(lambda0, color='black', linestyle='--',
                       linewidth=1.5, alpha=0.7)
        cb_res = plt.colorbar(im_res, ax=ax_res, shrink=0.85, aspect=20,
                              extend='neither')
        cb_res.set_label(r'$(F - \langle F\rangle_t)\,/\,\langle F\rangle$',
                         fontsize=15)
        cb_res.ax.tick_params(labelsize=12)
        ax_res.text(0.03, 0.93, line_label,
                    transform=ax_res.transAxes, fontsize=18,
                    color='black', ha='left', va='top',
                    bbox=dict(boxstyle='round,pad=0.25',
                              facecolor='white', alpha=0.6,
                              edgecolor='none'))

        for ax in (ax_raw, ax_res):
            ax.set_xlabel(r'$\lambda~[\rm \AA]$', fontsize=18)
            ax.set_ylabel(r'$\rm Time~[days]$', fontsize=18)
            ax.set_xlim(lambda_grid[0], lambda_grid[-1])
            ax.set_ylim(time_grid[0], time_grid[-1])
            ax.minorticks_on()
            ax.tick_params(top=False, right=True, axis='both', which='major',
                           length=8, width=1.5, direction='in', labelsize=14)
            ax.tick_params(top=False, right=True, axis='both', which='minor',
                           length=5, width=1, direction='in')

            # Top secondary axis: velocity (km/s) per row.
            lam2v = lambda lam, l0=lambda0: (C / KM_TO_CM) * (lam - l0) / l0
            v2lam = lambda v, l0=lambda0: l0 * (1.0 + v * KM_TO_CM / C)
            ax2 = ax.secondary_xaxis('top', functions=(lam2v, v2lam))
            ax2.set_xlabel(r'$v~[\rm km\,s^{-1}]$', fontsize=18, labelpad=8)
            v_lo = lam2v(lambda_grid[0])
            v_hi = lam2v(lambda_grid[-1])
            major_in_range = v_ticks_major[(v_ticks_major >= v_lo)
                                           & (v_ticks_major <= v_hi)]
            minor_in_range = v_ticks_minor[(v_ticks_minor >= v_lo)
                                           & (v_ticks_minor <= v_hi)]
            ax2.xaxis.set_major_locator(FixedLocator(major_in_range))
            ax2.xaxis.set_minor_locator(FixedLocator(minor_in_range))
            ax2.xaxis.set_minor_formatter(NullFormatter())
            ax2.tick_params(axis='both', which='major',
                            length=8, width=1.5, direction='in', labelsize=14)
            ax2.tick_params(axis='both', which='minor',
                            length=5, width=1, direction='in')

    plt.tight_layout()
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    print(f"Combined raw + residual maps saved to {filename}")
    plt.close()


def plot_time_evolving_lightcurves(lambda_grid_dict, time_grid, flux_map_dict,
                                    filename='velocity_time_cloudy_lightcurves.png'):
    """
    Wavelength-integrated lightcurves L(t) = int F(lambda, t) dlambda for each
    line, plotted as fractional variation (L(t) - <L>_t) / <L>_t so the lines
    are directly comparable. Reveals the line-dependent depth of the
    shadow-induced dip: CIV should drop the most because its emissivity
    responds steepest to log Phi.

    Parameters
    ----------
    lambda_grid_dict : dict[str, ndarray]
        Wavelength grid per line.
    time_grid : ndarray
        Time grid (days).
    flux_map_dict : dict[str, ndarray]
        Flux map per line, shape (ntime, nlambda).
    filename : str
        Output filename.
    """
    line_order = [k for k in ['C4', 'Mg2', 'Halpha'] if k in flux_map_dict]
    if not line_order:
        return

    line_labels = {'Halpha': r'$\rm H\alpha$',
                   'Mg2':    r'$\rm Mg\,II$',
                   'C4':     r'$\rm C\,IV$'}
    line_colors = {'Halpha': 'orangered',
                   'Mg2':    'forestgreen',
                   'C4':     'royalblue'}

    fig, ax = plt.subplots(1, 1, figsize=(11, 6.5))

    for line_key in line_order:
        lam = lambda_grid_dict[line_key]
        flux_map = flux_map_dict[line_key]
        # Wavelength-integrated lightcurve.
        L_t = np.trapezoid(flux_map, lam, axis=1)
        L_mean = float(np.mean(L_t))
        if not np.isfinite(L_mean) or L_mean == 0:
            continue
        frac = (L_t - L_mean) / L_mean
        ax.plot(time_grid, frac,
                color=line_colors.get(line_key, 'k'),
                linewidth=3, alpha=0.85,
                label=line_labels.get(line_key, line_key))

    ax.axhline(0.0, color='k', linewidth=1, alpha=0.4)
    ax.set_xlabel(r'$\rm Time~[days]$', fontsize=18)
    ax.set_ylabel(r'$(L(t) - \langle L\rangle_t)\,/\,\langle L\rangle_t$',
                  fontsize=18)
    ax.set_xlim(time_grid[0], time_grid[-1])
    ax.minorticks_on()
    ax.tick_params(top=True, right=True, axis='both', which='major',
                   length=8, width=1.5, direction='in', labelsize=14)
    ax.tick_params(top=True, right=True, axis='both', which='minor',
                   length=5, width=1, direction='in')
    ax.legend(fontsize=15, frameon=False, loc='best')

    plt.tight_layout()
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    print(f"Wavelength-integrated lightcurves saved to {filename}")
    plt.close()


def main(config_file='config_line.yaml', use_absolute_flux=True, plot_only=False):
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
                   'tx1', 'tirrad_tvisc_ratio', 'fcol', 'hlamp', 'dmpc', 'cosi', 'inclination', 'redshift',
                   'M_BH', 'r_isco_rg', 'f_trans', 'f_turb_line']
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

    # 3D geometry plots at different orbital phases for debugging foreshortening
    # The observer is FIXED (looking from +x direction at inclination i)
    # We show the pillar at different azimuthal positions (simulating rotation)
    # - near: φ=0 (pillar between observer and lamp) → observer sees shadow side
    # - quad1: φ=π/2 (pillar at +y) → observer sees pillar from side
    # - far: φ=π (pillar behind lamp from observer) → observer sees illuminated side
    # - quad2: φ=3π/2 (pillar at -y) → observer sees pillar from other side

    # Save original pillar positions
    original_pillars = [p.copy() for p in disk.pillars]

    # Get v_virial from config for orbital period calculation
    line_params = config.get('emission_line', {})
    if 'M_BH' in line_params:
        v_virial_geom = v_virial_from_mbh(float(line_params['M_BH']), disk.r0)
    else:
        v_virial_geom = float(line_params.get('v_virial', 1000.0))

    # Compute orbital period at pillar radius for time labels
    r_pillar = original_pillars[0]['r'] if len(original_pillars) > 0 else disk.r0
    T_orbital_pillar = compute_orbital_period(r_pillar, v_virial_geom, disk.r0)

    orbital_phases = [
        (0.0, 'phi0_near'),           # φ=0: near side (shadow visible)
        (np.pi/2, 'phi90_quad1'),     # φ=π/2: quadrature 1
        (np.pi, 'phi180_far'),        # φ=π: far side (illuminated visible)
        (3*np.pi/2, 'phi270_quad2'),  # φ=3π/2: quadrature 2
    ]

    for phi_offset, label in orbital_phases:
        # Temporarily move all pillars to the new azimuthal position
        for i, p in enumerate(disk.pillars):
            p['phi'] = phi_offset + (original_pillars[i]['phi'] - original_pillars[0]['phi'])

        # Compute time corresponding to this orbital phase
        t_days = phi_offset / (2 * np.pi) * T_orbital_pillar
        title = f't = {t_days:.0f} days'

        filename_3d = f'plots/geometry_{label}.png'
        os.makedirs('plots', exist_ok=True)
        # Lower resolution for faster plotting
        disk.plot_3d_geometry(show_light_rays=False, color_by='flux', filename=filename_3d,
                              title=title, n_phi_plot=180, n_r_plot=100)

    # Restore original pillar positions
    for i, p in enumerate(disk.pillars):
        p['phi'] = original_pillars[i]['phi']

    print(f"  Saved 3D geometry snapshots at different orbital phases:")
    print(f"    phi0_near: pillar at φ=0, observer sees shadow side")
    print(f"    phi180_far: pillar at φ=π, observer sees illuminated side")
    print(f"    phi90_quad1, phi270_quad2: pillar at quadratures")

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

    fig, axes = plt.subplots(2, 1, figsize=(7, 12), subplot_kw={'projection': 'polar'})

    for ax in axes:
        ax.grid(True, alpha=0.3)
        ax.tick_params(pad=5)
        angle_labels = [rf'$\bf{{{int(a)}}}$°' for a in np.degrees(ax.get_xticks())]
        ax.set_xticklabels(angle_labels)

    ax0 = axes[0]
    c0 = ax0.pcolormesh(phi_2d_plot, r_2d_plot, h_2d_plot, shading='auto', cmap='terrain')
    ax0.set_title(r'$\rm Disk~Height$', pad=15)
    cbar0 = fig.colorbar(c0, ax=ax0, pad=0.1)
    cbar0.set_label(r'$H~\rm [ld]$')

    ax1 = axes[1]
    c1 = ax1.pcolormesh(phi_2d_plot, r_2d_plot, log_phi_2d_plot, shading='auto', cmap='inferno',
                        vmin=15.0, vmax=21.5)
    ax1.set_title(r'$\rm Ionizing~Flux$', pad=15)
    cbar1 = fig.colorbar(c1, ax=ax1, pad=0.1)
    cbar1.set_label(r'$\log\Phi~\rm [photons~cm^{-2}~s^{-1}]$')

    plt.tight_layout()
    _geom_out = 'plots/geometry_t0_faceon.png'
    os.makedirs(os.path.dirname(_geom_out), exist_ok=True)
    plt.savefig(_geom_out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved face-on map: {_geom_out}")
    print("--- End geometry check ---\n")

    if plot_only:
        print("--plot-only flag set, skipping spectral computation.")
        return

    # Load Cloudy models
    cloudy_extension_file = None
    if use_absolute_flux:
        cloudy_file = '/Users/jiamuh/c23.01/my_models/loc_metal_flux/strong_LOC_varym_N25_v100_lineflux_LineList_BLR_Fe2_flux.txt'
        cloudy_extension_file = '/Users/jiamuh/c23.01/my_models/loc_metal_flux/strong_LOC_varym_N25_v100_lineflux_extlow_LineList_BLR_Fe2_flux.txt'
    else:
        cloudy_file = '/Users/jiamuh/c23.01/my_models/loc_metal/strong_LOC_varym_N25_LineList_BLR_Fe2.txt'

    print(f"\nLoading Cloudy models from: {cloudy_file}")
    cloudy_interp_dict, phi_grid, Z_grid = load_cloudy_models(
        cloudy_file, Z_target=1.0, extension_file=cloudy_extension_file)
    
    # Emission line parameters
    lambda0_dict = {
        'Halpha': 6562.8,   # Halpha rest wavelength
        'Mg2': 2798.0,      # MgII rest wavelength
        'C4': 1549.0        # CIV rest wavelength
    }
    
    v_virial = 1000.0  # Virial velocity (km/s) - default fallback

    line_params = config.get('emission_line', {})
    if line_params:
        # M_BH takes priority: compute v_virial from SMBH mass
        if 'M_BH' in line_params:
            M_BH = float(line_params['M_BH'])
            v_virial = v_virial_from_mbh(M_BH, disk.r0)
            print(f"  M_BH = {M_BH:.2e} Msun → v_virial = {v_virial:.1f} km/s at r0 = {disk.r0:.1f} ld")
        else:
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

    # Round tmax to integer orbital periods for exact periodicity, unless the
    # user has set ``round_tmax: false`` (e.g. when an orchestrator wants the
    # SAME tmax across multiple runs that have different mean pillar radii).
    round_tmax_flag = comp_params.get('round_tmax', True)
    if isinstance(round_tmax_flag, str):
        round_tmax_flag = round_tmax_flag.lower() not in ('false', '0', 'no')
    round_tmax_flag = bool(round_tmax_flag)

    if len(disk.pillars) > 0:
        pillar_radii = [p['r'] for p in disk.pillars]
        r_pillar_mean = np.mean(pillar_radii)
        T_orbital = compute_orbital_period(r_pillar_mean, v_virial, disk.r0)
        if tmax is not None:
            if round_tmax_flag:
                n_orbits = max(1, round(tmax / T_orbital))
                tmax_orig = tmax
                tmax = n_orbits * T_orbital
                print(f"Adjusted tmax: {tmax_orig:.1f} → {tmax:.1f} days "
                      f"({n_orbits} × T_orbital={T_orbital:.1f}d at r={r_pillar_mean:.1f})")
            else:
                print(f"Using tmax = {tmax:.1f} days verbatim "
                      f"(round_tmax disabled; T_orbital={T_orbital:.1f}d at r={r_pillar_mean:.1f})")
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
    # Plot main maps (with pillars). All default filenames are routed
    # through plots/ so intermediate files don't pollute the project root.
    plot_params = config.get('plotting', {})
    os.makedirs('plots', exist_ok=True)
    filename = plot_params.get('velocity_time_cloudy_filename',
                               'plots/velocity_time_cloudy.png')
    plot_time_evolving_cloudy_map(
        lambda_grid_dict, time_grid, flux_with_pillars, lambda0_dict,
        filename=filename, use_absolute_flux=use_absolute_flux,
        normalize_output=False
    )

    # Plot difference maps (with pillars - without pillars)
    diff_filename = plot_params.get('velocity_time_cloudy_diff_filename',
                                    'plots/velocity_time_cloudy_diff.png')
    plot_time_evolving_diff_map(
        lambda_grid_dict, time_grid,
        flux_with_pillars, flux_without_pillars,
        lambda0_dict,
        filename=diff_filename
    )

    # Plot true residual maps (flux - time-averaged flux)
    residual_filename = plot_params.get('velocity_time_cloudy_residual_filename',
                                        'plots/velocity_time_cloudy_residual.png')
    plot_time_evolving_residual_map(
        lambda_grid_dict, time_grid,
        flux_with_pillars, lambda0_dict,
        filename=residual_filename
    )

    # Combined 3x2 figure: rows = lines (CIV/MgII/Halpha), cols = (raw, residual)
    combined_filename = plot_params.get('velocity_time_cloudy_combined_filename',
                                         'plots/velocity_time_cloudy_combined.png')
    plot_time_evolving_raw_residual_combined(
        lambda_grid_dict, time_grid,
        flux_with_pillars, lambda0_dict,
        filename=combined_filename
    )

    # Wavelength-integrated lightcurves (velocity-integrated time series).
    lightcurve_filename = plot_params.get('velocity_time_cloudy_lightcurves_filename',
                                           'plots/velocity_time_cloudy_lightcurves.png')
    plot_time_evolving_lightcurves(
        lambda_grid_dict, time_grid,
        flux_with_pillars,
        filename=lightcurve_filename
    )
    print("Plotting complete.")
    
    print("\nSummary:")
    for line_key, lambda_grid in lambda_grid_dict.items():
        print(f"{line_key}:")
        print(f"  Wavelength range: {lambda_grid[0]:.2f} - {lambda_grid[-1]:.2f} Å")
        print(f"  Max flux (with pillars): {np.max(flux_with_pillars[line_key]):.2e}")
        print(f"  Max flux (without pillars): {np.max(flux_without_pillars[line_key]):.2e}")
    print(f"Time range: {time_grid[0]:.2f} - {time_grid[-1]:.2f} days")


def _load_cloudy_raw_data(cloudy_file, Z_target=1.0, extension_file=None,
                           ext_phi_grid=None):
    """
    Load Cloudy models and return raw collapsed data (picklable numpy arrays).

    Optionally also load an ``extension_file`` covering a lower
    ``log Phi`` range (e.g. 15 to 16.5) and concatenate it on the phi
    axis. ``ext_phi_grid`` defaults to ``arange(15, 17, 0.5)``.

    This avoids reloading the file in each worker process.
    """
    import pandas as pd

    def _load_one(path, phi_grid_arr):
        with open(path) as fh:
            header_line = fh.readline().lstrip('#').strip()
        col_names = header_line.split('\t')
        df = pd.read_csv(path, sep='\t', names=col_names, comment='#',
                         skiprows=0)
        df = df[df.iloc[:, 0] == 'iteration 1'].reset_index(drop=True)
        n_grid_loc = np.arange(9, 12.5, 0.5)
        n_Z_loc = len(df) // (len(phi_grid_arr) * len(n_grid_loc))
        Z_grid_loc = np.arange(1, 1 + n_Z_loc * 0.5, 0.5)
        expected = len(Z_grid_loc) * len(phi_grid_arr) * len(n_grid_loc)
        assert len(df) == expected, (
            f"Row count mismatch for {path}: got {len(df)}, expected "
            f"{expected}")
        shape = (len(n_grid_loc), len(phi_grid_arr), len(Z_grid_loc))
        intensity = {key: np.zeros(shape) for key in line_names}
        for j_phi, _ in enumerate(phi_grid_arr):
            for i_n, _ in enumerate(n_grid_loc):
                index = j_phi * len(n_grid_loc) + i_n
                start = index * len(Z_grid_loc)
                end = start + len(Z_grid_loc)
                df_slice = df.iloc[start:end]
                for key, colname in line_names.items():
                    intensity[key][i_n, j_phi, :] = df_slice[colname].values
        return intensity, Z_grid_loc

    print(f"Loading Cloudy models from {cloudy_file}...")

    # Grid definitions — phi and n are fixed; Z auto-detected from data
    phi_grid = np.arange(17, 21.5, 0.5)
    n_grid = np.arange(9, 12.5, 0.5)

    # Line names mapping
    line_names = {
        'Halpha': 'h  1 6562.80A',
        'HI':  'H  1 4861.32A',
        'Mg2': 'blnd 2798.00A',
        'C4':  'blnd 1549.00A',
    }

    intensity_data, Z_grid = _load_one(cloudy_file, phi_grid)

    # Optionally splice in the low-phi extension grid. The extension
    # file has only one Z column (Z = Z_sun); replicate it across all
    # Z values of the main grid so the resulting spline has the same
    # Z axis as the main file (otherwise it collapses to len(Z)=1 and
    # ky=2 spline construction fails with a dfitpack error).
    if extension_file is not None:
        if ext_phi_grid is None:
            ext_phi_grid = np.arange(15, 17, 0.5)
        ext_phi_grid = np.asarray(ext_phi_grid, dtype=float)
        print(f"  Splicing extension file {extension_file} "
              f"(phi {ext_phi_grid[0]} to {ext_phi_grid[-1]})...")
        intensity_ext, Z_grid_ext = _load_one(extension_file, ext_phi_grid)
        # Replicate the Z_sun column of the extension across all Z
        # values of the main grid (the paper only ever queries Z=Z_sun,
        # so this is fine in practice).
        for key in line_names:
            # intensity_ext[key] has shape (n_grid, len(ext_phi_grid), n_Z_ext);
            # tile along the Z axis to match Z_grid.
            ext_arr = intensity_ext[key]
            if ext_arr.shape[2] != len(Z_grid):
                ext_arr = np.tile(ext_arr[:, :, :1], (1, 1, len(Z_grid)))
            intensity_ext[key] = ext_arr
        # Concatenate along the phi axis.
        phi_grid = np.concatenate([ext_phi_grid, phi_grid])
        for key in line_names:
            intensity_data[key] = np.concatenate(
                [intensity_ext[key], intensity_data[key]], axis=1)

    # Collapse over density (average along axis=0)
    collapsed_data = {}
    for key in line_names:
        collapsed_data[key] = np.mean(intensity_data[key], axis=0)

    print(f"  Loaded lines: {list(line_names.keys())}")
    print(f"  phi range: {phi_grid[0]:.1f} to {phi_grid[-1]:.1f}")

    if os.environ.get('PILLAR_DEBUG'):
        Z_idx = np.argmin(np.abs(Z_grid - Z_target))
        print(f"\n  [DEBUG] Emissivity vs log_phi at Z={Z_grid[Z_idx]:.1f}:")
        print(f"  log_phi  |     CIV      |    MgII     |   Halpha")
        print(f"  ---------|--------------|-------------|------------")
        for i, phi_val in enumerate(phi_grid):
            c4_val = collapsed_data['C4'][i, Z_idx]
            mg2_val = collapsed_data['Mg2'][i, Z_idx]
            ha_val = collapsed_data['Halpha'][i, Z_idx]
            print(f"    {phi_val:.1f}   |  {c4_val:10.4e}  | {mg2_val:10.4e} | {ha_val:10.4e}")
        print()

    return collapsed_data, phi_grid, Z_grid


def _compute_vdmap_for_line_cloudy(args):
    """
    Compute Cloudy-weighted velocity-delay map for a single line (for parallel processing).

    This function is at module level so it can be pickled by multiprocessing.
    Receives pre-loaded Cloudy data (raw arrays) and recreates interpolators in worker.

    Tuple may have 12 elements (legacy, weighting='emissivity') or 13
    elements (with explicit ``weighting``) for backward compatibility.
    """
    if len(args) == 13:
        (line_key, lambda0, disk_params, pillars_list, v_virial, nlambda, ntau, taumax,
         cloudy_data_raw, phi_grid, Z_grid, Z_target, weighting) = args
    else:
        (line_key, lambda0, disk_params, pillars_list, v_virial, nlambda, ntau, taumax,
         cloudy_data_raw, phi_grid, Z_grid, Z_target) = args
        weighting = 'emissivity'

    # Recreate disk in worker process
    disk_worker = PillarDisk(**disk_params)
    for p in pillars_list:
        disk_worker.add_pillar(r_pillar=p['r'], phi_pillar=p['phi'],
                               height=p['height'], sigma_r=p['sigma_r'],
                               sigma_phi=p['sigma_phi'])

    # Recreate interpolator for this line from raw data (fast, no file I/O).
    # Wrap in try/except: scipy's dfitpack.error is not picklable across
    # the multiprocessing boundary, so if construction fails we re-raise
    # as a RuntimeError that the parent process can decode.
    from scipy.interpolate import RectBivariateSpline
    try:
        cloudy_interp_dict = {
            line_key: RectBivariateSpline(phi_grid, Z_grid,
                                          cloudy_data_raw[line_key],
                                          kx=2, ky=2)
        }
    except Exception as e:
        raise RuntimeError(
            f"RectBivariateSpline construction failed for line "
            f"{line_key!r}: {type(e).__name__}: {e}; "
            f"phi_grid shape={np.asarray(phi_grid).shape}, "
            f"Z_grid shape={np.asarray(Z_grid).shape}, "
            f"data shape={cloudy_data_raw[line_key].shape}"
        )

    # Compute Cloudy-weighted velocity-delay map
    lambda_grid, tau_grid, psi_map = compute_velocity_delay_map_cloudy(
        disk_worker, line_key, lambda0, v_virial,
        cloudy_interp_dict, Z_target=Z_target,
        nlambda=nlambda, ntau=ntau, taumax=taumax,
        weighting=weighting,
    )
    return line_key, lambda_grid, tau_grid, psi_map


def generate_movie_frames(config_file='config_line.yaml', n_frames=50, output_dir='movie_frames',
                          use_absolute_flux=True, skip_geometry=False,
                          weighting='emissivity'):
    """
    Generate movie frames showing velocity-delay maps and 3D geometry at each orbital phase.

    Each frame shows a 2x3 or 3x3 grid:
    - Top row: Halpha, Mg2, C4 velocity-delay maps Ψ(λ, τ)
    - Middle row: Difference maps ΔΨ for each line
    - Bottom row (optional): 3D disk geometry

    Parameters
    ----------
    skip_geometry : bool
        If True, omit the 3D geometry panel (produces a 2-row figure).
    """
    import os
    os.makedirs(output_dir, exist_ok=True)

    # Load config and create disk (same as main)
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)

    # Extract parameters (simplified version of main)
    disk_params = config.get('disk', {}).copy()
    temp_params = config.get('temperature', {}).copy()
    lamp_params = config.get('lamp', {}).copy()
    obs_params = config.get('observation', {}).copy()

    float_params = ['rin', 'rout', 'h1', 'r0', 'beta', 'tv1', 'alpha',
                   'tx1', 'tirrad_tvisc_ratio', 'fcol', 'hlamp', 'dmpc', 'cosi', 'inclination', 'redshift',
                   'M_BH', 'r_isco_rg', 'f_trans', 'f_turb_line']
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
                pass
        if key in int_params:
            return int(float(value))
        elif key in float_params:
            return float(value)
        return value

    for params_dict in [disk_params, temp_params, lamp_params, obs_params]:
        for key, value in params_dict.items():
            params_dict[key] = convert_value(key, value)

    if 'inclination' in obs_params and obs_params['inclination'] is not None:
        inclination_rad = np.radians(obs_params['inclination'])
        obs_params['cosi'] = np.cos(inclination_rad)
        del obs_params['inclination']

    all_disk_params = {**disk_params, **temp_params, **lamp_params, **obs_params}

    # Resolve rin='auto' to ISCO if needed
    resolve_rin(all_disk_params)

    disk = PillarDisk(**all_disk_params)

    # Add pillars
    pillars_cfg = config.get('pillars', {})

    def to_list(x):
        if isinstance(x, list):
            return x
        return [x]

    def parse_math_expr_local(value):
        if isinstance(value, (int, float)):
            return float(value)
        if isinstance(value, str):
            expr = value.replace('pi', str(np.pi))
            try:
                return float(eval(expr))
            except:
                return float(value)
        return float(value)

    # Check for make_many mode
    make_many = pillars_cfg.get('make_many', False)
    if isinstance(make_many, str):
        make_many = make_many.lower() in ('true', '1', 'yes')
    make_many = bool(make_many)

    if make_many:
        # Random pillar generation
        print("Generating random pillars...")
        N_pillar = int(pillars_cfg.get('N_pillar', 100))
        r_mean = float(pillars_cfg.get('r_mean', 25.0))
        sig_r = float(pillars_cfg.get('sig_r', 10.0))
        rmin = pillars_cfg.get('rmin', None)
        rmin = float(rmin) if rmin is not None else None

        def get_pillar_param(key, default):
            val = pillars_cfg.get(key, default)
            if isinstance(val, list):
                return val
            return [val]

        def expand_list(lst, n, default):
            if len(lst) == 0:
                return [default] * n
            if len(lst) == 1:
                return lst * n
            while len(lst) < n:
                lst.append(lst[-1])
            return lst[:n]

        h_pillar_list = expand_list(get_pillar_param('h_pillar', 0.12), N_pillar, 0.12)
        sigma_r_pillar_list = expand_list(get_pillar_param('sigma_r_pillar', 2.0), N_pillar, 2.0)
        sigma_phi_pillar_list = expand_list(get_pillar_param('sigma_phi_pillar', 0.2), N_pillar, 0.2)

        r_pillar_random = np.random.normal(r_mean, sig_r, N_pillar)
        rmin_eff = max(disk.rin, rmin) if rmin is not None else disk.rin
        r_pillar_random = np.clip(r_pillar_random, rmin_eff, disk.rout)
        phi_pillar_random = np.random.uniform(0, 2*np.pi, N_pillar)

        for i in range(N_pillar):
            disk.add_pillar(
                r_pillar=float(r_pillar_random[i]),
                phi_pillar=float(phi_pillar_random[i]),
                height=float(h_pillar_list[i]),
                sigma_r=float(sigma_r_pillar_list[i]),
                sigma_phi=float(sigma_phi_pillar_list[i]),
                modify_height=True,
                modify_temp=True,
                temp_factor=1.5
            )
        print(f"  Generated {N_pillar} random pillars")
    else:
        # Manual pillar placement
        r_pillar_list = to_list(pillars_cfg.get('r_pillar', [10.0]))
        phi_pillar_list = to_list(pillars_cfg.get('phi_pillar', [0.0]))
        height_list = to_list(pillars_cfg.get('height', [0.1]))
        sigma_r_list = to_list(pillars_cfg.get('sigma_r', [0.1]))
        sigma_phi_list = to_list(pillars_cfg.get('sigma_phi', [0.1]))
        modify_height_list = to_list(pillars_cfg.get('modify_height', [True]))
        modify_temp_list = to_list(pillars_cfg.get('modify_temp', [True]))
        temp_factor_list = to_list(pillars_cfg.get('temp_factor', [1.0]))

        n_pillars_cfg = len(r_pillar_list)
        for i in range(n_pillars_cfg):
            pillar = {
                'r_pillar': float(r_pillar_list[i] if i < len(r_pillar_list) else r_pillar_list[-1]),
                'phi_pillar': parse_math_expr_local(phi_pillar_list[i] if i < len(phi_pillar_list) else phi_pillar_list[-1]),
                'height': float(height_list[i] if i < len(height_list) else height_list[-1]),
                'sigma_r': float(sigma_r_list[i] if i < len(sigma_r_list) else sigma_r_list[-1]),
                'sigma_phi': parse_math_expr_local(sigma_phi_list[i] if i < len(sigma_phi_list) else sigma_phi_list[-1]),
                'modify_height': bool(modify_height_list[i] if i < len(modify_height_list) else modify_height_list[-1]),
                'modify_temp': bool(modify_temp_list[i] if i < len(modify_temp_list) else modify_temp_list[-1]),
                'temp_factor': float(temp_factor_list[i] if i < len(temp_factor_list) else temp_factor_list[-1]),
            }
            disk.add_pillar(**pillar)

    # Get orbital period and emission line parameters
    line_params = config.get('emission_line', {})
    if 'M_BH' in line_params:
        M_BH = float(line_params['M_BH'])
        v_virial = v_virial_from_mbh(M_BH, disk.r0)
        print(f"  M_BH = {M_BH:.2e} Msun → v_virial = {v_virial:.1f} km/s at r0 = {disk.r0:.1f} ld")
    else:
        v_virial = float(line_params.get('v_virial', 1000.0))
    # Use mean pillar radius for orbital period (determines movie duration)
    if len(disk.pillars) > 0:
        r_mean_pillars = np.mean([p['r'] for p in disk.pillars])
    else:
        r_mean_pillars = disk.r0
    T_orbital = compute_orbital_period(r_mean_pillars, v_virial, disk.r0)

    # Save original pillar positions
    original_pillars = [p.copy() for p in disk.pillars]
    original_phi = original_pillars[0]['phi'] if len(original_pillars) > 0 else 0

    # All pillars rotate together as disk rotates
    rotate_pillars = True

    # Emission line rest wavelengths. Hβ ("HI") is included so the
    # third-row colour-dependent diagnostic Δlog(Ψ_line/Ψ_Hβ) can be
    # built from the same set of Cloudy-weighted maps.
    lambda0_dict = {
        'Halpha': 6562.8,
        'Mg2': 2798.0,
        'C4': 1549.0,
        'HI': 4861.32,
    }

    # Update from config if provided
    if 'lambda0_Halpha' in line_params:
        lambda0_dict['Halpha'] = float(line_params['lambda0_Halpha'])
    if 'lambda0_Mg2' in line_params:
        lambda0_dict['Mg2'] = float(line_params['lambda0_Mg2'])
    if 'lambda0_C4' in line_params:
        lambda0_dict['C4'] = float(line_params['lambda0_C4'])

    # Velocity-delay map parameters
    comp_params = config.get('computation', {})
    nlambda = int(comp_params.get('nlambda', 200))
    ntau = int(comp_params.get('ntau_line', 100))
    taumax = float(comp_params.get('taumax', 150.0))

    # Cloudy model file - load ONCE here, pass raw data to workers
    cloudy_extension_file = None
    if use_absolute_flux:
        cloudy_file = '/Users/jiamuh/c23.01/my_models/loc_metal_flux/strong_LOC_varym_N25_v100_lineflux_LineList_BLR_Fe2_flux.txt'
        cloudy_extension_file = '/Users/jiamuh/c23.01/my_models/loc_metal_flux/strong_LOC_varym_N25_v100_lineflux_extlow_LineList_BLR_Fe2_flux.txt'
    else:
        cloudy_file = '/Users/jiamuh/c23.01/my_models/loc_metal/strong_LOC_varym_N25_LineList_BLR_Fe2.txt'
    Z_target = 1.0

    # Load Cloudy raw data ONCE (collapsed arrays, picklable). Splice the
    # low-phi extension grid in so the spline reaches down to log Phi ~ 15
    # and shadow cells are not clipped to the default [17, 21] floor.
    cloudy_data_raw, phi_grid_cloudy, Z_grid_cloudy = _load_cloudy_raw_data(
        cloudy_file, Z_target, extension_file=cloudy_extension_file)

    # Time grid for one orbital period
    time_grid = np.linspace(0, T_orbital, n_frames, endpoint=False)

    print(f"Generating {n_frames} movie frames with Cloudy-weighted velocity-delay maps...")
    print(f"  Mean pillar radius: {r_mean_pillars:.1f} ld, orbital period: {T_orbital:.1f} days")
    print(f"  Each pillar rotates at its own Keplerian velocity")
    print(f"  Lines: {list(lambda0_dict.keys())}")
    print(f"  Velocity-delay map: {nlambda} wavelength bins, {ntau} delay bins, taumax={taumax} days")
    print(f"  Output directory: {output_dir}/")

    # Serialize disk parameters for multiprocessing
    disk_params_serial = {
        'rin': disk.rin, 'rout': disk.rout, 'nr': disk.nr, 'nphi': disk.nphi,
        'h1': disk.h1, 'r0': disk.r0, 'beta': disk.beta,
        'tv1': disk.tv1, 'alpha': disk.alpha, 'hlamp': disk.hlamp,
        'tx1': disk.tx1, 'tirrad_tvisc_ratio': disk.tirrad_tvisc_ratio,
        'dmpc': disk.dmpc, 'cosi': disk.cosi, 'fcol': disk.fcol, 'redshift': disk.redshift,
        'flux_method': disk.flux_method, 'log_phi_inner': disk.log_phi_inner,
        'no_fluxfloor': disk.no_fluxfloor,
        'transparent': disk.transparent, 'f_trans': disk.f_trans,
        'f_turb_line': disk.f_turb_line,
    }

    # ========== Compute Cloudy-weighted velocity-delay maps WITHOUT pillars (only once) ==========
    print("\nComputing Cloudy-weighted velocity-delay maps without pillars (baseline)...")
    saved_pillars = disk.pillars.copy()
    disk.pillars = []

    # Parallel computation for baseline (no pillars)
    args_list = [(line_key, lambda0, disk_params_serial, [], v_virial, nlambda, ntau, taumax,
                  cloudy_data_raw, phi_grid_cloudy, Z_grid_cloudy, Z_target, weighting)
                 for line_key, lambda0 in lambda0_dict.items()]

    psi_no_pillars = {}
    lambda_grids = {}
    tau_grid = None

    with Pool(processes=min(3, cpu_count())) as pool:
        results = pool.map(_compute_vdmap_for_line_cloudy, args_list)
    for line_key, lambda_grid, tg, psi_map in results:
        psi_no_pillars[line_key] = psi_map
        lambda_grids[line_key] = lambda_grid
        tau_grid = tg

    disk.pillars = saved_pillars

    # Line labels for plot titles
    line_labels = {'Halpha': r'$\rm H\alpha$', 'Mg2': r'$\rm MgII$', 'C4': r'$\rm CIV$'}

    # Compute color scale limits from first frame
    print("\nEstimating color scale from first frame...")

    # Serialize pillars for parallel computation
    pillars_serial = [{'r': p['r'], 'phi': p['phi'], 'height': p['height'],
                       'sigma_r': p['sigma_r'], 'sigma_phi': p['sigma_phi']}
                      for p in disk.pillars]

    args_list = [(line_key, lambda0, disk_params_serial, pillars_serial, v_virial, nlambda, ntau, taumax,
                  cloudy_data_raw, phi_grid_cloudy, Z_grid_cloudy, Z_target, weighting)
                 for line_key, lambda0 in lambda0_dict.items()]

    with Pool(processes=min(3, cpu_count())) as pool:
        results = pool.map(_compute_vdmap_for_line_cloudy, args_list)

    # Compute SHARED color scales across all lines
    psi_max_global = 0
    diff_max_global = 0
    for line_key, _, _, psi_map in results:
        psi_max_global = max(psi_max_global, np.nanmax(psi_map), np.nanmax(psi_no_pillars[line_key]))
        diff_map = psi_map - psi_no_pillars[line_key]
        diff_max_global = max(diff_max_global, np.nanmax(np.abs(diff_map)))

    psi_vmax = psi_max_global * 1.1
    diff_vmax = diff_max_global * 1.1
    print(f"  Shared color scales: psi_vmax={psi_vmax:.2e}, diff_vmax={diff_vmax:.2e}")

    # ========== Generate frames ==========
    print(f"\nGenerating {n_frames} movie frames...")

    from matplotlib import cm
    from matplotlib.gridspec import GridSpec

    for iframe in range(n_frames):
        t = time_grid[iframe]

        # Rotate each pillar at its own Keplerian velocity: omega(r) = v_virial * sqrt(r0) / r^(3/2)
        # v(r) = v_virial * sqrt(r0/r), omega = v/r
        if rotate_pillars:
            for i, p in enumerate(disk.pillars):
                r_p = p['r']
                # Keplerian angular velocity at this radius
                v_kep = v_virial * np.sqrt(disk.r0 / r_p)  # km/s
                # Convert to rad/day: omega = v / r, with r in light-days and v in km/s
                # 1 light-day = 2.59e10 km, so omega = v [km/s] / (r [ld] * 2.59e10 [km/ld]) * 86400 [s/day]
                omega_rad_per_day = v_kep / (r_p * 2.59e10) * 86400.0
                phi_offset = omega_rad_per_day * t
                p['phi'] = original_pillars[i]['phi'] + phi_offset

        # ===== Compute Cloudy-weighted velocity-delay maps WITH pillars (parallel) =====
        pillars_serial = [{'r': p['r'], 'phi': p['phi'], 'height': p['height'],
                           'sigma_r': p['sigma_r'], 'sigma_phi': p['sigma_phi']}
                          for p in disk.pillars]

        args_list = [(line_key, lambda0, disk_params_serial, pillars_serial, v_virial, nlambda, ntau, taumax,
                      cloudy_data_raw, phi_grid_cloudy, Z_grid_cloudy, Z_target, weighting)
                     for line_key, lambda0 in lambda0_dict.items()]

        with Pool(processes=min(3, cpu_count())) as pool:
            results = pool.map(_compute_vdmap_for_line_cloudy, args_list)

        psi_with_pillars = {}
        for line_key, _, _, psi_map in results:
            psi_with_pillars[line_key] = psi_map

        # Compute difference maps
        psi_diff = {}
        for line_key in lambda0_dict.keys():
            psi_diff[line_key] = psi_with_pillars[line_key] - psi_no_pillars[line_key]

        # Layout depends on skip_geometry. With skip_geometry=True we
        # show: row 0 = ψ density, row 1 = ΔΨ/Ψ₀ residual,
        # row 2 = Δlog(Ψ_line/Ψ_Hβ) (colour-dependent diagnostic).
        if skip_geometry:
            fig = plt.figure(figsize=(15, 12))
            gs = GridSpec(3, 3, figure=fig, height_ratios=[1, 1, 1],
                          wspace=0.35, hspace=0.10,
                          top=0.94, bottom=0.07, left=0.06, right=0.94)
        else:
            fig = plt.figure(figsize=(15, 14))
            gs = GridSpec(3, 3, figure=fig, height_ratios=[1, 1, 1.2],
                          wspace=0.35, hspace=0.25,
                          top=0.95, bottom=0.05, left=0.06, right=0.94)

            # ===== 3D geometry at bottom (Row 2) =====
            ax_3d = fig.add_subplot(gs[2, :], projection='3d')

            n_phi_plot, n_r_plot = 180, 100
            r_plot = np.logspace(np.log10(disk.rin), np.log10(disk.rout), n_r_plot)
            phi_plot = np.linspace(0, 2*np.pi, n_phi_plot, endpoint=False)
            r_2d, phi_2d = np.meshgrid(r_plot, phi_plot, indexing='ij')
            h_2d = disk.get_height(r_2d, phi_2d)

            phi_2d_closed = np.concatenate([phi_2d, phi_2d[:, :1] + 2*np.pi], axis=1)
            r_2d_closed = np.concatenate([r_2d, r_2d[:, :1]], axis=1)
            h_2d_closed = np.concatenate([h_2d, h_2d[:, :1]], axis=1)

            x_2d = r_2d_closed * np.cos(phi_2d_closed)
            y_2d = r_2d_closed * np.sin(phi_2d_closed)
            z_2d = h_2d_closed

            log_phi_2d = disk.compute_log_ionizing_flux(r_2d_closed, phi_2d_closed)
            finite_lp = log_phi_2d[np.isfinite(log_phi_2d) & (log_phi_2d > 0)]
            phi_min_color = np.floor(np.min(finite_lp) * 2) / 2
            phi_max_color = np.ceil(np.max(finite_lp) * 2) / 2

            h_base = disk.h1 * (r_2d_closed / disk.r0) ** disk.beta
            pillar_mask = (h_2d_closed - h_base) > 0.001
            face_dot_obs = np.cos(phi_2d_closed) * disk.sini
            blend = 0.5 * (1.0 + np.tanh(face_dot_obs / 0.3))
            shadow_mask = disk._compute_shadow_mask(r_2d_closed, phi_2d_closed, h_2d_closed)
            shadowed_flux = np.where(shadow_mask < 0.5, log_phi_2d, phi_max_color)
            shadow_level = np.min(shadowed_flux)
            log_phi_view = np.where(pillar_mask,
                                    (1.0 - blend) * log_phi_2d + blend * shadow_level,
                                    log_phi_2d)
            log_phi_norm = np.clip((log_phi_view - phi_min_color) / (phi_max_color - phi_min_color), 0, 1)

            plot_colors = cm.inferno(log_phi_norm)
            plot_colors[:, :, 3] = 1.0

            ax_3d.plot_surface(x_2d, y_2d, z_2d, facecolors=plot_colors,
                              shade=False, antialiased=True, rstride=2, cstride=2)

            ax_3d.scatter([0], [0], [disk.hlamp], color='yellow', s=300, marker='*',
                          edgecolors='orange', linewidths=1.5, zorder=10)

            ax_3d.set_xlabel(r'$X$ [ld]', fontsize=12, labelpad=10)
            ax_3d.set_ylabel(r'$Y$ [ld]', fontsize=12, labelpad=10)
            ax_3d.set_zlabel(r'$Z$ [ld]', fontsize=12, labelpad=5)
            ax_3d.set_xlim([-disk.rout * 1.1, disk.rout * 1.1])
            ax_3d.set_ylim([-disk.rout * 1.1, disk.rout * 1.1])
            max_pillar_height = max([p['height'] for p in disk.pillars]) if len(disk.pillars) > 0 else 0.1
            ax_3d.set_zlim([0, max(0.15, max_pillar_height * 1.5)])

            ax_3d.tick_params(axis='x', labelsize=9, pad=1)
            ax_3d.tick_params(axis='y', labelsize=9, pad=1)
            ax_3d.tick_params(axis='z', labelsize=9, pad=1)

            inclination_deg = np.degrees(np.arccos(disk.cosi))
            ax_3d.view_init(elev=max(50, 90 - inclination_deg), azim=0)
            ax_3d.dist = 6
            ax_3d.grid(False)
            ax_3d.xaxis.pane.fill = False
            ax_3d.yaxis.pane.fill = False
            ax_3d.zaxis.pane.fill = False

        from matplotlib.ticker import AutoMinorLocator

        # Helper: add velocity twin axis with auto-placed round tick values
        def add_velocity_axis(ax_main, lambda_grid, lambda0, show_label=True):
            lam2v = lambda lam: (C / KM_TO_CM) * (lam - lambda0) / lambda0
            v2lam = lambda v: lambda0 * (1.0 + v * KM_TO_CM / C)
            ax2 = ax_main.secondary_xaxis('top', functions=(lam2v, v2lam))
            ax2.xaxis.set_major_locator(plt.MaxNLocator(nbins=5, integer=True, steps=[1, 2, 3, 5, 6, 10]))
            ax2.tick_params(axis='both', which='major', direction='in', length=8, width=1.5, labelsize=13)
            ax2.tick_params(axis='both', which='minor', direction='in', length=5, width=1.0)
            ax2.xaxis.set_minor_locator(AutoMinorLocator())
            if show_label:
                ax2.set_xlabel(r'$v~\rm [km\,s^{-1}]$', fontsize=14, labelpad=6)
            else:
                ax2.set_xlabel('')
                ax2.tick_params(labeltop=False)
            return ax2

        # Helper: style axis with inward ticks and minor ticks
        def style_axis(ax):
            ax.tick_params(axis='both', which='major', direction='in', length=8, width=1.5, labelsize=13)
            ax.tick_params(axis='both', which='minor', direction='in', length=5, width=1.0)
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())

        # ===== Per-line 2D probability density normalisation =====
        # We display Ψ as a 2D probability density ψ(λ, τ) = Ψ / ∬ Ψ dλ dτ
        # so that each line panel integrates to unity (units: AA^-1 day^-1).
        # The residual is computed from these normalised maps as well, so
        # the row-1 colorbar reflects the *shape change* in each line's
        # response under a fixed total normalisation, independent of the
        # very different absolute amplitudes between Hα/MgII/CIV.
        line_keys = ['C4', 'Mg2', 'Halpha']  # ordered by increasing rest wavelength

        def _to_density(psi_map, lambda_grid_loc, tau_grid_loc):
            dl = float(lambda_grid_loc[1] - lambda_grid_loc[0])
            dt = float(tau_grid_loc[1] - tau_grid_loc[0])
            tot = float(np.nansum(psi_map)) * dl * dt
            if tot <= 0 or not np.isfinite(tot):
                return np.zeros_like(psi_map)
            return psi_map / tot

        psi_with_density = {
            k: _to_density(psi_with_pillars[k], lambda_grids[k], tau_grid)
            for k in line_keys
        }
        psi_no_density = {
            k: _to_density(psi_no_pillars[k], lambda_grids[k], tau_grid)
            for k in line_keys
        }

        # ===== Row 0: Probability density ψ(λ, τ) =====
        # Shared, fixed vmax across all three lines so the panels are
        # directly comparable; CIV's peak will saturate (it's larger than
        # this fixed scale), making the bowl-envelope features in
        # Hα/MgII clearly visible.
        psi_vmax_fixed = 0.001  # AA^-1 d^-1
        for idx, line_key in enumerate(line_keys):
            ax = fig.add_subplot(gs[0, idx])

            lambda_grid = lambda_grids[line_key]
            psi_map = psi_with_density[line_key]
            lambda0 = lambda0_dict[line_key]

            im = ax.pcolormesh(lambda_grid, tau_grid, psi_map, shading='auto',
                              cmap=VDM_CMAP, vmin=0, vmax=psi_vmax_fixed)

            ax.axvline(lambda0, color='red', linestyle='--', linewidth=1.5, alpha=0.8)
            # Top row: no bottom x-axis label/ticks (shared with bottom row)
            ax.set_xticklabels([])
            ax.set_xlabel('')
            ax.set_ylabel(r'$\rm Delay~[days]$', fontsize=14)
            ax.set_title(line_labels[line_key], fontsize=14)
            ax.set_xlim(lambda_grid[0], lambda_grid[-1])
            ax.set_ylim(0, taumax)
            style_axis(ax)

            cbar = plt.colorbar(im, ax=ax, pad=0.02)
            cbar.set_label(r'$\psi~[\rm \AA^{-1}\,d^{-1}]$', fontsize=13)
            cbar.ax.tick_params(labelsize=12, direction='in', length=5, width=1.0)

            add_velocity_axis(ax, lambda_grid, lambda0, show_label=True)

        # ===== Row 1: Fractional difference maps =====
        # Computed from the ABSOLUTE Ψ (not the per-line densities), so
        # the residual reflects only local redistribution from the
        # pillars, not the global re-normalisation shift that would arise
        # from normalising each Ψ to its own total integral.
        for idx, line_key in enumerate(line_keys):
            ax = fig.add_subplot(gs[1, idx])

            lambda_grid = lambda_grids[line_key]
            psi_w_abs = psi_with_pillars[line_key]
            baseline = psi_no_pillars[line_key]
            lambda0 = lambda0_dict[line_key]

            # Δ log Ψ = log10(Ψ_with / Ψ_without). Symmetric in form
            # with row 3's Δlog(Ψ/Ψ_Hβ); both rows now read as
            # logarithmic differences. Mask cells where either map is
            # non-positive to avoid log(0) / log(negative).
            with np.errstate(divide='ignore', invalid='ignore'):
                delta_log_psi = np.log10(psi_w_abs / baseline)
            valid = (psi_w_abs > 0) & (baseline > 0)
            delta_log_psi = np.where(valid, delta_log_psi, np.nan)

            im = ax.pcolormesh(lambda_grid, tau_grid, delta_log_psi,
                               shading='auto', cmap='RdBu', vmin=-2, vmax=2)

            ax.axvline(lambda0, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
            # Middle row: no bottom x-axis label (shared with row 3 below)
            ax.set_xticklabels([])
            ax.set_xlabel('')
            ax.set_ylabel(r'$\rm Delay~[days]$', fontsize=14)
            ax.set_xlim(lambda_grid[0], lambda_grid[-1])
            ax.set_ylim(0, taumax)
            style_axis(ax)

            cbar = plt.colorbar(im, ax=ax, pad=0.02)
            cbar.set_label(r'$\Delta\log\Psi$', fontsize=13)
            cbar.ax.tick_params(labelsize=12, direction='in', length=5, width=1.0)

        # ===== Row 2: Δ log(Ψ_line / Ψ_Hβ) — colour-dependent diagnostic =====
        # For each of Hα, MgII, CIV we form the change in log(line / Hβ)
        # under the pillar perturbation. Hβ is a recombination reference,
        # so this map highlights the colour-dependent shadow response:
        # MgII gets boosted in shadow (low U → MgII enhanced relative to
        # recombination), CIV collapses in shadow.
        psi_Hb_w = psi_with_pillars.get('HI')
        psi_Hb_n = psi_no_pillars.get('HI')
        for idx, line_key in enumerate(line_keys):
            ax = fig.add_subplot(gs[2, idx])
            lambda_grid = lambda_grids[line_key]
            psi_line_w = psi_with_pillars[line_key]
            psi_line_n = psi_no_pillars[line_key]
            lambda0 = lambda0_dict[line_key]

            # Both ψ maps share the same (velocity, τ) grid even though
            # the absolute λ grids differ: the ratio is therefore well-
            # defined per (ibin, ilambda). Mask any cell where Hβ or the
            # line's baseline is non-positive.
            with np.errstate(divide='ignore', invalid='ignore'):
                log_ratio_w = np.log10(psi_line_w / psi_Hb_w)
                log_ratio_n = np.log10(psi_line_n / psi_Hb_n)
                delta_log = log_ratio_w - log_ratio_n
            valid = (psi_line_n > 0) & (psi_Hb_n > 0) & (psi_line_w > 0) & (psi_Hb_w > 0)
            delta_log = np.where(valid, delta_log, np.nan)

            im = ax.pcolormesh(lambda_grid, tau_grid, delta_log,
                               shading='auto', cmap='RdBu',
                               vmin=-1.0, vmax=1.0)

            ax.axvline(lambda0, color='black', linestyle='--',
                       linewidth=1.5, alpha=0.7)
            ax.set_xlabel(r'$\rm Wavelength~[\AA]$', fontsize=14)
            ax.set_ylabel(r'$\rm Delay~[days]$', fontsize=14)
            ax.set_xlim(lambda_grid[0], lambda_grid[-1])
            ax.set_ylim(0, taumax)
            style_axis(ax)

            cbar = plt.colorbar(im, ax=ax, pad=0.02)
            cbar.set_label(r'$\Delta\log(\Psi/\Psi_{\rm H\beta})$', fontsize=13)
            cbar.ax.tick_params(labelsize=12, direction='in', length=5, width=1.0)

            # Mirror the existing convention: top tick marks but no labels
            ax2 = ax.twiny()
            ax2.set_xlim(ax.get_xlim())
            ax2.set_xticklabels([])
            ax2.tick_params(axis='x', which='major', direction='in', length=8, width=1.5)
            ax2.tick_params(axis='x', which='minor', direction='in', length=5, width=1.0)
            ax2.xaxis.set_minor_locator(AutoMinorLocator())

        # Save frame
        frame_file = os.path.join(output_dir, f'frame_{iframe:04d}.png')
        plt.savefig(frame_file, dpi=200, bbox_inches='tight', pad_inches=0.2, facecolor='white')
        plt.close()

        if (iframe + 1) % 10 == 0:
            print(f"  Generated frame {iframe+1}/{n_frames}")

    # Restore original pillar positions
    for i, p in enumerate(disk.pillars):
        p['phi'] = original_pillars[i]['phi']

    print(f"\nMovie frames saved to {output_dir}/")
    print(f"To create movie, run:")
    print(f"  ffmpeg -framerate 10 -i {output_dir}/frame_%04d.png -c:v libx264 -pix_fmt yuv420p movie.mp4")


if __name__ == '__main__':
    import sys
    import argparse
    parser = argparse.ArgumentParser(description='Time-evolving emission line computation with Cloudy models')
    parser.add_argument('config', nargs='?', default='config_line.yaml', help='Config file (default: config_line.yaml)')
    parser.add_argument('--absolute-flux', action='store_true', help='Use absolute flux instead of normalized')
    parser.add_argument('--plot-only', action='store_true', help='Only plot geometry, skip spectral computation')
    parser.add_argument('--movie', action='store_true', help='Generate movie frames')
    parser.add_argument('--n-frames', type=int, default=50, help='Number of movie frames (default: 50)')
    args = parser.parse_args()

    if args.movie:
        generate_movie_frames(config_file=args.config, n_frames=args.n_frames,
                              use_absolute_flux=args.absolute_flux)
    else:
        main(config_file=args.config, use_absolute_flux=args.absolute_flux, plot_only=args.plot_only)
