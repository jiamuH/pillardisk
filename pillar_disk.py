#!/usr/bin/env python3
"""
pillar_disk.py - Compute time delay and SED for AGN accretion disk with Gaussian pillar bumps

This script places Gaussian bumps (pillars) in an axisymmetric accretion disk at fixed
r and phi positions, then predicts the expected lag time spectrum and SED.

The disk is still assumed to be axisymmetric overall, but the pillars are not.

Requirements:
    numpy
    matplotlib
    scipy (optional, for advanced features)
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple, Optional
import os
plt.rcParams.update({
    'text.usetex': True,
    'axes.linewidth': 2,
    'font.family': 'serif',
    'font.weight': 'heavy',
    'font.size': 20,
})

plt.rcParams['text.latex.preamble'] = r'\usepackage{bm} \boldmath'
# Try to import tqdm for progress bars
try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False
    # Dummy tqdm that does nothing
    def tqdm(iterable, *args, **kwargs):
        return iterable

# Try to import 3D plotting
try:
    from mpl_toolkits.mplot3d import Axes3D
    HAS_3D = True
except ImportError:
    HAS_3D = False

# Try to import multiprocessing for parallelization
try:
    from multiprocessing import Pool, cpu_count
    HAS_MP = True
except ImportError:
    HAS_MP = False

# Try to import YAML for config file
try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False
    print("Warning: PyYAML not installed. Install with 'pip install pyyaml' to use config files.")


def load_config(config_file='config.yaml'):
    """
    Load configuration from YAML file.
    
    Parameters:
    -----------
    config_file : str
        Path to configuration file
    
    Returns:
    --------
    config : dict
        Configuration dictionary
    """
    if not HAS_YAML:
        print("Warning: YAML not available. Using default parameters.")
        return None
    
    if not os.path.exists(config_file):
        print(f"Warning: Config file '{config_file}' not found. Using default parameters.")
        return None
    
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    
    return config


def resolve_rin(disk_params):
    """Compute rin from M_BH and r_isco_rg if rin='auto', then clean up extra keys."""
    rin_val = disk_params.get('rin', 0.1)
    if isinstance(rin_val, str) and rin_val.lower() == 'auto':
        M_BH = float(disk_params.get('M_BH', 7e7))
        r_isco_rg = float(disk_params.get('r_isco_rg', 6))
        G = 6.674e-11; M_sun = 1.989e30; c = 2.998e8; ld = c * 86400.0
        r_g = G * M_BH * M_sun / c**2 / ld  # light days
        disk_params['rin'] = r_isco_rg * r_g
        print(f"  r_g = {r_g:.5f} ld, r_ISCO = {r_isco_rg} r_g = {disk_params['rin']:.4f} ld")
    disk_params.pop('M_BH', None)
    disk_params.pop('r_isco_rg', None)


# Physical constants (cgs units)
C = 2.997925e10  # speed of light (cm/s)
H = 6.6262e-27   # Planck constant (erg s)
K = 1.3806e-16   # Boltzmann constant (erg/K)
SIGMA = 5.66956e-5  # Stefan-Boltzmann constant (erg/cm^2/s/K^4)
PC = 3.0857e18   # parsec (cm)
DAY = 24. * 3600.  # day (s)
ANGSTROM = 1e-8   # Angstrom (cm)

# Conversion factors
PC_TO_LD = 1e6 * PC / (C * DAY)  # parsec to light days
FNU_TO_MJY = 1e26  # erg/cm^2/s/Hz to mJy


class PillarDisk:
    """
    Accretion disk with Gaussian pillar bumps for AGN continuum modeling.
    """
    
    def __init__(self,
                 rin: float = 0.1,      # inner radius (light days)
                 rout: float = 100.0,   # outer radius (light days)
                 nr: int = 300,         # number of radial points
                 nphi: int = 360,       # number of azimuthal points
                 h1: float = 0.01,      # height at reference radius (light days)
                 r0: float = 10.0,      # reference radius (light days)
                 beta: float = 10.0,     # height power-law index h = h1 * (r/r0)^beta
                 tv1: float = 1e4,      # viscous temperature at r0 (K)
                 alpha: float = 0.75,   # temperature power-law index T = tv1 * (r0/r)^alpha
                 hlamp: float = 0.20,    # lamp height (light days)
                 tx1: float = 1e4,      # irradiation temperature at r0 (K) (deprecated, use tirrad_tvisc_ratio instead)
                 tirrad_tvisc_ratio: float = 10.0,  # ratio T_irrad/T_visc at r0 (both follow r^-3/4)
                 dmpc: float = 100.0,   # distance (Mpc)
                 cosi: float = 0.5,      # cos(inclination), 1=face-on
                 fcol: float = 1.0,     # color temperature factor Tcol = fcol * Teff
                 redshift: float = 0.0, # redshift
                 flux_method: str = 'scale_temp',  # 'scale_temp' or 'direct'
                 log_phi_inner: float = 24.0,  # target log10(Phi) at r_in for Q calibration
                 no_fluxfloor: bool = False,  # if True, shadow cells below Cloudy grid get zero emissivity
                 transparent: bool = False,  # if True, pillars transmit a fraction f_trans of lamp flux
                 f_trans: float = 0.0,       # transmission fraction (0=opaque, 1=fully transparent)
                 f_turb_line: float = 0.0):  # turbulent broadening: sigma_v = f_turb * v_Kep(r)
        """
        Initialize the disk model.
        """
        self.rin = rin
        self.rout = rout
        self.nr = nr
        self.nphi = nphi
        self.h1 = h1
        self.r0 = r0
        self.beta = beta
        self.tv1 = tv1
        self.alpha = alpha
        self.hlamp = hlamp
        self.tx1 = tx1
        self.tirrad_tvisc_ratio = tirrad_tvisc_ratio
        self.dmpc = dmpc
        self.cosi = cosi
        self.sini = np.sqrt(1.0 - cosi**2)
        self.fcol = fcol
        self.redshift = redshift
        self.flux_method = flux_method
        self.log_phi_inner = log_phi_inner
        self.no_fluxfloor = no_fluxfloor
        self.transparent = transparent
        self.f_trans = f_trans if transparent else 0.0
        self.f_turb_line = f_turb_line

        # Distance in light days
        self.d = dmpc * PC_TO_LD
        
        # Unit vector pointing to Earth (in x-z plane)
        self.ex = self.sini
        self.ey = 0.0
        self.ez = self.cosi
        
        # Pillar list: [(r_pillar, phi_pillar, height, sigma_r, sigma_phi), ...]
        self.pillars = []
        
        # Store parameters for easy access/modification
        self.params = {
            'rin': rin, 'rout': rout, 'nr': nr, 'nphi': nphi,
            'h1': h1, 'r0': r0, 'beta': beta,
            'tv1': tv1, 'alpha': alpha,
            'hlamp': hlamp, 'tx1': tx1, 'tirrad_tvisc_ratio': tirrad_tvisc_ratio,
            'dmpc': dmpc, 'cosi': cosi, 'fcol': fcol, 'redshift': redshift,
            'flux_method': flux_method, 'log_phi_inner': log_phi_inner,
            'no_fluxfloor': no_fluxfloor,
            'transparent': transparent, 'f_trans': f_trans,
            'f_turb_line': f_turb_line
        }
        
        # Set up radial and azimuthal grids
        self._setup_grids()
        
        # Compute base disk geometry and temperature
        self._compute_base_disk()
    
    def update_params(self, **kwargs):
        """
        Update disk parameters and recompute geometry.
        
        Example:
            disk.update_params(nr=200, nphi=180, beta=1.5)
        """
        # Update stored parameters
        for key, value in kwargs.items():
            if key in self.params:
                self.params[key] = value
                setattr(self, key, value)
            else:
                print(f"Warning: Unknown parameter '{key}'")
        
        # Recompute grids and geometry
        self._setup_grids()
        self._compute_base_disk()
    
    def _setup_grids(self):
        """Set up radial and azimuthal grids."""
        # Log-spaced radial grid
        rlog1 = np.log10(self.rin)
        rlog2 = np.log10(self.rout)
        self.r = np.logspace(rlog1, rlog2, self.nr)
        
        # Azimuthal grid (0 to 2pi)
        self.phi = np.linspace(0, 2*np.pi, self.nphi, endpoint=False)
        self.dphi = 2*np.pi / self.nphi
        
        # Radial spacing (for integration)
        self.dr = np.diff(self.r)
        self.dr = np.append(self.dr, self.dr[-1])  # extend last value
    
    def _compute_base_disk(self):
        """Compute base axisymmetric disk geometry and temperature."""
        # Base height profile (axisymmetric)
        x = self.r / self.r0
        self.h_base = self.h1 * x**self.beta
        
        # Base viscous temperature profile: T_visc ∝ r^(-3/4)
        # Use alpha = 0.75 for standard disk
        # Apply no-torque boundary condition at ISCO (rin)
        # Temperature should go to zero at rin due to no-torque condition
        # Use prescription: T_visc ∝ sqrt((r - rin) / r) * (r0/r)^0.75
        # This ensures T → 0 as r → rin
        r_factor = (self.r0 / self.r)**0.75
        # No-torque factor: sqrt((r - rin) / r) for r > rin, 0 for r <= rin
        no_torque_factor = np.sqrt(np.maximum(0.0, (self.r - self.rin) / self.r))
        self.tv_base = self.tv1 * r_factor * no_torque_factor
        
        # Irradiation temperature: T_irrad ∝ r^(-3/4), scaled by tirrad_tvisc_ratio
        # T_irrad = tirrad_tvisc_ratio * T_visc at r0
        # Also apply no-torque factor to prevent high temperatures at small r
        self.tx_base = self.tv1 * self.tirrad_tvisc_ratio * r_factor * no_torque_factor
        
        # Total effective temperature (T^4 = Tv^4 + Tx^4)
        self.t_base = (self.tv_base**4 + self.tx_base**4)**0.25
    
    def _compute_irradiation(self):
        """
        Compute irradiation temperature from lamp.
        
        Note: This method is now deprecated. Irradiation temperature is computed
        directly in _compute_base_disk() using tirrad_tvisc_ratio.
        This method is kept for backward compatibility but does nothing.
        """
        # Irradiation temperature is now computed in _compute_base_disk()
        # as T_irrad = tirrad_tvisc_ratio * T_visc, both following r^(-3/4)
        pass
    
    def add_pillar(self, r_pillar: float, phi_pillar: float, 
                   height: float = 0.01, sigma_r: float = 1.0, 
                   sigma_phi: float = 0.1, modify_height: bool = True,
                   modify_temp: bool = False, temp_factor: float = 1.5):
        """
        Add a Gaussian bump (pillar) to the disk.
        
        Parameters:
        -----------
        r_pillar : float
            Radial position of pillar center (light days)
        phi_pillar : float
            Azimuthal position of pillar center (radians)
        height : float
            Height of the 2D Gaussian bump at pillar center (light days)
        sigma_r : float
            Radial width of Gaussian (light days)
        sigma_phi : float
            Azimuthal width of Gaussian (radians)
        modify_height : bool
            If True, modify disk height at pillar location
        modify_temp : bool
            If True, modify disk temperature at pillar location
        temp_factor : float
            Temperature enhancement factor (deprecated - now computed from geometry/angle)
            Kept for backward compatibility but not used in calculation
        """
        self.pillars.append({
            'r': r_pillar,
            'phi': phi_pillar,
            'height': height,
            'sigma_r': sigma_r,
            'sigma_phi': sigma_phi,
            'modify_height': modify_height,
            'modify_temp': modify_temp,
            'temp_factor': temp_factor
        })
    
    def get_height(self, r: np.ndarray, phi: np.ndarray) -> np.ndarray:
        """
        Get disk height at given (r, phi) positions, including pillars.
        
        Parameters:
        -----------
        r : array
            Radial positions (light days), can be 1D or 2D
        phi : array
            Azimuthal positions (radians), can be 1D or 2D
        
        Returns:
        --------
        h : array
            Height above midplane (light days)
        """
        # Handle both 1D and 2D input arrays
        r = np.asarray(r)
        phi = np.asarray(phi)
        
        # If 1D arrays, create meshgrid
        if r.ndim == 1 and phi.ndim == 1:
            r_2d, phi_2d = np.meshgrid(r, phi, indexing='ij')
        elif r.ndim == 2 and phi.ndim == 2:
            r_2d = r
            phi_2d = phi
        else:
            # Try to broadcast
            r_2d, phi_2d = np.broadcast_arrays(r, phi)
        
        # Base axisymmetric height
        h = np.interp(r_2d.flatten(), self.r, self.h_base).reshape(r_2d.shape)
        
        # Add pillar contributions
        for pillar in self.pillars:
            # Distance from pillar center
            dr = r_2d - pillar['r']
            dphi = phi_2d - pillar['phi']
            # Wrap azimuthal difference to [-pi, pi]
            dphi = np.mod(dphi + np.pi, 2*np.pi) - np.pi
            
            # Gaussian bump with periodic azimuthal wrapping (sum over images)
            sig_phi = pillar['sigma_phi']
            gauss_phi = (np.exp(-0.5 * (dphi / sig_phi)**2)
                       + np.exp(-0.5 * ((dphi - 2*np.pi) / sig_phi)**2)
                       + np.exp(-0.5 * ((dphi + 2*np.pi) / sig_phi)**2))
            gauss = np.exp(-0.5 * (dr / pillar['sigma_r'])**2) * gauss_phi
            
            if pillar['modify_height']:
                # Add height directly (in light days) weighted by Gaussian
                h_bump = pillar['height'] * gauss
                h = h + h_bump
        
        return h
    
    def _compute_shadow_mask(self, r: np.ndarray, phi: np.ndarray, h: np.ndarray) -> np.ndarray:
        """
        Compute shadow mask: 1 if point is illuminated, 0 if in shadow.

        Since the lamp is at (0, 0, hlamp), a ray from the lamp to any disk
        point at (r, phi) travels radially outward at constant azimuth phi.
        At intermediate radius r_p, the ray height is:
            z_ray = hlamp + (r_p / r) * (h_disk - hlamp)

        A pillar at (r_p, phi_p) shadows the disk point if the pillar height
        at (r_p, phi) exceeds z_ray.  The pillar height profile at its own
        radius (dr=0) is just the azimuthal Gaussian:
            h_pillar(phi) = h_base(r_p) + height * exp(-0.5*((phi-phi_p)/sigma_phi)^2)

        Parameters:
        -----------
        r : array
            Radial positions (light days), 2D
        phi : array
            Azimuthal positions (radians), 2D
        h : array
            Heights at (r, phi), 2D

        Returns:
        --------
        shadow_mask : array
            Mask: 1.0 = illuminated, 0.0 = in shadow (2D, same shape as r)
        """
        shadow_mask = np.ones_like(r)

        if len(self.pillars) == 0:
            return shadow_mask

        for pillar in self.pillars:
            if not pillar['modify_height']:
                continue

            r_p = pillar['r']
            phi_p = pillar['phi']

            # Only points beyond the pillar radius can be shadowed
            behind = r > r_p
            if not np.any(behind):
                continue

            # Ray height at r_p for each disk point at (r, phi, h):
            #   z_ray = hlamp + (r_p / r) * (h - hlamp)
            # This is the height of the line-of-sight from the lamp at the
            # radial position of the pillar.
            z_ray = np.where(behind,
                             self.hlamp + (r_p / r) * (h - self.hlamp),
                             np.inf)

            # Pillar height at (r_p, phi_disk).  At r = r_p the radial
            # Gaussian term is 1, so only the azimuthal profile matters.
            dphi = phi - phi_p
            dphi = np.mod(dphi + np.pi, 2*np.pi) - np.pi
            sig_phi = pillar['sigma_phi']
            gauss_phi = (np.exp(-0.5 * (dphi / sig_phi)**2)
                       + np.exp(-0.5 * ((dphi - 2*np.pi) / sig_phi)**2)
                       + np.exp(-0.5 * ((dphi + 2*np.pi) / sig_phi)**2))

            h_base_rp = np.interp(r_p, self.r, self.h_base)
            h_pillar_at_phi = h_base_rp + pillar['height'] * gauss_phi

            # Shadow: where pillar is taller than the ray → fully opaque.
            # Use a steep sigmoid for a thin penumbra at the boundary
            # (scale ~ 1% of pillar height) for numerical smoothness.
            delta = h_pillar_at_phi - z_ray  # >0 means ray is blocked
            penumbra_scale = 0.01 * pillar['height'] + 1e-10
            shadow_strength = np.where(behind,
                                       np.clip(delta / penumbra_scale, 0.0, 1.0),
                                       0.0)

            # For transparent pillars, shadow floors at f_trans instead of 0
            shadow_mask *= (1.0 - shadow_strength * (1.0 - self.f_trans))

        shadow_mask = np.clip(shadow_mask, 0.0, 1.0)

        return shadow_mask
    
    def get_temperature(self, r: np.ndarray, phi: np.ndarray, compute_shadows: bool = True) -> np.ndarray:
        """
        Get disk temperature at given (r, phi) positions, including pillars and shadows.
        
        Parameters:
        -----------
        r : array
            Radial positions (light days), can be 1D or 2D
        phi : array
            Azimuthal positions (radians), can be 1D or 2D
        compute_shadows : bool
            If True, compute shadows cast by pillars (default: True)
        
        Returns:
        --------
        T : array
            Effective temperature (K)
        """
        # Handle both 1D and 2D input arrays
        r = np.asarray(r)
        phi = np.asarray(phi)
        
        # If 1D arrays, create meshgrid
        if r.ndim == 1 and phi.ndim == 1:
            r_2d, phi_2d = np.meshgrid(r, phi, indexing='ij')
        elif r.ndim == 2 and phi.ndim == 2:
            r_2d = r
            phi_2d = phi
        else:
            # Try to broadcast
            r_2d, phi_2d = np.broadcast_arrays(r, phi)
        
        # Get heights (includes pillar bumps)
        h_2d = self.get_height(r_2d, phi_2d)

        # Base viscous temperature (axisymmetric): T_visc ∝ r^(-3/4)
        tv_2d = np.interp(r_2d.flatten(), self.r, self.tv_base).reshape(r_2d.shape)

        # ── Irradiation from lamp geometry: T_irr^4 ∝ cos(θ) / d² ──
        # Compute for ALL disk elements (flat disk + pillar slopes alike).
        # Lamp position
        x_2d = r_2d * np.cos(phi_2d)
        y_2d = r_2d * np.sin(phi_2d)
        z_2d = h_2d

        # Vector from surface to lamp
        dx = -x_2d
        dy = -y_2d
        dz = self.hlamp - z_2d
        d2 = dx**2 + dy**2 + dz**2  # distance squared
        d = np.sqrt(d2)

        # Unit lamp direction (from surface toward lamp)
        lamp_x = dx / (d + 1e-10)
        lamp_y = dy / (d + 1e-10)
        lamp_z = dz / (d + 1e-10)

        # Surface normal including azimuthal gradient (pillar walls)
        nx, ny, nz = self._surface_normal(r_2d, phi_2d, h_2d)

        # cos(θ_inc) = dot(lamp_dir, normal)
        cos_theta = lamp_x * nx + lamp_y * ny + lamp_z * nz
        cos_theta = np.clip(cos_theta, 0.0, 1.0)

        # Irradiation flux ∝ cos(θ) / d²
        F_irr = cos_theta / (d2 + 1e-10)

        # Normalise so that T_irr at (r0, flat disk) matches tx_base(r0).
        # Reference: flat disk at r0, h≈0 → d0² = r0² + hlamp², cos θ0 = hlamp/d0.
        d0 = np.sqrt(self.r0**2 + self.hlamp**2)
        cos_theta0 = self.hlamp / d0
        F_irr_ref = cos_theta0 / d0**2
        tx_r0 = np.interp(self.r0, self.r, self.tx_base)

        # T_irr^4 = tx_r0^4 × (F_irr / F_irr_ref)
        tx_2d = tx_r0 * (F_irr / (F_irr_ref + 1e-30))**0.25

        # Apply temp_factor boost on pillar illuminated faces
        # (multiplicative on top of the geometric cos_theta enhancement)
        if len(self.pillars) > 0:
            h_base_2d = np.interp(r_2d.flatten(), self.r, self.h_base).reshape(r_2d.shape)
            for pillar in self.pillars:
                if not pillar.get('modify_temp', False):
                    continue
                tf = pillar.get('temp_factor', 1.0)
                if tf == 1.0:
                    continue
                dr = r_2d - pillar['r']
                dphi = phi_2d - pillar['phi']
                dphi = np.mod(dphi + np.pi, 2*np.pi) - np.pi
                sig_phi = pillar['sigma_phi']
                gauss = (np.exp(-0.5 * (dr / pillar['sigma_r'])**2)
                       * (np.exp(-0.5 * (dphi / sig_phi)**2)
                        + np.exp(-0.5 * ((dphi - 2*np.pi) / sig_phi)**2)
                        + np.exp(-0.5 * ((dphi + 2*np.pi) / sig_phi)**2)))
                # Boost factor smoothly weighted by pillar profile
                boost = 1.0 + (tf - 1.0) * gauss
                tx_2d = tx_2d * boost

        # Apply shadow mask
        if compute_shadows and len(self.pillars) > 0:
            shadow_mask = self._compute_shadow_mask(r_2d, phi_2d, h_2d)
            tx_2d = tx_2d * shadow_mask

        # Effective temperature: T_e^4 = T_v^4 + T_irr^4
        T = (tv_2d**4 + tx_2d**4)**0.25

        return T
    
    def get_temperature_profile(self, compute_shadows: bool = True) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Compute azimuthally-averaged temperature profile T(r) and components.
        
        Parameters:
        -----------
        compute_shadows : bool
            If True, compute shadows cast by pillars (default: True)
        
        Returns:
        --------
        r_profile : array
            Radial positions (light days)
        T_with_pillars : array
            Total temperature profile with pillars (K), averaged over phi
        T_without_pillars : array
            Total temperature profile without pillars (K), averaged over phi
        T_visc : array
            Viscous temperature profile (K), same with/without pillars
        T_irrad_with : array
            Irradiation temperature profile with pillars (K), averaged over phi
        T_irrad_without : array
            Irradiation temperature profile without pillars (K), averaged over phi
        """
        # Create 2D grid
        r_2d, phi_2d = np.meshgrid(self.r, self.phi, indexing='ij')
        
        # Get heights
        h_2d = self.get_height(r_2d, phi_2d)
        
        # Base viscous temperature (same with/without pillars)
        tv_2d = np.interp(r_2d.flatten(), self.r, self.tv_base).reshape(r_2d.shape)
        T_visc = np.mean(tv_2d, axis=1)
        
        # Base irradiation temperature (axisymmetric)
        tx_2d_base = np.interp(r_2d.flatten(), self.r, self.tx_base).reshape(r_2d.shape)
        
        # Compute temperature with pillars
        T_2d_with = self.get_temperature(r_2d, phi_2d, compute_shadows=compute_shadows)
        
        # Compute irradiation temperature with pillars (accounting for shadows)
        tx_2d_with = tx_2d_base.copy()
        if compute_shadows and len(self.pillars) > 0:
            shadow_mask = self._compute_shadow_mask(r_2d, phi_2d, h_2d)
            tx_2d_with = tx_2d_with * shadow_mask
        
        # Compute temperature without pillars
        pillars_backup = self.pillars.copy()
        self.pillars = []
        T_2d_without = self.get_temperature(r_2d, phi_2d, compute_shadows=False)
        # Without pillars, no shadows, so T_irrad is just the base
        tx_2d_without = tx_2d_base.copy()
        self.pillars = pillars_backup
        
        # Average over phi
        T_with_pillars = np.mean(T_2d_with, axis=1)
        T_without_pillars = np.mean(T_2d_without, axis=1)
        T_irrad_with = np.mean(tx_2d_with, axis=1)
        T_irrad_without = np.mean(tx_2d_without, axis=1)
        
        return self.r, T_with_pillars, T_without_pillars, T_visc, T_irrad_with, T_irrad_without
    
    def planck_function(self, wavelength: float, temperature: np.ndarray) -> np.ndarray:
        """
        Compute Planck function B_nu(wavelength, T) in erg/cm^2/s/Hz/ster.
        
        Parameters:
        -----------
        wavelength : float
            Wavelength in Angstroms
        temperature : array
            Temperature in Kelvin
        
        Returns:
        --------
        B_nu : array
            Blackbody intensity (erg/cm^2/s/Hz/ster)
        """
        # Constants
        c1 = 2.0 * np.pi * H * C**2  # 2*pi*h*c^2
        c2 = H * C / K  # h*c/k
        
        # Convert wavelength to cm
        wave_cm = wavelength * ANGSTROM
        
        # Frequency
        nu = C / wave_cm
        
        # Dimensionless frequency x = h*nu/(k*T)
        x = c2 / (wave_cm * temperature)
        
        # Avoid overflow/underflow
        x = np.clip(x, 1e-10, 700)
        
        # Planck function: B_nu = (2*h*nu^3/c^2) / (exp(h*nu/(k*T)) - 1)
        # In terms of wavelength: B_nu = (2*pi*h*c^2/lambda^5) / (exp(h*c/(k*T*lambda)) - 1)
        expx = np.exp(x)
        B_nu = (2.0 * np.pi * H * C**2) / (wave_cm**5 * (expx - 1.0))
        
        # Handle very small x (Rayleigh-Jeans limit)
        mask = x < 1e-4
        if np.any(mask):
            B_nu[mask] = (2.0 * K * temperature[mask] * nu**2) / (C**2)
        
        return B_nu
    
    def compute_sed(self, wavelengths: np.ndarray) -> np.ndarray:
        """
        Compute spectral energy distribution (SED) in mJy.
        
        Parameters:
        -----------
        wavelengths : array
            Wavelengths in Angstroms (rest frame)
        
        Returns:
        --------
        flux : array
            Flux in mJy
        """
        # Create 2D grids
        r_2d, phi_2d = np.meshgrid(self.r, self.phi, indexing='ij')
        
        # Get height and temperature
        h_2d = self.get_height(r_2d, phi_2d)
        T_2d = self.get_temperature(r_2d, phi_2d)
        
        # Color temperature
        Tcol_2d = T_2d * self.fcol
        f4 = self.fcol**4
        
        # Initialize flux array
        flux = np.zeros_like(wavelengths)
        
        # Units conversion
        ld_to_cm = C * DAY  # light days to cm
        d_cm = self.d * ld_to_cm  # distance in cm
        units = FNU_TO_MJY / (d_cm**2)  # conversion factor
        
        # Loop over wavelengths with progress bar
        for iw, w in enumerate(tqdm(wavelengths, desc="Computing SED", unit="wavelength")):
            # Planck function at this wavelength
            B_nu = self.planck_function(w, Tcol_2d) / f4
            
            # Sum over disk surface
            sum_flux = 0.0
            
            for ir in range(self.nr - 1):
                r_now = self.r[ir]
                h_now = h_2d[ir, :]
                B_now = B_nu[ir, :]
                
                # Next radius
                r_next = self.r[ir + 1]
                h_next = h_2d[ir + 1, :]
                
                # Surface element
                dr = r_next - r_now
                dh = h_next - h_now
                ds = np.sqrt(dr**2 + dh**2)
                
                # Tilt angle
                sintilt = dh / (ds + 1e-10)
                costilt = dr / (ds + 1e-10)
                
                # Loop over azimuth
                for iphi, phi in enumerate(self.phi):
                    # Normal vector
                    px = -sintilt[iphi] * np.cos(phi)
                    pz = costilt[iphi]
                    
                    # Foreshortening factor (dot product with Earth direction)
                    dot = self.ex * px + self.ez * pz
                    
                    # Only count if facing Earth
                    if dot > 0:
                        # Area element: ds * r * dphi
                        da = ds[iphi] * r_now * self.dphi
                        # Flux contribution: B_nu * area * foreshortening
                        sum_flux += B_now[iphi] * da * dot
            
            # Convert to mJy
            # F_nu = (1/D^2) * Int B_nu dOmega
            # B_nu is in erg/cm^2/s/Hz/ster, area in (light days)^2
            # Convert area to cm^2
            area_cm2 = sum_flux * (ld_to_cm**2)
            # Flux in erg/cm^2/s/Hz
            fnu_erg = area_cm2 / (d_cm**2)
            # Convert to mJy
            flux[iw] = fnu_erg * FNU_TO_MJY
        
        return flux
    
    def compute_lag_spectrum(self, wavelengths: np.ndarray, 
                            ntau: int = 1000, taumax: float = 200.0,
                            parallel: bool = False, n_jobs: int = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Compute lag time spectrum (mean delay vs wavelength).
        
        Parameters:
        -----------
        wavelengths : array
            Wavelengths in Angstroms (rest frame)
        ntau : int
            Number of delay bins
        taumax : float
            Maximum delay (days)
        parallel : bool
            If True, use multiprocessing to parallelize over wavelengths
        n_jobs : int
            Number of parallel jobs (default: number of CPU cores)
        
        Returns:
        --------
        tau_mean : array
            Mean delay for each wavelength (days)
        tau_grid : array
            Delay grid (days)
        psi : array
            Delay distribution for each wavelength (ntau, nw)
        """
        # Delay grid
        tau_grid = np.linspace(0, taumax, ntau)
        dtau = taumax / (ntau - 1)
        
        # Initialize arrays
        nw = len(wavelengths)
        tau_mean = np.zeros(nw)
        psi = np.zeros((ntau, nw))
        
        # Create 2D grids
        r_2d, phi_2d = np.meshgrid(self.r, self.phi, indexing='ij')
        
        # Get height and temperature
        h_2d = self.get_height(r_2d, phi_2d)
        T_2d = self.get_temperature(r_2d, phi_2d)
        
        # Color temperature
        Tcol_2d = T_2d * self.fcol
        f4 = self.fcol**4
        
        # Wien temperature for normalization
        T0 = H * C / (K * wavelengths * ANGSTROM)
        
        if parallel and HAS_MP:
            # Parallel computation
            if n_jobs is None:
                n_jobs = cpu_count()
            
            # Prepare arguments for parallel processing
            args_list = [(iw, w, r_2d, phi_2d, h_2d, T_2d, Tcol_2d, f4, T0, 
                         tau_grid, dtau, ntau, taumax) 
                        for iw, w in enumerate(wavelengths)]
            
            # Process in parallel
            with Pool(n_jobs) as pool:
                results = list(tqdm(pool.imap(self._compute_lag_single_wavelength, args_list),
                                   total=nw, desc="Computing lag spectrum (parallel)", 
                                   unit="wavelength"))
            
            # Unpack results
            for iw, (tau_mean_iw, psi_w) in enumerate(results):
                tau_mean[iw] = tau_mean_iw
                psi[:, iw] = psi_w
        else:
            # Sequential computation
            # Loop over wavelengths with progress bar
            for iw, w in enumerate(tqdm(wavelengths, desc="Computing lag spectrum", unit="wavelength")):
                # Planck function
                B_nu = self.planck_function(w, Tcol_2d) / f4
                
                # Temperature derivative dB/dT
                # dB/dT = (k/w^2) * x^2 / (cosh(x) - 1) where x = h*nu/(k*T)
                x = T0[iw] / Tcol_2d
                x = np.clip(x, 1e-3, 100)
                dBdT = (K / (w * ANGSTROM)**2) * (x**2 / (np.cosh(x) - 1.0 + 1e-10))
                
                # Response function: dI/dT * dB/dT
                # For lamp heating: dT^4/dL = (1-a)*cos(theta)/(4*pi*sigma*dlamp^2)
                # dT/dL = dT/d(T^4) * d(T^4)/dL = (1/(4*T^3)) * d(T^4)/dL
                
                # Initialize delay distribution
                psi_w = np.zeros(ntau)
                
                # Sum over disk surface
                for ir in range(self.nr - 1):
                    r_now = self.r[ir]
                    h_now = h_2d[ir, :]
                    T_now = T_2d[ir, :]
                    Tcol_now = Tcol_2d[ir, :]
                    B_now = B_nu[ir, :]
                    dBdT_now = dBdT[ir, :]
                    
                    # Next radius
                    r_next = self.r[ir + 1]
                    h_next = h_2d[ir + 1, :]
                    
                    # Surface element
                    dr = r_next - r_now
                    dh = h_next - h_now
                    ds = np.sqrt(dr**2 + dh**2)
                    da = ds * r_now * self.dphi
                    
                    # Tilt angle
                    sintilt = dh / (ds + 1e-10)
                    costilt = dr / (ds + 1e-10)
                    
                    # Loop over azimuth
                    for iphi, phi in enumerate(self.phi):
                        # Skip if temperature too low
                        if T_now[iphi] < T0[iw] / 200:
                            continue
                        
                        # Normal vector
                        px = -sintilt[iphi] * np.cos(phi)
                        pz = costilt[iphi]
                        
                        # Foreshortening factor
                        dot = self.ex * px + self.ez * pz
                        
                        if dot <= 0:
                            continue
                        
                        # Time delay
                        # Distance from lamp
                        dx = r_now * np.cos(phi)
                        dz = h_now[iphi] - self.hlamp
                        dlamp = np.sqrt(dx**2 + dz**2)
                        
                        # Delay = light travel time from lamp to disk element to Earth
                        # tau = dlamp - (r*cos(phi)*sin(i) + h*cos(i))
                        tau_delay = dlamp - (dx * self.sini + h_now[iphi] * self.cosi)
                        
                        # Response weight
                        # dI/dT * dB/dT * dOmega
                        # The response is proportional to dB/dT * B_nu * area
                        # For lamp heating, the response is weighted by dB/dT
                        weight = dBdT_now[iphi] * B_now[iphi] * da[iphi] * dot
                        
                        # Add to delay distribution (Gaussian blur to avoid delta function)
                        sigtau = dtau
                        if tau_delay >= 0:
                            # Find bin range for Gaussian spread
                            nsigma = 3
                            ibin_min = max(0, int((tau_delay - nsigma * sigtau) / dtau))
                            ibin_max = min(ntau - 1, int((tau_delay + nsigma * sigtau) / dtau))
                            
                            # Normalize Gaussian
                            gauss_norm = 0.0
                            for ibin in range(ibin_min, ibin_max + 1):
                                tau_rel = tau_delay - tau_grid[ibin]
                                gauss = np.exp(-0.5 * (tau_rel / sigtau)**2)
                                gauss_norm += gauss
                            
                            # Add weighted contribution
                            if gauss_norm > 0:
                                for ibin in range(ibin_min, ibin_max + 1):
                                    tau_rel = tau_delay - tau_grid[ibin]
                                    gauss = np.exp(-0.5 * (tau_rel / sigtau)**2) / gauss_norm
                                    psi_w[ibin] += weight * gauss
                
                # Normalize delay distribution
                norm = np.sum(psi_w) * dtau
                if norm > 0:
                    psi_w = psi_w / norm
                    # Mean delay
                    tau_mean[iw] = np.sum(psi_w * tau_grid) * dtau
                else:
                    tau_mean[iw] = 0.0
                
                psi[:, iw] = psi_w
        
        return tau_mean, tau_grid, psi
    
    def _compute_lag_single_wavelength(self, args):
        """
        Helper method for parallel computation of lag at a single wavelength.
        This is a static-like method that takes all needed parameters.
        """
        (iw, w, r_2d, phi_2d, h_2d, T_2d, Tcol_2d, f4, T0, 
         tau_grid, dtau, ntau, taumax) = args
        
        # Planck function
        B_nu = self.planck_function(w, Tcol_2d) / f4
        
        # Temperature derivative dB/dT
        x = T0[iw] / Tcol_2d
        x = np.clip(x, 1e-3, 100)
        dBdT = (K / (w * ANGSTROM)**2) * (x**2 / (np.cosh(x) - 1.0 + 1e-10))
        
        # Initialize delay distribution
        psi_w = np.zeros(ntau)
        
        # Sum over disk surface
        for ir in range(self.nr - 1):
            r_now = self.r[ir]
            h_now = h_2d[ir, :]
            T_now = T_2d[ir, :]
            Tcol_now = Tcol_2d[ir, :]
            B_now = B_nu[ir, :]
            dBdT_now = dBdT[ir, :]
            
            # Next radius
            r_next = self.r[ir + 1]
            h_next = h_2d[ir + 1, :]
            
            # Surface element
            dr = r_next - r_now
            dh = h_next - h_now
            ds = np.sqrt(dr**2 + dh**2)
            da = ds * r_now * self.dphi
            
            # Tilt angle
            sintilt = dh / (ds + 1e-10)
            costilt = dr / (ds + 1e-10)
            
            # Loop over azimuth
            for iphi, phi in enumerate(self.phi):
                # Skip if temperature too low
                if T_now[iphi] < T0[iw] / 200:
                    continue
                
                # Normal vector
                px = -sintilt[iphi] * np.cos(phi)
                pz = costilt[iphi]
                
                # Foreshortening factor
                dot = self.ex * px + self.ez * pz
                
                if dot <= 0:
                    continue
                
                # Time delay
                dx = r_now * np.cos(phi)
                dz = h_now[iphi] - self.hlamp
                dlamp = np.sqrt(dx**2 + dz**2)
                
                tau_delay = dlamp - (dx * self.sini + h_now[iphi] * self.cosi)
                
                # Response weight
                weight = dBdT_now[iphi] * B_now[iphi] * da[iphi] * dot
                
                # Add to delay distribution
                sigtau = dtau
                if tau_delay >= 0:
                    nsigma = 3
                    ibin_min = max(0, int((tau_delay - nsigma * sigtau) / dtau))
                    ibin_max = min(ntau - 1, int((tau_delay + nsigma * sigtau) / dtau))
                    
                    gauss_norm = 0.0
                    for ibin in range(ibin_min, ibin_max + 1):
                        tau_rel = tau_delay - tau_grid[ibin]
                        gauss = np.exp(-0.5 * (tau_rel / sigtau)**2)
                        gauss_norm += gauss
                    
                    if gauss_norm > 0:
                        for ibin in range(ibin_min, ibin_max + 1):
                            tau_rel = tau_delay - tau_grid[ibin]
                            gauss = np.exp(-0.5 * (tau_rel / sigtau)**2) / gauss_norm
                            psi_w[ibin] += weight * gauss
        
        # Normalize delay distribution
        norm = np.sum(psi_w) * dtau
        if norm > 0:
            psi_w = psi_w / norm
            tau_mean_iw = np.sum(psi_w * tau_grid) * dtau
        else:
            tau_mean_iw = 0.0
        
        return tau_mean_iw, psi_w
    
    def compute_sed_no_pillars(self, wavelengths: np.ndarray) -> np.ndarray:
        """
        Compute SED without pillars (temporarily disables pillars).
        """
        # Save current pillars
        saved_pillars = self.pillars.copy()
        # Temporarily remove pillars
        self.pillars = []
        # Compute SED
        flux = self.compute_sed(wavelengths)
        # Restore pillars
        self.pillars = saved_pillars
        return flux
    
    def compute_lag_spectrum_no_pillars(self, wavelengths: np.ndarray, 
                                       ntau: int = 1000, taumax: float = 200.0,
                                       parallel: bool = False, n_jobs: int = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Compute lag spectrum without pillars (temporarily disables pillars).
        """
        # Save current pillars
        saved_pillars = self.pillars.copy()
        # Temporarily remove pillars
        self.pillars = []
        # Compute lag spectrum
        result = self.compute_lag_spectrum(wavelengths, ntau, taumax, parallel, n_jobs)
        # Restore pillars
        self.pillars = saved_pillars
        return result
    
    def compute_log_ionizing_flux(self, r, phi):
        """Compute log ionizing photon flux.

        Two methods available (set via ``self.flux_method``):

        * ``'scale_temp'`` (legacy) -- rescale T_e / T_ref linearly to
          [17, 21.5].
        * ``'direct'`` -- geometric lamp flux + viscous blackbody Wien
          tail above 1 Ryd, clipped to the Cloudy grid [17, 21].

        Parameters
        ----------
        r, phi : array-like
            Radial (light-days) and azimuthal (radians) grids, same shape.

        Returns
        -------
        log_phi : ndarray
            log10 ionizing flux.
        """
        if self.flux_method == 'direct':
            return self._compute_direct_ionizing_flux(r, phi)
        else:
            return self._compute_scale_temp_ionizing_flux(r, phi)

    # ------------------------------------------------------------------
    @staticmethod
    def _surface_normal(r_2d, phi_2d, h_2d):
        """Compute outward surface normal for z = h(r, phi).

        Uses both radial (dh/dr) and azimuthal (dh/dphi) gradients so that
        pillar walls are correctly oriented.

        Tangent vectors on the surface:
            e_r   = (cos phi, sin phi, dh/dr)
            e_phi = (-r sin phi, r cos phi, dh/dphi)
        Normal = e_r x e_phi (unnormalised):
            n_x = sin(phi) dh/dphi - r cos(phi) dh/dr
            n_y = -r sin(phi) dh/dr - cos(phi) dh/dphi
            n_z = r

        Returns unit normal (nx, ny, nz).
        """
        # Radial gradient dh/dr (axis 0 = r)
        dr_grid = np.gradient(r_2d, axis=0)
        dh_r_grid = np.gradient(h_2d, axis=0)
        dh_dr = np.zeros_like(h_2d)
        ok_r = np.abs(dr_grid) > 1e-10
        dh_dr[ok_r] = dh_r_grid[ok_r] / dr_grid[ok_r]

        # Azimuthal gradient dh/dphi (axis 1 = phi)
        dphi_grid = np.gradient(phi_2d, axis=1)
        dh_p_grid = np.gradient(h_2d, axis=1)
        dh_dphi = np.zeros_like(h_2d)
        ok_p = np.abs(dphi_grid) > 1e-10
        dh_dphi[ok_p] = dh_p_grid[ok_p] / dphi_grid[ok_p]

        # Cross product
        nx_raw = np.sin(phi_2d) * dh_dphi - r_2d * np.cos(phi_2d) * dh_dr
        ny_raw = -r_2d * np.sin(phi_2d) * dh_dr - np.cos(phi_2d) * dh_dphi
        nz_raw = r_2d

        mag = np.sqrt(nx_raw**2 + ny_raw**2 + nz_raw**2) + 1e-30
        return nx_raw / mag, ny_raw / mag, nz_raw / mag

    # ------------------------------------------------------------------
    def _compute_scale_temp_ionizing_flux(self, r, phi):
        """Legacy method: linear rescaling of T_e/T_ref to [17, 21.5]."""
        T_2d = self.get_temperature(r, phi, compute_shadows=True)
        T_ref = self.tv1 * self.tirrad_tvisc_ratio
        T_ratio = np.clip(T_2d / T_ref, 0.1, 10.0)
        T_ratio_normalized = (T_ratio - 0.1) / (10.0 - 0.1)
        phi_min, phi_max = 17.0, 21.5
        return phi_min + T_ratio_normalized * (phi_max - phi_min)

    # ------------------------------------------------------------------
    def _compute_direct_ionizing_flux(self, r, phi):
        """Direct geometric lamp flux + viscous blackbody ionizing flux.

        Phi_total = S(r,phi) * Q * cos(theta) / (4 pi d^2) + Phi_visc(r)

        Q is calibrated so that Phi(r_in) = 10^log_phi_inner on the
        illuminated flat disk.  The result is clipped to [17, 21].
        """
        r = np.asarray(r, dtype=float)
        phi = np.asarray(phi, dtype=float)
        if r.ndim == 1 and phi.ndim == 1:
            r_2d, phi_2d = np.meshgrid(r, phi, indexing='ij')
        elif r.ndim == 2 and phi.ndim == 2:
            r_2d, phi_2d = r, phi
        else:
            r_2d, phi_2d = np.broadcast_arrays(r, phi)

        h_2d = self.get_height(r_2d, phi_2d)

        # ---- Lamp geometry ----
        x_2d = r_2d * np.cos(phi_2d)
        y_2d = r_2d * np.sin(phi_2d)
        z_2d = h_2d

        dx = -x_2d
        dy = -y_2d
        dz = self.hlamp - z_2d
        d2 = dx**2 + dy**2 + dz**2
        d = np.sqrt(d2)

        lamp_x = dx / (d + 1e-30)
        lamp_y = dy / (d + 1e-30)
        lamp_z = dz / (d + 1e-30)

        # Surface normal including azimuthal gradient (pillar walls)
        nx, ny, nz = self._surface_normal(r_2d, phi_2d, h_2d)

        cos_theta_raw = lamp_x * nx + lamp_y * ny + lamp_z * nz
        if self.f_trans > 0:
            # For transparent pillars, the pillar wall still receives transmitted
            # lamp flux. Floor at f_trans * cos_theta_flat, where cos_theta_flat
            # is what a flat disk element at that radius would see from the lamp
            # (using base disk height, not pillar-modified height).
            h_base = np.interp(r_2d.flatten(), self.r, self.h_base).reshape(r_2d.shape)
            d_flat = np.sqrt(r_2d**2 + (self.hlamp - h_base)**2)
            cos_theta_flat = (self.hlamp - h_base) / (d_flat + 1e-30)
            cos_theta_flat = np.clip(cos_theta_flat, 0.0, 1.0)
            cos_theta = np.maximum(cos_theta_raw, self.f_trans * cos_theta_flat)
        else:
            cos_theta = np.clip(cos_theta_raw, 0.0, 1.0)

        # ---- Calibrate Q at inner edge (flat disk) ----
        d_in = np.sqrt(self.rin**2 + self.hlamp**2)
        cos_theta_in = self.hlamp / d_in
        Q = 10.0**self.log_phi_inner * 4.0 * np.pi * d_in**2 / cos_theta_in

        # ---- Lamp flux ----
        Phi_lamp = Q * cos_theta / (4.0 * np.pi * d2)

        # Shadow mask
        if len(self.pillars) > 0:
            shadow_mask = self._compute_shadow_mask(r_2d, phi_2d, h_2d)
            Phi_lamp = Phi_lamp * shadow_mask

        # ---- Viscous blackbody ionizing flux (Wien tail above 1 Ryd) ----
        # Phi_visc = integral of B_nu(T) / (h*nu) dnu from nu_0 to inf
        # Wien-tail leading order:
        #   Phi_visc ≈ (2*pi*nu_0^2/c^2) * (k*T/h) * exp(-h*nu_0/(k*T))
        h_planck = 6.626e-27   # erg s
        k_B = 1.381e-16        # erg/K
        c_cgs = 2.998e10       # cm/s
        nu_0 = 13.6 * 1.602e-12 / h_planck  # 1 Ryd frequency (Hz)

        # Viscous temperature only (no irradiation)
        T_visc = np.interp(r_2d.flatten(), self.r, self.tv_base).reshape(r_2d.shape)

        A = 2.0 * np.pi * nu_0**2 / c_cgs**2 * (k_B / h_planck)
        x = h_planck * nu_0 / (k_B * np.maximum(T_visc, 100.0))
        # For numerical safety clip x; when x > 500 the exp is negligible
        x = np.clip(x, 0.0, 500.0)
        Phi_visc = A * T_visc * np.exp(-x)

        # ---- Total ----
        Phi_total = Phi_lamp + Phi_visc

        # Floor to avoid log10(0)
        Phi_total = np.maximum(Phi_total, 1e-30)
        log_phi = np.log10(Phi_total)

        # Clip upper bound to Cloudy grid max
        log_phi = np.clip(log_phi, None, 21.0)
        if not self.no_fluxfloor:
            # Clip lower bound to Cloudy grid min (legacy: shadows get min emissivity)
            log_phi = np.clip(log_phi, 17.0, None)
        return log_phi

    def plot_3d_geometry(self, n_phi_plot: int = 180, n_r_plot: int = 200,
                         show_light_rays: bool = True, show_shadows: bool = True,
                         show_pillar_markers: bool = True,
                         show_colorbar: bool = True,
                         color_by: str = 'shadow', filename: str = 'pillar_disk_3d.png',
                         title: str = None):
        """
        Create a 3D visualization of the disk geometry with pillars and light rays.

        The viewing angle is set to match the actual observer inclination:
        - Observer is at direction (sin(i), 0, cos(i)) in the x-z plane
        - Camera elevation = 90° - inclination
        - Camera azimuth = 0° (looking from +x direction)

        Parameters:
        -----------
        n_phi_plot : int
            Number of azimuthal points for plotting (increased default for better pillar resolution)
        n_r_plot : int
            Number of radial points for plotting (increased default for better pillar resolution)
        show_light_rays : bool
            Whether to show light rays from lamp to disk
        show_shadows : bool
            Whether to show shadow regions (illuminated regions in cyan, shadows in dark gray)
        color_by : str
            How to color the disk surface:
            - 'shadow': cyan for illuminated, gray for shadow (default)
            - 'flux' or 'temperature': color by ionizing flux (log Φ)
        filename : str
            Output filename for the plot
        """
        if not HAS_3D:
            print("Warning: 3D plotting not available. Install matplotlib with 3D support.")
            return
        
        fig = plt.figure(figsize=(14, 10))
        ax = fig.add_subplot(111, projection='3d')
        
        # Create plotting grid with higher resolution to capture pillars
        # Use finer sampling near pillar locations
        r_plot = np.logspace(np.log10(self.rin), np.log10(self.rout), n_r_plot)
        
        # If there are pillars, add extra points around them for better resolution
        if len(self.pillars) > 0:
            # Collect all pillar radii
            pillar_radii = [p['r'] for p in self.pillars]
            # Add extra points around each pillar
            extra_r_points = []
            for r_p in pillar_radii:
                # Add points around pillar (within 5*sigma_r)
                sigma_max = max([p.get('sigma_r', 3.0) for p in self.pillars if p['r'] == r_p])
                r_min = max(self.rin, r_p - 5 * sigma_max)
                r_max = min(self.rout, r_p + 5 * sigma_max)
                # Add finer grid around pillar
                n_extra = 20
                extra_r = np.linspace(r_min, r_max, n_extra)
                extra_r_points.extend(extra_r)
            
            # Combine and sort, remove duplicates
            if extra_r_points:
                all_r = np.concatenate([r_plot, extra_r_points])
                all_r = np.unique(np.sort(all_r))
                # Limit to reasonable number
                if len(all_r) > n_r_plot * 4:
                    # Downsample while keeping pillar regions. Cap raised
                    # from x2 to x4 because narrow pillars (sigma_r ~ 0.1
                    # ld) need finer adaptive refinement than the legacy
                    # cap allowed.
                    indices = np.linspace(0, len(all_r)-1, n_r_plot * 4, dtype=int)
                    all_r = all_r[indices]
                r_plot = all_r
        
        phi_plot = np.linspace(0, 2*np.pi, n_phi_plot, endpoint=False)
        r_2d, phi_2d = np.meshgrid(r_plot, phi_plot, indexing='ij')

        # Get height
        h_2d = self.get_height(r_2d, phi_2d)

        # Close the surface by appending the first column (φ=0) as the last column (φ=2π)
        phi_2d = np.concatenate([phi_2d, phi_2d[:, :1] + 2*np.pi], axis=1)
        r_2d = np.concatenate([r_2d, r_2d[:, :1]], axis=1)
        h_2d = np.concatenate([h_2d, h_2d[:, :1]], axis=1)

        # Convert to Cartesian coordinates
        x_2d = r_2d * np.cos(phi_2d)
        y_2d = r_2d * np.sin(phi_2d)
        z_2d = h_2d
        
        # Plot disk surface with coloring based on color_by parameter
        if color_by == 'temperature' or color_by == 'flux':
            log_phi_2d = self.compute_log_ionizing_flux(r_2d, phi_2d)
            finite = log_phi_2d[np.isfinite(log_phi_2d) & (log_phi_2d > 0)]
            phi_min = min(16.0, np.floor(np.min(finite) * 2) / 2)
            phi_max = np.ceil(np.max(finite) * 2) / 2

            # For pillar regions, apply view-dependent shading
            h_base = self.h1 * (r_2d / self.r0) ** self.beta
            pillar_mask = (h_2d - h_base) > 0.001

            # Which face does observer see?
            # face_dot_obs > 0: shadow side (near), face_dot_obs < 0: illuminated side (far)
            face_dot_obs = np.cos(phi_2d) * self.sini
            blend = 0.5 * (1.0 + np.tanh(face_dot_obs / 0.3))

            # Get the shadow flux from the disk shadow region (use minimum in shadowed areas)
            shadow_mask = self._compute_shadow_mask(r_2d, phi_2d, h_2d)
            # Where shadow_mask < 0.5, the disk is shadowed - find typical flux there
            shadowed_flux = np.where(shadow_mask < 0.5, log_phi_2d, phi_max)
            shadow_level = np.min(shadowed_flux)  # Darkest shadow on disk

            # Pillar shadow side should match disk shadow
            log_phi_view = np.where(pillar_mask,
                                    (1.0 - blend) * log_phi_2d + blend * shadow_level,
                                    log_phi_2d)

            # Normalize for colormap
            log_phi_norm = (log_phi_view - phi_min) / (phi_max - phi_min)
            log_phi_norm = np.clip(log_phi_norm, 0, 1)

            # Use 'inferno' colormap
            from matplotlib import cm
            cmap = cm.inferno
            plot_colors = cmap(log_phi_norm)
            plot_colors[:, :, 3] = 1.0  # Fully opaque

            # Hard back-face mask: any disc cell whose outward (concave-
            # side) normal points away from the observer is showing its
            # convex underside, which is never illuminated by the lamp.
            # Override its colour to pure black.
            #
            # We compute the normal in the *rendered* (display) frame, not
            # the physical frame: matplotlib auto-stretches the Z axis to
            # fill the same screen height as X/Y (here Z spans ~0.2 ld vs
            # X/Y spans 40 ld -> ~200x stretch), so the bowl that is
            # geometrically nearly flat (dh/dr ~ 0.025 at the rim) looks
            # like a deep cup on screen. A back-face mask computed in
            # physical coords would not fire on the bowl at all. We rescale
            # dh/dr and dh/dphi by the visual Z-stretch factor so the
            # mask matches what the viewer actually sees.
            r_vec = r_2d[:, 0]
            phi_vec = phi_2d[0, :]
            dh_dr = np.gradient(h_2d, r_vec, axis=0)
            dh_dphi = np.gradient(h_2d, phi_vec, axis=1)
            # z_stretch = full geometric ratio (X-span / Z-span ~ 200 for
            # the fiducial bowl) so the back-face mask matches what
            # matplotlib actually renders on screen (Z auto-stretched to
            # fill the same display height as X/Y). Smaller values let
            # too much of the visible underside escape the mask.
            x_range = 2.0 * float(self.rout)
            z_range = max(0.1, float(np.max(h_2d)) * 1.2)
            z_stretch = x_range / z_range
            dh_dr_vis = dh_dr * z_stretch
            dh_dphi_vis = dh_dphi * z_stretch
            cos_phi = np.cos(phi_2d)
            sin_phi = np.sin(phi_2d)
            dh_dx = dh_dr_vis * cos_phi - (dh_dphi_vis / r_2d) * sin_phi
            dh_dy = dh_dr_vis * sin_phi + (dh_dphi_vis / r_2d) * cos_phi
            n_x = -dh_dx
            n_y = -dh_dy
            n_z = np.ones_like(h_2d)
            n_norm = np.sqrt(n_x**2 + n_y**2 + n_z**2)
            n_x /= n_norm; n_y /= n_norm; n_z /= n_norm
            n_dot_o = n_x * float(self.sini) + n_z * float(self.cosi)
            back_mask = n_dot_o < 0
            plot_colors[back_mask, :3] = 0.0

            # Plot disk surface with ionizing flux coloring.
            # rstride/cstride=1 so every (r, phi) cell is rendered;
            # otherwise narrow pillars (sigma_phi ~ 0.05 rad) sub-sample
            # below the surface-mesh resolution and look like jagged
            # squares.
            surf = ax.plot_surface(x_2d, y_2d, z_2d, facecolors=plot_colors,
                           shade=False, antialiased=True,
                           rstride=1, cstride=1)

            # Add colorbar for ionizing flux (optional; the paper uses a
            # standalone horizontal cbar PNG generated separately so
            # each 3D panel can render without the right-side margin).
            if show_colorbar:
                from matplotlib.colors import Normalize
                norm = Normalize(vmin=phi_min, vmax=phi_max)
                mappable = cm.ScalarMappable(cmap=cmap, norm=norm)
                mappable.set_array(log_phi_view)
                cbar = plt.colorbar(mappable, ax=ax, shrink=0.5, pad=0.1,
                                    anchor=(0.5, 0.0))
                cbar.set_label(r'$\log\Phi~\rm [photons~cm^{-2}~s^{-1}]$')

        elif show_shadows and len(self.pillars) > 0:
            # Compute shadow mask
            shadow_mask = self._compute_shadow_mask(r_2d, phi_2d, h_2d)

            # Create colormap: illuminated regions are cyan, shadowed regions are dark gray
            # Use shadow_mask to create colors
            colors = np.zeros((r_2d.shape[0], r_2d.shape[1], 4))  # RGBA
            # Illuminated: cyan with transparency
            colors[:, :, 0] = 0.0  # R
            colors[:, :, 1] = 1.0  # G
            colors[:, :, 2] = 1.0  # B
            colors[:, :, 3] = 0.6 * shadow_mask + 0.2 * (1.0 - shadow_mask)  # Alpha: brighter when illuminated

            # Shadowed regions: dark gray
            shadow_colors = np.zeros_like(colors)
            shadow_colors[:, :, 0] = 0.3  # R
            shadow_colors[:, :, 1] = 0.3  # G
            shadow_colors[:, :, 2] = 0.3  # B
            shadow_colors[:, :, 3] = 0.3  # Alpha

            # Blend colors based on shadow mask
            plot_colors = colors * shadow_mask[:, :, np.newaxis] + shadow_colors * (1.0 - shadow_mask[:, :, np.newaxis])

            # Plot disk surface with shadow coloring
            ax.plot_surface(x_2d, y_2d, z_2d, facecolors=plot_colors,
                           edgecolor='none', shade=True, antialiased=True,
                           linewidth=0, rstride=1, cstride=1)
        else:
            # Plot disk surface without shadow coloring
            ax.plot_surface(x_2d, y_2d, z_2d, alpha=0.6, color='cyan',
                           edgecolor='none', shade=True, antialiased=True,
                           linewidth=0, rstride=1, cstride=1)
        
        # Plot pillars
        for i, pillar in enumerate(self.pillars):
            r_p = float(pillar['r'])
            phi_p = float(pillar['phi'])
            # Get height at pillar center - get_height returns 2D array, extract scalar
            h_p_arr = self.get_height(np.array([r_p]), np.array([phi_p]))
            # Extract scalar from 2D array
            if isinstance(h_p_arr, np.ndarray):
                h_p = float(np.asarray(h_p_arr).flat[0])
            else:
                h_p = float(h_p_arr)

            x_p = float(r_p * np.cos(phi_p))
            y_p = float(r_p * np.sin(phi_p))
            z_p = float(h_p)

            # Draw pillar marker
            if show_pillar_markers:
                ax.scatter([x_p], [y_p], [z_p], color='red', s=50,
                          marker='o', edgecolors='darkred', linewidths=1, zorder=5)

            # Draw light ray from lamp to pillar (only for first few if many)
            if show_light_rays and i < 5:
                ax.plot([0.0, x_p], [0.0, y_p], [float(self.hlamp), z_p],
                       'r--', linewidth=1, alpha=0.3)

        # Plot lamp post as vertical line with star on top
        ax.plot([0, 0], [0, 0], [0, self.hlamp], color='gold', linewidth=3)
        ax.plot([0], [0], [self.hlamp], color='yellow', marker='*', markersize=18,
                markeredgecolor='white', markeredgewidth=1, linestyle='none')

        # Draw light ray to Earth (viewing direction)
        # Earth is at infinity in direction (sini, 0, cosi)
        # Draw a long line in the viewing direction
        # view_length = float(self.rout * 2)
        # ax.plot([0.0, float(self.sini * view_length)], 
        #        [0.0, 0.0], 
        #        [float(self.hlamp), float(self.hlamp + self.cosi * view_length)],
        #        'g--', linewidth=2, alpha=0.7, label='Viewing Direction')
        
        # Set labels and title. Disable matplotlib's auto rotate_label on
        # the z axis: by default it tries to align the label with the
        # vertical axis as the view rotates and ends up rendering the
        # text upside down ("[pl] Z" instead of "Z [ld]").
        # X and Y labels: keep matplotlib's auto rotation so they tilt to
        # follow the axis direction (matches the tick labels and the
        # natural 3D reading orientation). Only the Z label needs the
        # rotate-disable hack, because auto-rotation makes it upside-down
        # at the current viewing angle (renders "Z [light days]" as
        # "[syad thgil] Z").
        ax.set_xlabel(r'$X~\rm [light~days]$', labelpad=10)
        ax.set_ylabel(r'$Y~\rm [light~days]$', labelpad=10)
        ax.zaxis.set_rotate_label(False)
        ax.set_zlabel(r'$Z~\rm [light~days]$', labelpad=20, rotation=90)

        # Turn off grid and panes for cleaner look
        ax.grid(False)
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False
        # Hide pane edges
        ax.xaxis.pane.set_edgecolor('none')
        ax.yaxis.pane.set_edgecolor('none')
        ax.zaxis.pane.set_edgecolor('none')
        # Hide axis lines
        ax.xaxis.line.set_color('none')
        ax.yaxis.line.set_color('none')
        ax.zaxis.line.set_color('none')

        # Set equal aspect ratio
        max_range = self.rout
        ax.set_xlim([-max_range, max_range])
        ax.set_ylim([-max_range, max_range])
        # Set zlim based on actual height range
        z_max = max(0.1, float(np.max(z_2d)) * 1.2)
        ax.set_zlim([0, z_max])

        # Set viewing angle to match actual observer inclination
        # Observer is at direction (sini, 0, cosi) - always in x-z plane
        # For matplotlib: elev = 90 - inclination, azim = 0
        inclination_deg = np.degrees(np.arccos(self.cosi))
        elev = 90 - inclination_deg  # elevation from horizontal
        azim = 0  # observer always in +x direction
        ax.view_init(elev=elev, azim=azim)

        # Add title if provided
        if title is not None:
            ax.set_title(title, fontsize=24, pad=20)

        plt.tight_layout()
        plt.savefig(filename, dpi=200, bbox_inches='tight')
        print(f"  -> {filename} (i={inclination_deg:.1f}°, elev={elev:.1f}°)")
        plt.close()


def main(config_file='config.yaml'):
    """
    Example usage: Create a disk with pillars and compute SED and lag spectrum.
    
    Parameters:
    -----------
    config_file : str
        Path to configuration YAML file (default: 'config.yaml')
    """
    # Load configuration from file
    config = load_config(config_file)
    
    if config is None:
        # Use default parameters if config file not available
        print("Using default parameters (config file not found or YAML not available)")
        disk_params = {
            'rin': 0.1, 'rout': 100.0, 'nr': 200, 'nphi': 180,
            'h1': 0.01, 'r0': 10.0, 'beta': 1.0,
            'tv1': 1e4, 'alpha': 0.75, 'hlamp': 5.0, 'tx1': 8e3,
            'dmpc': 100.0, 'cosi': 0.5, 'fcol': 1.0, 'redshift': 0.0
        }
        pillars_config = [{
            'r_pillar': 20.0, 'phi_pillar': 0.0, 'height': 0.05,
            'sigma_r': 3.0, 'sigma_phi': 0.3, 'modify_height': True,
            'modify_temp': True, 'temp_factor': 1.5
        }]
        comp_params = {
            'wlog_min': 2.0, 'wlog_max': 5.5, 'nwavelengths': 100,
            'ntau': 20, 'taumax': 150.0, 'use_parallel': True, 'n_jobs': None
        }
        plot_params = {
            'save_results': True, 'save_3d': True,
            'results_filename': 'pillar_disk_results.png',
            'plot3d_filename': 'pillar_disk_3d.png', 'dpi': 150
        }
    else:
        # Extract parameters from config
        disk_params = config.get('disk', {}).copy()
        temp_params = config.get('temperature', {}).copy()
        lamp_params = config.get('lamp', {}).copy()
        obs_params = config.get('observation', {}).copy()
        
        # Convert string values to appropriate types
        # Float parameters
        float_params = ['rin', 'rout', 'h1', 'r0', 'beta', 'tv1', 'alpha',
                       'tx1', 'fcol', 'hlamp', 'dmpc', 'cosi', 'inclination', 'redshift',
                       'M_BH', 'r_isco_rg', 'tirrad_tvisc_ratio', 'log_phi_inner', 'f_trans']
        # Integer parameters
        int_params = ['nr', 'nphi']

        def convert_value(key, value):
            """Convert YAML value to appropriate Python type."""
            if value is None:
                return None
            # Handle string representations of numbers (including scientific notation)
            if isinstance(value, str):
                # Try to convert string to number
                try:
                    if key in int_params:
                        return int(float(value))  # Handle scientific notation like "1e4"
                    elif key in float_params:
                        return float(value)
                except (ValueError, TypeError):
                    # If conversion fails, check for boolean strings
                    if value.lower() in ('true', '1', 'yes'):
                        return True
                    elif value.lower() in ('false', '0', 'no'):
                        return False
                    elif value.lower() in ('null', 'none'):
                        return None
            # Handle numeric types directly
            if key in int_params:
                return int(float(value))  # Handle scientific notation
            elif key in float_params:
                return float(value)
            return value

        # Convert all parameters
        for params_dict in [disk_params, temp_params, lamp_params, obs_params]:
            for key, value in params_dict.items():
                params_dict[key] = convert_value(key, value)

        # Convert inclination (degrees) to cosi if provided
        if 'inclination' in obs_params and obs_params['inclination'] is not None:
            import numpy as np
            inclination_rad = np.radians(obs_params['inclination'])
            obs_params['cosi'] = np.cos(inclination_rad)
            del obs_params['inclination']

        # Compute rin from M_BH and r_isco_rg if rin='auto'
        resolve_rin(disk_params)

        # Merge all disk parameters
        disk_params = {**disk_params, **temp_params, **lamp_params, **obs_params}
        
        # Helper function to parse mathematical expressions (e.g., 'pi', 'pi/2')
        def parse_math_expr(value):
            """Parse mathematical expressions like 'pi', 'pi/2', etc."""
            if isinstance(value, (int, float)):
                return float(value)
            if isinstance(value, str):
                # Replace common mathematical constants
                expr = value.lower().strip()
                expr = expr.replace('pi', str(np.pi))
                expr = expr.replace('π', str(np.pi))
                # Evaluate simple expressions
                try:
                    # Use eval with safe namespace (only math operations)
                    safe_dict = {
                        'pi': np.pi,
                        'π': np.pi,
                        'e': np.e,
                        'sin': np.sin,
                        'cos': np.cos,
                        'tan': np.tan,
                        'sqrt': np.sqrt,
                        'exp': np.exp,
                        'log': np.log,
                    }
                    result = eval(expr, {"__builtins__": {}}, safe_dict)
                    return float(result)
                except:
                    # If evaluation fails, try direct float conversion
                    return float(value)
            return float(value)
        
        # Process pillars config
        # Support both formats:
        # 1. List of dicts: pillars: [{r_pillar: 80, ...}, {r_pillar: 20, ...}]
        # 2. Dict with lists: pillars: {r_pillar: [80, 20], phi_pillar: [0, pi], ...}
        # 3. Random generation: pillars: {make_many: true, N_pillar: 100, ...}
        pillars_config = config.get('pillars', [])
        
        # Check if we should generate random pillars
        make_many = False
        if isinstance(pillars_config, dict):
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
            rmin_eff = max(disk_params['rin'], rmin) if rmin is not None else disk_params['rin']
            r_pillar_random = np.clip(r_pillar_random, rmin_eff, disk_params['rout'])
            
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
            print(f"  r: mean={r_mean:.1f}, sigma={sig_r:.1f} (clipped to [{disk_params['rin']:.1f}, {disk_params['rout']:.1f}])")
            print(f"  phi: uniform [0, 2π]")
        elif isinstance(pillars_config, dict):
            # New format: dict with lists - expand into individual pillars
            expanded_pillars = []
            
            # Get all parameter lists/values
            r_pillar_list = pillars_config.get('r_pillar', [])
            phi_pillar_list = pillars_config.get('phi_pillar', [])
            height_list = pillars_config.get('height', [0.01])
            sigma_r_list = pillars_config.get('sigma_r', [1.0])
            sigma_phi_list = pillars_config.get('sigma_phi', [0.1])
            temp_factor_list = pillars_config.get('temp_factor', [1.5])
            modify_height_list = pillars_config.get('modify_height', [True])
            modify_temp_list = pillars_config.get('modify_temp', [False])
            
            # Convert to lists if they're scalars
            def to_list(x):
                if isinstance(x, list):
                    return x
                return [x]
            
            r_pillar_list = to_list(r_pillar_list)
            phi_pillar_list = to_list(phi_pillar_list)
            height_list = to_list(height_list)
            sigma_r_list = to_list(sigma_r_list)
            sigma_phi_list = to_list(sigma_phi_list)
            temp_factor_list = to_list(temp_factor_list)
            modify_height_list = to_list(modify_height_list)
            modify_temp_list = to_list(modify_temp_list)
            
            # Determine number of pillars (max length of all lists)
            n_pillars = max(len(r_pillar_list), len(phi_pillar_list), len(height_list),
                          len(sigma_r_list), len(sigma_phi_list), len(temp_factor_list),
                          len(modify_height_list), len(modify_temp_list))
            
            # Expand to n_pillars, repeating last value if list is shorter
            def expand_list(lst, n, default):
                if len(lst) == 0:
                    return [default] * n
                if len(lst) == 1:
                    return lst * n
                # Pad with last value if needed
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
            
            # Create individual pillar dicts
            for i in range(n_pillars):
                pillar = {
                    'r_pillar': float(r_pillar_list[i]),
                    'phi_pillar': parse_math_expr(phi_pillar_list[i]),  # Support 'pi', 'pi/2', etc.
                    'height': float(height_list[i]),
                    'sigma_r': float(sigma_r_list[i]),
                    'sigma_phi': parse_math_expr(sigma_phi_list[i]),  # Support 'pi', etc.
                    'temp_factor': float(temp_factor_list[i]),
                    'modify_height': bool(modify_height_list[i]) if not isinstance(modify_height_list[i], str) 
                                   else modify_height_list[i].lower() in ('true', '1', 'yes'),
                    'modify_temp': bool(modify_temp_list[i]) if not isinstance(modify_temp_list[i], str)
                                 else modify_temp_list[i].lower() in ('true', '1', 'yes')
                }
                expanded_pillars.append(pillar)
            
            pillars_config = expanded_pillars
        else:
            # Old format: list of pillar dicts
            for pillar in pillars_config:
                for key in ['r_pillar', 'height', 'sigma_r', 'temp_factor']:
                    if key in pillar:
                        pillar[key] = float(pillar[key])
                # Special handling for phi_pillar and sigma_phi (may contain 'pi')
                for key in ['phi_pillar', 'sigma_phi']:
                    if key in pillar:
                        pillar[key] = parse_math_expr(pillar[key])
                for key in ['modify_height', 'modify_temp']:
                    if key in pillar:
                        if isinstance(pillar[key], str):
                            pillar[key] = pillar[key].lower() in ('true', '1', 'yes')
                        else:
                            pillar[key] = bool(pillar[key])
        
        # Process computation parameters
        comp_params = config.get('computation', {})
        for key in ['wlog_min', 'wlog_max', 'taumax']:
            if key in comp_params:
                comp_params[key] = float(comp_params[key])
        for key in ['nwavelengths', 'ntau']:
            if key in comp_params:
                comp_params[key] = int(comp_params[key])
        if 'use_parallel' in comp_params:
            if isinstance(comp_params['use_parallel'], str):
                comp_params['use_parallel'] = comp_params['use_parallel'].lower() in ('true', '1', 'yes')
            else:
                comp_params['use_parallel'] = bool(comp_params['use_parallel'])
        if 'n_jobs' in comp_params and comp_params['n_jobs'] is not None:
            comp_params['n_jobs'] = int(comp_params['n_jobs'])
        
        # Process plotting parameters
        plot_params = config.get('plotting', {})
        if 'save_results' in plot_params:
            if isinstance(plot_params['save_results'], str):
                plot_params['save_results'] = plot_params['save_results'].lower() in ('true', '1', 'yes')
            else:
                plot_params['save_results'] = bool(plot_params['save_results'])
        if 'save_3d' in plot_params:
            if isinstance(plot_params['save_3d'], str):
                plot_params['save_3d'] = plot_params['save_3d'].lower() in ('true', '1', 'yes')
            else:
                plot_params['save_3d'] = bool(plot_params['save_3d'])
        if 'dpi' in plot_params:
            plot_params['dpi'] = int(plot_params['dpi'])
    
    # Create disk model
    disk = PillarDisk(**disk_params)
    
    # Add pillars from config
    for pillar_cfg in pillars_config:
        # Convert YAML boolean to Python boolean if needed
        if 'modify_height' in pillar_cfg:
            pillar_cfg['modify_height'] = bool(pillar_cfg['modify_height'])
        if 'modify_temp' in pillar_cfg:
            pillar_cfg['modify_temp'] = bool(pillar_cfg['modify_temp'])
        
        disk.add_pillar(**pillar_cfg)
    
    # Wavelength grid (rest frame, Angstroms)
    wlog_min = comp_params.get('wlog_min', 2.0)
    wlog_max = comp_params.get('wlog_max', 5.5)
    nwavelengths = comp_params.get('nwavelengths', 100)
    wavelengths = np.logspace(wlog_min, wlog_max, nwavelengths)
    
    print("Computing SED with pillars...")
    flux_with = disk.compute_sed(wavelengths)
    
    print("Computing SED without pillars...")
    flux_without = disk.compute_sed_no_pillars(wavelengths)
    
    # Computation parameters
    ntau = comp_params.get('ntau', 20)
    taumax = comp_params.get('taumax', 150.0)
    use_parallel = comp_params.get('use_parallel', True)
    n_jobs = comp_params.get('n_jobs', None)
    
    print("Computing lag spectrum with pillars...")
    tau_mean_with, tau_grid, psi_with = disk.compute_lag_spectrum(
        wavelengths, ntau=ntau, taumax=taumax, parallel=use_parallel, n_jobs=n_jobs)
    
    print("Computing lag spectrum without pillars...")
    tau_mean_without, _, psi_without = disk.compute_lag_spectrum_no_pillars(
        wavelengths, ntau=ntau, taumax=taumax, parallel=use_parallel, n_jobs=n_jobs)
    
    # Plot results
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # SED - with and without pillars
    ax = axes[0, 0]
    ax.loglog(wavelengths, flux_without, 'b--', linewidth=2, label='Without pillars', alpha=0.7)
    ax.loglog(wavelengths, flux_with, 'r-', linewidth=2, label='With pillars')
    ax.set_xlabel('Wavelength (Angstroms)')
    ax.set_ylabel('Flux (mJy)')
    ax.set_title('SED')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Lag spectrum - with and without pillars
    ax = axes[0, 1]
    # Only plot where tau > 0
    mask_with = tau_mean_with > 0
    mask_without = tau_mean_without > 0
    ax.plot(wavelengths[mask_without], tau_mean_without[mask_without], 
             'b--', linewidth=2, label='Without pillars', alpha=0.7)
    ax.plot(wavelengths[mask_with], tau_mean_with[mask_with], 
             'r-', linewidth=2, label='With pillars')
    ax.set_xlabel('Wavelength (Angstroms)')
    ax.set_ylabel('Mean Delay (days)')
    ax.set_title('Lag Spectrum')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Delay distribution at a few wavelengths (with pillars)
    ax = axes[1, 0]
    w_indices = [10, 30, 50, 70]
    colors = ['b', 'g', 'r', 'm']
    for i, idx in enumerate(w_indices):
        if idx < len(wavelengths):
            ax.plot(tau_grid, psi_with[:, idx], color=colors[i], 
                   label=f'{wavelengths[idx]:.1f} Å', linewidth=1.5)
    ax.set_xlabel('Delay (days)')
    ax.set_ylabel('Response Function')
    ax.set_title('Delay Distributions (with pillars)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Disk geometry visualization - only show phi=0
    ax = axes[1, 1]
    ax.clear()  # Clear any previous plots
    # Remove any existing legend
    if ax.legend_ is not None:
        ax.legend_.remove()
    
    r_plot = np.logspace(np.log10(disk.rin), np.log10(disk.rout), 100)
    phi_val = 0.0  # Only show phi=0
    h_plot_with = disk.get_height(r_plot, np.full_like(r_plot, phi_val))
    # Also show without pillars for comparison
    saved_pillars = disk.pillars.copy()
    disk.pillars = []
    h_plot_without = disk.get_height(r_plot, np.full_like(r_plot, phi_val))
    disk.pillars = saved_pillars
    
    # Plot only once with clear labels - use handles to ensure single legend entry
    line1 = ax.plot(r_plot, h_plot_without, 'b--', linewidth=2, label='Without pillars', alpha=0.7)
    line2 = ax.plot(r_plot, h_plot_with, 'r-', linewidth=2, label=r'With pillars ($\phi=0$)')
    ax.set_xlabel('Radius (light days)')
    ax.set_ylabel('Height (light days)')
    ax.set_title(r'Disk Geometry ($\phi=0$)')
    ax.set_xscale('log')
    # Create legend explicitly with only the two lines
    ax.legend([line1[0], line2[0]], ['Without pillars', r'With pillars ($\phi=0$)'], 
             loc='best', frameon=True)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    results_filename = plot_params.get('results_filename', 'pillar_disk_results.png')
    dpi = plot_params.get('dpi', 150)
    if plot_params.get('save_results', True):
        plt.savefig(results_filename, dpi=dpi)
        print(f"Results saved to {results_filename}")
    else:
        plt.show()
    
    # Create 3D visualization
    if plot_params.get('save_3d', True):
        print("\nCreating 3D visualization...")
        plot3d_filename = plot_params.get('plot3d_filename', 'pillar_disk_3d.png')
        n_phi_plot = plot_params.get('n_phi_plot', 180)
        n_r_plot = plot_params.get('n_r_plot', 200)
        disk.plot_3d_geometry(show_light_rays=True, filename=plot3d_filename,
                             n_phi_plot=n_phi_plot, n_r_plot=n_r_plot)
    
    # Print summary
    print("\nSummary:")
    print(f"Number of pillars: {len(disk.pillars)}")
    for i, pillar in enumerate(disk.pillars):
        print(f"  Pillar {i+1}: r={pillar['r']:.1f} ld, φ={pillar['phi']:.2f} rad, "
              f"height={pillar['height']:.4f} ld")
    print(f"Wavelength range: {wavelengths[0]:.1f} - {wavelengths[-1]:.1f} Angstroms")
    mask = tau_mean_with > 0
    if np.any(mask):
        print(f"Lag range (with pillars): {tau_mean_with[mask].min():.2f} - {tau_mean_with[mask].max():.2f} days")
    mask = tau_mean_without > 0
    if np.any(mask):
        print(f"Lag range (without pillars): {tau_mean_without[mask].min():.2f} - {tau_mean_without[mask].max():.2f} days")
    print(f"Flux range (with pillars): {flux_with[flux_with>0].min():.2e} - {flux_with.max():.2e} mJy")
    print(f"Flux range (without pillars): {flux_without[flux_without>0].min():.2e} - {flux_without.max():.2e} mJy")


if __name__ == '__main__':
    import sys
    # Allow config file to be specified as command line argument
    config_file = sys.argv[1] if len(sys.argv) > 1 else 'config.yaml'
    main(config_file=config_file)

