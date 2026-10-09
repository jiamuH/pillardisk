#!/usr/bin/env python3
"""
transient_disk.py - PillarDisk subclass for the transient-SED experiment.

Adds two pieces of physics to the parent model (pillar_disk.PillarDisk),
without modifying the parent:

1. Heated pillar: a per-pillar absolute blackbody temperature `pillar_temp`
   (e.g. from a TDE of an embedded star). Added in quadrature over the
   pillar's Gaussian footprint: T^4 -> T^4 + pillar_temp^4 * footprint.

2. Observer line-of-sight occultation: a vectorized ray-march along the
   observer direction e = (sin i, 0, cos i). Disk cells whose sight line
   passes below the analytic surface (bowl + pillar Gaussians) are blocked,
   with a smooth penumbra. This is what suppresses the inner-disk UV when a
   tall, azimuthally extended pillar sits on the near side (phi = 0).

Also provides a vectorized compute_sed() that reproduces the parent's
discretization exactly (verified in test_transient.py) and applies the
occultation mask.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pillar_disk import PillarDisk, C, DAY, FNU_TO_MJY  # noqa: E402


def disks_from_config(cfg, nr=None, nphi=None):
    """Build (flare_disk, quiet_disk, dmpc) from a loaded config dict.

    The flare disk carries the pillar and, if transient.flare_tv1_scale is
    set, a brightened viscous temperature (bluer-when-brighter disk in the
    transient epoch). The quiet disk is the same disk with no pillar and
    the unscaled tv1.
    """
    obs = cfg['observation']
    dmpc = obs['dmpc']
    if dmpc == 'auto':
        from astropy.cosmology import FlatLambdaCDM
        dmpc = float(FlatLambdaCDM(H0=70, Om0=0.3)
                     .luminosity_distance(obs['redshift']).value)
    cosi = float(np.cos(np.radians(obs['inclination'])))
    occ = cfg['transient']['occultation']
    tr = cfg['transient']
    common = dict(
        rin=cfg['disk']['rin'], rout=cfg['disk']['rout'],
        nr=nr or cfg['disk']['nr'], nphi=nphi or cfg['disk']['nphi'],
        h1=cfg['disk']['h1'], r0=cfg['disk']['r0'], beta=cfg['disk']['beta'],
        tirrad_tvisc_ratio=cfg['temperature']['tirrad_tvisc_ratio'],
        fcol=cfg['temperature']['fcol'], hlamp=cfg['lamp']['hlamp'],
        dmpc=dmpc, cosi=cosi, redshift=obs['redshift'],
        occultation=occ['enabled'], f_trans_obs=occ.get('f_trans_obs', 0.0),
        n_ray=occ.get('n_ray', 300),
        penumbra_frac=occ.get('penumbra_frac', 0.05),
        wall_flux_conserve=tr.get('wall_flux_conserve', True),
        occult_mode=occ.get('mode', 'opaque'),
        tau_edge=occ.get('tau_edge', 5.0), n3_frac=occ.get('n3_frac', 0.0),
        abs_vsmear=occ.get('abs_vsmear', 15000.0),
        dust_law=occ.get('dust_law', 'smc'),
        dust_c=tuple(occ.get('dust_c', (27.0, 4.0, 5.5, 0.0))))
    tv1 = cfg['temperature']['tv1']
    s_hot = float(tr.get('flare_tv1_scale', 1.0))
    quiet = TransientPillarDisk(tv1=tv1, **common)
    flare = TransientPillarDisk(tv1=tv1 * s_hot, **common)
    p = tr['pillar']
    flare.add_pillar(p['r_pillar'], p['phi_pillar'], height=p['height'],
                     sigma_r=p['sigma_r'], sigma_phi=p['sigma_phi'],
                     pillar_temp=p.get('pillar_temp', 0.0),
                     heat_r=p.get('heat_r'),
                     heat_sigma_r=p.get('heat_sigma_r'),
                     heat_sigma_phi=p.get('heat_sigma_phi'),
                     heat_profile=p.get('heat_profile', 2.0),
                     spiral_shear=p.get('spiral_shear', 0.0),
                     arm_length=p.get('arm_length', 0.0),
                     heat_mode=p.get('heat_mode', 'blackbody'),
                     balmer_te=p.get('balmer_te', 1e4),
                     balmer_bb_frac=p.get('balmer_bb_frac', 0.0),
                     balmer_vsmear=p.get('balmer_vsmear', 0.0),
                     balmer_paschen=p.get('balmer_paschen', 0.12))
    return flare, quiet, dmpc


class TransientPillarDisk(PillarDisk):
    """PillarDisk with heated pillars and observer line-of-sight occultation."""

    def __init__(self, *args,
                 occultation: bool = True,
                 f_trans_obs: float = 0.0,
                 n_ray: int = 300,
                 penumbra_frac: float = 0.05,
                 wall_flux_conserve: bool = False,
                 occult_mode: str = 'opaque',
                 tau_edge: float = 5.0,
                 n3_frac: float = 0.0,
                 abs_vsmear: float = 15000.0,
                 dust_law: str = 'smc',
                 dust_c=(27.0, 4.0, 5.5, 0.0),
                 **kwargs):
        """
        Parameters (in addition to PillarDisk's):
        -----------------------------------------
        occultation : bool
            If True, apply observer line-of-sight occultation in compute_sed.
        f_trans_obs : float
            Transmission of the occulting wall (0 = opaque, 1 = transparent).
        n_ray : int
            Number of (log-spaced) samples along each observer ray.
        penumbra_frac : float
            Penumbra softening scale, as a fraction of the tallest pillar height.
        wall_flux_conserve : bool
            If True, conserve the ambient (disk) radiative flux on slanted
            pillar walls: the annulus luminosity is fixed, so spreading it
            over the larger slanted area lowers the ambient temperature by
            (dr/ds)^(1/4). The raised pillar then adds no spurious ambient
            luminosity; only the heated-pillar term and the occultation
            change the SED. The pillar_temp term is NOT flux-conserved
            (it is a specified surface temperature).
        """
        """
        occult_mode : str
            'opaque': geometric gray blocking (original behavior).
            'balmer_abs': the arm is a translucent gas screen; each sight
            line is attenuated by exp(-tau_edge * (N/N_ref) * opacity(lam))
            where N is the ray's path length through the arm and the
            opacity is the hydrogen n=2 bound-free cross-section,
            (lam/3646)^3 blueward of the Balmer edge and zero redward
            (plus an optional n=3 term scaled by n3_frac). The Balmer
            break then comes from the column density along the line of
            sight, a single absorption model.
        """
        super().__init__(*args, **kwargs)
        self.occultation = occultation
        self.f_trans_obs = f_trans_obs
        self.n_ray = n_ray
        self.penumbra_frac = penumbra_frac
        self.wall_flux_conserve = wall_flux_conserve
        self.occult_mode = occult_mode
        self.tau_edge = tau_edge
        self.n3_frac = n3_frac
        self.abs_vsmear = abs_vsmear
        self.dust_law = dust_law
        self.dust_c = tuple(dust_c)

    # ------------------------------------------------------------------
    # Heated pillar
    # ------------------------------------------------------------------
    def add_pillar(self, r_pillar, phi_pillar, height=0.01, sigma_r=1.0,
                   sigma_phi=0.1, modify_height=True, modify_temp=False,
                   temp_factor=1.5, pillar_temp=0.0,
                   heat_sigma_r=None, heat_sigma_phi=None, heat_r=None,
                   heat_profile=2.0, spiral_shear=0.0, arm_length=0.0,
                   heat_mode='blackbody', balmer_te=1e4, balmer_bb_frac=0.0,
                   balmer_vsmear=0.0, balmer_paschen=0.12):
        """Same as PillarDisk.add_pillar, plus:

        pillar_temp : float
            Absolute blackbody temperature of the heated pillar region (K);
            0 = not heated.
        heat_sigma_r, heat_sigma_phi : float or None
            Widths of the HEATED region (the TDE heats the part of the
            pillar near the embedded star). Default None = same as the
            geometric sigma_r / sigma_phi.
        heat_r : float or None
            Center radius of the heated region (the TDE debris orbits at
            the embedded star's radius, which need not coincide with the
            occulting wall). Default None = same as r_pillar.
        heat_profile : float
            Exponent of the heated-region profile exp(-0.5 |u|^p). p = 2 is
            a Gaussian; larger p gives a flat-topped, sharp-edged region
            radiating at nearly uniform temperature, which produces the
            narrowest (single-blackbody) bump.
        spiral_shear : float
            Keplerian shear winding parameter A = Omega(r_pillar) * t
            (radians). Debris stripped from a slowly disrupted embedded star
            shears into TWO spiral arms: the arm center at radius r sits at
            azimuth phi_pillar + A*((r/r_pillar)^(-3/2) - 1), so material
            inside the star's orbit leads (inner arm) and material outside
            trails (outer arm). sigma_r is then the radial envelope of the
            debris spread and sigma_phi the azimuthal thickness of the arm.
            0 = no shear (plain Gaussian pillar, parent behavior).
        arm_length : float
            Gaussian taper (radians) of each arm along the spiral, applied
            to the azimuthal drift |A*((r/r_p)^(-3/2)-1)|. Limits the
            debris stream to two finite arms instead of endless inner
            windings. 0 = no taper.
        heat_mode : str
            'blackbody': the heated region radiates as a pillar_temp
            blackbody added in quadrature (default).
            'balmer': the heated debris emits optically thin hydrogen
            recombination continuum (Balmer + Paschen free-bound) at
            electron temperature balmer_te. This has a sharp red cutoff at
            the 3646 A Balmer edge and an exponential blue decline, giving
            a much narrower bump than any blackbody. pillar_temp then sets
            the amplitude only (peak surface emittance equal to that of a
            pillar_temp blackbody at the Balmer edge).
        balmer_te : float
            Electron temperature (K) of the recombining gas in 'balmer'
            mode; controls the blue-side decline scale.
        balmer_bb_frac : float
            In 'balmer' mode, fraction of the arm emission that is
            thermalized (emitted as a balmer_te blackbody shape) rather
            than optically thin recombination continuum. Smears the sharp
            3646 A edge, mimicking a partially optically thick stream.
        balmer_vsmear : float
            Doppler smearing velocity (km/s) applied to the arm emission
            shape; Keplerian motion of the debris (~10^4 km/s at a few
            light days) rounds the Balmer edge over hundreds of Angstroms.
        balmer_paschen : float
            Paschen/Balmer edge emissivity ratio of the recombination
            continuum (Case B is ~0.1; 0 removes the Paschen term).
        """
        super().add_pillar(r_pillar, phi_pillar, height=height,
                           sigma_r=sigma_r, sigma_phi=sigma_phi,
                           modify_height=modify_height,
                           modify_temp=modify_temp, temp_factor=temp_factor)
        self.pillars[-1]['pillar_temp'] = pillar_temp
        self.pillars[-1]['heat_sigma_r'] = heat_sigma_r or sigma_r
        self.pillars[-1]['heat_sigma_phi'] = heat_sigma_phi or sigma_phi
        self.pillars[-1]['heat_r'] = heat_r or r_pillar
        self.pillars[-1]['heat_profile'] = heat_profile
        self.pillars[-1]['spiral_shear'] = spiral_shear
        self.pillars[-1]['arm_length'] = arm_length
        self.pillars[-1]['heat_mode'] = heat_mode
        self.pillars[-1]['balmer_te'] = balmer_te
        self.pillars[-1]['balmer_bb_frac'] = balmer_bb_frac
        self.pillars[-1]['balmer_vsmear'] = balmer_vsmear
        self.pillars[-1]['balmer_paschen'] = balmer_paschen

    @staticmethod
    def _pillar_footprint(r_2d, phi_2d, pillar, heat=False):
        """Periodic 2D Gaussian footprint of a pillar (same form as
        PillarDisk.get_height, 3-image sum in phi). With heat=True, use the
        heated-region widths instead of the geometric ones."""
        if heat:
            dr = r_2d - pillar.get('heat_r', pillar['r'])
            sig_r = pillar.get('heat_sigma_r', pillar['sigma_r'])
            sig_phi = pillar.get('heat_sigma_phi', pillar['sigma_phi'])
            p_exp = pillar.get('heat_profile', 2.0)
        else:
            dr = r_2d - pillar['r']
            sig_r = pillar['sigma_r']
            sig_phi = pillar['sigma_phi']
            p_exp = 2.0
        # two-armed spiral from Keplerian shear of the stripped debris:
        # arm center azimuth drifts as A*((r/r_p)^(-3/2) - 1)
        shear = pillar.get('spiral_shear', 0.0)
        phi_center = pillar['phi']
        taper = 1.0
        if shear != 0.0:
            rr = np.maximum(np.asarray(r_2d, dtype=float), 1e-6)
            drift = shear * ((rr / pillar['r']) ** -1.5 - 1.0)
            phi_center = phi_center + drift
            arm_len = pillar.get('arm_length', 0.0)
            if arm_len > 0.0:
                # finite debris stream: two arms, each tapering along the
                # spiral instead of winding indefinitely
                taper = np.exp(-0.5 * (drift / arm_len) ** 2)
        dphi = np.mod(phi_2d - phi_center + np.pi, 2 * np.pi) - np.pi
        gauss_phi = (np.exp(-0.5 * np.abs(dphi / sig_phi) ** p_exp)
                     + np.exp(-0.5 * np.abs((dphi - 2 * np.pi) / sig_phi) ** p_exp)
                     + np.exp(-0.5 * np.abs((dphi + 2 * np.pi) / sig_phi) ** p_exp))
        return np.exp(-0.5 * np.abs(dr / sig_r) ** p_exp) * gauss_phi * taper

    def get_height(self, r, phi):
        """Parent height with the subclass footprint (supports spiral arms).
        Note: the parent's lamp-shadow mask still approximates the pillar as
        an azimuthal Gaussian at r_pillar; with the weak lamp used here
        (small tirrad_tvisc_ratio) the error is negligible."""
        r = np.asarray(r)
        phi = np.asarray(phi)
        if r.ndim == 1 and phi.ndim == 1:
            r_2d, phi_2d = np.meshgrid(r, phi, indexing='ij')
        else:
            r_2d, phi_2d = np.broadcast_arrays(r, phi)
        h = np.interp(r_2d.flatten(), self.r,
                      self.h_base).reshape(r_2d.shape)
        for pillar in self.pillars:
            if pillar['modify_height']:
                h = h + pillar['height'] * self._pillar_footprint(
                    r_2d, phi_2d, pillar)
        return h

    def get_temperature(self, r, phi, compute_shadows=True):
        """Parent temperature plus heated-pillar term in quadrature."""
        T = super().get_temperature(r, phi, compute_shadows=compute_shadows)
        heated = [p for p in self.pillars
                  if p.get('pillar_temp', 0.0) > 0.0
                  and p.get('heat_mode', 'blackbody') == 'blackbody']
        if not heated:
            return T
        r = np.asarray(r)
        phi = np.asarray(phi)
        if r.ndim == 1 and phi.ndim == 1:
            r_2d, phi_2d = np.meshgrid(r, phi, indexing='ij')
        else:
            r_2d, phi_2d = np.broadcast_arrays(r, phi)
        T4 = T ** 4
        for p in heated:
            T4 = T4 + p['pillar_temp'] ** 4 * self._pillar_footprint(
                r_2d, phi_2d, p, heat=True)
        return T4 ** 0.25

    def surface_temperature_map(self, r, phi):
        """Effective surface temperature INCLUDING the TDE pillar heating,
        for visualization only (does not affect the SED). get_temperature
        omits the heating in 'balmer' emission mode, so the heated arm would
        otherwise look no hotter than the ambient disk; here pillar_temp is
        added over the heated footprint for every heated pillar."""
        T = self.get_temperature(r, phi)
        r = np.asarray(r)
        phi = np.asarray(phi)
        if r.ndim == 1 and phi.ndim == 1:
            r_2d, phi_2d = np.meshgrid(r, phi, indexing='ij')
        else:
            r_2d, phi_2d = np.broadcast_arrays(r, phi)
        T4 = T ** 4
        for p in self.pillars:
            Tp = p.get('pillar_temp', 0.0)
            # blackbody mode is already in get_temperature; add balmer here
            if Tp > 0.0 and p.get('heat_mode', 'blackbody') != 'blackbody':
                T4 = T4 + Tp ** 4 * self._pillar_footprint(
                    r_2d, phi_2d, p, heat=True)
        return T4 ** 0.25

    # ------------------------------------------------------------------
    # Observer line-of-sight occultation
    # ------------------------------------------------------------------
    def analytic_height(self, r, phi):
        """Closed-form surface height h(r, phi) = bowl + pillar Gaussians,
        valid at arbitrary off-grid positions (no interpolation)."""
        h = self.h1 * (np.asarray(r) / self.r0) ** self.beta
        for pillar in self.pillars:
            if pillar['modify_height']:
                h = h + pillar['height'] * self._pillar_footprint(r, phi, pillar)
        return h

    def compute_observer_occultation(self):
        """
        Visibility mask (nr, nphi): 1 = cell fully visible, f_trans_obs = cell
        fully behind an occulting wall.

        Ray-march from every cell toward the observer along
        e = (sin i, 0, cos i); the cell is blocked where the ray dips below
        the analytic surface. Log-spaced samples resolve nearby walls finely.
        """
        r_2d, phi_2d = np.meshgrid(self.r, self.phi, indexing='ij')
        if len(self.pillars) == 0:
            return np.ones_like(r_2d)

        h_2d = self.get_height(r_2d, phi_2d)
        x0 = r_2d * np.cos(phi_2d)
        y0 = r_2d * np.sin(phi_2d)
        z0 = h_2d

        sini = max(self.sini, 1e-6)
        s_max = 2.0 * self.rout / sini
        # Log-spaced samples: fine near the cell, coarse far out
        s_samples = np.logspace(-3, np.log10(s_max), self.n_ray)

        max_depth = np.full(r_2d.shape, -np.inf)
        for sk in s_samples:
            x = x0 + sk * self.sini
            z = z0 + sk * self.cosi
            rs = np.hypot(x, y0)
            H = self.analytic_height(rs, np.arctan2(y0, x))
            depth = np.where(rs < self.rout, H - z, -np.inf)
            np.maximum(max_depth, depth, out=max_depth)

        h_ref = max(p['height'] for p in self.pillars)
        penumbra = self.penumbra_frac * h_ref + 1e-12
        blocked = np.clip(max_depth / penumbra, 0.0, 1.0)
        return 1.0 - blocked * (1.0 - self.f_trans_obs)

    def compute_los_column(self):
        """
        Path length (light days) of each cell's sight line through the
        raised arm material (i.e. below the arm surface but above the base
        bowl). Used by occult_mode='balmer_abs': the ray's column density
        is proportional to this path length.
        """
        r_2d, phi_2d = np.meshgrid(self.r, self.phi, indexing='ij')
        if len(self.pillars) == 0:
            return np.zeros_like(r_2d)
        h_2d = self.get_height(r_2d, phi_2d)
        x0 = r_2d * np.cos(phi_2d)
        y0 = r_2d * np.sin(phi_2d)
        z0 = h_2d

        sini = max(self.sini, 1e-6)
        s_max = 2.0 * self.rout / sini
        s_samples = np.logspace(-3, np.log10(s_max), self.n_ray)
        ds = np.diff(np.concatenate([[0.0], s_samples]))

        h_ref = max(p['height'] for p in self.pillars)
        penumbra = self.penumbra_frac * h_ref + 1e-12
        column = np.zeros_like(r_2d)
        for sk, dsk in zip(s_samples, ds):
            x = x0 + sk * self.sini
            z = z0 + sk * self.cosi
            rs = np.hypot(x, y0)
            phis = np.arctan2(y0, x)
            H = self.analytic_height(rs, phis)
            h_bowl = self.h1 * (rs / self.r0) ** self.beta
            inside = np.clip((H - z) / penumbra, 0.0, 1.0) * (H > h_bowl + 1e-6)
            column += np.where(rs < self.rout, inside, 0.0) * dsk
        return column

    @staticmethod
    def smc_extinction_shape(wavelength):
        """SMC extinction curve xi(lambda) = A_lambda/A_B from the Pei (1992)
        parameterization (their Table 4), evaluated at rest wavelength in
        Angstroms. Smooth, steeply rising into the UV, no 2175 A bump."""
        lam_um = wavelength * 1e-4
        terms = [(185.0, 0.042, 90.0, 2.0),
                 (27.0, 0.08, 5.50, 4.0),
                 (0.005, 0.22, -1.95, 2.0),
                 (0.010, 9.7, -1.95, 2.0),
                 (0.012, 18.0, -1.80, 2.0),
                 (0.030, 25.0, 0.00, 2.0)]
        xi = 0.0
        for a, li, b, n in terms:
            xi += a / ((lam_um / li) ** n + (li / lam_um) ** n + b)
        return xi

    @staticmethod
    def smallgrain_extinction_shape(wavelength, c1, c2, c3, c4):
        """Small-grain extinction curve A_lambda/A_V (Ma et al. 2026, Eq. 1;
        a Pei-style sum of Drude terms with a free far-UV exponent c2). It is
        normalized to 1 at V (0.55 um) by construction; c2 sets the far-UV
        steepness, c4 the 2175 A bump strength. wavelength in Angstroms."""
        lam = wavelength * 1e-4  # micron
        t1 = c1 / ((lam / 0.08) ** c2 + (0.08 / lam) ** c2 + c3)
        norm2 = 233.0 * (1.0 - c1 / (6.88 ** c2 + 0.145 ** c2 + c3)
                         - c4 / 4.60)
        t2 = norm2 / ((lam / 0.046) ** 2 + (0.046 / lam) ** 2 + 90.0)
        t3 = c4 / ((lam / 0.2175) ** 2 + (0.2175 / lam) ** 2 - 1.95)
        return t1 + t2 + t3

    def _dust_opacity(self, wavelength):
        """Wavelength dependence of a DUSTY screen, normalized to 1 at 5500 A
        so tau_edge is the V-band optical depth of a central chord through the
        arm (occult_mode='dust_abs'). dust_law selects SMC or the small-grain
        (Ma et al. 2026) law."""
        if self.dust_law == 'smallgrain':
            return self.smallgrain_extinction_shape(wavelength, *self.dust_c)
        return (self.smc_extinction_shape(wavelength)
                / self.smc_extinction_shape(5500.0))

    def _balmer_opacity(self, wavelength):
        """Wavelength dependence of the screen opacity: hydrogen n=2
        bound-free, (lam/3646)^3 blueward of the Balmer edge, plus an
        optional n=3 (Paschen) term scaled by n3_frac. Each edge is
        Doppler-smeared by the arm's orbital motion (abs_vsmear, km/s) so
        the transmitted spectrum has no discontinuity at the edge."""
        from scipy.special import erfc
        su = self.abs_vsmear * 1e5 / C

        def edge_step(lam_edge):
            if su <= 0.0:
                return 1.0 if wavelength < lam_edge else 0.0
            sig = lam_edge * su
            return 0.5 * erfc((wavelength - lam_edge) / (np.sqrt(2.0) * sig))

        s = (wavelength / 3646.0) ** 3 * edge_step(3646.0)
        s += self.n3_frac * (wavelength / 8204.0) ** 3 * edge_step(8204.0)
        return s

    # ------------------------------------------------------------------
    # Vectorized SED
    # ------------------------------------------------------------------
    def _sed_weights(self, apply_occultation):
        """Wavelength-independent per-cell weight da * foreshortening * mask,
        shape (nr-1, nphi), reproducing the parent's discretization."""
        r_2d, phi_2d = np.meshgrid(self.r, self.phi, indexing='ij')
        h_2d = self.get_height(r_2d, phi_2d)

        dr = self.r[1:] - self.r[:-1]                       # (nr-1,)
        dh = h_2d[1:, :] - h_2d[:-1, :]                     # (nr-1, nphi)
        ds = np.sqrt(dr[:, None] ** 2 + dh ** 2)
        sintilt = dh / (ds + 1e-10)
        costilt = dr[:, None] / (ds + 1e-10)

        px = -sintilt * np.cos(self.phi)[None, :]
        pz = costilt
        dot = self.ex * px + self.ez * pz
        dot = np.where(dot > 0.0, dot, 0.0)                 # parent's dot>0 cut

        da = ds * self.r[:-1, None] * self.dphi
        weight = da * dot
        if apply_occultation:
            occ = self.compute_observer_occultation()
            weight = weight * occ[:-1, :]
        return weight, costilt

    def compute_sed(self, wavelengths):
        """Vectorized SED in mJy (same units/discretization as the parent),
        with observer occultation applied if self.occultation is True."""
        wavelengths = np.atleast_1d(np.asarray(wavelengths, dtype=float))
        r_2d, phi_2d = np.meshgrid(self.r, self.phi, indexing='ij')
        have_pillars = len(self.pillars) > 0
        absorbing = (self.occultation and have_pillars
                     and self.occult_mode in ('balmer_abs', 'dust_abs'))
        weight, costilt = self._sed_weights(self.occultation and have_pillars
                                            and not absorbing)
        col_scaled = None
        if absorbing:
            # column in units of a central chord through the arm
            col_ref = 2.0 * max(p['sigma_r'] for p in self.pillars)
            col_scaled = (self.compute_los_column()[:-1, :] / col_ref)

        # Ambient temperature (parent: viscous + irradiation + lamp shadows)
        T_amb = PillarDisk.get_temperature(self, r_2d, phi_2d)[:-1, :]
        if self.wall_flux_conserve:
            # Fixed annulus VISCOUS luminosity spread over the slanted wall
            # area: tv^4 -> tv^4 * (dr/ds). The irradiation part is NOT
            # conserved (intercepted lamp flux genuinely scales with the
            # wall area and tilt).
            tv = np.interp(r_2d[:-1, :].ravel(), self.r,
                           self.tv_base).reshape(costilt.shape)
            tx4 = np.clip(T_amb ** 4 - tv ** 4, 0.0, None)
            T_amb = (tv ** 4 * np.clip(costilt, 0.0, 1.0) + tx4) ** 0.25
        # Heated-pillar terms: blackbody pillars add in quadrature to T;
        # balmer pillars contribute additive optically thin recombination
        # continuum (precompute their footprint-weighted areas here)
        T4 = T_amb ** 4
        balmer_terms = []  # (weighted area sum, peak emittance, T_e)
        for p in self.pillars:
            Tp = p.get('pillar_temp', 0.0)
            if Tp <= 0.0:
                continue
            g = self._pillar_footprint(r_2d[:-1, :], phi_2d[:-1, :], p,
                                       heat=True)
            if p.get('heat_mode', 'blackbody') == 'balmer':
                amp = float(self.planck_function(
                    3640.0, np.array([Tp * self.fcol])) / self.fcol ** 4)
                balmer_terms.append((np.sum(g * weight), amp,
                                     p.get('balmer_te', 1e4),
                                     p.get('balmer_bb_frac', 0.0),
                                     p.get('balmer_vsmear', 0.0),
                                     p.get('balmer_paschen', 0.12)))
            else:
                T4 = T4 + Tp ** 4 * g
        Tcol = T4 ** 0.25 * self.fcol
        f4 = self.fcol ** 4

        ld_to_cm = C * DAY
        d_cm = self.d * ld_to_cm
        norm = ld_to_cm ** 2 / d_cm ** 2 * FNU_TO_MJY

        flux = np.empty(wavelengths.shape)
        for iw, w in enumerate(wavelengths):
            B_nu = self.planck_function(w, Tcol) / f4
            if col_scaled is not None:
                kappa = (self._dust_opacity(w)
                         if self.occult_mode == 'dust_abs'
                         else self._balmer_opacity(w))
                if kappa > 0.0:
                    B_nu = B_nu * np.exp(-self.tau_edge * kappa * col_scaled)
            total = np.sum(B_nu * weight)
            for area_w, amp, te, fbb, vsmear, paschen in balmer_terms:
                total += area_w * amp * self._arm_emission_shape(
                    w, te, fbb, vsmear, paschen)
            flux[iw] = total * norm
        return flux

    def _arm_emission_shape(self, w, te, fbb, vsmear, paschen_frac=0.12):
        """Recombination + thermalized arm emission shape at wavelength w,
        Doppler-smeared by the debris orbital motion (vsmear in km/s).
        The smeared one-sided exponential is evaluated in closed form
        (exponentially modified Gaussian), so the edge is exactly smooth."""
        if vsmear <= 0.0:
            shape = (1.0 - fbb) * self._balmer_shape(w, te, paschen_frac)
        else:
            from pillar_disk import H as HP, K as KB, ANGSTROM as ANG
            e_photon = HP * C / (w * ANG)
            kte = KB * te
            sig = e_photon * vsmear * 1e5 / C
            s = self._smeared_edge(e_photon, HP * C / (3646.0 * ANG),
                                   kte, sig)
            s += paschen_frac * self._smeared_edge(
                e_photon, HP * C / (8204.0 * ANG), kte, sig)
            e_b = HP * C / (3646.0 * ANG)
            e_p = HP * C / (8204.0 * ANG)
            s_edge = 1.0 + paschen_frac * np.exp(-(e_b - e_p) / kte)
            shape = (1.0 - fbb) * s / s_edge
        if fbb > 0.0:
            # thermalized part: balmer_te blackbody, edge-normalized
            shape += fbb * float(
                self.planck_function(w, np.array([te]))
                / self.planck_function(3646.0, np.array([te])))
        return shape

    @staticmethod
    def _smeared_edge(e_photon, e_edge, kte, sig):
        """One-sided exponential exp(-(E-E_edge)/kTe) for E > E_edge,
        convolved with a Gaussian of width sig: the exponentially modified
        Gaussian, computed in an overflow-safe form."""
        from scipy.special import erfc, erfcx
        de = e_photon - e_edge
        b = (sig / kte - de / sig) / np.sqrt(2.0)
        if b > 0.0:
            return 0.5 * erfcx(b) * np.exp(-de ** 2 / (2.0 * sig ** 2))
        a = sig ** 2 / (2.0 * kte ** 2) - de / kte
        return 0.5 * np.exp(a) * erfc(b)

    @staticmethod
    def _balmer_shape(wavelength, te, paschen_frac=0.12):
        """Optically thin hydrogen free-bound continuum shape (Balmer +
        Paschen), normalized to 1 just blueward of the 3646 A Balmer edge.
        The recombination continuum exists only BLUEWARD of each edge (sharp
        jump up at the edge, zero on its red side), declining toward the blue
        with scale k*T_e as exp(-(E - E_edge)/kT_e)."""
        from pillar_disk import H, C, K, ANGSTROM
        e_photon = H * C / (wavelength * ANGSTROM)
        e_balmer = H * C / (3646.0 * ANGSTROM)
        e_paschen = H * C / (8204.0 * ANGSTROM)
        kte = K * te
        s = 0.0
        if e_photon > e_balmer:
            s += np.exp(-(e_photon - e_balmer) / kte)
        if e_photon > e_paschen:
            s += paschen_frac * np.exp(-(e_photon - e_paschen) / kte)
        s_edge = 1.0 + paschen_frac * np.exp(-(e_balmer - e_paschen) / kte)
        return s / s_edge

    def arm_velocity_field(self, mbh_msun):
        """Line-of-sight velocity and emission weight of every disk cell
        for the HEATED arm footprint, for building the geometric Doppler
        kernel of the arm emission.

        Circular Keplerian motion v = v_K(r) phi_hat with
        v_K = c sqrt(r_g / r); observer at e = (sin i, 0, cos i); so
        v_los = -v_K(r) sin(i) sin(phi), positive = receding. The weight
        is the same per-cell factor compute_sed applies to the arm
        emission: heated footprint x area x foreshortening x occultation
        mask. Relativistic corrections are neglected.

        Returns (v_los [km/s], weight), both flattened over the heated
        pillars' cells; weights are NOT normalized. v_los scales exactly
        as sqrt(mbh_msun), so callers can rescale rather than recompute.
        """
        G_CGS = 6.6743e-8
        MSUN_G = 1.989e33
        r_2d, phi_2d = np.meshgrid(self.r, self.phi, indexing='ij')
        weight, _ = self._sed_weights(self.occultation
                                      and len(self.pillars) > 0)
        rg_ld = G_CGS * mbh_msun * MSUN_G / C ** 2 / (C * DAY)
        v_kep = C * 1e-5 * np.sqrt(rg_ld / r_2d[:-1, :])       # km/s
        v_los = -v_kep * self.sini * np.sin(phi_2d[:-1, :])
        v_all, w_all = [], []
        for p in self.pillars:
            if p.get('pillar_temp', 0.0) <= 0.0:
                continue
            g = self._pillar_footprint(r_2d[:-1, :], phi_2d[:-1, :], p,
                                       heat=True)
            v_all.append(v_los.ravel())
            w_all.append((g * weight).ravel())
        if not v_all:
            return np.array([]), np.array([])
        return np.concatenate(v_all), np.concatenate(w_all)

    def ring_velocity_field(self, mbh_msun, r_center, sigma_r):
        """Axisymmetric counterpart of arm_velocity_field: line-of-sight
        velocity and emission weight for a full Keplerian ring (Gaussian
        radial envelope around r_center, uniform in azimuth) — the arm
        gas BEFORE shearing/heating. Weights use the same per-cell area x
        foreshortening factors (no occultation: quiescent state has no
        raised arm). Returns (v_los [km/s], weight), v_los scaling as
        sqrt(mbh_msun)."""
        G_CGS = 6.6743e-8
        MSUN_G = 1.989e33
        r_2d, phi_2d = np.meshgrid(self.r, self.phi, indexing='ij')
        weight, _ = self._sed_weights(False)
        rg_ld = G_CGS * mbh_msun * MSUN_G / C ** 2 / (C * DAY)
        v_kep = C * 1e-5 * np.sqrt(rg_ld / r_2d[:-1, :])       # km/s
        v_los = -v_kep * self.sini * np.sin(phi_2d[:-1, :])
        g = np.exp(-0.5 * ((r_2d[:-1, :] - r_center) / sigma_r) ** 2)
        return v_los.ravel(), (g * weight).ravel()

    def compute_occulted_fraction(self, wavelengths):
        """Fraction of flux removed by occultation at each wavelength
        (same temperatures, mask on vs off)."""
        saved = self.occultation
        try:
            self.occultation = True
            f_on = self.compute_sed(wavelengths)
            self.occultation = False
            f_off = self.compute_sed(wavelengths)
        finally:
            self.occultation = saved
        return 1.0 - f_on / np.maximum(f_off, 1e-300)
