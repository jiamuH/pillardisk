# Pillar Disk Model

This Python script (`pillar_disk.py`) computes the time delay spectrum and spectral energy distribution (SED) for an AGN accretion disk with Gaussian pillar bumps.

## Overview

The script models an axisymmetric accretion disk with power-law temperature and height profiles, and allows you to add non-axisymmetric Gaussian bumps (pillars) at specific radial and azimuthal positions. The pillars can modify both the disk height and temperature.

## Requirements

```bash
pip install numpy matplotlib
pip install pyyaml  # Optional: for config file support
pip install tqdm    # Optional: for progress bars
```

## Configuration File

**NEW**: All parameters can now be configured via `config.yaml`! No need to edit the code.

Edit `config.yaml` to change:
- Disk geometry parameters (rin, rout, nr, nphi, etc.)
- Temperature parameters
- Lamp post parameters
- Pillar locations and properties
- Computation settings (wavelengths, parallel processing, etc.)
- Plotting options

Then run:
```bash
python3 pillar_disk.py
```

Or specify a different config file:
```bash
python3 pillar_disk.py my_config.yaml
```

## Usage

### Basic Example

```python
from pillar_disk import PillarDisk
import numpy as np

# Create a disk model
disk = PillarDisk(
    rin=0.1,           # inner radius (light days)
    rout=100.0,        # outer radius (light days)
    nr=200,            # number of radial points
    nphi=180,          # number of azimuthal points
    h1=0.01,           # height at reference radius (light days)
    r0=10.0,           # reference radius (light days)
    beta=1.0,          # height power-law index
    tv1=1e4,           # viscous temperature at r0 (K)
    alpha=0.75,        # temperature power-law index
    hlamp=5.0,         # lamp height (light days)
    tx1=8e3,           # irradiation temperature at r0 (K)
    dmpc=100.0,        # distance (Mpc)
    cosi=0.5,          # cos(inclination), 1=face-on
    fcol=1.0,          # color temperature factor
    redshift=0.0       # redshift
)

# Add a pillar at r=20 light days, phi=0
disk.add_pillar(
    r_pillar=20.0,
    phi_pillar=0.0,
    amplitude=0.2,      # 20% height increase
    sigma_r=2.0,        # radial width (light days)
    sigma_phi=0.2,      # azimuthal width (radians)
    modify_height=True,
    modify_temp=True,
    temp_factor=1.3    # 30% temperature increase
)

# Compute SED
wavelengths = np.logspace(2, 5.5, 100)  # 100 to 300,000 Angstroms
flux = disk.compute_sed(wavelengths)

# Compute lag spectrum
tau_mean, tau_grid, psi = disk.compute_lag_spectrum(
    wavelengths, 
    ntau=500, 
    taumax=150.0
)
```

### Running the Example

```bash
python3 pillar_disk.py
```

This will:
1. Create a disk with two pillars
2. Compute the SED and lag spectrum
3. Generate plots saved to `pillar_disk_results.png`

## Parameters

### Disk Parameters

- `rin`, `rout`: Inner and outer disk radii (light days)
- `nr`, `nphi`: Number of radial and azimuthal grid points
- `h1`: Disk height at reference radius (light days)
- `r0`: Reference radius (light days)
- `beta`: Height power-law index (h = h1 * (r/r0)^beta)
- `tv1`: Viscous temperature at reference radius (K)
- `alpha`: Temperature power-law index (T = tv1 * (r0/r)^alpha)
- `hlamp`: Height of the lamp post (light days)
- `tx1`: Irradiation temperature at reference radius (K)
- `dmpc`: Distance to AGN (Mpc)
- `cosi`: Cosine of inclination angle (1 = face-on, 0 = edge-on)
- `fcol`: Color temperature factor (Tcol = fcol * Teff)
- `redshift`: Redshift

### Pillar Parameters

- `r_pillar`: Radial position of pillar center (light days)
- `phi_pillar`: Azimuthal position of pillar center (radians)
- `amplitude`: Amplitude of height bump (fractional)
- `sigma_r`: Radial width of Gaussian (light days)
- `sigma_phi`: Azimuthal width of Gaussian (radians)
- `modify_height`: Whether to modify disk height
- `modify_temp`: Whether to modify disk temperature
- `temp_factor`: Temperature enhancement factor at pillar center

## Output

The script computes:

1. **SED (Spectral Energy Distribution)**: Flux in mJy as a function of wavelength
2. **Lag Spectrum**: Mean time delay (days) as a function of wavelength
3. **Delay Distributions**: Response functions showing the distribution of delays at different wavelengths

## Notes

- The disk is still assumed to be axisymmetric overall, but the pillars break this symmetry
- The computation includes both viscous heating and lamp post irradiation
- The time delay is computed as the light travel time from the lamp to the disk element to Earth
- The response function accounts for the temperature derivative of the Planck function

## References

Based on the FORTRAN code in `bowl.for`, `fnubowl.for`, and `tfbowl.for` which compute continuum lags and SEDs for lamp-post irradiated bowl-shaped accretion discs.

