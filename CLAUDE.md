# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is an AGN accretion disk modeling toolkit that computes time delay spectra and spectral energy distributions (SEDs) for accretion disks with Gaussian pillar bumps. The code models lamp-post irradiated bowl-shaped accretion discs for studying continuum lags in active galactic nuclei.

## Environment

All code runs in the conda environment `pypeit`:
```bash
conda activate pypeit
```

## Running the Code

```bash
# Main disk model (uses config.yaml by default)
python3 pillar_disk.py

# With custom config file
python3 pillar_disk.py my_config.yaml

# Emission line velocity-delay map
python3 pillar_line.py

# Time-evolving emission line (rotation effects)
python3 pillar_line_time.py

# Integration with Cloudy photoionization models
python3 pillar_line_cloudy.py
```

## Dependencies

Required: `numpy`, `matplotlib`
Optional: `pyyaml` (config files), `tqdm` (progress bars), `scipy` (for cloudy integration), `pandas` (for cloudy data)

## Architecture

### Core Module: `pillar_disk.py`
- `PillarDisk` class: Main disk model with radial/azimuthal grids
- Key methods:
  - `add_pillar()`: Add Gaussian bumps at (r, phi) positions
  - `compute_sed()`: Calculate flux vs wavelength
  - `compute_lag_spectrum()`: Calculate mean delay vs wavelength with response functions
  - `get_height()` / `get_temperature()`: 2D geometry including pillar modifications
  - `_compute_shadow_mask()`: Ray-tracing for pillar shadows
- Supports parallel processing via `multiprocessing` for wavelength loops

### Extended Modules
- `pillar_line.py`: Computes velocity-delay maps Psi(lambda, tau) for emission lines with Keplerian orbits
- `pillar_line_time.py`: Time-evolving spectra showing rotational barber-pole patterns
- `pillar_line_cloudy.py`: Integrates Cloudy photoionization model outputs with disk geometry

### Configuration: `config.yaml`
All parameters are YAML-configurable:
- Disk geometry: `rin`, `rout`, `nr`, `nphi`, `h1`, `r0`, `beta`
- Temperature: `tv1`, `alpha`, `tx1`, `fcol`, `tirrad_tvisc_ratio`
- Pillars: Manual placement (lists of r, phi, height, sigma) or random generation (`make_many: true`)
- Computation: `use_parallel`, `ntau`, `taumax`, wavelength grid
- Pillar phi values can use `pi` expressions (e.g., `phi_pillar: [0, pi, pi/2]`)

## Performance Tuning

Parameters affecting speed (most to least impact):
1. `nr` (radial points): 100-150 faster, 500+ more accurate
2. `nphi` (azimuthal points): 90-180 faster, 720 more accurate
3. `nwavelengths`: Linear scaling
4. `ntau` (delay bins): Minor impact

Enable parallel processing: `use_parallel: true` in config or `parallel=True` in `compute_lag_spectrum()`.

## Physical Units

- Distances: light days
- Temperatures: Kelvin
- Wavelengths: Angstroms
- Flux: mJy
- Delays: days
- Velocities: km/s

## Known Bugs / TODO

### Foreshortening Effect in Barber Pole Pattern (pillar_line_time_cloudy.py)
**Status**: Bug - needs investigation

The barber pole pattern should show asymmetry between near side and far side of the disk:
- **Near side (φ ≈ 0°)**: Observer sees shadow side of pillar → stronger MgII (low ionization), weaker CIV
- **Far side (φ ≈ 180°)**: Observer sees illuminated side of pillar → stronger CIV (high ionization), weaker MgII

Current implementation in `_process_time_step()` computes `illumination_visibility` based on the dot product of lamp-to-pillar direction with observer direction, but the effect is not clearly visible in the output plots.

**To debug**: Artificially exaggerate the `foreshortening_factor` and `line_asymmetry` parameters. The inclination is set via `cosi` in config (cosi=0.5 corresponds to i=60°, not 45°; for 45° use cosi≈0.707).

**Location**: `pillar_line_time_cloudy.py`, lines ~377-476 in `_process_time_step()`
