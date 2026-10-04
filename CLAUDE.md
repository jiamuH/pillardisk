# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## THE FIRST RULE — ASK WHEN YOU ARE NOT SURE

This rule outranks every other instruction in this file and every habit you may have to "act decisively."

If you are not certain what I want — about the scope, the layout, which file to touch, which variant of a parameter, which of several reasonable interpretations to take — **stop and ask one short question first**. Do not pick a "reasonable" interpretation and act. Do not start editing or running anything until I have answered.

Concrete cases where this is mandatory in this project:

- **Figure layout, panel counts, included files.** Words like "top-bottom layout", "1-column plot", "combined figure" are ambiguous. Ask whether I mean the *panel arrangement inside one figure*, the *page-width class* (`figure` vs `figure*`), or *which configurations to include*.
- **Paper-text JH notes.** Many of these are personal shorthand. If the action is not unambiguous, ask before editing.
- **Reusing a config, deck, or directory that already exists.** Confirm scope: which inputs, what is varied, the minimal version that answers the question. Do not default to whatever is already there.
- **Anything destructive or hard to reverse** (overwriting figures, modifying shared configs, moving files into the paper directory).

When asking, give 1 to 3 concrete options. Do not preamble. Wait for the answer.

This rule has been the source of nearly every avoidable mistake in this project. Re-read it before each non-trivial action.

## Project Overview

This is an AGN accretion disk modeling toolkit that computes time delay spectra and spectral energy distributions (SEDs) for accretion disks with Gaussian pillar bumps. The code models lamp-post irradiated bowl-shaped accretion discs for studying continuum lags in active galactic nuclei.

## Environment

All code runs in the conda environment `pypeit`:
```bash
conda activate pypeit
```

## Architecture

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
