# Speed Optimization Tips for Lag Spectrum Computation

## Parameters that affect computation speed:

1. **`nr`** (number of radial points): Most impactful
   - Default: 300
   - Faster: 100-150
   - Slower but more accurate: 500+

2. **`nphi`** (number of azimuthal points): Second most impactful
   - Default: 360
   - Faster: 90-180
   - Slower but more accurate: 720

3. **Number of wavelengths**: Linear scaling
   - Default: 100
   - Faster: 50
   - Slower but more accurate: 200+

4. **`ntau`** (number of delay bins): Minor impact
   - Default: 20-50
   - Faster: 10-20
   - Slower but more accurate: 100+

## Parallel Processing

The code now supports parallel processing over wavelengths using multiprocessing:

```python
# Enable parallel processing (uses all CPU cores by default)
tau_mean, tau_grid, psi = disk.compute_lag_spectrum(
    wavelengths, parallel=True)

# Or specify number of cores
tau_mean, tau_grid, psi = disk.compute_lag_spectrum(
    wavelengths, parallel=True, n_jobs=4)
```

## Example: Fast computation setup

```python
# Create disk with fewer grid points
disk = PillarDisk(
    nr=150,      # Reduced from 300
    nphi=180,    # Reduced from 360
    # ... other parameters
)

# Use fewer wavelengths
wavelengths = np.logspace(2, 5.5, 50)  # Reduced from 100

# Enable parallel processing
tau_mean, tau_grid, psi = disk.compute_lag_spectrum(
    wavelengths, ntau=20, taumax=150.0, parallel=True)
```

This can speed up computation by 5-10x depending on your CPU.

