# Answers to Your Questions

## 1. How to Speed Up Lag Spectrum Computation

### Parameters that affect speed (in order of impact):

1. **`nr`** (number of radial points) - **BIGGEST IMPACT**
   - Current: 300
   - Faster: 100-150 (2-3x speedup)
   - Example: `disk = PillarDisk(nr=150, ...)`

2. **`nphi`** (number of azimuthal points) - **SECOND BIGGEST**
   - Current: 360
   - Faster: 90-180 (2-4x speedup)
   - Example: `disk = PillarDisk(nphi=180, ...)`

3. **Number of wavelengths** - **LINEAR SCALING**
   - Current: 100
   - Faster: 50 (2x speedup)
   - Example: `wavelengths = np.logspace(2, 5.5, 50)`

4. **`ntau`** (delay bins) - **MINOR IMPACT**
   - Current: 20
   - Faster: 10 (small speedup)
   - Example: `compute_lag_spectrum(..., ntau=10)`

### Parallel Processing

The code now supports parallel processing! Enable it in `main()`:

```python
# In main() function, change:
use_parallel = True  # Enable parallel processing

tau_mean_with, tau_grid, psi_with = disk.compute_lag_spectrum(
    wavelengths, ntau=20, taumax=150.0, parallel=use_parallel)
```

This can give 4-8x speedup on multi-core CPUs.

### Combined Speedup Example:

```python
# Fast setup (10-20x faster than default)
disk = PillarDisk(
    nr=150,      # Reduced from 300
    nphi=180,    # Reduced from 360
    # ... other params
)

wavelengths = np.logspace(2, 5.5, 50)  # Reduced from 100

tau_mean, tau_grid, psi = disk.compute_lag_spectrum(
    wavelengths, ntau=10, taumax=150.0, parallel=True)
```

## 2. How to Change Basic Parameters

### Option 1: Change in `__init__` default values (lines 56-71)

Edit the default values directly in the `__init__` method:

```python
def __init__(self, 
             rin: float = 0.1,      # Change this
             rout: float = 100.0,   # Change this
             nr: int = 300,         # Change this
             # ... etc
```

### Option 2: Pass parameters when creating disk (RECOMMENDED)

```python
# In main() function, modify the disk creation:
disk = PillarDisk(
    rin=0.1,
    rout=100.0,
    nr=200,          # Changed from 300
    nphi=180,        # Changed from 360
    h1=0.01,
    r0=10.0,
    beta=1.0,        # Changed from 10.0
    tv1=1e4,
    alpha=0.75,
    hlamp=5.0,       # Changed from 0.20
    tx1=8e3,
    dmpc=100.0,
    cosi=0.5,
    fcol=1.0,
    redshift=0.0
)
```

### Option 3: Use `update_params()` method (NEW)

After creating the disk, you can update parameters:

```python
disk = PillarDisk(...)  # Create with defaults

# Later, update parameters
disk.update_params(nr=200, nphi=180, beta=1.0, hlamp=5.0)
```

## 3. Why "Without pillars" appears multiple times

"Without pillars" appears in **3 different panels** (which is correct):

1. **Top Left (SED)**: Shows flux with/without pillars
2. **Top Right (Lag Spectrum)**: Shows delay with/without pillars  
3. **Bottom Right (Disk Geometry)**: Shows height profile with/without pillars

Each panel has its own legend, so you see "Without pillars" 3 times - once per panel. This is **expected behavior** for comparison plots.

If you want to remove it from specific panels, you can modify the plotting code in `main()`:

```python
# To remove from SED plot:
ax.loglog(wavelengths, flux_without, 'b--', linewidth=2, alpha=0.7)  # Remove label='...'

# To remove from lag spectrum:
ax.loglog(..., label='')  # Empty label

# To remove from geometry plot:
ax.plot(..., label='')  # Empty label
```

Or if you only want to show "With pillars" results, simply don't compute the "without" versions.

