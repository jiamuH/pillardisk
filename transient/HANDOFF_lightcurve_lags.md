# Handoff: ATLAS + ZTF light curves → inter-band time-delay (lag) measurement

## Goal
Download the public ATLAS and ZTF light curves for our transient quasar,
measure the **inter-band continuum time delays** (reverberation lags), and
compare them to the pillar-disk model's predicted lag spectrum. The paper
(compiled.pdf, Liu et al.) shows these light curves; we want the lags
ourselves.

## Object (Liu et al. "Main Target")
- R.A. = 25.1172, Dec. = +13.3891 (deg, J2000);  z = 0.494;  V ~ 17 mag.
- Point-like, no host complication. The transient bump is present in the
  2019 eBOSS epoch (MJD 58523); quiescent in 2000 (SDSS, MJD 51882) and
  2021 (DESI, MJD 59518). The optical brightening should bracket ~2019.

## Data to download
- **ZTF** (g, r, i): forced photometry via the IRSA ZTF Forced-Photometry
  Service (batch request at the object coords), or alert-broker light
  curves from Lasair / ALeRCE / ANTARES. Prefer forced photometry for
  clean continuum sampling. ZTF starts ~2018.
- **ATLAS** (cyan `c`, orange `o`): forced photometry at
  https://fallingstar-data.com/forcedphot/ (free account, then submit the
  coords; returns a table). ATLAS starts ~2015-2017.
- Save raw tables and cleaned per-band light curves (MJD, mag/flux, err)
  under `transient/data/` (e.g. `lc_ztf_g.npz`, `lc_atlas_o.npz`, ...).
  Convert mags to flux; clip bad/low-SNR points; remove flags.

## Measuring the lag
- Bands to cross-correlate: ZTF g vs r (and vs i), ATLAS c vs o, and
  cross-survey (ATLAS c vs ZTF g are similar bands). The bluer band should
  **lead** the redder band (disk reverberation: tau ∝ lambda^{4/3}).
- Method: interpolated cross-correlation (ICCF / PyCCF) for the CCF
  centroid + peak, with flux-randomization / random-subset-selection for
  the lag uncertainty; cross-check with JAVELIN (damped-random-walk) or a
  Gaussian-process approach. Watch for the ~1 yr seasonal-gap aliases.
- Deliverable: tau(lambda) between the bands, with errors.

## Connection to the model (why we want this)
- `pillar_disk.py` already predicts lag spectra: `compute_lag_spectrum()`
  (mean delay vs wavelength via response functions). The pillar/arm and
  occultation distort the lag spectrum vs a smooth bowl disk.
- Compare the **measured** ZTF/ATLAS inter-band lags to the model's
  predicted tau(lambda) for the adopted geometry (config_transient.yaml:
  r_p = 2 ld, i ~ 55-70 deg, etc.). This is an independent test of the
  disk size / geometry beyond the SED fitting.

## Repo state / where the SED work left off
- `transient/config_transient.yaml` = current model (first-principles
  blackbody heated arm, `heat_mode: blackbody`; occultation `opaque`;
  normalization anchored red-side 5500-6050 A). The `balmer` emission mode
  is kept as an optional switch, not deleted.
- Key SED scripts: make_transient_figure.py (inclination sweep),
  plot_sed_tempsweep.py / plot_sed_phisweep.py (T and azimuth sweeps),
  plot_bestfit.py, scan_ipT.py (i/phi/T chi^2 grid), fit_components.py
  (NNLS decomposition: hot BB + Balmer b-f + occultation, chi2/dof ~ 133).
- Best SED fit so far: T_peak ~ 9.5 kK, weak occultation; a single BB is
  too broad, so the bump needs Balmer b-f (red-truncation) + occultation.
  Data show a smooth, Fe II-free continuum bump peaking ~3300 A.

## First steps for the new session
1. Get ATLAS + ZTF forced photometry at the coords above; clean and save
   per-band light curves in `transient/data/`.
2. Plot the multi-band light curves; identify the 2019 brightening.
3. Run ICCF/PyCCF between bands → measured lags + errors.
4. Compute the model lag spectrum with `compute_lag_spectrum()` for the
   adopted geometry and overlay.
