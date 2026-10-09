# Handoff: modeling the transient UV bump with the spiral-arm disk model

## 0. Read first

- Read both CLAUDE.md files (global `~/.claude/CLAUDE.md` + the rules it
  includes, and the repo `CLAUDE.md`). THE FIRST RULE: if scope, file, or
  parameter variant is ambiguous, ask one short question (1-3 options)
  before editing or running anything.
- Plain single-line `python3` commands; the transient scripts run in the
  user's usual conda env (pypeit nominally). The `pyroa` env is only for
  PyROA lag fitting, which is NOT the focus now.

## 1. The goal

Reproduce the transient UV/optical bump of the Liu et al. "Main Target"
(R.A. 25.1172, Dec +13.3891, z = 0.494) with our spiral-arm disk model.
The 2019 eBOSS epoch (MJD 58523) shows the bump; SDSS 2000 (51882) and
DESI 2021 (59518) are quiescent power laws. The model idea: a slow TDE
heats a pillar on the lamp-irradiated bowl disk at r_p = 2 light days;
differential rotation shears it into two finite spiral arms; the raised
arm also occults the inner disk along the line of sight at high
inclination, suppressing the UV.

## 2. What the data show (the target to fit)

The 2019-minus-quiescent difference spectrum is a smooth, Fe II-free
continuum bump peaking near rest 3300 A, with a sharp blue edge
(2620 -> 3050 A), a red side that is ~zero by rest 5000 A, and
L ~ 1e45 erg/s. It is NARROWER than any single-temperature blackbody.
That narrowness is the central modeling challenge; the paper's own
single blackbody (10453 K) and modified blackbody fail the same way.

## 3. The model machinery

`transient/transient_disk.py` defines `TransientPillarDisk(PillarDisk)`:

- Sheared arm geometry: `spiral_shear` A = 2 rad, `arm_length` 2.5 rad
  taper, arm thickness `sigma_phi` = 0.9, envelope `sigma_r` = 1.5 ld,
  height 1.2 ld at r_p = 2 ld.
- TDE heating of the arm: `heat_mode: blackbody` (quadrature T^4 add,
  `pillar_temp` ~ 10.5 kK) — the ADOPTED mode. `heat_mode: balmer`
  (hydrogen recombination continuum with Doppler-smeared edges) is kept
  as a switch; its edges must use the closed erfc/erfcx forms, discrete
  node sums alias into staircases.
- Occultation: `opaque` geometric blocking (ADOPTED) or `balmer_abs`
  (translucent Balmer-opacity screen from the line-of-sight column).
- `disks_from_config(cfg)` -> (flare, quiet, dmpc).
  Config: `transient/config_transient.yaml` (i = 80 deg, normalization
  anchored on the red side 5500-6050 A observed).

## 4. Established results (chi2 vs the 2019 difference spectrum)

- Pure single-T blackbody bump saturates at chi2 ~ 850; its red tail
  cannot drop below ~35-40% of peak (data < 5%). Planck-width floor.
- Balmer-opacity screen + arm recombination emission: chi2 = 272
  (thick arm) / 307 (thin arm, i = 70). Absorption without arm emission
  is rejected (chi2 >= 3371).
- Geometric occultation + blackbody-heated arm (the adopted figure
  configuration): chi2 = 397 at i = 80 — competitive, better red side,
  overshoots the blue edge.
- `fit_components.py` NNLS decomposition (hot BB + Balmer bound-free +
  occultation) reaches chi2/dof ~ 133: the bump wants Balmer
  red-truncation plus occultation, not a pure blackbody.
- Viewing-angle selection effect: UV suppression switches on only for
  i above ~arctan(r_p/h) (~60-70 deg). The SED figure draws the whole
  inclination family with a colorbar.

## 5. Open threads for the new session

- The red-side residual: the model difference spectrum flattens at ~10
  (1e-17 cgs) redward of 4500 A where the data go to zero. Knobs:
  `balmer_bb_frac`, `n3_frac`, the hardcoded 0.12 Paschen/Balmer ratio.
- Can the adopted blackbody-mode arm be made to reproduce the sharp
  blue edge, or does that force the Balmer-opacity mode back in?
- Ask the user what to attack first; do not assume.

## 6. Key scripts and files (all under `transient/`)

- `make_transient_figure.py` — the main 3-panel SED figure
  (`plots/transient_sed_figure.png`), inclination family + colorbar.
- `plot_bestfit.py` (BEST = i 58, phi 90, T 9500), `scan_ipT.py`,
  `fit_components.py`, `plot_sed_tempsweep.py`, `plot_sed_phisweep.py`.
- Geometry/appearance maps: `plot_pillar_xy.py`, `plot_height_xy.py`,
  `plot_los_view.py`. Tests: `test_transient.py`.
- Spectra: `data/epoch_{sdss2000,eboss2019,desi2021}.npz`
  (keys wave_obs, flam [1e-17 cgs], ivar, mjd, label).
- UNIT GOTCHA: the parent `compute_sed` returns rest-frame f_lambda
  x 1e26, NOT mJy, despite labels elsewhere.

## 7. Lag work (done; deprioritized for now)

ZTF g/r/i light curves are downloaded and lags measured: ICCF
(`lag_iccf.py`; light curves MUST be detrended or the CCF is
trend-dominated) and PyROA (`run_pyroa.py`, pyroa env, gridsize=1000;
diagnostics via `plot_pyroa_diagnostics.py`, which wraps PyROA built-ins
PG0844-style). Result: r lags g by +3.0 +/- 0.7 d, i by +7.5 -1.6/+1.5 d;
the model lag spectrum (`compute_model_lags.py`, `plot_lag_overlay.py`)
under-predicts these by a factor ~6-8 and the arm barely affects it.
Details in `HANDOFF_lightcurve_lags.md`. Park this unless the user asks.
