# Note: the per-ray Cloudy emulator (`ray_emulator.py`)

Context: `src/pillardisk/transient/ray_emulator.py` (console command `emu`),
emulating the per-ray Cloudy runs of `cloudy_rays.py` (`rays`) so that ray
sets for any lamp height can be made without Cloudy. Written 2026-11-02;
add to it rather than re-deriving.

## What is emulated, and what is not (settled 2026-11-02)

- The emulator predicts **per-ray Cloudy output only**, never maps. For each
  ray (from its deck alone: entry flux, n_H along the ray) and each line
  (Ha, Hb, C IV, Mg II, He II), it predicts:
  - `lI`, the total;
  - `lq005..lq995`, the column quantiles of the emission (where along the
    ray it comes out);
  - `lTe`, T_e at the emission centroid;
  - for C IV only, the C3+ layer: `lnl`, `lNab`, `lNtot`, `lthk`.
- It does **not** emulate the ionisation structure itself (ion fractions
  against depth). The Stromgren quantities below are *inputs* computed from
  the deck, not predictions.
- `synth_ray` turns the predictions into a synthetic `.ems` per ray, saved
  in `<run>/emulated.npz`. `load_rays`, `civ_layers` and `line_layers` read
  that file in place of Cloudy's output.
- **Everything after the rays is the same code for Cloudy and emulated
  runs:**
  - face-on emission maps (`faceon_surface`, `faceon_map`);
  - observed maps (`civ_emission_map`);
  - line profiles (`line_profile_data`).
- **Line-of-sight turbulence is free, with no retraining.** Escape toward
  the observer is a post-processing step (`civ_escape` / `line_escape`):
  - f = beta(tau_obs) / <beta>, at `LOSG_INCL` = 0, 45 and 60 deg.
  - Variants: `noshear`, `sim` (Sobolev from the sim's velocities),
    `v_t,z` (vertical turbulence `CIV_VTZ`), and `kep` / `kep0`.
  - See `NOTE_line_escape.md`.
- **Not free: the turbulence inside Cloudy along the lamp ray.** The
  training rays used `VTURB = 10 km/s` and `iterate 2`, so each ray's total
  and its trapping are the 10 km/s ones. A different value needs new Cloudy
  training runs.
- `faceon_compare.png` (`--val`) compares *all-direction emission*
  (`faceon_surface`, no escape weight). It is not the observed map.

## Training and validation

- **Training runs** (`HEIGHT_RUNS`):
  - `train/h0.0`, `train/h0.2`, `edd10_h0.5_adapt`, `train/h1.0`, under
    `/data2/jhuang/runs/cloudy_rays`.
  - Models are cached per training set in `train/emulator_models`.
  - Tables are cached in `<run>/emulator/table.npz` and rebuilt whenever
    `TABLE_VERSION` changes.
- **Hold-one-height-out check:** `emu emulate --val h0.0` (or `h0.2`,
  `h0.5`, `h1.0`). It trains without that height, emulates its ray plan
  into `emulated/<val>_val`, and plots the Cloudy / emulated / ratio maps.
- `synth` takes a running median of each target over 5 rays along alpha
  within each phi column (`SMOOTH_N`). Tested: it costs almost nothing in
  placement accuracy.

## History (so we do not redo it)

1. **TABLE_VERSION <= 3: one Gaussian layer per line** (centroid, 5% and
   95% column).
   - **Failed at h = 0 in the fan (phi 240-300 deg).** There Cloudy lights
     the rim, the ray crosses thin gas, and Cloudy lights a second layer
     far out. One layer cannot hold both, so the fan came out dark in
     He II and C IV.
2. **v4 (commit 4c713d9): 12 column quantiles per line** (0.5% to 99.5%).
   - `column_cdf` puts a PCHIP spline through F(log N). `synth_ray` places
     dF/dN along the ray by the deck's column, so a curve that is
     continuous in column splits into separate layers in depth where the
     gas is thin.
   - **Given the true quantiles it reproduces the far layer.**
   - With the predicted quantiles the fan stayed dark. The far layer sits
     in the last ~0.05 dex of column, and the models were biased short by
     ~0.1 dex there (He II `lq995` predicted 24.39 vs 24.51 true).
     Physically, that edge is an ionisation-front breakthrough threshold.
   - The per-ray totals `lI` were unchanged from v3. h0.0_val map totals
     (Ha, He II, C IV, Mg II): +0.076, +0.044, +0.047, +0.027 dex.
3. **Stromgren as a target offset** (predict the quantile minus a
   Stromgren column): **worse**, with scatter of 0.2-0.4 vs 0.1-0.2 dex.
   Not used.
4. **v5 (commit 22b0e8d): Stromgren columns as extra features.**
   - recombinations rec(s) = integral of alpha_B 1.1 n_H^2 s^2 ds, with
     alpha_B = 2.6e-13;
   - photon supply need = phi s0^2;
   - `lNS<f>`: the column where rec = f x need, for f = 1, 0.3, 0.1, 0.03,
     0.01 (log column at the end + 0.5 if never reached);
   - `lrec`: log rec_end / need.
   - Held-out h0 fan, median quantile bias:

     | quantile | v4 features | v5 features |
     |---|---|---|
     | He II 99.5% | -0.093 | -0.027 |
     | He II 95% | -0.076 | -0.016 |
     | C IV 99.5% | -0.133 | -0.032 |
     | Ha 95% | -0.081 | -0.008 |

   - All-ray errors: ~0.02 -> <0.01 dex.
   - h0.0_val map totals: Ha +0.042, He II +0.004, C IV -0.006,
     Mg II -0.006 dex. **The fan is lit.**

## State after v5 (h0.0_val, fan phi 240-300 deg)

- **Cell by cell the ratio map looks bad, but that is misleading.** The
  median cell is +0.24 to +0.46 dex while the region sums are within
  ~0.1 dex, because many faint cells hold little light. Weighted by light,
  the share of emulated fan light in cells where Cloudy is less than half
  as bright is: He II 3%, Mg II 11%, Ha 18%, C IV 19%. None lands in cells
  Cloudy has dark.
- **The placement test points at the predictions, not the synthesis.** The
  test covered 77 fan rays and measured the fraction of each ray's light
  in the wrong 0.01 r0 depth bin (median [75%]):

  | line | true quantiles | predicted | predicted + smoothing |
  |---|---|---|---|
  | Ha | 0.08 [0.10] | 0.16 [0.21] | 0.18 [0.22] |
  | Hb | 0.04 [0.10] | 0.10 [0.20] | 0.11 [0.21] |
  | C IV | 0.03 [0.06] | 0.06 [0.40] | 0.07 [0.40] |
  | Mg II | 0.11 [0.17] | 0.31 [0.48] | 0.32 [0.49] |
  | He II | 0.02 [0.05] | 0.08 [0.18] | 0.08 [0.19] |

  `synth_ray` is fine. The remaining error is in the predictions: Mg II
  overall, and the C IV far layer on about a quarter of the rays.
- The upper-left wedge (phi 110-160 deg) of C IV is ~2x too faint
  (median cell -0.34 dex).

## Open

- **Mg II placement** and the **C IV far layer**: better targets or
  features.
  - Option: emulate the H, He and C ionisation fronts as column quantiles
    and tie each line to them. This is bigger; not started.
- **Observed-map validation:** run `civ_emission_map` / `line_escape` on
  both Cloudy and emulated rays and compare the ratio. Not done yet; only
  the all-direction maps have been compared.
- **h0.05 streaks:** not diagnosed.
- **Out of date since v5:** `emu redo --run h0.5_n30000`,
  `emu redo --run h0.05_n30000` and `emu emulate --val h0.5` all need
  rerunning.

## RGB maps across lamp heights (`emu rgb`, added 2026-11-02)

- `emu rgb` (default `--hs 0 0.05 0.2 0.5 1`, in ld) uses the dense
  emulated runs `emulated/h<h>_n30000`:
  - makes any that are missing;
  - re-emulates any older than `TABLE_VERSION` (`version` in
    `emulated_pred.npz`);
  - caches each run's face-on maps in `faceon_rgb.npz`.
- Output: `emulated/plots/rgb_heights.png` and
  `emulated/rgb_heights.npz`.
- Colours: R = Ha, G = Mg II, B = C IV. All-direction emission. One
  stretch shared by every panel (each line divided by its largest cell
  over all heights; hue over 1.5 dex, brightness over 4 dex, as in
  `faceon_map`).
- Heights in run names are in ld: h0.5 = 0.25 r0. The training runs
  cover 0-1 ld.
