# Note: line escape, turbulence and velocities in the per-ray Cloudy runs

Context: `transient/cloudy_rays.py`, one Cloudy c23.01 model per lamp ray
through the Athena++ snapshot (`disk.out1.00012.athdf`), compared with
MOCASSIN (`/data2/jhuang/runs/mocassin/moc_edd10`). Written 2026-10-17.

## Sim velocity structure (`cloudy_rays vmaps`)

- Command: `python3 -m transient.cloudy_rays vmaps --out
  /data2/jhuang/runs/cloudy_rays/edd10`; writes `plots/vel_{grad,shear,
  disp,normal}.png`, meridional cuts through the star (phi 335, index
  238), the arm (phi 70, index 49, strongest Sigma excess away from the
  star, 3.6x the phi-median at R = 1) and the quietest column (phi 199,
  index 141). Black contours: n_H = 1e9, 1e10, 1e11 cm^-3 (thin to
  thick). The lit skin is between the 1e10 and 1e11 contours.
- Fields (`velocity_fields`):
  - `grad`: dv_r/dr, along an in-plane lamp ray [km/s per 0.01 r0].
  - `shear`: |sigma| = sqrt(2 sigma_ij sigma_ij), traceless strain rate;
    dominated by Keplerian r dOmega/dr; excludes the divergence.
  - `disp`: 1D dispersion of (v_r, v_th, v_ph) over 3x3x3 cells [km/s].
    Not sub-grid turbulence: the resolved gradients over ~1e14 cm.
  - `normal`: n.e.n, gradient of the normal velocity along the surface
    normal n = -grad n_H / |grad n_H| (full strain rate, trace
    included). Red = layer stretched (outflow base), blue = squeezed.
- Grid: 128 r x 64 th x 256 ph; dr = 0.0117 r0, r dth = 0.005-0.02 r0,
  r dph = 0.012-0.05 r0. "Per 0.01 r0" is ~1 pixel at R ~ 1.
- Lit skin numbers (16/50/84%): |dv_r/dr| 0.3/2.8/26, |n.e.n| 0.15/1.4/21
  km/s per 0.01 r0; |n.e.n| x 1e12 cm = 0.003/0.03/0.4 km/s; skin normal
  is vertical (|n_r| ~ 0.04). Away from the star the skin is nearly
  hydrostatic (v_z ~ -1 +- 30 km/s). Above the star: squeezed (n.e.n ~
  -9). Envelope around the star: X pattern, rising (v_z ~ +14), T ~
  6500 K. Hot layer (n < 1e9): +-300 km/s streams; n is ill defined
  there (checkerboard), so only coherent patches mean anything.
- Conclusion: resolved velocities shift a photon by << the thermal width
  (Mg II 1.7, H 13 km/s) across a ~1e12 cm emitting layer; only the
  star envelope / arm reach ~0.5-2 km/s. The sim says nothing about
  sub-cell (< 1e14 cm) motions.

## Turbulence in Cloudy

- Stock c23.01: `turbulence` (constant), `vlaw alpha` (power law in
  radius, alpha must be < 0), `turbulence ... dissipate` (exp decay, adds
  wave heating). No table form. `wind` solves its own dynamics (LVG), so
  an imposed v(s) is impossible.
- A `turbulence table depth` patch would be easy: `vlaw` resets
  `DoppVel.TurbVel` in `ZoneStart` (`source/zone_startend.cpp:55`) and
  every line width comes from it via `GetDopplerWidth`. Dropped: the sim
  velocities give vturb ~ 0 anyway.
- The current 100 km/s is supersonic and not supported by the sim;
  sub-cell turbulence, if any, is plausibly subsonic (c_s ~ 10 km/s at
  1e4 K). MOCASSIN has no turbulence (thermal widths only).

## Line escape geometry in Cloudy (important)

- Escape is angle-integrated, in a plane-parallel (or spherical-shell)
  geometry whose layers are perpendicular to the lamp ray: P = [P(tau_in)
  + P(T - tau_in)] / 2 (`docs/latex/hazy2/linedetl.tex:238`), each term a
  K_2 function containing E_2 (exp(-tau/mu) integrated over all
  directions mu, `linedetl.tex:161`), out of the lamp face and the back
  face. A photon leaving at angle theta to the ray sees tau / cos theta,
  so the cheapest exits are along the ray. (Corrected 2026-10-18: earlier
  versions of this note called it "1D two-stream along the ray".)
- `save lines emissivity` (what the maps and spectra use) is the net
  escaping emission 4 pi J = n_u A_ul P_ul h nu (`hazy1/conout.tex`,
  "Save lines, emissivity"; intrinsic by default, `emergent` would add
  absorption on the way out), so the face-on maps already show escaping
  light, in Cloudy's geometry.
- Tilted skin: the lamp ray grazes the skin, whose normal is nearly
  vertical. The real shortest exit is through the skin surface, sideways
  to the ray, with column ~ tau sin(grazing angle) << tau; Cloudy's slab,
  with its normal along the ray, has no such path. For thick lines (Mg
  II) Cloudy's escape is then likely too LOW, which moves it toward
  MOCASSIN. The deeper direction's velocity gradient is dv_r/dr, not
  n.e.n; across the layer n.e.n gives 0.03-0.4 km/s, negligible next to
  vturb and the thermal widths, so velocity does not fix this.
- Deep route is closed physically: Sobolev length along the ray ~3e13 cm
  in the skin, Mg+ column ~1e20 over it, Mg II tau0 ~ 1e8 even then.
- OUR DECKS RUN ONE ITERATION (no `iterate`; Cloudy warns "This is the
  last iteration ... Use the ITERATE command"). On iteration 1, T is
  unknown and starts ~0, so tau_out = |T - tau_in| ~ tau_in
  (`tau_from_here`, `source/rt_escprob.cpp:47`): the deep side is as
  open as the lamp side. Only H Lya gets a guessed T (from `stop neutral
  column density`, `source/rt_tau_init.cpp:325`). Effect: thick lines
  (Mg II, Ha, C IV) escape up to 2x too easily.

## Horne & Marsh (1986, MNRAS 218, 761): shear and U vs V profiles

- Shear broadening across an emitting layer of thickness H:
  dv = dv_th [1 + Q^2 sin^2 i tan^2 i sin^2 2theta]^1/2,
  Q = (3/4)(V_K / V_th)(H/R) (eqs. 20-21); theta = disk azimuth from the
  line of sight; largest at theta = 45, 135 deg.
- Optically thin lines: isotropic, sharp peaks, shallow U-shaped valley.
  Optically thick lines: intensity ~ S_L dv (8 ln tau0)^1/2 (eq. 25),
  so they leak out along the high-shear 45 deg directions -> broader
  peaks, deeper V-shaped valley (once shear > thermal, sin i tan i >~ 1).
- Here (i = 45, sin i tan i = 0.71; R = 1, V_K ~ 5060 km/s; layer 1e12-
  1e13 cm): Q(Mg II) ~ 0.5-4.5 (dv up to ~3x thermal), Q(Ha) ~ 0.06-0.6
  (<~10%), Q with vturb = 100 km/s ~ 0.01-0.08 (no effect). Neither
  MOCASSIN (isotropic) nor our profile code includes this.
- The sim's real velocities belong here, in the profile synthesis, not
  in Cloudy: weight each emitting element by its line-of-sight width
  from |o.e.o| (strain tensor) x path through the layer; the arm and
  star then get their actual non-Keplerian shear.

## The Mg II gap (per-ray Cloudy ~100x below MOCASSIN)

- More turbulence makes thick-line flux larger (eq. 25), so vturb = 100
  already flatters Cloudy's Mg II; vturb ~ 0 makes it weaker. Turbulence
  cannot close the gap.
- Per-ray Cloudy Mg II/Ha ~ 0.4 is normal for AGN (Mg II/Hb ~ 1-2,
  Ha/Hb ~ 3-4). MOCASSIN gives ~4-5 (peak heights of its profiles, per-A
  units corrected). Its excess comes from the deep neutral midplane
  layer (|z| < 0.04, R 0.7-1.2, x(H+) ~ 0-0.01, 4000-8000 K), where
  MOCASSIN treats Mg II as optically thin (no trapping/thermalization).
- Working hypothesis: MOCASSIN, not Cloudy, is off for Mg II.

## Step 1 result: MOCASSIN line totals by x(H+) (2026-10-17)

- All-sky MOCASSIN totals (active cells): Mg II 2796+2804 = 3.0e42, Ha
  3.7e41, C IV 1.25e42 erg/s -> Mg II/Ha = 8.2 (per-ray Cloudy: 2.0e40,
  5.1e40, 4.1e41; Mg II/Ha = 0.4).
- Fraction kept when dropping cells below an x(H+) cut:

  | cut x(H+) >= | Mg II | Ha   | C IV | He II | Mg II/Ha |
  |--------------|-------|------|------|-------|----------|
  | 0.01         | 1.00  | 0.99 | 1.00 | 1.00  | 8.2      |
  | 0.1          | 0.51  | 0.76 | 1.00 | 0.98  | 5.5      |
  | 0.5          | 0.17  | 0.47 | 1.00 | 0.91  | 2.9      |

- So the truly neutral midplane (x < 0.01) emits ~nothing; half of
  MOCASSIN's Mg II is from 0.01 < x < 0.1 cells and a third from 0.1-0.5,
  i.e. partially ionized cells. Removing them does not close the gap
  (still ~75x at x >= 0.1).
- New suspicion besides optically thin Mg II: resolution. A MOCASSIN
  cell (~1e14 cm) straddling the ionization front is given one mean
  ionization, so its whole volume is "partially ionized" and emits Mg
  II, whereas in Cloudy the partially ionized zone is ~1e12 cm thick.
  This volume excess is of the order of the gap (cf. the N(C3+) layer
  being thinner than a cell).

## Plan (in order)

1. DONE (above): MOCASSIN line totals excluding cells with x(H+) < 0.1.
2. ~10 test rays (rim, ray 16, skin, star, arm, quiet) with `iterate`, at
   vturb = 0, 10, 100 km/s (~30 runs, user runs); compare with the
   1-iteration 100 km/s outputs; record tau0 per line.
   Decks: `python3 -m transient.cloudy_rays vtest --out
   /data2/jhuang/runs/cloudy_rays/vtest` (VTEST_RAYS: t13 and t31 at phi
   index 0, 48 arm, 140 quiet, 240 star; plus t16, t24 at 0; `iterate to
   convergence 4`; saves `.tau`). Run with the `run` step on that dir.
   RESULT (2026-10-17, all 30 exited OK; sum over the 10 rays relative to
   the 1-iteration 100 km/s production):

   | run                 | Ha   | Hb   | C IV | Mg II | Mg II/Ha |
   |---------------------|------|------|------|-------|----------|
   | production (1 it)   | 1    | 1    | 1    | 1     | 0.47     |
   | 100 km/s + iterate  | 0.91 | 0.90 | 0.92 | 0.87  | 0.45     |
   | 10 km/s + iterate   | 0.46 | 0.33 | 0.41 | 0.22  | 0.23     |
   | 0 km/s + iterate    | 0.44 | 0.30 | 0.32 | 0.13  | 0.14     |

   - `iterate` alone: lines down 3-18% (the deep side is now opaque).
   - Turbulence dominates: 100 -> 0 km/s cuts Mg II ~8x, C IV ~3x, Ha
     ~2x (thick-line flux ~ line width, Horne & Marsh eq. 25). He II
     1640 changes little (0.85-1.4 per ray).
   - tau0 at 0 / 100 km/s: Lya ~1e11 / 1e10, Ha ~1e6 / 1e5, C IV 1548
     ~2e7 / 1e6, Mg II 2796 ~1e8 / 2e6: all very thick.
   - Mg II/Ha falls to 0.14 at 0 km/s: further from MOCASSIN (8), so the
     gap is on the MOCASSIN side (optically thin Mg II, cell-size
     partially ionized volume), not Cloudy's turbulence.
   - Cost: every ray hit 4 iterations ("did not converge in 4", Cloudy's
     all-lines 5% criterion); CPU 3-5x production (t13 rays 1850-3560 s).
     The 4 rim rays t31 at phi 0/48/140/240 give near-identical results.
   - Decision (2026-10-18): production decks use `iterate 2` (DECK,
     ITERATE in cloudy_rays.py); 4 passes are unaffordable for ~3000
     rays. The vtest decks saved only the last pass (`print last`, save
     ... `last`), so the pass-by-pass history is not available; the
     `vtest --it2` check was skipped. Production vturb = 10 km/s
     (VTURB; subsonic ceiling for sub-cell motions).
   - Cloudy cannot take the sim's velocity field (checked c23.01 source):
     static/negative-velocity geometries add the plain cumulated opacity
     (`rt_line_one_tauinc.cpp:80`, no velocity shift); the only Sobolev /
     LVG path is the ballistic `wind`, whose dv/dr = |a_rad - g| / v
     comes from Cloudy's own momentum equation (`:99`), with density from
     continuity; `wind velocity` sets only the launch speed. A patch
     imposing a table dv/dr would still be a gradient along the lamp ray
     only (Sobolev tau there stays ~1e8) and isotropic. The
     real velocities go into the profile synthesis (step 3).
3. Horne & Marsh anisotropic weighting in `ray_spectra` for thick lines,
   using the sim strain tensor; needs the vturb ~ 0 outputs from step 2.
